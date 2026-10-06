#pragma once

// Global forward fit: generated-bin coefficients contribute to all reached
// reconstructed rows. Free overflow coefficients propagate their uncertainty
// through the full covariance in the default unconstrained fit. Optional
// positivity is imposed inside each fit, never by clipping fitted points.
#include "xsec_analysis.h"
#include "xsec_positive_solver.h"
#include "xsec_proxy_fit.h"
#include <set>
#include <TMatrixDSymEigen.h>

inline void ExclPi0XSecAnalysis::compute_ratios_and_xsec() {
    // Hydrogen yield reduction is external to MC and data fitting. Its common
    // uncertainty is exported separately, never added as independent row noise.
    const double factor = cfg.tgt_contam;
    for (size_t b = 0; b < slices.size(); ++b) {
        auto &s = slices[b];
        const auto &m = truth_moments[b];
        s.has_model_xsec = m.weight > 0;
        s.truth_response_sum = m.weight;
        if (s.has_model_xsec) {
            s.mean_q2_vertex_sim = m.q2 / m.weight;
            s.mean_xb_vertex_sim = m.xb / m.weight;
            s.mean_t_sim = m.t / m.weight;
            s.mean_tprime_vertex_sim = m.tprime / m.weight;
            s.epsilon = m.epsilon / m.weight;
            s.gamma_flux =
                virtual_photon_flux(cfg.ebeam, s.mean_q2_vertex_sim, s.mean_xb_vertex_sim, cfg.mp, s.epsilon);
        } else {
            s.mean_q2_vertex_sim = s.mean_xb_vertex_sim = s.mean_t_sim =
                s.mean_tprime_vertex_sim = s.epsilon = s.gamma_flux = std::numeric_limits<double>::quiet_NaN();
        }
        s.sumw_data /= factor;
        s.sumw2_data /= factor * factor;
        for (auto &p : s.phi) {
            p.data /= factor;
            p.data_sumw2 /= factor * factor;
            p.weights.divide(factor);
            p.data_plus /= factor;
            p.data_minus /= factor;
            p.data_plus_sumw2 /= factor * factor;
            p.data_minus_sumw2 /= factor * factor;
        }
    }
    for (auto& entry : run_weight_moments) entry.second.divide(factor);
}

inline void ExclPi0XSecAnalysis::fit_slices() {
    fit_attempts.clear();
    successful_fit_groups = 0;
    fit_fallback = false;
    const int ngroups = cfg.n_q2 * cfg.n_xb;
    retained_fit_groups.assign(ngroups, true);
    std::vector<std::string> exclusions(ngroups);
    auto attempt = [&](const std::vector<bool>& groups) {
        FitAttempt result;
        result.groups = groups;
        invalidate_fit_results();
        log("Global fit attempt " + std::to_string(fit_attempts.size() + 1) + ": " +
            std::to_string(std::count(groups.begin(), groups.end(), true)) + " Q2/xB groups retained");
        try {
            fit_global_subset(groups);
            result.ok = true;
        } catch (const std::runtime_error& error) {
            // A failed solve may already have filled some slice records.
            // Never publish a partially converged or partially populated fit.
            invalidate_fit_results();
            result.reason = error.what();
            migration_fit = nps_xsec::LinearSolution{};
            fit_curvature_inverse.clear();
            positivity_boundary_tolerances.clear();
            positivity_feasibility_tolerances.clear();
            positivity_boundary_active = false;
            positivity_iterations = 0;
            fit_rows.clear(); fit_variance.clear(); scaled_rows.clear(); mc_converged = false;
        }
        fit_attempts.push_back(std::move(result));
        if (!fit_attempts.back().ok)
            warn("Global fit attempt " + std::to_string(fit_attempts.size()) +
                 " failed: " + fit_attempts.back().reason);
        return fit_attempts.back().ok;
    };
    if (attempt(retained_fit_groups)) {
        successful_fit_groups = ngroups;
        return;
    }
    fit_fallback = true;
    invalidate_fit_results();
    // These failures identify a Q2/xB group unambiguously. Exclude its entire
    // reconstructed group (all t'/phi rows), not just the troublesome row.
    for (size_t b = 0; b < slices.size(); ++b) {
        const size_t g = b % ngroups;
        if (!(truth_moments[b].weight > 0))
            exclusions[g] = "No response for generated bin " + std::to_string(b);
        for (int ip = 0; ip < cfg.n_phi; ++ip) {
            const size_t r = b * cfg.n_phi + ip;
            bool support = false;
            for (const auto& cell : migration_response[r])
                support = support || std::any_of(cell.basis.begin(), cell.basis.end(), [](double v) { return v != 0; });
            const auto& p = slices[b].phi[ip];
            if (!support && (p.data != 0 || p.data_sumw2 > 0) && exclusions[g].empty())
                exclusions[g] = "Data outside MC support in row " + std::to_string(r) +
                                " (it=" + std::to_string(b / ngroups) + ", ip=" + std::to_string(ip) + ")";
        }
        if (!exclusions[g].empty()) retained_fit_groups[g] = false;
    }
    std::vector<int> candidates;
    for (int g = 0; g < ngroups; ++g)
        if (retained_fit_groups[g]) candidates.push_back(g);

    // Rank/convergence failures need not have a unique offending bin. Try
    // global subsets in decreasing size; ties use increasing excluded index.
    // Never select on chi2 or silently regularize an unidentifiable problem.
    bool recovered = false;
    if (static_cast<int>(candidates.size()) < ngroups && !candidates.empty())
        recovered = attempt(retained_fit_groups);
    for (size_t drop = 1; !recovered && drop < candidates.size(); ++drop) {
        auto search = [&](auto&& self, size_t start, size_t left) -> bool {
            if (left == 0) return attempt(retained_fit_groups);
            for (size_t i = start; i + left <= candidates.size(); ++i) {
                retained_fit_groups[candidates[i]] = false;
                if (self(self, i + 1, left - 1)) return true;
                retained_fit_groups[candidates[i]] = true;
            }
            return false;
        };
        recovered = search(search, 0, drop);
    }
    if (!recovered) {
        retained_fit_groups.assign(ngroups, false);
        active_truth_blocks.clear();
        response_design.assign(migration_response.size(), {});
        nonphysical_truth_bins = 0;
    }
    successful_fit_groups = static_cast<int>(std::count(retained_fit_groups.begin(), retained_fit_groups.end(), true));
    for (int g = 0; g < ngroups; ++g) {
        if (retained_fit_groups[g]) continue;
        if (exclusions[g].empty())
            exclusions[g] = recovered ? "Excluded to recover global fit; see fit_attempts.csv" :
                                       "No viable global subset: " + fit_attempts.back().reason;
        for (size_t b = g; b < slices.size(); b += ngroups)
            slices[b].fit_failure_reason = exclusions[g];
        warn("Excluded Q2/xB bin iq=" + std::to_string(g / cfg.n_xb) + ", ix=" +
             std::to_string(g % cfg.n_xb) + ": " + exclusions[g]);
    }
    warn("Global fit summary: " + std::to_string(successful_fit_groups) + " Q2/xB bins retained, " +
         std::to_string(ngroups - successful_fit_groups) + " excluded; see fit_status.csv");
}

inline void ExclPi0XSecAnalysis::invalidate_fit_results() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (size_t b = 0; b < slices.size(); ++b) {
        auto& s = slices[b];
        s.fit_scope = fit_fallback ? "global_migration_subset" : "global_migration";
        if (cfg.fit_objective == "scaled-poisson") s.fit_scope += "_scaled_poisson";
        if (model_fit_mode) s.fit_scope += "_sigparam2021_pi0";
        if (cfg.positive_xsec) s.fit_scope += "_positive";
        s.fit_failure_reason.clear();
        s.fit_rank = s.fit_mc_iterations = 0;
        s.fit_condition = nan;
        auto& f = s.fit_xsec;
        f.ok = f.absolute_xsec_fit = false;
        f.p.clear(); f.perr.clear();
        f.chi2 = f.ndf = nan;
        f.sigmaU = f.sigmaU_err = f.sigmaTL = f.sigmaTL_err = f.sigmaTT = f.sigmaTT_err = nan;
        for (auto& p : s.phi) {
            p.sim = p.sim_sumw2 = p.ratio = p.ratio_err = nan;
            p.xsec = p.xsec_err = p.xsec_sys_tgt = nan;
            p.mean_q2_xsec = p.mean_xb_xsec = p.mean_tprime_xsec = nan;
        }
    }
}

inline void ExclPi0XSecAnalysis::fit_global_subset(const std::vector<bool>& groups) {
    active_truth_blocks.clear(); fixed_truth_blocks.clear(); fit_rows.clear(); fit_variance.clear(); scaled_rows.clear();
    scaled_minuit_status=scaled_covariance_status=-1;
    scaled_edm=std::numeric_limits<double>::quiet_NaN();
    scaled_calls=0;
    migration_fit = nps_xsec::LinearSolution{};
    fit_curvature_inverse.clear();
    positivity_boundary_tolerances.clear();
    positivity_feasibility_tolerances.clear();
    positivity_boundary_active = false;
    positivity_iterations = 0;
    mc_iterations = omitted_zero_variance_rows = nonphysical_truth_bins = 0;
    mc_converged = false;
    const auto selected = [&](size_t b) {
        return groups[b % (cfg.n_q2 * cfg.n_xb)];
    };
    // All published bins must be identifiable. Only low-tprime exterior
    // feed-in receives fitted columns; other exterior events remain fixed.
    for (size_t b = 0; b < truth_moments.size(); ++b) {
        bool contributes = false;
        for (size_t r = 0; r < migration_response.size(); ++r)
            if (selected(r / cfg.n_phi)) {
                const auto& basis = migration_response[r][b].basis;
                contributes = contributes || std::any_of(basis.begin(), basis.end(), [](double v) { return v != 0; });
            }
        if ((b < slices.size() && selected(b)) || contributes) {
            if(!event_model() && nps_xsec::is_fixed_model_feedin(static_cast<int>(b),static_cast<int>(slices.size())))
                fixed_truth_blocks.push_back(static_cast<int>(b));
            else active_truth_blocks.push_back(static_cast<int>(b));
        }
    }
    const size_t npar = 3 * active_truth_blocks.size();
    response_design.assign(migration_response.size(), std::vector<double>(npar, 0));
    if (cfg.fit_objective == "scaled-poisson") {
        for (size_t r=0;r<migration_response.size();++r)
            for (size_t b=0;b<active_truth_blocks.size();++b)
                for (int a=0;a<3;++a)
                    response_design[r][3*b+a]=migration_response[r][active_truth_blocks[b]].basis[a];
        fit_scaled_poisson_subset(groups);
        return;
    }
    std::vector<std::vector<double>> X;
    std::vector<double> data, data_variance;
    for (size_t r = 0; r < migration_response.size(); ++r) {
        auto &row = response_design[r];
        for (size_t b = 0; b < active_truth_blocks.size(); ++b)
            for (int a = 0; a < 3; ++a)
                row[3 * b + a] = migration_response[r][active_truth_blocks[b]].basis[a];
        if (!selected(r / cfg.n_phi)) continue;
        const auto &p = slices[r / cfg.n_phi].phi[r % cfg.n_phi];
        const bool support = std::any_of(row.begin(), row.end(), [](double v) { return v != 0; }) ||
                             (!event_model() && fixed_feedin_prediction[r]!=0.);
        if (!support && (p.data != 0 || p.data_sumw2 > 0))
            die("Data outside MC support in row " + std::to_string(r) +
                " (it=" + std::to_string(r / cfg.n_phi / (cfg.n_q2 * cfg.n_xb)) +
                ", iq=" + std::to_string((r / cfg.n_phi / cfg.n_xb) % cfg.n_q2) +
                ", ix=" + std::to_string((r / cfg.n_phi) % cfg.n_xb) +
                ", ip=" + std::to_string(r % cfg.n_phi) + ")");
        if (p.data_sumw2 <= 0) {
            // Weighted Gaussian fits cannot determine an observed variance in
            // empty rows. Do not invent pseudocounts: export and flag exclusions.
            // Sparse-bin bias requires independent toys/coarser-bin validation.
            if (support)
                ++omitted_zero_variance_rows;
            continue;
        }
        if (!support)
            continue;
        X.push_back(row);
        data.push_back(p.data-(event_model()?0.:fixed_feedin_prediction[r]));
        data_variance.push_back(p.data_sumw2);
        fit_rows.push_back(static_cast<int>(r));
    }
    if (omitted_zero_variance_rows)
        warn(std::to_string(omitted_zero_variance_rows) +
             " MC-supported rows have zero observed variance; excluded and recorded");

    if (model_fit_mode) { fit_proxy_subset(groups); return; }

    // Every fitted block is constrained, including low-tprime feed-in and
    // excluded published bins retained as independent parameters. At fixed
    // binwise U/LT/TT the angular minimum decreases with epsilon, so checking
    // the largest response-event epsilon covers its full observed envelope.
    std::vector<double> epsilon_max;
    for (int b : active_truth_blocks) epsilon_max.push_back(truth_moments[b].epsilon_max);
    const auto solve = [&](const std::vector<double>& variance) {
        if (!cfg.positive_xsec)
            return nps_xsec::solve_weighted_response(X, data, variance, cfg.rank_tolerance);
        auto constrained = nps_xsec::solve_positive_response(X, data, variance, epsilon_max, cfg.rank_tolerance);
        positivity_boundary_active = constrained.boundary_active;
        positivity_boundary_tolerances = constrained.boundary_tolerances;
        positivity_feasibility_tolerances = constrained.feasibility_tolerances;
        positivity_iterations = constrained.iterations;
        return constrained.fit;
    };
    // Feasible GLS recomputes finite-MC variance from the constrained trial
    // coefficients as well. Constraining only the final unconstrained fit
    // would leave those response variances inconsistent with the estimate.
    fit_variance = data_variance;
    migration_fit = solve(fit_variance);
    mc_converged = cfg.fit_variance_mode == "data";
    for (mc_iterations = 0; cfg.fit_variance_mode == "finite-mc" && mc_iterations < cfg.mc_max_iterations;) {
        ++mc_iterations;
        std::vector<double> next_variance = data_variance;
        for (size_t i = 0; i < fit_rows.size(); ++i)
            next_variance[i] += nps_xsec::mc_prediction_variance(
                migration_response[fit_rows[i]], active_truth_blocks, migration_fit.parameters)+
                fixed_feedin_mc_variance[fit_rows[i]];
        auto next = solve(next_variance);
        double parameter_change = 0, variance_change = 0;
        // Error-scaled convergence avoids relative division by a coefficient
        // crossing zero. Both coefficients and row variances must settle.
        for (size_t j = 0; j < npar; ++j) {
            const double scale =
                std::max(std::abs(next.parameters[j]), std::sqrt(next.covariance[j * npar + j]));
            parameter_change = std::max(parameter_change,
                                        std::abs(next.parameters[j] - migration_fit.parameters[j]) / scale);
        }
        for (size_t i = 0; i < fit_rows.size(); ++i)
            variance_change =
                std::max(variance_change, std::abs(next_variance[i] - fit_variance[i]) / next_variance[i]);
        migration_fit = std::move(next);
        fit_variance = std::move(next_variance);
        if (std::max(parameter_change, variance_change) < cfg.mc_fit_tolerance) {
            mc_converged = true;
            break;
        }
    }
    if (!mc_converged)
        die("Finite-MC iteration did not converge; inspect MC statistics or increase --mc-max-iterations");

    fit_curvature_inverse = migration_fit.covariance;
    if (auto* statistics=dynamic_cast<TTree*>(f_data->Get("analysis_sigma_covariance"))) {
        if(cfg.fit_objective!="gaussian" || cfg.fit_variance_mode!="data" || cfg.positive_xsec ||
           !f_data->Get("analysis_reco_yields"))
            die("External event-bootstrap covariance requires the matching unconstrained data-only point estimator");
        int i=0,j=0,truth_i=0,truth_j=0;
        double cdata=0,cmc=0,nominal_i=0,nominal_j=0;
        bind_branch(statistics,"i",&i);bind_branch(statistics,"j",&j);
        bind_branch(statistics,"truth_i",&truth_i);bind_branch(statistics,"truth_j",&truth_j);
        bind_branch(statistics,"data_covariance",&cdata);bind_branch(statistics,"mc_covariance",&cmc);
        bind_branch(statistics,"nominal_i",&nominal_i);bind_branch(statistics,"nominal_j",&nominal_j);
        const size_t np=migration_fit.parameters.size();
        if(statistics->GetEntries()!=static_cast<Long64_t>(np*np)) die("Bootstrap covariance dimension mismatch");
        std::set<std::pair<int,int>> seen;
        std::vector<double> covariance(np*np);
        const auto same_point=[](double a,double b) { return std::abs(a-b)<=1e-10*std::max(std::abs(a),std::abs(b))+1e-20; };
        for(Long64_t row=0;row<statistics->GetEntries();++row) {
            statistics->GetEntry(row);
            if(i<0 || j<0 || i>=static_cast<int>(np) || j>=static_cast<int>(np) || !seen.emplace(i,j).second ||
               truth_i!=active_truth_blocks[i/3] || truth_j!=active_truth_blocks[j/3] ||
               !same_point(nominal_i,migration_fit.parameters[i]) || !same_point(nominal_j,migration_fit.parameters[j]) ||
               !std::isfinite(cdata) || !std::isfinite(cmc) || (i==j && (cdata<0 || cmc<0)))
                die("Bootstrap covariance does not match this nominal estimator");
            covariance[i*np+j]=cdata+cmc;
        }
        statistics->ResetBranchAddresses();
        for(size_t i=0;i<np;++i) for(size_t j=0;j<np;++j)
            if(std::abs(covariance[i*np+j]-covariance[j*np+i])>
               1e-10*std::sqrt(covariance[i*np+i]*covariance[j*np+j])+1e-30)
                die("Asymmetric bootstrap covariance");
        TMatrixDSym correlation(np);
        for(size_t i=0;i<np;++i) {
            if(!(covariance[i*np+i]>0)) die("Bootstrap covariance has no variance for an active coefficient");
            for(size_t j=0;j<np;++j)
                correlation(i,j)=covariance[i*np+j]/std::sqrt(covariance[i*np+i]*covariance[j*np+j]);
        }
        const auto eigenvalues=TMatrixDSymEigen(correlation).GetEigenValues();
        for(size_t i=0;i<np;++i) if(eigenvalues[i]<-1e-10*np) die("Bootstrap covariance is not positive semidefinite");
        migration_fit.covariance=std::move(covariance);
        log("Reported covariance: full event data bootstrap plus supplied independent finite-MC estimate; nominal coefficients unchanged");
    }
    if (positivity_boundary_active) {
        // A boundary changes the sampling distribution and invalidates the
        // unconstrained hat-matrix prediction errors. Do not report zero or
        // symmetric Gaussian errors as confidence intervals. NaN means not
        // evaluated; constrained profiles/toys are a separate calculation.
        std::fill(migration_fit.covariance.begin(), migration_fit.covariance.end(),
                  std::numeric_limits<double>::quiet_NaN());
        warn("Positivity boundary active: ordinary covariance/errors unavailable; "
             "conditional refit-toy plot spread is computed separately. "
             "Unconstrained inverse curvature is a diagnostic only.");
    }

    finalize_fit_subset(groups);
}

inline void ExclPi0XSecAnalysis::finalize_fit_subset(const std::vector<bool>& groups) {
    const size_t npar=migration_fit.parameters.size();
    const auto selected=[&](size_t b) { return groups[b%(cfg.n_q2*cfg.n_xb)]; };
    // Per-bin records expose marginal covariance blocks. Complete cross-bin
    // and nuisance covariance is separately saved for downstream fitting.
    for (size_t b = 0; b < slices.size(); ++b) {
        if (!selected(b)) continue;
        const size_t column = 3 * static_cast<size_t>(
            std::find(active_truth_blocks.begin(), active_truth_blocks.end(), static_cast<int>(b)) - active_truth_blocks.begin());
        auto &s = slices[b];
        s.fit_rank = static_cast<int>(migration_fit.rank);
        s.fit_condition = migration_fit.condition;
        s.fit_mc_iterations = mc_iterations;
        auto &f = s.fit_xsec;
        f.ok = true;
        f.absolute_xsec_fit = true;
        f.p.resize(3);
        f.perr.resize(3);
        f.cov.ResizeTo(3, 3);
        f.chi2 = migration_fit.chi2;
        f.ndf = migration_fit.ndf; // global quantities
        for (int i = 0; i < 3; ++i) {
            f.p[i] = migration_fit.parameters[column + i];
            f.perr[i] = std::sqrt(migration_fit.covariance[(column + i) * npar + column + i]);
            for (int j = 0; j < 3; ++j)
                f.cov(i, j) = migration_fit.covariance[(column + i) * npar + column + j];
        }
        f.sigmaU = f.p[0];
        f.sigmaU_err = f.perr[0];
        f.sigmaTL = f.p[1];
        f.sigmaTL_err = f.perr[1];
        f.sigmaTT = f.p[2];
        f.sigmaTT_err = f.perr[2];
        const double check_epsilon = cfg.positive_xsec ? truth_moments[b].epsilon_max : s.epsilon;
        const double positivity_roundoff = cfg.positive_xsec ? positivity_feasibility_tolerances[column/3] : 0.0;
        if (nps_xsec::minimum_response(f.p[0], f.p[1], f.p[2], check_epsilon) < -positivity_roundoff)
            ++nonphysical_truth_bins;
    }
    for (size_t r = 0; r < response_design.size(); ++r) {
        if (!selected(r / cfg.n_phi)) continue;
        auto &s = slices[r / cfg.n_phi];
        auto &p = s.phi[r % cfg.n_phi];
        const auto &row = response_design[r];
        p.sim = fixed_feedin_prediction[r]+std::inner_product(row.begin(), row.end(), migration_fit.parameters.begin(), 0.0);
        double variance = 0;
        for (size_t i = 0; i < npar; ++i)
            for (size_t j = 0; j < npar; ++j)
                variance += row[i] * migration_fit.covariance[i * npar + j] * row[j];
        if(event_model()){p.sim=model_rows[r].prediction;variance=model_row_variance(r);}
        // The response MC also determined the fitted coefficients. Its effect
        // on this same-fit prediction is therefore correlated with the fit.
        // With H = X (X' V^-1 X)^-1 X' V^-1, first-order propagation gives
        // Var(pred_r) = (X C X')_rr + Vmc_r - 2 H_rr Vmc_r.
        // For an excluded row H_rr=0: its MC is independent of fitted rows.
        // This is conditional on final GLS weights; it is not a posterior
        // predictive uncertainty or a treatment of systematic correlations.
        const double mc_variance = (cfg.fit_objective == "scaled-poisson" || cfg.fit_variance_mode == "data") ? 0.0 : event_model()?model_rows[r].mc_variance:
            nps_xsec::mc_prediction_variance(
                migration_response[r], active_truth_blocks, migration_fit.parameters)+fixed_feedin_mc_variance[r];
        const auto fitted = std::find(fit_rows.begin(), fit_rows.end(), static_cast<int>(r));
        const double leverage = (cfg.fit_objective == "scaled-poisson" || fitted == fit_rows.end())
                                    ? 0
                                    : variance / fit_variance[static_cast<size_t>(fitted - fit_rows.begin())];
        const double prediction_variance = variance + mc_variance * (1 - 2 * leverage);
        if (!positivity_boundary_active &&
            (!std::isfinite(prediction_variance) || prediction_variance < -1e-10 * (variance + mc_variance)))
            die("Invalid fitted-prediction variance in reconstructed row " + std::to_string(r));
        p.sim_sumw2 = positivity_boundary_active ? std::numeric_limits<double>::quiet_NaN() :
                         std::max(0.0, prediction_variance);
        // Fitted same-data yield ratio is a residual diagnostic, not an
        // independent acceptance correction or normalization measurement.
        const double undefined = std::numeric_limits<double>::quiet_NaN();
        p.ratio = p.sim > 0 ? p.data / p.sim : undefined;
        p.ratio_err = p.sim > 0 ? std::sqrt(p.data_sumw2) / p.sim : undefined;
        const int ip = static_cast<int>(r % cfg.n_phi);
        const auto basis =
            sigma_model_basis(.5 * (phi_edges[ip] + phi_edges[ip + 1]), s.epsilon, s.gamma_flux, 0);
        p.xsec = 0;
        double xvariance = 0;
        for (int i = 0; i < 3; ++i) {
            p.xsec += basis[i] * s.fit_xsec.p[i];
            for (int j = 0; j < 3; ++j)
                xvariance += basis[i] * s.fit_xsec.cov(i, j) * basis[j];
        }
        p.xsec_err = positivity_boundary_active ? std::numeric_limits<double>::quiet_NaN() :
                        std::sqrt(std::max(0.0, xvariance));
        p.xsec_sys_tgt = std::abs(p.xsec) * cfg.tgt_contam_err / cfg.tgt_contam;
        p.mean_q2_xsec = s.mean_q2_vertex_sim;
        p.mean_xb_xsec = s.mean_xb_vertex_sim;
        p.mean_tprime_xsec = s.mean_tprime_vertex_sim;
    }
    if (nonphysical_truth_bins)
        warn(std::to_string(nonphysical_truth_bins) +
             " generated bins have negative fitted cross section at some phi; inspect covariance/closure, no "
             "clipping applied");
    std::cout << "Global migration fit: rows=" << fit_rows.size() << ", parameters=" << (model_fit_mode ? proxy_result.parameters.size() : npar)
              << ", rank=" << migration_fit.rank << ", condition=" << migration_fit.condition
              << ", " << (cfg.fit_objective == "scaled-poisson" ? "deviance/nominal_ndf=" : "chi2/ndf=")
              << migration_fit.chi2 << "/" << migration_fit.ndf
              << ", MC iterations=" << mc_iterations << std::endl;
    if (cfg.fit_objective=="scaled-poisson") {
        for (size_t b=0;b<slices.size();++b) if (slices[b].fit_xsec.ok)
            std::cout << "Truth block " << b << ": U=" << slices[b].fit_xsec.sigmaU
                      << ", LT=" << slices[b].fit_xsec.sigmaTL
                      << ", TT=" << slices[b].fit_xsec.sigmaTT << std::endl;
    }
    if (cfg.fit_objective=="scaled-poisson")
        std::cout << "Minuit2: status=" << scaled_minuit_status
                  << ", EDM=" << scaled_edm << ", calls=" << scaled_calls
                  << ", covariance_status=" << scaled_covariance_status
                  << ", response=fixed, weight_correction=fixed" << std::endl;
    if (cfg.positive_xsec)
        std::cout << "Positivity: full phi, all active truth blocks, observed epsilon envelope; "
                  << "boundary=" << positivity_boundary_active << ", final solver iterations="
                  << positivity_iterations << "; nominal ndf is descriptive only" << std::endl;
}
