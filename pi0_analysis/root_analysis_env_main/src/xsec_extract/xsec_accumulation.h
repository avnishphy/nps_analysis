#pragma once

// Book owned histograms, read event trees, accumulate weighted data and Monte Carlo response.
#include "xsec_analysis.h"
#include "xsec_vertex_epsilon.h"

inline void ExclPi0XSecAnalysis::init_storage() {
    slices.assign(cfg.n_tprime * cfg.n_q2 * cfg.n_xb, SliceResult{});
    truth_moments.resize(slices.size()+6);
    migration_response.assign(slices.size()*cfg.n_phi,
        std::vector<nps_xsec::ResponseCell>(truth_moments.size()));
    truth_phi_response.resize(migration_response.size());
    truth_phi_moments.resize(truth_moments.size()*cfg.n_phi);
    for (auto& s : slices) s.phi.resize(cfg.n_phi);
    if (cfg.diagnostics) {
        const double qspan = cfg.q2_max - cfg.q2_min;
        const double xspan = cfg.xb_max - cfg.xb_min;
        auto coverage_hist = [&](const char* name, const char* title) {
            auto h = std::make_unique<TH2D>(name, title, 100,
                std::max(0.0, cfg.xb_min - 0.5 * xspan), std::min(1.0, cfg.xb_max + 0.5 * xspan),
                100, std::max(0.0, cfg.q2_min - 0.5 * qspan), cfg.q2_max + 0.5 * qspan);
            h->SetDirectory(nullptr);
            // Outside-range truth can feed selected reconstructed rows. Extend
            // diagnostic axes to retain it instead of hiding it in overflow.
            h->SetCanExtend(TH1::kAllAxes);
            return h;
        };
        h_migration_vertex_q2_xb = coverage_hist("vertex_q2_xb_selected",
            "Generated coordinates of selected MC;x_{B,gen};Q^{2}_{gen} [GeV^{2}];Selected MC entries");
        h_migration_reco_q2_xb = coverage_hist("reco_q2_xb_selected",
            "Reconstructed coordinates of the same MC;x_{B,reco};Q^{2}_{reco} [GeV^{2}];Selected MC entries");
    }
    if (cfg.mmiss_select != "window") {
        auto book_mass = [](const char* name) {
            // Use coarse display bins: sparse selected data remain visible
            // while the MCD/ellipse boundary is drawn at full precision.
            auto h = std::make_unique<TH2D>(name,
                ";m_{#pi^{0}} [GeV];M_{miss} [GeV];Events",
                60, 0.10, 0.17, 80, 0.5, 1.5);
            h->SetDirectory(nullptr);
            return h;
        };
        h_mass_data_all = book_mass("h_mass_data_all");
        h_mass_data_selected = book_mass("h_mass_data_selected");
        h_mass_sim_all = book_mass("h_mass_sim_all");
        h_mass_sim_selected = book_mass("h_mass_sim_selected");
        // Corrected Mx follows the workflow's unweighted event-count convention.
        h_mmiss_corr_data_all = std::make_unique<TH1D>(
            "h_mmiss_corr_data_all", ";mmiss_all_corr [GeV];Events / 20 MeV",
            125, 0.0, 2.5);
        h_mmiss_corr_data_selected = std::make_unique<TH1D>(
            "h_mmiss_corr_data_selected", ";mmiss_all_corr [GeV];Events / 20 MeV",
            125, 0.0, 2.5);
        h_mmiss_corr_data_all->SetDirectory(nullptr);
        h_mmiss_corr_data_selected->SetDirectory(nullptr);
        h_mmiss_reco_data_all = std::make_unique<TH1D>(
            "h_mmiss_reco_data_all", ";mmiss_all [GeV];Events / 20 MeV",
            125, 0.0, 2.5);
        h_mmiss_reco_data_selected = std::make_unique<TH1D>(
            "h_mmiss_reco_data_selected", ";mmiss_all [GeV];Events / 20 MeV",
            125, 0.0, 2.5);
        h_mmiss_reco_data_all->SetDirectory(nullptr);
        h_mmiss_reco_data_selected->SetDirectory(nullptr);
        h_mmiss_sim_all = std::make_unique<TH1D>(
            "h_mmiss_sim_all", ";SIMC mmiss [GeV];Events / 20 MeV",
            125, 0.0, 2.5);
        h_mmiss_sim_selected = std::make_unique<TH1D>(
            "h_mmiss_sim_selected", ";SIMC mmiss [GeV];Events / 20 MeV",
            125, 0.0, 2.5);
        h_mmiss_sim_all->SetDirectory(nullptr);
        h_mmiss_sim_selected->SetDirectory(nullptr);
    }
    mmiss_slices.resize(slices.size());
    for (size_t i = 0; i < mmiss_slices.size(); ++i) {
        auto create_mass_hist = [&](const std::string& channel) {
            const std::string name = "h_mmiss_" + channel + "_slice" + std::to_string(i);
            auto h = std::make_unique<TH1D>(name.c_str(),
                ";Reconstructed missing mass [GeV];Weighted yield / 25 MeV", 100, 0.0, 2.5);
            h->SetDirectory(nullptr);
            h->Sumw2();
            return h;
        };
        mmiss_slices[i].data = create_mass_hist("data");
        mmiss_slices[i].exclusive = create_mass_hist("exclusive");
    }

    h_q2_data.reset(new TH1D("h_q2_data", "Data;Q^{2} [GeV^{2}];Weighted counts", 100, cfg.q2_min, cfg.q2_max));
    h_q2_sim.reset(new TH1D("h_q2_sim", "SIMC;Q^{2} [GeV^{2}];Weighted counts", 100, cfg.q2_min, cfg.q2_max));
    h_xb_data.reset(new TH1D("h_xb_data", "Data;x_{B};Weighted counts", 100, cfg.xb_min, cfg.xb_max));
    h_xb_sim.reset(new TH1D("h_xb_sim", "SIMC;x_{B};Weighted counts", 100, cfg.xb_min, cfg.xb_max));
    h_tprime_data.reset(new TH1D("h_tprime_data", "Data;t' [GeV^{2}];Weighted counts", 120, cfg.tprime_min, cfg.tprime_max));
    h_tprime_sim.reset(new TH1D("h_tprime_sim", "SIMC;t' [GeV^{2}];Weighted counts", 120, cfg.tprime_min, cfg.tprime_max));
    h_phi_data.reset(new TH1D("h_phi_data", "Data;#phi [rad];Weighted counts", cfg.n_phi, phi_edges.data()));
    h_phi_sim.reset(new TH1D("h_phi_sim", "SIMC;#phi [rad];Weighted counts", cfg.n_phi, phi_edges.data()));
    h_q2_xb_data.reset(new TH2D("h_q2_xb_data", "Data;Q^{2} [GeV^{2}];x_{B}", 100, cfg.q2_min, cfg.q2_max, 100, cfg.xb_min, cfg.xb_max));
    h_q2_xb_sim.reset(new TH2D("h_q2_xb_sim", "SIMC;Q^{2} [GeV^{2}];x_{B}", 100, cfg.q2_min, cfg.q2_max, 100, cfg.xb_min, cfg.xb_max));
    h_tprime_phi_data.reset(new TH2D("h_tprime_phi_data", "Data;t' [GeV^{2}];#phi [rad]", 100, cfg.tprime_min, cfg.tprime_max, cfg.n_phi, phi_edges.data()));
    h_tprime_phi_sim.reset(new TH2D("h_tprime_phi_sim", "SIMC;t' [GeV^{2}];#phi [rad]", 100, cfg.tprime_min, cfg.tprime_max, cfg.n_phi, phi_edges.data()));
}

inline void ExclPi0XSecAnalysis::fill_mmiss_data_diagnostic(double q2, double t, double tmin,
        double xb, double phi, double mmiss, double pi0_weight, float scale,
        double charge_uC, double total_charge_uC) {
    // Do not apply the cross-section's high-side mass cut. Full reconstructed
    // range must remain visible when locating data/exclusive-SIMC divergence.
    const double tp = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tp) || !std::isfinite(phi) ||
        !std::isfinite(mmiss) || mmiss <= 0.0 || mmiss >= 2.5 ||
        !std::isfinite(pi0_weight) || !std::isfinite(scale)) return;
    const int it = find_bin(tprime_edges, tp, false);
    const int iq = find_bin(q2_edges, q2, false);
    if (it < 0 || iq < 0 || iq >= cfg.n_q2 || find_bin(phi_edges, phi, true) < 0) return;
    const int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], xb, false);
    if (ix < 0) return;
    double charge_fraction = 1.0;
    if (std::isfinite(charge_uC) && charge_uC > 0.0 &&
        std::isfinite(total_charge_uC) && total_charge_uC > 0.0)
        charge_fraction = charge_uC / total_charge_uC;
    const double w = pi0_weight * static_cast<double>(scale) * charge_fraction;
    mmiss_slices[static_cast<size_t>(slice_index(it, iq, ix))].data->Fill(mmiss, w);
}

inline void ExclPi0XSecAnalysis::fill_mmiss_sim_diagnostic(float q2, float t, float tmin,
        float xb, float phi, float mmiss, float full_weight,
        int exclusive) {
    // Use generated-exclusive events only. full_weight retains sigcm and the
    // production normalization, so it describes the reconstructed exclusive
    // yield shape; full_weight/sigcm is reserved for acceptance demodelling.
    if (exclusive == 0) return;
    const double tp = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tp) || !std::isfinite(phi) ||
        !std::isfinite(mmiss) || mmiss <= 0.0f || mmiss >= 2.5f ||
        !std::isfinite(full_weight)) return;
    const int it = find_bin(tprime_edges, tp, false);
    const int iq = find_bin(q2_edges, q2, false);
    if (it < 0 || iq < 0 || iq >= cfg.n_q2 || find_bin(phi_edges, phi, true) < 0) return;
    const int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], xb, false);
    if (ix < 0) return;
    mmiss_slices[static_cast<size_t>(slice_index(it, iq, ix))].exclusive->Fill(mmiss, full_weight);
}

inline void ExclPi0XSecAnalysis::accumulate_global_histograms(float q2, float xb, double tprime, double phi, double weight_data, double weight_sim) {
    h_q2_data->Fill(q2, weight_data);
    h_xb_data->Fill(xb, weight_data);
    h_tprime_data->Fill(tprime, weight_data);
    h_phi_data->Fill(phi, weight_data);
    h_q2_xb_data->Fill(q2, xb, weight_data);
    h_tprime_phi_data->Fill(tprime, phi, weight_data);

    h_q2_sim->Fill(q2, weight_sim);
    h_xb_sim->Fill(xb, weight_sim);
    h_tprime_sim->Fill(tprime, weight_sim);
    h_phi_sim->Fill(phi, weight_sim);
    h_q2_xb_sim->Fill(q2, xb, weight_sim);
    h_tprime_phi_sim->Fill(tprime, phi, weight_sim);
}

inline void ExclPi0XSecAnalysis::fill_data_event(double q2, double t, double tmin, double xb, double phi, double pi0_weight, float scale, double charge_uC, double total_charge_uC, double mmiss_all, double mpi0_all, int exclusive_flag, int helicity, bool use_helicity, double W, int run_number) {
    cutflow.n_data_total++;
    if (!selected_data(mmiss_all, mpi0_all, exclusive_flag)) return;
    double tprime = calc_tprime(t, tmin);
    if (cfg.prepare_forward_inputs) {
        if (!xsec_inside_diamond(cfg, xb, q2)) return;
        const double w = pi0_weight * static_cast<double>(scale) * charge_uC /
                         total_charge_uC / cfg.tgt_contam;
        for (double v : {q2, xb, tprime, phi, w})
            if (!std::isfinite(v)) die("Nonfinite selected forward data event");
        if (w < 0) die("Forward scaled-Poisson input requires nonnegative event weights; signed background subtraction needs a separate likelihood.");
        forward_data_stream << forward_event_id << ',' << run_number << ',' << q2 << ',' << xb
                            << ',' << tprime << ',' << wrap_phi(phi) << ',' << w << '\n';
        ++forward_data_count;
        return;
    }
    if (!slice_passes_kin(q2, xb, tprime)) return;
    if (!std::isfinite(phi)) return;
    if (!std::isfinite(pi0_weight) || !std::isfinite(scale)) die("Nonfinite selected data weight/scale");

    cutflow.n_data_pass++;

    int it = find_bin(tprime_edges, tprime, false);
    int iq = find_bin(q2_edges, q2, false);
    if (iq < 0 || iq >= cfg.n_q2) return;
    int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], xb, false);
    int ip = find_bin(phi_edges, phi, true);
    if (it < 0 || ix < 0 || ip < 0) return;
    cutflow.n_data_inrange++;

    // src/analysis/combine_analysis_branches.py supplies the run's prescale,
    // livetime, efficiency and inverse charge (mC) in scale. Multiplication by
    // Qrun_uC/Qtotal_uC replaces the run denominator with the total exposure:
    // sum_r [PS_r/(LT_r*eff_r)] * signal_counts_r / Qtotal_mC.
    // The charge fraction is dimensionless; do not insert another 1000 here.
    double charge_fraction = 1.0;
    if (std::isfinite(charge_uC) && charge_uC > 0.0 && std::isfinite(total_charge_uC) && total_charge_uC > 0.0) {
        charge_fraction = charge_uC / total_charge_uC;
    }
    // pi0_weight is the producer's mass-based signal estimate, not a binary
    // exclusivity gate. An upper missing-mass cut can still leave non-pi0
    // contributions; compare their tails in the companion notebook.
    double w = pi0_weight * static_cast<double>(scale) * charge_fraction;

    PhiBin& pb = slice(it, iq, ix).phi[ip];
    pb.data += w;
    // Conditional weighted-count variance only. pi0_weight was itself inferred
    // from mass/background fits upstream; their shared uncertainty is not
    // contained in sum(w^2), nor in the extraction's finite-MC covariance.
    pb.data_sumw2 += w * w;
    pb.n_data += 1;
    pb.weights.add(w);
    run_weight_moments[run_number].add(w);
    pb.mean_q2_data += q2 * w;
    pb.mean_xb_data += xb * w;
    pb.mean_tprime_data += tprime * w;
    if (use_helicity && helicity > 0) {
        pb.data_plus += w;
        pb.data_plus_sumw2 += w * w;
    } else if (use_helicity && helicity < 0) {
        pb.data_minus += w;
        pb.data_minus_sumw2 += w * w;
    }

    SliceResult& s = slice(it, iq, ix);
    s.sumw_data += w;
    s.sumw2_data += w * w;
    s.mean_q2_data += q2 * w;
    s.mean_xb_data += xb * w;
    s.mean_tprime_data += tprime * w;
    accumulate_global_histograms(q2, xb, tprime, wrap_phi(phi), w, 0.0);
}

inline void ExclPi0XSecAnalysis::fill_sim_event(float q2,
                                         float t,
                                         float tmin,
                                         float xb,
                                         float phi,
                                         float full_weight,
                                         int is_exclusive,
                                         float mmiss,
                                         float mpi0,
                                         float model_xsec,
                                         float vertex_q2,
                                         float vertex_W,
                                         float vertex_t,
                                         float vertex_phi,
                                         float vertex_hsxptari,
                                         float vertex_hsyptari,
                                         int helicity,
                                         bool use_helicity,
                                         float W) {
    (void)W;
    (void)helicity;
    (void)use_helicity; // This fit uses unpolarized totals only.
    cutflow.n_sim_total++;
    if (!selected_sim(is_exclusive, mmiss, mpi0)) return;
    double tprime = calc_tprime(t, tmin);
    if (cfg.prepare_forward_inputs) {
        if (!xsec_inside_diamond(cfg, xb, q2)) return;
        if (!std::isfinite(model_xsec) || std::fabs(model_xsec) < 1e-20f)
            die("Forward event cache: absent generator sigcm support");
        const double base = static_cast<double>(full_weight) / model_xsec;
        if (!(has_vertex_kinematics && vertex_q2 > 0 && vertex_W > cfg.mp && vertex_t > 0))
            die("Forward event cache: invalid matched vertex kinematics");
        const double tx = vertex_q2 / (static_cast<double>(vertex_W)*vertex_W - cfg.mp*cfg.mp + vertex_q2);
        const double tp = -static_cast<double>(vertex_t) - nps_xsec::forward_t(vertex_q2, vertex_W, cfg.mp, cfg.mpi0);
        const double eps = nps_xsec::vertex_epsilon_from_exclusive_simc(vertex_q2, vertex_W,
            vertex_hsxptari, vertex_hsyptari, cfg.hms_theta_deg, cfg.mp);
        for (double v : {static_cast<double>(q2), static_cast<double>(xb), tprime,
                         static_cast<double>(phi), static_cast<double>(vertex_q2), tx, tp,
                         static_cast<double>(vertex_phi), eps, base})
            if (!std::isfinite(v)) die("Forward event cache: nonfinite MC coordinates/weight");
        if (base < 0 || !(tx > 0 && tx < 1) || eps < 0 || eps > 1)
            die("Forward event cache: invalid MC weight/xB/epsilon");
        forward_mc_stream << forward_event_id << ',' << q2 << ',' << xb << ',' << tprime << ','
                          << wrap_phi(phi) << ',' << vertex_q2 << ',' << tx << ',' << tp << ','
                          << wrap_phi(vertex_phi) << ',' << eps << ',' << base << '\n';
        ++forward_mc_count;
        return;
    }
    if (!slice_passes_kin(q2, xb, tprime)) return;
    if (!std::isfinite(phi)) return;
    if (!std::isfinite(full_weight) || !std::isfinite(model_xsec)) die("Nonfinite selected SIMC weight/sigcm");
    if (std::fabs(model_xsec) < 1e-20f) {
        // A model zero is a hole in generator support: 0/0 cannot create a
        // unit-cross-section acceptance event. Report its frequency explicitly.
        ++cutflow.n_sim_zero_sigcm;
        return;
    }

    cutflow.n_sim_pass++;

    int it = find_bin(tprime_edges, tprime, false);
    int iq = find_bin(q2_edges, q2, false);
    if (iq < 0 || iq >= cfg.n_q2) return;
    int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], xb, false);
    int ip = find_bin(phi_edges, phi, true);
    if (it < 0 || ix < 0 || ip < 0) return;
    cutflow.n_sim_inrange++;

    // Divide out only the vertex model. The remaining integration weight keeps
    // all flux/Jacobian/radiation factors needed for a unit hadronic response.
    // It is an integrated yield-per-cross-section factor, not an acceptance
    // probability. Multiplying by the basis and fitted coefficients restores
    // yield/mC; another Q2, xB, t' or phi bin-width divisor would change units.
    const double base_w = static_cast<double>(full_weight) / static_cast<double>(model_xsec);
    if (!(std::isfinite(base_w) && base_w >= 0)) die("Invalid or negative demodeled SIMC weight");

    const double phiw = wrap_phi(phi);
    // SIMC sigcm is evaluated at the interaction vertex. Using event-level
    // epsilon and phi inside the response integral avoids replacing
    // sum(w*epsilon*cos(2phi)) by epsilon(mean)*sum(w*cos(2phi)).
    const bool valid_vertex = has_vertex_kinematics && std::isfinite(vertex_q2) &&
                              std::isfinite(vertex_W) && std::isfinite(vertex_t) &&
                              std::isfinite(vertex_phi) &&
                              vertex_q2 > 0.0f && vertex_W > cfg.mp && vertex_t > 0.0f;
    if (!valid_vertex) die("Selected SIMC event lacks valid vertex kinematics; migration requires truth");
    double physics_q2 = vertex_q2;
    double physics_W = vertex_W;
    double physics_phi = vertex_phi;
    double physics_xb = 0;
    const double xb_den = physics_W * physics_W - cfg.mp * cfg.mp + physics_q2;
    if (xb_den > 0.0) physics_xb = physics_q2 / xb_den;
    if (!(physics_xb>0 && physics_xb<1)) die("Invalid generated xB in selected SIMC event");
    const double eps_event = nps_xsec::vertex_epsilon_from_exclusive_simc(
        physics_q2, physics_W, vertex_hsxptari, vertex_hsyptari,
        cfg.hms_theta_deg, cfg.mp);
    const double physics_t=-static_cast<double>(vertex_t);
    const double physics_tprime=physics_t-nps_xsec::forward_t(physics_q2,physics_W,cfg.mp,cfg.mpi0);
    const int origin=nps_xsec::truth_block(physics_q2,physics_xb,physics_tprime,
                                         q2_edges,xb_edges_by_q2,tprime_edges);

    const double k_lt_event = std::sqrt(std::max(0.0, 2.0 * eps_event * (1.0 + eps_event)));
    // SIMC phipqi uses y=q x k and x=y x q. The data plane normals k x k'
    // and q x p_pi give the same cos(phi), cos(2phi) convention. Adding pi
    // would reverse the fitted LT sign while leaving U and TT unchanged.
    const double physics_phiw = wrap_phi(physics_phi);
    const int truth_ip = find_bin(phi_edges, physics_phiw, true);
    if (truth_ip < 0) die("Selected SIMC event has invalid vertex phi bin");

    PhiBin& pb = slice(it, iq, ix).phi[ip];
    pb.sim_base += base_w;
    // Use hard-vertex epsilon with generated phipqi for each event.
    const std::array<double, 3> event_basis = {
        base_w / (2.0 * TMath::Pi()),
        base_w * k_lt_event * std::cos(physics_phiw) / (2.0 * TMath::Pi()),
        base_w * eps_event * std::cos(2.0 * physics_phiw) / (2.0 * TMath::Pi())
    };
    // Off-diagonal response: the row is reconstruction, the block is origin.
    // Reconstructed phi determines the row; generated phi and vertex epsilon
    // determine the Fourier integrand. A phi migration needs no fitted phi bins.
    const int row=slice_index(it,iq,ix)*cfg.n_phi+ip;
    migration_response[row][origin].add(event_basis,q2,xb,tprime,base_w*eps_event);
    truth_phi_response[row][origin*cfg.n_phi+truth_ip].add(event_basis,q2,xb,tprime,
                                                           base_w*eps_event);
    truth_moments[origin].add(base_w,physics_q2,physics_xb,physics_t,physics_tprime,eps_event);
    truth_phi_moments[origin*cfg.n_phi+truth_ip].add(base_w,physics_q2,physics_xb,physics_t,physics_tprime,eps_event);
    if (cfg.diagnostics) {
        // Count once per response event, with no model/cross-section weighting.
        // The denominator excludes detector-lost/rejected events; these maps
        // show conditional coverage of the sample that actually enters A.
        h_migration_vertex_q2_xb->Fill(physics_xb, physics_q2);
        h_migration_reco_q2_xb->Fill(xb, q2);
    }

    pb.n_sim += 1;
    pb.mean_q2_sim += q2 * base_w;
    pb.mean_xb_sim += xb * base_w;
    pb.mean_tprime_sim += tprime * base_w;

    SliceResult& s = slice(it, iq, ix);
    s.sumw_sim += base_w;
    s.sumw2_sim += base_w * base_w;
    s.mean_q2_sim += q2 * base_w;
    s.mean_xb_sim += xb * base_w;
    s.mean_tprime_sim += tprime * base_w;
    // Generated means are accumulated by origin above, not by this row.
    accumulate_global_histograms(q2, xb, tprime, phiw, 0.0, base_w);
}

inline void ExclPi0XSecAnalysis::fill_from_trees() {
    // Data loop
    double dq2 = 0, dt = 0, dtmin = 0, dxb = 0, dphi = 0, dpi0_weight = 0, dW = 0;
    float dcharge_uC = 0;
    float dscale = 0;
    int dhelicity = 0;
    double dmmiss_all = 0, dmpi0_all = 0, dmmiss_all_corr = 0;
    int dexclusive_flag = 0;
    int drun_number = 0;
    bind_branch(t_data,"Q2", &dq2);
    bind_branch(t_data,"t", &dt);
    bind_branch(t_data,"tmin", &dtmin);
    bind_branch(t_data,"xB", &dxb);
    bind_branch(t_data,"phi", &dphi);
    bind_branch(t_data,"pi0_weight", &dpi0_weight);
    bind_branch(t_data,"scale", &dscale);
    bind_branch(t_data,"mmiss_all", &dmmiss_all);
    if (cfg.mmiss_select != "window") {
        bind_branch(t_data,"mpi0_all", &dmpi0_all);
        bind_branch(t_data,"mmiss_all_corr", &dmmiss_all_corr);
        const char* flag = cfg.mmiss_select == "mcd" ?
            "is_exclusive_mcd_combined" : "is_exclusive_ellipse_combined";
        bind_branch(t_data,flag, &dexclusive_flag);
    }
    bind_branch(t_data,"charge_uC", &dcharge_uC);
    bind_branch(t_data,"run_number", &drun_number);
    bind_branch(t_data,"W", &dW);
    if (has_helicity) bind_branch(t_data,"helicity", &dhelicity);

    // Sum each run's exposure exactly once. Inconsistent or invalid exposure
    // must stop an absolute extraction rather than use a neutral fallback.
    const Long64_t ndata = t_data->GetEntries();
    std::map<int,std::pair<float,float>> exposure_by_run;
    double total_charge_uC = 0.0;

    for (Long64_t i = 0; i < ndata; ++i) {
        t_data->GetEntry(i);
        if(!(std::isfinite(dcharge_uC) && dcharge_uC>0 && std::isfinite(dscale) && dscale>0))
            die("Invalid run charge/scale in data entry "+std::to_string(i));
        const auto inserted=exposure_by_run.emplace(drun_number,std::make_pair(dcharge_uC,dscale));
        if (inserted.second) {
            total_charge_uC += dcharge_uC;
        } else {
            const auto reference=inserted.first->second;
            if(std::abs(dcharge_uC-reference.first)>1e-6*reference.first ||
               std::abs(dscale-reference.second)>1e-6*reference.second)
                die("Inconsistent charge/scale within run "+std::to_string(drun_number));
        }
    }

    if (!(total_charge_uC > 0.0)) {
        die("Combined data have no valid total run charge");
    }

    // Second pass: fill events using total_charge_uC
    for (Long64_t i = 0; i < ndata; ++i) {
        t_data->GetEntry(i);
        forward_event_id = static_cast<unsigned long long>(i);
        if ((cfg.mmiss_select == "mcd" || cfg.mmiss_select == "ellipse") &&
            std::isfinite(dpi0_weight * dscale) && dpi0_weight * dscale > 0 &&
            mass_geometry.contains(dmpi0_all, dmmiss_all) != (dexclusive_flag != 0))
            die("Combined " + cfg.mmiss_select + " geometry disagrees with stored data flag at entry " +
                std::to_string(i) + "; use metadata from the matching data file.");
        fill_mmiss_data_diagnostic(dq2, dt, dtmin, dxb, dphi, dmmiss_all,
                                   dpi0_weight, dscale, dcharge_uC, total_charge_uC);
        if (h_mmiss_reco_data_all && std::isfinite(dmmiss_all)) {
            h_mmiss_reco_data_all->Fill(dmmiss_all);
            if (selected_data(dmmiss_all, dmpi0_all, dexclusive_flag))
                h_mmiss_reco_data_selected->Fill(dmmiss_all);
        }
        if (h_mmiss_corr_data_all && std::isfinite(dmmiss_all_corr)) {
            h_mmiss_corr_data_all->Fill(dmmiss_all_corr);
            if (selected_data(dmmiss_all, dmpi0_all, dexclusive_flag))
                h_mmiss_corr_data_selected->Fill(dmmiss_all_corr);
        }
        if (h_mass_data_all && std::isfinite(dmpi0_all) && std::isfinite(dmmiss_all)) {
            h_mass_data_all->Fill(dmpi0_all, dmmiss_all);
            if (selected_data(dmmiss_all, dmpi0_all, dexclusive_flag))
                h_mass_data_selected->Fill(dmpi0_all, dmmiss_all);
        }
        fill_data_event(dq2, dt, dtmin, dxb, dphi, dpi0_weight, dscale, dcharge_uC, total_charge_uC, dmmiss_all, dmpi0_all, dexclusive_flag, dhelicity, has_helicity, dW, drun_number);
    }

    // SIMC loop
    float sim_q2 = 0, sim_t = 0, sim_tmin = 0, sim_xb = 0, sim_phi = 0, full_weight = 0, sim_model_xsec = 0, sim_W = 0;
    int sim_is_exclusive = 0, sim_helicity = 0;
    float sim_mmiss = 0, sim_mpi0 = 0, sim_vertex_q2 = 0, sim_vertex_W = 0, sim_vertex_t = 0, sim_vertex_phi = 0;
    ULong64_t sim_event_id = 0;
    Float_t raw_q2i = 0, raw_Wi = 0, raw_ti = 0, raw_phipqi = 0, raw_sigcm = 0;
    Float_t raw_hsxptari = 0, raw_hsyptari = 0;
    bind_branch(t_sim,"Q2", &sim_q2);
    bind_branch(t_sim,"t", &sim_t);
    bind_branch(t_sim,"tmin", &sim_tmin);
    bind_branch(t_sim,"xB", &sim_xb);
    bind_branch(t_sim,"phi", &sim_phi);
    bind_branch(t_sim,"full_weight", &full_weight);
    bind_branch(t_sim,"is_exclusive", &sim_is_exclusive);
    bind_branch(t_sim,"mmiss", &sim_mmiss);
    if (cfg.mmiss_select != "window") bind_branch(t_sim,"mpi0", &sim_mpi0);
    bind_branch(t_sim,"event_id", &sim_event_id);
    t_vertex->SetBranchStatus("*", 0);
    for (const char* branch : {"Q2i", "Wi", "ti", "phipqi", "sigcm", "hsxptari", "hsyptari"})
        t_vertex->SetBranchStatus(branch, 1);
    bind_branch(t_vertex,"Q2i", &raw_q2i);
    bind_branch(t_vertex,"Wi", &raw_Wi);
    bind_branch(t_vertex,"ti", &raw_ti);
    bind_branch(t_vertex,"phipqi", &raw_phipqi);
    bind_branch(t_vertex,"sigcm", &raw_sigcm);
    bind_branch(t_vertex,"hsxptari", &raw_hsxptari);
    bind_branch(t_vertex,"hsyptari", &raw_hsyptari);
    bind_branch(t_sim,"W", &sim_W);
    bind_branch(t_sim,model_xsec_branch.c_str(), &sim_model_xsec);
    if (has_sim_helicity) bind_branch(t_sim,"helicity", &sim_helicity);

    const Long64_t nsim = t_sim->GetEntries();
    for (Long64_t i = 0; i < nsim; ++i) {
        t_sim->GetEntry(i);
        forward_event_id = sim_event_id;
        fill_mmiss_sim_diagnostic(sim_q2, sim_t, sim_tmin, sim_xb, sim_phi,
                                  sim_mmiss, full_weight, sim_is_exclusive);
        if (h_mmiss_sim_all && sim_is_exclusive && std::isfinite(sim_mmiss)) {
            h_mmiss_sim_all->Fill(sim_mmiss);
            if (selected_sim(sim_is_exclusive, sim_mmiss, sim_mpi0))
                h_mmiss_sim_selected->Fill(sim_mmiss);
        }
        if (h_mass_sim_all && sim_is_exclusive && std::isfinite(sim_mpi0) && std::isfinite(sim_mmiss)) {
            h_mass_sim_all->Fill(sim_mpi0, sim_mmiss);
            if (selected_sim(sim_is_exclusive, sim_mmiss, sim_mpi0))
                h_mass_sim_selected->Fill(sim_mpi0, sim_mmiss);
        }
        if (selected_sim(sim_is_exclusive, sim_mmiss, sim_mpi0)) {
            if (sim_event_id >= static_cast<ULong64_t>(t_vertex->GetEntries()) ||
                t_vertex->GetEntry(static_cast<Long64_t>(sim_event_id)) <= 0)
                die("Raw h10 entry unavailable for smeared event_id " + std::to_string(sim_event_id));
            // Both values are copied from the same generated event in this
            // production. An event-wise match guards against a different
            // SIMC file or altered GEANT entry ordering.
            const double tolerance = 1e-6 * std::max(std::fabs(static_cast<double>(raw_sigcm)),
                                                     std::fabs(static_cast<double>(sim_model_xsec))) + 1e-20;
            if (!std::isfinite(raw_sigcm) || !std::isfinite(sim_model_xsec) ||
                std::fabs(static_cast<double>(raw_sigcm) - static_cast<double>(sim_model_xsec)) > tolerance)
                die("Raw SIMC sigcm mismatch for smeared event_id " + std::to_string(sim_event_id) +
                    "; wrong SIMC production or event ordering.");
            sim_vertex_q2 = raw_q2i; sim_vertex_W = raw_Wi;
            sim_vertex_t = raw_ti; sim_vertex_phi = raw_phipqi;
            ++n_vertex_matched;
        }
        fill_sim_event(sim_q2,
                       sim_t,
                       sim_tmin,
                       sim_xb,
                       sim_phi,
                       full_weight,
                       sim_is_exclusive,
                       sim_mmiss,
                       sim_mpi0,
                       sim_model_xsec,
                       sim_vertex_q2,
                       sim_vertex_W,
                       sim_vertex_t,
                       sim_vertex_phi,
                       raw_hsxptari,
                       raw_hsyptari,
                       sim_helicity,
                       has_sim_helicity,
                       sim_W);
    }
    if (cutflow.n_sim_zero_sigcm > 0)
        die("Selected exclusive SIMC events with near-zero sigcm: " +
             std::to_string(cutflow.n_sim_zero_sigcm) +
             ". Their response is absent; inspect generator-model support before interpreting an absolute cross section.");
    if (has_vertex_kinematics && cutflow.n_sim_vertex_fallback > 0)
        warn("Selected SIMC events with invalid vertex values, using reconstructed harmonics: " +
             std::to_string(cutflow.n_sim_vertex_fallback));
    if (vertex_from_raw_simc)
        std::cout << "Raw SIMC vertex sigcm matches: " << n_vertex_matched
                  << " selected exclusive smeared events" << std::endl;

    // The integrated response keeps its physical SIMC normalization. Never
    // scale it to the data yield. Target correction acts on data in the fit.

    // finalize all weighted means and per-phi errors
    for (auto& s : slices) {
        if (s.sumw_data > 0.0) {
            s.mean_q2_data /= s.sumw_data;
            s.mean_xb_data /= s.sumw_data;
            s.mean_tprime_data /= s.sumw_data;
        }
        if (s.sumw_sim > 0.0) {
            s.mean_q2_sim /= s.sumw_sim;
            s.mean_xb_sim /= s.sumw_sim;
            s.mean_tprime_sim /= s.sumw_sim;
        }
        for (auto& pb : s.phi) {
            // Signed signal weights can have a negative total. A weighted
            // reconstructed mean is defined whenever its sum is nonzero.
            if (pb.data != 0.0) {
                pb.mean_q2_data /= pb.data;
                pb.mean_xb_data /= pb.data;
                pb.mean_tprime_data /= pb.data;
            }
            if (pb.sim_base > 0.0) {
                pb.mean_q2_sim /= pb.sim_base;
                pb.mean_xb_sim /= pb.sim_base;
                pb.mean_tprime_sim /= pb.sim_base;
            }
        }
    }
}

inline void ExclPi0XSecAnalysis::compute_mmiss_area_scales() {
    // Same area-matching operation used by nps_sim_smearing_new.C previews:
    //   scale = Integral(data) / Integral(simulation)
    // Here simulation means the smeared generated-exclusive spectrum, and the
    // integral uses exactly the configured 0.6-1.1 GeV extraction window.
    for (auto& m : mmiss_slices) {
        double data_integral = 0.0, sim_integral = 0.0;
        for (int bin = 1; bin <= m.data->GetNbinsX(); ++bin) {
            const double center = m.data->GetBinCenter(bin);
            if (center >= cfg.mmiss_lower_gev && center <= cfg.mmiss_upper_gev) {
                data_integral += m.data->GetBinContent(bin);
                sim_integral += m.exclusive->GetBinContent(bin);
            }
        }
        if (!(data_integral > 0.0 && sim_integral > 0.0)) continue;
        m.exclusive_normalization_scale = data_integral / sim_integral;
    }
}
