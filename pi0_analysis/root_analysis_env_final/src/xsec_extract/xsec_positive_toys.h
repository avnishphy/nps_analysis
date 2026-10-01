#pragma once

// Conditional parametric refits for plots when the physical positivity
// boundary makes the ordinary inverse-information covariance inapplicable.
// Draw the observed weighted yields at the fitted forward mean using the
// converged row variances, and repeat the complete constrained variance fit.
// Response moments, binning, target factor and variance estimates are held
// fixed. These are sampling standard deviations, not coverage intervals.
#include "xsec_experimental_points.h"

inline void ExclPi0XSecAnalysis::compute_positive_toy_errors() {
    positive_toy_parameter_covariance.clear();
    positive_toy_point_covariance.clear();
    positive_toys_successful = 0;
    if (cfg.fit_objective=="scaled-poisson" || !positivity_boundary_active || migration_fit.parameters.empty()) return;

    constexpr int requested = 256;
    constexpr unsigned long long seed = 20260924ULL;
    const size_t np = migration_fit.parameters.size(), nr = experimental_points.size();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> design, parameter_samples, point_samples   ;
    std::vector<double> data_variance, generating_mean;
    std::vector<int> active_index(truth_moments.size(), -1);
    for (size_t i=0;i<active_truth_blocks.size();++i)
        active_index[active_truth_blocks[i]]=static_cast<int>(i);
    std::vector<double> epsilon_max;
    for (int b:active_truth_blocks) epsilon_max.push_back(truth_moments[b].epsilon_max);
    for (int r:fit_rows) {
        design.push_back(response_design[r]);
        data_variance.push_back(slices[r/cfg.n_phi].phi[r%cfg.n_phi].data_sumw2);
        generating_mean.push_back(std::inner_product(response_design[r].begin(),
            response_design[r].end(),migration_fit.parameters.begin(),0.));
    }
    std::mt19937_64 engine(seed);
    std::normal_distribution<double> normal(0.,1.);
    for (int toy=0;toy<requested;++toy) {
        std::vector<double> y(fit_rows.size());
        for (size_t i=0;i<y.size();++i)
            y[i]=generating_mean[i]+std::sqrt(fit_variance[i])*normal(engine);
        try {
            auto solve=[&](const std::vector<double>& variance) {
                return nps_xsec::solve_positive_response(design,y,variance,epsilon_max,
                                                           cfg.rank_tolerance).fit;
            };
            std::vector<double> variance=data_variance;
            auto fitted=solve(variance);
            bool converged=cfg.fit_variance_mode=="data";
            for (int iteration=0;!converged && iteration<cfg.mc_max_iterations;++iteration) {
                std::vector<double> next_variance=data_variance;
                for (size_t i=0;i<fit_rows.size();++i)
                    next_variance[i]+=nps_xsec::mc_prediction_variance(
                        migration_response[fit_rows[i]],active_truth_blocks,fitted.parameters);
                auto next=solve(next_variance);
                double parameter_change=0.,variance_change=0.;
                for (size_t j=0;j<np;++j) {
                    const double scale=std::max(std::abs(next.parameters[j]),
                        std::sqrt(next.covariance[j*np+j]));
                    parameter_change=std::max(parameter_change,
                        std::abs(next.parameters[j]-fitted.parameters[j])/scale);
                }
                for (size_t i=0;i<variance.size();++i)
                    variance_change=std::max(variance_change,
                        std::abs(next_variance[i]-variance[i])/next_variance[i]);
                fitted=std::move(next);
                variance=std::move(next_variance);
                converged=std::max(parameter_change,variance_change)<cfg.mc_fit_tolerance;
            }
            if (!converged) continue;
            std::vector<double> points(nr,nan);
            for (size_t i=0;i<fit_rows.size();++i) {
                const int r=fit_rows[i], b=r/cfg.n_phi, a=active_index[b];
                if (experimental_points[r].status!="central_only_boundary" || a<0) continue;
                const auto found=truth_phi_response[r].find(r);
                if (found==truth_phi_response[r].end()) continue;
                const auto basis=sigma_model_basis(.5*(phi_edges[r%cfg.n_phi]+
                    phi_edges[r%cfg.n_phi+1]),experimental_points[r].epsilon_reference,0.,0);
                double m=0.,f=0.,absolute=0.;
                for (int term=0;term<3;++term) {
                    m+=found->second.basis[term]*fitted.parameters[3*a+term];
                    f+=basis[term]*fitted.parameters[3*a+term];
                }
                const auto& row=response_design[r];
                const double mu=std::inner_product(row.begin(),row.end(),fitted.parameters.begin(),0.);
                for (size_t j=0;j<np;++j) absolute+=std::abs(row[j]*fitted.parameters[j]);
                if (!std::isfinite(m) || std::abs(m)<=1e-12*absolute || m==0.) continue;
                points[r]=(y[i]-mu+m)*f/m;
            }
            parameter_samples.push_back(std::move(fitted.parameters));
            point_samples.push_back(std::move(points));
        } catch (const std::runtime_error&) {
            // Singular or nonconvergent toys are counted below, never replaced
            // with the central estimate or a zero uncertainty.
        }
    }
    positive_toys_successful=static_cast<int>(parameter_samples.size());
    if (positive_toys_successful<requested*9/10) {
        warn("Positivity refit toys: too many failures ("+std::to_string(positive_toys_successful)+
             "/"+std::to_string(requested)+"); toy spread unavailable");
        positive_toys_successful=0;
        return;
    }
    std::vector<double> mean(np,0.);
    for (const auto& sample:parameter_samples)
        for (size_t i=0;i<np;++i) mean[i]+=sample[i]/positive_toys_successful;
    positive_toy_parameter_covariance.assign(np*np,0.);
    for (const auto& sample:parameter_samples)
        for (size_t i=0;i<np;++i) for (size_t j=0;j<np;++j)
            positive_toy_parameter_covariance[i*np+j]+=
                (sample[i]-mean[i])*(sample[j]-mean[j])/(positive_toys_successful-1);
    positive_toy_point_covariance.assign(nr*nr,nan);
    for (size_t r=0;r<nr;++r) for (size_t s=0;s<=r;++s) {
        double mr=0.,ms=0.,cross=0.; int count=0;
        for (const auto& sample:point_samples) if (std::isfinite(sample[r]) && std::isfinite(sample[s])) {
            mr+=sample[r];ms+=sample[s];cross+=sample[r]*sample[s];++count;
        }
        if (count<positive_toys_successful*9/10 || count<2) continue;
        const double covariance=(cross-mr*ms/count)/(count-1);
        positive_toy_point_covariance[r*nr+s]=covariance;
        positive_toy_point_covariance[s*nr+r]=covariance;
    }
    std::ofstream out(fs::path(cfg.out_dir)/"positivity_refit_toy_covariance.csv");
    out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"kind,index_i,index_j,covariance,unit,successful_toys,requested_toys,seed\n";
    for(size_t i=0;i<np;++i) for(size_t j=0;j<np;++j)
        out<<"parameter,"<<i<<','<<j<<','<<positive_toy_parameter_covariance[i*np+j]
           <<",(ub/MeV2)^2,"<<positive_toys_successful<<','<<requested<<','<<seed<<'\n';
    for(size_t i=0;i<nr;++i) for(size_t j=0;j<nr;++j)
        if(std::isfinite(positive_toy_point_covariance[i*nr+j]))
            out<<"experimental_point,"<<i<<','<<j<<','<<positive_toy_point_covariance[i*nr+j]
               <<",(ub/MeV2/rad)^2,"<<positive_toys_successful<<','<<requested<<','<<seed<<'\n';
    log("Positivity conditional refit toys: "+std::to_string(positive_toys_successful)+
        "/"+std::to_string(requested)+" converged; plotting sampling standard deviations");
}

inline double ExclPi0XSecAnalysis::plot_parameter_variance(const std::vector<double>& gradient) const {
    const auto& covariance=positivity_boundary_active ? positive_toy_parameter_covariance : migration_fit.covariance;
    const size_t np=migration_fit.parameters.size();
    if (covariance.size()!=np*np || gradient.size()!=np)
        return std::numeric_limits<double>::quiet_NaN();
    double value=0.;
    for(size_t i=0;i<np;++i) for(size_t j=0;j<np;++j)
        value+=gradient[i]*covariance[i*np+j]*gradient[j];
    return std::isfinite(value) && value>=-1e-10*std::abs(value) ? std::max(0.,value) :
        std::numeric_limits<double>::quiet_NaN();
}

inline double ExclPi0XSecAnalysis::plot_coefficient_error(int truth_block,int term) const {
    const auto found=std::find(active_truth_blocks.begin(),active_truth_blocks.end(),truth_block);
    if(found==active_truth_blocks.end() || term<0 || term>=3)
        return std::numeric_limits<double>::quiet_NaN();
    std::vector<double> gradient(migration_fit.parameters.size(),0.);
    gradient[3*(found-active_truth_blocks.begin())+term]=1.;
    return std::sqrt(plot_parameter_variance(gradient));
}

inline double ExclPi0XSecAnalysis::plot_point_error(int row) const {
    if(!positivity_boundary_active) return experimental_points[row].sigma_exp_err;
    const size_t nr=experimental_points.size();
    if(positive_toy_point_covariance.size()!=nr*nr)
        return std::numeric_limits<double>::quiet_NaN();
    const double variance=positive_toy_point_covariance[static_cast<size_t>(row)*nr+row];
    return std::isfinite(variance) && variance>=0. ? std::sqrt(variance) :
        std::numeric_limits<double>::quiet_NaN();
}
