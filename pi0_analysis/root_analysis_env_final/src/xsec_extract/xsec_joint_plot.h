#pragma once

// Import one setting's marginal joint-fit coefficients/covariance. Rebuild
// the event diagnostics, but never call the single-setting solver or toys.
#include "xsec_analysis.h"

inline void ExclPi0XSecAnalysis::load_joint_plot_input() {
    std::ifstream input(cfg.joint_plot_input);
    if (!input) die("Cannot read joint plotting input: " + cfg.joint_plot_input);
    std::string magic;
    size_t blocks=0, nr=0;
    int boundary=0;
    input >> magic >> blocks >> nr >> boundary;
    if (magic != "joint_plot_v1" || !blocks || nr != migration_response.size())
        die("Invalid joint plotting input dimensions/version");
    invalidate_fit_results();
    const size_t np=3*blocks;
    migration_fit.parameters.resize(np);
    migration_fit.covariance.resize(np*np);
    input >> migration_fit.chi2 >> migration_fit.ndf >> migration_fit.rank
          >> migration_fit.condition >> mc_iterations;
    active_truth_blocks.resize(blocks);
    positivity_feasibility_tolerances.resize(blocks);
    for (size_t i=0;i<blocks;++i) {
        input >> active_truth_blocks[i] >> positivity_feasibility_tolerances[i];
        if (active_truth_blocks[i]<0 || active_truth_blocks[i]>=static_cast<int>(truth_moments.size()) ||
            (i && active_truth_blocks[i]<=active_truth_blocks[i-1])) die("Invalid joint plotting block map");
        for (int a=0;a<3;++a) input >> migration_fit.parameters[3*i+a];
    }
    // std::istream's floating extraction is not portable for CSV NaN tokens.
    const auto number=[&]() { std::string token; input>>token; return std::stod(token); };
    for (double& value:migration_fit.covariance) value=number();
    response_design.assign(nr,std::vector<double>(np));
    fit_rows.clear(); fit_variance.clear();
    const auto agree=[&](double actual,double expected,const char* kind,size_t row) {
        if (!std::isfinite(actual) || !std::isfinite(expected) ||
            std::abs(actual-expected)>1e-8*std::max({std::abs(actual),std::abs(expected),1e-100}))
            die(std::string("Joint plotting ")+kind+" differs from saved fit in row "+std::to_string(row));
    };
    std::vector<double> predictions(nr);
    for (size_t r=0;r<nr;++r) {
        int included=0;
        input >> included;
        const double data=number(),variance=number(),used=number();
        predictions[r]=number();
        agree(slices[r/cfg.n_phi].phi[r%cfg.n_phi].data,data,"data",r);
        agree(slices[r/cfg.n_phi].phi[r%cfg.n_phi].data_sumw2,variance,"data variance",r);
        for (size_t i=0;i<blocks;++i) for (int a=0;a<3;++a) {
            const double expected=number();
            const double actual=migration_response[r][active_truth_blocks[i]].basis[a];
            agree(actual,expected,"response",r);
            response_design[r][3*i+a]=actual;
        }
        if (included) { fit_rows.push_back(static_cast<int>(r)); fit_variance.push_back(used); }
    }
    if (!input) die("Truncated joint plotting input");
    positivity_boundary_active=boundary!=0;
    mc_converged=true;
    retained_fit_groups.assign(cfg.n_q2*cfg.n_xb,true);
    successful_fit_groups=cfg.n_q2*cfg.n_xb;
    for (size_t b=0;b<slices.size();++b)
        if (std::find(active_truth_blocks.begin(),active_truth_blocks.end(),b)==active_truth_blocks.end())
            die("Joint plotting input is missing a published truth block");
    finalize_fit_subset(retained_fit_groups);
    for (size_t r=0;r<nr;++r)
        agree(slices[r/cfg.n_phi].phi[r%cfg.n_phi].sim,predictions[r],"prediction",r);
    for (auto& slice:slices) slice.fit_scope="joint_shared_LT_TT_independent_U";
    log("Imported joint-fit coefficients and covariance for plotting; no individual fit.");
}
