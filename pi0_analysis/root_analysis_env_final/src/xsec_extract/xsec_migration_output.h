#pragma once

// Preserve the complete linear problem, including nuisance parameters, so an
// independent program can reconstruct the fit without rerunning event loops.
// CSV uses max_digits10 precision; ROOT matrices retain native doubles. All
// cross-section coefficients use ub/MeV2 and covariance uses their square.
#include "xsec_analysis.h"

inline void ExclPi0XSecAnalysis::write_migration_results() {
    fout->cd();
    const int nr=static_cast<int>(response_design.size());
    const int np=static_cast<int>(migration_fit.parameters.size());
    const int nb=static_cast<int>(truth_moments.size());
    const int published=static_cast<int>(slices.size());
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const double target_fraction=cfg.tgt_contam_err/cfg.tgt_contam;
    auto csv=[&](const char* name) {
        std::ofstream out(fs::path(cfg.out_dir)/name);
        if(!out) die(std::string("Cannot create migration CSV: ")+name);
        out<<std::setprecision(std::numeric_limits<double>::max_digits10);
        return out;
    };
    auto label=[&](int b) {
        return b<published?std::string("published"):std::string(nps_xsec::guard_name(b-published));
    };
    auto indices=[&](int b) {
        return b<published?std::array<int,3>{b/(cfg.n_q2*cfg.n_xb),(b/cfg.n_xb)%cfg.n_q2,b%cfg.n_xb}:
                           std::array<int,3>{-1,-1,-1};
    };

    auto status_csv = csv("fit_status.csv");
    status_csv << "iq,ix,q2_lo,q2_hi,xb_lo,xb_hi,fit_ok,fit_scope,failure_reason\n";
    TTree status_tree("fit_status", "Q2/xB groups retained in the final global fit");
    int status_iq = 0, status_ix = 0, status_ok = 0;
    std::string status_scope, status_reason;
    status_tree.Branch("iq", &status_iq); status_tree.Branch("ix", &status_ix);
    status_tree.Branch("fit_ok", &status_ok); status_tree.Branch("fit_scope", &status_scope);
    status_tree.Branch("failure_reason", &status_reason);
    for (status_iq = 0; status_iq < cfg.n_q2; ++status_iq)
        for (status_ix = 0; status_ix < cfg.n_xb; ++status_ix) {
            const auto& s = slice(0, status_iq, status_ix);
            status_ok = s.fit_xsec.ok; status_scope = s.fit_scope; status_reason = s.fit_failure_reason;
            status_tree.Fill();
            status_csv << status_iq << ',' << status_ix << ',' << q2_edges[status_iq] << ','
                       << q2_edges[status_iq + 1] << ',' << xb_edges_by_q2[status_iq][status_ix] << ','
                       << xb_edges_by_q2[status_iq][status_ix + 1] << ',' << status_ok << ','
                       << status_scope << ',' << xsec_csv_quote(status_reason) << '\n';
        }
    status_tree.Write();
    auto attempts_csv = csv("fit_attempts.csv");
    attempts_csv << "attempt,retained_q2_xb_groups,fit_ok,failure_reason\n";
    TTree attempts_tree("fit_attempts", "Global subset attempts in deterministic search order");
    int attempt_index = 0, attempt_ok = 0;
    std::string attempt_groups, attempt_reason;
    attempts_tree.Branch("attempt", &attempt_index); attempts_tree.Branch("fit_ok", &attempt_ok);
    attempts_tree.Branch("retained_q2_xb_groups", &attempt_groups);
    attempts_tree.Branch("failure_reason", &attempt_reason);
    for (const auto& a : fit_attempts) {
        ++attempt_index; attempt_ok = a.ok; attempt_reason = a.reason; attempt_groups.clear();
        for (size_t g = 0; g < a.groups.size(); ++g) if (a.groups[g]) {
            if (!attempt_groups.empty()) attempt_groups += ';';
            attempt_groups += std::to_string(g / cfg.n_xb) + ':' + std::to_string(g % cfg.n_xb);
        }
        attempts_tree.Fill();
        attempts_csv << attempt_index << ',' << xsec_csv_quote(attempt_groups) << ',' << attempt_ok
                     << ',' << xsec_csv_quote(attempt_reason) << '\n';
    }
    attempts_tree.Write();
    TParameter<int>("fit_retained_q2_xb_groups", successful_fit_groups).Write();
    TParameter<int>("fit_used_subset_recovery", fit_fallback ? 1 : 0).Write();
    TParameter<int>("fit_positive_xsec", cfg.positive_xsec ? 1 : 0).Write();
    TObjString((cfg.fit_objective=="scaled-poisson" ? "ignored_fixed_response" : cfg.fit_variance_mode.c_str())).Write("fit_variance_mode");
    TObjString(cfg.fit_objective.c_str()).Write("fit_objective");
    TParameter<int>("fit_positivity_boundary_active", positivity_boundary_active ? 1 : 0).Write();
    TParameter<int>("fit_positivity_iterations", positivity_iterations).Write();
    TParameter<int>("fit_positive_refit_toys_successful",positive_toys_successful).Write();
    TObjString(np == 0 ? "unavailable: fit failed" : positivity_boundary_active ?
        "unavailable: constrained boundary; use dedicated profile/toy intervals" :
        (cfg.fit_objective=="scaled-poisson" ? "conditional Minuit2 Hessian transformed to physical coefficients" :
        "conditional known-variance inverse information")).Write("migration_covariance_status");

    // This is the curvature of the UNCONSTRAINED quadratic objective at the
    // final GLS weights, not a covariance of an estimate on the boundary.
    auto curvature_csv = csv("migration_curvature_inverse.csv");
    curvature_csv << "parameter_i,parameter_j,"
                  << (cfg.fit_objective=="scaled-poisson" ? "physical_Hessian_covariance_diagnostic" : "unconstrained_curvature_inverse") << "\n";
    TMatrixD curvature(np, np);
    for (int i = 0; i < np; ++i) for (int j = 0; j < np; ++j) {
        curvature(i,j) = fit_curvature_inverse[static_cast<size_t>(i)*np+j];
        curvature_csv << i << ',' << j << ',' << curvature(i,j) << '\n';
    }
    if (np > 0) curvature.Write(cfg.fit_objective=="scaled-poisson" ?
        "migration_Hessian_covariance_diagnostic" : "migration_unconstrained_curvature_inverse");
    if (cfg.fit_objective=="scaled-poisson" && np>0) {
        auto hcorr_csv=csv("scaled_poisson_hessian_correlation_diagnostic.csv");
        hcorr_csv<<"parameter_i,parameter_j,Hessian_correlation_diagnostic_not_boundary_interval\n";
        TMatrixD hcorr(np,np);
        for (int i=0;i<np;++i) for (int j=0;j<np;++j) {
            const double vi=curvature(i,i),vj=curvature(j,j);
            hcorr(i,j)=vi>0. && vj>0. ? curvature(i,j)/std::sqrt(vi*vj) : nan;
            hcorr_csv<<i<<','<<j<<','<<hcorr(i,j)<<'\n';
        }
        hcorr.Write("migration_Hessian_correlation_diagnostic");
    }
    if (positive_toy_parameter_covariance.size()==static_cast<size_t>(np)*np) {
        TMatrixD toy_covariance(np,np),toy_correlation(np,np);
        for(int i=0;i<np;++i) for(int j=0;j<np;++j) {
            const double v=positive_toy_parameter_covariance[static_cast<size_t>(i)*np+j];
            const double vi=positive_toy_parameter_covariance[static_cast<size_t>(i)*np+i];
            const double vj=positive_toy_parameter_covariance[static_cast<size_t>(j)*np+j];
            toy_covariance(i,j)=v;
            toy_correlation(i,j)=vi>0. && vj>0. ? v/std::sqrt(vi*vj) :
                std::numeric_limits<double>::quiet_NaN();
        }
        toy_covariance.Write("migration_positive_refit_toy_covariance");
        toy_correlation.Write("migration_positive_refit_toy_correlation");
    }

    auto positivity_csv = csv("positivity_diagnostics.csv");
    positivity_csv << "active_block_index,truth_block,epsilon_max,minimum_response_bracket,at_cos_phi,positivity_enabled,boundary_active,feasibility_tolerance,boundary_tolerance\n";
    TTree positivity_tree("positivity_diagnostics", "Full-phi minimum bracket (2pi times angular cross section) at largest observed epsilon");
    int pos_active = 0, pos_block = 0, pos_enabled = cfg.positive_xsec, pos_boundary = 0;
    double pos_eps = 0, pos_min = 0, pos_cos = 0;
    double pos_feasibility_tolerance = 0, pos_boundary_tolerance = 0;
    positivity_tree.Branch("active_block_index", &pos_active);
    positivity_tree.Branch("truth_block", &pos_block);
    positivity_tree.Branch("epsilon_max", &pos_eps);
    positivity_tree.Branch("minimum_response_bracket", &pos_min);
    positivity_tree.Branch("at_cos_phi", &pos_cos);
    positivity_tree.Branch("positivity_enabled", &pos_enabled);
    positivity_tree.Branch("boundary_active", &pos_boundary);
    positivity_tree.Branch("feasibility_tolerance", &pos_feasibility_tolerance);
    positivity_tree.Branch("boundary_tolerance", &pos_boundary_tolerance);
    for (pos_active = 0; pos_active < np/3; ++pos_active) {
        pos_block = active_truth_blocks[pos_active];
        pos_eps = truth_moments[pos_block].epsilon_max;
        const double u = migration_fit.parameters[3*pos_active];
        const double lt = migration_fit.parameters[3*pos_active+1];
        const double tt = migration_fit.parameters[3*pos_active+2];
        const double a = 2*pos_eps*tt, b = std::sqrt(2*pos_eps*(1+pos_eps))*lt;
        pos_cos = b > 0 ? -1 : 1;
        if (a > 0 && std::abs(b/(2*a)) <= 1) pos_cos = -b/(2*a);
        pos_min = nps_xsec::minimum_response(u, lt, tt, pos_eps);
        pos_feasibility_tolerance = cfg.positive_xsec ? positivity_feasibility_tolerances[pos_active] : 0.0;
        pos_boundary_tolerance = cfg.positive_xsec ? positivity_boundary_tolerances[pos_active] : 0.0;
        pos_boundary = cfg.positive_xsec && pos_min <= pos_boundary_tolerance;
        positivity_tree.Fill();
        positivity_csv << pos_active << ',' << pos_block << ',' << pos_eps << ',' << pos_min << ','
                       << pos_cos << ',' << pos_enabled << ',' << pos_boundary << ','
                       << pos_feasibility_tolerance << ',' << pos_boundary_tolerance << '\n';
    }
    positivity_tree.Write();

    // Column index j is 3*active_block_index + component, with component
    // 0=U, 1=LT (legacy CSV calls it TL), 2=TT. No omitted guard is assigned
    // a zero-valued fitted parameter; inactive blocks are listed separately.
    TMatrixD design(nr,np),covariance(np,np),correlation(np,np),target_covariance(np,np);
    TVectorD parameters(np),singular(migration_fit.singular_values.size());
    auto response_csv=csv("migration_design.csv");
    response_csv<<"reco_row,parameter_index,response\n";
    for(int r=0;r<nr;++r) for(int j=0;j<np;++j) {
        design(r,j)=response_design[r][j];
        response_csv<<r<<','<<j<<','<<design(r,j)<<'\n';
    }
    auto covariance_csv=csv("migration_covariance.csv");
    covariance_csv<<"parameter_i,parameter_j,"
                  <<(cfg.fit_objective=="scaled-poisson" ? "conditional_Hessian_covariance" : "stat_plus_mc_covariance")
                  <<",target_correlated_covariance\n";
    auto correlation_csv=csv("migration_correlation.csv");
    correlation_csv<<"parameter_i,parameter_j,correlation\n";
    for(int i=0;i<np;++i) {
        parameters[i]=migration_fit.parameters[i];
        for(int j=0;j<np;++j) {
            covariance(i,j)=migration_fit.covariance[static_cast<size_t>(i)*np+j];
            const double vi=migration_fit.covariance[static_cast<size_t>(i)*np+i];
            const double vj=migration_fit.covariance[static_cast<size_t>(j)*np+j];
            correlation(i,j)=(!positivity_boundary_active && vi>0 && vj>0) ?
                covariance(i,j)/std::sqrt(vi*vj) : nan;
            correlation_csv<<i<<','<<j<<','<<correlation(i,j)<<'\n';
            // One shared divisor moves every coefficient coherently. Keep
            // this rank-one systematic separate from statistical/MC errors.
            target_covariance(i,j)=migration_fit.parameters[i]*migration_fit.parameters[j]*target_fraction*target_fraction;
            covariance_csv<<i<<','<<j<<','<<covariance(i,j)<<','<<target_covariance(i,j)<<'\n';
        }
    }
    auto singular_csv=csv("migration_singular_values.csv");
    singular_csv<<"index,whitened_column_normalized_singular_value\n";
    for(size_t i=0;i<migration_fit.singular_values.size();++i) {
        singular[i]=migration_fit.singular_values[i];singular_csv<<i<<','<<singular[i]<<'\n';
    }
    if (np > 0) {
        design.Write("migration_design_all_rows");
        parameters.Write("migration_parameters");
        covariance.Write(cfg.fit_objective=="scaled-poisson" ?
            "migration_covariance_scaled_poisson" : "migration_covariance_stat_plus_mc");
        correlation.Write("migration_parameter_correlation");
        target_covariance.Write("migration_covariance_target_correlated");
        singular.Write("migration_singular_values_scaled");
    }

    auto parameter_csv=csv("migration_parameters.csv");
    parameter_csv<<"parameter_index,active_block_index,truth_block,region,it,iq,ix,component,value,"
                 <<(cfg.fit_objective=="scaled-poisson" ? "error_conditional_Hessian" : "error_stat_plus_mc")
                 <<",error_target,is_nuisance\n";
    TTree parameter_tree("migration_parameter_index","Column mapping and global-fit parameters including free guard coefficients");
    int column=0,active=0,block=0,it=0,iq=0,ix=0,component=0;
    std::string region,component_name;
    int is_nuisance=0;
    double value=0,error=0,target_error=0;
    parameter_tree.Branch("parameter_index",&column);parameter_tree.Branch("active_block_index",&active);
    parameter_tree.Branch("truth_block",&block);parameter_tree.Branch("region",&region);
    parameter_tree.Branch("it",&it);parameter_tree.Branch("iq",&iq);parameter_tree.Branch("ix",&ix);
    parameter_tree.Branch("component",&component);parameter_tree.Branch("component_name",&component_name);
    parameter_tree.Branch("value",&value);parameter_tree.Branch("error_stat_plus_mc",&error);
    parameter_tree.Branch("error_target",&target_error);
    parameter_tree.Branch("is_nuisance",&is_nuisance);
    const char* components[]={"U","LT","TT"};
    for(column=0;column<np;++column) {
        active=column/3;component=column%3;block=active_truth_blocks[active];
        region=label(block);component_name=components[component];
        const auto index=indices(block);it=index[0];iq=index[1];ix=index[2];
        is_nuisance=block>=published || !retained_fit_groups[block%(cfg.n_q2*cfg.n_xb)];
        value=parameters[column];error=std::sqrt(covariance(column,column));target_error=std::abs(value)*target_fraction;
        parameter_tree.Fill();
        parameter_csv<<column<<','<<active<<','<<block<<','<<region<<','<<it<<','<<iq<<','<<ix<<','
                     <<component_name<<','<<value<<','<<error<<','<<target_error<<','<<is_nuisance<<'\n';
    }
    parameter_tree.Write();

    // Save raw weighted moments AND their means, including all six exterior
    // regions. Means describe the response-weighted generated events reaching
    // accepted reconstructed bins, not a bin-centering correction.
    auto truth_csv=csv("migration_truth_blocks.csv");
    truth_csv<<"truth_block,active_block_index,region,it,iq,ix,events,response_weight,sum_q2,sum_xb,sum_t,sum_tprime,sum_epsilon,mean_q2,mean_xb,mean_t,mean_tprime,mean_epsilon,epsilon_max\n";
    TTree truth_tree("migration_truth_blocks","Generated-origin weighted moments; inactive guard means are NaN");
    Long64_t events=0;double weight=0,sumq=0,sumx=0,sumt=0,sumtp=0,sume=0,mq=0,mx=0,mt=0,mtp=0,me=0;
    truth_tree.Branch("truth_block",&block);truth_tree.Branch("active_block_index",&active);
    truth_tree.Branch("region",&region);truth_tree.Branch("it",&it);truth_tree.Branch("iq",&iq);truth_tree.Branch("ix",&ix);
    truth_tree.Branch("events",&events);truth_tree.Branch("response_weight",&weight);
    truth_tree.Branch("sum_q2",&sumq);truth_tree.Branch("sum_xb",&sumx);truth_tree.Branch("sum_t",&sumt);
    truth_tree.Branch("sum_tprime",&sumtp);truth_tree.Branch("sum_epsilon",&sume);
    truth_tree.Branch("mean_q2",&mq);truth_tree.Branch("mean_xb",&mx);truth_tree.Branch("mean_t",&mt);
    truth_tree.Branch("mean_tprime",&mtp);truth_tree.Branch("mean_epsilon",&me);
    double max_eps=0;truth_tree.Branch("epsilon_max",&max_eps);
    for(block=0;block<nb;++block) {
        const auto found=std::find(active_truth_blocks.begin(),active_truth_blocks.end(),block);
        active=found==active_truth_blocks.end()?-1:static_cast<int>(found-active_truth_blocks.begin());
        region=label(block);const auto index=indices(block);it=index[0];iq=index[1];ix=index[2];
        const auto& m=truth_moments[block];events=m.events;weight=m.weight;
        sumq=m.q2;sumx=m.xb;sumt=m.t;sumtp=m.tprime;sume=m.epsilon;
        mq=weight>0?sumq/weight:nan;mx=weight>0?sumx/weight:nan;mt=weight>0?sumt/weight:nan;
        mtp=weight>0?sumtp/weight:nan;me=weight>0?sume/weight:nan;max_eps=m.epsilon_max;truth_tree.Fill();
        truth_csv<<block<<','<<active<<','<<region<<','<<it<<','<<iq<<','<<ix<<','<<events<<','<<weight<<','
                 <<sumq<<','<<sumx<<','<<sumt<<','<<sumtp<<','<<sume<<','<<mq<<','<<mx<<','<<mt<<','<<mtp<<','<<me<<','<<max_eps<<'\n';
    }
    truth_tree.Write();

    // A cell holds sum(event Fourier vectors) and sum(their outer products).
    // Keeping all 3x3 terms is essential: independent U/LT/TT MC errors would
    // destroy their event-level correlations. Under the Poissonized MC model,
    // distinct reconstructed-row/truth-block cells are independent.
    auto cell_csv=csv("migration_response_cells.csv");
    cell_csv<<"reco_row,truth_block,events,basis_U,basis_LT,basis_TT,cov_U_U,cov_U_LT,cov_U_TT,cov_LT_U,cov_LT_LT,cov_LT_TT,cov_TT_U,cov_TT_LT,cov_TT_TT\n";
    TTree cell_tree("migration_response_cells","Integrated basis and Poissonized event outer products for every row/truth cell");
    int row=0;double basis[3]{},cell_covariance[9]{};
    cell_tree.Branch("reco_row",&row);cell_tree.Branch("truth_block",&block);cell_tree.Branch("events",&events);
    cell_tree.Branch("basis",basis,"basis[3]/D");cell_tree.Branch("covariance",cell_covariance,"covariance[9]/D");
    for(row=0;row<nr;++row) for(block=0;block<nb;++block) {
        const auto& cell=migration_response[row][block];events=cell.events;
        std::copy(cell.basis.begin(),cell.basis.end(),basis);std::copy(cell.covariance.begin(),cell.covariance.end(),cell_covariance);
        cell_tree.Fill();cell_csv<<row<<','<<block<<','<<events;
        for(double v:basis) cell_csv<<','<<v;
        for(double v:cell_covariance) cell_csv<<','<<v;
        cell_csv<<'\n';
    }
    cell_tree.Write();

    // Four-term response for a multi-kinematic Rosenbluth fit. The first
    // column is T, the second is L and carries each event's vertex epsilon.
    // Store the full T/L/LT/TT outer product so the L column's finite-MC
    // covariance with the other Fourier terms is not discarded.
    auto joint_csv=csv("migration_joint_response_cells.csv");
    joint_csv<<"reco_row,truth_block,events,basis_T,basis_L,basis_LT,basis_TT";
    const std::array<const char*,4> term_names{{"T","L","LT","TT"}};
    for(const char* first:term_names) for(const char* second:term_names)
        joint_csv<<",cov_"<<first<<'_'<<second;
    joint_csv<<'\n';
    const double inv_twopi=1.0/(2.0*TMath::Pi());
    for(int jr=0;jr<nr;++jr) for(int jb=0;jb<nb;++jb) {
        const auto& cell=migration_response[jr][jb];
        const std::array<double,4> basis4{{cell.basis[0],cell.epsilon_weight*inv_twopi,
                                             cell.basis[1],cell.basis[2]}};
        const std::array<double,16> cov4{{
            cell.covariance[0],cell.basis_epsilon_covariance[0]*inv_twopi,cell.covariance[1],cell.covariance[2],
            cell.basis_epsilon_covariance[0]*inv_twopi,cell.epsilon_weight_covariance*inv_twopi*inv_twopi,
                cell.basis_epsilon_covariance[1]*inv_twopi,cell.basis_epsilon_covariance[2]*inv_twopi,
            cell.covariance[3],cell.basis_epsilon_covariance[1]*inv_twopi,cell.covariance[4],cell.covariance[5],
            cell.covariance[6],cell.basis_epsilon_covariance[2]*inv_twopi,cell.covariance[7],cell.covariance[8]
        }};
        joint_csv<<jr<<','<<jb<<','<<cell.events;
        for(double value:basis4) joint_csv<<','<<value;
        for(double value:cov4) joint_csv<<','<<value;
        joint_csv<<'\n';
    }

    // Exact last-solve variance is exported separately from the MC variance
    // recomputed at the final coefficients (equal within convergence tolerance).
    // Excluded rows have fit_index=-1 and variance_used=NaN, never pseudocounts.
    auto reco_csv=csv("migration_reco_rows.csv");
    reco_csv<<"reco_row,fit_index,it,iq,ix,ip,phi_lo,phi_hi,data,data_variance,variance_used,mc_variance_used,mc_variance_at_final,prediction,parameter_prediction_variance,residual,pull_conditional,has_response,exclusion\n";
    TTree reco_tree("migration_reco_rows","All reconstructed rows; exact final-solve inputs and exclusion mapping");
    int fit_index=0,ip=0,support=0;std::string exclusion;
    double plo=0,phi=0,data=0,data_variance=0,variance_used=0,mc_used=0,mc_final=0,prediction=0,parameter_variance=0,residual=0,pull=0;
    reco_tree.Branch("reco_row",&row);reco_tree.Branch("fit_index",&fit_index);
    reco_tree.Branch("it",&it);reco_tree.Branch("iq",&iq);reco_tree.Branch("ix",&ix);reco_tree.Branch("ip",&ip);
    reco_tree.Branch("phi_lo",&plo);reco_tree.Branch("phi_hi",&phi);
    reco_tree.Branch("data",&data);reco_tree.Branch("data_variance",&data_variance);
    reco_tree.Branch("variance_used",&variance_used);reco_tree.Branch("mc_variance_used",&mc_used);
    reco_tree.Branch("mc_variance_at_final",&mc_final);reco_tree.Branch("prediction",&prediction);
    reco_tree.Branch("parameter_prediction_variance",&parameter_variance);reco_tree.Branch("residual",&residual);
    reco_tree.Branch("pull_conditional",&pull);reco_tree.Branch("has_response",&support);reco_tree.Branch("exclusion",&exclusion);
    std::vector<int> fit_index_by_row(nr,-1);
    TVectorD used_rows(fit_rows.size()),used_data(fit_rows.size()),used_variance(fit_rows.size());
    for(size_t i=0;i<fit_rows.size();++i) {fit_index_by_row[fit_rows[i]]=static_cast<int>(i);used_rows[i]=fit_rows[i];used_variance[i]=fit_variance[i];}
    for(row=0;row<nr;++row) {
        const auto index=indices(row/cfg.n_phi);it=index[0];iq=index[1];ix=index[2];ip=row%cfg.n_phi;
        plo=phi_edges[ip];phi=phi_edges[ip+1];const auto& p=slices[row/cfg.n_phi].phi[ip];
        data=p.data;data_variance=p.data_sumw2;fit_index=fit_index_by_row[row];
        support=0;
        for (const auto& cell : migration_response[row])
            support=support || std::any_of(cell.basis.begin(),cell.basis.end(),[](double v){return v!=0;});
        const bool retained=retained_fit_groups[(row/cfg.n_phi)%(cfg.n_q2*cfg.n_xb)];
        exclusion=!retained?"excluded_q2_xb":(fit_index>=0?"included":(support?"zero_observed_variance":"no_response_and_no_data"));
        if (cfg.fit_objective=="scaled-poisson" && row<static_cast<int>(scaled_rows.size()))
            exclusion=scaled_rows[row].included ? "included" : scaled_rows[row].exclusion_reason;
        variance_used=(cfg.fit_objective=="scaled-poisson" || fit_index<0)?nan:fit_variance[fit_index];
        mc_used=(cfg.fit_objective=="scaled-poisson" || fit_index<0)?nan:variance_used-data_variance;
        mc_final=prediction=parameter_variance=nan;
        if (retained && np>0) {
            mc_final=nps_xsec::mc_prediction_variance(migration_response[row],active_truth_blocks,migration_fit.parameters);
            prediction=std::inner_product(response_design[row].begin(),response_design[row].end(),migration_fit.parameters.begin(),0.0);
            parameter_variance=0;
            for(int i=0;i<np;++i) for(int j=0;j<np;++j) parameter_variance+=design(row,i)*covariance(i,j)*design(row,j);
        }
        residual=data-prediction;
        pull=cfg.fit_objective=="scaled-poisson" ?
            (row<static_cast<int>(scaled_rows.size()) ? scaled_rows[row].residual : nan) :
            (fit_index>=0?residual/std::sqrt(variance_used):nan);
        if(fit_index>=0) used_data[fit_index]=data;
        reco_tree.Fill();
        reco_csv<<row<<','<<fit_index<<','<<it<<','<<iq<<','<<ix<<','<<ip<<','<<plo<<','<<phi<<','<<data<<','
                <<data_variance<<','<<variance_used<<','<<mc_used<<','<<mc_final<<','<<prediction<<','<<parameter_variance<<','
                <<residual<<','<<pull<<','<<support<<','<<exclusion<<'\n';
    }
    reco_tree.Write();used_rows.Write("migration_fit_row_indices");used_data.Write("migration_fit_data");used_variance.Write("migration_fit_variance");
    TParameter<int>("migration_mc_iterations",mc_iterations).Write();
    TParameter<int>("migration_mc_converged",(cfg.fit_objective=="scaled-poisson" || cfg.fit_variance_mode=="data")?-1:(mc_converged?1:0)).Write();
    TParameter<int>("migration_excluded_supported_zero_variance_rows",omitted_zero_variance_rows).Write();
    TParameter<int>("migration_nonphysical_truth_bins",nonphysical_truth_bins).Write();
    TParameter<int>("migration_rank",static_cast<int>(migration_fit.rank)).Write();
    TParameter<int>("migration_ndf",migration_fit.ndf).Write();
    TParameter<double>("migration_chi2",cfg.fit_objective=="scaled-poisson" ? nan : migration_fit.chi2).Write();
    TParameter<double>("migration_deviance",cfg.fit_objective=="scaled-poisson" ? migration_fit.chi2 : nan).Write();
    TParameter<int>("scaled_minuit_status",scaled_minuit_status).Write();
    TParameter<double>("scaled_minuit_edm",scaled_edm).Write();
    TParameter<int>("scaled_minuit_covariance_status",scaled_covariance_status).Write();
    TParameter<double>("migration_scaled_condition",migration_fit.condition).Write();
    // Detect full filesystem/write errors before reporting a successful export.
    for(auto* stream:{&status_csv,&attempts_csv,&response_csv,&covariance_csv,&correlation_csv,&singular_csv,&parameter_csv,&truth_csv,&cell_csv,&reco_csv,&curvature_csv,&positivity_csv}) {
        stream->flush();if(!*stream) die("Failed writing migration CSV output");
    }
}
