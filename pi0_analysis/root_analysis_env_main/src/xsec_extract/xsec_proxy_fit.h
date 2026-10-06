#pragma once
#include "xsec_analysis.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include "xsec_sigparam_diagnostics.h"

inline void ExclPi0XSecAnalysis::fit_proxy_subset(const std::vector<bool>& groups) {
    nps_xsec::ProxyProblem problem;
    problem.poisson=cfg.fit_objective=="scaled-poisson";
    problem.positive=cfg.positive_xsec;
    problem.blocks=active_truth_blocks;
    for(int b:active_truth_blocks) {
        const auto& m=truth_moments[b];
        const bool physical=b<int(slices.size()) && groups[b%(cfg.n_q2*cfg.n_xb)];
        const bool fixed=!physical && !nps_xsec::is_fitted_tprime_feedin(b,static_cast<int>(slices.size()));
        if(physical && !(m.weight>0)) die("No response for generated bin "+std::to_string(b));
        problem.physical.push_back(physical);
        problem.fixed.push_back(fixed);
        problem.tprime.push_back(m.weight>0?m.tprime/m.weight:0.);
        problem.epsilon.push_back(m.epsilon_max);
        if(physical) std::cout<<"[PROXY_POINT] block="<<b<<" tprime="<<problem.tprime.back()
                             <<" tau="<<-problem.tprime.back()<<'\n';
    }
    for(int r:fit_rows) {
        const auto& p=slices[r/cfg.n_phi].phi[r%cfg.n_phi];
        problem.design.push_back(response_design[r]);problem.y.push_back(p.data);
        problem.variance.push_back(p.data_sumw2);
        problem.scale.push_back(problem.poisson?scaled_rows[r].s_used:0.);
    }
    const auto data_variance=problem.variance;
    if(event_model()) {
        model_all_rows=problem;
        model_all_rows.design=response_design;
        model_all_rows.events.resize(response_design.size());
        model_all_rows.reporting.resize(problem.blocks.size());
        for(const auto& e:model_events) {
            const auto it=std::find(problem.blocks.begin(),problem.blocks.end(),e.block);
            if(it==problem.blocks.end())continue;
            model_all_rows.events[e.row].push_back(&e);
            model_all_rows.reporting[it-problem.blocks.begin()].push_back(&e);
        }
        model_all_rows.set_pivot();problem.tau0=model_all_rows.tau0;
        std::cout<<std::setprecision(17)<<"[MODEL_PIVOT] tau0="<<problem.tau0<<" GeV2; positive physical response weighted, fixed\n";
        problem.reporting=model_all_rows.reporting;
        for(int r:fit_rows)problem.events.push_back(model_all_rows.events[r]);
    }
    const auto seed=nps_xsec::proxy_seed(problem,proxy_options);
    std::cout<<"[PROXY_INITIAL]";for(double v:seed)std::cout<<' '<<v;std::cout<<'\n';
    std::ofstream starts(fs::path(cfg.out_dir)/"model_starts.csv");
    starts<<std::setprecision(17)<<"start,converged,status,covariance_status,edm,objective,mc_iterations,mc_converged";
    for(const auto& n:problem.names())starts<<','<<n;starts<<'\n';
    nps_xsec::ProxyResult best;std::vector<double> best_variance;int best_iteration=0;
    const bool staged=proxy_options.fit_strategy=="staged_feasible";
    if(staged) {
        proxy_options.staged_trace=(fs::path(cfg.out_dir)/"staged_solver_history.csv").string();
        std::ofstream strategy(fs::path(cfg.out_dir)/"model_fit_strategy.txt");
        strategy<<"staged_feasible: Gaussian central-fit diagnostic. Status 0 means constrained solver acceptance, not Minuit status. No HESSE errors or covariance.\n";
    }
    auto minimize=[&](const std::vector<double>& initial) {
        return staged?nps_xsec::minimize_staged_proxy(problem,proxy_options,initial):nps_xsec::minimize_proxy(problem,proxy_options,initial);
    };
    for(int attempt=0;attempt<proxy_options.starts;++attempt) {
        auto initial=seed;
        if(attempt>0) {
            initial[2]*=attempt%2?-1.:1.;
            initial[3]*=(attempt/2)%2?-1.:1.;
            if(attempt>=4) {
                initial[0]*=attempt%2?1.5:.7;
                if(!proxy_options.fix_u_slope)initial[1]+=attempt%2?-.5:.5;
            }
        }
        problem.variance=data_variance;
        auto result=minimize(initial);
        bool settled=problem.poisson || cfg.fit_variance_mode=="data";
        int iteration=0;
        while(!settled && result.converged && result.status==0 && iteration<cfg.mc_max_iterations) {
            ++iteration;
            auto next_variance=data_variance;
            for(size_t i=0;i<fit_rows.size();++i)
                next_variance[i]+=problem.event_row(i,result.parameters).mc_variance;
            const auto previous_variance=problem.variance;problem.variance=next_variance;
            auto next=minimize(result.parameters);
            double change=0.;
            for(size_t j=0;j<result.parameters.size();++j) {
                // With no boundary confidence errors, use the stricter
                // sufficient parameter test; the variance rule is unchanged.
                const double scale=staged?std::max(std::abs(next.parameters[j]),1e-30):std::max({std::abs(next.parameters[j]),next.errors[j],1e-30});
                change=std::max(change,std::abs(next.parameters[j]-result.parameters[j])/scale);
            }
            for(size_t i=0;i<fit_rows.size();++i)
                change=std::max(change,std::abs(next_variance[i]-previous_variance[i])/next_variance[i]);
            result=std::move(next);settled=change<cfg.mc_fit_tolerance;
            if(iteration==1 || iteration%20==0)std::cout<<"[PROXY_MC] start="<<attempt<<" iteration="<<iteration<<" max_change="<<change<<'\n';
        }
        starts<<attempt<<','<<result.converged<<','<<result.status<<','<<result.covariance_status<<','
              <<result.edm<<','<<result.objective<<','<<iteration<<','<<settled;
        for(double p:result.parameters)starts<<','<<p;starts<<'\n';starts.flush();
        std::cout<<"[PROXY_START] "<<attempt<<" status="<<result.status<<" cov="<<result.covariance_status
                 <<" objective="<<result.objective<<" EDM="<<result.edm<<" MC="<<iteration<<" settled="<<settled<<'\n';
        if(result.converged && result.status==0 && settled && result.objective<best.objective) {
            best=result;best_variance=problem.variance;best_iteration=iteration;
        }
    }
    if(best.parameters.empty()) die("No converged proxy start; see model_starts.csv");
    proxy_result=best;problem.variance=best_variance;fit_variance=best_variance;
    if(event_model()) {
        model_rows.clear();
        for(size_t r=0;r<response_design.size();++r)model_rows.push_back(model_all_rows.event_row(r,best.parameters));
    }
    mc_iterations=best_iteration;mc_converged=true;
    scaled_minuit_status=best.status;scaled_covariance_status=best.covariance_status;
    scaled_edm=best.edm;scaled_calls=best.calls;
    std::vector<std::vector<double>> jac;
    migration_fit.parameters=problem.coefficients(best.parameters,&jac);
    migration_fit.covariance=nps_xsec::proxy_covariance(jac,best.covariance);
    fit_curvature_inverse=migration_fit.covariance;
    migration_fit.chi2=best.objective;
    migration_fit.ndf=int(fit_rows.size())-int(best.parameters.size())+int(proxy_options.fix_u_slope);
    // Rank/condition of the actual model-yield Jacobian, not the larger
    // independent-bin design. This does not regularize weak directions.
    std::vector<std::vector<double>> yield_jac(problem.design.size(),std::vector<double>(best.parameters.size()));
    for(size_t r=0;r<yield_jac.size();++r)yield_jac[r]=problem.event_row(r,best.parameters).jacobian;
    if(proxy_options.fix_u_slope)for(auto& row:yield_jac)row.erase(row.begin()+1);
    auto rank_variance=best_variance;
    const auto mu=problem.prediction(best.parameters);
    if(problem.poisson) for(size_t i=0;i<mu.size();++i)rank_variance[i]=problem.scale[i]*mu[i];
    try {
        const auto rank=nps_xsec::solve_weighted_response(yield_jac,problem.y,rank_variance,cfg.rank_tolerance);
        migration_fit.rank=rank.rank;migration_fit.condition=rank.condition;migration_fit.singular_values=rank.singular_values;
    } catch(const std::exception& e) {
        warn(std::string("Proxy identifiability: ")+e.what());best.covariance_status=0;proxy_result.covariance_status=0;
    }
    positivity_boundary_active=best.boundary || best.covariance_status!=3;
    positivity_boundary_tolerances.assign(active_truth_blocks.size(),0.);
    positivity_feasibility_tolerances.assign(active_truth_blocks.size(),0.);
    if(positivity_boundary_active) {
        warn("Proxy covariance unavailable as confidence errors: parameter/angular boundary or invalid covariance; raw Hessian saved only as diagnostic.");
        std::fill(migration_fit.covariance.begin(),migration_fit.covariance.end(),std::numeric_limits<double>::quiet_NaN());
    }
    if(problem.poisson) {
        fit_variance.clear();
        for(size_t i=0;i<fit_rows.size();++i) {
            auto& row=scaled_rows[fit_rows[i]];
            row.deviance=nps_xsec::scaled_poisson_deviance(problem.y[i],mu[i],row.s_used);
            row.residual=(problem.y[i]>=mu[i]?1.:-1.)*std::sqrt(row.deviance);
            fit_variance.push_back(row.s_used*std::max(problem.y[i],mu[i]));
        }
    }
    scaled_covariance_status=proxy_result.covariance_status;
    write_proxy_diagnostics(problem);
    finalize_fit_subset(groups);
}

inline void ExclPi0XSecAnalysis::write_proxy_diagnostics(const nps_xsec::ProxyProblem& problem) {
    constexpr size_t nphysics=nps_xsec::kModelParameters;
    const auto& result=proxy_result;const size_t np=result.parameters.size();
    const double nan=std::numeric_limits<double>::quiet_NaN();
    auto csv=[&](const char* name) {std::ofstream f(fs::path(cfg.out_dir)/name);f<<std::setprecision(17);return f;};
    auto pars=csv("model_parameters.csv");
    pars<<"index,name,role,value,error,raw_hessian_error,at_lower_limit,poorly_constrained,lower_bound,upper_bound,at_bound,relative_error,fixed,unit\n";
    const auto names=problem.names();
    for(size_t i=0;i<np;++i) {
        const auto& definition=nps_xsec::xsec_model();
        const double lower=i<nphysics?definition.lower[i]:-std::numeric_limits<double>::infinity();
        const double upper=i<nphysics?definition.upper[i]:std::numeric_limits<double>::infinity();
        const bool at_lower=result.parameters[i]-lower<1e-6*(i==0?result.initial[0]:1.);
        const bool at_limit=at_lower || upper-result.parameters[i]<1e-6;
        const bool weak=result.errors[i]>=std::abs(result.parameters[i]);
        pars<<i<<','<<names[i]<<','<<(i<nphysics?"model":"migration_nuisance")<<','<<result.parameters[i]<<','
            <<(positivity_boundary_active?nan:result.errors[i])<<','<<result.errors[i]<<','<<at_lower<<','<<weak<<','
            <<lower<<','<<upper<<','<<at_limit<<','
            <<(positivity_boundary_active||result.parameters[i]==0?nan:result.errors[i]/std::abs(result.parameters[i]))<<','
            <<(i==1 && proxy_options.fix_u_slope)<<','<<(i==1?"GeV^-2":i<nphysics?"dimensionless":"microbarn/MeV2")<<'\n';
        std::cout<<"[PROXY_PARAMETER] "<<names[i]<<'='<<result.parameters[i]<<" error="
                 <<(positivity_boundary_active?nan:result.errors[i])<<" raw_Hesse="<<result.errors[i]
                 <<" lower_limit="<<at_limit<<" weak="<<weak<<'\n';
    }
    TMatrixDSym cov(np),corr(np);
    auto covcsv=csv("model_covariance.csv");covcsv<<"i,j,role,covariance,correlation,raw_hessian_covariance\n";
    for(size_t i=0;i<np;++i)for(size_t j=0;j<np;++j) {
        cov(i,j)=result.covariance[i*np+j];
        const double den=std::sqrt(std::abs(result.covariance[i*np+i]*result.covariance[j*np+j]));
        corr(i,j)=den>0?cov(i,j)/den:(proxy_options.fix_u_slope && (i==1 || j==1)?0.:nan);
        covcsv<<i<<','<<j<<','<<(i<nphysics?(j<nphysics?"model":"model_nuisance"):(j<nphysics?"model_nuisance":"nuisance"))<<','
              <<(positivity_boundary_active?nan:cov(i,j))<<','<<corr(i,j)<<','<<cov(i,j)<<'\n';
        if(i<j && std::abs(corr(i,j))>0.9)warn("Strong proxy correlation: "+names[i]+" / "+names[j]+" = "+std::to_string(corr(i,j)));
    }
    TMatrixDSym active_corr(np-int(proxy_options.fix_u_slope));
    for(size_t i=0,ii=0;i<np;++i)if(!(proxy_options.fix_u_slope && i==1)) {
        for(size_t j=0,jj=0;j<np;++j)if(!(proxy_options.fix_u_slope && j==1))active_corr(ii,jj++)=corr(i,j);
        ++ii;
    }
    double smallest=nan,largest=nan;
    bool finite_correlation=true;
    for(int i=0;i<active_corr.GetNrows();++i)for(int j=0;j<active_corr.GetNcols();++j)
        finite_correlation=finite_correlation&&std::isfinite(active_corr(i,j));
    if(finite_correlation) {
        TMatrixDSymEigen eigen(active_corr);const auto ev=eigen.GetEigenValues();
        smallest=largest=ev[0];for(int i=0;i<ev.GetNrows();++i){smallest=std::min(smallest,ev[i]);largest=std::max(largest,ev[i]);}
    }
    auto status=csv("model_fit_status.csv");
    status<<"converged,status,covariance_status,edm,objective,rows,parameters,nominal_dof,objective_per_dof,mc_iterations,mc_converged,boundary_or_invalid_covariance,jacobian_condition,correlation_min_eigenvalue,correlation_condition,model_identifier,statistical_objective,physics_parameters,nuisance_parameters\n"
          <<result.converged<<','<<result.status<<','<<result.covariance_status<<','<<result.edm<<','<<result.objective<<','
          <<fit_rows.size()<<','<<np<<','<<migration_fit.ndf<<','<<result.objective/migration_fit.ndf<<','<<mc_iterations<<','
          <<mc_converged<<','<<positivity_boundary_active<<','<<migration_fit.condition<<','<<smallest<<','<<largest/smallest
          <<','<<nps_xsec::xsec_model().id<<','<<cfg.fit_objective<<','<<nphysics<<','<<np-nphysics<<'\n';
    if(!finite_correlation || !(smallest>0) || largest/smallest>1e8)warn("Proxy correlation matrix has near-singular or invalid directions");
    int negative_rows=0,negative_blocks=0;
    for(const auto& row:response_design)if(std::inner_product(row.begin(),row.end(),migration_fit.parameters.begin(),0.)<0)++negative_rows;
    if(event_model()){negative_rows=0;for(const auto& row:model_rows)if(row.prediction<0)++negative_rows;}
    for(size_t b=0;b<problem.blocks.size();++b)if(!problem.is_fixed(b) && nps_xsec::minimum_response(migration_fit.parameters[3*b],migration_fit.parameters[3*b+1],migration_fit.parameters[3*b+2],problem.epsilon[b])<0)++negative_blocks;
    if(event_model()) {
        negative_blocks=0;
        for(size_t b=0;b<problem.blocks.size();++b) {
            bool negative=false;
            for(const auto* e:problem.reporting[b]) {
                const auto f=problem.event_value(*e,result.parameters);
                const double eps=problem.is_physics(b)?e->kinematics.epsilon:problem.epsilon[b];
                negative=negative || (!problem.is_fixed(b) && nps_xsec::minimum_response(f.value[0],f.value[1],f.value[2],eps)<0);
            }
            if(negative)++negative_blocks;
        }
    }
    auto validity=csv("model_prediction_validity.csv");
    validity<<"negative_reconstructed_rows,negative_angular_blocks,positivity_required\n"<<negative_rows<<','<<negative_blocks<<','<<(problem.positive||problem.poisson)<<'\n';
    if(negative_rows || negative_blocks)warn("Gaussian proxy has "+std::to_string(negative_rows)+" negative reconstructed predictions and "+std::to_string(negative_blocks)+" angular-negative fitted truth blocks; no clipping applied.");
    auto points=csv("model_structure_functions.csv");
    points<<"truth_block,tprime,tau,sigma_U,sigma_U_error,sigma_LT,sigma_LT_error,sigma_TT,sigma_TT_error,it,iq,ix\n";
    for(size_t b=0;b<problem.blocks.size();++b)if(problem.physical[b]) {
        points<<problem.blocks[b]<<','<<problem.tprime[b]<<','<<-problem.tprime[b];
        for(int k=0;k<3;++k)points<<','<<migration_fit.parameters[3*b+k]<<','<<
            (positivity_boundary_active?nan:std::sqrt(migration_fit.covariance[(3*b+k)*migration_fit.parameters.size()+3*b+k]));
        const int block=problem.blocks[b];
        points<<','<<block/(cfg.n_q2*cfg.n_xb)<<','<<(block/cfg.n_xb)%cfg.n_q2<<','<<block%cfg.n_xb<<'\n';
    }
    auto folded=csv("model_reconstructed_yields.csv");
    auto fixed_feed=csv("model_fixed_feedin.csv");
    fixed_feed<<"row,prediction,poissonized_mc_variance,events\n";
    for(size_t r=0;r<model_all_rows.events.size();++r) {
        double value=0.,variance=0.;size_t events=0;
        for(const auto* e:model_all_rows.events[r]) {
            const auto found=std::find(problem.blocks.begin(),problem.blocks.end(),e->block);
            if(found==problem.blocks.end() || !problem.is_fixed(found-problem.blocks.begin()))continue;
            double event=0.;const auto baseline=e->baseline.unseparated();
            for(int k=0;k<3;++k)event+=e->basis[k]*baseline[k];
            value+=event;variance+=event*event;++events;
        }
        fixed_feed<<r<<','<<value<<','<<variance<<','<<events<<'\n';
    }
    auto derived=csv("model_structure_covariance.csv");
    derived<<"truth_block_i,component_i,truth_block_j,component_j,covariance,raw_hessian_propagation\n";
    const size_t nc=migration_fit.parameters.size();
    for(size_t a=0;a<problem.blocks.size();++a)if(problem.physical[a])
        for(size_t b=0;b<problem.blocks.size();++b)if(problem.physical[b])
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)
                derived<<problem.blocks[a]<<','<<i<<','<<problem.blocks[b]<<','<<j<<','
                       <<migration_fit.covariance[(3*a+i)*nc+3*b+j]<<','<<fit_curvature_inverse[(3*a+i)*nc+3*b+j]<<'\n';
    folded<<"row,included,data,data_sumw2,prediction,residual,objective_pull,objective_variance,it,iq,ix,ip,phi_lo,phi_hi,component_U,component_LT,component_TT,prediction_error\n";
    TGraphErrors data,prediction;TGraph pull;
    for(size_t r=0;r<response_design.size();++r) {
        const auto& p=slices[r/cfg.n_phi].phi[r%cfg.n_phi];
        const double mu=event_model()?model_rows[r].prediction:std::inner_product(response_design[r].begin(),response_design[r].end(),migration_fit.parameters.begin(),0.);
        const auto it=std::find(fit_rows.begin(),fit_rows.end(),int(r));const bool included=it!=fit_rows.end();
        const size_t i=it-fit_rows.begin();
        const double residual=included?(problem.poisson?scaled_rows[r].residual:(p.data-mu)/std::sqrt(fit_variance[i])):nan;
        const int slice_id=r/cfg.n_phi, ip=r%cfg.n_phi;
        std::array<double,3> component{{0.,0.,0.}};
        for(size_t j=0;j<nc;++j)component[j%3]+=response_design[r][j]*migration_fit.parameters[j];
        double prediction_variance=0.;
        for(size_t a=0;a<nc;++a)for(size_t b=0;b<nc;++b)
            prediction_variance+=response_design[r][a]*migration_fit.covariance[a*nc+b]*response_design[r][b];
        if(event_model()){component=model_rows[r].components;prediction_variance=model_row_variance(r);}
        folded<<r<<','<<included<<','<<p.data<<','<<p.data_sumw2<<','<<mu<<','<<p.data-mu<<','<<residual<<','<<(included?fit_variance[i]:nan)
              <<','<<slice_id/(cfg.n_q2*cfg.n_xb)<<','<<(slice_id/cfg.n_xb)%cfg.n_q2<<','<<slice_id%cfg.n_xb<<','<<ip
              <<','<<phi_edges[ip]<<','<<phi_edges[ip+1]<<','<<component[0]<<','<<component[1]<<','<<component[2]<<','
              <<(positivity_boundary_active||prediction_variance<0?nan:std::sqrt(prediction_variance))<<'\n';
        data.SetPoint(r,r,p.data);data.SetPointError(r,0,std::sqrt(p.data_sumw2));prediction.SetPoint(r,r,mu);
        if(included)pull.SetPoint(pull.GetN(),r,residual);
    }
    // Export model evaluations once, using the same model/Jacobian as the fit.
    // Plotting consumers need no knowledge of the parameterization.
    auto curves=csv("model_structure_curves.csv");
    curves<<"tprime,tau,component,value,error\n";
    for(int i=0;i<=200;++i) {
        const double tp=cfg.tprime_min+(cfg.tprime_max-cfg.tprime_min)*i/200.;
        const auto f=model_curve(tp);
        for(int k=0;k<3;++k) {
            double variance=0.;
            for(size_t a=0;a<nphysics;++a)for(size_t b=0;b<nphysics;++b)variance+=f.jacobian[k][a]*cov(a,b)*f.jacobian[k][b];
            curves<<tp<<','<<-tp<<','<<(k==0?"U":k==1?"LT":"TT")<<','<<f.value[k]<<','
                  <<(positivity_boundary_active||variance<0?nan:std::sqrt(variance))<<'\n';
        }
    }
    if(event_model())write_sigparam_diagnostics();
    const size_t fixed_blocks=std::count(problem.fixed.begin(),problem.fixed.end(),true);
    std::cout<<"[MODEL_SUMMARY] model="<<nps_xsec::xsec_model().id<<" physics_parameters="<<nphysics<<" nuisance_parameters="<<np-nphysics
             <<" fixed_feedin_blocks="<<fixed_blocks<<" total_fit_parameters="<<np
             <<" supported_rows="<<fit_rows.size()<<" objective="<<cfg.fit_objective<<'\n';
    // Separate ROOT file retains raw covariance and an explicit validity flag.
    TFile output((fs::path(cfg.out_dir)/"model_fit.root").c_str(),"RECREATE");
    cov.Write("raw_parameter_hessian_covariance");corr.Write("parameter_correlation");
    TParameter<int>("covariance_trustworthy",!positivity_boundary_active).Write();
    TParameter<int>("minuit_covariance_status",result.covariance_status).Write();
    TParameter<double>("tau0_GeV2",problem.tau0).Write();
    TParameter<int>("U_slope_fixed",proxy_options.fix_u_slope).Write();
    TVectorD fitted(np);for(size_t i=0;i<np;++i)fitted[i]=result.parameters[i];fitted.Write("parameters");
    auto save=[&](TCanvas& canvas,const char* name){
        canvas.Write(name);
        const size_t pages_before=combined_pdf_pages.size();
        write_canvas_pdf_png(&canvas,(fs::path(cfg.out_dir)/name).string());
        // Keep legacy standalone previews and ROOT canvases. The production
        // report supplies richer versions; avoid duplicate combined-PDF pages.
        combined_pdf_pages.resize(pages_before);
    };
    TCanvas yields("proxy_yields","Detector-level closure",1100,800);yields.Divide(1,2);
    yields.cd(1);data.SetTitle("Forward-folded model and measured yield;reconstructed row;yield / mC");data.SetMarkerStyle(20);data.Draw("AP");
    prediction.SetLineColor(kRed);prediction.Draw("L SAME");
    TLegend legend(.65,.73,.89,.89);legend.AddEntry(&data,"Measured weighted yield","pe");legend.AddEntry(&prediction,"Forward-folded proxy","l");legend.Draw();
    yields.cd(2);pull.SetTitle("Objective residuals;reconstructed row;pull / signed sqrt(deviance)");pull.SetMarkerStyle(20);pull.Draw("AP");save(yields,"model_yield_closure");
    for(int k=0;k<3;++k) {
        TGraphErrors curve,reported;
        for(int i=0;i<=200;++i) {
            const double tp=cfg.tprime_min+(cfg.tprime_max-cfg.tprime_min)*i/200.;
            const auto f=model_curve(tp);double variance=0.;
            for(size_t a=0;a<nphysics;++a)for(size_t b=0;b<nphysics;++b)variance+=f.jacobian[k][a]*cov(a,b)*f.jacobian[k][b];
            curve.SetPoint(i,tp,1e9*f.value[k]);curve.SetPointError(i,0,positivity_boundary_active?0:1e9*std::sqrt(std::max(0.,variance)));
        }
        for(size_t b=0;b<problem.blocks.size();++b)if(problem.physical[b]) {
            int i=reported.GetN();reported.SetPoint(i,problem.tprime[b],1e9*migration_fit.parameters[3*b+k]);
            reported.SetPointError(i,0,positivity_boundary_active?0:1e9*std::sqrt(migration_fit.covariance[(3*b+k)*migration_fit.parameters.size()+3*b+k]));
        }
        const std::string name=std::string("model_sigma_")+(k==0?"U":k==1?"LT":"TT");
        TCanvas c(name.c_str(),name.c_str(),900,650);reported.SetTitle((name+(positivity_boundary_active?" (errors unavailable)":"")+";signed t' [GeV^{2}];#sigma [nb/GeV^{2}]").c_str());
        reported.SetMarkerStyle(20);reported.Draw("APZ");curve.SetLineColor(kBlue);curve.SetLineStyle(2);curve.SetLineWidth(1);curve.Draw("L SAME");
        TLegend terms(.15,.73,.85,.88);terms.AddEntry(&reported,"Event-averaged fitted model","pe");terms.AddEntry(&curve,"Reference kinematics - diagnostic only","l");terms.Draw();save(c,name.c_str());
    }
    output.Close();
}
