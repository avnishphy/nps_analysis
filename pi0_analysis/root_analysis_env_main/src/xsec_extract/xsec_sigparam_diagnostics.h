#pragma once
#include "xsec_analysis.h"

// Curves describe a declared fixed context; the detector fit always integrates
// actual events. Bin markers are response-weighted averages, not curve samples.
inline nps_xsec::XSecModelKinematics sigparam_curve_context(const nps_xsec::ProxyProblem& q,
        double tp,double mp,double mpi) {
    double wgt=0,q2=0,w2=0,eps=0;
    for(size_t b=0;b<q.blocks.size();++b)if(q.physical[b])for(const auto* e:q.reporting[b]) {
        wgt+=e->weight;q2+=e->weight*e->kinematics.Q2;
        w2+=e->weight*e->kinematics.W2;eps+=e->weight*e->kinematics.epsilon;
    }
    if(!(wgt>0))throw std::runtime_error("No physical curve context");
    q2/=wgt;w2/=wgt;eps/=wgt;
    const double tf=nps_xsec::forward_t(q2,std::sqrt(w2),mp,mpi);
    return nps_xsec::sigparam2021::kinematics(q2,std::sqrt(w2),tf+tp,tp,eps,mp,mpi);
}
inline nps_xsec::ModelEvaluation ExclPi0XSecAnalysis::model_curve(double tp) const {
    const auto x=sigparam_curve_context(model_all_rows,tp,cfg.mp,cfg.mpi0);
    return nps_xsec::evaluate_cached_model(nps_xsec::xsec_model().baseline(x).unseparated(),proxy_result.parameters.data(),x.tau,model_all_rows.tau0);
}
inline void ExclPi0XSecAnalysis::write_sigparam_diagnostics() {
    using namespace nps_xsec;
    auto csv=[&](const char* name){std::ofstream f(fs::path(cfg.out_dir)/name);f<<std::setprecision(17);return f;};
    auto metadata=csv("model_context.csv");
    metadata<<"tau0_GeV2,pivot_weighting,U_slope_fixed,free_physics_parameters,synthetic,internal_structure_unit,display_structure_unit,display_scale\n"
            <<model_all_rows.tau0<<",positive_physical_response,"<<proxy_options.fix_u_slope<<','
            <<kModelParameters-int(proxy_options.fix_u_slope)<<','<<synthetic_validation
            <<",microbarn/MeV2,nb/GeV2,1000000000\n";
    auto events=csv("model_event_cache.csv");
    auto spread=csv("model_charged_spread.csv");
    spread<<"truth_block,tprime,Q2,W2,epsilon,component,plus,minus,D,status\n";
    events<<"row,truth_block,truth_phi,physical,treatment,response_weight,Q2,W2,t,tprime,tau,theta_cm,epsilon,phi,baseline_T,baseline_L,baseline_U,baseline_LT,baseline_TT,basis_U,basis_LT,basis_TT\n";
    const std::array<const char*,9> names{"Q2","W","W2","abs_t","tprime","tau","theta_cm","epsilon","f_L"};
    std::array<double,9> lo,hi;lo.fill(std::numeric_limits<double>::infinity());hi.fill(-std::numeric_limits<double>::infinity());
    double lsum=0,usum=0,wgt=0;size_t physical_count=0,dominant=0;
    for(size_t b=0;b<model_all_rows.blocks.size();++b)for(const auto* e:model_all_rows.reporting[b]) {
        const auto& x=e->kinematics;const auto& f=e->baseline;
        const char* treatment=model_all_rows.is_physics(b)?"physics_model":model_all_rows.is_nuisance(b)?"fitted_tprime_feedin":"fixed_model_feedin";
        events<<e->row<<','<<e->block<<','<<e->truth_phi<<','<<model_all_rows.physical[b]<<','<<treatment<<','<<e->weight<<','
              <<x.Q2<<','<<x.W2<<','<<x.t<<','<<x.tprime<<','<<x.tau<<','<<x.theta_cm<<','<<x.epsilon<<','<<e->phi<<','
              <<f.sigma_T<<','<<f.sigma_L<<','<<f.sigma_U<<','<<f.sigma_LT<<','<<f.sigma_TT;
        for(double v:e->basis)events<<','<<v;events<<'\n';
        if(!model_all_rows.physical[b])continue;
        const auto plus=sigparam2021::components(x,sigparam2021::pp).unseparated();
        const auto minus=sigparam2021::components(x,sigparam2021::pm).unseparated();
        for(int k=0;k<3;++k) {
            const double den=.5*(std::abs(plus[k])+std::abs(minus[k]));
            const bool stable=den>1e-12*std::max({std::abs(plus[0]),std::abs(minus[0]),1e-30});
            spread<<e->block<<','<<x.tprime<<','<<x.Q2<<','<<x.W2<<','<<x.epsilon<<','<<(k==0?"U":k==1?"LT":"TT")<<','
                <<plus[k]<<','<<minus[k]<<','<<(stable?(plus[k]-minus[k])/den:std::numeric_limits<double>::quiet_NaN())<<','
                <<(stable?"ok":"near_zero_denominator")<<'\n';
        }
        ++physical_count;const double frac=x.epsilon*f.sigma_L/f.sigma_U;if(frac>.5)++dominant;
        const std::array<double,9> values{x.Q2,std::sqrt(x.W2),x.W2,std::abs(x.t),x.tprime,x.tau,x.theta_cm,x.epsilon,frac};
        for(int i=0;i<9;++i){lo[i]=std::min(lo[i],values[i]);hi[i]=std::max(hi[i],values[i]);}
        lsum+=e->weight*x.epsilon*f.sigma_L;usum+=e->weight*f.sigma_U;wgt+=e->weight;
    }
    auto ranges=csv("model_kinematic_ranges.csv");ranges<<"quantity,min,max\n";
    for(int i=0;i<9;++i){ranges<<names[i]<<','<<lo[i]<<','<<hi[i]<<'\n';std::cout<<"[MODEL_RANGE] "<<names[i]<<' '<<lo[i]<<' '<<hi[i]<<'\n';}
    auto longitudinal=csv("model_longitudinal_summary.csv");
    longitudinal<<"physical_events,events_L_fraction_above_half,fraction_above_half,response_weighted_L_fraction,T_L_status\n"
        <<physical_count<<','<<dominant<<','<<double(dominant)/physical_count<<','<<lsum/usum<<",model_assumed_not_Rosenbluth_separated\n";
    if(dominant)warn("PROVISIONAL pi0 longitudinal ansatz dominates U (f_L>0.5) at "+std::to_string(dominant)+" physical events; requires physics scrutiny");
    auto positivity=csv("model_event_positivity.csv");
    positivity<<"stage,points,negative_points,negative_fraction,minimum_response,dense_phi_minimum,Q2,W2,t,tprime,theta_cm,epsilon,phi_at_grid_minimum\n";
    for(bool fitted:{false,true}) {
        auto p=proxy_result.parameters;if(!fitted)std::copy(xsec_model().initial.begin(),xsec_model().initial.end(),p.begin());
        double worst=std::numeric_limits<double>::infinity(),gridworst=worst,gridphi=0;const ModelEvent* bad=nullptr;size_t neg=0;
        for(size_t b=0;b<model_all_rows.blocks.size();++b)if(model_all_rows.physical[b])for(const auto* e:model_all_rows.reporting[b]) {
            const auto f=model_all_rows.event_value(*e,p);const double eps=e->kinematics.epsilon;
            const double minimum=minimum_response(f.value[0],f.value[1],f.value[2],eps);
            if(minimum<0)++neg;if(minimum<worst){worst=minimum;bad=e;}
            for(int j=0;j<720;++j) {
                const double phi=2*TMath::Pi()*j/720.;
                const double v=f.value[0]+std::sqrt(2*eps*(1+eps))*std::cos(phi)*f.value[1]+eps*std::cos(2*phi)*f.value[2];
                if(v<gridworst){gridworst=v;gridphi=phi;}
            }
        }
        const auto& x=bad->kinematics;
        positivity<<(fitted?"fitted":"baseline")<<','<<physical_count<<','<<neg<<','<<double(neg)/physical_count<<','<<worst<<','<<gridworst<<','
            <<x.Q2<<','<<x.W2<<','<<x.t<<','<<x.tprime<<','<<x.theta_cm<<','<<x.epsilon<<','<<gridphi<<'\n';
    }
    auto grid=csv("model_baseline_diagnostics.csv");
    grid<<"tprime,tau,abs_t,Q2,W2,epsilon,theta_cm,kernel_charged_plus,kernel_neutral_plus,kernel_charged_minus,kernel_neutral_minus,charged_L_plus,charged_L_minus,old_pi0_L";
    for(const char* c:{"T","L","U","LT","TT"})grid<<','<<c<<"_plus,"<<c<<"_minus,"<<c<<"_baseline,"<<c<<"_fitted";
    grid<<",f_L\n";
    for(int i=0;i<=200;++i) {
        const double tp=cfg.tprime_min+(cfg.tprime_max-cfg.tprime_min)*i/200.;
        const auto x=sigparam_curve_context(model_all_rows,tp,cfg.mp,cfg.mpi0);
        const auto a=sigparam2021::components(x,sigparam2021::pp),b=sigparam2021::components(x,sigparam2021::pm),f=sigparam2021::baseline(x);
        const auto fit=model_curve(tp);const double factor=proxy_result.parameters[0]*std::exp(-proxy_result.parameters[1]*(x.tau-model_all_rows.tau0));
        grid<<tp<<','<<x.tau<<','<<std::abs(x.t)<<','<<x.Q2<<','<<x.W2<<','<<x.epsilon<<','<<x.theta_cm;
        for(const auto& params:{sigparam2021::pp,sigparam2021::pm})for(bool charged:{true,false})grid<<','<<sigparam2021::kernel(x.Q2,std::abs(x.t),params,charged);
        grid<<','<<sigparam2021::components(x,sigparam2021::pp,true).sigma_L<<','<<sigparam2021::components(x,sigparam2021::pm,true).sigma_L<<",0";
        const std::array<double,5> av{a.sigma_T,a.sigma_L,a.sigma_U,a.sigma_LT,a.sigma_TT},bv{b.sigma_T,b.sigma_L,b.sigma_U,b.sigma_LT,b.sigma_TT},fv{f.sigma_T,f.sigma_L,f.sigma_U,f.sigma_LT,f.sigma_TT},cv{factor*f.sigma_T,factor*f.sigma_L,fit.value[0],fit.value[1],fit.value[2]};
        for(int k=0;k<5;++k)grid<<','<<av[k]<<','<<bv[k]<<','<<fv[k]<<','<<cv[k];
        grid<<','<<x.epsilon*f.sigma_L/f.sigma_U<<'\n';
    }
    auto jac=csv("model_row_jacobian.csv");jac<<"row,parameter,derivative\n";
    for(size_t r=0;r<model_rows.size();++r)for(size_t j=0;j<proxy_result.parameters.size();++j)jac<<r<<','<<j<<','<<model_rows[r].jacobian[j]<<'\n';
}
