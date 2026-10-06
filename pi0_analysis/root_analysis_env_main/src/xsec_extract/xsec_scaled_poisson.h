#pragma once

// Scaled-Poisson approximation to a fixed-weight compound-Poisson yield.
// The upstream pi0_weight correction is accepted as fixed. For a populated
// row s=sum(w^2)/sum(w); zero rows borrow s from independently populated
// neighboring rows. No fractional weight is treated as an integer count.
// Near mu=Y, D=(mu-Y)^2/(s*Y)+O((mu-Y)^3), and s*Y=sum(w^2).
#include "xsec_analysis.h"
#include "xsec_scaled_poisson_stat.h"
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>

inline void ExclPi0XSecAnalysis::fit_scaled_poisson_subset(const std::vector<bool>& groups) {
    const size_t nr=response_design.size(), np=3*active_truth_blocks.size();
    const int ngroups=cfg.n_q2*cfg.n_xb;
    scaled_rows.assign(nr,XsecScaledRow{});
    fit_rows.clear(); fit_variance.clear();
    std::vector<XsecWeightMoments> blocks(slices.size()), slices_pooled(ngroups);
    XsecWeightMoments global;
    for (size_t b=0;b<slices.size();++b)
        for (const auto& p:slices[b].phi)
            blocks[b].merge(p.weights);
    for (size_t b=0;b<slices.size();++b)
        slices_pooled[b%ngroups].merge(blocks[b]);
    for (const auto& block:blocks) global.merge(block);
    if (!(global.scale()>0.)) die("No positive selected-data weight scale");
    if (global.min<0.) die("Scaled-Poisson cannot use signed selected event weights");
    for (size_t b=0;b<blocks.size();++b)
        if (blocks[b].sum2>0. && blocks[b].max2/blocks[b].sum2>0.10)
            warn("Scaled-Poisson block "+std::to_string(b)+
                 " has one event contributing >10% of sumw2; approximation sensitivity should be checked");
    // Twenty effective events is a deliberately conservative minimum for a
    // borrowed scale. Populated rows always use their own observed moments.
    constexpr double min_reference_neff=20.;
    for (size_t r=0;r<nr;++r) {
        auto& record=scaled_rows[r];
        const size_t b=r/cfg.n_phi;
        const auto& p=slices[b].phi[r%cfg.n_phi];
        record.support=std::any_of(response_design[r].begin(),response_design[r].end(),
                                   [](double v){ return v!=0.; }) ||
                       (!event_model() && fixed_feedin_prediction[r]!=0.);
        record.s_observed=p.data>0. ? p.data_sumw2/p.data :
            std::numeric_limits<double>::quiet_NaN();
        if (!groups[b%ngroups]) { record.exclusion_reason="excluded_kinematic_group"; continue; }
        if (!record.support) {
            record.exclusion_reason=p.data>0. ? "data_without_SIMC_support" : "zero_data_zero_SIMC_support";
            if (p.data>0.) die("Data outside MC support in row "+std::to_string(r));
            continue;
        }
        if (!(std::isfinite(p.data) && p.data>=0. && std::isfinite(p.data_sumw2) && p.data_sumw2>=0.))
            die("Scaled-Poisson requires nonnegative finite selected weighted yields");
        if (p.data>0.) {
            if (!(record.s_observed>0.)) die("Positive weighted yield has no sumw2 in row "+std::to_string(r));
            record.s_used=record.s_observed;
            record.scale_source="row";
            const double neff=p.data/record.s_used;
            if (std::abs(neff*record.s_used-p.data)>1e-10*std::max(1.,p.data) ||
                std::abs(neff-p.data*p.data/p.data_sumw2)>1e-10*std::max(1.,neff))
                die("Effective-count identity failed in row "+std::to_string(r));
        } else {
            const auto& block=blocks[b];
            const auto& slice=slices_pooled[b%ngroups];
            if (cfg.scaled_empty_scale=="auto" && block.effective_n()>=min_reference_neff) {
                record.s_used=block.scale(); record.scale_source="same-block-phi-pooled";
            } else if (cfg.scaled_empty_scale!="global" && slice.effective_n()>=min_reference_neff) {
                record.s_used=slice.scale(); record.scale_source=cfg.scaled_empty_scale=="auto" ? "slice-fallback" : "slice-sensitivity";
            } else {
                record.s_used=global.scale(); record.scale_source=cfg.scaled_empty_scale=="auto" ? "global-fallback" : "global-sensitivity";
            }
        }
        if (!(std::isfinite(record.s_used) && record.s_used>0.))
            die("Invalid reference weight scale in row "+std::to_string(r));
        record.included=true;
        fit_rows.push_back(static_cast<int>(r));
    }
    if (model_fit_mode) { fit_proxy_subset(groups); return; }
    if (fit_rows.size()<=np) die("Scaled-Poisson fit has no positive degrees of freedom");

    // The constrained Gaussian solve supplies a feasible seed only; its
    // covariance is never used as the scaled-Poisson fit covariance.
    std::vector<std::vector<double>> xseed;
    std::vector<double> yseed,vseed,epsilon_max;
    for (int b:active_truth_blocks) epsilon_max.push_back(truth_moments[b].epsilon_max);
    for (int r:fit_rows) {
        const auto& p=slices[r/cfg.n_phi].phi[r%cfg.n_phi];
        const double s=scaled_rows[r].s_used;
        xseed.push_back(response_design[r]);
        yseed.push_back(p.data-fixed_feedin_prediction[r]);
        vseed.push_back(s*std::max(p.data,s));
    }
    auto seed=nps_xsec::solve_positive_response(xseed,yseed,vseed,epsilon_max,cfg.rank_tolerance).fit;
    std::vector<double> unit(np,1.), start(np,0.);
    for (size_t b=0;b<active_truth_blocks.size();++b) {
        for (int a=0;a<3;++a) {
            const size_t j=3*b+a;
            unit[j]=std::max({std::abs(seed.parameters[j]),
                              std::sqrt(std::max(0.,seed.covariance[j*np+j])),1e-30});
        }
        const double lt=seed.parameters[3*b+1],tt=seed.parameters[3*b+2];
        const double required=std::max(0.,-nps_xsec::minimum_response(0.,lt,tt,epsilon_max[b]));
        start[3*b]=std::sqrt(std::max(1e-8,(seed.parameters[3*b]-required)/unit[3*b]));
        start[3*b+1]=lt/unit[3*b+1];
        start[3*b+2]=tt/unit[3*b+2];
    }
    auto physical=[&](const double* z) {
        std::vector<double> p(np);
        for (size_t b=0;b<active_truth_blocks.size();++b) {
            const size_t j=3*b;
            p[j+1]=unit[j+1]*z[j+1];
            p[j+2]=unit[j+2]*z[j+2];
            const double required=std::max(0.,-nps_xsec::minimum_response(
                0.,p[j+1],p[j+2],epsilon_max[b]));
            p[j]=required+unit[j]*z[j]*z[j];
        }
        return p;
    };
    auto objective=[&](const double* z) {
        const auto p=physical(z);
        double total=0.;
        for (int r:fit_rows) {
            const auto& row=response_design[r];
            const double mu=fixed_feedin_prediction[r]+std::inner_product(row.begin(),row.end(),p.begin(),0.);
            if (!(std::isfinite(mu) && mu>0.)) return 1e100;
            const double y=slices[r/cfg.n_phi].phi[r%cfg.n_phi].data;
            total+=nps_xsec::scaled_poisson_deviance(y,mu,scaled_rows[r].s_used);
        }
        return std::isfinite(total) ? total : 1e100;
    };
    auto minimizer=std::unique_ptr<ROOT::Math::Minimizer>(
        ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
    if (!minimizer) die("Minuit2 minimizer unavailable");
    ROOT::Math::Functor functor(objective,static_cast<unsigned int>(np));
    minimizer->SetFunction(functor);
    minimizer->SetErrorDef(1.);
    minimizer->SetMaxFunctionCalls(250000);
    minimizer->SetMaxIterations(50000);
    // Minuit2 uses EDM < 0.001*tolerance*ErrorDef; 0.01 requests
    // an absolute EDM below 1e-5 for this deviance (ErrorDef=1).
    minimizer->SetTolerance(0.01);
    for (size_t j=0;j<np;++j) {
        const std::string name=std::string(j%3==0?"gap_":j%3==1?"LT_":"TT_")+
                               std::to_string(j/3);
        if (j%3==0)
            minimizer->SetLowerLimitedVariable(j,name,start[j],0.03,0.);
        else
            minimizer->SetVariable(j,name,start[j],0.03);
    }
    const bool ok=minimizer->Minimize();
    minimizer->Hesse();
    scaled_minuit_status=minimizer->Status();
    scaled_covariance_status=minimizer->CovMatrixStatus();
    scaled_edm=minimizer->Edm();
    scaled_calls=minimizer->NCalls();
    // Status 2 means forced positive-definite, not an accurate confidence
    // covariance. Do not publish it as a successful uncertainty calculation.
    if (!ok || scaled_minuit_status!=0 || scaled_covariance_status!=3 || !std::isfinite(minimizer->MinValue()))
        die("Scaled-Poisson Minuit2 failed: status "+std::to_string(scaled_minuit_status)+
            " ok="+std::to_string(ok)+" cov="+std::to_string(scaled_covariance_status)+
            " edm="+std::to_string(scaled_edm)+" value="+std::to_string(minimizer->MinValue())+
            " calls="+std::to_string(scaled_calls));
    const double* z=minimizer->X();
    migration_fit.parameters=physical(z);
    migration_fit.chi2=objective(z);
    migration_fit.ndf=static_cast<int>(fit_rows.size()-np);
    migration_fit.rank=seed.rank;
    migration_fit.condition=seed.condition;
    migration_fit.singular_values=seed.singular_values;

    // Transform Minuit's internal Hessian covariance to physical U,LT,TT.
    // The envelope derivative gives dUmin/dLT and dUmin/dTT at the angular
    // minimum. At a boundary its symmetric Hessian errors are not intervals.
    std::vector<double> jac(np*np,0.);
    positivity_boundary_tolerances.resize(active_truth_blocks.size());
    positivity_feasibility_tolerances.resize(active_truth_blocks.size());
    positivity_boundary_active=false;
    for (size_t b=0;b<active_truth_blocks.size();++b) {
        const size_t j=3*b;
        const double lt=migration_fit.parameters[j+1],tt=migration_fit.parameters[j+2];
        const double eps=epsilon_max[b];
        const double cosine=static_cast<double>(nps_xsec::positive_detail::minimum(
            0.L,lt,tt,eps).second);
        jac[j*np+j]=2.*unit[j]*z[j];
        jac[j*np+j+1]=-std::sqrt(2.*eps*(1.+eps))*cosine*unit[j+1];
        jac[j*np+j+2]=-eps*(2.*cosine*cosine-1.)*unit[j+2];
        jac[(j+1)*np+j+1]=unit[j+1];
        jac[(j+2)*np+j+2]=unit[j+2];
        const double scale=std::max({std::abs(migration_fit.parameters[j]),
                                     std::sqrt(2.*eps*(1.+eps))*std::abs(lt),
                                     eps*std::abs(tt),unit[j]});
        positivity_boundary_tolerances[b]=2e-8*scale;
        positivity_feasibility_tolerances[b]=2e-10*scale;
        const double minimum=nps_xsec::minimum_response(migration_fit.parameters[j],lt,tt,eps);
        if (minimum < -positivity_feasibility_tolerances[b])
            die("Scaled-Poisson parameterization violated angular positivity");
        if (minimum<=positivity_boundary_tolerances[b])
            positivity_boundary_active=true;
    }
    migration_fit.covariance.assign(np*np,0.);
    for (size_t i=0;i<np;++i) for (size_t j=0;j<np;++j)
        for (size_t k=0;k<np;++k) for (size_t l=0;l<np;++l)
            migration_fit.covariance[i*np+j]+=jac[i*np+k]*
                minimizer->CovMatrix(k,l)*jac[j*np+l];
    for (size_t i=0;i<np;++i)
        if (!(std::isfinite(migration_fit.covariance[i*np+i]) &&
              migration_fit.covariance[i*np+i]>=0.))
            die("Invalid Minuit physical Hessian covariance");
    fit_curvature_inverse=migration_fit.covariance;
    if (positivity_boundary_active)
        std::fill(migration_fit.covariance.begin(),migration_fit.covariance.end(),
                  std::numeric_limits<double>::quiet_NaN());

    // MINOS profiles physical LT and TT directly because they are linear
    // internal coordinates. U mixes the gap and interference coefficients;
    // profiling the gap alone would not be a U interval, so label it absent.
    std::ofstream profiles(fs::path(cfg.out_dir)/"scaled_poisson_profile_intervals.csv");
    profiles<<std::setprecision(std::numeric_limits<double>::max_digits10);
    profiles<<"truth_block,component,value,lower_error,upper_error,status\n";
    for (size_t b=0;b<active_truth_blocks.size();++b) {
        const int block=active_truth_blocks[b];
        if (block<0 || block>=static_cast<int>(slices.size())) continue;
        profiles<<block<<",U,"<<migration_fit.parameters[3*b]
                <<",nan,nan,nonlinear_boundary_profile_not_computed\n";
        for (int a=1;a<3;++a) {
            const size_t j=3*b+a;
            double low=0.,high=0.;
            const bool good=minimizer->GetMinosError(j,low,high);
            profiles<<block<<','<<(a==1?"LT":"TT")<<','<<migration_fit.parameters[j]
                    <<','<<(good?low*unit[j]:std::numeric_limits<double>::quiet_NaN())
                    <<','<<(good?high*unit[j]:std::numeric_limits<double>::quiet_NaN())
                    <<','<<(good?"MINOS_delta_deviance_1":"MINOS_failed")<<'\n';
        }
    }
    fit_variance.reserve(fit_rows.size());
    for (int r:fit_rows) {
        auto& record=scaled_rows[r];
        const auto& row=response_design[r];
        const double mu=fixed_feedin_prediction[r]+std::inner_product(row.begin(),row.end(),
                                           migration_fit.parameters.begin(),0.);
        if (!(mu>0.)) die("Scaled-Poisson final prediction is nonpositive");
        const double y=slices[r/cfg.n_phi].phi[r%cfg.n_phi].data;
        record.deviance=nps_xsec::scaled_poisson_deviance(y,mu,record.s_used);
        record.residual=(y>=mu?1.:-1.)*std::sqrt(record.deviance);
        fit_variance.push_back(record.s_used*std::max(y,mu)); // display only
    }
    mc_converged=true;
    mc_iterations=0;
    omitted_zero_variance_rows=0;
    positivity_iterations=static_cast<int>(scaled_calls);
    finalize_fit_subset(groups);
}

inline void ExclPi0XSecAnalysis::write_scaled_poisson_diagnostics() {
    if (cfg.fit_objective!="scaled-poisson") return;
    const auto nan=std::numeric_limits<double>::quiet_NaN();
    std::ofstream out(fs::path(cfg.out_dir)/"scaled_poisson_rows.csv");
    out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"row,it,iq,ix,ip,phi_low,phi_high,phi_center,n_data,SIMC_response_events,sumw,sumw2,s_observed,s_used,scale_source,n_eff,mu,deviance,deviance_residual,legacy_chi2_at_current_prediction,SIMC_support,included,zero_bin,exclusion_reason\n";
    for (size_t r=0;r<scaled_rows.size();++r) {
        const size_t b=r/cfg.n_phi,ip=r%cfg.n_phi;
        const auto& p=slices[b].phi[ip]; const auto& d=scaled_rows[r];
        long long response_events=0;
        for (const auto& cell:migration_response[r]) response_events+=cell.events;
        const double mu=slices[b].fit_xsec.ok ? p.sim : nan;
        const double legacy=p.data_sumw2>0. && std::isfinite(mu) ?
            (p.data-mu)*(p.data-mu)/p.data_sumw2 : nan;
        const double neff=p.data_sumw2>0. ? p.data*p.data/p.data_sumw2 : 0.;
        out<<r<<','<<b/(cfg.n_q2*cfg.n_xb)<<','<<(b/cfg.n_xb)%cfg.n_q2<<','
           <<b%cfg.n_xb<<','<<ip<<','<<phi_edges[ip]<<','<<phi_edges[ip+1]<<','
           <<.5*(phi_edges[ip]+phi_edges[ip+1])<<','<<p.n_data<<','<<response_events<<','<<p.data<<','
           <<p.data_sumw2<<','<<d.s_observed<<','<<d.s_used<<','
           <<d.scale_source<<','<<neff<<','<<mu<<','<<d.deviance<<','
           <<d.residual<<','<<legacy<<','<<d.support<<','<<d.included<<','
           <<(p.data==0.)<<','<<d.exclusion_reason<<'\n';
        if (d.included && p.data==0.)
            std::cout<<"Scaled-Poisson zero row "<<r<<": raw_n="<<p.n_data
                     <<", SIMC_events="<<response_events<<", s="<<d.s_used
                     <<", mu="<<mu<<", D="<<d.deviance
                     <<", included=1"<<std::endl;
    }
    std::ofstream weights(fs::path(cfg.out_dir)/"scaled_poisson_weight_groups.csv");
    weights<<std::setprecision(std::numeric_limits<double>::max_digits10);
    weights<<"group,index,n,mean,rms,min,max,sumw,sumw2,scale,n_eff,max_sumw2_fraction\n";
    auto write=[&](const std::string& group,int index,const XsecWeightMoments& m) {
        weights<<group<<','<<index<<','<<m.n<<','
               <<(m.n?m.sum/m.n:nan)<<','
               <<(m.n?std::sqrt(m.sum2/m.n):nan)<<','
               <<(m.n?m.min:nan)<<','<<(m.n?m.max:nan)<<','
               <<m.sum<<','<<m.sum2<<','<<m.scale()<<','<<m.effective_n()<<','
               <<(m.sum2>0.?m.max2/m.sum2:nan)<<'\n';
    };
    XsecWeightMoments global;
    for (size_t b=0;b<slices.size();++b) {
        XsecWeightMoments block;
        for (const auto& p:slices[b].phi) block.merge(p.weights);
        write("block",static_cast<int>(b),block); global.merge(block);
    }
    std::vector<XsecWeightMoments> phi(cfg.n_phi), slice(cfg.n_q2*cfg.n_xb);
    for (size_t b=0;b<slices.size();++b)
        for (int ip=0;ip<cfg.n_phi;++ip) {
            phi[ip].merge(slices[b].phi[ip].weights);
            slice[b%(cfg.n_q2*cfg.n_xb)].merge(slices[b].phi[ip].weights);
        }
    for (int ip=0;ip<cfg.n_phi;++ip) write("phi",ip,phi[ip]);
    for (size_t g=0;g<slice.size();++g) write("slice",static_cast<int>(g),slice[g]);
    for (const auto& entry:run_weight_moments) write("run",entry.first,entry.second);
    write("global",0,global);

    TGraph residual;
    residual.SetName("scaled_poisson_deviance_residual");
    residual.SetTitle("Signed scaled-Poisson deviance residual;reconstructed row;signed #sqrt{D_{r}}");
    for (size_t r=0;r<scaled_rows.size();++r)
        if (scaled_rows[r].included)
            residual.SetPoint(residual.GetN(),static_cast<double>(r),scaled_rows[r].residual);
    fout->cd();
    residual.Write();
    if (cfg.diagnostics && (cfg.write_png || cfg.write_pdf) && residual.GetN()>0) {
        TCanvas canvas("c_scaled_deviance_residual","Scaled-Poisson residual",1100,600);
        residual.SetMarkerStyle(20);
        residual.SetMarkerSize(0.8);
        residual.Draw("AP");
        TLine zero(0.,0.,static_cast<double>(scaled_rows.size()-1),0.);
        zero.SetLineStyle(2); zero.Draw();
        canvas.Update();
        if (cfg.write_png) canvas.SaveAs((fs::path(cfg.out_dir)/"scaled_poisson_residual.png").c_str());
        if (cfg.write_pdf) canvas.SaveAs((fs::path(cfg.out_dir)/"scaled_poisson_residual.pdf").c_str());
    }
}
