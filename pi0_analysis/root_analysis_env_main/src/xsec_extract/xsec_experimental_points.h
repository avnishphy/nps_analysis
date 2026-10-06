#pragma once

// Residual-corrected Eq. 5.30-5.31 points. A truth cell v=(block,vertex phi
// bin) inherits its block's U/LT/TT coefficients. The corresponding row r
// uses the same t'/Q2/xB/phi indices at reconstruction. Exterior and other
// truth cells are subtracted through the complete forward prediction.
#include "xsec_analysis.h"

inline void ExclPi0XSecAnalysis::compute_experimental_points() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int nr = static_cast<int>(migration_response.size());
    const int np = static_cast<int>(migration_fit.parameters.size());
    experimental_points.assign(nr, ExperimentalPoint{});
    experimental_point_covariance.assign(static_cast<size_t>(nr)*nr, nan);
    for(int r=0;r<nr;++r) { experimental_points[r].row=r; experimental_points[r].truth_cell=r; }
    if (np == 0) return;
    if(event_model()) {
        // The independent-bin residual-correction formula assumes constant
        // coefficients within a truth cell. It is not an event-model observable.
        // Actual folded predictions and generated-average SF remain available.
        for(auto& p:experimental_points)p.status="unavailable_event_model_bin_correction";
        return;
    }
    std::vector<int> active_index(truth_moments.size(), -1), fit_index(nr, -1);
    for (size_t i=0; i<active_truth_blocks.size(); ++i) active_index[active_truth_blocks[i]]=static_cast<int>(i);
    for (size_t i=0; i<fit_rows.size(); ++i) fit_index[fit_rows[i]]=static_cast<int>(i);
    std::vector<double> prediction(nr, nan);
    for (int r=0; r<nr; ++r)
        prediction[r]=fixed_feedin_prediction[r]+std::inner_product(response_design[r].begin(),response_design[r].end(),
                                         migration_fit.parameters.begin(),0.0);

    for (int r=0; r<nr; ++r) {
        auto& p=experimental_points[r];
        p.row=r; p.truth_cell=r;
        const int b=r/cfg.n_phi, ip=r%cfg.n_phi;
        if (!slices[b].fit_xsec.ok) { p.status="excluded_or_failed_fit"; continue; }
        const int a=active_index[b];
        if (a<0) { p.status="inactive_truth_block"; continue; }
        if (r>=static_cast<int>(truth_phi_response.size())) { p.status="no_vertex_phi_bookkeeping"; continue; }
        const auto found=truth_phi_response[r].find(r);
        if (found==truth_phi_response[r].end() || found->second.events==0) {
            p.status="no_diagonal_response"; continue;
        }
        const auto& tm=truth_phi_moments[r];
        if (!(tm.weight>0)) { p.status="no_truth_reference"; continue; }
        p.q2_reference=tm.q2/tm.weight;
        p.xb_reference=tm.xb/tm.weight;
        p.tprime_reference=tm.tprime/tm.weight;
        p.epsilon_reference=tm.epsilon/tm.weight;
        const double phi_ref=.5*(phi_edges[ip]+phi_edges[ip+1]);
        const auto gb=sigma_model_basis(phi_ref,p.epsilon_reference,0.,0);
        double m=0., f=0.;
        for (int n=0;n<3;++n) {
            m+=found->second.basis[n]*migration_fit.parameters[3*a+n];
            f+=gb[n]*migration_fit.parameters[3*a+n];
        }
        double abs_sum=0.;
        for (int k=0;k<np;++k)
            abs_sum+=std::abs(response_design[r][k]*migration_fit.parameters[k]);
        p.contribution=m;
        p.sigma_reference=f;
        if (!std::isfinite(m) || std::abs(m)<=1e-12*abs_sum || m==0.) {
            p.status="unstable_denominator"; continue;
        }
        const double y=slices[b].phi[ip].data;
        const double residual=y-prediction[r];
        p.subtracted_yield=residual+m;
        p.correction=p.subtracted_yield/m;
        p.sigma_exp=p.correction*f;
        if (!std::isfinite(p.sigma_exp)) { p.status="nonfinite_point"; continue; }
        // Forward-predicted reconstructed means use the same accepted MC
        // events and fitted angular weights as the reconstructed yield.
        const double fitted_prediction=prediction[r]-fixed_feedin_prediction[r];
        if (std::isfinite(fitted_prediction) && std::abs(fitted_prediction)>1e-12*abs_sum) {
            double q=0.,x=0.,tp=0.;
            for (int j=0;j<static_cast<int>(active_truth_blocks.size());++j) {
                const auto& cell=migration_response[r][active_truth_blocks[j]];
                for (int n=0;n<3;++n) {
                    const double v=migration_fit.parameters[3*j+n];
                    q+=v*cell.reco_q2[n]; x+=v*cell.reco_xb[n];
                    tp+=v*cell.reco_tprime[n];
                }
            }
            // These moments cover fitted-response components only; fixed
            // feed-in remains in the detector prediction but is not assigned
            // fictitious U/LT/TT coefficients for moment reconstruction.
            p.reco_prediction_q2=q/fitted_prediction;
            p.reco_prediction_xb=x/fitted_prediction;
            p.reco_prediction_tprime=tp/fitted_prediction;
        }
        p.status=model_fit_mode ? "central_only_proxy_diagnostic" : !cfg.joint_plot_input.empty() ? "central_only_joint" :
            cfg.fit_objective=="scaled-poisson" ? "central_only_scaled_poisson" :
            (positivity_boundary_active ? "central_only_boundary" : "ok");
    }
    // The single-setting influence sums below omit the other joint settings.
    // Keep these residual-corrected points explicitly central-only in joint
    // plots; coefficient/curve errors use the full imported joint covariance.
    if (model_fit_mode || !cfg.joint_plot_input.empty() || positivity_boundary_active || cfg.fit_objective=="scaled-poisson") return;

    // For fixed D,V: Xhat=K y, K=C D' V^-1. For E=f*(1+z/m), with
    // z=y_r-D_r Xhat, f=g_v Xhat, m=d_rv Xhat, the derivative with
    // respect to Xhat is h=(1+z/m)g-(f/m)D_r-(f*z/m^2)d_rv.
    // Hence dE/dy_s=(f/m)delta_rs+h'K_s. This keeps the same-data
    // numerator, denominator and fitted curve correlations exactly at
    // first order for the data-only reference mode.
    std::vector<std::vector<double>> data_gradient(nr,std::vector<double>(nr,0.));
    std::vector<std::vector<double>> ch(nr,std::vector<double>(np,0.));
    std::vector<double> point_f(nr,0.),point_m(nr,0.),point_z(nr,0.);
    for (int r=0;r<nr;++r) {
        if (experimental_points[r].status!="ok") continue;
        const auto& p=experimental_points[r];
        const int a=active_index[r/cfg.n_phi];
        const auto& cell=truth_phi_response[r].at(r);
        const auto gb=sigma_model_basis(.5*(phi_edges[r%cfg.n_phi]+phi_edges[r%cfg.n_phi+1]),
                                        p.epsilon_reference,0.,0);
        const double f=p.sigma_reference,m=p.contribution;
        const double z=slices[r/cfg.n_phi].phi[r%cfg.n_phi].data-prediction[r];
        point_f[r]=f;point_m[r]=m;point_z[r]=z;
        std::vector<double> h(np,0.);
        for (int k=0;k<np;++k) h[k]=-(f/m)*response_design[r][k];
        for (int n=0;n<3;++n)
            h[3*a+n]+=(1.+z/m)*gb[n]-(f*z/(m*m))*cell.basis[n];
        for (int j=0;j<np;++j)
            for (int k=0;k<np;++k)
                ch[r][j]+=migration_fit.covariance[static_cast<size_t>(j)*np+k]*h[k];
        data_gradient[r][r]+=f/m;
        for (size_t i=0;i<fit_rows.size();++i) {
            const int t=fit_rows[i];
            const double dot=std::inner_product(ch[r].begin(),ch[r].end(),
                                                response_design[t].begin(),0.);
            data_gradient[r][t]+=dot/fit_variance[i];
        }
    }
    for (int r=0;r<nr;++r) {
        if (experimental_points[r].status!="ok") continue;
        for (int s=0;s<=r;++s) {
            if (experimental_points[s].status!="ok") continue;
            double c=0.;
            for (int t=0;t<nr;++t)
                c+=data_gradient[r][t]*data_gradient[s][t]*
                   slices[t/cfg.n_phi].phi[t%cfg.n_phi].data_sumw2;
            experimental_point_covariance[static_cast<size_t>(r)*nr+s]=c;
            experimental_point_covariance[static_cast<size_t>(s)*nr+r]=c;
        }
    }
    if (cfg.fit_variance_mode=="finite-mc") {
        // First-order independent Poissonized MC-cell extension, conditional
        // on the converged row variances and frozen binning. Different vertex
        // phi cells have independent event sums; their 3x3 outer products keep
        // U/LT/TT event correlations. dX/dD_{t,b,n} =
        // C[e_{b,n} z_t/V_t - D_t' X_{b,n}/V_t] for fitted t.
        for (int t=0;t<nr;++t) for (const auto& entry:truth_phi_response[t]) {
            const int v=entry.first,b=v/cfg.n_phi,a=active_index[b];
            if (a<0 || entry.second.events==0) continue;
            std::vector<std::array<double,4>> gradient(nr);
            const int fi=fit_index[t];
            const double residual=fi>=0 ?
                slices[t/cfg.n_phi].phi[t%cfg.n_phi].data-prediction[t] : 0.;
            for (int r=0;r<nr;++r) {
                if (experimental_points[r].status!="ok") continue;
                const double f=point_f[r],m=point_m[r],z=point_z[r];
                const double projection=fi>=0 ?
                    std::inner_product(ch[r].begin(),ch[r].end(),response_design[t].begin(),0.)/fit_variance[fi] : 0.;
                for (int n=0;n<3;++n) {
                    double g=fi>=0 ? ch[r][3*a+n]*residual/fit_variance[fi]-
                                         projection*migration_fit.parameters[3*a+n] : 0.;
                    if (t==r) {
                        g-=f/m*migration_fit.parameters[3*a+n];
                        if (v==r) g-=f*z/(m*m)*migration_fit.parameters[3*a+n];
                    }
                    gradient[r][n]=g;
                }
                gradient[r][3]=0.;
                if (v==r) {
                    // The reference epsilon is sum(base_w*eps)/sum(base_w)
                    // over every accepted reconstructed destination of v.
                    // Its numerator shares MC events with all three basis
                    // terms, so retain their event-level cross moments.
                    const auto& p=experimental_points[r];
                    const double phi=.5*(phi_edges[r%cfg.n_phi]+phi_edges[r%cfg.n_phi+1]);
                    const double eps=p.epsilon_reference;
                    const double k=std::sqrt(2.*eps*(1.+eps));
                    const double df_deps=((1.+2.*eps)/k*std::cos(phi)*
                                          migration_fit.parameters[3*a+1]+ 
                                          std::cos(2.*phi)*migration_fit.parameters[3*a+2])/
                                         (2.*TMath::Pi());
                    const double ge=p.correction*df_deps/truth_phi_moments[v].weight;
                    gradient[r][0]-=ge*(2.*TMath::Pi())*eps;
                    gradient[r][3]=ge;
                }
            }
            for (int r=0;r<nr;++r) {
                if (experimental_points[r].status!="ok") continue;
                for (int s=0;s<=r;++s) {
                    if (experimental_points[s].status!="ok") continue;
                    double c=0.;
                    for (int i=0;i<3;++i) for (int j=0;j<3;++j)
                        c+=gradient[r][i]*entry.second.covariance[3*i+j]*gradient[s][j];
                    for (int i=0;i<3;++i)
                        c+=(gradient[r][i]*gradient[s][3]+gradient[r][3]*gradient[s][i])*
                           entry.second.basis_epsilon_covariance[i];
                    c+=gradient[r][3]*gradient[s][3]*entry.second.epsilon_weight_covariance;
                    experimental_point_covariance[static_cast<size_t>(r)*nr+s]+=c;
                    if (r!=s) experimental_point_covariance[static_cast<size_t>(s)*nr+r]+=c;
                }
            }
        }
    }
    for (int r=0;r<nr;++r) if (experimental_points[r].status=="ok") {
        const double v=experimental_point_covariance[static_cast<size_t>(r)*nr+r];
        if (!std::isfinite(v) || v<0.) die("Invalid experimental-point variance in row "+std::to_string(r));
        experimental_points[r].sigma_exp_err=std::sqrt(v);
    }
}

inline void ExclPi0XSecAnalysis::write_experimental_points() {
    fout->cd();
    const int nr=static_cast<int>(experimental_points.size());
    std::ofstream out(fs::path(cfg.out_dir)/"experimental_points.csv");
    out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"reco_row,truth_cell,truth_block,reco_phi_bin,vertex_phi_bin,status,"
          "phi_reference,q2_reference,xb_reference,tprime_reference,epsilon_reference,"
          "data_reco_q2,data_reco_xb,data_reco_tprime,predicted_reco_q2,predicted_reco_xb,predicted_reco_tprime,"
          "diagonal_contribution,subtracted_yield,correction,sigma_fit_reference,sigma_exp,sigma_exp_error,sigma_exp_target_error\n";
    TTree tree("experimental_points","Residual-corrected virtual-photon cross section at truth-cell reference");
    int row=0,cell=0,block=0,ip=0; std::string status;
    double phi=0,q2=0,xb=0,tp=0,eps=0,y=0,corr=0,fit=0,sigma=0,error=0;
    tree.Branch("reco_row",&row);tree.Branch("truth_cell",&cell);tree.Branch("truth_block",&block);
    tree.Branch("phi_bin",&ip);tree.Branch("status",&status);
    tree.Branch("phi_reference",&phi);tree.Branch("q2_reference",&q2);
    tree.Branch("xb_reference",&xb);tree.Branch("tprime_reference",&tp);
    tree.Branch("epsilon_reference",&eps);tree.Branch("subtracted_yield",&y);
    tree.Branch("correction",&corr);tree.Branch("sigma_fit_reference",&fit);
    tree.Branch("sigma_exp",&sigma);tree.Branch("sigma_exp_error",&error);
    for (const auto& p:experimental_points) {
        row=p.row;cell=p.truth_cell;block=cell/cfg.n_phi;ip=cell%cfg.n_phi;status=p.status;
        phi=.5*(phi_edges[ip]+phi_edges[ip+1]);q2=p.q2_reference;xb=p.xb_reference;
        tp=p.tprime_reference;eps=p.epsilon_reference;y=p.subtracted_yield;
        corr=p.correction;fit=p.sigma_reference;sigma=p.sigma_exp;error=p.sigma_exp_err;
        tree.Fill();
        const auto& pb=slices[row/cfg.n_phi].phi[row%cfg.n_phi];
        const double data_q=pb.n_data>0 && pb.data!=0. ? pb.mean_q2_data :
            std::numeric_limits<double>::quiet_NaN();
        const double data_x=pb.n_data>0 && pb.data!=0. ? pb.mean_xb_data :
            std::numeric_limits<double>::quiet_NaN();
        const double data_tp=pb.n_data>0 && pb.data!=0. ? pb.mean_tprime_data :
            std::numeric_limits<double>::quiet_NaN();
        out<<row<<','<<cell<<','<<block<<','<<ip<<','<<ip<<','<<status<<','
           <<phi<<','<<q2<<','<<xb<<','<<tp<<','<<eps<<','
           <<data_q<<','<<data_x<<','<<data_tp<<','
           <<p.reco_prediction_q2<<','<<p.reco_prediction_xb<<','<<p.reco_prediction_tprime<<','
           <<p.contribution<<','<<y<<','<<corr<<','<<fit<<','<<sigma<<','<<error<<','
           <<std::abs(sigma)*cfg.tgt_contam_err/cfg.tgt_contam<<'\n';
    }
    tree.Write();
    TMatrixD covariance(nr,nr),correlation(nr,nr),target_covariance(nr,nr);
    std::ofstream covout(fs::path(cfg.out_dir)/"experimental_point_covariance.csv");
    std::ofstream corrout(fs::path(cfg.out_dir)/"experimental_point_correlation.csv");
    covout<<std::setprecision(std::numeric_limits<double>::max_digits10);
    corrout<<std::setprecision(std::numeric_limits<double>::max_digits10);
    covout<<"reco_row_i,reco_row_j,statistical_covariance,target_correlated_covariance\n";
    corrout<<"reco_row_i,reco_row_j,statistical_correlation\n";
    for(int i=0;i<nr;++i) for(int j=0;j<nr;++j) {
        covariance(i,j)=experimental_point_covariance[static_cast<size_t>(i)*nr+j];
        const double vi=experimental_point_covariance[static_cast<size_t>(i)*nr+i];
        const double vj=experimental_point_covariance[static_cast<size_t>(j)*nr+j];
        correlation(i,j)=(vi>0 && vj>0 && std::isfinite(covariance(i,j))) ?
            covariance(i,j)/std::sqrt(vi*vj) : std::numeric_limits<double>::quiet_NaN();
        target_covariance(i,j)=experimental_points[i].sigma_exp*experimental_points[j].sigma_exp*
            std::pow(cfg.tgt_contam_err/cfg.tgt_contam,2);
        covout<<i<<','<<j<<','<<covariance(i,j)<<','<<target_covariance(i,j)<<'\n';
        corrout<<i<<','<<j<<','<<correlation(i,j)<<'\n';
    }
    covariance.Write("experimental_point_covariance");
    correlation.Write("experimental_point_correlation");
    target_covariance.Write("experimental_point_target_covariance");
    if (positive_toy_point_covariance.size()==static_cast<size_t>(nr)*nr) {
        TMatrixD toy_covariance(nr,nr),toy_correlation(nr,nr);
        for(int i=0;i<nr;++i) for(int j=0;j<nr;++j) {
            const double v=positive_toy_point_covariance[static_cast<size_t>(i)*nr+j];
            const double vi=positive_toy_point_covariance[static_cast<size_t>(i)*nr+i];
            const double vj=positive_toy_point_covariance[static_cast<size_t>(j)*nr+j];
            toy_covariance(i,j)=v;
            toy_correlation(i,j)=vi>0. && vj>0. ? v/std::sqrt(vi*vj) :
                std::numeric_limits<double>::quiet_NaN();
        }
        toy_covariance.Write("experimental_point_positive_refit_toy_covariance");
        toy_correlation.Write("experimental_point_positive_refit_toy_correlation");
    }
    std::ofstream cells(fs::path(cfg.out_dir)/"truth_phi_response_cells.csv");
    cells<<std::setprecision(std::numeric_limits<double>::max_digits10);
    cells<<"reco_row,truth_cell,truth_block,vertex_phi_bin,events,basis_U,basis_LT,basis_TT,cov_U_U,cov_U_LT,cov_U_TT,cov_LT_U,cov_LT_LT,cov_LT_TT,cov_TT_U,cov_TT_LT,cov_TT_TT,reco_q2_U,reco_q2_LT,reco_q2_TT,reco_xb_U,reco_xb_LT,reco_xb_TT,reco_tprime_U,reco_tprime_LT,reco_tprime_TT,epsilon_weight,cov_U_epsilon_weight,cov_LT_epsilon_weight,cov_TT_epsilon_weight,cov_epsilon_weight_epsilon_weight\n";
    for(int r=0;r<nr;++r) for(const auto& entry:truth_phi_response[r]) {
        cells<<r<<','<<entry.first<<','<<entry.first/cfg.n_phi<<','<<entry.first%cfg.n_phi<<','
             <<entry.second.events;
        for(double v:entry.second.basis) cells<<','<<v;
        for(double v:entry.second.covariance) cells<<','<<v;
        for(double v:entry.second.reco_q2) cells<<','<<v;
        for(double v:entry.second.reco_xb) cells<<','<<v;
        for(double v:entry.second.reco_tprime) cells<<','<<v;
        cells<<','<<entry.second.epsilon_weight;
        for(double v:entry.second.basis_epsilon_covariance) cells<<','<<v;
        cells<<','<<entry.second.epsilon_weight_covariance;
        cells<<'\n';
    }
    out.flush();covout.flush();corrout.flush();cells.flush();
    if(!out || !covout || !corrout || !cells) die("Failed writing experimental-point output");
}
