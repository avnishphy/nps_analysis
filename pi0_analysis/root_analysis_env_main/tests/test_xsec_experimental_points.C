// Deterministic response fixture; no production files or generated artifacts.
#include "../src/xsec_extract/xsec_root.h"
#include "../src/xsec_extract/xsec_response.h"
#define private public
#include "../src/xsec_extract/xsec_analysis.h"
#undef private
#define main xsec_production_main
#include "../src/xsec_extract/excl_xsec_pi0_analysis_no_simc_model.C"
#undef main

static void check(bool pass,const char* what) {
    if(!pass) throw std::runtime_error(what);
}
static void close(double a,double b,double tol,const char* what) {
    if(!(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)})))
        throw std::runtime_error(what);
}
static AnalysisConfig config(const char* mode) {
    AnalysisConfig c;c.n_tprime=2;c.n_q2=1;c.n_xb=1;c.n_phi=8;
    c.fit_variance_mode=mode;c.tgt_contam=1.;c.tgt_contam_err=0.;
    c.verbose=false;c.write_pdf=c.write_png=c.diagnostics=false;
    return c;
}
static void fixture(ExclPi0XSecAnalysis& a,bool migration,bool mc=true) {
    const int nphi=8,nrow=16,nblock=8;
    a.slices.resize(2);a.truth_moments.resize(nblock);
    a.migration_response.assign(nrow,std::vector<nps_xsec::ResponseCell>(nblock));
    a.truth_phi_response.resize(nrow);a.truth_phi_moments.resize(nblock*nphi);
    a.phi_edges.resize(nphi+1);a.q2_edges={3.,5.};a.tprime_edges={-.4,-.2,0.};
    a.xb_edges_by_q2={{.25,.45}};
    for(int i=0;i<=nphi;++i) a.phi_edges[i]=2*TMath::Pi()*i/nphi;
    for(auto& s:a.slices) s.phi.resize(nphi);
    for(auto& m:a.truth_moments) {
        m.weight=1.;m.q2=4.;m.xb=.35;m.t=-.3;m.tprime=-.1;m.epsilon=.7;m.epsilon_max=.7;
    }
    for(auto& m:a.truth_phi_moments) {
        m.weight=1.;m.q2=4.;m.xb=.35;m.t=-.3;m.tprime=-.1;m.epsilon=.7;m.epsilon_max=.7;
    }
    const double truth[9]={10.,.6,-.3,12.,-.4,.2,3.,.1,-.2};
    auto event=[&](int row,int block,int vertex_phi,const std::array<double,3>& basis) {
        a.migration_response[row][block].add(basis,4.1,.34,-.18);
        a.truth_phi_response[row][block*nphi+vertex_phi].add(basis,4.1,.34,-.18);
    };
    for(int row=0;row<nrow;++row) {
        const int block=row/nphi,ip=row%nphi,other=1-block;
        const double phi=2*TMath::Pi()*(ip+.5)/nphi;
        event(row,block,ip,{1.,std::cos(phi),std::cos(2*phi)});
        if(migration) {
            const int vp=(ip+1)%nphi;
            const double shifted=2*TMath::Pi()*(vp+.5)/nphi;
            event(row,block,vp,{.25,.25*std::cos(shifted),.25*std::cos(2*shifted)});
            event(row,other,(ip+2)%nphi,{.12,.10*std::sin(phi),.08*std::cos(3*phi)});
            event(row,2,(ip+3)%nphi,{.09*(1+.3*std::sin(3*phi)),.07*std::sin(phi+.2),.06*std::cos(3*phi)});
        }
        double y=0.;
        for(int b=0;b<(migration?3:2);++b)
            for(int j=0;j<3;++j) y+=a.migration_response[row][b].basis[j]*truth[3*b+j];
        a.slices[block].phi[ip].data=y;
        a.slices[block].phi[ip].data_sumw2=1.+.1*row;
    }
    if(!mc) for(auto& row:a.truth_phi_response) for(auto& kv:row)
        kv.second.covariance.fill(0.);
    if(!mc) for(auto& row:a.migration_response) for(auto& cell:row)
        cell.covariance.fill(0.);
}
int main(int argc,char** argv) {
    gROOT->SetBatch(kTRUE);
    ExclPi0XSecAnalysis a(config("data"));
    fixture(a,true);
    a.fit_slices();a.compute_experimental_points();
    check(a.successful_fit_groups==1,"off-diagonal fixture did not fit");
    for(int b=0;b<2;++b) {
        close(a.slices[b].fit_xsec.sigmaU,b?12.:10.,1e-9,"injected U");
        close(a.slices[b].fit_xsec.sigmaTL,b?-.4:.6,1e-9,"injected LT");
    }
    for(int r=0;r<16;++r) {
        const auto& p=a.experimental_points[r];
        check(p.status=="ok","closure point unavailable");
        close(p.correction,1.,1e-10,"exact closure correction");
        close(p.sigma_exp,p.sigma_reference,1e-10,"exact closure point");
    }
    // A changed row must subtract all other truth cells, including guard 2.
    a.slices[0].phi[0].data+=.7;
    a.fit_slices();a.compute_experimental_points();
    const auto& p=a.experimental_points[0];
    const double mu=a.slices[0].phi[0].sim;
    close(p.subtracted_yield,a.slices[0].phi[0].data-(mu-p.contribution),1e-12,"Eq. 5.31 subtraction");
    check(a.truth_phi_response[0].count(2*8+3)>0,"exterior feed-in absent");
    check(a.truth_phi_response[0].count(1)>0,"vertex phi migration absent");
    check(std::abs(p.correction-a.slices[0].phi[0].data/mu)>1e-5,
          "incorrect total-yield correction");

    // Independent centered differences reproduce cross-bin point covariance.
    const auto cov=a.experimental_point_covariance;
    double J[2][16]{};
    for(int t=0;t<16;++t) {
        auto& y=a.slices[t/8].phi[t%8].data; const double original=y;
        y=original+1e-5;a.fit_slices();a.compute_experimental_points();
        const double plus0=a.experimental_points[0].sigma_exp,plus1=a.experimental_points[1].sigma_exp;
        y=original-1e-5;a.fit_slices();a.compute_experimental_points();
        J[0][t]=(plus0-a.experimental_points[0].sigma_exp)/2e-5;
        J[1][t]=(plus1-a.experimental_points[1].sigma_exp)/2e-5;
        y=original;
    }
    double independent=0.;
    for(int t=0;t<16;++t)
        independent+=J[0][t]*J[1][t]*a.slices[t/8].phi[t%8].data_sumw2;
    close(cov[1],independent,1e-5,"cross-bin point covariance");

    ExclPi0XSecAnalysis diagonal(config("data"));fixture(diagonal,false);
    diagonal.fit_slices();diagonal.compute_experimental_points();
    for(int r=0;r<16;++r)
        close(diagonal.experimental_points[r].correction,
              diagonal.slices[r/8].phi[r%8].data/diagonal.slices[r/8].phi[r%8].sim,
              1e-12,"no-migration limit");
    diagonal.truth_phi_response[0].erase(0);
    diagonal.compute_experimental_points();
    check(diagonal.experimental_points[0].status=="no_diagonal_response","missing diagonal status");
    diagonal.truth_phi_response[1].at(1).basis={0.,0.,0.};
    diagonal.compute_experimental_points();
    check(diagonal.experimental_points[1].status=="unstable_denominator","zero denominator status");

    ExclPi0XSecAnalysis d(config("data"));fixture(d,true,false);d.fit_slices();d.compute_experimental_points();
    ExclPi0XSecAnalysis f(config("finite-mc"));fixture(f,true,false);f.fit_slices();f.compute_experimental_points();
    check(f.mc_converged,"finite-MC data-only limit did not converge");
    for(int i=0;i<16;++i)
        for(int j=0;j<16;++j)
            close(f.experimental_point_covariance[i*16+j],d.experimental_point_covariance[i*16+j],
                  1e-8,"finite-MC data-only covariance limit");
    const double weight[]={.3,2.,5.},sigcm[]={.1,1.5,4.};
    for(int i=0;i<3;++i) {
        const double factor=1.+.7*i;
        close((weight[i]*factor)/(sigcm[i]*factor),weight[i]/sigcm[i],1e-14,
              "full_weight/sigcm event-scale invariance");
    }
    ExclPi0XSecAnalysis failed(config("data"));fixture(failed,false);
    for(auto& s:failed.slices) for(auto& pb:s.phi) pb.data_sumw2=0.;
    failed.fit_slices();failed.compute_experimental_points();
    check(failed.successful_fit_groups==0 && failed.experimental_points[0].row==0 &&
          failed.experimental_points[0].status=="unavailable_fit","failed-fit point status");
    const fs::path output=argc>1?argv[1]:"/tmp/nps_xsec_points_failed";
    fs::create_directories(output);failed.cfg.out_dir=output.string();
    failed.fout=TFile::Open((output/"failed.root").c_str(),"RECREATE");
    failed.write_experimental_points();
    check(fs::exists(output/"experimental_points.csv"),"failed-fit point export");
    std::cout<<"PASS experimental points: migration, vertex phi, exterior, closure, covariance, support, MC limit, weight ratio\n";
}
