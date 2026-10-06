// Small C ABI for the Python event-resampling driver. All nonlinear fits and
// Gaussian response solves use the production C++ implementations.
#include "nps_comb_bg_pepsi.h"
#include "../xsec_extract/xsec_linear_solver.h"
#include <TROOT.h>
#include <streambuf>
namespace {
struct NullBuffer : std::streambuf { int overflow(int c) override { return c; } };
struct Quiet {
    NullBuffer buffer;
    std::streambuf *out, *err;
    Quiet():out(std::cout.rdbuf(&buffer)),err(std::cerr.rdbuf(&buffer)) {}
    ~Quiet() { std::cout.rdbuf(out); std::cerr.rdbuf(err); }
};
}
extern "C" int nps_fit_spectrum(const double* y,const double* variance,
        double* final,double* info) {
    try {
        Quiet quiet; gROOT->SetBatch(true);
        TH1D h("bootstrap_mass","",200,0,.4); h.SetDirectory(nullptr); h.Sumw2();
        for(int i=0;i<200;++i) { h.SetBinContent(i+1,y[i]); h.SetBinError(i+1,std::sqrt(variance[i])); }
        auto fit=nps::FitCombinatorialBGAndSubtract(&h,"",-1,4,.01,.11,.15,.4,false);
        info[0]=fit.zero_background; info[1]=fit.minimizer_status; info[2]=fit.covariance_status;
        info[3]=fit.chi2; info[4]=fit.ndf; info[5]=fit.amplitude;
        info[6]=fit.background_integral; info[7]=fit.boundary_score_upper;
        if (!fit.success || !fit.h_final) return 1;
        for(int i=0;i<200;++i) final[i]=fit.h_final->GetBinContent(i+1);
        delete fit.h_final; return 0;
    } catch(...) { return 2; }
}
extern "C" int nps_solve_response(int nr,int np,const double* design,
        const double* y,const double* variance,double* parameters,double* covariance) {
    try {
        std::vector<std::vector<double>> a(nr,std::vector<double>(np));
        for(int r=0;r<nr;++r) for(int p=0;p<np;++p) a[r][p]=design[r*np+p];
        const auto fit=nps_xsec::solve_weighted_response(a,
            std::vector<double>(y,y+nr),std::vector<double>(variance,variance+nr));
        std::copy(fit.parameters.begin(),fit.parameters.end(),parameters);
        std::copy(fit.covariance.begin(),fit.covariance.end(),covariance);
        return 0;
    } catch(const std::exception& e) { std::cerr<<"[bootstrap solve] "<<e.what()<<'\n'; return 1; }
}
