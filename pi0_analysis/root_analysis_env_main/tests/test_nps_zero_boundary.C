#include "../src/analysis/nps_comb_bg_pepsi.h"
#include <cassert>
void test_nps_zero_boundary() {
    TH1D zero("zero_boundary","",200,0,.4); zero.SetDirectory(nullptr);
    for(int b=1;b<=200;++b) {
        const double x=zero.GetBinCenter(b);
        zero.SetBinContent(b,(x>.11 && x<.15)?10:-.1);
        zero.SetBinError(b,1);
    }
    auto fit=nps::FitCombinatorialBGAndSubtract(&zero,"",1,4,.01,.11,.15,.4,false);
    if(!fit.success || !fit.zero_background || fit.amplitude!=0 || fit.background_integral!=0)
        throw std::runtime_error("Zero background was not accepted exactly");
    for(int b=1;b<=200;++b) if(fit.h_final->GetBinContent(b)!=zero.GetBinContent(b))
        throw std::runtime_error("Zero background changed the input residual");
    delete fit.h_final;
    // A small amplitude that improves chi2 above its numerical precision is
    // not clipped. Numerical classification must be invariant to yield units.
    TGraphErrors graph(20);
    for(int i=0;i<20;++i) { graph.SetPoint(i,.01+.005*i,1e-12); graph.SetPointError(i,0,1e-12); }
    double upper,chi0;
    if(nps::CertifyPepsiZero(graph,.11,.22,.001,.1,upper,chi0))
        throw std::runtime_error("A numerically significant positive amplitude was clipped");
    std::cout<<"PASS certified zero prediction; tiny positive score is not clipped\n";
}
