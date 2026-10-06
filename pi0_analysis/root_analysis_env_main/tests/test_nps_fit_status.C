#include "../src/analysis/nps_comb_bg_pepsi.h"
#include <cassert>
#include <TSystem.h>
void test_nps_fit_status(const char* good_fixture="", int run=6418) {
    TH1D insufficient("insufficient","",2,0,.4);
    insufficient.SetBinContent(1,1); insufficient.SetBinError(1,1);
    auto bad=nps::FitCombinatorialBGAndSubtract(&insufficient,"",1,4,.01,.11,.15,.4,false);
    if(bad.success || bad.h_final) throw std::runtime_error("Invalid fit produced a signal histogram");
    std::cout << "PASS invalid fit returns no subtraction histogram\n";
    if (std::string(good_fixture).empty()) return;
    TFile f(good_fixture);
    auto input=dynamic_cast<TH1D*>(f.Get(Form("run_%d/h_coin_bgsub_input",run)));
    auto expected=dynamic_cast<TH1D*>(f.Get(Form("run_%d/h_bgsub_final",run)));
    if(!input || !expected) throw std::runtime_error("Missing healthy fit fixture");
    auto good=nps::FitCombinatorialBGAndSubtract(input,"",run,4,.01,.11,.15,.4,false);
    if(!good.success || !good.h_final || good.attempts!=1) throw std::runtime_error("Healthy fit changed fit path");
    double maximum_difference=0;
    for (int i=1;i<=expected->GetNbinsX();++i)
        maximum_difference=std::max(maximum_difference,std::abs(expected->GetBinContent(i)-good.h_final->GetBinContent(i)));
    if(maximum_difference>=1e-10) throw std::runtime_error("Healthy fit changed nominal yield");
    delete good.h_final;
    std::cout << "PASS valid fit unchanged; max bin difference=" << maximum_difference << '\n';
}
