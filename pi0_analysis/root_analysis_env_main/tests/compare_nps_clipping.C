// Controlled comparison: same corrected timing input, fit model and events.
#include "../src/analysis/nps_comb_bg_pepsi.h"
#include <TTree.h>
#include <fstream>
#include <stdexcept>
void compare_nps_clipping(const char* input, const char* output, int run=6418) {
    TFile file(input);
    auto timing=dynamic_cast<TH1D*>(file.Get("h_pi0_coin_bgsub"));
    auto all=dynamic_cast<TH1D*>(file.Get(Form("h_mpi0_all_run%d",run)));
    auto tree=dynamic_cast<TTree*>(file.Get("physics"));
    if (!timing || !all || !tree) throw std::runtime_error("Missing comparison inputs");
    auto clipped=static_cast<TH1D*>(timing->Clone("clipped_timing"));
    clipped->SetDirectory(nullptr);
    for (int i=1;i<=clipped->GetNbinsX();++i)
        clipped->SetBinContent(i,std::max(0.,clipped->GetBinContent(i)));
    auto fit=nps::FitCombinatorialBGAndSubtract(clipped,"",run,4,.01,.11,.15,.4,false);
    if (!fit.success) throw std::runtime_error("Clipped comparison fit failed");
    double mass=0,weight=0,mmiss=0;
    tree->SetBranchAddress("mpi0_all",&mass);
    tree->SetBranchAddress("pi0_weight",&weight);
    tree->SetBranchAddress("mmiss_all",&mmiss);
    double old_y[2]={},new_y[2]={},old_v[2]={},new_v[2]={};
    for (Long64_t i=0;i<tree->GetEntries();++i) {
        tree->GetEntry(i);
        const int bin=all->FindBin(mass);
        const double n=all->GetBinContent(bin);
        const double old_w=(bin>=1 && bin<=all->GetNbinsX() && n>0) ?
            std::max(0.,fit.h_final->GetBinContent(bin)/n):0.;
        for (int selected=0;selected<2;++selected) {
            if (selected && !(mmiss>=.8 && mmiss<=1.1)) continue;
            old_y[selected]+=old_w;new_y[selected]+=weight;
            old_v[selected]+=old_w*old_w;new_v[selected]+=weight*weight;
        }
    }
    std::ofstream out(output);
    out << "run,selection,Y_clipped,Y_signed,delta,relative_delta,frozen_sigma_clipped,frozen_sigma_signed\n" << std::setprecision(17);
    for (int i=0;i<2;++i) out << run << ',' << (i?"mmiss_0.8_1.1_only":"all_accepted")
        << ',' << old_y[i] << ',' << new_y[i] << ',' << new_y[i]-old_y[i]
        << ',' << (new_y[i]-old_y[i])/old_y[i] << ',' << std::sqrt(old_v[i])
        << ',' << std::sqrt(new_v[i]) << '\n';
    tree->ResetBranchAddresses();
    delete fit.h_final;delete clipped;
}
