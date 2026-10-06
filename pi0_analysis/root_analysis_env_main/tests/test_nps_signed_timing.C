#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TLine.h>
#include <TBox.h>
#include <cassert>
#include "../src/analysis/nps_helper.h"
#include "../src/analysis/nps_time_bg.h"
void test_nps_signed_timing() {
    TH2D timing("timing","",200,140,160,200,140,160);
    TH1D coin("coin","",200,0,.4), diag("diag","",200,0,.4);
    coin.Sumw2(); diag.Sumw2();
    coin.Fill(.135);
    for(int i=0;i<12;++i) { diag.Fill(.135); timing.Fill(142,142); }
    timing.Fill(150,150);
    auto bg=nps::estimate_coincidence_background_default(&timing);
    std::vector<TH1D*> diagonals(6,nullptr); diagonals[0]=&diag;
    auto result=nps::make_and_subtract_accidentals_data_driven(&coin,bg,
        nullptr,nullptr,{}, {},diagonals,nullptr,&timing,
        nps::default_diag_windows(),nps::default_side_windows(),"/tmp",999);
    if(!result || std::abs(result->GetBinContent(result->FindBin(.135))+1.)>=1e-12)
        throw std::runtime_error("Signed timing residual failed");
    delete result;
    std::cout << "PASS signed timing subtraction: C=1 T=2 residual=-1\n";
    timing.Reset();
    timing.Fill(142,154);
    timing.Fill(154,142); timing.Fill(154,142);
    auto asymmetric=nps::estimate_coincidence_background_default(&timing);
    if(asymmetric.n_full1_raw!=1 || asymmetric.n_full2_raw!=2)
        throw std::runtime_error("Full accidental boxes use the same quadrant");
    std::cout << "PASS full accidental boxes count opposite timing quadrants\n";
}
