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
void test_nps_timing_geometry() {
    for (bool shifted : {false,true}) {
        TH2D h("timing","",shifted?220:200,shifted?139:140,shifted?161:160,
              shifted?220:200,shifted?139:140,shifted?161:160);
        auto windows=shifted?std::vector<std::pair<double,double>>{
            {139,141},{141,143},{143,145},{155,157},{157,159},{159,161}}:
            nps::default_diag_windows();
        auto low=shifted?std::make_pair(139.,145.):std::make_pair(141.,147.);
        auto high=shifted?std::make_pair(155.,161.):std::make_pair(153.,159.);
        auto bg=nps::estimate_coincidence_background_default(&h,{149,151},
            windows,windows,high,low,high,low);
        assert(bg.area_diag==24 && bg.area_full1==36);
        assert(std::abs(h.GetXaxis()->GetBinWidth(1)-.1)<1e-12);
    }
    TH2D old("old","",200,140,160,200,140,160);
    bool rejected=false;
    try { nps::integral_and_area_TH2(&old,139,141,139,141); }
    catch (const std::runtime_error&) { rejected=true; }
    assert(rejected);
    std::cout << "PASS timing geometry in both modes and out-of-range rejection\n";
}
