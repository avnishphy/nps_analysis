// Exercise the production accumulator directly with a controlled signed event.
#include "../src/xsec_extract/xsec_physics.h"
#include "../src/xsec_extract/xsec_response.h"
#include "../src/xsec_extract/xsec_mass_cut.h"
#include "../src/xsec_extract/xsec_vertex_epsilon.h"
#define private public
#include "../src/xsec_extract/xsec_accumulation.h"
#include "../src/xsec_extract/xsec_binning.h"
#undef private
#include "../src/xsec_extract/xsec_plot_global.h"
#include "../src/xsec_extract/xsec_output.h"
#include <cassert>

int main() {
    gROOT->SetBatch(true);
    AnalysisConfig cfg;
    cfg.mmiss_select="window";
    cfg.fit_objective="gaussian";
    ExclPi0XSecAnalysis gaussian(cfg);
    gaussian.build_binning(); gaussian.init_storage();
    auto fill=[](ExclPi0XSecAnalysis& a) {
        a.fill_data_event(4.,-.4,-.1,.36,.5,-.25,2.f,100.,100.,.95,.135,1,0,false,2.6,11);
    };
    fill(gaussian);
    double sum=0,sum2=0;
    for (const auto& s:gaussian.slices) for (const auto& p:s.phi) {
        sum+=p.data;sum2+=p.data_sumw2;
    }
    assert(sum==-.5 && sum2==.25);
    cfg.fit_objective="scaled-poisson";
    ExclPi0XSecAnalysis poisson(cfg);
    poisson.build_binning(); poisson.init_storage();
    bool rejected=false;
    try {fill(poisson);} catch(const std::runtime_error& e) {
        rejected=std::string(e.what()).find("require --fit-objective gaussian")!=std::string::npos;
    }
    assert(rejected);
    std::cout << "PASS Gaussian signed accumulation and explicit scaled-Poisson rejection\n";
}
