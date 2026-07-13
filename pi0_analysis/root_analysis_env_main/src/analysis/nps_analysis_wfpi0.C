// ============================================================================
// File: nps_analysis_wfpi0.C
// Purpose: Deprecated compatibility wrapper for unified nps_analysis_main.C
// ============================================================================

#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <iostream>

void nps_analysis_wfpi0(const TString &skimDir_in = "/lustre24/expphy/volatile/hallc/nps/hhuang/farmFile/Production/DVCS",
                        const TString &outBase_in = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output",
                        const TString &runlistFile = "",
                        const double Ebeam = -1.0)
{
    std::cerr << "[WARN] nps_analysis_wfpi0.C is deprecated. Forwarding to nps_analysis_main.C in waveform mode.\n";

    if (!skimDir_in.IsNull()) gSystem->Setenv("NPS_INPUT_DIR", skimDir_in.Data());
    if (!outBase_in.IsNull()) gSystem->Setenv("NPS_OUTPUT_BASE", outBase_in.Data());
    gSystem->Setenv("NPS_MODE", "waveform");

    if (Ebeam > 0.0) {
        gSystem->Setenv("NPS_EBEAM", Form("%.9f", Ebeam));
    }
    if (!runlistFile.IsNull() && runlistFile.Length() > 0) {
        std::cerr << "[WARN] runlist argument is ignored in unified workflow; run selection comes from config CSV + NPS_KIN.\n";
    }

    TString macro_path = (gSystem->AccessPathName("src/analysis/nps_analysis_main.C") == 0)
        ? "src/analysis/nps_analysis_main.C"
        : "nps_analysis_main.C";

    if (gROOT->ProcessLine(Form(".L %s", macro_path.Data())) != 0) {
        std::cerr << "[ERROR] Failed to load nps_analysis_main.C\n";
        return;
    }
    gROOT->ProcessLine("nps_analysis_main();");
}
