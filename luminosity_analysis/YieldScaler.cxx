#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void plot_yield_scaler(run, scaler_tree){

    int run = run;
    std::string filename = Form("/lustre24/expphy/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%d_1_-1.root", run);
    TFile *data_file = new TFile(filename.c_str(), "READ");

    TFile *scaler_tree = scaler_tree; //should be "TSH"

    // defining the branches to be used
    double_t H_BCM4A_scalerCharge, H_BCM2_scalerCharge, H_BCM4A_scalerCurrent, H_BCM2_scalerCurrent,
    H_EDTM_scaler, H_hTRIG4_scaler, H_1MHz_scalerTime;

    scaler_tree->SetBranchAddress("H.BCM4A.scalerCharge", &H_BCM4A_scalerCharge);
    scaler_tree->SetBranchAddress("H.BCM2.scalerCharge", &H_BCM2_scalerCharge);
    scaler_tree->SetBranchAddress("H.BCM4A.scalerCurrent", &H_BCM4A_scalerCurrent);
    scaler_tree->SetBranchAddress("H.BCM2.scalerCurrent", &H_BCM2_scalerCurrent);
    scaler_tree->SetBranchAddress("H.EDTM.scaler", &H_EDTM_scaler);
    scaler_tree->SetBranchAddress("H.hTRIG4.scaler", &H_hTRIG4_scaler);
    scaler_tree->SetBranchAddess("H.1MHz.scalerTime", H_1MHz_scalerTime);

    

}