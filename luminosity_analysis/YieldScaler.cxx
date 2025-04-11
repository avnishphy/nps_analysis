#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void plot_yield_scaler(run, current, correction){

    int run = run;
    int current = current;
    int correction = correction;
    std::string filename = Form("/lustre24/expphy/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%d_1_-1.root", run);
    TFile *data_file = new TFile(filename.c_str(), "READ");

    TFile *scaler_tree = "TSH"; //should be "TSH"

    // defining the branches to be used
    double_t H_BCM4A_scalerCharge, H_BCM2_scalerCharge, H_BCM4A_scalerCurrent, H_BCM2_scalerCurrent,
    H_EDTM_scaler, H_hTRIG4_scaler, H_1MHz_scalerTime;

    double_t charge_accumulated;

    scaler_tree->SetBranchAddress("H.BCM4A.scalerCharge", &H_BCM4A_scalerCharge);
    scaler_tree->SetBranchAddress("H.BCM2.scalerCharge", &H_BCM2_scalerCharge);
    scaler_tree->SetBranchAddress("H.BCM4A.scalerCurrent", &H_BCM4A_scalerCurrent);
    scaler_tree->SetBranchAddress("H.BCM2.scalerCurrent", &H_BCM2_scalerCurrent);
    scaler_tree->SetBranchAddress("H.EDTM.scaler", &H_EDTM_scaler);
    scaler_tree->SetBranchAddress("H.hTRIG4.scaler", &H_hTRIG4_scaler);
    scaler_tree->SetBranchAddess("H.1MHz.scalerTime", H_1MHz_scalerTime);

    for (int i = 0; scaler_tree->GetEntries(); i++){
        scaler_tree->GetEntry(i);
        if (abs(H_BCM4A_scalerCurrent - current)<2.5){
            double_t curr_corr = (H_BCM4A_current + correction)/H_BCM4A_current;
            double_t final_current = current*curr_corr;
            accumulated_charge = accumulated_charge + final_current*(H_1MHz_scalerTime - H_1MHz_scalerTime[i-1]);
        }
        else{
            accumulated_edtm_scaler = H_EDTM_scaler[i-1];
            accumulated_hTRIG4_scaler = H_hTRIG4_scaler[i-1];
            continue;
        }
    }

    yield_charge_normalized = (accumulated_hTRIG4_scaler - accumulated_edtm_scaler)/accumulated_charge;

    return yield_charge_normalized;

}

double_t main(run, current, correction){
    int run = run;
    int current = current;
    int correction = correction;

    //plot yield_charge_normalized vs final_current for the runs in the given run range.
    // the input run range for main function would look like 1523-1530 but 
    // each run in the given range could have a different current so the code should be able to ask as input that current value somehow too 
}