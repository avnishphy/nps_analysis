#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

double plot_yield_scaler(int run, double current, double correction) {
    std::string filename = Form("/lustre24/expphy/cache/hallc/c-nps/analysis/pass1/replays/skim/nps_hms_skim_%d_1_-1.root", run);
    TFile *data_file = new TFile(filename.c_str(), "READ");

    if (!data_file || data_file->IsZombie()) {
        std::cerr << "Error: Couldn't open file for run " << run << std::endl;
        return -1;
    }

    TTree *scaler_tree = (TTree*) data_file->Get("TSH");
    if (!scaler_tree) {
        std::cerr << "Error: Couldn't find scaler tree 'TSH' in file " << filename << std::endl;
        return -1;
    }

    // Variables to read from tree
    double H_BCM4A_scalerCharge, H_BCM4A_scalerCurrent;
    double H_EDTM_scaler, H_hTRIG4_scaler, H_1MHz_scalerTime;

    scaler_tree->SetBranchAddress("H.BCM4A.scalerCharge", &H_BCM4A_scalerCharge);
    scaler_tree->SetBranchAddress("H.BCM4A.scalerCurrent", &H_BCM4A_scalerCurrent);
    scaler_tree->SetBranchAddress("H.EDTM.scaler", &H_EDTM_scaler);
    scaler_tree->SetBranchAddress("H.hTRIG4.scaler", &H_hTRIG4_scaler);
    scaler_tree->SetBranchAddress("H.1MHz.scalerTime", &H_1MHz_scalerTime);

    Long64_t nentries = scaler_tree->GetEntries();
    if (nentries < 2) {
        std::cerr << "Error: Not enough scaler entries in run " << run << std::endl;
        return -1;
    }

    double accumulated_charge = 0;
    double accumulated_edtm = 0;
    double accumulated_hTRIG4 = 0;
    double prev_time = 0;
    double prev_EDTM = 0;
    double prev_hTRIG4 = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        scaler_tree->GetEntry(i);

        if (i == 0) {
            prev_time = H_1MHz_scalerTime;
            prev_EDTM = H_EDTM_scaler;
            prev_hTRIG4 = H_hTRIG4_scaler;
            continue;
        }

        double delta_time = H_1MHz_scalerTime - prev_time;

        if (fabs(H_BCM4A_scalerCurrent - current) < 2.5) {
            double corr_factor = (H_BCM4A_scalerCurrent + correction) / H_BCM4A_scalerCurrent;
            double corrected_current = current * corr_factor;
            accumulated_charge += corrected_current * delta_time;
        }

        accumulated_edtm += (H_EDTM_scaler - prev_EDTM);
        accumulated_hTRIG4 += (H_hTRIG4_scaler - prev_hTRIG4);

        prev_time = H_1MHz_scalerTime;
        prev_EDTM = H_EDTM_scaler;
        prev_hTRIG4 = H_hTRIG4_scaler;
    }

    data_file->Close();
    delete data_file;

    if (accumulated_charge <= 0) {
        std::cerr << "Warning: Accumulated charge is zero or negative for run " << run << std::endl;
        return -1;
    }

    double yield_norm = (accumulated_hTRIG4 - accumulated_edtm) / accumulated_charge;
    return yield_norm;
}

void analyze_run_range() {
    // Example input: run numbers and corresponding currents
    std::vector<int> runs = {1523, 1524, 1525};
    std::vector<double> currents = {60.0, 65.0, 60.0}; // You may prompt or input manually
    double correction = 0.0; // Or any value you'd like

    for (size_t i = 0; i < runs.size(); ++i) {
        double yield = plot_yield_scaler(runs[i], currents[i], correction);
        if (yield > 0)
            std::cout << "Run " << runs[i] << ": Yield/Charge = " << yield << std::endl;
    }
}
