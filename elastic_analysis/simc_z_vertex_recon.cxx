#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

void plot_z_vertex() {
    // Open the ROOT file
    TFile *data_file = new TFile("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_6834_0_1_-1.root", "READ");

    TFile *simc_file = new TFile("/u/group/nps/singhav/simc_gfortran/worksim/eep_hydrogen_6828_6841_HMS4042.root", "READ");
    if (!simc_file || simc_file->IsZombie()) {
        std::cerr << "Error: Unable to open file!" << std::endl;
        return;
    }

    // Get the tree from the file
    TTree *data_tree = (TTree*) data_file->Get("T");
    TTree *simc_tree = (TTree*) simc_file->Get("h10");
    if (!simc_tree) {
        std::cerr << "Error: Unable to get tree!" << std::endl;
        simc_file->Close();
        return;
    }

    // Define histogram for z_vertex
    TH1F *hist_simc_z_vertex = new TH1F("hist_simc_z_vertex", "Reconstructed z Vertex; z (cm); Counts", 100, -10, 10);
    TH1F *hist_data_z_vertex = new TH1F("hist_data_z_vertex", "Reconstructed z Vertex; z (cm); Counts", 100, -10, 10);

    // Define cut conditions for Data
    double_t H_kin_W, H_kin_x_bj, H_cal_etottracknorm, H_cer_npeSum, H_gtr_dp, H_gtr_th, H_gtr_ph, H_react_z, H_gtr_y, H_kin_Q2;
    double_t H_hod_goodscinhit, H_kin_scat_ang_deg, H_kin_scat_energy;
    data_tree->SetBranchAddress("H.kin.W", &H_kin_W);
    data_tree->SetBranchAddress("H.kin.x_bj", &H_kin_x_bj);
    data_tree->SetBranchAddress("H.cal.etottracknorm", &H_cal_etottracknorm);
    data_tree->SetBranchAddress("H.hod.goodscinhit", &H_hod_goodscinhit);
    data_tree->SetBranchAddress("H.cer.npeSum", &H_cer_npeSum);
    data_tree->SetBranchAddress("H.gtr.dp", &H_gtr_dp);
    data_tree->SetBranchAddress("H.gtr.th", &H_gtr_th);
    data_tree->SetBranchAddress("H.gtr.ph", &H_gtr_ph);
    data_tree->SetBranchAddress("H.gtr.y", &H_gtr_y);
    data_tree->SetBranchAddress("H.react.z", &H_react_z);
    data_tree->SetBranchAddress("H.kin.scat_ang_deg", &H_kin_scat_ang_deg);
    data_tree->SetBranchAddress("H.kin.Q2", &H_kin_Q2);

    // Fill Data Histogram with cuts
    Long64_t nentries_data = data_tree->GetEntries();
    for (Long64_t i = 0; i < nentries_data; i++) {
        data_tree->GetEntry(i);
        if (H_kin_W > 0.88 && H_kin_W < 1 && 
            H_cal_etottracknorm > 0.9 &&
            H_hod_goodscinhit == 1 && H_cer_npeSum > 0.5 && H_gtr_dp<8.5 && H_gtr_dp>-8.5 && H_gtr_th<0.09 && H_gtr_th>-0.09 && H_gtr_ph<0.055
            && H_gtr_ph>-0.055) {
            hist_data_z_vertex->Fill(H_react_z);
        }
    }


    // Define SIMC branches
    float_t W, Weight, fhsdelta, fhsyptar, fhsxptar, fQ2, fnu, fhsytar; 

    simc_tree->SetBranchAddress("hsdelta", &fhsdelta);  
    simc_tree->SetBranchAddress("hsyptar", &fhsyptar);
    simc_tree->SetBranchAddress("hsxptar", &fhsxptar);
    simc_tree->SetBranchAddress("hsytar", &fhsytar);
    simc_tree->SetBranchAddress("Q2", &fQ2);
    simc_tree->SetBranchAddress("nu", &fnu);
    simc_tree->SetBranchAddress("W", &W);
    simc_tree->SetBranchAddress("Weight", &Weight);
    

    // Define constants
    float x_beam = 0.0;  // Set the correct value
    float theta_scat_deg = 24.86; // Set the correct value
    float theta_scat_rad = theta_scat_deg * (3.141592653589793 / 180.0); // Convert to radians
    float y_mispoint = 0.0; // Set if applicable

    double_t normfac = 0.160082E+08;
    double_t nevents = 100000;
    double_t weight_factor = normfac / nevents;
    double_t full_weight;

    double_t effective_charge_data = 46.9562;

    // Loop over the entries in the tree
    Long64_t nentries = simc_tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        simc_tree->GetEntry(i);
        full_weight = Weight*weight_factor;
        if (W > 0.88 && W < 1 && fhsdelta < 8.5 && fhsdelta > -8.5 && fhsxptar > -0.09 && fhsxptar < 0.09 && fhsyptar > -0.055 && fhsyptar < 0.055) {

            // Compute reconstructed z_vertex
            float z_vertex = (fhsytar + y_mispoint - x_beam * (cos(theta_scat_rad) - fhsyptar * sin(theta_scat_rad))) / 
                             (sin(theta_scat_rad) + fhsyptar * cos(theta_scat_rad));

            hist_simc_z_vertex->Fill(z_vertex, full_weight);
        }
    }

    // Disable default stat box
    gStyle->SetOptStat(0);

    // Draw the histogram
    TCanvas *c1 = new TCanvas("c1", "Reconstructed z Vertex", 800, 600);
    hist_data_z_vertex->Scale(1.0 / effective_charge_data);
    hist_data_z_vertex->SetLineColor(kBlue);
    hist_simc_z_vertex->SetLineColor(kRed);
    hist_simc_z_vertex->Draw("hist");
    hist_data_z_vertex->Draw("hist same");

    c1->Update();
    
    // Save the canvas as an image (optional)
    // c1->SaveAs("z_vertex_histogram.png");

    // Close the file
    // simc_file->Close();
}
