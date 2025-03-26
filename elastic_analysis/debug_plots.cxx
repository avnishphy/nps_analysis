#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"

void plot_diagnostic_2D() {
    // Load Data and SIMC File
    TFile *simc_file = new TFile("/u/group/nps/singhav/simc_gfortran/worksim/eep_hydrogen_6828_6841_HMS4042.root", "READ");
    int run = 6834;
    std::string filename = Form("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_%d_0_1_-1.root", run);
    TFile *data_file = new TFile(filename.c_str(), "READ");

    // TChain *data_tree = new TChain("T");
    // for (int i = 6828; i <= 6840; i++) {
    //     data_tree->Add(Form("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_%d_0_1_-1.root", i));
    // }
    TTree *simc_tree = (TTree*) simc_file->Get("h10");
    TTree *data_tree = (TTree*) data_file->Get("T");

    if (!data_tree || !simc_tree) {
        std::cerr << "Error: Unable to get trees!" << std::endl;
        return;
    }

    // Weighting for simc
    double_t normfac = 0.159733E+08;
    double_t nevents = 100000;
    double_t weight_factor = normfac / nevents;
    double_t full_weight;

    // Define branches for Data
    double H_dc_x_fp, H_dc_y_fp, H_dc_xp_fp, H_dc_yp_fp, H_gtr_dp;
    data_tree->SetBranchAddress("H.dc.x_fp", &H_dc_x_fp);
    data_tree->SetBranchAddress("H.dc.y_fp", &H_dc_y_fp);
    data_tree->SetBranchAddress("H.dc.xp_fp", &H_dc_xp_fp);
    data_tree->SetBranchAddress("H.dc.yp_fp", &H_dc_yp_fp);
    data_tree->SetBranchAddress("H.gtr.dp", &H_gtr_dp);

    // Define branches for SIMC
    float fhsxfp, fhsyfp, fhsxpfp, fhsypfp, fhsdelta;
    float Weight;
    simc_tree->SetBranchAddress("hsxfp", &fhsxfp);
    simc_tree->SetBranchAddress("hsyfp", &fhsyfp);
    simc_tree->SetBranchAddress("hsxpfp", &fhsxpfp); 
    simc_tree->SetBranchAddress("hsypfp", &fhsypfp); 
    simc_tree->SetBranchAddress("hsdelta", &fhsdelta);
    simc_tree->SetBranchAddress("Weight", &Weight);
    

    // Define histograms
    TH2F *hist_data_xfp_yfp = new TH2F("data_xfp_yfp", "Data x_fp vs y_fp; x_fp; y_fp", 100, -60, 30, 100, -30, 30);
    TH2F *hist_simc_xfp_yfp = new TH2F("simc_xfp_yfp", "SIMC x_fp vs y_fp; x_fp; y_fp", 100, -60, 30, 100, -30, 30);
    TH2F *hist_data_xpfp_ypfp = new TH2F("data_xpfp_ypfp", "Data x_fp vs y_fp; x_fp; y_fp", 100, -0.07, 0.04, 100, -0.04, 0.04);
    TH2F *hist_simc_xpfp_ypfp = new TH2F("simc_xpfp_ypfp", "SIMC x_fp vs y_fp; x_fp; y_fp", 100, -0.07, 0.04, 100, -0.04, 0.04);
    TH2F *hist_xfp_corr = new TH2F("xfp_corr", "Data x_fp vs SIMC x_fp; SIMC x_fp; Data x_fp", 100, -60, 60, 100, -60, 60);
    TH2F *hist_yfp_corr = new TH2F("yfp_corr", "Data y_fp vs SIMC y_fp; SIMC y_fp; Data y_fp", 100, -30, 30, 100, -30, 30);
    TH2F *hist_delta_corr = new TH2F("delta_corr", "Data delta vs SIMC delta; SIMC delta; Data delta", 100, -15, 5, 100, -15, 5);

    // Loop over Data Entries (No Cuts)
    Long64_t nentries_data = data_tree->GetEntries();
    for (Long64_t i = 0; i < nentries_data; i++) {
        data_tree->GetEntry(i);
        hist_data_xfp_yfp->Fill(H_dc_x_fp, H_dc_y_fp);
        hist_data_xpfp_ypfp->Fill(H_dc_xp_fp, H_dc_yp_fp);
    }

    // Loop over SIMC Entries (No Cuts, Apply Weights)
    Long64_t nentries_simc = simc_tree->GetEntries();
    for (Long64_t i = 0; i < nentries_simc; i++) {
        simc_tree->GetEntry(i);
        full_weight = Weight*weight_factor;
        hist_simc_xfp_yfp->Fill(fhsxfp, fhsyfp, full_weight);
        hist_simc_xpfp_ypfp->Fill(fhsxpfp, fhsypfp, full_weight);
    }

    // Define histograms with cuts
    TH2F *hist_data_xfp_yfp_cut = new TH2F("data_xfp_yfp_cut", "Data x_fp vs y_fp (Cut); x_fp; y_fp", 100, -60, 30, 100, -30, 30);
    TH2F *hist_simc_xfp_yfp_cut = new TH2F("simc_xfp_yfp_cut", "SIMC x_fp vs y_fp (Cut); x_fp; y_fp", 100, -60, 30, 100, -30, 30);
    TH2F *hist_data_xpfp_ypfp_cut = new TH2F("data_xpfp_ypfp_cut", "Data x_fp vs y_fp (Cut); x_fp; y_fp", 100, -0.07, 0.04, 100, -0.04, 0.04);
    TH2F *hist_simc_xpfp_ypfp_cut = new TH2F("simc_xpfp_ypfp_cut", "SIMC x_fp vs y_fp (Cut); x_fp; y_fp", 100, -0.07, 0.04, 100, -0.04, 0.04);
    TH2F *hist_xfp_corr_cut = new TH2F("xfp_corr_cut", "Data x_fp vs SIMC x_fp (Cut); SIMC x_fp; Data x_fp", 100, -60, 60, 100, -60, 60);
    TH2F *hist_yfp_corr_cut = new TH2F("yfp_corr_cut", "Data y_fp vs SIMC y_fp (Cut); SIMC y_fp; Data y_fp", 100, -30, 30, 100, -30, 30);
    TH2F *hist_delta_corr_cut = new TH2F("delta_corr_cut", "Data delta vs SIMC delta (Cut); SIMC delta; Data delta", 100, -15, 5, 100, -15, 5);

    // Apply cuts and fill histograms
    for (Long64_t i = 0; i < nentries_data; i++) {
        data_tree->GetEntry(i);
        if (H_gtr_dp > -8.5 && H_gtr_dp < 8.5) {
            hist_data_xfp_yfp_cut->Fill(H_dc_x_fp, H_dc_y_fp);
            hist_data_xpfp_ypfp_cut->Fill(H_dc_xp_fp, H_dc_yp_fp);
        }
    }

    for (Long64_t i = 0; i < nentries_simc; i++) {
        simc_tree->GetEntry(i);
        if (fhsdelta > -8.5 && fhsdelta < 8.5) {
            full_weight = Weight*weight_factor;
            hist_simc_xfp_yfp_cut->Fill(fhsxfp, fhsyfp, full_weight);
            hist_simc_xpfp_ypfp_cut->Fill(fhsxpfp, fhsypfp, full_weight);
        }
    }

    // Correlation plots
    for (Long64_t i = 0; i < nentries_data; i++) {
        data_tree->GetEntry(i);
        for (Long64_t j = 0; j < nentries_simc; j++) {
            simc_tree->GetEntry(j);
            hist_xfp_corr->Fill(fhsxfp, H_dc_x_fp);
            hist_yfp_corr->Fill(fhsyfp, H_dc_y_fp);
            hist_delta_corr->Fill(fhsdelta, H_gtr_dp);

            if (H_gtr_dp > -8.5 && H_gtr_dp < 8.5 && fhsdelta > -8.5 && fhsdelta < 8.5) {
                hist_xfp_corr_cut->Fill(fhsxfp, H_dc_x_fp);
                hist_yfp_corr_cut->Fill(fhsyfp, H_dc_y_fp);
                hist_delta_corr_cut->Fill(fhsdelta, H_gtr_dp);
            }
        }
    }

    // Save plots
    TCanvas *c1 = new TCanvas("c1", "Diagnostics", 1200, 800);
    c1->Divide(2, 4);
    c1->cd(1); hist_data_xfp_yfp->Draw("COLZ");
    c1->cd(2); hist_simc_xfp_yfp->Draw("COLZ");
    c1->cd(3); hist_data_xpfp_ypfp->Draw("COLZ");
    c1->cd(4); hist_simc_xpfp_ypfp->Draw("COLZ");
    c1->cd(5); hist_xfp_corr->Draw("COLZ");
    c1->cd(6); hist_yfp_corr->Draw("COLZ");
    c1->cd(7); hist_delta_corr->Draw("COLZ");
    c1->Update();
    c1->SaveAs("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/elastic_plots/diagnostic_plots/diagnostics.png");

    TCanvas *c2 = new TCanvas("c2", "Diagnostics with Cuts", 1200, 800);
    c2->Divide(2, 4);
    c2->cd(1); hist_data_xfp_yfp_cut->Draw("COLZ");
    c2->cd(2); hist_simc_xfp_yfp_cut->Draw("COLZ");
    c2->cd(3); hist_data_xpfp_ypfp_cut->Draw("COLZ");
    c2->cd(4); hist_simc_xpfp_ypfp_cut->Draw("COLZ");
    c2->cd(5); hist_xfp_corr_cut->Draw("COLZ");
    c2->cd(6); hist_yfp_corr_cut->Draw("COLZ");
    c2->cd(7); hist_delta_corr_cut->Draw("COLZ");
    c2->Update();
    c2->SaveAs("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/elastic_plots/diagnostic_plots/diagnostics_cut.png");
}
