#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TList.h"
#include "TString.h"

void analyze_tree(TString ana_type="") {

    // Declare common variables
    TFile *input_file;
    TTree *tree;
    Long64_t nentries;

    // Declare histograms
    TH1F *H_W, *H_Q2, *H_theta_e, *H_theta_p, *H_Pmx, *H_Pmy, *H_Pmz, *H_h_delta, *H_e_delta;
    TH1F *H_h_xptar, *H_e_xptar, *H_h_yptar, *H_e_yptar, *H_h_ytar, *H_e_ytar, *H_Em;
    TH1F *H_W_data, *H_Q2_data, *H_Pmx_data, *H_Pmy_data, *H_Pmz_data, *H_e_delta_data, *H_h_delta_data, *H_h_xptar_data, *H_e_xptar_data, *H_h_yptar_data, *H_e_yptar_data, *H_Em_data;
    // Create a list to manage histograms
    TList *kin_HList = new TList();

    // Check analysis type and initialize accordingly
    if (ana_type == "simc") {
        std::cout << "Analyzing SIMC data..." << std::endl;

        // SIMC-specific input file
        TString simc_fname = "d2_heep_scan_rad_0.root";
        input_file = new TFile(simc_fname.Data(), "READ");
        tree = (TTree*)input_file->Get("SNT");
        nentries = tree->GetEntries();

        // Define SIMC-specific variables
        Double_t W, Q2, theta_e, theta_p, Pmx, Pmy, Pmz, h_delta, e_delta;
        Double_t h_xptar, e_xptar, h_yptar, e_yptar, h_ytar, e_ytar;
        Double_t Normfac, Weight, FullWeight, Em;

        // Set branch addresses
        tree->SetBranchAddress("W", &W);
        tree->SetBranchAddress("Q2", &Q2);
        tree->SetBranchAddress("theta_e", &theta_e);
        tree->SetBranchAddress("theta_p", &theta_p);
        tree->SetBranchAddress("Pmx", &Pmx);
        tree->SetBranchAddress("Pmy", &Pmy);
        tree->SetBranchAddress("Pmz", &Pmz);
        tree->SetBranchAddress("h_delta", &h_delta);
        tree->SetBranchAddress("e_delta", &e_delta);
        tree->SetBranchAddress("h_xptar", &h_xptar);
        tree->SetBranchAddress("e_xptar", &e_xptar);
        tree->SetBranchAddress("h_yptar", &h_yptar);
        tree->SetBranchAddress("e_yptar", &e_yptar);
        tree->SetBranchAddress("h_ytar", &h_ytar);
        tree->SetBranchAddress("e_ytar", &e_ytar);
        tree->SetBranchAddress("Normfac", &Normfac);
        tree->SetBranchAddress("Weight", &Weight);
        tree->SetBranchAddress("Em", &Em);

        // Define histograms
        H_W = new TH1F("W", "W", 100, 0.7, 1.3);
        H_Q2 = new TH1F("Q2", "Q2", 100, 2.8, 4.5);
        H_theta_e = new TH1F("theta_e", "theta_e", 100, 0.18, 0.23);
        H_theta_p = new TH1F("theta_p", "theta_p", 100, 0.6, 0.75);
        H_Pmx = new TH1F("Pmx", "Pmx", 100, -0.1, 0.1);
        H_Pmy = new TH1F("Pmy", "Pmy", 100, -0.1, 0.1);
        H_Pmz = new TH1F("Pmz", "Pmz", 100, -0.1, 0.2);
        H_h_delta = new TH1F("h_delta", "h_delta", 100, -20, 15);
        H_e_delta = new TH1F("e_delta", "e_delta", 100, -17, 4);
        H_h_xptar = new TH1F("h_xptar", "h_xptar", 100, -0.2, 0.2);
        H_e_xptar = new TH1F("e_xptar", "e_xptar", 100, -0.03, 0.03);
        H_h_yptar = new TH1F("h_yptar", "h_yptar", 100, -0.06, 0.06);
        H_e_yptar = new TH1F("e_yptar", "e_yptar", 100, -0.02, 0.02);
        H_h_ytar = new TH1F("h_ytar", "h_ytar", 100, -15, 15);
        H_e_ytar = new TH1F("e_ytar", "e_ytar", 100, -2, 2);
        H_Em = new TH1F("Em","Em", 100, -0.1, 0.3);
        // Add histograms to the list
        kin_HList->Add(H_W);
        kin_HList->Add(H_Q2);
        kin_HList->Add(H_theta_e);
        kin_HList->Add(H_theta_p);
        kin_HList->Add(H_Pmx);
        kin_HList->Add(H_Pmy);
        kin_HList->Add(H_Pmz);
        kin_HList->Add(H_h_delta);
        kin_HList->Add(H_e_delta);
        kin_HList->Add(H_h_xptar);
        kin_HList->Add(H_e_xptar);
        kin_HList->Add(H_h_yptar);
        kin_HList->Add(H_e_yptar);
        kin_HList->Add(H_h_ytar);
        kin_HList->Add(H_e_ytar);
        kin_HList->Add(H_Em);
        // Event loop for SIMC
        for (Long64_t ientry = 0; ientry < nentries; ientry++) {
            tree->GetEntry(ientry);
            FullWeight = Normfac * Weight / nentries;
            H_W->Fill(W, FullWeight);
            H_Q2->Fill(Q2, FullWeight);
            H_theta_e->Fill(theta_e, FullWeight);
            H_theta_p->Fill(theta_p, FullWeight);
            H_Pmx->Fill(Pmx, FullWeight);
            H_Pmy->Fill(Pmy, FullWeight);
            H_Pmz->Fill(Pmz, FullWeight);
            H_h_delta->Fill(h_delta, FullWeight);
            H_e_delta->Fill(e_delta, FullWeight);
            H_h_xptar->Fill(h_xptar, FullWeight);
            H_e_xptar->Fill(e_xptar, FullWeight);
            H_h_yptar->Fill(h_yptar, FullWeight);
            H_e_yptar->Fill(e_yptar, FullWeight);
            H_h_ytar->Fill(h_ytar, FullWeight);
            H_e_ytar->Fill(e_ytar, FullWeight);
            H_Em->Fill(Em, FullWeight);
            std::cout << "SIMC Events Completed: " << std::setprecision(2)
                      << double(ientry) / nentries * 100. << " % " << std::flush << "\r";
        }

        TString simc_OutputFileName = "d2_heep_0_output.root";
        TFile *outROOT = new TFile(simc_OutputFileName, "RECREATE");
        outROOT->cd();
        kin_HList->Write();
        outROOT->Close();

    } else if (ana_type == "data") {
        std::cout << "Analyzing DATA..." << std::endl;

        // DATA-specific input file
        TString data_fname = "deut_replay_prod_20851_-1.root";
        input_file = new TFile(data_fname.Data(), "READ");
        tree = (TTree*)input_file->Get("T");
        nentries = tree->GetEntries();

        // Define DATA-specific variables
        Double_t W_data, Q2_data, h_delta_data, e_delta_data, Pmx_data, Pmy_data, Pmz_data,h_xptar_data, h_yptar_data, e_xptar_data, e_yptar_data, Em_data ;

        // Set branch addresses
        tree->SetBranchAddress("P.kin.primary.W", &W_data);
        tree->SetBranchAddress("P.kin.primary.Q2", &Q2_data);
        tree->SetBranchAddress("H.gtr.dp", &h_delta_data);
        tree->SetBranchAddress("P.gtr.dp", &e_delta_data);
        tree->SetBranchAddress("H.kin.secondary.Prec_x", &Pmx_data);
        tree->SetBranchAddress("H.kin.secondary.Prec_y", &Pmy_data);
        tree->SetBranchAddress("H.kin.secondary.Prec_z", &Pmz_data);
        tree->SetBranchAddress("H.gtr.th", &h_xptar_data);
        tree->SetBranchAddress("H.gtr.ph", &h_yptar_data);
        tree->SetBranchAddress("P.gtr.th", &e_xptar_data);
        tree->SetBranchAddress("P.gtr.ph", &e_yptar_data);
        tree->SetBranchAddress("H.kin.secondary.emiss", &Em_data);
        // Define histograms
        H_W_data = new TH1F("W_data", "W_data", 100, 0.5, 1.3);
        H_Q2_data = new TH1F("Q2_data", "Q2_data", 100, 2.5, 5.5);
        H_h_delta_data = new TH1F("h_delta_data", "h_delta_data", 100, -15, 15);
        H_e_delta_data = new TH1F("e_delta_data", "e_delta_data", 100, -3, 3);
        H_Pmx_data = new TH1F("Pmx_data", "Pmx_data", 100, -0.1, 0.1);
        H_Pmy_data = new TH1F("Pmy_data", "Pmy_data", 100, -0.1, 0.1);
        H_Pmz_data = new TH1F("Pmz_data", "Pmz_data", 100, -0.1, 0.1);
        H_h_xptar_data = new TH1F("h_xptar_data", "h_xptar_data", 100, -0.2, 0.2);
        H_h_yptar_data = new TH1F("h_yptar_data", "h_yptar_data", 100, -0.06, 0.06);
        H_e_xptar_data = new TH1F("e_xptar_data", "e_xptar_data", 100, -0.03, 0.03);
        H_e_yptar_data = new TH1F("e_yptar_data", "e_yptar_data", 100, -0.02, 0.02);
        H_Em_data = new TH1F("Em_data", "Em_data", 100, -0.1, 0.8);
        // Add histograms to the list
        kin_HList->Add(H_W_data);
        kin_HList->Add(H_Q2_data);
        kin_HList->Add(H_h_delta_data);
        kin_HList->Add(H_e_delta_data);
        kin_HList->Add(H_Pmx_data);
        kin_HList->Add(H_Pmy_data);
        kin_HList->Add(H_Pmz_data);
        kin_HList->Add(H_h_xptar_data);
        kin_HList->Add(H_h_yptar_data);
        kin_HList->Add(H_e_xptar_data);
        kin_HList->Add(H_e_yptar_data);
        kin_HList->Add(H_Em_data);
        // Event loop for DATA
        for (Long64_t ientry = 0; ientry < nentries; ientry++) {
            tree->GetEntry(ientry);
            H_W_data->Fill(W_data);
            H_Q2_data->Fill(Q2_data);
            H_h_delta_data->Fill(h_delta_data);
            H_e_delta_data->Fill(e_delta_data);
            H_Pmx_data->Fill(Pmx_data);
            H_Pmy_data->Fill(Pmy_data);
            H_Pmz_data->Fill(Pmz_data);
            H_e_xptar_data->Fill(e_xptar_data);
            H_e_yptar_data->Fill(e_yptar_data);
            H_h_xptar_data->Fill(h_xptar_data);
            H_h_yptar_data->Fill(h_yptar_data);
            H_Em_data->Fill(Em_data);
//Set Line color for Histograms
            H_W_data->SetLineColor(kRed);
            H_Q2_data->SetLineColor(kRed);
            H_h_delta_data->SetLineColor(kRed);
            H_e_delta_data->SetLineColor(kRed);
            H_Pmx_data->SetLineColor(kRed);
            H_Pmy_data->SetLineColor(kRed);
            H_Pmz_data->SetLineColor(kRed);
            H_e_xptar_data->SetLineColor(kRed);
            H_e_yptar_data->SetLineColor(kRed);
            H_h_xptar_data->SetLineColor(kRed);
            H_h_yptar_data->SetLineColor(kRed);
            H_Em_data->SetLineColor(kRed);

            std::cout << "DATA Events Completed: " << std::setprecision(2)
                      << double(ientry) / nentries * 100. << " % " << std::flush << "\r";
        }

        TString data_OutputFileName = "deut_20851_-1_output.root";
        TFile *outROOT = new TFile(data_OutputFileName, "RECREATE");
        outROOT->cd();
        kin_HList->Write();
        outROOT->Close();

    } else {
        std::cerr << "Invalid analysis type. Use 'data' or 'simc'." << std::endl;
    }

    // Clean up
    delete input_file;
    delete kin_HList;
}
