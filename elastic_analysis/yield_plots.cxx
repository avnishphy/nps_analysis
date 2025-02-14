#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void plot_W_distribution() {
    // Load Data and SIMC Files
    TFile *data_file = new TFile("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_6834_0_1_-1.root", "READ");
    TFile *simc_file = new TFile("/u/group/nps/singhav/simc_gfortran/worksim/eep_hydrogen_6828_6841_HMS4042.root", "READ");
    TFile *dummy_file = new TFile("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_6705_0_1_-1.root", "READ");

    if (!data_file || data_file->IsZombie() || !simc_file || simc_file->IsZombie()) {
        std::cerr << "Error: Unable to open files!" << std::endl;
        return;
    }

    TTree *data_tree = (TTree*) data_file->Get("T");
    TTree *simc_tree = (TTree*) simc_file->Get("h10");
    TTree *dummy_tree = (TTree*) dummy_file->Get("T");

    if (!data_tree || !simc_tree) {
        std::cerr << "Error: Unable to get trees!" << std::endl;
        return;
    }

    double_t m_proton = 0.938; //GeV

    // Define histograms
    double_t nbins = 100;
    double_t xmin = 0.6, xmax = 1.4;

    TH1F *hist_data = new TH1F("hist_data", "W Distribution;W [GeV];Counts", nbins, xmin, xmax);
    TH1F *hist_simc = new TH1F("hist_simc", "W Distribution;W [GeV];Counts", nbins, xmin, xmax);
    TH1F *hist_simc_weighted = new TH1F("hist_simc_weighted", "W Distribution (Weighted);W [GeV];Counts", nbins, xmin, xmax);
    TH1F *hist_exp = new TH1F("hist_exp", "W Distribution;W [GeV];Counts", nbins, xmin, xmax);
    TH1F *hist_dummy = new TH1F("hist_dummy", "W Distribution;W [GeV];Counts", nbins, xmin, xmax);

    // Define cut conditions for Data
    double_t H_kin_W, H_kin_x_bj, H_cal_etottracknorm, H_cer_npeSum, H_gtr_dp, H_gtr_th, H_gtr_ph;
    double_t H_hod_goodscinhit;
    data_tree->SetBranchAddress("H.kin.W", &H_kin_W);
    data_tree->SetBranchAddress("H.kin.x_bj", &H_kin_x_bj);
    data_tree->SetBranchAddress("H.cal.etottracknorm", &H_cal_etottracknorm);
    data_tree->SetBranchAddress("H.hod.goodscinhit", &H_hod_goodscinhit);
    data_tree->SetBranchAddress("H.cer.npeSum", &H_cer_npeSum);
    data_tree->SetBranchAddress("H.gtr.dp", &H_gtr_dp);
    data_tree->SetBranchAddress("H.gtr.th", &H_gtr_th);
    data_tree->SetBranchAddress("H.gtr.ph", &H_gtr_ph);

    // Fill Data Histogram with cuts
    Long64_t nentries_data = data_tree->GetEntries();
    for (Long64_t i = 0; i < nentries_data; i++) {
        data_tree->GetEntry(i);
        if (H_kin_x_bj > 0.96 && H_kin_x_bj < 1.04 && 
            H_cal_etottracknorm > 0.96 && H_cal_etottracknorm < 1.04 &&
            H_hod_goodscinhit == 1 && H_cer_npeSum > 0.5 && H_gtr_dp<8.5 && H_gtr_dp>-8.5 && H_gtr_th<0.09 && H_gtr_th>-0.09 && H_gtr_ph<0.055
            && H_gtr_ph>-0.055) {
            hist_data->Fill(H_kin_W);
        }
    }

    // Define cut conditions for Data
    double_t dH_kin_W, dH_kin_x_bj, dH_cal_etottracknorm, dH_cer_npeSum, dH_gtr_dp, dH_gtr_th, dH_gtr_ph;
    double_t dH_hod_goodscinhit;
    dummy_tree->SetBranchAddress("H.kin.W", &dH_kin_W);
    dummy_tree->SetBranchAddress("H.kin.x_bj", &dH_kin_x_bj);
    dummy_tree->SetBranchAddress("H.cal.etottracknorm", &dH_cal_etottracknorm);
    dummy_tree->SetBranchAddress("H.hod.goodscinhit", &dH_hod_goodscinhit);
    dummy_tree->SetBranchAddress("H.cer.npeSum", &dH_cer_npeSum);
    dummy_tree->SetBranchAddress("H.gtr.dp", &dH_gtr_dp);
    dummy_tree->SetBranchAddress("H.gtr.th", &dH_gtr_th);
    dummy_tree->SetBranchAddress("H.gtr.ph", &dH_gtr_ph);

    // Fill Dummy Histogram with cuts
    Long64_t nentries_dummy = dummy_tree->GetEntries();
    for (Long64_t i = 0; i < nentries_dummy; i++) {
        dummy_tree->GetEntry(i);
        if (dH_kin_x_bj > 0.96 && dH_kin_x_bj < 1.04 && 
            dH_cal_etottracknorm > 0.96 && dH_cal_etottracknorm < 1.04 &&
            dH_hod_goodscinhit == 1 && dH_cer_npeSum > 0.5 && dH_gtr_dp<8.5 && dH_gtr_dp>-8.5 && dH_gtr_th<0.09 && dH_gtr_th>-0.09 && dH_gtr_ph<0.055
            && dH_gtr_ph>-0.055) {
            hist_dummy->Fill(dH_kin_W);
        }
    }

    // Define SIMC branches
    float_t W, Weight;
    simc_tree->SetBranchAddress("W", &W);
    simc_tree->SetBranchAddress("Weight", &Weight);

    float_t fhsdelta, fhsyptar, fhsxptar, fQ2, fnu; 

    simc_tree->SetBranchAddress("hsdelta", &fhsdelta);
    simc_tree->SetBranchAddress("hsyptar", &fhsyptar);
    simc_tree->SetBranchAddress("hsxptar", &fhsxptar);
    simc_tree->SetBranchAddress("Q2", &fQ2);
    simc_tree->SetBranchAddress("nu", &fnu);


    // Fill SIMC Histogram (Unweighted)
    Long64_t nentries_simc = simc_tree->GetEntries();
    // for (Long64_t i = 0; i < nentries_simc; i++) {
    //     simc_tree->GetEntry(i);
    //     if (W>0.96 && W<1.04 && fhsdelta<8.5 && fhsdelta>-8.5 && fhsxptar>-0.09 && fhsxptar<0.09 && fhsyptar>-0.055 && fhsyptar<0.055 ){
    //         hist_simc->Fill(W);
    //     }
    // }

    // Apply Weight and Normalize SIMC Histogram
    // double_t normfac = 0.698679E+07; HMS 3910.20
    // double_t normfac = 0.692370E+07; //HMS 4042
    // double_t normfac = 0.685014E+07; //HMS 3910 20% extra e arm limits
    // double_t normfac = 0.677792E+07; //HMS 4042 10% extra e arm limits
    // double_t normfac =  0.670731E+07; //HMS 4042 20% extra e arm limits
    // double_t normfac = 0.670773E+07; // HMS 4042 run 6706
    double_t normfac = 0.671016E+07; // HMS 4042 run 6834 Ebeam 6.3967
    double_t nevents = 100000;
    double_t weight_factor = normfac / nevents;
    double_t full_weight;
    double_t effective_charge_data;
    double_t effective_charge_dummy;
    double_t scale_factor_dummy = 0.864;
    double_t mass_p = 0.938272; //mass of proton in GeV
    double_t f_xbj; 

    for (Long64_t i = 0; i < nentries_simc; i++) {
        simc_tree->GetEntry(i);
        full_weight = Weight*weight_factor;
        f_xbj = (fQ2)/(2*mass_p*fnu);
        if (f_xbj > 0.96 && f_xbj < 1.04 && fhsdelta < 8.5 && fhsdelta > -8.5 && fhsxptar > -0.09 && fhsxptar < 0.09 && fhsyptar > -0.055 && fhsyptar < 0.055) {
            hist_simc_weighted->Fill(W, full_weight);
        }
        
    }

    // Normalize histograms
    // if (hist_data->Integral() > 0) hist_data->Scale(1.0 / hist_data->Integral());
    // effective_charge_data = 46266.080/1000.0; //BCM4A charge in mC run 6834
    effective_charge_data = 49676.807/1000; //run 6706
    effective_charge_dummy = 22103.178/1000.0;
    hist_data->Scale(1.0 / effective_charge_data);
    hist_dummy->Scale(scale_factor_dummy / effective_charge_dummy);
    hist_exp->Add(hist_data, hist_dummy, 1.0, -1.0); // hist_exp = hist_data - hist_dummy

    // // Get mean and standard deviation
    // double mean_data = hist_data->GetMean();
    // double stddev_data = hist_data->GetStdDev();
    // double mean_simc = hist_simc_weighted->GetMean();
    // double stddev_simc = hist_simc_weighted->GetStdDev();
    // double mean_exp = hist_exp->GetMean();
    // double stddev_exp = hist_exp->GetStdDev();

    // Fit histograms with a Gaussian
    TF1 *fit_data = new TF1("fit_data", "gaus", 0.88, 0.96);
    TF1 *fit_simc_weighted = new TF1("fit_simc_weighted", "gaus", 0.88, 0.96);
    TF1 *fit_exp = new TF1("fit_exp", "gaus", 0.88, 0.96);

    hist_data->Fit(fit_data, "RQ");
    hist_simc_weighted->Fit(fit_simc_weighted, "RQ");
    hist_exp->Fit(fit_exp, "RQ");

    // Get fit parameters
    double mean_data = fit_data->GetParameter(1);
    double sigma_data = fit_data->GetParameter(2);
    double mean_simc = fit_simc_weighted->GetParameter(1);
    double sigma_simc = fit_simc_weighted->GetParameter(2);
    double mean_exp = fit_exp->GetParameter(1);
    double sigma_exp = fit_exp->GetParameter(2);

    // Disable default stat box
    gStyle->SetOptStat(0);

    // Create Canvas
    TCanvas *canvas = new TCanvas("canvas", "W Distribution", 800, 600);

    // Draw Histograms
    hist_data->SetLineColor(kBlue);
    hist_simc_weighted->SetLineColor(kRed);
    hist_exp->SetLineColor(kBlack);
    hist_dummy->SetLineColor(kGreen);
    hist_exp->SetLineWidth(2);

    hist_data->Draw("hist");
    hist_simc_weighted->Draw("hist same");
    hist_exp->Draw("hist same");
    hist_dummy->Draw("hist same");

    double_t yield_exp = hist_exp->Integral();
    double_t yield_simc_weighted = hist_simc_weighted->Integral();
    double_t yield_ratio = (yield_simc_weighted != 0) ? (yield_exp / yield_simc_weighted) : 0.0; // Prevent division by zero
    
    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->AddEntry(hist_data, Form("Data: Mean = %.5f, #sigma = %.5f", mean_data, sigma_data), "l");
    leg->AddEntry(hist_simc_weighted, Form("SIMC: Mean = %.5f, #sigma = %.5f", mean_simc, sigma_simc), "l");
    leg->AddEntry(hist_exp, Form("Exp: Mean = %.5f, #sigma = %.5f", mean_exp, sigma_exp), "l");
    leg->AddEntry((TObject*)0, Form("Yield Ratio = Yield_exp/Yield_simc: %.5f", yield_ratio), ""); // Display yield ratio safely
    leg->Draw();
    

    // Show plot interactively
    canvas->Update();

}
