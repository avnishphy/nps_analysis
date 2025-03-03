#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void plot_W_distribution() {
    // Load Data and SIMC File
    TFile *data_file = new TFile("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_6705_0_1_-1_24_85.root", "READ");
    TFile *simc_file = new TFile("/u/group/nps/singhav/simc_gfortran/worksim/eep_hydrogen_6828_6841_HMS4042.root", "READ");
    TFile *dummy_file = new TFile("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_6705_0_1_-1_24_85.root", "READ");

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

    // Define histograms for W
    double_t nbins = 100;
    double_t xmin = 0.6, xmax = 1.4;

    TH1F *hist_data_w = new TH1F("hist_data_w", "W Distribution;W [GeV];Counts", nbins, xmin, xmax);
    TH1F *hist_simc_w = new TH1F("hist_simc_w", "W Distribution;W [GeV];Counts", nbins, xmin, xmax);
    TH1F *hist_simc_weighted_w = new TH1F("hist_simc_weighted_w", "W Distribution (Weighted);W [GeV];Counts", nbins, xmin, xmax);
    TH1F *hist_exp_w = new TH1F("hist_exp_w", "W Distribution;W [GeV];Counts", nbins, xmin, xmax);
    TH1F *hist_dummy_w = new TH1F("hist_dummy_w", "W Distribution;W [GeV];Counts", nbins, xmin, xmax);
    TH1F* hist_dummy_w_up = new TH1F("hist_dummy_w_up", "W distribution (dummy, ytar > 0)", nbins, xmin, xmax);
    TH1F* hist_dummy_w_down = new TH1F("hist_dummy_w_down", "W distribution (dummy, ytar < 0)", nbins, xmin, xmax);
    //Define histograms for delta, xptar, yptar, ytar
    TH1F *hist_data_delta = new TH1F("hist_data_delta", "Delta Distribution;W [GeV];Counts", nbins, -20, 10);
    TH1F *hist_data_xptar = new TH1F("hist_data_xptar", "xptar Distribution;mili radians;Counts", nbins, -0.1, 0.1);
    TH1F *hist_data_yptar = new TH1F("hist_data_yptar", "yptar Distribution;miliradians;Counts", nbins, -0.05, 0.05);
    TH1F *hist_data_ytar = new TH1F("hist_data_ytar", "ytar Distribution;ytar;Counts", nbins, -7, 7);
    TH1F *hist_data_q2 = new TH1F("hist_data_q2", "q2 Distribution;Q2 [GeV^2];Counts", nbins, -7, 7);
    TH1F *hist_data_scat_ang = new TH1F("hist_data_scat_angle", "Scattering angle Distribution;angle [degrees];Counts", nbins, 20, 30);

    TH1F *hist_dummy_delta = new TH1F("hist_dummy_delta", "Delta Distribution;W [GeV];Counts", nbins, -20, 10);
    TH1F *hist_dummy_xptar = new TH1F("hist_dummy_xptar", "xptar Distribution;mili radians;Counts", nbins, -0.1, 0.1);
    TH1F *hist_dummy_yptar = new TH1F("hist_dummy_yptar", "yptar Distribution;miliradians;Counts", nbins, -0.05, 0.05);
    TH1F *hist_dummy_ytar = new TH1F("hist_dummy_ytar", "ytar Distribution;ytar;Counts", nbins, -7, 7);
    TH1F *hist_dummy_delta_up = new TH1F("hist_dummy_delta_up", "Delta Distribution (dummy, ytar > 0);W [GeV];Counts", nbins, -20, 10);
    TH1F *hist_dummy_xptar_up = new TH1F("hist_dummy_xptar_up", "xptar Distribution (dummy, ytar > 0);mili radians;Counts", nbins, -0.1, 0.1);
    TH1F *hist_dummy_yptar_up = new TH1F("hist_dummy_yptar_up", "yptar Distribution (dummy, ytar > 0);miliradians;Counts", nbins, -0.05, 0.05);
    TH1F *hist_dummy_ytar_up = new TH1F("hist_dummy_ytar_up", "ytar Distribution (dummy, ytar > 0);ytar;Counts", nbins, -7, 7);
    TH1F *hist_dummy_delta_down = new TH1F("hist_dummy_delta_down", "Delta Distribution (dummy, ytar < 0);W [GeV];Counts", nbins, -20, 10);
    TH1F *hist_dummy_xptar_down = new TH1F("hist_dummy_xptar_down", "xptar Distribution (dummy, ytar < 0);mili radians;Counts", nbins, -0.1, 0.1);
    TH1F *hist_dummy_yptar_down = new TH1F("hist_dummy_yptar_down", "yptar Distribution (dummy, ytar < 0);miliradians;Counts", nbins, -0.05, 0.05);
    TH1F *hist_dummy_ytar_down = new TH1F("hist_dummy_ytar_down", "ytar Distribution;ytar (dummy, ytar < 0);Counts", nbins, -7, 7);

    //All simc histograms shall be weighted with the fullweight
    TH1F *hist_simc_delta = new TH1F("hist_simc_delta", "Delta Distribution;W [GeV];Counts", nbins, -20, 10);
    TH1F *hist_simc_xptar = new TH1F("hist_simc_xptar", "xptar Distribution;mili radians;Counts", nbins, -0.1, 0.1);
    TH1F *hist_simc_yptar = new TH1F("hist_simc_yptar", "yptar Distribution;miliradians;Counts", nbins, -0.05, 0.05);
    TH1F *hist_simc_ytar = new TH1F("hist_simc_ytar", "ytar Distribution;ytar;Counts", nbins, -7, 7);
    TH1F *hist_simc_q2 = new TH1F("hist_simc_q2", "q2 Distribution;Q2 [GeV^2];Counts", nbins, -7, 7);
    TH1F *hist_simc_scat_ang = new TH1F("simc_data_scat_angle", "Scattering angle Distribution;angle [degrees];Counts", nbins, 20, 30);


    // Define cut conditions for Data
    double_t H_kin_W, H_kin_x_bj, H_cal_etottracknorm, H_cer_npeSum, H_gtr_dp, H_gtr_th, H_gtr_ph, H_react_z, H_gtr_y;
    double_t H_hod_goodscinhit, H_kin_scat_ang_deg;
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

    // Fill Data Histogram with cuts
    Long64_t nentries_data = data_tree->GetEntries();
    for (Long64_t i = 0; i < nentries_data; i++) {
        data_tree->GetEntry(i);
        if (H_kin_x_bj > 0.96 && H_kin_x_bj < 1.04 && 
            H_cal_etottracknorm > 0.96 && H_cal_etottracknorm < 1.04 &&
            H_hod_goodscinhit == 1 && H_cer_npeSum > 0.5 && H_gtr_dp<8.5 && H_gtr_dp>-8.5 && H_gtr_th<0.09 && H_gtr_th>-0.09 && H_gtr_ph<0.055
            && H_gtr_ph>-0.055) {
            hist_data_w->Fill(H_kin_W);
            hist_data_delta->Fill(H_gtr_dp);
            hist_data_xptar->Fill(H_gtr_th);
            hist_data_yptar->Fill(H_gtr_ph);
            hist_data_ytar->Fill(H_gtr_y);
            hist_data_scat_ang->Fill(H_kin_scat_ang_deg);

        }
    }

    // Define cut conditions for dummy
    double_t dH_kin_W, dH_kin_x_bj, dH_cal_etottracknorm, dH_cer_npeSum, dH_gtr_dp, dH_gtr_th, dH_gtr_ph, dH_gtr_y, dH_react_z;
    double_t dH_hod_goodscinhit;
    dummy_tree->SetBranchAddress("H.kin.W", &dH_kin_W);
    dummy_tree->SetBranchAddress("H.kin.x_bj", &dH_kin_x_bj);
    dummy_tree->SetBranchAddress("H.cal.etottracknorm", &dH_cal_etottracknorm);
    dummy_tree->SetBranchAddress("H.hod.goodscinhit", &dH_hod_goodscinhit);
    dummy_tree->SetBranchAddress("H.cer.npeSum", &dH_cer_npeSum);
    dummy_tree->SetBranchAddress("H.gtr.dp", &dH_gtr_dp);
    dummy_tree->SetBranchAddress("H.gtr.th", &dH_gtr_th);
    dummy_tree->SetBranchAddress("H.gtr.ph", &dH_gtr_ph);
    dummy_tree->SetBranchAddress("H.gtr.y", &dH_gtr_y);
    dummy_tree->SetBranchAddress("H.react.z", &dH_react_z);

    // Fill Dummy Histogram with cuts
    // Long64_t nentries_dummy = dummy_tree->GetEntries();
    // for (Long64_t i = 0; i < nentries_dummy; i++) {
    //     dummy_tree->GetEntry(i);
    //     if (dH_kin_x_bj > 0.96 && dH_kin_x_bj < 1.04 && 
    //         dH_cal_etottracknorm > 0.96 && dH_cal_etottracknorm < 1.04 &&
    //         dH_hod_goodscinhit == 1 && dH_cer_npeSum > 0.5 && dH_gtr_dp<8.5 && dH_gtr_dp>-8.5 && dH_gtr_th<0.09 && dH_gtr_th>-0.09 && dH_gtr_ph<0.055
    //         && dH_gtr_ph>-0.055) {
    //         hist_dummy_w->Fill(dH_kin_W);
    //         hist_dummy_delta->Fill(dH_kin_W);
    //         hist_dummy_xptar->Fill(dH_gtr_th);
    //         hist_dummy_yptar->Fill(dH_gtr_ph);
    //         hist_dummy_ytar->Fill(dH_gtr_y);
    //     }
    // }

    // Loop over dummy data tree
    for (int i = 0; i < dummy_tree->GetEntries(); i++) {
        dummy_tree->GetEntry(i);
    
        if (dH_kin_x_bj > 0.96 && dH_kin_x_bj < 1.04 && 
            dH_cal_etottracknorm > 0.96 && dH_cal_etottracknorm < 1.04 &&
            dH_hod_goodscinhit == 1 && dH_cer_npeSum > 0.5 && dH_gtr_dp<8.5 && dH_gtr_dp>-8.5 && dH_gtr_th<0.09 && dH_gtr_th>-0.09 && dH_gtr_ph<0.055
            && dH_gtr_ph>-0.055 && dH_gtr_y > 0) {
            hist_dummy_w_down->Fill(dH_kin_W);
            hist_dummy_delta_down->Fill(dH_gtr_dp);
            hist_dummy_xptar_down->Fill(dH_gtr_th);
            hist_dummy_yptar_down->Fill(dH_gtr_ph);
            hist_dummy_ytar_down->Fill(dH_gtr_y);
        } 
        else if (dH_kin_x_bj > 0.96 && dH_kin_x_bj < 1.04 && 
            dH_cal_etottracknorm > 0.96 && dH_cal_etottracknorm < 1.04 &&
            dH_hod_goodscinhit == 1 && dH_cer_npeSum > 0.5 && dH_gtr_dp<8.5 && dH_gtr_dp>-8.5 && dH_gtr_th<0.09 && dH_gtr_th>-0.09 && dH_gtr_ph<0.055
            && dH_gtr_ph>-0.055 && dH_gtr_y < 0) {
            hist_dummy_w_up->Fill(dH_kin_W);
            hist_dummy_delta_up->Fill(dH_gtr_dp);
            hist_dummy_xptar_up->Fill(dH_gtr_th);
            hist_dummy_yptar_up->Fill(dH_gtr_ph);
            hist_dummy_ytar_up->Fill(dH_gtr_y);
        }
    }

    // Define SIMC branches
    float_t W, Weight;
    simc_tree->SetBranchAddress("W", &W);
    simc_tree->SetBranchAddress("Weight", &Weight);

    float_t fhsdelta, fhsyptar, fhsxptar, fQ2, fnu, fhsytar, fhs; 

    simc_tree->SetBranchAddress("hsdelta", &fhsdelta);
    simc_tree->SetBranchAddress("hsyptar", &fhsyptar);
    simc_tree->SetBranchAddress("hsxptar", &fhsxptar);
    simc_tree->SetBranchAddress("hsytar", &fhsytar);
    simc_tree->SetBranchAddress("Q2", &fQ2);
    simc_tree->SetBranchAddress("nu", &fnu);

    // Apply Weight and Normalize SIMC Histogram

    // double_t normfac = 0.671016E+07; // HMS 4042 run 6834 Ebeam 6.3967
    double_t normfac = 0.159800E+08;
    double_t nevents = 100000;
    double_t weight_factor = normfac / nevents;
    double_t full_weight;
    double_t hms_tracking_effic, hgc_cerenkov_effic, calo_effic, hod_effic;
    double_t effective_charge_data;
    double_t effective_charge_dummy;
    double_t scale_factor_dummy_upstream = 8.467;
    double_t scale_factor_dummy_downstream = 4.256; //using tip
    // double_t scale_factor_dummy_downstream = 3.711; //using wall
    double_t mass_p = 0.938272; //mass of proton in GeV
    double_t f_xbj; 
    double_t live_time= 0.999931;//Pre-Scaled Ps4 HMS Computer Live Time : 99.9931 %

    Long64_t nentries_simc = simc_tree->GetEntries();

    for (Long64_t i = 0; i < nentries_simc; i++) {
        simc_tree->GetEntry(i);
        full_weight = Weight*weight_factor;
        f_xbj = (fQ2)/(2*mass_p*fnu);
        if (f_xbj > 0.96 && f_xbj < 1.04 && fhsdelta < 8.5 && fhsdelta > -8.5 && fhsxptar > -0.09 && fhsxptar < 0.09 && fhsyptar > -0.055 && fhsyptar < 0.055) {
            hist_simc_weighted_w->Fill(W, full_weight);
            hist_simc_delta->Fill(fhsdelta, full_weight);
            hist_simc_xptar->Fill(fhsxptar, full_weight);
            hist_simc_yptar->Fill(fhsyptar, full_weight);
            hist_simc_ytar->Fill(fhsytar, full_weight);
        }
        
    }

    // Normalize histograms
    // if (hist_data->Integral() > 0) hist_data->Scale(1.0 / hist_data->Integral());
    // effective_charge_data = 46266.080/1000.0; //BCM4A charge in mC run 6834
    
    // efficiencies from the REPORT_OUTPUT_pass1
    hms_tracking_effic = 0.9955;
    hgc_cerenkov_effic = 0.990136;
    hod_effic = 0.999618;
    
    effective_charge_data = 46266.080/(1000*hms_tracking_effic*hgc_cerenkov_effic*hod_effic); 
    effective_charge_dummy = 22103.178/1000.0;

    hist_data_w->Scale(1.0 / effective_charge_data);

    // hist_dummy_w->Scale(scale_factor_dummy_upstream / effective_charge_dummy);
    hist_dummy_w_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_w_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));
    hist_dummy_delta_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_delta_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));
    hist_dummy_xptar_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_xptar_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));
    hist_dummy_yptar_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_yptar_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));
    hist_dummy_ytar_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_ytar_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));

    // Combine into final dummy histogram
    hist_dummy_w->Add(hist_dummy_w_up);
    hist_dummy_w->Add(hist_dummy_w_down);
    hist_dummy_delta->Add(hist_dummy_delta_up);
    hist_dummy_delta->Add(hist_dummy_delta_down);
    hist_dummy_xptar->Add(hist_dummy_xptar_up);
    hist_dummy_xptar->Add(hist_dummy_xptar_down);
    hist_dummy_yptar->Add(hist_dummy_yptar_up);
    hist_dummy_yptar->Add(hist_dummy_yptar_down);
    hist_dummy_ytar->Add(hist_dummy_ytar_up);
    hist_dummy_ytar->Add(hist_dummy_ytar_down);
    

    hist_exp_w->Add(hist_data_w, hist_dummy_w, 1.0, -1.0); // hist_exp = hist_data - hist_dummy

    // // Get mean and standard deviation
    // double mean_data = hist_data->GetMean();
    // double stddev_data = hist_data->GetStdDev();
    // double mean_simc = hist_simc_weighted->GetMean();
    // double stddev_simc = hist_simc_weighted->GetStdDev();
    // double mean_exp = hist_exp->GetMean();
    // double stddev_exp = hist_exp->GetStdDev();

    // Fit histograms with a Gaussian
    TF1 *fit_data = new TF1("fit_data", "gaus", 0.82, 0.96);
    TF1 *fit_simc_weighted = new TF1("fit_simc_weighted", "gaus", 0.92, 0.96);
    TF1 *fit_exp = new TF1("fit_exp", "gaus", 0.88, 0.96);

    hist_data_w->Fit(fit_data, "RQ");
    hist_simc_weighted_w->Fit(fit_simc_weighted, "RQ");
    hist_exp_w->Fit(fit_exp, "RQ");

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
    hist_data_w->SetLineColor(kBlue);
    // hist_data->SetLineWidth(1.5);
    hist_simc_weighted_w->SetLineColor(kRed);
    hist_exp_w->SetLineColor(kBlack);
    hist_dummy_w->SetLineColor(kGreen);
    // hist_exp->SetLineWidth(2);

    double maxW = std::max({hist_simc_weighted_w->GetMaximum(), hist_data_w->GetMaximum(), hist_exp_w->GetMaximum(), hist_dummy_w->GetMaximum()});

    hist_simc_weighted_w->SetMaximum(1.2 * maxW);  // Add some margin

    hist_simc_weighted_w->Draw("hist");
    hist_data_w->Draw("hist same");
    hist_exp_w->Draw("hist same");
    hist_dummy_w->Draw("hist same");

    // double_t yield_data = hist_data_w->Integral();
    double_t yield_exp = hist_exp_w->Integral();
    double_t yield_simc_weighted = hist_simc_weighted_w->Integral();
    double_t yield_ratio = (yield_simc_weighted != 0) ? (yield_exp / yield_simc_weighted) : 0.0; // Prevent division by zero
    
    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->AddEntry(hist_data_w, Form("Data: Mean = %.5f, #sigma = %.5f", mean_data, sigma_data), "l");
    leg->AddEntry(hist_simc_weighted_w, Form("SIMC: Mean = %.5f, #sigma = %.5f", mean_simc, sigma_simc), "l");
    leg->AddEntry(hist_exp_w, Form("Exp = Data - Dummy: Mean = %.5f, #sigma = %.5f", mean_exp, sigma_exp), "l");
    leg->AddEntry((TObject*)0, Form("Yield Ratio = Yield_exp/Yield_simc: %.5f", yield_ratio), ""); // Display yield ratio safely
    leg->Draw();
    

    // Show plot interactively
    canvas->Update();

    TCanvas *canvas_delta = new TCanvas("canvas_delta", "Delta Distribution", 800, 600);
    hist_data_delta->SetLineColor(kBlue);
    hist_simc_delta->SetLineColor(kRed);
    hist_data_delta->Scale(1.0 / effective_charge_data);
    hist_data_delta->Add(hist_dummy_delta, -1);
    double maxDelta = std::max(hist_simc_delta->GetMaximum(), hist_data_delta->GetMaximum());
    hist_simc_delta->SetMaximum(1.2 * maxDelta);
    hist_simc_delta->Draw("hist");
    hist_data_delta->Draw("hist same");
    // hist_simc_delta->Draw("hist same");
    canvas_delta->Update();

    TCanvas *canvas_xptar = new TCanvas("canvas_xptar", "xptar Distribution", 800, 600);
    hist_data_xptar->SetLineColor(kBlue);
    hist_simc_xptar->SetLineColor(kRed);
    hist_data_xptar->Scale(1.0 / effective_charge_data);
    hist_data_xptar->Add(hist_dummy_xptar, -1);
    double maxXptar = std::max(hist_simc_xptar->GetMaximum(), hist_data_xptar->GetMaximum());
    hist_simc_xptar->SetMaximum(1.2 * maxXptar);
    hist_simc_xptar->Draw("hist");
    hist_data_xptar->Draw("hist same");
    // hist_simc_xptar->Draw("hist same");
    canvas_xptar->Update();
    
    TCanvas *canvas_yptar = new TCanvas("canvas_yptar", "yptar Distribution", 800, 600);
    hist_data_yptar->SetLineColor(kBlue);
    hist_simc_yptar->SetLineColor(kRed);
    hist_data_yptar->Scale(1.0 / effective_charge_data);
    hist_data_yptar->Add(hist_dummy_yptar, -1);
    double maxYptar = std::max(hist_simc_yptar->GetMaximum(), hist_data_yptar->GetMaximum());
    hist_simc_yptar->SetMaximum(1.2 * maxYptar);
    hist_simc_yptar->Draw("hist");
    hist_data_yptar->Draw("hist same");
    // hist_simc_yptar->Draw("hist same");
    canvas_yptar->Update();

    TCanvas *canvas_ytar = new TCanvas("canvas_ytar", "ytar Distribution", 800, 600);
    hist_data_ytar->SetLineColor(kBlue);
    hist_simc_ytar->SetLineColor(kRed);
    hist_dummy_ytar->SetLineColor(kGreen);
    hist_data_ytar->Scale(1.0 / effective_charge_data);
    hist_data_ytar->Add(hist_dummy_ytar, -1);
    double maxYtar = std::max({hist_simc_ytar->GetMaximum(), hist_data_ytar->GetMaximum(),hist_dummy_ytar->GetMaximum()});
    hist_simc_ytar->SetMaximum(1.2 * maxYtar);
    hist_simc_ytar->Draw("hist");
    hist_data_ytar->Draw("hist same");
    hist_dummy_ytar->Draw("hist same");
    canvas_ytar->Update();


}
