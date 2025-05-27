#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void plot_W_distribution(){
    int run = 6705;
    TFile *dummy_file = new TFile("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_6705_0_1_-1.root", "READ");
    TTree *dummy_tree = (TTree*) dummy_file->Get("T");

    double_t nbins = 100;
    double_t xmin = 0.6, xmax = 1.4;

    TH1F *hist_dummy_w = new TH1F("hist_dummy_w", "W Distribution;W [GeV];Counts", nbins, xmin, xmax);
    TH1F* hist_dummy_w_up = new TH1F("hist_dummy_w_up", "W distribution (dummy, ytar > 0)", nbins, xmin, xmax);
    TH1F* hist_dummy_w_down = new TH1F("hist_dummy_w_down", "W distribution (dummy, ytar < 0)", nbins, xmin, xmax);

    TH1F *hist_dummy_delta = new TH1F("hist_dummy_delta", "Delta Distribution;W [GeV];Counts", nbins, -20, 10);
    TH1F *hist_dummy_xptar = new TH1F("hist_dummy_xptar", "xptar Distribution;mili radians;Counts", nbins, -0.1, 0.1);
    TH1F *hist_dummy_yptar = new TH1F("hist_dummy_yptar", "yptar Distribution;miliradians;Counts", nbins, -0.05, 0.05);
    TH1F *hist_dummy_ytar = new TH1F("hist_dummy_ytar", "ytar Distribution;ytar;Counts", nbins, -7, 7);
    TH1F *hist_dummy_q2 = new TH1F("hist_dummy_q2", "Q2 Distribution (dummy);ytar;Counts", nbins, 3, 6);
    TH1F *hist_dummy_x_fp = new TH1F("hist_dummy_x_fp", "x focal plane (dummy)", nbins, -60,30);
    TH1F *hist_dummy_y_fp = new TH1F("hist_dummy_y_fp", "y focal plane (dummy)", nbins, -30 , 30);
    TH1F *hist_dummy_xp_fp = new TH1F("hist_dummy_xp_fp", "xp focal plane (dummy)", nbins, -0.06, 0.06);
    TH1F *hist_dummy_yp_fp = new TH1F("hist_dummy_yp_fp", "yp focal plane (dummy)", nbins, -0.02 , 0.02);

    TH1F *hist_dummy_delta_up = new TH1F("hist_dummy_delta_up", "Delta Distribution (dummy, ytar > 0);W [GeV];Counts", nbins, -20, 10);
    TH1F *hist_dummy_xptar_up = new TH1F("hist_dummy_xptar_up", "xptar Distribution (dummy, ytar > 0);mili radians;Counts", nbins, -0.1, 0.1);
    TH1F *hist_dummy_yptar_up = new TH1F("hist_dummy_yptar_up", "yptar Distribution (dummy, ytar > 0);miliradians;Counts", nbins, -0.05, 0.05);
    TH1F *hist_dummy_ytar_up = new TH1F("hist_dummy_ytar_up", "ytar Distribution (dummy, ytar > 0);ytar;Counts", nbins, -7, 7);
    TH1F *hist_dummy_q2_up = new TH1F("hist_dummy_q2_up", "Q2 Distribution (dummy, ytar > 0);ytar;Counts", nbins, 3, 6);
    TH1F *hist_dummy_x_fp_up = new TH1F("hist_dummy_x_fp_up", "x focal plane (dummy, ytar > 0)", nbins, -60,30);
    TH1F *hist_dummy_y_fp_up = new TH1F("hist_dummy_y_fp_up", "y focal plane (dummy, ytar > 0)", nbins, -30 , 30);
    TH1F *hist_dummy_xp_fp_up = new TH1F("hist_dummy_xp_fp_up", "xp focal plane (dummy, ytar > 0)", nbins, -0.06, 0.06);
    TH1F *hist_dummy_yp_fp_up = new TH1F("hist_dummy_yp_fp_up", "yp focal plane (dummy, ytar > 0)", nbins, -0.02 , 0.02);
    
    TH1F *hist_dummy_delta_down = new TH1F("hist_dummy_delta_down", "Delta Distribution (dummy, ytar < 0);W [GeV];Counts", nbins, -20, 10);
    TH1F *hist_dummy_xptar_down = new TH1F("hist_dummy_xptar_down", "xptar Distribution (dummy, ytar < 0);mili radians;Counts", nbins, -0.1, 0.1);
    TH1F *hist_dummy_yptar_down = new TH1F("hist_dummy_yptar_down", "yptar Distribution (dummy, ytar < 0);miliradians;Counts", nbins, -0.05, 0.05);
    TH1F *hist_dummy_ytar_down = new TH1F("hist_dummy_ytar_down", "ytar Distribution;ytar (dummy, ytar < 0);Counts", nbins, -7, 7);
    TH1F *hist_dummy_q2_down = new TH1F("hist_dummy_q2_down", "Q2 Distribution (dummy, ytar > 0);ytar;Counts", nbins, 3, 6);
    TH1F *hist_dummy_x_fp_down = new TH1F("hist_dummy_x_fp_down", "x focal plane (dummy, ytar < 0)", nbins, -60,30);
    TH1F *hist_dummy_y_fp_down = new TH1F("hist_dummy_y_fp_down", "y focal plane (dummy, ytar < 0)", nbins, -30 , 30);
    TH1F *hist_dummy_xp_fp_down = new TH1F("hist_dummy_xp_fp_down", "xp focal plane (dummy, ytar < 0)", nbins, -0.06, 0.06);
    TH1F *hist_dummy_yp_fp_down = new TH1F("hist_dummy_yp_fp_down", "yp focal plane (dummy, ytar < 0)", nbins, -0.02 , 0.02);

    double_t dH_kin_W, dH_kin_x_bj, dH_cal_etottracknorm, dH_cer_npeSum, dH_gtr_dp, dH_gtr_th, dH_gtr_ph, dH_gtr_y, dH_react_z, dH_kin_Q2, dH_dc_x_fp, dH_dc_y_fp, dH_dc_xp_fp, dH_dc_yp_fp;
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
    dummy_tree->SetBranchAddress("H.kin.Q2", &dH_kin_Q2);
    dummy_tree->SetBranchAddress("H.dc.x_fp", &dH_dc_x_fp);
    dummy_tree->SetBranchAddress("H.dc.y_fp", &dH_dc_y_fp);
    dummy_tree->SetBranchAddress("H.dc.xp_fp", &dH_dc_xp_fp);
    dummy_tree->SetBranchAddress("H.dc.yp_fp", &dH_dc_yp_fp);

    // Loop over dummy data tree
    for (int i = 0; i < dummy_tree->GetEntries(); i++) {
        dummy_tree->GetEntry(i);
    
        if (dH_kin_W > 0.92 && dH_kin_W < 0.96 && 
            dH_cal_etottracknorm > 0.9 &&
            dH_cer_npeSum > 0.5 && dH_gtr_dp<8.5 && dH_gtr_dp>-8.5 && dH_gtr_th<0.09 && dH_gtr_th>-0.09 && dH_gtr_ph<0.055
            && dH_gtr_ph>-0.055 && dH_gtr_y > 0) {
            hist_dummy_w_down->Fill(dH_kin_W);
            hist_dummy_delta_down->Fill(dH_gtr_dp);
            hist_dummy_xptar_down->Fill(dH_gtr_th);
            hist_dummy_yptar_down->Fill(dH_gtr_ph);
            hist_dummy_ytar_down->Fill(dH_gtr_y);
            hist_dummy_q2_down->Fill(dH_kin_Q2);
            hist_dummy_x_fp_down->Fill(dH_dc_x_fp);
            hist_dummy_y_fp_down->Fill(dH_dc_y_fp);
            hist_dummy_xp_fp_down->Fill(dH_dc_xp_fp);
            hist_dummy_yp_fp_down->Fill(dH_dc_yp_fp);
        } 
        else if (dH_kin_W > 0.92 && dH_kin_W < .96 && 
            dH_cal_etottracknorm > 0.9 &&
            dH_cer_npeSum > 0.5 && dH_gtr_dp<8.5 && dH_gtr_dp>-8.5 && dH_gtr_th<0.09 && dH_gtr_th>-0.09 && dH_gtr_ph<0.055
            && dH_gtr_ph>-0.055 && dH_gtr_y < 0) {
            hist_dummy_w_up->Fill(dH_kin_W);
            hist_dummy_delta_up->Fill(dH_gtr_dp);
            hist_dummy_xptar_up->Fill(dH_gtr_th);
            hist_dummy_yptar_up->Fill(dH_gtr_ph);
            hist_dummy_ytar_up->Fill(dH_gtr_y);
            hist_dummy_q2_up->Fill(dH_kin_Q2);
            hist_dummy_x_fp_up->Fill(dH_dc_x_fp);
            hist_dummy_y_fp_up->Fill(dH_dc_y_fp);
            hist_dummy_xp_fp_up->Fill(dH_dc_xp_fp);
            hist_dummy_yp_fp_up->Fill(dH_dc_yp_fp);
        }
    }

    double_t effective_charge_dummy = 22103.178/1000.0;
    // double_t scale_factor_dummy_upstream = 8.467;
    // double_t scale_factor_dummy_downstream = 4.256;

    double_t scale_factor_dummy_upstream = 1.0; //for unscaled plots
    double_t scale_factor_dummy_downstream = 1.0; //for unscaled plots

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
    hist_dummy_q2_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_q2_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));
    hist_dummy_x_fp_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_x_fp_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));
    hist_dummy_y_fp_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_y_fp_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));
    hist_dummy_xp_fp_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_xp_fp_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));
    hist_dummy_yp_fp_up->Scale(1.0 /(scale_factor_dummy_upstream*effective_charge_dummy));
    hist_dummy_yp_fp_down->Scale(1.0 /(scale_factor_dummy_downstream*effective_charge_dummy));

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
    hist_dummy_q2->Add(hist_dummy_q2_up);
    hist_dummy_q2->Add(hist_dummy_q2_down);
    hist_dummy_x_fp->Add(hist_dummy_x_fp_up);
    hist_dummy_x_fp->Add(hist_dummy_x_fp_down);
    hist_dummy_y_fp->Add(hist_dummy_y_fp_up);
    hist_dummy_y_fp->Add(hist_dummy_y_fp_down);
    hist_dummy_xp_fp->Add(hist_dummy_xp_fp_up);
    hist_dummy_xp_fp->Add(hist_dummy_xp_fp_down);
    hist_dummy_yp_fp->Add(hist_dummy_yp_fp_up);
    hist_dummy_yp_fp->Add(hist_dummy_yp_fp_down);

    // Disable default stat box
    gStyle->SetOptStat(0);

    // Create Canvas
    TCanvas *canvas = new TCanvas("canvas", "W Distribution", 800, 600);

    hist_dummy_w->Draw("hist");


    // Show plot interactively
    canvas->Update();
    std::string save_filename1 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/W_%d.png", run);
    canvas->SaveAs(save_filename1.c_str());

    TCanvas *canvas_delta = new TCanvas("canvas_delta", "Delta Distribution", 800, 600);
    hist_dummy_delta->Draw("hist");
    canvas_delta->Update();
    std::string save_filename2 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/delta_%d.png", run);
    canvas_delta->SaveAs(save_filename2.c_str());

    TCanvas *canvas_xptar = new TCanvas("canvas_xptar", "xptar Distribution", 800, 600);
    hist_dummy_xptar->Draw("hist");
    // hist_simc_xptar->Draw("hist same");
    canvas_xptar->Update();
    std::string save_filename3 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/xptar_%d.png", run);
    canvas_xptar->SaveAs(save_filename3.c_str());
    
    TCanvas *canvas_yptar = new TCanvas("canvas_yptar", "yptar Distribution", 800, 600);
    hist_dummy_yptar->Draw("hist");
    // hist_simc_yptar->Draw("hist same");
    canvas_yptar->Update();
    std::string save_filename4 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/yptar_%d.png", run);
    canvas_yptar->SaveAs(save_filename4.c_str());

    TCanvas *canvas_ytar = new TCanvas("canvas_ytar", "ytar Distribution", 800, 600);
    hist_dummy_ytar->SetLineWidth(2);
    hist_dummy_ytar->Draw("hist");
    hist_dummy_ytar->SetTitle("ytar Distribution;ytar (cm);Counts");
    canvas_ytar->SetGrid();
    canvas_ytar->Update();
    std::string save_filename5 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/ytar_%d.png", run);
    canvas_ytar->SaveAs(save_filename5.c_str());

    TCanvas *canvas_q2 = new TCanvas("canvas_q2", "Q2 Distribution", 800, 600);
    hist_dummy_q2->Draw("hist");
    canvas_q2->Update();
    std::string save_filename6 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/q2_%d.png", run);
    canvas_q2->SaveAs(save_filename6.c_str());

    TCanvas *canvas_x_fp = new TCanvas("x_fp_delta", "x focal plane Distribution", 800, 600);
    hist_dummy_x_fp->Draw("hist");
    // hist_simc_delta->Draw("hist same");
    canvas_x_fp->Update();
    std::string save_filename7 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/xfp_%d.png", run);
    canvas_x_fp->SaveAs(save_filename7.c_str());

    TCanvas *canvas_y_fp = new TCanvas("y_fp_delta", "y focal plane Distribution", 800, 600);
    hist_dummy_y_fp->Draw("hist");
    // hist_simc_delta->Draw("hist same");
    canvas_y_fp->Update();
    std::string save_filename8 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/yfp_%d.png", run);
    canvas_y_fp->SaveAs(save_filename8.c_str());

    TCanvas *canvas_xp_fp = new TCanvas("xp_fp_delta", "xp focal plane Distribution", 800, 600);
    hist_dummy_xp_fp->Draw("hist");
    // hist_simc_delta->Draw("hist same");
    canvas_xp_fp->Update();
    std::string save_filename9 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/xpfp_%d.png", run);
    canvas_xp_fp->SaveAs(save_filename9.c_str());

    TCanvas *canvas_yp_fp = new TCanvas("yp_fp_delta", "yp focal plane Distribution", 800, 600);
    hist_dummy_yp_fp->Draw("hist");
    // hist_simc_delta->Draw("hist same");
    canvas_yp_fp->Update();
    std::string save_filename10 = Form("/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/dummy_plots_unscaled/ypfp_%d.png", run);
    canvas_yp_fp->SaveAs(save_filename10.c_str());
}