#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TStyle.h"

void compare_data_simc() {
    // Turn off stat boxes
    gStyle->SetOptStat(0);

    // Load Data and SIMC File
    int run = 6834;
    std::string filename = Form("/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/nps_hms_coin_%d_0_1_-1.root", run);
    TFile *data_file = new TFile(filename.c_str(), "READ");
    TFile *simc_file = new TFile("/u/group/nps/singhav/simc_gfortran/worksim/eep_hydrogen_6828_6841_HMS4042.root", "READ");

    if (!data_file || data_file->IsZombie() || !simc_file || simc_file->IsZombie()) {
        std::cerr << "Error: Unable to open files!" << std::endl;
        return;
    }

    TTree *data_tree = (TTree*) data_file->Get("T");
    TTree *simc_tree = (TTree*) simc_file->Get("h10");

    if (!data_tree || !simc_tree) {
        std::cerr << "Error: Unable to get trees!" << std::endl;
        return;
    }

    // Create Canvas with two pads
    TCanvas *c1 = new TCanvas("c1", "Data vs SIMC", 1200, 600);
    c1->Divide(2, 1); // Two side-by-side pads

    // Define histograms
    TH2F *h_simc = new TH2F("h_simc", "SIMC: hsdelta vs hsyptar;hsyptar (rad); hsdelta (%)", 
                            100, -0.05, 0.05, 100, -10, 10);
    TH2F *h_data = new TH2F("h_data", "Data: H.gtr.dp vs H.gtr.ph;H.gtr.ph (rad); H.gtr.dp (%)", 
                            100, -0.05, 0.05, 100, -10, 10);

    // Fill histograms
    simc_tree->Draw("hsdelta:hsyptar>>h_simc", 
                    "W>0.88 && W<1 && hsdelta<8.5 && hsdelta>-8.5 && hsxptar<0.09 && hsyptar>-0.055", "goff");
    data_tree->Draw("H.gtr.dp:H.gtr.ph>>h_data", 
                    "H.kin.W>0.88 && H.kin.W<1 && H.cal.etottracknorm>0.9 && H.hod.goodscinhit==1 && H.cer.npeSum>0.5 && H.gtr.dp<8.5 && H.gtr.dp>-8.5 && H.gtr.th<0.09 && H.gtr.th>-0.09 && H.gtr.ph<0.055 && H.gtr.ph>-0.055", 
                    "goff");

    // --- SIMC FIT ---
    c1->cd(1);
    h_simc->Draw("COLZ");

    // Fit SIMC histogram
    TF1 *fit_simc = new TF1("fit_simc", "pol1", -10, 10);
    h_simc->Fit("fit_simc", "R");
    fit_simc->SetLineColor(kRed);

    // Extract fit parameters
    double m_simc = fit_simc->GetParameter(1);
    double b_simc = fit_simc->GetParameter(0);

    // Draw fit manually using TGraph
    TGraph *fit_graph_simc = new TGraph();
    for (double x = -10; x <= 10; x += 0.1) {
        fit_graph_simc->SetPoint(fit_graph_simc->GetN(), x, m_simc * x + b_simc);
    }
    fit_graph_simc->SetLineColor(kRed);
    fit_graph_simc->SetLineWidth(2);
    fit_graph_simc->Draw("L SAME");

    // Display fit equation for SIMC
    TLatex latex_simc;
    latex_simc.SetNDC();
    latex_simc.SetTextSize(0.04);
    latex_simc.SetTextColor(kRed);
    latex_simc.DrawLatex(0.15, 0.85, Form("y = %.3fx + %.3f", m_simc, b_simc));

    // --- DATA FIT ---
    c1->cd(2);
    h_data->Draw("COLZ");

    // Fit Data histogram
    TF1 *fit_data = new TF1("fit_data", "pol1", -10, 10);
    h_data->Fit("fit_data", "R");
    fit_data->SetLineColor(kBlue);

    // Extract fit parameters
    double m_data = fit_data->GetParameter(1);
    double b_data = fit_data->GetParameter(0);

    // Draw fit manually using TGraph
    TGraph *fit_graph_data = new TGraph();
    for (double x = -10; x <= 10; x += 0.1) {
        fit_graph_data->SetPoint(fit_graph_data->GetN(), x, m_data * x + b_data);
    }
    fit_graph_data->SetLineColor(kBlue);
    fit_graph_data->SetLineWidth(2);
    fit_graph_data->Draw("L SAME");

    // Display fit equation for Data
    TLatex latex_data;
    latex_data.SetNDC();
    latex_data.SetTextSize(0.04);
    latex_data.SetTextColor(kBlue);
    latex_data.DrawLatex(0.15, 0.85, Form("y = %.3fx + %.3f", m_data, b_data));

    // Show canvas
    c1->Update();
}
