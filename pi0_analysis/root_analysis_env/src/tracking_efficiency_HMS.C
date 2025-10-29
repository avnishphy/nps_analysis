#include "utils.C"  // logmsg(), trim(), etc.

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

// ----------------------
// Runlist reader
// ----------------------
vector<int> readRunList(const string &fname) {
    vector<int> runs;
    ifstream in(fname);
    if (!in.is_open()) { logmsg(WARN, "Cannot open runlist: " + fname); return runs; }
    string line;
    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        try { runs.push_back(stoi(line)); } catch (...) { logmsg(WARN, "Skipping invalid runlist line: " + line); }
    }
    return runs;
}

// ============================================================
// HMS Tracking Efficiency Analysis
// ============================================================
void tracking_efficiency_HMS(const TString &skimDir_in="output/skimmed/",
                             const TString &outPlotDir_in="output/plots/efficiency/",
                             const TString &runlistFile="config/runlist_x60_4b.txt") {

    // --------------------------- Initialize
    TStopwatch sw_total;
    sw_total.Start();
    logmsg(INFO, "=========== HMS Tracking Efficiency ===========");

    TString skimDir = skimDir_in.EndsWith("/") ? skimDir_in : skimDir_in + "/";
    TString outPlotDir = outPlotDir_in.EndsWith("/") ? outPlotDir_in : outPlotDir_in + "/";
    gSystem->mkdir(outPlotDir, true);

    vector<int> runs = readRunList(runlistFile.Data());
    if (runs.empty()) { logmsg(ERROR, "No runs found!"); return; }

    // --------------------------- Global histograms
    TH1D *h_beta_notrack = new TH1D("h_beta_notrack", "HMS Beta (no track);#beta;Counts", 100, 0.4, 1.4);
    TH1D *h_etotnorm    = new TH1D("h_etotnorm", "E_{tot}/p;E_{tot}/p;Counts", 100, 0.4, 1.4);
    TH1D *h_npe         = new TH1D("h_npe", "HMS Cherenkov NPE;NPE;Counts", 100, 0, 15);

    map<int, pair<double,double>> runEffMap; // run -> (eff, eff_unc)
    double total_Nshould = 0;
    double total_Ndid = 0;

    // ============================================================ Loop over runs
    for (int run : runs) {
        TString infile = Form("%sskim_run%d.root", skimDir.Data(), run);
        if (gSystem->AccessPathName(infile)) {
            logmsg(WARN, Form("Skipping run %d: file not found.", run));
            continue;
        }

        logmsg(INFO, Form("Processing run %d", run));
        TFile *f = TFile::Open(infile, "READ");
        if (!f || f->IsZombie()) { logmsg(ERROR, Form("Error opening file for run %d", run)); continue; }

        TTree *T = (TTree*)f->Get("T");
        if (!T) { logmsg(ERROR, Form("Tree not found in run %d", run)); f->Close(); continue; }

        // --------------------------- Branch setup
        double_t H_hod_betanotrack, H_cal_etotnorm, H_cer_npeSum, H_hod_goodscinhit, H_dc_ntrack;

        T->SetBranchStatus("*", 0);
        T->SetBranchStatus("H.hod.betanotrack", 1);  T->SetBranchAddress("H.hod.betanotrack", &H_hod_betanotrack);
        T->SetBranchStatus("H.cal.etotnorm", 1);     T->SetBranchAddress("H.cal.etotnorm", &H_cal_etotnorm);
        T->SetBranchStatus("H.cer.npeSum", 1);       T->SetBranchAddress("H.cer.npeSum", &H_cer_npeSum);
        T->SetBranchStatus("H.hod.goodscinhit", 1);  T->SetBranchAddress("H.hod.goodscinhit", &H_hod_goodscinhit);
        T->SetBranchStatus("H.dc.ntrack", 1);        T->SetBranchAddress("H.dc.ntrack", &H_dc_ntrack);

        Long64_t nentries = T->GetEntries();
        Long64_t Nshould = 0, Ndid = 0;

        // --------------------------- Event loop
        for (Long64_t i=0; i<nentries; i++) {
            T->GetEntry(i);

            bool should = (H_hod_goodscinhit == 1 &&
                           H_hod_betanotrack > 0.5 && H_hod_betanotrack < 1.5 &&
                           H_cal_etotnorm > 0.6 &&
                           H_cer_npeSum > 0.6);

            if (should) {
                Nshould++;
                h_beta_notrack->Fill(H_hod_betanotrack);
                h_etotnorm->Fill(H_cal_etotnorm);
                h_npe->Fill(H_cer_npeSum);

                if (H_dc_ntrack > 0) Ndid++;
            }
        }

        f->Close();

        // --------------------------- Compute run efficiency
        double eff = (Nshould>0) ? (double)Ndid/Nshould : 0.0;
        double eff_unc = (Nshould>0) ? sqrt(Nshould-Ndid)/Nshould : 0.0;
        runEffMap[run] = {eff, eff_unc};
        logmsg(INFO, Form("Run %d: Nshould=%lld, Ndid=%lld, Eff=%.4f ± %.4f",
                          run, Nshould, Ndid, eff, eff_unc));

        total_Nshould += Nshould;
        total_Ndid += Ndid;

        // --------------------------- Save run-wise histograms (clone to avoid deleting global hist)
        TH1D *h_beta_run = (TH1D*)h_beta_notrack->Clone(Form("h_beta_run%d", run));
        TH1D *h_etot_run  = (TH1D*)h_etotnorm->Clone(Form("h_etot_run%d", run));
        TH1D *h_npe_run   = (TH1D*)h_npe->Clone(Form("h_npe_run%d", run));

        TFile fout_run(Form("%shms_tracking_eff_run%d.root", outPlotDir.Data(), run), "RECREATE");
        h_beta_run->Write();
        h_etot_run->Write();
        h_npe_run->Write();
        fout_run.Close();

        delete h_beta_run;
        delete h_etot_run;
        delete h_npe_run;
    }

    // ============================================================ Efficiency vs run
    vector<double> run_nums, eff_vals, eff_errs;
    for (auto &[run, vals] : runEffMap) {
        run_nums.push_back(run);
        eff_vals.push_back(vals.first);
        eff_errs.push_back(vals.second);
    }

    TGraphErrors *g_eff = new TGraphErrors(run_nums.size());
    for (size_t i=0; i<run_nums.size(); i++) {
        g_eff->SetPoint(i, run_nums[i], eff_vals[i]);
        g_eff->SetPointError(i, 0, eff_errs[i]);
    }
    g_eff->SetTitle("HMS Tracking Efficiency vs Run;Run Number;Tracking Efficiency");
    g_eff->SetMarkerStyle(20); g_eff->SetMarkerSize(1.0);
    g_eff->SetLineWidth(2); g_eff->SetLineColor(kBlue+1); g_eff->SetMarkerColor(kBlue+2);

    TCanvas *c_eff = new TCanvas("c_eff","Tracking Efficiency",800,600);
    c_eff->SetGrid();
    g_eff->Draw("AP");
    c_eff->SaveAs(outPlotDir + "hms_tracking_eff_vs_run.pdf");

    // ============================================================ Overall efficiency
    double total_eff = (total_Nshould>0) ? total_Ndid/total_Nshould : 0.0;
    double total_eff_unc = (total_Nshould>0) ? sqrt(total_Nshould-total_Ndid)/total_Nshould : 0.0;
    logmsg(INFO, "----------------------------------------------");
    logmsg(INFO, Form("Overall HMS Tracking Efficiency = %.4f ± %.4f", total_eff, total_eff_unc));
    logmsg(INFO, "----------------------------------------------");

    // ============================================================ Save CSV
    TString csvFile = outPlotDir + "hms_tracking_efficiency.csv";
    ofstream csv(csvFile.Data());
    csv << "Run,Efficiency,Uncertainty\n";
    for (auto &[run, vals] : runEffMap)
        csv << run << "," << vals.first << "," << vals.second << "\n";
    csv.close();
    logmsg(INFO, Form("Efficiency CSV saved to %s", csvFile.Data()));

    // ============================================================ Save global ROOT
    TString outfile = outPlotDir + "hms_tracking_efficiency.root";
    TFile *fout = new TFile(outfile, "RECREATE");
    h_beta_notrack->Write();
    h_etotnorm->Write();
    h_npe->Write();
    g_eff->Write("g_tracking_eff_vs_run");
    fout->Close();
    logmsg(INFO, Form("Global histograms saved to %s", outfile.Data()));

    sw_total.Stop();
    logmsg(INFO, Form("Total time: %.2f s", sw_total.RealTime()));
}