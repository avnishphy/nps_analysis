#include "utils.C"  // logmsg(), trim(), etc.

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <utility>

using namespace std;

// ============================================================
// Runlist reader
// ============================================================
vector<int> readRunList(const string &fname) {
    vector<int> runs;
    ifstream in(fname);
    if (!in.is_open()) {
        logmsg(WARN, "Cannot open runlist: " + fname);
        return runs;
    }
    string line;
    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        try { runs.push_back(stoi(line)); }
        catch (...) { logmsg(WARN, "Skipping invalid runlist line: " + line); }
    }
    return runs;
}

// ============================================================
// Robust Gaussian peak fitter
// ============================================================
pair<double,double> fit_peak(TH1D* h, const string &label) {
    if (!h || h->GetEntries() < 10) {
        int ib = h ? h->GetMaximumBin() : 1;
        return {h ? h->GetBinCenter(ib) : 0.0, 0.0};
    }

    double mean = h->GetMean();
    double rms  = h->GetRMS();
    int ibmax   = h->GetMaximumBin();
    double xpeak = h->GetBinCenter(ibmax);
    double binwidth = h->GetBinWidth(1);
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    cout << "\n==== Fit Debug Info for " << label << " ====\n";
    cout << "Entries: " << h->GetEntries() << "\n";
    cout << "Mean: " << mean << "  RMS: " << rms << "\n";
    cout << "xpeak (max bin center): " << xpeak << "  BinWidth: " << binwidth << "\n";
    cout << "X Range: [" << xmin << ", " << xmax << "]\n";

    if (!isfinite(mean) || rms <= 0)
        mean = xpeak, rms = std::max(binwidth, 0.05 * std::fabs(mean == 0 ? 1.0 : mean));

    double fit_xmin = std::max(xmin, mean - 1 * rms);
    double fit_xmax = std::min(xmax, mean + 1 * rms);
    if (fit_xmax - fit_xmin < 5 * binwidth)
        fit_xmin = std::max(xmin, xpeak - 5 * binwidth),
        fit_xmax = std::min(xmax, xpeak + 5 * binwidth);

    TF1 *f = new TF1(Form("gaus_%s_%p", label.c_str(), (void*)h), "gaus", fit_xmin, fit_xmax);
    f->SetParameters(h->GetMaximum(), mean, rms);
    f->SetParLimits(2, 0.2 * binwidth, xmax - xmin);

    int fitStatus = h->Fit(f, "QRN");
    double peak = (fitStatus == 0) ? f->GetParameter(1) : mean;
    double peak_err = (fitStatus == 0) ? f->GetParError(1) : 0.0;

    delete f;
    return {peak, peak_err};
}

// ============================================================
// HMS Tracking Efficiency + Scaler Rates
// ============================================================
void tracking_efficiency_HMS1(const TString &skimDir_in="output/skimmed/",
                             const TString &outPlotDir_in="output/plots/efficiency/",
                             const TString &runlistFile="config/runlist_x60_4b.txt") {

    TStopwatch sw_total;
    sw_total.Start();
    logmsg(INFO, "=========== HMS Tracking Efficiency ===========");

    TString skimDir = skimDir_in.EndsWith("/") ? skimDir_in : skimDir_in + "/";
    TString outPlotDir = outPlotDir_in.EndsWith("/") ? outPlotDir_in : outPlotDir_in + "/";
    gSystem->mkdir(outPlotDir, true);

    vector<int> runs = readRunList(runlistFile.Data());
    if (runs.empty()) {
        logmsg(ERROR, "No runs found!");
        return;
    }

    // Histograms for sanity checks
    TH1D *h_beta_notrack = new TH1D("h_beta_notrack", "HMS Beta (no track);#beta;Counts", 100, 0.4, 1.4);
    TH1D *h_etotnorm    = new TH1D("h_etotnorm", "E_{tot}/p;E_{tot}/p;Counts", 100, 0.4, 1.4);
    TH1D *h_npe         = new TH1D("h_npe", "HMS Cherenkov NPE;NPE;Counts", 100, 0, 15);

    map<int, pair<double,double>> runEffMap;
    map<int, map<string,double>> runRateMap;
    map<int, map<string,double>> runRateErr;

    double total_Nshould = 0;
    double total_Ndid = 0;

    // ============================================================ Open output ROOT file
    TString outfile = outPlotDir + "hms_tracking_efficiency_extended.root";
    TFile *fout = new TFile(outfile, "RECREATE");
    fout->mkdir("scaler_fits");  // subdirectory for per-run scaler fits

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

        // --- Branch variables ---
        double_t H_hod_betanotrack, H_cal_etotnorm, H_cer_npeSum, H_hod_goodscinhit, H_dc_ntrack;
        double_t s1x_rate, s1y_rate, s2x_rate, s2y_rate;

        T->SetBranchStatus("*", 0);
        T->SetBranchStatus("H.hod.betanotrack", 1);  T->SetBranchAddress("H.hod.betanotrack", &H_hod_betanotrack);
        T->SetBranchStatus("H.cal.etotnorm", 1);     T->SetBranchAddress("H.cal.etotnorm", &H_cal_etotnorm);
        T->SetBranchStatus("H.cer.npeSum", 1);       T->SetBranchAddress("H.cer.npeSum", &H_cer_npeSum);
        T->SetBranchStatus("H.hod.goodscinhit", 1);  T->SetBranchAddress("H.hod.goodscinhit", &H_hod_goodscinhit);
        T->SetBranchStatus("H.dc.ntrack", 1);        T->SetBranchAddress("H.dc.ntrack", &H_dc_ntrack);
        T->SetBranchStatus("H.S1X.scalerRate", 1);   T->SetBranchAddress("H.S1X.scalerRate", &s1x_rate);
        T->SetBranchStatus("H.S1Y.scalerRate", 1);   T->SetBranchAddress("H.S1Y.scalerRate", &s1y_rate);
        T->SetBranchStatus("H.S2X.scalerRate", 1);   T->SetBranchAddress("H.S2X.scalerRate", &s2x_rate);
        T->SetBranchStatus("H.S2Y.scalerRate", 1);   T->SetBranchAddress("H.S2Y.scalerRate", &s2y_rate);

        Long64_t nentries = T->GetEntries();
        Long64_t Nshould = 0, Ndid = 0;

        TH1D *h_s1x = new TH1D("h_s1x","H.S1X.scalerRate;Rate (Hz);Counts",500,0,200e3);
        TH1D *h_s1y = new TH1D("h_s1y","H.S1Y.scalerRate;Rate (Hz);Counts",500,0,200e3);
        TH1D *h_s2x = new TH1D("h_s2x","H.S2X.scalerRate;Rate (Hz);Counts",500,0,200e3);
        TH1D *h_s2y = new TH1D("h_s2y","H.S2Y.scalerRate;Rate (Hz);Counts",500,0,200e3);

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
            h_s1x->Fill(s1x_rate);
            h_s1y->Fill(s1y_rate);
            h_s2x->Fill(s2x_rate);
            h_s2y->Fill(s2y_rate);
        }

        double eff = (Nshould>0) ? (double)Ndid/Nshould : 0.0;
        double eff_unc = (Nshould>0) ? sqrt(eff*(1-eff)/Nshould) : 0.0;
        runEffMap[run] = {eff, eff_unc};

        auto [s1x_peak, s1x_err] = fit_peak(h_s1x, "s1x");
        auto [s1y_peak, s1y_err] = fit_peak(h_s1y, "s1y");
        auto [s2x_peak, s2x_err] = fit_peak(h_s2x, "s2x");
        auto [s2y_peak, s2y_err] = fit_peak(h_s2y, "s2y");

        runRateMap[run] = { {"S1X", s1x_peak}, {"S1Y", s1y_peak}, {"S2X", s2x_peak}, {"S2Y", s2y_peak} };
        runRateErr[run] = { {"S1X", s1x_err}, {"S1Y", s1y_err}, {"S2X", s2x_err}, {"S2Y", s2y_err} };

        logmsg(INFO, Form("Run %d: Eff=%.4f ± %.4f | S1X=%.1f | S1Y=%.1f | S2X=%.1f | S2Y=%.1f",
                          run, eff, eff_unc, s1x_peak, s1y_peak, s2x_peak, s2y_peak));

        total_Nshould += Nshould;
        total_Ndid += Ndid;

        // ------------------- Save fitted scaler histograms -------------------
        TDirectory *scalerDir = fout->GetDirectory("scaler_fits");
        scalerDir->cd();
        TDirectory *thisRunDir = scalerDir->mkdir(Form("run_%d", run));
        thisRunDir->cd();

        auto make_fit = [&](TH1D *h, const char* tag) {
            TF1 *fgaus = new TF1(Form("fgaus_%s_%d", tag, run), "gaus",
                                 h->GetMean()-3*h->GetRMS(), h->GetMean()+3*h->GetRMS());
            fgaus->SetLineColor(kRed);
            h->Fit(fgaus, "QR");
            h->Write(Form("%s_hist", tag));
            fgaus->Write(Form("%s_fit", tag));
            delete fgaus;
        };

        make_fit(h_s1x, "S1X");
        make_fit(h_s1y, "S1Y");
        make_fit(h_s2x, "S2X");
        make_fit(h_s2y, "S2Y");

        delete h_s1x; delete h_s1y; delete h_s2x; delete h_s2y;
        f->Close();
    }

    // ============================================================ Summary plots
    vector<double> run_nums, eff_vals, eff_errs;
    vector<double> s1x_rates, s1y_rates, s2x_rates, s2y_rates;
    for (auto &[run, vals] : runEffMap) {
        run_nums.push_back(run);
        eff_vals.push_back(vals.first);
        eff_errs.push_back(vals.second);
        s1x_rates.push_back(runRateMap[run]["S1X"]);
        s1y_rates.push_back(runRateMap[run]["S1Y"]);
        s2x_rates.push_back(runRateMap[run]["S2X"]);
        s2y_rates.push_back(runRateMap[run]["S2Y"]);
    }

    TCanvas *c_eff = new TCanvas("c_eff","Tracking Efficiency vs Run",800,600);
    TGraphErrors *g_eff = new TGraphErrors(run_nums.size());
    for (size_t i=0; i<run_nums.size(); i++) {
        g_eff->SetPoint(i, run_nums[i], eff_vals[i]);
        g_eff->SetPointError(i, 0, eff_errs[i]);
    }
    g_eff->SetTitle("HMS Tracking Efficiency vs Run;Run Number;Tracking Efficiency");
    g_eff->SetMarkerStyle(20); g_eff->SetMarkerSize(1.0);
    g_eff->SetLineWidth(2); g_eff->SetLineColor(kBlue+1); g_eff->SetMarkerColor(kBlue+2);
    g_eff->Draw("AP");
    c_eff->SetGrid();
    c_eff->SaveAs(outPlotDir + "hms_tracking_eff_vs_run.pdf");

    TCanvas *c_rates = new TCanvas("c_rates","Efficiency vs Scaler Rates",1000,800);
    c_rates->Divide(2,2);
    auto make_graph = [&](int pad, const vector<double>& rate, const char* name, int color){
        c_rates->cd(pad);
        TGraphErrors *g = new TGraphErrors(rate.size());
        for (size_t i=0; i<rate.size(); i++){
            g->SetPoint(i, rate[i], eff_vals[i]);
            g->SetPointError(i, 0, eff_errs[i]);
        }
        g->SetTitle(Form("Tracking Efficiency vs %s; %s Rate (Hz); Efficiency", name, name));
        g->SetMarkerStyle(20);
        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->Draw("AP");
    };
    make_graph(1, s1x_rates, "S1X", kRed+1);
    make_graph(2, s1y_rates, "S1Y", kBlue+1);
    make_graph(3, s2x_rates, "S2X", kGreen+2);
    make_graph(4, s2y_rates, "S2Y", kMagenta+1);
    c_rates->SaveAs(outPlotDir + "hms_tracking_eff_vs_scalerRates.pdf");

    // ============================================================ CSV output
    TString csvFile = outPlotDir + "hms_tracking_efficiency.csv";
    ofstream csv(csvFile.Data());
    csv << "Run,Efficiency,Uncertainty,"
        << "S1X_rate,S1X_err,S1Y_rate,S1Y_err,"
        << "S2X_rate,S2X_err,S2Y_rate,S2Y_err\n";
    for (auto &[run, vals] : runEffMap)
        csv << run << "," << vals.first << "," << vals.second << ","
            << runRateMap[run]["S1X"] << "," << runRateErr[run]["S1X"] << ","
            << runRateMap[run]["S1Y"] << "," << runRateErr[run]["S1Y"] << ","
            << runRateMap[run]["S2X"] << "," << runRateErr[run]["S2X"] << ","
            << runRateMap[run]["S2Y"] << "," << runRateErr[run]["S2Y"] << "\n";
    csv.close();
    logmsg(INFO, Form("CSV saved to %s", csvFile.Data()));

    // ============================================================ Write summary histograms
    fout->cd();
    h_beta_notrack->Write();
    h_etotnorm->Write();
    h_npe->Write();
    g_eff->Write("g_eff_vs_run");
    c_eff->Write("c_eff_vs_run");
    c_rates->Write("c_eff_vs_scalerRates");

    fout->Close();

    logmsg(INFO, Form("Extended output saved to %s", outfile.Data()));
    sw_total.Stop();
    logmsg(INFO, Form("Total time: %.2f s", sw_total.RealTime()));
}
