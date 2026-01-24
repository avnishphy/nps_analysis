// nps_HMS_track_eff.h
#ifndef NPS_HMS_TRACK_EFF_H
#define NPS_HMS_TRACK_EFF_H

// nps_HMS_track_eff.h
// Safe, ACLiC-friendly per-run HMS efficiency utilities.
// - HMSTrackingEffCounter collects per-event scalers/summaries (no persistent ROOT objects).
// - finalizeAndWrite(fout) creates per-run histograms inside fout and writes scalars.
// - HMSEffAggregator stores only POD results and writes a global CSV + plots at the end.

// Note: include this header from your macro (nps_analysis.C) and instantiate HMSTrackingEffCounter
// inside each run before the event loop, call processEvent(...) per event, then call
// finalizeAndWrite(fout) after the event loop. After all runs call
// HMSEffAggregator::Instance().WriteGlobalSummary(outPlotDir).

#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TDirectory.h>
#include <TSystem.h>

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <utility>
#include <algorithm>

namespace nps { namespace hms_track_eff {

// ------------------------------------------------------------------
// Local helper: robust gaussian peak fitter for a TH1D
// Uses a unique TF1 name to avoid name collisions in ACLiC.
// ------------------------------------------------------------------
static inline std::pair<double,double> fit_peak_local(TH1D* h, const std::string &label = "") {
    if (!h || h->GetEntries() < 10) {
        int ib = h ? h->GetMaximumBin() : 1;
        return { h ? h->GetBinCenter(ib) : 0.0, 0.0 };
    }

    double mean = h->GetMean();
    double rms  = h->GetRMS();
    int ibmax   = h->GetMaximumBin();
    double xpeak = h->GetBinCenter(ibmax);
    double binwidth = h->GetBinWidth(1);
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    if (!std::isfinite(mean) || rms <= 0) {
        mean = xpeak;
        rms = std::max(binwidth, 0.05 * std::fabs(mean == 0 ? 1.0 : mean));
    }

    double fit_xmin = std::max(xmin, mean - 1 * rms);
    double fit_xmax = std::min(xmax, mean + 1 * rms);
    if (fit_xmax - fit_xmin < 5 * binwidth) {
        fit_xmin = std::max(xmin, xpeak - 5 * binwidth);
        fit_xmax = std::min(xmax, xpeak + 5 * binwidth);
    }

    TString fname = TString::Format("tmp_gaus_%s_%p", label.c_str(), (void*)h);
    TF1 *f = new TF1(fname.Data(), "gaus", fit_xmin, fit_xmax);
    f->SetParameters(h->GetMaximum(), mean, rms);
    f->SetParLimits(2, 0.2 * binwidth, xmax - xmin);

    int status = h->Fit(f, "QRN"); // quiet, range, no draw
    double peak = (status == 0) ? f->GetParameter(1) : mean;
    double peak_err = (status == 0) ? f->GetParError(1) : 0.0;

    delete f;
    return {peak, peak_err};
}

// ------------------------------------------------------------------
// HMSEffAggregator: stores results across runs and writes summary plots/CSV
// Singleton; stores only POD data (no ROOT pointers stored).
// ------------------------------------------------------------------
class HMSEffAggregator {
public:
    struct RunEntry {
        int run = -1;
        double eff = 0.0, eff_err = 0.0;
        double s1x = 0.0, s1x_err = 0.0;
        double s1y = 0.0, s1y_err = 0.0;
        double s2x = 0.0, s2x_err = 0.0;
        double s2y = 0.0, s2y_err = 0.0;
    };

private:
    std::vector<RunEntry> runs;
    HMSEffAggregator() {}
public:
    static HMSEffAggregator & Instance() {
        static HMSEffAggregator inst;
        return inst;
    }

    inline void AddRunResult(int run, double eff, double eff_err,
                             double s1x, double s1x_err,
                             double s1y, double s1y_err,
                             double s2x, double s2x_err,
                             double s2y, double s2y_err)
    {
        RunEntry e;
        e.run = run;
        e.eff = eff; e.eff_err = eff_err;
        e.s1x = s1x; e.s1x_err = s1x_err;
        e.s1y = s1y; e.s1y_err = s1y_err;
        e.s2x = s2x; e.s2x_err = s2x_err;
        e.s2y = s2y; e.s2y_err = s2y_err;
        runs.push_back(e);
    }

    // Write CSV and simple ROOT plots to outPlotDir (call once after all runs).
    inline void WriteGlobalSummary(const TString &outPlotDir_in = "output/plots/efficiency/") {
        TString outPlotDir = outPlotDir_in.EndsWith("/") ? outPlotDir_in : outPlotDir_in + "/";
        gSystem->mkdir(outPlotDir, true);

        if (runs.empty()) {
            std::cout << "[HMSEffAggregator] No runs to write.\n";
            return;
        }

        std::sort(runs.begin(), runs.end(), [](const RunEntry &a, const RunEntry &b){ return a.run < b.run; });

        // CSV
        TString csvFile = outPlotDir + "hms_tracking_efficiency_summary.csv";
        std::ofstream csv(csvFile.Data());
        csv << "Run,Eff,Eff_err,S1X_peak,S1X_err,S1Y_peak,S1Y_err,S2X_peak,S2X_err,S2Y_peak,S2Y_err\n";
        for (const auto &r : runs) {
            csv << r.run << "," << std::setprecision(6) << r.eff << "," << r.eff_err << ","
                << r.s1x << "," << r.s1x_err << ","
                << r.s1y << "," << r.s1y_err << ","
                << r.s2x << "," << r.s2x_err << ","
                << r.s2y << "," << r.s2y_err << "\n";
        }
        csv.close();
        std::cout << "[HMSEffAggregator] Wrote CSV: " << csvFile.Data() << "\n";

        // Simple ROOT graph: Eff vs Run (created once at the end)
        std::vector<double> runnums, effs, efferrs;
        for (const auto &r : runs) { runnums.push_back(r.run); effs.push_back(r.eff); efferrs.push_back(r.eff_err); }

        TCanvas *c_eff = new TCanvas("c_eff_hms", "HMS tracking eff vs run", 900,600);
        TGraphErrors *g_eff = new TGraphErrors((int)runnums.size());
        for (size_t i=0;i<runnums.size();++i) {
            g_eff->SetPoint((int)i, runnums[i], effs[i]);
            g_eff->SetPointError((int)i, 0, efferrs[i]);
        }
        g_eff->SetTitle("HMS Tracking Efficiency vs Run;Run;Efficiency");
        g_eff->SetMarkerStyle(20);
        g_eff->Draw("AP");
        c_eff->SaveAs(outPlotDir + "hms_track_eff_vs_run.png");
        c_eff->SaveAs(outPlotDir + "hms_track_eff_vs_run.pdf");

        // Leave ROOT objects to be owned by ROOT. Do not store them in the aggregator.
        std::cout << "[HMSEffAggregator] Wrote global plots to: " << outPlotDir.Data() << "\n";
    }
};

// ------------------------------------------------------------------
// HMSTrackingEffCounter
// Per-run collector: stores per-event scalar samples (vectors) not histograms.
// This avoids long-lived ROOT ownership issues when running multiple runs.
// ------------------------------------------------------------------
class HMSTrackingEffCounter {
public:
    const int run;
    Long64_t Nshould = 0;
    Long64_t Ndid = 0;

    // store raw samples per-event (lightweight doubles)
    std::vector<double> s1x_samples;
    std::vector<double> s1y_samples;
    std::vector<double> s2x_samples;
    std::vector<double> s2y_samples;

    // store diagnostics only for events that satisfy 'should' (like previous logic)
    std::vector<double> beta_notrack_samples;
    std::vector<double> etotnorm_samples;
    std::vector<double> npe_samples;

    HMSTrackingEffCounter(int runnum = -1) : run(runnum) {
        s1x_samples.reserve(16384);
        s1y_samples.reserve(16384);
        s2x_samples.reserve(16384);
        s2y_samples.reserve(16384);
        beta_notrack_samples.reserve(4096);
        etotnorm_samples.reserve(4096);
        npe_samples.reserve(4096);
    }

    // processEvent: call once per event; pass branch values already read in your macro
    inline void processEvent(double H_hod_betanotrack,
                             double H_cal_etotnorm,
                             double H_cer_npeSum,
                             double H_hod_goodscinhit,
                             double H_dc_ntrack,
                             double s1x_rate = 0.0,
                             double s1y_rate = 0.0,
                             double s2x_rate = 0.0,
                             double s2y_rate = 0.0)
    {
        bool should = (H_hod_goodscinhit == 1 &&
                       H_hod_betanotrack > 0.5 && H_hod_betanotrack < 1.5 &&
                       H_cal_etotnorm > 0.6 &&
                       H_cer_npeSum > 0.6);
        if (should) {
            ++Nshould;
            beta_notrack_samples.push_back(H_hod_betanotrack);
            etotnorm_samples.push_back(H_cal_etotnorm);
            npe_samples.push_back(H_cer_npeSum);
            if (H_dc_ntrack > 0) ++Ndid;
        }
        s1x_samples.push_back(s1x_rate);
        s1y_samples.push_back(s1y_rate);
        s2x_samples.push_back(s2x_rate);
        s2y_samples.push_back(s2y_rate);
    }

    // finalizeAndWrite: create per-run histograms inside fout and write scalars; return {eff, eff_unc}
    inline std::pair<double,double> finalizeAndWrite(TFile *fout = nullptr, const TString &outPlotDir = "") {
        double eff = 0.0, eff_unc = 0.0;
        if (Nshould > 0) {
            eff = static_cast<double>(Ndid) / static_cast<double>(Nshould);
            eff_unc = std::sqrt(eff * (1.0 - eff) / static_cast<double>(std::max<Long64_t>(1, Nshould)));
        }

        // Prepare per-run output directory in the provided fout
        if (fout && !fout->IsZombie()) {
            TString dirName = TString::Format("hms_track_eff/run_%d", run);
            TDirectory *base = fout->GetDirectory("");
            TDirectory *d = base->mkdir(dirName.Data());
            d->cd();

            // Create & fill scaler histograms (owned by the TFile/directory)
            TH1D *h_s1x = new TH1D(Form("h_s1x_run%d", run), Form("H.S1X.scalerRate_run%d;Rate (Hz);Counts", run), 500, 0, 200e3);
            TH1D *h_s1y = new TH1D(Form("h_s1y_run%d", run), Form("H.S1Y.scalerRate_run%d;Rate (Hz);Counts", run), 500, 0, 200e3);
            TH1D *h_s2x = new TH1D(Form("h_s2x_run%d", run), Form("H.S2X.scalerRate_run%d;Rate (Hz);Counts", run), 500, 0, 200e3);
            TH1D *h_s2y = new TH1D(Form("h_s2y_run%d", run), Form("H.S2Y.scalerRate_run%d;Rate (Hz);Counts", run), 500, 0, 200e3);

            for (double v : s1x_samples) h_s1x->Fill(v);
            for (double v : s1y_samples) h_s1y->Fill(v);
            for (double v : s2x_samples) h_s2x->Fill(v);
            for (double v : s2y_samples) h_s2y->Fill(v);

            // Diagnostic histograms (only for 'should' events)
            TH1D *h_beta = new TH1D(Form("h_beta_notrack_run%d", run), Form("HMS Beta (no track) run%d;#beta;Counts", run), 100, 0.4, 1.4);
            TH1D *h_etot = new TH1D(Form("h_etotnorm_run%d", run), Form("E_{tot}/p run%d;E_{tot}/p;Counts", run), 100, 0.4, 1.4);
            TH1D *h_npe  = new TH1D(Form("h_npe_run%d", run), Form("HMS NPE run%d;NPE;Counts", run), 100, 0, 15);

            for (double v : beta_notrack_samples) h_beta->Fill(v);
            for (double v : etotnorm_samples) h_etot->Fill(v);
            for (double v : npe_samples) h_npe->Fill(v);

            // Fit peaks (use the local fitter)
            auto p1 = fit_peak_local(h_s1x, TString::Format("s1x_run%d", run).Data());
            auto p2 = fit_peak_local(h_s1y, TString::Format("s1y_run%d", run).Data());
            auto p3 = fit_peak_local(h_s2x, TString::Format("s2x_run%d", run).Data());
            auto p4 = fit_peak_local(h_s2y, TString::Format("s2y_run%d", run).Data());

            // Write histograms & scalar TParameters
            h_s1x->Write("h_s1x", TObject::kOverwrite);
            h_s1y->Write("h_s1y", TObject::kOverwrite);
            h_s2x->Write("h_s2x", TObject::kOverwrite);
            h_s2y->Write("h_s2y", TObject::kOverwrite);

            h_beta->Write("h_beta_notrack", TObject::kOverwrite);
            h_etot->Write("h_etotnorm", TObject::kOverwrite);
            h_npe->Write("h_npe", TObject::kOverwrite);

            // scalars as TParameter
            TParameter<double> (Form("hms_track_eff_run_%d_eff", run), eff).Write();
            TParameter<double> (Form("hms_track_eff_run_%d_eff_unc", run), eff_unc).Write();
            TParameter<double> (Form("hms_track_eff_run_%d_Nshould", run), static_cast<double>(Nshould)).Write();
            TParameter<double> (Form("hms_track_eff_run_%d_Ndid", run), static_cast<double>(Ndid)).Write();

            TParameter<double> (Form("hms_track_eff_run_%d_S1X_peak", run), p1.first).Write();
            TParameter<double> (Form("hms_track_eff_run_%d_S1Y_peak", run), p2.first).Write();
            TParameter<double> (Form("hms_track_eff_run_%d_S2X_peak", run), p3.first).Write();
            TParameter<double> (Form("hms_track_eff_run_%d_S2Y_peak", run), p4.first).Write();

            // Note: do NOT delete these histograms; they are owned by the TDirectory (fout) now.
            // They will be cleaned up when fout->Close() is called in the outer code.
        } // end if fout valid

        // Register summary in aggregator for global plots
        // We pass peaks computed if fout existed, or zero otherwise.
        // If fout not provided, estimate peaks by simple statistics from the vectors:
        double s1x_peak=0, s1x_err=0, s1y_peak=0, s1y_err=0, s2x_peak=0, s2x_err=0, s2y_peak=0, s2y_err=0;
        // If samples exist and we did not compute above (no fout), compute rough mean + RMS
        auto compute_mean_rms = [](const std::vector<double> &v)->std::pair<double,double>{
            if (v.empty()) return {0.0, 0.0};
            double sum = 0.0;
            for (double x : v) sum += x;
            double mean = sum / v.size();
            double var = 0.0;
            for (double x : v) var += (x-mean)*(x-mean);
            var /= v.size();
            return {mean, std::sqrt(var)};
        };
        if (!s1x_samples.empty()) { auto pr = compute_mean_rms(s1x_samples); s1x_peak = pr.first; s1x_err = pr.second; }
        if (!s1y_samples.empty()) { auto pr = compute_mean_rms(s1y_samples); s1y_peak = pr.first; s1y_err = pr.second; }
        if (!s2x_samples.empty()) { auto pr = compute_mean_rms(s2x_samples); s2x_peak = pr.first; s2x_err = pr.second; }
        if (!s2y_samples.empty()) { auto pr = compute_mean_rms(s2y_samples); s2y_peak = pr.first; s2y_err = pr.second; }

        HMSEffAggregator::Instance().AddRunResult(run, eff, eff_unc,
                                                  s1x_peak, s1x_err,
                                                  s1y_peak, s1y_err,
                                                  s2x_peak, s2x_err,
                                                  s2y_peak, s2y_err);

        return {eff, eff_unc};
    }

    // no destructor needed; we do not own ROOT objects
    ~HMSTrackingEffCounter() = default;
};

} } // namespace nps::hms_track_eff

#endif // NPS_HMS_TRACK_EFF_H
