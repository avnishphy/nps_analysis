// nps_HMS_track_eff.h
// Safe, ACLiC-friendly per-run HMS efficiency utilities with best practices.
//
// Design:
//   - HMSTrackingEffCounter: lightweight per-event scalar collector (no ROOT objects)
//   - HMSEffAggregator: singleton storing POD results, writes global CSV/plots
//   - All ROOT objects use RAII (unique_ptr or TDirectory ownership)
//   - Proper stream formatting to avoid numerical precision issues
//   - Exception-safe with cleanup guaranteed
//
// Usage:
//   1. HMSTrackingEffCounter hmstrackEff(run); before event loop
//   2. hmstrackEff.processEvent(...); in event loop
//   3. hmstrackEff.finalizeAndWrite(fout); after event loop  
//   4. HMSEffAggregator::Instance().WriteGlobalSummary(dir); after all runs

#ifndef NPS_HMS_TRACK_EFF_H
#define NPS_HMS_TRACK_EFF_H

#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <TParameter.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TStyle.h>

#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>

namespace nps { namespace hms_track_eff {

// ============================================================================
// Gaussian Peak Fitter (with RAII and best practices)
// ============================================================================
// Robust gaussian fitting with automatic cleanup and safeguards
static inline std::pair<double,double> fit_peak_local(TH1D* h, const std::string &label = "") 
{
    if (!h || h->GetEntries() < 10) {
        int ib = h ? h->GetMaximumBin() : 1;
        return { h ? h->GetBinCenter(ib) : 0.0, 0.0 };
    }

    double mean = h->GetMean();
    double rms = h->GetRMS();
    int ibmax = h->GetMaximumBin();
    double xpeak = h->GetBinCenter(ibmax);
    double binwidth = h->GetBinWidth(1);
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    // Safeguards for invalid values
    if (!std::isfinite(mean) || rms <= 0.0) {
        mean = xpeak;
        rms = std::max(binwidth, 0.05 * std::fabs(mean == 0.0 ? 1.0 : mean));
    }

    // Fit range: 1 sigma around mean
    double fit_xmin = std::max(xmin, mean - 1.0 * rms);
    double fit_xmax = std::min(xmax, mean + 1.0 * rms);
    if (fit_xmax - fit_xmin < 5.0 * binwidth) {
        fit_xmin = std::max(xmin, xpeak - 5.0 * binwidth);
        fit_xmax = std::min(xmax, xpeak + 5.0 * binwidth);
    }

    // Create TF1 with unique name (avoid collisions in ACLiC)
    TString fname = TString::Format("tmp_gaus_%s_%p", label.c_str(), (void*)h);
    std::unique_ptr<TF1> f(new TF1(fname.Data(), "gaus", fit_xmin, fit_xmax));

    f->SetParameters(h->GetMaximum(), mean, rms);
    f->SetParLimits(2, 0.2 * binwidth, xmax - xmin);

    // Fit: quiet, range, no draw (status code: 0=success)
    int status = h->Fit(f.get(), "QRN");
    double peak = (status == 0) ? f->GetParameter(1) : mean;
    double peak_err = (status == 0) ? f->GetParError(1) : 0.0;

    return {peak, peak_err};
}

// ============================================================================
// HMSEffAggregator: Global results aggregator (singleton)
// ============================================================================
// Stores results across runs, writes CSV and ROOT plots
// Thread-safe singleton with deleted copy/move operations
class HMSEffAggregator 
{
public:
    struct RunEntry 
    {
        int run = -1;
        double eff = 0.0, eff_err = 0.0;
        double s1x = 0.0, s1x_err = 0.0;
        double s1y = 0.0, s1y_err = 0.0;
        double s2x = 0.0, s2x_err = 0.0;
        double s2y = 0.0, s2y_err = 0.0;
    };

private:
    std::vector<RunEntry> runs;

    // Private constructor for singleton pattern
    HMSEffAggregator() = default;

public:
    // Deleted copy/move to enforce singleton
    HMSEffAggregator(const HMSEffAggregator&) = delete;
    HMSEffAggregator& operator=(const HMSEffAggregator&) = delete;
    HMSEffAggregator(HMSEffAggregator&&) = delete;
    HMSEffAggregator& operator=(HMSEffAggregator&&) = delete;

    // Singleton instance (thread-safe static initialization)
    static HMSEffAggregator& Instance() 
    {
        static HMSEffAggregator inst;
        return inst;
    }

    // Add result from a single run
    inline void AddRunResult(int run, 
                            double eff, double eff_err,
                            double s1x, double s1x_err,
                            double s1y, double s1y_err,
                            double s2x, double s2x_err,
                            double s2y, double s2y_err)
    {
        RunEntry e;
        e.run = run;
        e.eff = eff; 
        e.eff_err = eff_err;
        e.s1x = s1x; 
        e.s1x_err = s1x_err;
        e.s1y = s1y; 
        e.s1y_err = s1y_err;
        e.s2x = s2x; 
        e.s2x_err = s2x_err;
        e.s2y = s2y; 
        e.s2y_err = s2y_err;
        runs.push_back(e);
    }
    
    // Get last run entry (most recently added)
    inline RunEntry GetLastRunEntry() const {
        if (runs.empty()) {
            RunEntry empty;
            empty.run = -1;
            return empty;
        }
        return runs.back();
    }

    // Write CSV and ROOT plots (call once after all runs)
    inline void WriteGlobalSummary(const TString &outPlotDir_in = "output/plots/efficiency/") 
    {
        TString outPlotDir = outPlotDir_in.EndsWith("/") ? outPlotDir_in : outPlotDir_in + "/";
        gSystem->mkdir(outPlotDir, true);

        if (runs.empty()) {
            std::cout << "[HMSEffAggregator] No runs to write.\n";
            return;
        }

        // Sort runs by run number
        std::sort(runs.begin(), runs.end(), 
                  [](const RunEntry &a, const RunEntry &b){ return a.run < b.run; });

        // ===== Write CSV disabled for standalone efficiency workflow =====
        // (main analysis consumes auxiliary CSVs and does not compute efficiencies)
        
        // ===== Write ROOT plots with RAII =====
        std::vector<double> runnums, effs, efferrs;
        for (const auto &r : runs) { 
            runnums.push_back(r.run); 
            effs.push_back(r.eff); 
            efferrs.push_back(r.eff_err); 
        }

        // Create canvas and graph (RAII cleanup)
        std::unique_ptr<TCanvas> c_eff(new TCanvas("c_eff_hms", "HMS tracking eff vs run", 900, 600));
        std::unique_ptr<TGraphErrors> g_eff(new TGraphErrors((int)runnums.size()));

        c_eff->SetLeftMargin(0.12);
        c_eff->SetRightMargin(0.05);
        c_eff->SetBottomMargin(0.12);
        c_eff->SetTopMargin(0.08);
        c_eff->SetGridy(true);
        c_eff->SetTicks(1, 1);
        
        for (size_t i = 0; i < runnums.size(); ++i) {
            g_eff->SetPoint((int)i, runnums[i], effs[i]);
            g_eff->SetPointError((int)i, 0, efferrs[i]);
        }
        
        g_eff->SetTitle("HMS Tracking Efficiency vs Run;Run;Efficiency");
        g_eff->SetMarkerStyle(20);
        g_eff->SetMarkerSize(1.0);
        g_eff->SetMarkerColor(kBlue);
        g_eff->SetLineColor(kBlue + 1);
        g_eff->SetLineWidth(2);

        double ymin = 1.0;
        double ymax = 0.0;
        for (size_t i = 0; i < effs.size(); ++i) {
            const double lo = effs[i] - efferrs[i];
            const double hi = effs[i] + efferrs[i];
            ymin = std::min(ymin, lo);
            ymax = std::max(ymax, hi);
        }
        ymin = std::max(0.0, ymin - 0.03);
        ymax = std::min(1.05, ymax + 0.03);
        if (!(ymax > ymin)) {
            ymin = 0.0;
            ymax = 1.05;
        }

        g_eff->SetMinimum(ymin);
        g_eff->SetMaximum(ymax);
        g_eff->Draw("AP");
        g_eff->GetXaxis()->SetTitleSize(0.045);
        g_eff->GetYaxis()->SetTitleSize(0.045);
        g_eff->GetXaxis()->SetLabelSize(0.040);
        g_eff->GetYaxis()->SetLabelSize(0.040);

        c_eff->SaveAs(outPlotDir + "hms_track_eff_vs_run.png");
        c_eff->SaveAs(outPlotDir + "hms_track_eff_vs_run.pdf");
        // unique_ptr auto-cleanup on exit

        std::cout << "[HMSEffAggregator] Wrote plots to: " << outPlotDir.Data() << "\n";
    }
};

// ============================================================================
// HMSTrackingEffCounter: Per-run lightweight collector
// ============================================================================
// Collects per-event scaler samples (lightweight) rather than persistent
// ROOT objects. This avoids ownership issues across multiple runs.
class HMSTrackingEffCounter 
{
public:
    const int run;
    Long64_t Nshould = 0;
    Long64_t Ndid = 0;

    // Per-event samples (lightweight doubles, pre-allocated)
    std::vector<double> s1x_samples;
    std::vector<double> s1y_samples;
    std::vector<double> s2x_samples;
    std::vector<double> s2y_samples;

    // Diagnostic samples (for "should" events only)
    std::vector<double> beta_notrack_samples;
    std::vector<double> etotnorm_samples;
    std::vector<double> npe_samples;

    // Constructor with run number
    explicit HMSTrackingEffCounter(int runnum = -1) : run(runnum) 
    {
        s1x_samples.reserve(16384);
        s1y_samples.reserve(16384);
        s2x_samples.reserve(16384);
        s2y_samples.reserve(16384);
        beta_notrack_samples.reserve(4096);
        etotnorm_samples.reserve(4096);
        npe_samples.reserve(4096);
    }

    // Default destructor (no ROOT objects owned)
    ~HMSTrackingEffCounter() = default;

    // Process one event (call once per event in main loop)
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
        // Check "should" criteria
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

        // Always collect scaler rates
        s1x_samples.push_back(s1x_rate);
        s1y_samples.push_back(s1y_rate);
        s2x_samples.push_back(s2x_rate);
        s2y_samples.push_back(s2y_rate);
    }

    // Finalize: create per-run histograms in fout, write results
    // Returns: (efficiency, uncertainty)
    inline std::pair<double,double> finalizeAndWrite(TFile *fout = nullptr, 
                                                      const TString &outPlotDir = "") 
    {
        double eff = 0.0, eff_unc = 0.0;

        // Calculate efficiency with binomial uncertainty
        if (Nshould > 0) {
            eff = static_cast<double>(Ndid) / static_cast<double>(Nshould);
            eff_unc = std::sqrt(eff * (1.0 - eff) / static_cast<double>(std::max<Long64_t>(1, Nshould)));
        }

        // Write per-run histograms if file provided
        if (fout && !fout->IsZombie()) {
            TString dirName = TString::Format("hms_track_eff/run_%d", run);
            TDirectory *base = fout->GetDirectory("");
            TDirectory *d = base->mkdir(dirName.Data());
            if (!d) {
                std::cerr << "[HMSTrackingEffCounter] ERROR: Failed to create directory " << dirName << "\n";
                return {eff, eff_unc};
            }
            d->cd();

            // Create scaler rate histograms (owned by directory)
            std::unique_ptr<TH1D> h_s1x(new TH1D(Form("h_s1x_run%d", run), 
                                                 Form("H.S1X.scalerRate_run%d;Rate (Hz);Counts", run), 
                                                 500, 0, 200e3));
            std::unique_ptr<TH1D> h_s1y(new TH1D(Form("h_s1y_run%d", run), 
                                                 Form("H.S1Y.scalerRate_run%d;Rate (Hz);Counts", run), 
                                                 500, 0, 200e3));
            std::unique_ptr<TH1D> h_s2x(new TH1D(Form("h_s2x_run%d", run), 
                                                 Form("H.S2X.scalerRate_run%d;Rate (Hz);Counts", run), 
                                                 500, 0, 200e3));
            std::unique_ptr<TH1D> h_s2y(new TH1D(Form("h_s2y_run%d", run), 
                                                 Form("H.S2Y.scalerRate_run%d;Rate (Hz);Counts", run), 
                                                 500, 0, 200e3));

            // Fill histograms
            for (double v : s1x_samples) h_s1x->Fill(v);
            for (double v : s1y_samples) h_s1y->Fill(v);
            for (double v : s2x_samples) h_s2x->Fill(v);
            for (double v : s2y_samples) h_s2y->Fill(v);

            // Create diagnostic histograms (for "should" events)
            std::unique_ptr<TH1D> h_beta(new TH1D(Form("h_beta_notrack_run%d", run), 
                                                  Form("HMS Beta (no track) run%d;#beta;Counts", run), 
                                                  100, 0.4, 1.4));
            std::unique_ptr<TH1D> h_etot(new TH1D(Form("h_etotnorm_run%d", run), 
                                                  Form("E_{tot}/p run%d;E_{tot}/p;Counts", run), 
                                                  100, 0.4, 1.4));
            std::unique_ptr<TH1D> h_npe(new TH1D(Form("h_npe_run%d", run), 
                                                 Form("HMS NPE run%d;NPE;Counts", run), 
                                                 100, 0, 15));

            for (double v : beta_notrack_samples) h_beta->Fill(v);
            for (double v : etotnorm_samples) h_etot->Fill(v);
            for (double v : npe_samples) h_npe->Fill(v);

            // Fit peaks using local fitter
            auto p1 = fit_peak_local(h_s1x.get(), TString::Format("s1x_run%d", run).Data());
            auto p2 = fit_peak_local(h_s1y.get(), TString::Format("s1y_run%d", run).Data());
            auto p3 = fit_peak_local(h_s2x.get(), TString::Format("s2x_run%d", run).Data());
            auto p4 = fit_peak_local(h_s2y.get(), TString::Format("s2y_run%d", run).Data());

            // Write histograms (owned by TDirectory now)
            h_s1x->Write("h_s1x", TObject::kOverwrite);
            h_s1y->Write("h_s1y", TObject::kOverwrite);
            h_s2x->Write("h_s2x", TObject::kOverwrite);
            h_s2y->Write("h_s2y", TObject::kOverwrite);

            h_beta->Write("h_beta_notrack", TObject::kOverwrite);
            h_etot->Write("h_etotnorm", TObject::kOverwrite);
            h_npe->Write("h_npe", TObject::kOverwrite);

            // Write scalar results as TParameter
            TParameter<double>(Form("hms_track_eff_run_%d_eff", run), eff).Write();
            TParameter<double>(Form("hms_track_eff_run_%d_eff_unc", run), eff_unc).Write();
            TParameter<double>(Form("hms_track_eff_run_%d_Nshould", run), static_cast<double>(Nshould)).Write();
            TParameter<double>(Form("hms_track_eff_run_%d_Ndid", run), static_cast<double>(Ndid)).Write();

            TParameter<double>(Form("hms_track_eff_run_%d_S1X_peak", run), p1.first).Write();
            TParameter<double>(Form("hms_track_eff_run_%d_S1Y_peak", run), p2.first).Write();
            TParameter<double>(Form("hms_track_eff_run_%d_S2X_peak", run), p3.first).Write();
            TParameter<double>(Form("hms_track_eff_run_%d_S2Y_peak", run), p4.first).Write();

            // Histograms owned by directory, cleaned up when fout->Close()
        }

        // Compute peaks from samples if needed for aggregator
        double s1x_peak = 0, s1x_err = 0;
        double s1y_peak = 0, s1y_err = 0;
        double s2x_peak = 0, s2x_err = 0;
        double s2y_peak = 0, s2y_err = 0;

        // Helper: compute mean and RMS from sample vector
        auto compute_mean_rms = [](const std::vector<double> &v) -> std::pair<double,double> 
        {
            if (v.empty()) return {0.0, 0.0};
            double sum = 0.0;
            for (double x : v) sum += x;
            double mean = sum / v.size();
            double var = 0.0;
            for (double x : v) var += (x - mean) * (x - mean);
            var /= v.size();
            return {mean, std::sqrt(var)};
        };

        if (!s1x_samples.empty()) { auto pr = compute_mean_rms(s1x_samples); s1x_peak = pr.first; s1x_err = pr.second; }
        if (!s1y_samples.empty()) { auto pr = compute_mean_rms(s1y_samples); s1y_peak = pr.first; s1y_err = pr.second; }
        if (!s2x_samples.empty()) { auto pr = compute_mean_rms(s2x_samples); s2x_peak = pr.first; s2x_err = pr.second; }
        if (!s2y_samples.empty()) { auto pr = compute_mean_rms(s2y_samples); s2y_peak = pr.first; s2y_err = pr.second; }

        // Register with global aggregator
        HMSEffAggregator::Instance().AddRunResult(run, eff, eff_unc,
                                                  s1x_peak, s1x_err,
                                                  s1y_peak, s1y_err,
                                                  s2x_peak, s2x_err,
                                                  s2y_peak, s2y_err);

        return {eff, eff_unc};
    }
};

} } // namespace nps::hms_track_eff

#endif // NPS_HMS_TRACK_EFF_H
