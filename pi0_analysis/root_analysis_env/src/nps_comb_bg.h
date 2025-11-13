// nps_comb_bg_plots.h
// Enhanced combinatorial BG fit + subtraction with rich diagnostic plotting.
// Usage:
//   #include "nps_comb_bg_plots.h"
//   nps::BGSubtractionResult res = nps::FitCombinatorialBGAndSubtract(h_coin_bgsub, "output/plots", run, 2);
// Returns: BGSubtractionResult containing the bg-subtracted histogram (caller owns h_final) and fit summary.

#ifndef NPS_COMB_BG_H
#define NPS_COMB_BG_H

#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TMatrixDSym.h>
#include <TFitResultPtr.h>
#include <TFile.h>
#include <TLatex.h>
#include <TMath.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TAttFill.h>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>

namespace nps {

struct BGSubtractionResult {
    TH1D* h_final = nullptr;
    double chi2_ndf = -1.0;
    double mu_MeV = 0.0;
    double sigma_MeV = 0.0;
    double signal_counts = 0.0;
};

// Build a TGraphErrors from histogram using only bins inside the two sideband windows
inline TGraphErrors* MakeSidebandGraph(const TH1D* h, double left_lo, double left_hi, double right_lo, double right_hi) {
    if (!h) return nullptr;
    int nb = h->GetNbinsX();
    std::vector<double> xs, ys, ex, ey;
    xs.reserve(nb); ys.reserve(nb); ex.reserve(nb); ey.reserve(nb);
    for (int b = 1; b <= nb; ++b) {
        double x = h->GetXaxis()->GetBinCenter(b);
        if ((x >= left_lo && x <= left_hi) || (x >= right_lo && x <= right_hi)) {
            double y = h->GetBinContent(b);
            double err = h->GetBinError(b);
            if (err <= 0.0) err = (y > 0.0 ? std::sqrt(y) : 1.0);
            xs.push_back(x); ys.push_back(y); ex.push_back(0.0); ey.push_back(err);
        }
    }
    if (xs.empty()) return nullptr;
    TGraphErrors* g = new TGraphErrors((int)xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {
        g->SetPoint((int)i, xs[i], ys[i]);
        g->SetPointError((int)i, ex[i], ey[i]);
    }
    return g;
}

// Helper: build polynomial formula string for given order (e.g. order=2 -> "[0] + [1]*x + [2]*x*x")
inline std::string PolyFormula(int order) {
    std::ostringstream os;
    for (int i = 0; i <= order; ++i) {
        if (i > 0) os << " + ";
        os << "[" << i << "]";
        if (i > 0) {
            os << "*";
            for (int p = 0; p < i; ++p) {
                if (p > 0) os << "*";
                os << "x";
            }
        }
    }
    return os.str();
}

// Compute background error at x using covariance matrix cov (parameters are polynomial coefficients par[0..order])
// For polynomial f(x) = sum_i p_i * x^i, gradient wrt p_i is x^i.
inline double PolyErrorAtX(const TMatrixDSym &cov, double x) {
    int n = cov.GetNrows();
    if (n == 0) return 0.0;
    std::vector<double> v(n);
    for (int i = 0; i < n; ++i) v[i] = std::pow(x, i);
    double var = 0.0;
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) var += v[i] * cov(i,j) * v[j];
    return (var > 0.0) ? std::sqrt(var) : 0.0;
}

// Nicely format a double with precision for TLatex
inline std::string fmt(double x, int prec=3) {
    std::ostringstream os; os<<std::fixed<<std::setprecision(prec)<<x; return os.str();
}

// Main helper: enhanced plotting
// NOTE: returns BGSubtractionResult. Caller takes ownership of result.h_final and should delete it when done.
inline BGSubtractionResult FitCombinatorialBGAndSubtract(TH1D *h_coin_bgsub,
                                                         const char *outDir = "",
                                                         int run = -1,
                                                         int poly_order = 2,
                                                         double left_lo = 0.01, double left_hi = 0.1,
                                                         double right_lo = 0.15, double right_hi = 0.40,
                                                         bool draw = true)
{
    BGSubtractionResult result;

    if (!h_coin_bgsub) {
        std::cerr << "[nps::FitCombinatorialBGAndSubtract] null histogram\n";
        return result;
    }

    // guard Sumw2
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    if (h_coin_bgsub->GetSumw2N() == 0) h_coin_bgsub->Sumw2();
    #else
    h_coin_bgsub->Sumw2();
    #endif

    if (poly_order < 0) poly_order = 2;

    // Set a clean, modern style for the canvas
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Build sideband graph
    TGraphErrors* g_side = MakeSidebandGraph(h_coin_bgsub, left_lo, left_hi, right_lo, right_hi);
    if (!g_side) {
        std::cerr << "[nps::FitCombinatorialBGAndSubtract] ERROR: no sideband points found for given windows.\n";
        return result;
    }
    g_side->SetMarkerStyle(21);
    g_side->SetMarkerSize(0.9);
    g_side->SetMarkerColor(kBlue+2);

    // Build polynomial TF1 covering full histogram range
    std::string formula = PolyFormula(poly_order);
    double x_min = h_coin_bgsub->GetXaxis()->GetXmin();
    double x_max = h_coin_bgsub->GetXaxis()->GetXmax();
    TF1 *f_bg = new TF1("f_bg", formula.c_str(), x_min, x_max);
    // init parameters: p0 = average sideband height
    double sumy = 0.0; int ny = 0;
    for (int i = 0; i < g_side->GetN(); ++i) {
        double x, y; g_side->GetPoint(i, x, y);
        sumy += y; ny++;
    }
    double avg = (ny > 0) ? (sumy / ny) : 0.0;
    for (int i = 0; i <= poly_order; ++i) f_bg->SetParameter(i, (i==0)? avg : 0.0);

    // Fit sideband TGraphErrors with TF1 and request the fit result ("S") so we can extract covariance.
    TFitResultPtr fitres = g_side->Fit(f_bg, "RQS"); // R:range, Q:quiet, S:return result
    if (!fitres.Get()) std::cerr << "[nps::FitCombinatorialBGAndSubtract] Warning: fit returned null TFitResultPtr\n";

    // Extract covariance matrix for polynomial params (fallback zeros)
    TMatrixDSym cov_all;
    if (fitres.Get() && fitres->IsValid()) {
        cov_all = fitres->GetCovarianceMatrix();
    } else {
        cov_all.ResizeTo(poly_order+1, poly_order+1);
        for (int i=0;i<=poly_order;++i) for (int j=0;j<=poly_order;++j) cov_all(i,j) = 0.0;
        std::cerr << "[nps::FitCombinatorialBGAndSubtract] Warning: no valid covariance matrix from fit; background errors will be zero.\n";
    }

    // Prefer chi2/ndf from fitres if available
    double chi2 = 0.0; int ndf = 0;
    if (fitres.Get() && fitres->IsValid()) { chi2 = fitres->Chi2(); ndf = fitres->Ndf(); }
    else { chi2 = f_bg->GetChisquare(); ndf = f_bg->GetNDF(); }

    // Print fit params
    std::cout << "[nps::FitCombinatorialBGAndSubtract] run="<<run<<" poly_order="<<poly_order
              <<"  chi2="<<chi2<<" ndf="<<ndf<<" chi2/ndf="<<(ndf>0?chi2/ndf:-1)<<"\n";
    for (int i=0;i<=poly_order;++i) std::cout << "  p["<<i<<"] = "<< f_bg->GetParameter(i) << " ± " << f_bg->GetParError(i) << "\n";

    // Build bg-subtracted histogram
    TString name_final = TString::Format("%s_bgsub_fit_sub_run%d", h_coin_bgsub->GetName(), run);
    TH1D *h_final = (TH1D*)h_coin_bgsub->Clone(name_final);
    h_final->SetDirectory(nullptr);
    h_final->SetTitle(TString::Format("%s (bg fit-subtracted)", h_coin_bgsub->GetTitle()).Data());
    #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
    if (h_final->GetSumw2N() == 0) h_final->Sumw2();
    #else
    h_final->Sumw2();
    #endif

    // central exclusion window
    double excl_lo = left_hi;
    double excl_hi = right_lo;
    int bin_lo = h_coin_bgsub->FindBin(excl_lo);
    int bin_hi = h_coin_bgsub->FindBin(excl_hi);
    double data_integral_excl = h_coin_bgsub->Integral(bin_lo, bin_hi);

    // build cov_poly
    TMatrixDSym cov_poly(poly_order+1);
    for (int i=0;i<=poly_order;++i) for (int j=0;j<=poly_order;++j) {
        if (i < cov_all.GetNrows() && j < cov_all.GetNcols()) cov_poly(i,j) = cov_all(i,j);
        else cov_poly(i,j) = 0.0;
    }

    // subtract with propagated error
    int nbins = h_coin_bgsub->GetNbinsX();
    for (int b=1; b<=nbins; ++b) {
        double x = h_coin_bgsub->GetXaxis()->GetBinCenter(b);
        double data = h_coin_bgsub->GetBinContent(b);
        double data_err = h_coin_bgsub->GetBinError(b);
        if (data_err <= 0.0) data_err = (data>0.0? std::sqrt(data) : 1.0);

        double bgval = f_bg->Eval(x);
        double bg_err = PolyErrorAtX(cov_poly, x);

        double newc = data - bgval;
        double newerr = std::sqrt(data_err*data_err + bg_err*bg_err);

        h_final->SetBinContent(b, newc);
        h_final->SetBinError(b, newerr);
    }

    // Fit Gaussian to h_final in central window (prefit)
    int imax = h_final->GetMaximumBin();
    double mu_guess = h_final->GetBinCenter(imax);
    double A_guess = h_final->GetBinContent(imax);
    double rms_guess = h_final->GetRMS();
    if (rms_guess <= 0.0) rms_guess = (excl_hi - excl_lo) / 6.0;
    double sigma_guess = std::min(rms_guess, 0.02);

    double signal_lo = std::max(h_final->GetXaxis()->GetXmin(), mu_guess - 3.0*sigma_guess);
    double signal_hi = std::min(h_final->GetXaxis()->GetXmax(), mu_guess + 3.0*sigma_guess);

    TF1 f_sig("f_sig", "gaus", signal_lo, signal_hi);
    f_sig.SetParameters(A_guess, mu_guess, sigma_guess);
    TFitResultPtr sigres = h_final->Fit(&f_sig, "RQ");

    double sig_A = f_sig.GetParameter(0);
    double sig_mu = f_sig.GetParameter(1);
    double sig_sigma = std::abs(f_sig.GetParameter(2));
    double sig_A_err = f_sig.GetParError(0);
    double sig_mu_err = f_sig.GetParError(1);
    double sig_sigma_err = f_sig.GetParError(2);
    double gaus_area_cont = f_sig.Integral(signal_lo, signal_hi);

    // compute yields (discrete) in the signal window (bins inclusive)
    int sig_bin_lo = h_final->FindBin(signal_lo);
    int sig_bin_hi = h_final->FindBin(signal_hi);
    double data_counts_win = h_coin_bgsub->Integral(sig_bin_lo, sig_bin_hi);
    double bg_counts_win = 0.0;
    for (int b=sig_bin_lo; b<=sig_bin_hi; ++b) bg_counts_win += f_bg->Eval(h_coin_bgsub->GetXaxis()->GetBinCenter(b));
    double final_counts_win = h_final->Integral(sig_bin_lo, sig_bin_hi);

    // Purity (approx) - still available if needed but we won't display it
    double purity = (final_counts_win > 0.0) ? (final_counts_win / (final_counts_win + bg_counts_win)) : 0.0;

    // Start plotting: large canvas with 4 pads:
    if (draw) {
        // prepare output subdir consistent with workflow
        if (outDir && std::string(outDir).size()>0) gSystem->mkdir(outDir, true);
        TString runDir = (run>=0) ? TString::Format("%s/run_%d", outDir, run) : TString::Format("%s/run_all", outDir);
        gSystem->mkdir(runDir, true);

        TString cname = (run >= 0) ? TString::Format("c_combbg_run%d", run) : TString("c_combbg");
        TCanvas *c = new TCanvas(cname, "Combinatorial BG fit & subtraction", 1400, 1000);
        c->Divide(2,2);

        // -------------------------
        // (1) Top-left: full-range data with fitted background and sidebands
        // -------------------------
        c->cd(1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        // Colors: data = black, fitted bg = orange, sideband boxes = soft cyan, sideband points = azure
        h_coin_bgsub->SetLineColor(kBlack);
        h_coin_bgsub->SetLineWidth(1);
        h_coin_bgsub->SetMarkerStyle(20);
        h_coin_bgsub->SetMarkerSize(0.9);
        h_coin_bgsub->Draw("E"); // draw data with errors

        // fitted background (f_bg is TF1*)
        f_bg->SetLineColor(kOrange+1);
        f_bg->SetLineWidth(2);
        f_bg->SetLineStyle(1);
        f_bg->Draw("SAME");

        // translucent sideband boxes
        TBox bL(left_lo, gPad->GetUymin(), left_hi, gPad->GetUymax());
        bL.SetFillStyle(1001); bL.SetFillColorAlpha(kCyan-9, 0.12); bL.SetLineColor(kCyan-9); bL.Draw("SAME");
        TBox bR(right_lo, gPad->GetUymin(), right_hi, gPad->GetUymax());
        bR.SetFillStyle(1001); bR.SetFillColorAlpha(kCyan-9, 0.12); bR.SetLineColor(kCyan-9); bR.Draw("SAME");

        // central excluded region lines
        TLine linL(excl_lo, gPad->GetUymin(), excl_lo, gPad->GetUymax());
        linL.SetLineStyle(2); linL.SetLineColor(kGray+2); linL.SetLineWidth(1); linL.Draw("Same");
        TLine linR(excl_hi, gPad->GetUymin(), excl_hi, gPad->GetUymax());
        linR.SetLineStyle(2); linR.SetLineColor(kGray+2); linR.SetLineWidth(1); linR.Draw("Same");

        // sideband points on top (TGraphErrors)
        g_side->SetMarkerStyle(21); g_side->SetMarkerSize(0.9); g_side->SetMarkerColor(kAzure+7);
        g_side->Draw("P SAME");

        // legend & compact parameter box (only concise items)
        TLegend leg1(0.58, 0.65, 0.92, 0.88);
        leg1.SetBorderSize(0); leg1.SetFillColor(0);
        leg1.AddEntry(h_coin_bgsub, "Data (coin bgsub)", "lep");
        leg1.AddEntry(g_side, "Sideband points", "p");
        leg1.AddEntry(f_bg, Form("Fitted background (order %d)", poly_order), "l");
        leg1.Draw();

        // -------------------------
        // (2) Top-right: before vs after overlay (line histograms only)
        // -------------------------
        c->cd(2);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        // Draw h_coin_bgsub and h_final as lines (HIST). Keep same axis scale.
        TH1D *h_before = (TH1D*)h_coin_bgsub->Clone(TString::Format("%s_before_draw", h_coin_bgsub->GetName()));
        h_before->SetDirectory(nullptr);
        h_before->SetMarkerStyle(0);
        h_before->SetLineColor(kGray+2);
        h_before->SetLineWidth(2);
        h_before->Draw("HIST");

        h_final->SetLineColor(kBlue-7); // calm deep blue
        h_final->SetLineWidth(2);
        h_final->SetMarkerStyle(0);
        h_final->Draw("HIST SAME");

        // Add legend: clear, only show before/after
        TLegend leg2(0.55, 0.72, 0.92, 0.88);
        leg2.SetBorderSize(0); leg2.SetFillColor(0);
        leg2.AddEntry(h_before, "Before (data)", "l");
        leg2.AddEntry(h_final, "After (data - fitted bg)", "l");
        leg2.Draw();

        // -------------------------
        // (3) Bottom-left: zoom around peak, x axis range 0.08-0.17, μ/σ in MeV
        // -------------------------
        c->cd(3);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        // set fixed x-range in GeV
        const double zoom_xlo = 0.08;
        const double zoom_xhi = 0.17;
        h_final->GetXaxis()->SetRangeUser(zoom_xlo, zoom_xhi);

        // improve axis label / tick sizes for readability
        h_final->GetXaxis()->SetTitleSize(0.045);
        h_final->GetXaxis()->SetLabelSize(0.040);
        h_final->GetYaxis()->SetTitleSize(0.045);
        h_final->GetYaxis()->SetLabelSize(0.040);

        // draw final histogram as line (no markers) for clarity
        h_final->SetLineColor(kBlue-7);
        h_final->SetLineWidth(2);
        h_final->SetMarkerStyle(0);
        h_final->Draw("HIST E"); // "HIST E" keeps errors if present but draws as line

        // draw fitted gaussian
        f_sig.SetLineColor(kRed);
        f_sig.SetLineWidth(2);
        f_sig.SetLineStyle(1);
        f_sig.Draw("Same");

        // Legend: keep it inside the pad, increase text size and remove amplitude/purity
        TLegend leg3(0.60, 0.66, 0.88, 0.86);
        leg3.SetBorderSize(0);
        leg3.SetFillColor(0);
        leg3.SetTextSize(0.040);
        leg3.SetTextFont(42);
        leg3.SetMargin(0.12);

        // Format mu and sigma in MeV and show only those in legend
        const double mu_MeV = sig_mu * 1000.0;
        const double muErr_MeV = sig_mu_err * 1000.0;
        const double sigma_MeV = sig_sigma * 1000.0;
        const double sigmaErr_MeV = sig_sigma_err * 1000.0;
        TString gaussInfo = TString::Format("μ = %4.1f ± %4.1f MeV; σ = %4.1f ± %4.1f MeV",
                                            mu_MeV, muErr_MeV, sigma_MeV, sigmaErr_MeV);

        leg3.AddEntry(h_final, "Final (bg-subtracted)", "l");
        leg3.AddEntry(&f_sig, gaussInfo.Data(), "l");
        leg3.Draw("same");

        // Counts text: print the discrete counts and the exact signal window range used (signal_lo/signal_hi are in GeV)
        TLatex tx3;
        tx3.SetNDC();
        tx3.SetTextSize(0.038);
        tx3.SetTextFont(42);
        tx3.DrawLatex(0.12, 0.92, Form("Run %d", run));
        tx3.DrawLatex(0.12, 0.88, Form("Counts in [%.3f, %.3f] GeV = %.1f", signal_lo, signal_hi, final_counts_win));

        // -------------------------
        // (4) Bottom-right: simplified residual/pull plot (easier to follow)
        // -------------------------
        c->cd(4);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        gPad->SetGridy(true);

        TH1D *h_pull = (TH1D*)h_coin_bgsub->Clone(TString::Format("%s_pull_run%d", h_coin_bgsub->GetName(), run));
        h_pull->Reset();
        h_pull->SetDirectory(nullptr);

        // compute pulls using combined uncertainty (data err + bg err)
        for (int b=1; b<=nbins; ++b) {
            double x = h_coin_bgsub->GetXaxis()->GetBinCenter(b);
            double yi = h_coin_bgsub->GetBinContent(b);
            double ei = h_coin_bgsub->GetBinError(b); if (ei <= 0.0) ei = (yi>0 ? std::sqrt(yi) : 1.0);
            double mval = f_bg->Eval(x);
            double bg_err = PolyErrorAtX(cov_poly, x);
            double tot_err = std::sqrt(ei*ei + bg_err*bg_err);
            double pull = (tot_err > 0.0) ? (yi - mval) / tot_err : 0.0;
            h_pull->SetBinContent(b, pull);
            h_pull->SetBinError(b, 0.0);
        }

        // styling pull plot
        h_pull->SetTitle(";M_{#gamma#gamma} [GeV];(data - bg)/#sigma_{tot}");
        h_pull->GetYaxis()->SetRangeUser(-5.5, 5.5);
        h_pull->SetMarkerStyle(20);
        h_pull->SetMarkerSize(0.8);
        h_pull->SetLineColor(kBlack);
        h_pull->Draw("E");

        // zero line and RMS annotation
        TLine l0(x_min, 0.0, x_max, 0.0);
        l0.SetLineColor(kRed); l0.SetLineWidth(1); l0.Draw("SAME");

        double pull_rms = h_pull->GetRMS();
        TLatex txPR; txPR.SetNDC(); txPR.SetTextSize(0.028);
        txPR.DrawLatex(0.12, 0.92, Form("Pull RMS = %.3f", pull_rms));

        // -------------------------
        // Save canvas and write diagnostics (robust)
        // -------------------------
        c->Update();
        TString png = TString::Format("%s/combbg_run%d_order%d_enhanced.png", runDir.Data(), run, poly_order);
        c->SaveAs(png);
        std::cout << Form("[nps::FitCombinatorialBGAndSubtract] wrote PNG: %s\n", png.Data());

        // write ROOT diagnostics into run directory (clone objects, do not SetDirectory on TF1/TGraph)
        TString rootout = TString::Format("%s/combbg_run%d_order%d_results.root", runDir.Data(), run, poly_order);
        TFile fout(rootout, "RECREATE");
        if (!fout.IsOpen() || fout.IsZombie()) {
            std::cerr << "[nps::FitCombinatorialBGAndSubtract] ERROR: cannot open output root file " << rootout << " for writing. Check permissions.\n";
        } else {
            fout.mkdir(TString::Format("run_%d", run));
            fout.cd(TString::Format("run_%d", run));

            // write histograms (clones)
            if (h_coin_bgsub) { TH1D *h_in = (TH1D*)h_coin_bgsub->Clone("h_coin_bgsub_input"); h_in->Write(); delete h_in; }
            if (h_final)      { TH1D *h_fin = (TH1D*)h_final->Clone("h_bgsub_final"); h_fin->Write(); delete h_fin; }

            // write sideband points and fits (clone & write)
            if (g_side) { TGraphErrors *g_clone = (TGraphErrors*)g_side->Clone("g_sideband_points"); g_clone->Write(); delete g_clone; }
            if (f_bg)   { TF1 *fbg_clone = (TF1*)f_bg->Clone("fitted_bg"); fbg_clone->Write(); delete fbg_clone; }
            // f_sig is TF1 object in this code; copy-construct and write
            TF1 f_sig_clone = f_sig;
            f_sig_clone.SetName("fitted_gaus_on_bgsub");
            f_sig_clone.Write();

            // write pull histogram
            if (h_pull) { TH1D *h_pclone = (TH1D*)h_pull->Clone("h_pull"); h_pclone->Write(); delete h_pclone; }

            fout.Close();
            std::cout << Form("[nps::FitCombinatorialBGAndSubtract] wrote diagnostics to %s\n", rootout.Data());
        }

        // cleanup clones / temporaries created in pads
        if (h_before) { delete h_before; h_before = nullptr; }
        if (h_pull)   { delete h_pull;   h_pull = nullptr; }
        // delete original objects created earlier
    } // end draw

    // clean up objects that were allocated at function start
    if (g_side) { delete g_side; g_side = nullptr; }
    if (f_bg)   { delete f_bg;   f_bg = nullptr; }

    // Return results in struct
    result.h_final = h_final;
    result.chi2_ndf = (ndf>0) ? (chi2 / double(ndf)) : -1.0;
    result.mu_MeV = sig_mu * 1000.0;
    result.sigma_MeV = sig_sigma * 1000.0;
    result.signal_counts = final_counts_win;

    return result;
}

} // namespace nps

#endif // NPS_COMB_BG_H
