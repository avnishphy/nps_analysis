// ============================================================================
// File: nps_helper.h  (coincidence accidental estimator)
// Author: Avnish Singh (helpers refined + BG estimator by ChatGPT)
// Purpose: Helper and utility functions for NPS π0 and DVCS analysis
// ============================================================================

#ifndef NPS_TIME_BG
#define NPS_TIME_BG

#include <array>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <utility> // pair
#include <string>

// Forward-declare ROOT types here; header using this file that calls the functions
// must include <TH2D.h>, <TH1D.h> etc.
class TH2D;
class TH1D;

namespace nps{
// --------------------
// A small struct to hold the accidental estimate result with errors & components
// --------------------
struct CoincidenceBGResult {
    // raw counts (from histograms) for core/coincidence and sideband boxes
    double n_coin_raw = 0.0;
    double n_diag_raw = 0.0;
    double n_hor_raw = 0.0;
    double n_ver_raw = 0.0;
    double n_full1_raw = 0.0;
    double n_full2_raw = 0.0;

    // areas (ns^2) for each box
    double area_coin = 0.0;
    double area_diag = 0.0;
    double area_hor = 0.0;
    double area_ver = 0.0;
    double area_full1 = 0.0;
    double area_full2 = 0.0;

    // normalized contributions (scaled to coincidence box area)
    double diag_norm = 0.0;
    double hor_norm = 0.0;
    double ver_norm = 0.0;
    double full1_norm = 0.0;
    double full2_norm = 0.0;

    // final estimated accidentals in the coincidence box and propagated error
    double n_accidentals = 0.0;
    double n_accidentals_err = 0.0;

    // convenience: total background per-component as printed strings
    std::string summary() const {
        char buf[512];
        sprintf(buf,
                "coin_raw=%.1f diag_raw=%.1f hor_raw=%.1f ver_raw=%.1f full1_raw=%.1f full2_raw=%.1f\n"
                "areas (ns^2): coin=%.3f diag=%.3f hor=%.3f ver=%.3f full1=%.3f full2=%.3f\n"
                "normalized (to coin area): diag=%.3f hor=%.3f ver=%.3f full1=%.3f full2=%.3f\n"
                "estimated accidentals = %.3f +- %.3f",
                n_coin_raw, n_diag_raw, n_hor_raw, n_ver_raw, n_full1_raw, n_full2_raw,
                area_coin, area_diag, area_hor, area_ver, area_full1, area_full2,
                diag_norm, hor_norm, ver_norm, full1_norm, full2_norm,
                n_accidentals, n_accidentals_err);
        return std::string(buf);
    }
};

// --------------------
// Helper: count entries in a TH2D between two axis ranges (x_lo,x_hi) x (y_lo,y_hi)
// Returns raw counts (number of entries) and area (x range * y range).
// NOTE: caller must include <TH2D.h> and pass a valid TH2D*.
// --------------------
inline std::pair<double,double> integral_and_area_TH2(TH2D *h2,
                                                      double x_lo, double x_hi,
                                                      double y_lo, double y_hi)
{
    if (!h2) return {0.0, 0.0};

    // Adjust high edges to match < x_hi logical selection
    double eps = 1e-9;
    int bx1 = h2->GetXaxis()->FindBin(x_lo);
    int bx2 = h2->GetXaxis()->FindBin(x_hi - eps);
    int by1 = h2->GetYaxis()->FindBin(y_lo);
    int by2 = h2->GetYaxis()->FindBin(y_hi - eps);

    bx1 = std::max(1, bx1);
    bx2 = std::min(h2->GetNbinsX(), bx2);
    by1 = std::max(1, by1);
    by2 = std::min(h2->GetNbinsY(), by2);

    double raw = 0.0;
    for (int ix = bx1; ix <= bx2; ++ix)
        for (int iy = by1; iy <= by2; ++iy)
            raw += h2->GetBinContent(ix, iy);

    double area = (x_hi - x_lo) * (y_hi - y_lo);
    return {raw, area};
}



// --------------------
// Default windows based on your comment. These are used by the convenience wrapper
// and can be tuned by the user. All ranges are [lo,hi] in ns.
// --------------------
inline std::vector<std::pair<double,double>> default_diag_windows() {
    return {
        {141.0, 143.0},
        {143.0, 145.0},
        {145.0, 147.0},
        {153.0, 155.0},
        {155.0, 157.0},
        {157.0, 159.0}
    };
}
inline std::vector<std::pair<double,double>> default_side_windows() {
    // same as diag windows above; used for horizontal / vertical combos
    return default_diag_windows();
}
inline std::pair<double,double> default_coin_window() {
    // default 2-ns coincidence window around 150 ns
    return {149.0, 151.0};
}
inline std::pair<double,double> default_full_acc1_t1() { return {153.0, 159.0}; }
inline std::pair<double,double> default_full_acc1_t2() { return {141.0, 147.0}; }
inline std::pair<double,double> default_full_acc2_t2() { return {153.0, 159.0}; }
inline std::pair<double,double> default_full_acc2_t1() { return {141.0, 147.0}; }

// --------------------
// The main estimation function (convenience wrapper using the defaults described in your comment).
// Inputs:
//   h2: TH2D* with t2 on X axis and t1 on Y axis (user requested convention t1 = y, t2 = x).
//   coin_win: pair (t_lo,t_hi) for both t1 and t2 coincidence window; default {149,151}.
//   diag_windows: vector of small diagonal sideband windows (same window for t1 & t2 each).
//   side_windows: same as diag_windows by default (used to build horizontal/vertical boxes).
//   full_acc windows: two large boxes used in subtraction to correct double counting.
// Returns CoincidenceBGResult populated with raws, normalized contributions and final estimate.
// --------------------
inline CoincidenceBGResult estimate_coincidence_background_default(TH2D *h2,
                                                                   std::pair<double,double> coin_win = default_coin_window(),
                                                                   std::vector<std::pair<double,double>> diag_windows = default_diag_windows(),
                                                                   std::vector<std::pair<double,double>> side_windows = default_side_windows(),
                                                                   std::pair<double,double> full1_t1 = default_full_acc1_t1(),
                                                                   std::pair<double,double> full1_t2 = default_full_acc1_t2(),
                                                                   std::pair<double,double> full2_t2 = default_full_acc2_t2(),
                                                                   std::pair<double,double> full2_t1 = default_full_acc2_t1())
{
    CoincidenceBGResult R;

    if (!h2) {
        std::cerr << "[nps::estimate_coincidence_background_default] ERROR: null TH2D* passed\n";
        return R;
    }

    const double coin_x_lo = coin_win.first;
    const double coin_x_hi = coin_win.second;
    const double coin_y_lo = coin_win.first;
    const double coin_y_hi = coin_win.second;

    // 1) coin raw
    auto coin_pair = integral_and_area_TH2(h2, coin_x_lo, coin_x_hi, coin_y_lo, coin_y_hi);
    R.n_coin_raw = coin_pair.first;
    R.area_coin = coin_pair.second;

    // 2) diagonal contributions: sum boxes where both t1 and t2 are in same small window
    double diag_raw_sum = 0.0;
    double diag_area_sum = 0.0;
    for (auto &w : diag_windows) {
        auto pr = integral_and_area_TH2(h2, w.first, w.second, w.first, w.second);
        diag_raw_sum += pr.first;
        diag_area_sum += pr.second;
    }
    R.n_diag_raw = diag_raw_sum;
    R.area_diag = diag_area_sum;
    // normalized to coin area:
    if (R.area_diag > 0) R.diag_norm = R.n_diag_raw * (R.area_coin / R.area_diag);
    else R.diag_norm = 0.0;

    // 3) horizontal contributions: t1 in coin window, t2 in side windows
    double hor_raw_sum = 0.0;
    double hor_area_sum = 0.0;
    for (auto &w : side_windows) {
        auto pr = integral_and_area_TH2(h2, w.first, w.second, coin_y_lo, coin_y_hi); // x: side, y: coin
        hor_raw_sum += pr.first;
        hor_area_sum += pr.second;
    }
    R.n_hor_raw = hor_raw_sum;
    R.area_hor = hor_area_sum;
    if (R.area_hor > 0) R.hor_norm = R.n_hor_raw * (R.area_coin / R.area_hor);
    else R.hor_norm = 0.0;

    // 4) vertical contributions: t2 in coin window, t1 in side windows
    double ver_raw_sum = 0.0;
    double ver_area_sum = 0.0;
    for (auto &w : side_windows) {
        auto pr = integral_and_area_TH2(h2, coin_x_lo, coin_x_hi, w.first, w.second); // x: coin, y: side
        ver_raw_sum += pr.first;
        ver_area_sum += pr.second;
    }
    R.n_ver_raw = ver_raw_sum;
    R.area_ver = ver_area_sum;
    if (R.area_ver > 0) R.ver_norm = R.n_ver_raw * (R.area_coin / R.area_ver);
    else R.ver_norm = 0.0;

    // 5) "complete accidental" boxes used for subtraction (two rectangles)
    auto full1_pr = integral_and_area_TH2(h2, full1_t2.first, full1_t2.second, full1_t1.first, full1_t1.second);
    auto full2_pr = integral_and_area_TH2(h2, full2_t1.first, full2_t1.second, full2_t2.first, full2_t2.second);
    R.n_full1_raw = full1_pr.first;
    R.area_full1 = full1_pr.second;
    R.n_full2_raw = full2_pr.first;
    R.area_full2 = full2_pr.second;
    if (R.area_full1 > 0) R.full1_norm = R.n_full1_raw * (R.area_coin / R.area_full1);
    else R.full1_norm = 0.0;
    if (R.area_full2 > 0) R.full2_norm = R.n_full2_raw * (R.area_coin / R.area_full2);
    else R.full2_norm = 0.0;

    // 6) Combine using the formula you specified:
    // n_accidentals = diag_contri + 0.5*(vertical_contri+horizontal_contri) - 0.5*(complete_acc_contri1+complete_acc_contri2)
    R.n_accidentals = R.diag_norm + 0.5 * (R.ver_norm + R.hor_norm) - 0.5 * (R.full1_norm + R.full2_norm);

    // 7) Propagate errors (Poisson on raw counts, scale by same scaling factors used for normalization).
    // Errors on normalized quantities:
    //   err(diag_norm) = sqrt(diag_raw_sum) * (area_coin/area_diag_sum)
    double err_diag = (R.area_diag > 0) ? std::sqrt(std::max(0.0, R.n_diag_raw)) * (R.area_coin / R.area_diag) : 0.0;
    double err_hor  = (R.area_hor  > 0) ? std::sqrt(std::max(0.0, R.n_hor_raw))  * (R.area_coin / R.area_hor)  : 0.0;
    double err_ver  = (R.area_ver  > 0) ? std::sqrt(std::max(0.0, R.n_ver_raw))  * (R.area_coin / R.area_ver)  : 0.0;
    double err_full1= (R.area_full1> 0) ? std::sqrt(std::max(0.0, R.n_full1_raw)) * (R.area_coin / R.area_full1) : 0.0;
    double err_full2= (R.area_full2> 0) ? std::sqrt(std::max(0.0, R.n_full2_raw)) * (R.area_coin / R.area_full2) : 0.0;

    // n_accidentals = diag + 0.5*(ver+hor) - 0.5*(full1+full2)
    // errors add in quadrature with scaling factors:
    R.n_accidentals_err = std::sqrt( nps::sqr(err_diag) + nps::sqr(0.5*err_ver) + nps::sqr(0.5*err_hor) + nps::sqr(0.5*err_full1) + nps::sqr(0.5*err_full2) );

    // // 8) Print a verbose summary to stdout (also available via R.summary())
    std::cout << "[nps::estimate_coincidence_background_default] Summary:\n";
    std::cout << " Coincidence box: [" << coin_x_lo << "," << coin_x_hi << "] x [" << coin_y_lo << "," << coin_y_hi << "]\n";
    std::cout << "  raw coin counts = " << R.n_coin_raw << "   area (ns^2) = " << R.area_coin << "\n";
    std::cout << " Diagonal sideband raw sum = " << R.n_diag_raw << "   total area = " << R.area_diag
              << "  normalized -> " << R.diag_norm << "\n";
    std::cout << " Horizontal sideband raw sum = " << R.n_hor_raw << "   total area = " << R.area_hor
              << "  normalized -> " << R.hor_norm << "\n";
    std::cout << " Vertical sideband raw sum = " << R.n_ver_raw << "   total area = " << R.area_ver
              << "  normalized -> " << R.ver_norm << "\n";
    std::cout << " Full accidental box1 raw = " << R.n_full1_raw << "  area = " << R.area_full1 << "  norm -> " << R.full1_norm << "\n";
    std::cout << " Full accidental box2 raw = " << R.n_full2_raw << "  area = " << R.area_full2 << "  norm -> " << R.full2_norm << "\n";
    std::cout << " Final estimated accidental counts in coin box = " << R.n_accidentals
              << " +/- " << R.n_accidentals_err << "\n";

    // return the struct
    return R;
}

inline TH1D* make_and_subtract_accidentals_data_driven(
                TH1D *h_m_pi0_coin,
                const nps::CoincidenceBGResult &bg,
                TH1D *h_mgg_full1, TH1D *h_mgg_full2,
                const std::vector<TH1D*> &h_mgg_hor,
                const std::vector<TH1D*> &h_mgg_ver,
                const std::vector<TH1D*> &h_mgg_diag,
                TH1D *h_m_pi0_acc, // fallback optional
                TH2D *h_t1_t2,
                const std::vector<std::pair<double,double>> &diag_windows,
                const std::vector<std::pair<double,double>> &side_windows,
                const TString &outPlotDir = "", int run = -1)
{
    if (!h_m_pi0_coin) { std::cerr << "[make_and_subtract_accidentals_data_driven] null h_m_pi0_coin\n"; return nullptr; }
    h_m_pi0_coin->Sumw2();

    // Helper: safe clone with Sumw2 and detached directory
    auto safe_clone = [&](TH1D *h, const char *newname)->TH1D* {
        if (!h) return (TH1D*)nullptr;
        TH1D *c = (TH1D*)h->Clone(newname);
        c->SetDirectory(nullptr);
        c->Sumw2();
        c->SetTitle(TString::Format("%s (scaled template)", h->GetTitle()).Data());
        return c;
    };

    // Helper to compute area of a box using integral_and_area_TH2
    auto area_of_box = [&](double xlo, double xhi, double ylo, double yhi)->double {
        if (!h_t1_t2) return 0.0;
        auto pr = nps::integral_and_area_TH2(h_t1_t2, xlo, xhi, ylo, yhi);
        return pr.second;
    };

    const double coin_area = bg.area_coin;
    const double norm_unc_frac = (bg.n_accidentals > 0.0) ? (bg.n_accidentals_err / bg.n_accidentals) : 0.0;
    std::vector<TH1D*> tmp_clones; tmp_clones.reserve(16);

    
    // Create empty template that will hold the algebraic combination:
    // h_template = diag_sum * 1.0 + 0.5*(hor_sum + ver_sum) - 0.5*(full1 + full2)
    TH1D *h_template = (TH1D*)h_m_pi0_coin->Clone(TString::Format("h_template_run%d", run).Data());
    h_template->Reset();
    h_template->Sumw2();

    // ---------- 1) DIAGONAL: sum per-window histograms, sum areas, then scale by coin_area/total_area ----------
    double area_d_sum = 0.0;
    TH1D *diag_comb = nullptr;
    for (size_t i = 0; i < h_mgg_diag.size() && i < diag_windows.size(); ++i) {
        TH1D *hwin = h_mgg_diag[i];
        double xlo = diag_windows[i].first, xhi = diag_windows[i].second;
        double area_d = area_of_box(xlo, xhi, xlo, xhi);
        area_d_sum += area_d;

        if (!hwin) {
            printf("[DEBUG] diag[%zu]: input hist is null, skipping\n", i);
            continue;
        }

        if (!diag_comb) {
            // clone first non-null window as the accumulator (detached, Sumw2 on)
            diag_comb = safe_clone(hwin, TString::Format("%s_diag_comb_run%d", hwin->GetName(), run).Data());
            tmp_clones.push_back(diag_comb);
        } else {
            // basic binning compatibility check
            if (diag_comb->GetNbinsX() == hwin->GetNbinsX()) {
                diag_comb->Add(hwin);
            } else {
                printf("[WARNING] diag[%zu]: histogram binning mismatch (skipping add)\n", i);
            }
        }

        printf("[DEBUG] diag[%zu]: window=[%.3f,%.3f] area=%.6g  hist_int=%g\n",
               i, xlo, xhi, area_d, hwin ? hwin->Integral() : 0.0);
    }

    if (diag_comb && area_d_sum > 0.0) {
        double scale = coin_area / area_d_sum;
        printf("[DEBUG] diag: total_area=%.6g  coin_area=%.6g  applying scale=%.6g\n", area_d_sum, coin_area, scale);
        diag_comb->Scale(scale);
        printf("[DEBUG] diag: combined HIST integral AFTER scaling = %g (before scale was %g)\n",
               diag_comb->Integral(), diag_comb->Integral() / fabs(scale > 0 ? scale : 1.0));
        h_template->Add(diag_comb, 1.0); // weight +1
    } else {
        if (!diag_comb) printf("[DEBUG] diag: no diag histograms were available to combine\n");
        else printf("[DEBUG] diag: total area sum <= 0 (%.6g) -> skipping diag contribution\n", area_d_sum);
    }

    // ---------- 2) HORIZONTAL: sum per-window histograms (t2 side, t1 coin), total-area scaling, weight=0.5 ----------
    double area_h_sum = 0.0;
    TH1D *hor_comb = nullptr;
    const double coin_y_lo = nps::default_coin_window().first;
    const double coin_y_hi = nps::default_coin_window().second;
    for (size_t i = 0; i < h_mgg_hor.size() && i < side_windows.size(); ++i) {
        TH1D *hwin = h_mgg_hor[i];
        double xlo = side_windows[i].first, xhi = side_windows[i].second;
        auto pr = nps::integral_and_area_TH2(h_t1_t2, xlo, xhi, coin_y_lo, coin_y_hi);
        double area_h = pr.second;
        area_h_sum += area_h;

        if (!hwin) {
            printf("[DEBUG] hor[%zu]: input hist is null, skipping\n", i);
            continue;
        }
        if (!hor_comb) {
            hor_comb = safe_clone(hwin, TString::Format("%s_hor_comb_run%d", hwin->GetName(), run).Data());
            tmp_clones.push_back(hor_comb);
        } else {
            if (hor_comb->GetNbinsX() == hwin->GetNbinsX()) hor_comb->Add(hwin);
            else printf("[WARNING] hor[%zu]: histogram binning mismatch (skipping add)\n", i);
        }
        printf("[DEBUG] hor[%zu]: x=[%.3f,%.3f], y=coin[%.3f,%.3f] area=%.6g hist_int=%g\n",
               i, xlo, xhi, coin_y_lo, coin_y_hi, area_h, hwin ? hwin->Integral() : 0.0);
    }

    if (hor_comb && area_h_sum > 0.0) {
        double scale = coin_area / area_h_sum;
        printf("[DEBUG] hor: total_area=%.6g  coin_area=%.6g  applying scale=%.6g (weight 0.5)\n", area_h_sum, coin_area, scale);
        hor_comb->Scale(scale);
        printf("[DEBUG] hor: combined HIST integral AFTER scaling = %g\n", hor_comb->Integral());
        h_template->Add(hor_comb, 0.5);
    } else {
        if (!hor_comb) printf("[DEBUG] hor: no hor histograms provided\n");
        else printf("[DEBUG] hor: total area sum <= 0 (%.6g) -> skipping hor contribution\n", area_h_sum);
    }

    // ---------- 3) VERTICAL: sum per-window histograms (t2 coin, t1 side), total-area scaling, weight=0.5 ----------
    double area_v_sum = 0.0;
    TH1D *ver_comb = nullptr;
    const double coin_x_lo = nps::default_coin_window().first;
    const double coin_x_hi = nps::default_coin_window().second;
    for (size_t i = 0; i < h_mgg_ver.size() && i < side_windows.size(); ++i) {
        TH1D *hwin = h_mgg_ver[i];
        double ylo = side_windows[i].first, yhi = side_windows[i].second;
        auto pr = nps::integral_and_area_TH2(h_t1_t2, coin_x_lo, coin_x_hi, ylo, yhi);
        double area_v = pr.second;
        area_v_sum += area_v;

        if (!hwin) {
            printf("[DEBUG] ver[%zu]: input hist is null, skipping\n", i);
            continue;
        }
        if (!ver_comb) {
            ver_comb = safe_clone(hwin, TString::Format("%s_ver_comb_run%d", hwin->GetName(), run).Data());
            tmp_clones.push_back(ver_comb);
        } else {
            if (ver_comb->GetNbinsX() == hwin->GetNbinsX()) ver_comb->Add(hwin);
            else printf("[WARNING] ver[%zu]: histogram binning mismatch (skipping add)\n", i);
        }
        printf("[DEBUG] ver[%zu]: x=coin[%.3f,%.3f], y=[%.3f,%.3f] area=%.6g hist_int=%g\n",
               i, coin_x_lo, coin_x_hi, ylo, yhi, area_v, hwin ? hwin->Integral() : 0.0);
    }

    if (ver_comb && area_v_sum > 0.0) {
        double scale = coin_area / area_v_sum;
        printf("[DEBUG] ver: total_area=%.6g  coin_area=%.6g applying scale=%.6g (weight 0.5)\n", area_v_sum, coin_area, scale);
        ver_comb->Scale(scale);
        printf("[DEBUG] ver: combined HIST integral AFTER scaling = %g\n", ver_comb->Integral());
        h_template->Add(ver_comb, 0.5);
    } else {
        if (!ver_comb) printf("[DEBUG] ver: no ver histograms provided\n");
        else printf("[DEBUG] ver: total area sum <= 0 (%.6g) -> skipping ver contribution\n", area_v_sum);
    }

    // ---------- 4) Full boxes ----------
    double full1_contribution = 0.0, full2_contribution = 0.0;
    if (h_mgg_full1 && bg.area_full1 > 0.0) {
        TH1D *tmp = safe_clone(h_mgg_full1, TString::Format("%s_full1_tmp_run%d", h_mgg_full1->GetName(), run).Data());
        tmp_clones.push_back(tmp);
        double before = tmp->Integral();
        tmp->Scale( coin_area / bg.area_full1 );
        double after = tmp->Integral();
        printf("[DEBUG] full1: bg.area=%.6g scale=%.6g before=%g after=%g\n", bg.area_full1, coin_area / bg.area_full1, before, after);
        h_template->Add(tmp, -0.5);
        full1_contribution = -0.5 * after;
    } else {
        if (!h_mgg_full1) printf("[DEBUG] full1: none provided\n");
        else printf("[DEBUG] full1: bg.area_full1<=0 (%.6g), skipping\n", bg.area_full1);
    }

    if (h_mgg_full2 && bg.area_full2 > 0.0) {
        TH1D *tmp = safe_clone(h_mgg_full2, TString::Format("%s_full2_tmp_run%d", h_mgg_full2->GetName(), run).Data());
        tmp_clones.push_back(tmp);
        double before = tmp->Integral();
        tmp->Scale( coin_area / bg.area_full2 );
        double after = tmp->Integral();
        printf("[DEBUG] full2: bg.area=%.6g scale=%.6g before=%g after=%g\n", bg.area_full2, coin_area / bg.area_full2, before, after);
        h_template->Add(tmp, -0.5);
        full2_contribution = -0.5 * after;
    } else {
        if (!h_mgg_full2) printf("[DEBUG] full2: none provided\n");
        else printf("[DEBUG] full2: bg.area_full2<=0 (%.6g), skipping\n", bg.area_full2);
    }


    // If we added nothing (no sidebands/full), fallback to h_m_pi0_acc scaled to bg.n_accidentals (if available)
    bool have_template_content = false;
    for (int b=1; b<=h_template->GetNbinsX(); ++b) {
        if (h_template->GetBinContent(b) != 0.0) { have_template_content = true; break; }
    }
    if (!have_template_content) {
        if (h_m_pi0_acc && h_m_pi0_acc->Integral() > 0.0) {
            delete h_template;
            h_template = safe_clone(h_m_pi0_acc, TString::Format("%s_acc_fallback_run%d", h_m_pi0_acc->GetName(), run).Data());
            // normalize to coin area per-bin? We'll just scale to bg.n_accidentals below.
        } else {
            std::cerr << "[make_and_subtract_accidentals_data_driven] ERROR: no templates available and fallback not possible.\n";
            for (auto *p : tmp_clones) delete p;
            return nullptr;
        }
    }

    // Compute template integral & error (before scaling to bg.n_accidentals)
    double templ_int_err = 0.0;
    int ib_lo = 1;
    int ib_hi = h_template->GetNbinsX();
    double templ_integral = h_template->IntegralAndError(ib_lo, ib_hi, templ_int_err);

    // double templ_int = h_template->IntegralAndError(1, h_template->GetNbinsX(), templ_int_err);
    std::cout << "constructed template integral = " << templ_integral << "\n";
    std::cout << "timing-plane estimate (bg.n_accidentals) = " << bg.n_accidentals << " +- " << bg.n_accidentals_err << "\n";

    // If zero integral and fallback exists, try fallback logic already handled above; if still zero, abort.
    if (templ_integral <= 0.0) {
        std::cerr << "[make_and_subtract_accidentals_data_driven] WARNING: template integral <= 0.\n";
        for (auto *p : tmp_clones) delete p;
        delete h_template;
        return nullptr;
    }

    // 2) Subtract template from h_m_pi0_coin (bin-by-bin) with improved error propagation
    TH1D *h_coin_bgsub = (TH1D*)h_m_pi0_coin->Clone(TString::Format("%s_bgsub_run%d", h_m_pi0_coin->GetName(), run).Data());
    h_coin_bgsub->SetDirectory(nullptr);
    h_coin_bgsub->SetTitle(TString::Format("%s (acc-subtracted)", h_m_pi0_coin->GetTitle()).Data());
    h_coin_bgsub->Sumw2();

    for (int b = 1; b <= h_m_pi0_coin->GetNbinsX(); ++b) {
        double orig = h_m_pi0_coin->GetBinContent(b);
        double orig_err = h_m_pi0_coin->GetBinError(b);
        double sub = h_template->GetBinContent(b);
        double sub_err_stat = h_template->GetBinError(b);                    // template statistical/bin error
        double sub_err_norm = fabs(sub) * norm_unc_frac;                    // normalization uncertainty propagated to bin
        double total_sub_err = std::sqrt(sub_err_stat*sub_err_stat + sub_err_norm*sub_err_norm);

        double newc = orig - sub;
        if (newc < 0.0) newc = 0.0;

        double newerr = std::sqrt(orig_err*orig_err + total_sub_err*total_sub_err);
        h_coin_bgsub->SetBinContent(b, newc);
        h_coin_bgsub->SetBinError(b, newerr);
    }

// ------------------------
// Single-plot: before, template (scaled), after; no statbox; legend+info top-right
// ------------------------
{
    // disable default statbox
    gStyle->SetOptStat(0);

    // clones for drawing (keep originals intact)
    TH1D *h_before = h_m_pi0_coin;
    TH1D *h_tpl    = (TH1D*)h_template->Clone(TString::Format("%s_draw_clone_run%d", h_template->GetName(), run).Data());
    h_tpl->SetDirectory(nullptr);
    TH1D *h_after  = h_coin_bgsub;

    // style lines (clear and distinct)
    h_before->SetLineColor(kBlack);
    h_before->SetLineWidth(2);
    h_before->SetLineStyle(1);
    h_before->SetMarkerStyle(0);

    h_tpl->SetLineColor(kRed);
    h_tpl->SetLineWidth(2);
    h_tpl->SetLineStyle(7); // dashed
    h_tpl->SetMarkerStyle(0);

    h_after->SetLineColor(kBlue+2);
    h_after->SetLineWidth(2);
    h_after->SetLineStyle(1);
    h_after->SetMarkerStyle(0);

    // compute y-range to show all three comfortably
    double max_before = (h_before->GetEntries()>0) ? h_before->GetMaximum() : 0.0;
    double max_tpl    = (h_tpl   ->GetEntries()>0) ? h_tpl   ->GetMaximum() : 0.0;
    double max_after  = (h_after ->GetEntries()>0) ? h_after ->GetMaximum() : 0.0;
    double ymax = std::max({max_before, max_tpl, max_after});
    if (ymax <= 0) ymax = 1.0;
    h_before->SetMaximum(ymax * 1.20); // small headroom

    // canvas
    TString cname = (run >= 0) ? TString::Format("c_pi0_acc_sub_run%d_oneplot", run) : TString("c_pi0_acc_sub_oneplot");
    TCanvas *c = new TCanvas(cname, "Pi0 accidental subtraction", 1000, 700);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.12);
    c->cd();

    // draw: base histogram defines axes
    h_before->Draw("HIST");          // before (raw) drawn first so axis fits
    h_tpl->Draw("HIST SAME");        // template (already scaled to bg.n_accidentals)
    h_after->Draw("HIST SAME");      // after subtraction

    // Legend and info box on top-right, arranged to not overlap
   // === Info box (tight, top-right) ===
    TPaveText *info = new TPaveText(0.53, 0.56, 0.87, 0.77, "NDC");
    info->SetFillColor(0);
    info->SetBorderSize(0);
    info->SetTextFont(42);
    info->SetTextSize(0.030);
    info->AddText(Form("Run: %d", run));
    info->AddText(Form("Coin integral: %.2f", h_before->Integral()));
    info->AddText(Form("Timing accidental model: %.3f",
                       templ_integral));
    info->AddText(Form("Final (coin - model): %.3f", h_after->Integral()));
    info->Draw("same");

    // === Legend (directly below info box, fully separated) ===
    TLegend *leg = new TLegend(0.53, 0.32, 0.87, 0.46, "", "NDC");
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);

    leg->AddEntry(h_before, "All coincidence window", "l");
    leg->AddEntry(h_tpl,    "Timing accidental model", "l");
    leg->AddEntry(h_after,  "Coin - model (final)", "l");

    leg->Draw("same");

    gPad->Modified();
    gPad->Update();

    // optional: draw a vertical line or shaded region around pi0 peak if you want focus (commented)
    // TLine *peakL = new TLine(m_peak - width, 0, m_peak - width, ymax*1.2); peakL->SetLineStyle(2); peakL->Draw("same");

    // Save if requested
    if (outPlotDir != "") {
        TString png = TString::Format("%s/pi0_acc_template_sub_run%d_oneplot.png", outPlotDir.Data(), run);
        c->SaveAs(png);
    }

    // cleanup
    delete h_tpl; // safe: it was cloned and detached
    // keep c/h_before/h_after for caller if they want; clones removed
}




    // cleanup temporary clones
    for (auto *p : tmp_clones) { if (p) { delete p; } }
    // keep h_template and h_coin_bgsub in memory for caller to use/write
    return h_coin_bgsub;
}


} // namespace nps

#endif // NPS_HELPER_H