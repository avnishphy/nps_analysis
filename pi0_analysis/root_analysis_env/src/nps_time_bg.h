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

namespace nps_bg{
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
// Helper: integrate a TH2D between two axis ranges (x_lo,x_hi) x (y_lo,y_hi)
// Returns raw counts and area (x range * y range).
// NOTE: caller must include <TH2D.h> and pass a valid TH2D*.
// --------------------
inline std::pair<double,double> integral_and_area_TH2(TH2D *h2,
                                                      double x_lo, double x_hi,
                                                      double y_lo, double y_hi)
{
    if (!h2) return {0.0, 0.0};
    int bx1 = h2->GetXaxis()->FindBin(x_lo);
    int bx2 = h2->GetXaxis()->FindBin(x_hi);
    int by1 = h2->GetYaxis()->FindBin(y_lo);
    int by2 = h2->GetYaxis()->FindBin(y_hi);
    // Integral is inclusive of bin indices
    double raw = h2->Integral(bx1, bx2, by1, by2);
    double area = (x_hi - x_lo) * (y_hi - y_lo); // ns^2
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

    // 8) Print a verbose summary to stdout (also available via R.summary())
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

} // namespace nps_bg

#endif // NPS_HELPER_H