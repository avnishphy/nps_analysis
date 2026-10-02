#include <cassert>
#include <cmath>
#include <iostream>

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TString.h>

namespace nps { inline double sqr(double value) { return value * value; } }
#include "../src/analysis/nps_time_bg.h"

void test_legacy_timing_window() {
    TH2D histogram("legacy_window", "", 220, 139.0, 161.0,
                   220, 139.0, 161.0);

    // These shifted-sideband entries were under/overflow in the old
    // 140--160 ns histogram.
    histogram.Fill(139.5, 139.5);
    histogram.Fill(160.5, 160.5);
    histogram.Fill(139.5, 150.0);
    histogram.Fill(150.0, 160.5);
    histogram.Fill(139.5, 160.5);
    histogram.Fill(160.5, 139.5);

    assert(histogram.GetNbinsX() == 220);
    assert(std::abs(histogram.GetXaxis()->GetXmin() - 139.0) < 1.0e-12);
    assert(std::abs(histogram.GetXaxis()->GetXmax() - 161.0) < 1.0e-12);
    assert(std::abs(histogram.GetXaxis()->GetBinWidth(1) - 0.1) < 1.0e-12);

    const std::vector<std::pair<double, double>> side_windows = {
        {139.0, 141.0}, {141.0, 143.0}, {143.0, 145.0},
        {155.0, 157.0}, {157.0, 159.0}, {159.0, 161.0},
    };
    const auto result = nps::estimate_coincidence_background_default(
        &histogram, {149.0, 151.0}, side_windows, side_windows,
        {155.0, 161.0}, {139.0, 145.0},
        {139.0, 145.0}, {155.0, 161.0});

    assert(result.n_diag_raw == 2.0);
    assert(result.n_hor_raw == 1.0);
    assert(result.n_ver_raw == 1.0);
    assert(result.n_full1_raw == 1.0);
    assert(result.n_full2_raw == 1.0);
    std::cout << "legacy 139--161 ns timing window test: PASS\n";
}
