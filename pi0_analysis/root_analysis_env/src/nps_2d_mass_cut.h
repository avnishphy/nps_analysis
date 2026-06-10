#ifndef NPS_2D_MASS_CUT_H
#define NPS_2D_MASS_CUT_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <limits>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TPolyLine.h"
#include "TStyle.h"
#include "TSystem.h"

namespace nps2d {

struct Point {
    Long64_t id = 0;
    double mpi0 = 0.0;
    double mmiss = 0.0;
    double weight = 0.0;
};

struct Config {
    double peak_fraction = -1.0;
    int n_mpi0_bins = 160;
    double mpi0_min = 0.11;
    double mpi0_max = 0.15;
    int n_mmiss_bins = 200;
    double mmiss_min = 0.6;
    double mmiss_max = 1.5;
    double ellipse_core_quantile = 0.995;
    double ellipse_padding = 1.05;
    double auto_peak_min = 0.05;
    double auto_peak_max = 0.60;
    double auto_peak_step = 0.005;
    double auto_min_core_total_fraction = 0.05;
    int auto_min_core_bins = 30;
    double mcd_keep_total_fraction = 0.30;
    int mcd_iterations = 8;
    double mcd_ellipse_quantile = 0.995;
    double mcd_padding = 1.05;
    bool write_debug = true;
    std::string output_dir = ".";
    std::string tag = "mass_cut";
};

struct Params {
    bool valid = false;
    bool auto_peak_fraction = false;
    double peak_fraction = 0.0;
    double auto_jump_ratio = 1.0;
    double peak_weight = 0.0;
    double threshold_weight = 0.0;
    double total_weight = 0.0;
    int peak_ix = 0;
    int peak_iy = 0;
    double peak_mpi0 = 0.0;
    double peak_mmiss = 0.0;
    int core_bins = 0;
    double core_weight = 0.0;
    double core_total_fraction = 0.0;

    double mean_mpi0 = 0.0;
    double mean_mmiss = 0.0;
    double cov_mpi0_mpi0 = 0.0;
    double cov_mpi0_mmiss = 0.0;
    double cov_mmiss_mmiss = 0.0;
    double cov_det = 0.0;
    double ellipse_d2_cut = 0.0;
    double ellipse_weight = 0.0;
    double ellipse_total_fraction = 0.0;

    double mcd_mean_mpi0 = 0.0;
    double mcd_mean_mmiss = 0.0;
    double mcd_cov_mpi0_mpi0 = 0.0;
    double mcd_cov_mpi0_mmiss = 0.0;
    double mcd_cov_mmiss_mmiss = 0.0;
    double mcd_det = 0.0;
    double mcd_d2_cut = 0.0;
    double mcd_subset_weight = 0.0;
    int mcd_subset_bins = 0;
    double mcd_weight = 0.0;
    double mcd_total_fraction = 0.0;
};

struct Result {
    Params params;
    std::vector<int> pass_ellipse;
    std::vector<int> pass_mcd;
};

namespace detail {

struct WeightedValue {
    double value = 0.0;
    double weight = 0.0;
};

struct WeightedPoint {
    double x = 0.0;
    double y = 0.0;
    double weight = 0.0;
    int global_bin = 0;
};

struct CovModel {
    double mean_x = 0.0;
    double mean_y = 0.0;
    double cov_xx = 0.0;
    double cov_xy = 0.0;
    double cov_yy = 0.0;
    double det = 0.0;
};

struct CoreStats {
    double peak_fraction = 0.0;
    double weight = 0.0;
    double total_fraction = 0.0;
    int bins = 0;
};

struct IndexedDistance {
    double d2 = 0.0;
    double weight = 0.0;
    std::size_t index = 0;
};

inline std::string replace_suffix(const std::string& path, const std::string& suffix) {
    const std::size_t dot = path.find_last_of('.');
    if (dot == std::string::npos) return path + suffix;
    return path.substr(0, dot) + suffix;
}

inline double weighted_quantile(std::vector<WeightedValue>& values, double quantile) {
    if (values.empty()) return 0.0;
    if (quantile <= 0.0) return values.front().value;
    if (quantile > 1.0) quantile = 1.0;

    std::sort(values.begin(), values.end(),
              [](const WeightedValue& a, const WeightedValue& b) { return a.value < b.value; });

    double total = 0.0;
    for (const auto& v : values) total += v.weight;
    if (total <= 0.0) return values.back().value;

    const double target = quantile * total;
    double running = 0.0;
    for (const auto& v : values) {
        running += v.weight;
        if (running >= target) return v.value;
    }
    return values.back().value;
}

inline double covariance_d2(const CovModel& model, double x, double y) {
    const double dx = x - model.mean_x;
    const double dy = y - model.mean_y;
    return (model.cov_yy * dx * dx - 2.0 * model.cov_xy * dx * dy +
            model.cov_xx * dy * dy) / model.det;
}

inline bool compute_cov_model(const std::vector<WeightedPoint>& points,
                              const std::vector<unsigned char>& mask,
                              CovModel& model,
                              double& weight_sum,
                              int& n_bins)
{
    weight_sum = 0.0;
    n_bins = 0;
    model = CovModel();

    for (const auto& p : points) {
        if (!mask[p.global_bin]) continue;
        weight_sum += p.weight;
        model.mean_x += p.weight * p.x;
        model.mean_y += p.weight * p.y;
        ++n_bins;
    }
    if (weight_sum <= 0.0 || n_bins < 3) return false;

    model.mean_x /= weight_sum;
    model.mean_y /= weight_sum;

    for (const auto& p : points) {
        if (!mask[p.global_bin]) continue;
        const double dx = p.x - model.mean_x;
        const double dy = p.y - model.mean_y;
        model.cov_xx += p.weight * dx * dx;
        model.cov_xy += p.weight * dx * dy;
        model.cov_yy += p.weight * dy * dy;
    }
    model.cov_xx /= weight_sum;
    model.cov_xy /= weight_sum;
    model.cov_yy /= weight_sum;
    model.det = model.cov_xx * model.cov_yy - model.cov_xy * model.cov_xy;
    return model.det > 0.0;
}

inline TPolyLine* make_ellipse_line(const CovModel& model, double d2_cut,
                                    Color_t color, Style_t style, Width_t width)
{
    const double trace = model.cov_xx + model.cov_yy;
    const double diff = model.cov_xx - model.cov_yy;
    const double root = std::sqrt(0.25 * diff * diff + model.cov_xy * model.cov_xy);
    const double lambda1 = std::max(0.0, 0.5 * trace + root);
    const double lambda2 = std::max(0.0, 0.5 * trace - root);
    const double angle = 0.5 * std::atan2(2.0 * model.cov_xy, diff);
    const double ca = std::cos(angle);
    const double sa = std::sin(angle);
    const double axis1 = std::sqrt(std::max(0.0, d2_cut * lambda1));
    const double axis2 = std::sqrt(std::max(0.0, d2_cut * lambda2));

    const int n = 240;
    TPolyLine* line = new TPolyLine(n + 1);
    for (int i = 0; i <= n; ++i) {
        const double t = 2.0 * TMath::Pi() * static_cast<double>(i) / static_cast<double>(n);
        const double u = axis1 * std::cos(t);
        const double v = axis2 * std::sin(t);
        line->SetPoint(i, model.mean_x + u * ca - v * sa,
                       model.mean_y + u * sa + v * ca);
    }
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    return line;
}

} // namespace detail

inline Result apply_mass_cuts(const std::vector<Point>& input_points,
                              const Params& params,
                              const Config& cfg = Config())
{
    Result result;
    result.params = params;
    result.pass_ellipse.assign(input_points.size(), 0);
    result.pass_mcd.assign(input_points.size(), 0);
    if (!params.valid || params.cov_det <= 0.0 || params.mcd_det <= 0.0) return result;

    detail::CovModel model;
    model.mean_x = params.mean_mpi0;
    model.mean_y = params.mean_mmiss;
    model.cov_xx = params.cov_mpi0_mpi0;
    model.cov_xy = params.cov_mpi0_mmiss;
    model.cov_yy = params.cov_mmiss_mmiss;
    model.det = params.cov_det;

    detail::CovModel mcd_model;
    mcd_model.mean_x = params.mcd_mean_mpi0;
    mcd_model.mean_y = params.mcd_mean_mmiss;
    mcd_model.cov_xx = params.mcd_cov_mpi0_mpi0;
    mcd_model.cov_xy = params.mcd_cov_mpi0_mmiss;
    mcd_model.cov_yy = params.mcd_cov_mmiss_mmiss;
    mcd_model.det = params.mcd_det;

    double total_weight = 0.0;
    double ellipse_weight = 0.0;
    double mcd_weight = 0.0;
    for (std::size_t i = 0; i < input_points.size(); ++i) {
        const auto& p = input_points[i];
        if (!std::isfinite(p.mpi0) || !std::isfinite(p.mmiss) ||
            !std::isfinite(p.weight) || p.weight <= 0.0) continue;
        if (p.mpi0 < cfg.mpi0_min || p.mpi0 >= cfg.mpi0_max ||
            p.mmiss < cfg.mmiss_min || p.mmiss >= cfg.mmiss_max) continue;
        total_weight += p.weight;
        if (detail::covariance_d2(model, p.mpi0, p.mmiss) <= params.ellipse_d2_cut) {
            result.pass_ellipse[i] = 1;
            ellipse_weight += p.weight;
        }
        if (detail::covariance_d2(mcd_model, p.mpi0, p.mmiss) <= params.mcd_d2_cut) {
            result.pass_mcd[i] = 1;
            mcd_weight += p.weight;
        }
    }

    result.params.total_weight = total_weight;
    result.params.ellipse_weight = ellipse_weight;
    result.params.mcd_weight = mcd_weight;
    result.params.ellipse_total_fraction = (total_weight > 0.0) ? ellipse_weight / total_weight : 0.0;
    result.params.mcd_total_fraction = (total_weight > 0.0) ? mcd_weight / total_weight : 0.0;
    return result;
}

inline Result evaluate_mass_cuts(const std::vector<Point>& input_points, const Config& cfg) {
    using namespace detail;

    Result result;
    result.pass_ellipse.assign(input_points.size(), 0);
    result.pass_mcd.assign(input_points.size(), 0);

    if (input_points.empty()) return result;

    TH2D* h2 = new TH2D((cfg.tag + "_h_mmiss_vs_mpi0_weighted").c_str(),
                        "Weighted mmiss_all:mpi0_all;mpi0_all [GeV];mmiss_all [GeV]",
                        cfg.n_mpi0_bins, cfg.mpi0_min, cfg.mpi0_max,
                        cfg.n_mmiss_bins, cfg.mmiss_min, cfg.mmiss_max);
    h2->SetDirectory(nullptr);

    for (const auto& p : input_points) {
        if (!std::isfinite(p.mpi0) || !std::isfinite(p.mmiss) || !std::isfinite(p.weight)) continue;
        if (p.weight <= 0.0) continue;
        if (p.mpi0 < cfg.mpi0_min || p.mpi0 >= cfg.mpi0_max) continue;
        if (p.mmiss < cfg.mmiss_min || p.mmiss >= cfg.mmiss_max) continue;
        h2->Fill(p.mpi0, p.mmiss, p.weight);
    }

    const int nx = h2->GetNbinsX();
    const int ny = h2->GetNbinsY();
    const int n_cells = h2->GetNcells();
    const double total_weight = h2->Integral(1, nx, 1, ny);
    if (total_weight <= 0.0) {
        delete h2;
        return result;
    }

    std::vector<WeightedPoint> points;
    points.reserve(nx * ny);
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const double w = h2->GetBinContent(ix, iy);
            if (w <= 0.0) continue;
            points.push_back({h2->GetXaxis()->GetBinCenter(ix),
                              h2->GetYaxis()->GetBinCenter(iy),
                              w,
                              h2->GetBin(ix, iy)});
        }
    }

    int peak_global = h2->GetMaximumBin();
    int peak_ix = 0, peak_iy = 0, peak_iz = 0;
    h2->GetBinXYZ(peak_global, peak_ix, peak_iy, peak_iz);
    const double peak_weight = h2->GetBinContent(peak_global);
    const double peak_tolerance = std::max(1e-12, std::abs(peak_weight) * 1e-9);

    auto build_core = [&](double peak_fraction, std::vector<unsigned char>& core_mask) {
        core_mask.assign(n_cells + 1, 0);
        std::queue<std::pair<int, int> > frontier;
        const double threshold_weight = peak_fraction * peak_weight;

        auto push = [&](int ix, int iy) {
            if (ix < 1 || ix > nx || iy < 1 || iy > ny) return;
            const int global = h2->GetBin(ix, iy);
            if (core_mask[global]) return;
            const double w = h2->GetBinContent(global);
            if (w + 1e-12 < threshold_weight) return;
            core_mask[global] = 1;
            frontier.push(std::make_pair(ix, iy));
        };

        for (int ix = 1; ix <= nx; ++ix) {
            for (int iy = 1; iy <= ny; ++iy) {
                const double w = h2->GetBinContent(ix, iy);
                if (std::abs(w - peak_weight) <= peak_tolerance) push(ix, iy);
            }
        }

        while (!frontier.empty()) {
            const auto bin = frontier.front();
            frontier.pop();
            const int ix = bin.first;
            const int iy = bin.second;
            push(ix - 1, iy); push(ix + 1, iy); push(ix, iy - 1); push(ix, iy + 1);
            push(ix - 1, iy - 1); push(ix - 1, iy + 1);
            push(ix + 1, iy - 1); push(ix + 1, iy + 1);
        }

        CoreStats stats;
        stats.peak_fraction = peak_fraction;
        for (int ix = 1; ix <= nx; ++ix) {
            for (int iy = 1; iy <= ny; ++iy) {
                const int global = h2->GetBin(ix, iy);
                if (!core_mask[global]) continue;
                stats.weight += h2->GetBinContent(ix, iy);
                ++stats.bins;
            }
        }
        stats.total_fraction = stats.weight / total_weight;
        return stats;
    };

    std::vector<CoreStats> scan_stats;
    double peak_fraction = cfg.peak_fraction;
    double auto_jump_ratio = 1.0;
    const bool auto_peak_fraction = (peak_fraction <= 0.0);
    if (auto_peak_fraction) {
        std::vector<unsigned char> scan_mask;
        CoreStats previous;
        bool have_previous = false;
        double best_jump = 0.0;
        double chosen_peak_fraction = cfg.auto_peak_max;
        const int n_scan = static_cast<int>(std::floor((cfg.auto_peak_max - cfg.auto_peak_min) / cfg.auto_peak_step + 0.5)) + 1;
        for (int i = 0; i < n_scan; ++i) {
            double candidate = cfg.auto_peak_max - i * cfg.auto_peak_step;
            if (candidate < cfg.auto_peak_min) candidate = cfg.auto_peak_min;
            CoreStats current = build_core(candidate, scan_mask);
            scan_stats.push_back(current);
            if (have_previous && previous.bins > 0 && previous.weight > 0.0) {
                const double jump = std::max(current.weight / previous.weight,
                                             static_cast<double>(current.bins) / static_cast<double>(previous.bins));
                const bool previous_core_is_real =
                    previous.bins >= cfg.auto_min_core_bins &&
                    previous.total_fraction >= cfg.auto_min_core_total_fraction;
                if (previous_core_is_real && jump > best_jump) {
                    best_jump = jump;
                    chosen_peak_fraction = previous.peak_fraction;
                }
            }
            previous = current;
            have_previous = true;
        }
        peak_fraction = chosen_peak_fraction;
        auto_jump_ratio = best_jump;
    }

    std::vector<unsigned char> selected(n_cells + 1, 0);
    CoreStats core_stats = build_core(peak_fraction, selected);
    if (core_stats.bins < 3 || core_stats.weight <= 0.0) {
        delete h2;
        return result;
    }

    CovModel model;
    double model_weight = 0.0;
    int model_bins = 0;
    if (!compute_cov_model(points, selected, model, model_weight, model_bins)) {
        delete h2;
        return result;
    }

    std::vector<WeightedValue> core_d2_values;
    for (const auto& p : points) {
        if (!selected[p.global_bin]) continue;
        core_d2_values.push_back({covariance_d2(model, p.x, p.y), p.weight});
    }
    const double ellipse_d2_cut = weighted_quantile(core_d2_values, cfg.ellipse_core_quantile) *
                                  cfg.ellipse_padding * cfg.ellipse_padding;

    CovModel mcd_model = model;
    CovModel best_mcd_model = model;
    std::vector<unsigned char> current_mcd_subset = selected;
    std::vector<unsigned char> best_mcd_subset = selected;
    double best_mcd_det = std::numeric_limits<double>::infinity();
    double best_mcd_subset_weight = model_weight;
    int best_mcd_subset_bins = model_bins;
    const double mcd_target_weight = cfg.mcd_keep_total_fraction * total_weight;

    for (int iter = 0; iter < cfg.mcd_iterations; ++iter) {
        std::vector<IndexedDistance> distances;
        distances.reserve(points.size());
        for (std::size_t ip = 0; ip < points.size(); ++ip) {
            distances.push_back({covariance_d2(mcd_model, points[ip].x, points[ip].y),
                                 points[ip].weight, ip});
        }
        std::sort(distances.begin(), distances.end(),
                  [](const IndexedDistance& a, const IndexedDistance& b) { return a.d2 < b.d2; });

        current_mcd_subset.assign(n_cells + 1, 0);
        double kept_weight = 0.0;
        for (const auto& d : distances) {
            const auto& p = points[d.index];
            current_mcd_subset[p.global_bin] = 1;
            kept_weight += p.weight;
            if (kept_weight >= mcd_target_weight) break;
        }

        CovModel candidate_model;
        double candidate_weight = 0.0;
        int candidate_bins = 0;
        if (!compute_cov_model(points, current_mcd_subset, candidate_model, candidate_weight, candidate_bins)) break;
        mcd_model = candidate_model;
        if (candidate_model.det < best_mcd_det) {
            best_mcd_det = candidate_model.det;
            best_mcd_model = candidate_model;
            best_mcd_subset = current_mcd_subset;
            best_mcd_subset_weight = candidate_weight;
            best_mcd_subset_bins = candidate_bins;
        }
    }
    if (!std::isfinite(best_mcd_det)) {
        best_mcd_det = model.det;
        best_mcd_model = model;
        best_mcd_subset = selected;
        best_mcd_subset_weight = model_weight;
        best_mcd_subset_bins = model_bins;
    }

    std::vector<WeightedValue> mcd_d2_values;
    for (const auto& p : points) {
        if (!best_mcd_subset[p.global_bin]) continue;
        mcd_d2_values.push_back({covariance_d2(best_mcd_model, p.x, p.y), p.weight});
    }
    const double mcd_d2_cut = weighted_quantile(mcd_d2_values, cfg.mcd_ellipse_quantile) *
                              cfg.mcd_padding * cfg.mcd_padding;

    double ellipse_weight = 0.0;
    double mcd_weight = 0.0;
    for (std::size_t i = 0; i < input_points.size(); ++i) {
        const auto& p = input_points[i];
        if (!std::isfinite(p.mpi0) || !std::isfinite(p.mmiss) || !std::isfinite(p.weight) || p.weight <= 0.0) continue;
        if (p.mpi0 < cfg.mpi0_min || p.mpi0 >= cfg.mpi0_max || p.mmiss < cfg.mmiss_min || p.mmiss >= cfg.mmiss_max) continue;
        const double d2 = covariance_d2(model, p.mpi0, p.mmiss);
        const double mcd_d2 = covariance_d2(best_mcd_model, p.mpi0, p.mmiss);
        if (d2 <= ellipse_d2_cut) {
            result.pass_ellipse[i] = 1;
            ellipse_weight += p.weight;
        }
        if (mcd_d2 <= mcd_d2_cut) {
            result.pass_mcd[i] = 1;
            mcd_weight += p.weight;
        }
    }

    result.params.valid = true;
    result.params.auto_peak_fraction = auto_peak_fraction;
    result.params.peak_fraction = peak_fraction;
    result.params.auto_jump_ratio = auto_jump_ratio;
    result.params.peak_weight = peak_weight;
    result.params.threshold_weight = peak_fraction * peak_weight;
    result.params.total_weight = total_weight;
    result.params.peak_ix = peak_ix;
    result.params.peak_iy = peak_iy;
    result.params.peak_mpi0 = h2->GetXaxis()->GetBinCenter(peak_ix);
    result.params.peak_mmiss = h2->GetYaxis()->GetBinCenter(peak_iy);
    result.params.core_bins = core_stats.bins;
    result.params.core_weight = core_stats.weight;
    result.params.core_total_fraction = core_stats.total_fraction;
    result.params.mean_mpi0 = model.mean_x;
    result.params.mean_mmiss = model.mean_y;
    result.params.cov_mpi0_mpi0 = model.cov_xx;
    result.params.cov_mpi0_mmiss = model.cov_xy;
    result.params.cov_mmiss_mmiss = model.cov_yy;
    result.params.cov_det = model.det;
    result.params.ellipse_d2_cut = ellipse_d2_cut;
    result.params.ellipse_weight = ellipse_weight;
    result.params.ellipse_total_fraction = ellipse_weight / total_weight;
    result.params.mcd_mean_mpi0 = best_mcd_model.mean_x;
    result.params.mcd_mean_mmiss = best_mcd_model.mean_y;
    result.params.mcd_cov_mpi0_mpi0 = best_mcd_model.cov_xx;
    result.params.mcd_cov_mpi0_mmiss = best_mcd_model.cov_xy;
    result.params.mcd_cov_mmiss_mmiss = best_mcd_model.cov_yy;
    result.params.mcd_det = best_mcd_model.det;
    result.params.mcd_d2_cut = mcd_d2_cut;
    result.params.mcd_subset_weight = best_mcd_subset_weight;
    result.params.mcd_subset_bins = best_mcd_subset_bins;
    result.params.mcd_weight = mcd_weight;
    result.params.mcd_total_fraction = mcd_weight / total_weight;

    if (cfg.write_debug) {
        gSystem->mkdir(cfg.output_dir.c_str(), kTRUE);
        const std::string base = cfg.output_dir + "/" + cfg.tag;
        TFile* fdebug = TFile::Open((base + ".root").c_str(), "RECREATE");
        TH2D* h_ellipse = dynamic_cast<TH2D*>(h2->Clone((cfg.tag + "_ellipse_selected").c_str()));
        TH2D* h_mcd = dynamic_cast<TH2D*>(h2->Clone((cfg.tag + "_mcd_selected").c_str()));
        for (int ix = 1; ix <= nx; ++ix) {
            for (int iy = 1; iy <= ny; ++iy) {
                const double x = h2->GetXaxis()->GetBinCenter(ix);
                const double y = h2->GetYaxis()->GetBinCenter(iy);
                if (covariance_d2(model, x, y) > ellipse_d2_cut) h_ellipse->SetBinContent(ix, iy, 0.0);
                if (covariance_d2(best_mcd_model, x, y) > mcd_d2_cut) h_mcd->SetBinContent(ix, iy, 0.0);
            }
        }

        TCanvas* c = new TCanvas((cfg.tag + "_canvas").c_str(), "2D mass cut debug", 1600, 1200);
        c->Divide(2, 2);
        TPolyLine* ellipse_line = make_ellipse_line(model, ellipse_d2_cut, kMagenta + 2, 1, 3);
        TPolyLine* mcd_line = make_ellipse_line(best_mcd_model, mcd_d2_cut, kGreen + 2, 1, 3);
        TPolyLine* ellipse_line2 = make_ellipse_line(model, ellipse_d2_cut, kMagenta + 2, 1, 3);
        TPolyLine* mcd_line2 = make_ellipse_line(best_mcd_model, mcd_d2_cut, kGreen + 2, 1, 3);
        c->cd(1); gPad->SetRightMargin(0.14); h2->Draw("COLZ"); ellipse_line->Draw("SAME"); mcd_line->Draw("SAME");
        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.030);
        latex.DrawLatex(0.13, 0.93, Form("peak_fraction %.3f", peak_fraction));
        latex.DrawLatex(0.13, 0.88, Form("ellipse %.2f%%, MCD %.2f%%",
                                         100.0 * result.params.ellipse_total_fraction,
                                         100.0 * result.params.mcd_total_fraction));
        c->cd(2); gPad->SetRightMargin(0.14); h_mcd->Draw("COLZ"); ellipse_line2->Draw("SAME"); mcd_line2->Draw("SAME");
        c->cd(3);
        TH1D* hx_all = h2->ProjectionX((cfg.tag + "_mpi0_all").c_str(), 1, ny);
        TH1D* hx_ellipse = h_ellipse->ProjectionX((cfg.tag + "_mpi0_ellipse").c_str(), 1, ny);
        TH1D* hx_mcd = h_mcd->ProjectionX((cfg.tag + "_mpi0_mcd").c_str(), 1, ny);
        hx_all->SetLineColor(kBlack); hx_ellipse->SetLineColor(kMagenta + 2); hx_mcd->SetLineColor(kGreen + 2);
        hx_all->Draw("HIST"); hx_ellipse->Draw("HIST SAME"); hx_mcd->Draw("HIST SAME");
        TLegend* legx = new TLegend(0.58, 0.72, 0.88, 0.88);
        legx->AddEntry(hx_all, "all", "l"); legx->AddEntry(hx_ellipse, "ellipse", "l"); legx->AddEntry(hx_mcd, "MCD", "l"); legx->Draw();
        c->cd(4);
        TH1D* hy_all = h2->ProjectionY((cfg.tag + "_mmiss_all").c_str(), 1, nx);
        TH1D* hy_ellipse = h_ellipse->ProjectionY((cfg.tag + "_mmiss_ellipse").c_str(), 1, nx);
        TH1D* hy_mcd = h_mcd->ProjectionY((cfg.tag + "_mmiss_mcd").c_str(), 1, nx);
        hy_all->SetLineColor(kBlack); hy_ellipse->SetLineColor(kMagenta + 2); hy_mcd->SetLineColor(kGreen + 2);
        hy_all->Draw("HIST"); hy_ellipse->Draw("HIST SAME"); hy_mcd->Draw("HIST SAME");
        TLegend* legy = new TLegend(0.58, 0.72, 0.88, 0.88);
        legy->AddEntry(hy_all, "all", "l"); legy->AddEntry(hy_ellipse, "ellipse", "l"); legy->AddEntry(hy_mcd, "MCD", "l"); legy->Draw();
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());

        std::ofstream params((base + "_params.csv").c_str());
        params << "parameter,value\n";
        params << "peak_fraction," << result.params.peak_fraction << "\n";
        params << "auto_peak_fraction," << (result.params.auto_peak_fraction ? 1 : 0) << "\n";
        params << "auto_jump_ratio," << result.params.auto_jump_ratio << "\n";
        params << "auto_min_core_total_fraction," << cfg.auto_min_core_total_fraction << "\n";
        params << "auto_min_core_bins," << cfg.auto_min_core_bins << "\n";
        params << "mean_mpi0," << result.params.mean_mpi0 << "\n";
        params << "mean_mmiss," << result.params.mean_mmiss << "\n";
        params << "cov_mpi0_mpi0," << result.params.cov_mpi0_mpi0 << "\n";
        params << "cov_mpi0_mmiss," << result.params.cov_mpi0_mmiss << "\n";
        params << "cov_mmiss_mmiss," << result.params.cov_mmiss_mmiss << "\n";
        params << "ellipse_d2_cut," << result.params.ellipse_d2_cut << "\n";
        params << "ellipse_total_fraction," << result.params.ellipse_total_fraction << "\n";
        params << "mcd_mean_mpi0," << result.params.mcd_mean_mpi0 << "\n";
        params << "mcd_mean_mmiss," << result.params.mcd_mean_mmiss << "\n";
        params << "mcd_cov_mpi0_mpi0," << result.params.mcd_cov_mpi0_mpi0 << "\n";
        params << "mcd_cov_mpi0_mmiss," << result.params.mcd_cov_mpi0_mmiss << "\n";
        params << "mcd_cov_mmiss_mmiss," << result.params.mcd_cov_mmiss_mmiss << "\n";
        params << "mcd_d2_cut," << result.params.mcd_d2_cut << "\n";
        params << "mcd_total_fraction," << result.params.mcd_total_fraction << "\n";
        params.close();

        std::ofstream scan((base + "_peak_scan.csv").c_str());
        scan << "peak_fraction,selected_bins,selected_weight,selected_total_fraction,weight_ratio_to_previous,bin_ratio_to_previous\n";
        for (std::size_t i = 0; i < scan_stats.size(); ++i) {
            double wr = 0.0, br = 0.0;
            if (i > 0 && scan_stats[i - 1].weight > 0.0 && scan_stats[i - 1].bins > 0) {
                wr = scan_stats[i].weight / scan_stats[i - 1].weight;
                br = static_cast<double>(scan_stats[i].bins) / static_cast<double>(scan_stats[i - 1].bins);
            }
            scan << scan_stats[i].peak_fraction << "," << scan_stats[i].bins << ","
                 << scan_stats[i].weight << "," << scan_stats[i].total_fraction << ","
                 << wr << "," << br << "\n";
        }
        scan.close();

        if (fdebug && !fdebug->IsZombie()) {
            fdebug->cd();
            h2->Write("h_mmiss_vs_mpi0_weighted");
            h_ellipse->Write("h_mmiss_vs_mpi0_ellipse");
            h_mcd->Write("h_mmiss_vs_mpi0_mcd");
            c->Write("c_2d_mass_cut_debug");
            fdebug->Close();
        }
        delete c;
        if (fdebug) delete fdebug;
    }

    delete h2;
    result.params.valid = true;
    return result;
}

} // namespace nps2d

#endif
