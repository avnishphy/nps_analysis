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
    double seed_mpi0 = 0.13498;
    double seed_mpi0_half_width = 0.010;
    double seed_mmiss = 0.938;
    double seed_mmiss_half_width = 0.030;
    int smoothing_radius_bins = 1;
    double core_quantile = 0.6827;
    double ellipse_quantile = 0.990;
    int ellipse_grow_iterations = 20;
    double ellipse_padding = 1.05;
    double auto_peak_min = 0.05;
    double auto_peak_max = 0.60;
    double auto_peak_step = 0.005;
    double auto_min_core_total_fraction = 0.005;
    double auto_max_core_total_fraction = 0.30;
    int auto_min_core_bins = 8;
    double max_model_mpi0_offset = 0.008;
    double max_model_mmiss_offset = 0.100;
    double mcd_candidate_mpi0_half_width = 0.015;
    double mcd_candidate_mmiss_half_width = 0.180;
    double mcd_keep_candidate_fraction = 0.50;
    int mcd_iterations = 8;
    double mcd_ellipse_quantile = 0.975;
    double mcd_padding = 1.05;
    double covariance_regularization_bins = 0.50;
    bool write_debug = true;
    std::string output_dir = ".";
    std::string tag = "mass_cut";
};

struct Params {
    bool valid = false;
    bool ellipse_valid = false;
    bool mcd_valid = false;
    bool auto_peak_fraction = false;
    bool core_leak_rejected = false;
    double peak_fraction = 0.0;
    double auto_jump_ratio = 1.0;
    double peak_weight = 0.0;
    double threshold_weight = 0.0;
    double total_weight = 0.0;
    int peak_ix = 0;
    int peak_iy = 0;
    double peak_mpi0 = 0.0;
    double peak_mmiss = 0.0;
    double smoothed_peak_weight = 0.0;
    int core_bins = 0;
    double core_weight = 0.0;
    double core_total_fraction = 0.0;
    int fit_subset_bins = 0;
    double fit_subset_weight = 0.0;
    double fit_subset_total_fraction = 0.0;
    int seed_core_bins = 0;
    double seed_core_weight = 0.0;
    double seed_core_total_fraction = 0.0;
    int ellipse_growth_steps = 0;

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
    double mean_x = 0.0;
    double mean_y = 0.0;
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
                              int& n_bins,
                              double regularization_x = 0.0,
                              double regularization_y = 0.0)
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
    model.cov_xx += regularization_x * regularization_x;
    model.cov_yy += regularization_y * regularization_y;
    const double max_abs_cov_xy = 0.995 * std::sqrt(model.cov_xx * model.cov_yy);
    model.cov_xy = std::max(-max_abs_cov_xy, std::min(max_abs_cov_xy, model.cov_xy));
    model.det = model.cov_xx * model.cov_yy - model.cov_xy * model.cov_xy;
    return model.det > 0.0;
}

inline bool model_is_anchored(const CovModel& model, const Config& cfg) {
    return std::abs(model.mean_x - cfg.seed_mpi0) <= cfg.max_model_mpi0_offset &&
           std::abs(model.mean_y - cfg.seed_mmiss) <= cfg.max_model_mmiss_offset;
}

inline double chi2_quantile_df2(double quantile) {
    const double q = std::max(0.0, std::min(1.0 - 1e-12, quantile));
    return -2.0 * std::log(1.0 - q);
}

inline double mcd_consistency_scale_df2(double keep_fraction) {
    const double h = std::max(1e-6, std::min(1.0 - 1e-6, keep_fraction));
    const double cutoff = chi2_quantile_df2(h);
    const double truncated_variance = 1.0 - cutoff * (1.0 - h) / (2.0 * h);
    return (truncated_variance > 1e-6) ? 1.0 / truncated_variance : 1.0;
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
    if (!params.valid) return result;

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
        if (params.ellipse_valid && params.cov_det > 0.0 &&
            detail::covariance_d2(model, p.mpi0, p.mmiss) <= params.ellipse_d2_cut) {
            result.pass_ellipse[i] = 1;
            ellipse_weight += p.weight;
        }
        if (params.mcd_valid && params.mcd_det > 0.0 &&
            detail::covariance_d2(mcd_model, p.mpi0, p.mmiss) <= params.mcd_d2_cut) {
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

    // Smooth only the density used for peak finding and connectivity. Final
    // covariance and accepted weights always use the unsmoothed histogram.
    std::vector<double> density(n_cells + 1, 0.0);
    const int smooth_radius = std::max(0, cfg.smoothing_radius_bins);
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double weighted_sum = 0.0;
            double kernel_sum = 0.0;
            for (int dx = -smooth_radius; dx <= smooth_radius; ++dx) {
                for (int dy = -smooth_radius; dy <= smooth_radius; ++dy) {
                    const int jx = ix + dx;
                    const int jy = iy + dy;
                    if (jx < 1 || jx > nx || jy < 1 || jy > ny) continue;
                    const double kernel_weight =
                        static_cast<double>((smooth_radius + 1 - std::abs(dx)) *
                                            (smooth_radius + 1 - std::abs(dy)));
                    weighted_sum += kernel_weight * h2->GetBinContent(jx, jy);
                    kernel_sum += kernel_weight;
                }
            }
            density[h2->GetBin(ix, iy)] = (kernel_sum > 0.0) ? weighted_sum / kernel_sum : 0.0;
        }
    }

    int peak_ix = 0;
    int peak_iy = 0;
    double smoothed_peak_weight = 0.0;
    double best_anchor_distance = std::numeric_limits<double>::infinity();
    for (int ix = 1; ix <= nx; ++ix) {
        const double x = h2->GetXaxis()->GetBinCenter(ix);
        if (std::abs(x - cfg.seed_mpi0) > cfg.seed_mpi0_half_width) continue;
        for (int iy = 1; iy <= ny; ++iy) {
            const double y = h2->GetYaxis()->GetBinCenter(iy);
            if (std::abs(y - cfg.seed_mmiss) > cfg.seed_mmiss_half_width) continue;
            const double value = density[h2->GetBin(ix, iy)];
            const double anchor_distance =
                std::pow((x - cfg.seed_mpi0) / cfg.seed_mpi0_half_width, 2) +
                std::pow((y - cfg.seed_mmiss) / cfg.seed_mmiss_half_width, 2);
            const double tolerance = std::max(1e-12, std::abs(smoothed_peak_weight) * 1e-9);
            if (value > smoothed_peak_weight + tolerance ||
                (std::abs(value - smoothed_peak_weight) <= tolerance &&
                 anchor_distance < best_anchor_distance)) {
                smoothed_peak_weight = value;
                best_anchor_distance = anchor_distance;
                peak_ix = ix;
                peak_iy = iy;
            }
        }
    }
    if (peak_ix == 0 || peak_iy == 0 || smoothed_peak_weight <= 0.0) {
        delete h2;
        return result;
    }
    const double peak_weight = h2->GetBinContent(peak_ix, peak_iy);

    auto build_core = [&](double peak_fraction, std::vector<unsigned char>& core_mask) {
        core_mask.assign(n_cells + 1, 0);
        std::queue<std::pair<int, int> > frontier;
        const double threshold_weight = peak_fraction * smoothed_peak_weight;

        auto push = [&](int ix, int iy) {
            if (ix < 1 || ix > nx || iy < 1 || iy > ny) return;
            const int global = h2->GetBin(ix, iy);
            if (core_mask[global]) return;
            const double w = density[global];
            if (w + 1e-12 < threshold_weight) return;
            core_mask[global] = 1;
            frontier.push(std::make_pair(ix, iy));
        };

        // One anchored seed only. Seeding every tied maximum makes sparse,
        // disconnected background islands part of the same "core".
        push(peak_ix, peak_iy);

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
                const double w = h2->GetBinContent(ix, iy);
                if (w <= 0.0) continue;
                stats.weight += w;
                stats.mean_x += w * h2->GetXaxis()->GetBinCenter(ix);
                stats.mean_y += w * h2->GetYaxis()->GetBinCenter(iy);
                ++stats.bins;
            }
        }
        stats.total_fraction = stats.weight / total_weight;
        if (stats.weight > 0.0) {
            stats.mean_x /= stats.weight;
            stats.mean_y /= stats.weight;
        }
        return stats;
    };

    std::vector<CoreStats> scan_stats;
    double peak_fraction = cfg.peak_fraction;
    double auto_jump_ratio = 1.0;
    const bool auto_peak_fraction = (peak_fraction <= 0.0);
    std::vector<unsigned char> selected(n_cells + 1, 0);
    if (auto_peak_fraction) {
        std::vector<unsigned char> scan_mask;
        std::vector<unsigned char> last_safe_mask;
        std::vector<unsigned char> last_qualified_mask;
        CoreStats previous;
        bool have_previous = false;
        bool have_safe = false;
        bool have_qualified = false;
        double best_jump = 0.0;
        double chosen_peak_fraction = 0.0;
        const int n_scan = static_cast<int>(std::floor((cfg.auto_peak_max - cfg.auto_peak_min) / cfg.auto_peak_step + 0.5)) + 1;
        for (int i = 0; i < n_scan; ++i) {
            double candidate = cfg.auto_peak_max - i * cfg.auto_peak_step;
            if (candidate < cfg.auto_peak_min) candidate = cfg.auto_peak_min;
            CoreStats current = build_core(candidate, scan_mask);
            scan_stats.push_back(current);
            const bool anchored =
                std::abs(current.mean_x - cfg.seed_mpi0) <= cfg.max_model_mpi0_offset &&
                std::abs(current.mean_y - cfg.seed_mmiss) <= cfg.max_model_mmiss_offset;
            const bool safe = current.bins >= 3 && current.weight > 0.0 && anchored &&
                              current.total_fraction <= cfg.auto_max_core_total_fraction;
            const bool qualified = safe && current.bins >= cfg.auto_min_core_bins &&
                                   current.total_fraction >= cfg.auto_min_core_total_fraction;
            if (have_previous && previous.bins > 0 && previous.weight > 0.0) {
                const double jump = std::max(current.weight / previous.weight,
                                             static_cast<double>(current.bins) / static_cast<double>(previous.bins));
                best_jump = std::max(best_jump, jump);
            }
            if (safe) {
                have_safe = true;
                last_safe_mask = scan_mask;
                chosen_peak_fraction = current.peak_fraction;
                if (qualified) {
                    have_qualified = true;
                    last_qualified_mask = scan_mask;
                }
            } else if (have_safe &&
                       (current.total_fraction > cfg.auto_max_core_total_fraction || !anchored)) {
                result.params.core_leak_rejected = true;
                break;
            }
            previous = current;
            have_previous = true;
        }
        if (have_qualified) {
            selected = last_qualified_mask;
        } else if (have_safe) {
            selected = last_safe_mask;
        } else {
            delete h2;
            return result;
        }
        peak_fraction = chosen_peak_fraction;
        auto_jump_ratio = best_jump;
    }

    CoreStats core_stats;
    if (auto_peak_fraction) {
        core_stats = build_core(peak_fraction, selected);
    } else {
        core_stats = build_core(peak_fraction, selected);
        const bool safe = core_stats.bins >= 3 && core_stats.weight > 0.0 &&
            core_stats.total_fraction <= cfg.auto_max_core_total_fraction &&
            std::abs(core_stats.mean_x - cfg.seed_mpi0) <= cfg.max_model_mpi0_offset &&
            std::abs(core_stats.mean_y - cfg.seed_mmiss) <= cfg.max_model_mmiss_offset;
        if (!safe) {
            delete h2;
            return result;
        }
    }
    if (core_stats.bins < 3 || core_stats.weight <= 0.0) {
        delete h2;
        return result;
    }

    CovModel model;
    double model_weight = 0.0;
    int model_bins = 0;
    const double regularization_x = cfg.covariance_regularization_bins * h2->GetXaxis()->GetBinWidth(1);
    const double regularization_y = cfg.covariance_regularization_bins * h2->GetYaxis()->GetBinWidth(1);
    if (!compute_cov_model(points, selected, model, model_weight, model_bins,
                           regularization_x, regularization_y) ||
        !model_is_anchored(model, cfg)) {
        delete h2;
        return result;
    }
    const CovModel seed_core_model = model;
    const CoreStats seed_core_stats = core_stats;

    // Grow the tiny density seed into the local, anti-correlated event band by
    // iterative Mahalanobis inclusion. No rectangular gate is applied: the
    // ellipse itself defines smooth membership. Anchor and weight checks stop
    // growth before it reaches the remote global background.
    const double ellipse_grow_d2 = chi2_quantile_df2(cfg.ellipse_quantile);
    int ellipse_growth_steps = 0;
    for (int iter = 0; iter < cfg.ellipse_grow_iterations && points.size() >= 3; ++iter) {
        std::vector<unsigned char> grown_subset(n_cells + 1, 0);
        for (const auto& p : points) {
            if (covariance_d2(model, p.x, p.y) <= ellipse_grow_d2)
                grown_subset[p.global_bin] = 1;
        }
        CovModel grown_model;
        double grown_weight = 0.0;
        int grown_bins = 0;
        if (!compute_cov_model(points, grown_subset, grown_model, grown_weight, grown_bins,
                               regularization_x, regularization_y) ||
            !model_is_anchored(grown_model, cfg) ||
            grown_weight / total_weight > cfg.auto_max_core_total_fraction) break;
        const bool converged = (grown_subset == selected);
        selected = grown_subset;
        model = grown_model;
        model_weight = grown_weight;
        model_bins = grown_bins;
        ++ellipse_growth_steps;
        if (converged) break;
    }
    // The iteratively grown population is a fitting subset, not a physical
    // core cut. Its stopping edge is controlled by the contamination guard.
    CoreStats fit_subset_stats;
    fit_subset_stats.weight = model_weight;
    fit_subset_stats.total_fraction = model_weight / total_weight;
    fit_subset_stats.bins = model_bins;
    fit_subset_stats.mean_x = model.mean_x;
    fit_subset_stats.mean_y = model.mean_y;

    // Define the reported core from the final covariance. This makes it a
    // smooth, nested inner contour rather than exposing the fitting subset's
    // algorithmic stopping boundary.
    const double core_d2_cut = chi2_quantile_df2(cfg.core_quantile);
    core_stats = CoreStats{};
    for (const auto& p : points) {
        if (covariance_d2(model, p.x, p.y) > core_d2_cut) continue;
        core_stats.weight += p.weight;
        ++core_stats.bins;
    }
    core_stats.total_fraction = core_stats.weight / total_weight;
    core_stats.mean_x = model.mean_x;
    core_stats.mean_y = model.mean_y;

    const double ellipse_d2_cut = chi2_quantile_df2(cfg.ellipse_quantile) *
                                  cfg.ellipse_padding * cfg.ellipse_padding;

    CovModel mcd_model = seed_core_model;
    CovModel best_mcd_model = seed_core_model;
    std::vector<unsigned char> current_mcd_subset(n_cells + 1, 0);
    std::vector<unsigned char> best_mcd_subset(n_cells + 1, 0);
    double best_mcd_det = std::numeric_limits<double>::infinity();
    double best_mcd_subset_weight = 0.0;
    int best_mcd_subset_bins = 0;

    std::vector<unsigned char> mcd_candidates(n_cells + 1, 0);
    double mcd_candidate_weight = 0.0;
    int mcd_candidate_bins = 0;
    for (const auto& p : points) {
        if (std::abs(p.x - cfg.seed_mpi0) > cfg.mcd_candidate_mpi0_half_width ||
            std::abs(p.y - cfg.seed_mmiss) > cfg.mcd_candidate_mmiss_half_width) continue;
        mcd_candidates[p.global_bin] = 1;
        mcd_candidate_weight += p.weight;
        ++mcd_candidate_bins;
    }
    const double mcd_keep_fraction =
        std::max(0.05, std::min(0.95, cfg.mcd_keep_candidate_fraction));
    const double mcd_target_weight = mcd_keep_fraction * mcd_candidate_weight;

    for (int iter = 0; iter < cfg.mcd_iterations && mcd_candidate_bins >= 3; ++iter) {
        std::vector<IndexedDistance> distances;
        distances.reserve(mcd_candidate_bins);
        for (std::size_t ip = 0; ip < points.size(); ++ip) {
            if (!mcd_candidates[points[ip].global_bin]) continue;
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
        if (!compute_cov_model(points, current_mcd_subset, candidate_model,
                               candidate_weight, candidate_bins,
                               regularization_x, regularization_y) ||
            !model_is_anchored(candidate_model, cfg)) break;
        mcd_model = candidate_model;
        if (candidate_model.det < best_mcd_det) {
            best_mcd_det = candidate_model.det;
            best_mcd_model = candidate_model;
            best_mcd_subset = current_mcd_subset;
            best_mcd_subset_weight = candidate_weight;
            best_mcd_subset_bins = candidate_bins;
        }
    }
    const bool mcd_valid = std::isfinite(best_mcd_det) && best_mcd_subset_bins >= 3;
    if (mcd_valid) {
        // MCD covariance comes from a central truncated sample. Correct its
        // Gaussian scale before using a chi-square contour for final flags.
        const double consistency_scale = mcd_consistency_scale_df2(mcd_keep_fraction);
        best_mcd_model.cov_xx *= consistency_scale;
        best_mcd_model.cov_xy *= consistency_scale;
        best_mcd_model.cov_yy *= consistency_scale;
        best_mcd_model.det = best_mcd_model.cov_xx * best_mcd_model.cov_yy -
                             best_mcd_model.cov_xy * best_mcd_model.cov_xy;
    }
    const double mcd_d2_cut = chi2_quantile_df2(cfg.mcd_ellipse_quantile) *
                              cfg.mcd_padding * cfg.mcd_padding;

    double ellipse_weight = 0.0;
    double mcd_weight = 0.0;
    for (std::size_t i = 0; i < input_points.size(); ++i) {
        const auto& p = input_points[i];
        if (!std::isfinite(p.mpi0) || !std::isfinite(p.mmiss) || !std::isfinite(p.weight) || p.weight <= 0.0) continue;
        if (p.mpi0 < cfg.mpi0_min || p.mpi0 >= cfg.mpi0_max || p.mmiss < cfg.mmiss_min || p.mmiss >= cfg.mmiss_max) continue;
        const double d2 = covariance_d2(model, p.mpi0, p.mmiss);
        if (d2 <= ellipse_d2_cut) {
            result.pass_ellipse[i] = 1;
            ellipse_weight += p.weight;
        }
        if (mcd_valid && covariance_d2(best_mcd_model, p.mpi0, p.mmiss) <= mcd_d2_cut) {
            result.pass_mcd[i] = 1;
            mcd_weight += p.weight;
        }
    }

    result.params.valid = true;
    result.params.ellipse_valid = true;
    result.params.mcd_valid = mcd_valid;
    result.params.auto_peak_fraction = auto_peak_fraction;
    result.params.peak_fraction = peak_fraction;
    result.params.auto_jump_ratio = auto_jump_ratio;
    result.params.peak_weight = peak_weight;
    result.params.smoothed_peak_weight = smoothed_peak_weight;
    result.params.threshold_weight = peak_fraction * smoothed_peak_weight;
    result.params.total_weight = total_weight;
    result.params.peak_ix = peak_ix;
    result.params.peak_iy = peak_iy;
    result.params.peak_mpi0 = h2->GetXaxis()->GetBinCenter(peak_ix);
    result.params.peak_mmiss = h2->GetYaxis()->GetBinCenter(peak_iy);
    result.params.core_bins = core_stats.bins;
    result.params.core_weight = core_stats.weight;
    result.params.core_total_fraction = core_stats.total_fraction;
    result.params.fit_subset_bins = fit_subset_stats.bins;
    result.params.fit_subset_weight = fit_subset_stats.weight;
    result.params.fit_subset_total_fraction = fit_subset_stats.total_fraction;
    result.params.seed_core_bins = seed_core_stats.bins;
    result.params.seed_core_weight = seed_core_stats.weight;
    result.params.seed_core_total_fraction = seed_core_stats.total_fraction;
    result.params.ellipse_growth_steps = ellipse_growth_steps;
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
                if (!mcd_valid || covariance_d2(best_mcd_model, x, y) > mcd_d2_cut)
                    h_mcd->SetBinContent(ix, iy, 0.0);
            }
        }

        TCanvas* c = new TCanvas((cfg.tag + "_canvas").c_str(), "2D mass cut debug", 1600, 1200);
        c->Divide(2, 2);
        TPolyLine* ellipse_line = make_ellipse_line(model, ellipse_d2_cut, kMagenta + 2, 1, 3);
        TPolyLine* mcd_line = mcd_valid ? make_ellipse_line(best_mcd_model, mcd_d2_cut, kGreen + 2, 1, 3) : nullptr;
        TPolyLine* ellipse_line2 = make_ellipse_line(model, ellipse_d2_cut, kMagenta + 2, 1, 3);
        TPolyLine* mcd_line2 = mcd_valid ? make_ellipse_line(best_mcd_model, mcd_d2_cut, kGreen + 2, 1, 3) : nullptr;
        c->cd(1); gPad->SetRightMargin(0.14); h2->Draw("COLZ"); ellipse_line->Draw("SAME"); if (mcd_line) mcd_line->Draw("SAME");
        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.030);
        latex.DrawLatex(0.13, 0.93, Form("peak_fraction %.3f", peak_fraction));
        latex.DrawLatex(0.13, 0.88, Form("ellipse %.2f%%, MCD %.2f%%",
                                         100.0 * result.params.ellipse_total_fraction,
                                         100.0 * result.params.mcd_total_fraction));
        c->cd(2); gPad->SetRightMargin(0.14); h_mcd->Draw("COLZ"); ellipse_line2->Draw("SAME"); if (mcd_line2) mcd_line2->Draw("SAME");
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
        params << "ellipse_valid," << (result.params.ellipse_valid ? 1 : 0) << "\n";
        params << "mcd_valid," << (result.params.mcd_valid ? 1 : 0) << "\n";
        params << "core_leak_rejected," << (result.params.core_leak_rejected ? 1 : 0) << "\n";
        params << "peak_fraction," << result.params.peak_fraction << "\n";
        params << "auto_peak_fraction," << (result.params.auto_peak_fraction ? 1 : 0) << "\n";
        params << "auto_jump_ratio," << result.params.auto_jump_ratio << "\n";
        params << "auto_min_core_total_fraction," << cfg.auto_min_core_total_fraction << "\n";
        params << "auto_min_core_bins," << cfg.auto_min_core_bins << "\n";
        params << "auto_max_core_total_fraction," << cfg.auto_max_core_total_fraction << "\n";
        params << "core_quantile," << cfg.core_quantile << "\n";
        params << "seed_mpi0," << cfg.seed_mpi0 << "\n";
        params << "seed_mpi0_half_width," << cfg.seed_mpi0_half_width << "\n";
        params << "seed_mmiss," << cfg.seed_mmiss << "\n";
        params << "seed_mmiss_half_width," << cfg.seed_mmiss_half_width << "\n";
        params << "peak_mpi0," << result.params.peak_mpi0 << "\n";
        params << "peak_mmiss," << result.params.peak_mmiss << "\n";
        params << "peak_weight," << result.params.peak_weight << "\n";
        params << "smoothed_peak_weight," << result.params.smoothed_peak_weight << "\n";
        params << "core_bins," << result.params.core_bins << "\n";
        params << "core_total_fraction," << result.params.core_total_fraction << "\n";
        params << "fit_subset_bins," << result.params.fit_subset_bins << "\n";
        params << "fit_subset_total_fraction," << result.params.fit_subset_total_fraction << "\n";
        params << "seed_core_bins," << result.params.seed_core_bins << "\n";
        params << "seed_core_total_fraction," << result.params.seed_core_total_fraction << "\n";
        params << "ellipse_growth_steps," << result.params.ellipse_growth_steps << "\n";
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
        params << "mcd_keep_candidate_fraction," << cfg.mcd_keep_candidate_fraction << "\n";
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
