// Test 2D mass cut on weighted mmiss_all:mpi0_all.
//
// Run:
//   root -l -b -q 'scripts/test_2D_mass_cut.C()'
//
// Method:
//   1. Fill TH2D(mpi0_all, mmiss_all) with pi0_weight.
//   2. Find highest-count bin.
//   3. Select connected region around the peak where each bin has at least
//      peak_fraction * peak_bin_weight.
//   4. Fit a covariance ellipse to that core and use it as a smooth geometric cut.
//   5. Run a weighted minimum-covariance-distance (MCD) iteration as a robust comparison.
//   6. Save masks, selected-only histograms, and diagnostic plots.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TTree.h"
#include "TString.h"

namespace {

struct WeightedValue {
    double value;
    double weight;
};

struct CoreStats {
    double peak_fraction = 0.0;
    double weight = 0.0;
    double total_fraction = 0.0;
    int bins = 0;
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

struct IndexedDistance {
    double d2 = 0.0;
    double weight = 0.0;
    std::size_t index = 0;
};

std::string replace_suffix(const std::string& path, const std::string& suffix) {
    const std::size_t dot = path.find_last_of('.');
    if (dot == std::string::npos) return path + suffix;
    return path.substr(0, dot) + suffix;
}

std::string dirname_of(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return ".";
    return path.substr(0, slash);
}

bool require_branch(TTree* tree, const char* name) {
    if (tree->GetBranch(name)) return true;
    std::cerr << "ERROR: missing branch '" << name << "'" << std::endl;
    return false;
}

void draw_peak_lines(TH2D* h, int ix, int iy) {
    const double x = h->GetXaxis()->GetBinCenter(ix);
    const double y = h->GetYaxis()->GetBinCenter(iy);

    TLine* lx = new TLine(x, h->GetYaxis()->GetXmin(), x, h->GetYaxis()->GetXmax());
    lx->SetLineColor(kBlack);
    lx->SetLineStyle(2);
    lx->SetLineWidth(2);
    lx->Draw();

    TLine* ly = new TLine(h->GetXaxis()->GetXmin(), y, h->GetXaxis()->GetXmax(), y);
    ly->SetLineColor(kBlack);
    ly->SetLineStyle(2);
    ly->SetLineWidth(2);
    ly->Draw();
}

double weighted_quantile(std::vector<WeightedValue>& values, double quantile) {
    if (values.empty()) return 0.0;
    if (quantile <= 0.0) return values.front().value;
    if (quantile >= 1.0) quantile = 1.0;

    std::sort(values.begin(), values.end(),
              [](const WeightedValue& a, const WeightedValue& b) {
                  return a.value < b.value;
              });

    double total = 0.0;
    for (const WeightedValue& v : values) total += v.weight;
    if (total <= 0.0) return values.back().value;

    const double target = quantile * total;
    double running = 0.0;
    for (const WeightedValue& v : values) {
        running += v.weight;
        if (running >= target) return v.value;
    }
    return values.back().value;
}

double covariance_d2(const CovModel& model, double x, double y) {
    const double dx = x - model.mean_x;
    const double dy = y - model.mean_y;
    return (model.cov_yy * dx * dx - 2.0 * model.cov_xy * dx * dy +
            model.cov_xx * dy * dy) / model.det;
}

bool compute_cov_model(const std::vector<WeightedPoint>& points,
                       const std::vector<unsigned char>& mask,
                       CovModel& model,
                       double& weight_sum,
                       int& n_bins)
{
    weight_sum = 0.0;
    n_bins = 0;
    model = CovModel();

    for (const WeightedPoint& p : points) {
        if (!mask[p.global_bin]) continue;
        weight_sum += p.weight;
        model.mean_x += p.weight * p.x;
        model.mean_y += p.weight * p.y;
        ++n_bins;
    }
    if (weight_sum <= 0.0 || n_bins < 3) return false;

    model.mean_x /= weight_sum;
    model.mean_y /= weight_sum;

    for (const WeightedPoint& p : points) {
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

TPolyLine* make_ellipse_line(double mean_x, double mean_y,
                             double cov_xx, double cov_xy, double cov_yy,
                             double d2_cut,
                             Color_t color, Style_t style, Width_t width)
{
    const double trace = cov_xx + cov_yy;
    const double diff = cov_xx - cov_yy;
    const double root = std::sqrt(0.25 * diff * diff + cov_xy * cov_xy);
    const double lambda1 = std::max(0.0, 0.5 * trace + root);
    const double lambda2 = std::max(0.0, 0.5 * trace - root);
    const double angle = 0.5 * std::atan2(2.0 * cov_xy, diff);
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
        const double x = mean_x + u * ca - v * sa;
        const double y = mean_y + u * sa + v * ca;
        line->SetPoint(i, x, y);
    }
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    return line;
}

} // namespace

void test_2D_mass_cut(
    const char* input_path = "output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root",
    const char* tree_name = "physics",
    // double peak_fraction = -1.0,
    double peak_fraction = 0.10,

    int n_mpi0_bins = 160,
    double mpi0_min = 0.11,
    double mpi0_max = 0.15,
    int n_mmiss_bins = 200,
    double mmiss_min = 0.6,
    double mmiss_max = 1.5,
    double ellipse_core_quantile = 0.995,
    double ellipse_padding = 1.05,
    double auto_peak_min = 0.05,
    double auto_peak_max = 0.60,
    double auto_peak_step = 0.005,
    double auto_min_core_total_fraction = 0.05,
    int auto_min_core_bins = 30,
    double mcd_keep_total_fraction = 0.30,
    int mcd_iterations = 8,
    double mcd_ellipse_quantile = 0.995,
    double mcd_padding = 1.05,
    const char* output_root_path = "output/plots/x60_4b/production_wfpi0/test_2D_mass_cut.root")
    // const char* output_root_path = "output/plots/x60_4b/diagnostics_run4398_test_2D_mass_cut.root")
{
    const bool auto_peak_fraction = (peak_fraction <= 0.0);
    if (!auto_peak_fraction && peak_fraction > 1.0) {
        std::cerr << "ERROR: peak_fraction must be in (0, 1], or <=0 for auto. Got "
                  << peak_fraction << std::endl;
        return;
    }
    if (ellipse_core_quantile <= 0.0 || ellipse_core_quantile > 1.0) {
        std::cerr << "ERROR: ellipse_core_quantile must be in (0, 1]. Got "
                  << ellipse_core_quantile << std::endl;
        return;
    }
    if (ellipse_padding <= 0.0) {
        std::cerr << "ERROR: ellipse_padding must be positive. Got "
                  << ellipse_padding << std::endl;
        return;
    }
    if (auto_peak_fraction &&
        (auto_peak_min <= 0.0 || auto_peak_max > 1.0 ||
         auto_peak_min >= auto_peak_max || auto_peak_step <= 0.0)) {
        std::cerr << "ERROR: invalid auto peak scan settings." << std::endl;
        return;
    }
    if (mcd_keep_total_fraction <= 0.0 || mcd_keep_total_fraction > 1.0 ||
        mcd_iterations < 1 ||
        mcd_ellipse_quantile <= 0.0 || mcd_ellipse_quantile > 1.0 ||
        mcd_padding <= 0.0) {
        std::cerr << "ERROR: invalid MCD settings." << std::endl;
        return;
    }

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    TFile* input_file = TFile::Open(input_path, "READ");
    if (!input_file || input_file->IsZombie()) {
        std::cerr << "ERROR: cannot open input file: " << input_path << std::endl;
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(input_file->Get(tree_name));
    if (!tree) {
        std::cerr << "ERROR: cannot find tree '" << tree_name
                  << "' in " << input_path << std::endl;
        input_file->Close();
        return;
    }

    if (!require_branch(tree, "mpi0_all") ||
        !require_branch(tree, "mmiss_all") ||
        !require_branch(tree, "pi0_weight")) {
        input_file->Close();
        return;
    }

    Double_t mpi0_all = 0.0;
    Double_t mmiss_all = 0.0;
    Double_t pi0_weight = 0.0;

    tree->SetBranchAddress("mpi0_all", &mpi0_all);
    tree->SetBranchAddress("mmiss_all", &mmiss_all);
    tree->SetBranchAddress("pi0_weight", &pi0_weight);

    TH2D* h2 = new TH2D(
        "h_mmiss_vs_mpi0_weighted",
        "Weighted mmiss_all:mpi0_all;mpi0_all [GeV];mmiss_all [GeV]",
        n_mpi0_bins, mpi0_min, mpi0_max,
        n_mmiss_bins, mmiss_min, mmiss_max);
    h2->SetDirectory(nullptr);

    Long64_t n_read = 0;
    Long64_t n_filled = 0;
    double raw_weight_sum = 0.0;

    const Long64_t n_entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < n_entries; ++entry) {
        tree->GetEntry(entry);
        ++n_read;

        if (!std::isfinite(mpi0_all) ||
            !std::isfinite(mmiss_all) ||
            !std::isfinite(pi0_weight)) {
            continue;
        }
        if (pi0_weight <= 0.0) continue;
        if (mpi0_all < mpi0_min || mpi0_all >= mpi0_max) continue;
        if (mmiss_all < mmiss_min || mmiss_all >= mmiss_max) continue;

        h2->Fill(mpi0_all, mmiss_all, pi0_weight);
        raw_weight_sum += pi0_weight;
        ++n_filled;
    }

    const int nx = h2->GetNbinsX();
    const int ny = h2->GetNbinsY();
    const int n_cells = h2->GetNcells();
    const double total_weight = h2->Integral(1, nx, 1, ny);
    if (total_weight <= 0.0) {
        std::cerr << "ERROR: weighted histogram is empty after cuts." << std::endl;
        input_file->Close();
        return;
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
    int peak_ix = 0;
    int peak_iy = 0;
    int peak_iz = 0;
    h2->GetBinXYZ(peak_global, peak_ix, peak_iy, peak_iz);
    const double peak_weight = h2->GetBinContent(peak_global);
    const double peak_tolerance = std::max(1e-12, std::abs(peak_weight) * 1e-9);

    auto build_core = [&](double candidate_peak_fraction,
                          std::vector<unsigned char>& core_mask) {
        core_mask.assign(n_cells + 1, 0);
        std::queue<std::pair<int, int> > core_frontier;
        const double candidate_threshold_weight = candidate_peak_fraction * peak_weight;

        auto push_if_selected = [&](int ix, int iy) {
            if (ix < 1 || ix > nx || iy < 1 || iy > ny) return;
            const int global = h2->GetBin(ix, iy);
            if (core_mask[global]) return;
            const double w = h2->GetBinContent(global);
            if (w + 1e-12 < candidate_threshold_weight) return;
            core_mask[global] = 1;
            core_frontier.push(std::make_pair(ix, iy));
        };

        for (int ix = 1; ix <= nx; ++ix) {
            for (int iy = 1; iy <= ny; ++iy) {
                const double w = h2->GetBinContent(ix, iy);
                if (std::abs(w - peak_weight) <= peak_tolerance) {
                    push_if_selected(ix, iy);
                }
            }
        }

        while (!core_frontier.empty()) {
            const std::pair<int, int> cand = core_frontier.front();
            core_frontier.pop();
            const int ix = cand.first;
            const int iy = cand.second;
            push_if_selected(ix - 1, iy);
            push_if_selected(ix + 1, iy);
            push_if_selected(ix, iy - 1);
            push_if_selected(ix, iy + 1);
            push_if_selected(ix - 1, iy - 1);
            push_if_selected(ix - 1, iy + 1);
            push_if_selected(ix + 1, iy - 1);
            push_if_selected(ix + 1, iy + 1);
        }

        CoreStats stats;
        stats.peak_fraction = candidate_peak_fraction;
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
    double auto_jump_ratio = 1.0;
    if (auto_peak_fraction) {
        std::vector<unsigned char> scan_mask;
        CoreStats previous;
        bool have_previous = false;
        double best_jump = 0.0;
        double chosen_peak_fraction = auto_peak_max;

        const int n_scan = static_cast<int>(std::floor((auto_peak_max - auto_peak_min) / auto_peak_step + 0.5)) + 1;
        for (int i = 0; i < n_scan; ++i) {
            double candidate = auto_peak_max - i * auto_peak_step;
            if (candidate < auto_peak_min) candidate = auto_peak_min;
            CoreStats current = build_core(candidate, scan_mask);
            scan_stats.push_back(current);

            if (have_previous && previous.bins > 0 && previous.weight > 0.0) {
                const double weight_ratio = current.weight / previous.weight;
                const double bin_ratio = static_cast<double>(current.bins) / static_cast<double>(previous.bins);
                const double jump = std::max(weight_ratio, bin_ratio);
                const bool previous_core_is_real =
                    previous.bins >= auto_min_core_bins &&
                    previous.total_fraction >= auto_min_core_total_fraction;
                if (previous_core_is_real && jump > best_jump) {
                    best_jump = jump;
                    chosen_peak_fraction = previous.peak_fraction;
                }
            }
            previous = current;
            have_previous = true;
        }

        auto_jump_ratio = best_jump;
        peak_fraction = chosen_peak_fraction;
    }

    const double threshold_weight = peak_fraction * peak_weight;
    std::vector<unsigned char> selected(n_cells + 1, 0);
    build_core(peak_fraction, selected);

    double selected_weight = 0.0;
    int selected_bins = 0;
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int global = h2->GetBin(ix, iy);
            if (!selected[global]) continue;
            selected_weight += h2->GetBinContent(ix, iy);
            ++selected_bins;
        }
    }
    const double selected_total_fraction = selected_weight / total_weight;
    if (selected_bins < 3 || selected_weight <= 0.0) {
        std::cerr << "ERROR: too few selected core bins to fit ellipse." << std::endl;
        input_file->Close();
        return;
    }

    double mean_x = 0.0;
    double mean_y = 0.0;
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int global = h2->GetBin(ix, iy);
            if (!selected[global]) continue;
            const double w = h2->GetBinContent(ix, iy);
            mean_x += w * h2->GetXaxis()->GetBinCenter(ix);
            mean_y += w * h2->GetYaxis()->GetBinCenter(iy);
        }
    }
    mean_x /= selected_weight;
    mean_y /= selected_weight;

    double cov_xx = 0.0;
    double cov_xy = 0.0;
    double cov_yy = 0.0;
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int global = h2->GetBin(ix, iy);
            if (!selected[global]) continue;
            const double w = h2->GetBinContent(ix, iy);
            const double dx = h2->GetXaxis()->GetBinCenter(ix) - mean_x;
            const double dy = h2->GetYaxis()->GetBinCenter(iy) - mean_y;
            cov_xx += w * dx * dx;
            cov_xy += w * dx * dy;
            cov_yy += w * dy * dy;
        }
    }
    cov_xx /= selected_weight;
    cov_xy /= selected_weight;
    cov_yy /= selected_weight;

    double cov_det = cov_xx * cov_yy - cov_xy * cov_xy;
    if (cov_det <= 0.0) {
        std::cerr << "ERROR: ellipse covariance is singular." << std::endl;
        input_file->Close();
        return;
    }

    auto mahalanobis_d2 = [&](double x, double y) {
        const double dx = x - mean_x;
        const double dy = y - mean_y;
        return (cov_yy * dx * dx - 2.0 * cov_xy * dx * dy + cov_xx * dy * dy) / cov_det;
    };

    std::vector<WeightedValue> core_d2_values;
    core_d2_values.reserve(selected_bins);
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int global = h2->GetBin(ix, iy);
            if (!selected[global]) continue;
            const double w = h2->GetBinContent(ix, iy);
            const double x = h2->GetXaxis()->GetBinCenter(ix);
            const double y = h2->GetYaxis()->GetBinCenter(iy);
            core_d2_values.push_back({mahalanobis_d2(x, y), w});
        }
    }

    const double core_d2_quantile = weighted_quantile(core_d2_values, ellipse_core_quantile);
    const double ellipse_d2_cut = core_d2_quantile * ellipse_padding * ellipse_padding;

    std::vector<unsigned char> ellipse_selected(n_cells + 1, 0);
    double ellipse_weight = 0.0;
    int ellipse_bins = 0;
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const double x = h2->GetXaxis()->GetBinCenter(ix);
            const double y = h2->GetYaxis()->GetBinCenter(iy);
            if (mahalanobis_d2(x, y) > ellipse_d2_cut) continue;
            const int global = h2->GetBin(ix, iy);
            ellipse_selected[global] = 1;
            const double w = h2->GetBinContent(ix, iy);
            ellipse_weight += w;
            if (w > 0.0) ++ellipse_bins;
        }
    }
    const double ellipse_total_fraction = ellipse_weight / total_weight;

    CovModel mcd_model;
    double mcd_subset_weight = 0.0;
    int mcd_subset_bins = 0;
    if (!compute_cov_model(points, selected, mcd_model, mcd_subset_weight, mcd_subset_bins)) {
        std::cerr << "ERROR: cannot initialize MCD covariance." << std::endl;
        input_file->Close();
        return;
    }

    CovModel best_mcd_model = mcd_model;
    std::vector<unsigned char> best_mcd_subset(n_cells + 1, 0);
    std::vector<unsigned char> current_mcd_subset = selected;
    best_mcd_subset = selected;
    double best_mcd_det = std::numeric_limits<double>::infinity();
    double best_mcd_subset_weight = mcd_subset_weight;
    int best_mcd_subset_bins = mcd_subset_bins;
    const double mcd_target_weight = mcd_keep_total_fraction * total_weight;

    for (int iter = 0; iter < mcd_iterations; ++iter) {
        std::vector<IndexedDistance> distances;
        distances.reserve(points.size());
        for (std::size_t ip = 0; ip < points.size(); ++ip) {
            distances.push_back({covariance_d2(mcd_model, points[ip].x, points[ip].y),
                                 points[ip].weight,
                                 ip});
        }
        std::sort(distances.begin(), distances.end(),
                  [](const IndexedDistance& a, const IndexedDistance& b) {
                      return a.d2 < b.d2;
                  });

        current_mcd_subset.assign(n_cells + 1, 0);
        double kept_weight = 0.0;
        for (const IndexedDistance& d : distances) {
            const WeightedPoint& p = points[d.index];
            current_mcd_subset[p.global_bin] = 1;
            kept_weight += p.weight;
            if (kept_weight >= mcd_target_weight) break;
        }

        CovModel candidate_model;
        double candidate_weight = 0.0;
        int candidate_bins = 0;
        if (!compute_cov_model(points, current_mcd_subset, candidate_model,
                               candidate_weight, candidate_bins)) {
            break;
        }
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
        best_mcd_det = mcd_model.det;
        best_mcd_model = mcd_model;
        best_mcd_subset = selected;
        best_mcd_subset_weight = mcd_subset_weight;
        best_mcd_subset_bins = mcd_subset_bins;
    }

    std::vector<WeightedValue> mcd_d2_values;
    mcd_d2_values.reserve(best_mcd_subset_bins);
    for (const WeightedPoint& p : points) {
        if (!best_mcd_subset[p.global_bin]) continue;
        mcd_d2_values.push_back({covariance_d2(best_mcd_model, p.x, p.y), p.weight});
    }

    const double mcd_core_d2_quantile = weighted_quantile(mcd_d2_values, mcd_ellipse_quantile);
    const double mcd_d2_cut = mcd_core_d2_quantile * mcd_padding * mcd_padding;
    std::vector<unsigned char> mcd_selected(n_cells + 1, 0);
    double mcd_weight = 0.0;
    int mcd_bins = 0;
    for (const WeightedPoint& p : points) {
        if (covariance_d2(best_mcd_model, p.x, p.y) > mcd_d2_cut) continue;
        mcd_selected[p.global_bin] = 1;
        mcd_weight += p.weight;
        ++mcd_bins;
    }
    const double mcd_total_fraction = mcd_weight / total_weight;

    TH2D* h2_selected = dynamic_cast<TH2D*>(h2->Clone("h_mmiss_vs_mpi0_selected"));
    h2_selected->SetTitle("Selected weighted mmiss_all:mpi0_all;mpi0_all [GeV];mmiss_all [GeV]");
    h2_selected->SetDirectory(nullptr);

    TH2D* h2_ellipse = dynamic_cast<TH2D*>(h2->Clone("h_mmiss_vs_mpi0_ellipse"));
    h2_ellipse->SetTitle("Ellipse selected weighted mmiss_all:mpi0_all;mpi0_all [GeV];mmiss_all [GeV]");
    h2_ellipse->SetDirectory(nullptr);

    TH2D* h2_mcd = dynamic_cast<TH2D*>(h2->Clone("h_mmiss_vs_mpi0_mcd"));
    h2_mcd->SetTitle("MCD selected weighted mmiss_all:mpi0_all;mpi0_all [GeV];mmiss_all [GeV]");
    h2_mcd->SetDirectory(nullptr);

    TH2D* h_mask = dynamic_cast<TH2D*>(h2->Clone("h_2d_mass_cut_mask"));
    h_mask->SetTitle("2D mass cut mask;mpi0_all [GeV];mmiss_all [GeV]");
    h_mask->Reset("ICES");
    h_mask->SetDirectory(nullptr);

    TH2D* h_ellipse_mask = dynamic_cast<TH2D*>(h2->Clone("h_ellipse_mass_cut_mask"));
    h_ellipse_mask->SetTitle("Ellipse mass cut mask;mpi0_all [GeV];mmiss_all [GeV]");
    h_ellipse_mask->Reset("ICES");
    h_ellipse_mask->SetDirectory(nullptr);

    TH2D* h_mcd_mask = dynamic_cast<TH2D*>(h2->Clone("h_mcd_mass_cut_mask"));
    h_mcd_mask->SetTitle("MCD mass cut mask;mpi0_all [GeV];mmiss_all [GeV]");
    h_mcd_mask->Reset("ICES");
    h_mcd_mask->SetDirectory(nullptr);

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int global = h2->GetBin(ix, iy);
            if (selected[global]) {
                h_mask->SetBinContent(ix, iy, 1.0);
            } else {
                h2_selected->SetBinContent(ix, iy, 0.0);
            }
            if (ellipse_selected[global]) {
                h_ellipse_mask->SetBinContent(ix, iy, 1.0);
            } else {
                h2_ellipse->SetBinContent(ix, iy, 0.0);
            }
            if (mcd_selected[global]) {
                h_mcd_mask->SetBinContent(ix, iy, 1.0);
            } else {
                h2_mcd->SetBinContent(ix, iy, 0.0);
            }
        }
    }

    TH1D* h_mpi0_all = h2->ProjectionX("h_mpi0_all_weighted", 1, ny);
    TH1D* h_mpi0_selected = h2_selected->ProjectionX("h_mpi0_selected_weighted", 1, ny);
    TH1D* h_mpi0_ellipse = h2_ellipse->ProjectionX("h_mpi0_ellipse_weighted", 1, ny);
    TH1D* h_mpi0_mcd = h2_mcd->ProjectionX("h_mpi0_mcd_weighted", 1, ny);
    TH1D* h_mmiss_all = h2->ProjectionY("h_mmiss_all_weighted", 1, nx);
    TH1D* h_mmiss_selected = h2_selected->ProjectionY("h_mmiss_selected_weighted", 1, nx);
    TH1D* h_mmiss_ellipse = h2_ellipse->ProjectionY("h_mmiss_ellipse_weighted", 1, nx);
    TH1D* h_mmiss_mcd = h2_mcd->ProjectionY("h_mmiss_mcd_weighted", 1, nx);

    h_mpi0_all->SetTitle("mpi0_all projection;mpi0_all [GeV];Weighted counts");
    h_mpi0_selected->SetTitle("mpi0_all projection;mpi0_all [GeV];Weighted counts");
    h_mpi0_ellipse->SetTitle("mpi0_all projection;mpi0_all [GeV];Weighted counts");
    h_mpi0_mcd->SetTitle("mpi0_all projection;mpi0_all [GeV];Weighted counts");
    h_mmiss_all->SetTitle("mmiss_all projection;mmiss_all [GeV];Weighted counts");
    h_mmiss_selected->SetTitle("mmiss_all projection;mmiss_all [GeV];Weighted counts");
    h_mmiss_ellipse->SetTitle("mmiss_all projection;mmiss_all [GeV];Weighted counts");
    h_mmiss_mcd->SetTitle("mmiss_all projection;mmiss_all [GeV];Weighted counts");

    h_mpi0_all->SetLineColor(kBlack);
    h_mpi0_all->SetLineWidth(2);
    h_mpi0_selected->SetLineColor(kRed + 1);
    h_mpi0_selected->SetLineWidth(2);
    h_mpi0_ellipse->SetLineColor(kMagenta + 2);
    h_mpi0_ellipse->SetLineWidth(2);
    h_mpi0_mcd->SetLineColor(kGreen + 2);
    h_mpi0_mcd->SetLineWidth(2);
    h_mmiss_all->SetLineColor(kBlack);
    h_mmiss_all->SetLineWidth(2);
    h_mmiss_selected->SetLineColor(kRed + 1);
    h_mmiss_selected->SetLineWidth(2);
    h_mmiss_ellipse->SetLineColor(kMagenta + 2);
    h_mmiss_ellipse->SetLineWidth(2);
    h_mmiss_mcd->SetLineColor(kGreen + 2);
    h_mmiss_mcd->SetLineWidth(2);

    const std::string output_root = output_root_path;
    const std::string output_dir = dirname_of(output_root);
    gSystem->mkdir(output_dir.c_str(), kTRUE);

    const std::string output_pdf = replace_suffix(output_root, ".pdf");
    const std::string output_png = replace_suffix(output_root, ".png");
    const std::string output_csv = replace_suffix(output_root, "_mask_bins.csv");
    const std::string output_ellipse_csv = replace_suffix(output_root, "_ellipse_bins.csv");
    const std::string output_mcd_csv = replace_suffix(output_root, "_mcd_bins.csv");
    const std::string output_params_csv = replace_suffix(output_root, "_ellipse_params.csv");
    const std::string output_scan_csv = replace_suffix(output_root, "_peak_scan.csv");

    std::ofstream csv(output_csv.c_str());
    csv << "ix,iy,mpi0_low,mpi0_high,mmiss_low,mmiss_high,bin_weight,bin_over_peak\n";
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int global = h2->GetBin(ix, iy);
            if (!selected[global]) continue;
            csv << ix << "," << iy << ","
                << h2->GetXaxis()->GetBinLowEdge(ix) << ","
                << h2->GetXaxis()->GetBinUpEdge(ix) << ","
                << h2->GetYaxis()->GetBinLowEdge(iy) << ","
                << h2->GetYaxis()->GetBinUpEdge(iy) << ","
                << h2->GetBinContent(ix, iy) << ","
                << h2->GetBinContent(ix, iy) / peak_weight << "\n";
        }
    }
    csv.close();

    std::ofstream ellipse_csv(output_ellipse_csv.c_str());
    ellipse_csv << "ix,iy,mpi0_low,mpi0_high,mmiss_low,mmiss_high,bin_weight,d2\n";
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int global = h2->GetBin(ix, iy);
            if (!ellipse_selected[global]) continue;
            const double x = h2->GetXaxis()->GetBinCenter(ix);
            const double y = h2->GetYaxis()->GetBinCenter(iy);
            ellipse_csv << ix << "," << iy << ","
                        << h2->GetXaxis()->GetBinLowEdge(ix) << ","
                        << h2->GetXaxis()->GetBinUpEdge(ix) << ","
                        << h2->GetYaxis()->GetBinLowEdge(iy) << ","
                        << h2->GetYaxis()->GetBinUpEdge(iy) << ","
                        << h2->GetBinContent(ix, iy) << ","
                        << mahalanobis_d2(x, y) << "\n";
        }
    }
    ellipse_csv.close();

    std::ofstream mcd_csv(output_mcd_csv.c_str());
    mcd_csv << "ix,iy,mpi0_low,mpi0_high,mmiss_low,mmiss_high,bin_weight,mcd_d2\n";
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            const int global = h2->GetBin(ix, iy);
            if (!mcd_selected[global]) continue;
            const double x = h2->GetXaxis()->GetBinCenter(ix);
            const double y = h2->GetYaxis()->GetBinCenter(iy);
            mcd_csv << ix << "," << iy << ","
                    << h2->GetXaxis()->GetBinLowEdge(ix) << ","
                    << h2->GetXaxis()->GetBinUpEdge(ix) << ","
                    << h2->GetYaxis()->GetBinLowEdge(iy) << ","
                    << h2->GetYaxis()->GetBinUpEdge(iy) << ","
                    << h2->GetBinContent(ix, iy) << ","
                    << covariance_d2(best_mcd_model, x, y) << "\n";
        }
    }
    mcd_csv.close();

    const double trace = cov_xx + cov_yy;
    const double diff = cov_xx - cov_yy;
    const double root = std::sqrt(0.25 * diff * diff + cov_xy * cov_xy);
    const double lambda1 = std::max(0.0, 0.5 * trace + root);
    const double lambda2 = std::max(0.0, 0.5 * trace - root);
    const double ellipse_angle_rad = 0.5 * std::atan2(2.0 * cov_xy, diff);
    const double ellipse_axis1 = std::sqrt(std::max(0.0, ellipse_d2_cut * lambda1));
    const double ellipse_axis2 = std::sqrt(std::max(0.0, ellipse_d2_cut * lambda2));

    std::ofstream params_csv(output_params_csv.c_str());
    params_csv << "parameter,value\n";
    params_csv << "mean_mpi0," << mean_x << "\n";
    params_csv << "mean_mmiss," << mean_y << "\n";
    params_csv << "cov_mpi0_mpi0," << cov_xx << "\n";
    params_csv << "cov_mpi0_mmiss," << cov_xy << "\n";
    params_csv << "cov_mmiss_mmiss," << cov_yy << "\n";
    params_csv << "ellipse_d2_cut," << ellipse_d2_cut << "\n";
    params_csv << "ellipse_axis1," << ellipse_axis1 << "\n";
    params_csv << "ellipse_axis2," << ellipse_axis2 << "\n";
    params_csv << "ellipse_angle_rad," << ellipse_angle_rad << "\n";
    params_csv << "ellipse_angle_deg," << ellipse_angle_rad * 180.0 / TMath::Pi() << "\n";
    params_csv << "peak_fraction," << peak_fraction << "\n";
    params_csv << "auto_peak_fraction," << (auto_peak_fraction ? 1 : 0) << "\n";
    params_csv << "auto_peak_min," << auto_peak_min << "\n";
    params_csv << "auto_peak_max," << auto_peak_max << "\n";
    params_csv << "auto_peak_step," << auto_peak_step << "\n";
    params_csv << "auto_min_core_total_fraction," << auto_min_core_total_fraction << "\n";
    params_csv << "auto_min_core_bins," << auto_min_core_bins << "\n";
    params_csv << "auto_jump_ratio," << auto_jump_ratio << "\n";
    params_csv << "ellipse_core_quantile," << ellipse_core_quantile << "\n";
    params_csv << "ellipse_padding," << ellipse_padding << "\n";
    params_csv << "mcd_mean_mpi0," << best_mcd_model.mean_x << "\n";
    params_csv << "mcd_mean_mmiss," << best_mcd_model.mean_y << "\n";
    params_csv << "mcd_cov_mpi0_mpi0," << best_mcd_model.cov_xx << "\n";
    params_csv << "mcd_cov_mpi0_mmiss," << best_mcd_model.cov_xy << "\n";
    params_csv << "mcd_cov_mmiss_mmiss," << best_mcd_model.cov_yy << "\n";
    params_csv << "mcd_det," << best_mcd_model.det << "\n";
    params_csv << "mcd_d2_cut," << mcd_d2_cut << "\n";
    params_csv << "mcd_keep_total_fraction," << mcd_keep_total_fraction << "\n";
    params_csv << "mcd_target_weight," << mcd_target_weight << "\n";
    params_csv << "mcd_subset_weight," << best_mcd_subset_weight << "\n";
    params_csv << "mcd_subset_bins," << best_mcd_subset_bins << "\n";
    params_csv << "mcd_iterations," << mcd_iterations << "\n";
    params_csv << "mcd_ellipse_quantile," << mcd_ellipse_quantile << "\n";
    params_csv << "mcd_padding," << mcd_padding << "\n";
    params_csv.close();

    std::ofstream scan_csv(output_scan_csv.c_str());
    scan_csv << "peak_fraction,selected_bins,selected_weight,selected_total_fraction,weight_ratio_to_previous,bin_ratio_to_previous\n";
    for (std::size_t i = 0; i < scan_stats.size(); ++i) {
        double weight_ratio = 0.0;
        double bin_ratio = 0.0;
        if (i > 0 && scan_stats[i - 1].weight > 0.0 && scan_stats[i - 1].bins > 0) {
            weight_ratio = scan_stats[i].weight / scan_stats[i - 1].weight;
            bin_ratio = static_cast<double>(scan_stats[i].bins) / static_cast<double>(scan_stats[i - 1].bins);
        }
        scan_csv << scan_stats[i].peak_fraction << ","
                 << scan_stats[i].bins << ","
                 << scan_stats[i].weight << ","
                 << scan_stats[i].total_fraction << ","
                 << weight_ratio << ","
                 << bin_ratio << "\n";
    }
    scan_csv.close();

    TCanvas* c = new TCanvas("c_test_2d_mass_cut", "2D mass cut test", 1600, 1200);
    c->Divide(2, 2);

    double contour_level[1] = {0.5};
    h_mask->SetContour(1, contour_level);
    h_mask->SetLineColor(kRed + 1);
    h_mask->SetLineWidth(3);
    h_ellipse_mask->SetContour(1, contour_level);
    h_ellipse_mask->SetLineColor(kMagenta + 2);
    h_ellipse_mask->SetLineWidth(3);
    h_mcd_mask->SetContour(1, contour_level);
    h_mcd_mask->SetLineColor(kGreen + 2);
    h_mcd_mask->SetLineWidth(3);
    TPolyLine* ellipse_line = make_ellipse_line(mean_x, mean_y, cov_xx, cov_xy, cov_yy,
                                                ellipse_d2_cut, kMagenta + 2, 1, 3);
    TPolyLine* ellipse_line2 = make_ellipse_line(mean_x, mean_y, cov_xx, cov_xy, cov_yy,
                                                 ellipse_d2_cut, kMagenta + 2, 1, 3);
    TPolyLine* mcd_line = make_ellipse_line(best_mcd_model.mean_x, best_mcd_model.mean_y,
                                            best_mcd_model.cov_xx, best_mcd_model.cov_xy,
                                            best_mcd_model.cov_yy, mcd_d2_cut,
                                            kGreen + 2, 1, 3);
    TPolyLine* mcd_line2 = make_ellipse_line(best_mcd_model.mean_x, best_mcd_model.mean_y,
                                             best_mcd_model.cov_xx, best_mcd_model.cov_xy,
                                             best_mcd_model.cov_yy, mcd_d2_cut,
                                             kGreen + 2, 1, 3);

    c->cd(1);
    gPad->SetRightMargin(0.14);
    h2->Draw("COLZ");
    h_mask->Draw("CONT3 SAME");
    ellipse_line->Draw("SAME");
    mcd_line->Draw("SAME");
    draw_peak_lines(h2, peak_ix, peak_iy);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.032);
    latex.DrawLatex(0.13, 0.93, Form("Threshold %.1f%% of peak bin",
                                      100.0 * peak_fraction));
    latex.DrawLatex(0.13, 0.88, Form("Peak: mpi0=%.5f, mmiss=%.5f",
                                      h2->GetXaxis()->GetBinCenter(peak_ix),
                                      h2->GetYaxis()->GetBinCenter(peak_iy)));
    latex.DrawLatex(0.13, 0.83, Form("Selected total weight: %.2f%%",
                                      100.0 * selected_total_fraction));
    latex.DrawLatex(0.13, 0.78, Form("Ellipse total weight: %.2f%%",
                                      100.0 * ellipse_total_fraction));
    latex.DrawLatex(0.13, 0.73, Form("MCD total weight: %.2f%%",
                                      100.0 * mcd_total_fraction));

    c->cd(2);
    gPad->SetRightMargin(0.14);
    h2_mcd->Draw("COLZ");
    h_mask->Draw("CONT3 SAME");
    ellipse_line2->Draw("SAME");
    mcd_line2->Draw("SAME");

    c->cd(3);
    h_mpi0_all->Draw("HIST");
    h_mpi0_selected->Draw("HIST SAME");
    h_mpi0_ellipse->Draw("HIST SAME");
    h_mpi0_mcd->Draw("HIST SAME");
    TLegend* leg_mpi0 = new TLegend(0.56, 0.74, 0.88, 0.88);
    leg_mpi0->AddEntry(h_mpi0_all, "all weighted", "l");
    leg_mpi0->AddEntry(h_mpi0_selected, "core threshold", "l");
    leg_mpi0->AddEntry(h_mpi0_ellipse, "ellipse cut", "l");
    leg_mpi0->AddEntry(h_mpi0_mcd, "MCD cut", "l");
    leg_mpi0->Draw();

    c->cd(4);
    h_mmiss_all->Draw("HIST");
    h_mmiss_selected->Draw("HIST SAME");
    h_mmiss_ellipse->Draw("HIST SAME");
    h_mmiss_mcd->Draw("HIST SAME");
    TLegend* leg_mmiss = new TLegend(0.56, 0.74, 0.88, 0.88);
    leg_mmiss->AddEntry(h_mmiss_all, "all weighted", "l");
    leg_mmiss->AddEntry(h_mmiss_selected, "core threshold", "l");
    leg_mmiss->AddEntry(h_mmiss_ellipse, "ellipse cut", "l");
    leg_mmiss->AddEntry(h_mmiss_mcd, "MCD cut", "l");
    leg_mmiss->Draw();

    c->SaveAs(output_pdf.c_str());
    c->SaveAs(output_png.c_str());

    TFile* output_file = TFile::Open(output_root.c_str(), "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        std::cerr << "ERROR: cannot create output file: " << output_root << std::endl;
        input_file->Close();
        return;
    }

    h2->Write();
    h2_selected->Write();
    h2_ellipse->Write();
    h2_mcd->Write();
    h_mask->Write();
    h_ellipse_mask->Write();
    h_mcd_mask->Write();
    h_mpi0_all->Write();
    h_mpi0_selected->Write();
    h_mpi0_ellipse->Write();
    h_mpi0_mcd->Write();
    h_mmiss_all->Write();
    h_mmiss_selected->Write();
    h_mmiss_ellipse->Write();
    h_mmiss_mcd->Write();
    c->Write();

    output_file->Close();
    input_file->Close();

    std::cout << "[2D_MASS_CUT_TEST]" << std::endl;
    std::cout << "  input: " << input_path << std::endl;
    std::cout << "  tree: " << tree_name << std::endl;
    std::cout << "  entries read: " << n_read << std::endl;
    std::cout << "  entries filled: " << n_filled << std::endl;
    std::cout << "  raw in-range pi0_weight sum: " << raw_weight_sum << std::endl;
    std::cout << "  histogram total weight: " << total_weight << std::endl;
    std::cout << "  peak bin weight: " << peak_weight << std::endl;
    std::cout << "  auto peak fraction: " << (auto_peak_fraction ? "true" : "false") << std::endl;
    std::cout << "  auto jump ratio: " << auto_jump_ratio << std::endl;
    std::cout << "  peak fraction threshold: " << peak_fraction << std::endl;
    std::cout << "  bin weight threshold: " << threshold_weight << std::endl;
    std::cout << "  peak bin ix,iy: " << peak_ix << "," << peak_iy << std::endl;
    std::cout << "  peak mpi0,mmiss: "
              << h2->GetXaxis()->GetBinCenter(peak_ix) << ","
              << h2->GetYaxis()->GetBinCenter(peak_iy) << std::endl;
    std::cout << "  selected bins: " << selected_bins << std::endl;
    std::cout << "  selected weight: " << selected_weight << std::endl;
    std::cout << "  selected total fraction: " << selected_total_fraction << std::endl;
    std::cout << "  ellipse center mpi0,mmiss: " << mean_x << "," << mean_y << std::endl;
    std::cout << "  ellipse axes: " << ellipse_axis1 << "," << ellipse_axis2 << std::endl;
    std::cout << "  ellipse angle deg: " << ellipse_angle_rad * 180.0 / TMath::Pi() << std::endl;
    std::cout << "  ellipse d2 cut: " << ellipse_d2_cut << std::endl;
    std::cout << "  ellipse bins with weight: " << ellipse_bins << std::endl;
    std::cout << "  ellipse weight: " << ellipse_weight << std::endl;
    std::cout << "  ellipse total fraction: " << ellipse_total_fraction << std::endl;
    std::cout << "  MCD center mpi0,mmiss: "
              << best_mcd_model.mean_x << "," << best_mcd_model.mean_y << std::endl;
    std::cout << "  MCD det: " << best_mcd_model.det << std::endl;
    std::cout << "  MCD d2 cut: " << mcd_d2_cut << std::endl;
    std::cout << "  MCD subset bins: " << best_mcd_subset_bins << std::endl;
    std::cout << "  MCD subset weight: " << best_mcd_subset_weight << std::endl;
    std::cout << "  MCD bins with weight: " << mcd_bins << std::endl;
    std::cout << "  MCD weight: " << mcd_weight << std::endl;
    std::cout << "  MCD total fraction: " << mcd_total_fraction << std::endl;
    std::cout << "  output root: " << output_root << std::endl;
    std::cout << "  output pdf: " << output_pdf << std::endl;
    std::cout << "  output png: " << output_png << std::endl;
    std::cout << "  output csv: " << output_csv << std::endl;
    std::cout << "  output ellipse csv: " << output_ellipse_csv << std::endl;
    std::cout << "  output MCD csv: " << output_mcd_csv << std::endl;
    std::cout << "  output params csv: " << output_params_csv << std::endl;
    std::cout << "  output peak scan csv: " << output_scan_csv << std::endl;
}
