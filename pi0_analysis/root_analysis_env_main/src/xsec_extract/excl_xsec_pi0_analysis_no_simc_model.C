// g++ -O2 -std=c++17 -o excl_xsec_no_simc_model scripts/excl_xsec_pi0_analysis_no_simc_model.C `root-config --cflags --libs`

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGaxis.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TDecompLU.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TMath.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

struct AnalysisConfig {
    std::string simc_file = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/simc_pi0_analysis_output_smeared.root";
    std::string data_file = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root";
    std::string simc_tree = "simulation";
    std::string data_tree = "physics";

    std::string out_root = "excl_xsec_pi0_analysis_no_simc_model_output.root";
    std::string out_csv  = "excl_xsec_pi0_analysis_no_simc_model_summary.csv";
    
    std::string out_slice_csv = "excl_xsec_pi0_analysis_no_simc_model_slice_summary.csv";
    std::string out_dir  = "output_pi0_xsec_no_simc_model";
    std::string out_all_plots_pdf = "all_generated_plots_no_simc_model.pdf";

    double ebeam = 10.538; // GeV, set from Hall C beam energy used in the run period
    double mp    = 0.9382720813;
    double mpi0  = 0.1349768;

    int n_phi    = 12;
    int n_tprime = 5;
    int n_q2     = 2; // quantile bins in Q2
    int n_xb     = 2; // xB bins are built per-Q2 bin for balanced 2D occupancy

    double phi_min    = 0.0;
    double phi_max    = 2.0 * TMath::Pi();
    double tprime_min = -1.0;
    double tprime_max = 0.0;
    double q2_min     = 5.1;
    double q2_max     = 6.4;
    double xb_min     = 0.45;
    double xb_max     = 0.65;

    double tgt_contam     = 0.584;
    double tgt_contam_err  = 0.014;

    bool verbose = true;
    bool diagnostics = true;
    bool write_png = true;
    bool write_pdf = true;

    // Event-level SIMC CM model branch used for de-modeling weight:
    //   w_base = full_weight / sigcm
    // SIMC has siglab = sigcm * davejac * gtpr * fac. Dividing by sigcm
    // removes only the model cross section while preserving the pion
    // lab/CM Jacobian, SIMC flux convention, and Fermi-motion correction.
    // Candidate branch names are searched in this order.
    std::vector<std::string> model_xsec_candidates = {
        "sigcm"
    };
};

struct PhiBin {
    double data = 0.0, data_sumw2 = 0.0;
    double sim  = 0.0, sim_sumw2  = 0.0;

    // De-modeled SIMC basis sums, built with weight = full_weight/sigcm.
    double sim_base = 0.0;
    double sim_base_cos1 = 0.0;
    double sim_base_cos2 = 0.0;
    double sim_base_plus = 0.0;
    double sim_base_minus = 0.0;
    double sim_base_cos1_plus = 0.0;
    double sim_base_cos1_minus = 0.0;
    double sim_base_cos2_plus = 0.0;
    double sim_base_cos2_minus = 0.0;
    double sim_base_sin1_plus = 0.0;
    double sim_base_sin1_minus = 0.0;

    double data_plus = 0.0, data_plus_sumw2 = 0.0;
    double data_minus = 0.0, data_minus_sumw2 = 0.0;
    double ratio = 0.0, ratio_err = 0.0;
    double xsec = 0.0, xsec_err = 0.0;
    double xsec_sys_tgt = 0.0;
    double mean_q2_data = 0.0, mean_xb_data = 0.0, mean_tprime_data = 0.0;
    double mean_q2_sim = 0.0, mean_xb_sim = 0.0, mean_tprime_sim = 0.0;
    double mean_q2_xsec = 0.0, mean_xb_xsec = 0.0, mean_tprime_xsec = 0.0;
    int n_data = 0, n_sim = 0;
};

struct FourierFit {
    bool ok = false;
    bool absolute_xsec_fit = false;
    bool helicity_asymmetry_fit = false;
    std::vector<double> p;
    std::vector<double> perr;
    TMatrixD cov;
    double chi2 = 0.0;
    double ndf  = 0.0;
    double sigmaU = 0.0, sigmaU_err = 0.0;
    double sigmaTL = 0.0, sigmaTL_err = 0.0;
    double sigmaTT = 0.0, sigmaTT_err = 0.0;
    double sigmaTLp = 0.0, sigmaTLp_err = 0.0;

    // Custom copy constructor
    FourierFit(const FourierFit& other)
        : ok(other.ok), absolute_xsec_fit(other.absolute_xsec_fit), helicity_asymmetry_fit(other.helicity_asymmetry_fit),
          p(other.p), perr(other.perr), chi2(other.chi2), ndf(other.ndf),
          sigmaU(other.sigmaU), sigmaU_err(other.sigmaU_err),
          sigmaTL(other.sigmaTL), sigmaTL_err(other.sigmaTL_err),
          sigmaTT(other.sigmaTT), sigmaTT_err(other.sigmaTT_err),
          sigmaTLp(other.sigmaTLp), sigmaTLp_err(other.sigmaTLp_err)
    {
        cov.ResizeTo(other.cov.GetNrows(), other.cov.GetNcols());
        cov = other.cov;
    }

    // Custom assignment operator
    FourierFit& operator=(const FourierFit& other) {
        if (this != &other) {
            ok = other.ok;
            absolute_xsec_fit = other.absolute_xsec_fit;
            helicity_asymmetry_fit = other.helicity_asymmetry_fit;
            p = other.p;
            perr = other.perr;
            chi2 = other.chi2;
            ndf = other.ndf;
            sigmaU = other.sigmaU;
            sigmaU_err = other.sigmaU_err;
            sigmaTL = other.sigmaTL;
            sigmaTL_err = other.sigmaTL_err;
            sigmaTT = other.sigmaTT;
            sigmaTT_err = other.sigmaTT_err;
            sigmaTLp = other.sigmaTLp;
            sigmaTLp_err = other.sigmaTLp_err;
            cov.ResizeTo(other.cov.GetNrows(), other.cov.GetNcols());
            cov = other.cov;
        }
        return *this;
    }

    FourierFit() = default;
};

struct SliceResult {
    std::vector<PhiBin> phi;
    double sumw_data = 0.0, sumw2_data = 0.0;
    double sumw_sim  = 0.0, sumw2_sim  = 0.0;
    double mean_q2_data = 0.0, mean_xb_data = 0.0, mean_tprime_data = 0.0;
    double mean_q2_sim  = 0.0, mean_xb_sim  = 0.0, mean_tprime_sim  = 0.0;
    double mean_q2_abs  = 0.0, mean_xb_abs  = 0.0, mean_tprime_abs  = 0.0;
    bool has_model_xsec = false;
    double epsilon = 0.0;
    double gamma_flux = 0.0;
    FourierFit fit_ratio;
    FourierFit fit_xsec;
    FourierFit fit_asym; // used for helicity-odd diagnostic / TL' extraction
};

struct CutFlow {
    long long n_data_total = 0, n_data_pass = 0, n_data_inrange = 0;
    long long n_sim_total = 0, n_sim_pass = 0, n_sim_inrange = 0;
};

static double wrap_phi(double x) {
    double y = std::fmod(x, 2.0 * TMath::Pi());
    if (y < 0) y += 2.0 * TMath::Pi();
    // map exact 2pi to 0
    if (y >= 2.0 * TMath::Pi()) y = 0.0;
    return y;
}

static double clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

static double kallen_lambda_sqrt(double a, double b, double c) {
    double arg = (a - (b + c)) * (a - (b + c)) - 4.0 * b * c;
    return (arg > 0.0) ? std::sqrt(arg) : 0.0;
}

static double q2_xb_to_w2(double q2, double xb, double mp) {
    // W^2 = M^2 + Q^2(1/xB - 1)
    return mp * mp + q2 * (1.0 / xb - 1.0);
}

static double epsilon_virtual(double ebeam, double q2, double xb, double mp) {
    // Using y = nu/E = Q2/(2 M xB E)
    // epsilon = [1 - y - Q2/(4E^2)] / [1 - y + y^2/2 + Q2/(4E^2)]
    // This is the standard electron-scattering form for negligible electron mass.
    if (ebeam <= 0.0 || q2 <= 0.0 || xb <= 0.0) return 0.0;
    double y = q2 / (2.0 * mp * xb * ebeam);
    if (y <= 0.0) return 0.0;
    double e2 = ebeam * ebeam;
    double term = q2 / (4.0 * e2);
    double num = 1.0 - y - term;
    double den = 1.0 - y + 0.5 * y * y + term;
    if (den <= 0.0) return 0.0;
    return clamp(num / den, 0.0, 1.0);
}

static double virtual_photon_flux(double ebeam, double q2, double xb, double mp, double epsilon) {
    // Gamma_{gamma*} = alpha/(8pi) * Q2/(M^2 E^2) * (1-xB)/xB^3 * 1/(1-epsilon).
    // This is a flux factor, not an additional SIMC generation Jacobian.
    constexpr double alpha = 1.0 / 137.035999084;
    if (ebeam <= 0.0 || q2 <= 0.0 || xb <= 0.0 || mp <= 0.0) return 0.0;
    const double eps = clamp(epsilon, 0.0, 1.0 - 1e-12);
    return (alpha / (8.0 * TMath::Pi())) *
           (q2 / (mp * mp * ebeam * ebeam)) *
           ((1.0 - xb) / (xb * xb * xb)) *
           (1.0 / (1.0 - eps));
}

static std::vector<double> quantile_edges(std::vector<double> values, int nbins, double lo, double hi) {
    std::vector<double> edges(nbins + 1);
    edges.front() = lo;
    edges.back()  = hi;

    std::vector<double> clipped;
    clipped.reserve(values.size());
    for (double v : values) {
        if (std::isfinite(v) && v >= lo && v <= hi) clipped.push_back(v);
    }

    if (clipped.size() < static_cast<size_t>(nbins * 4)) {
        for (int i = 1; i < nbins; ++i) edges[i] = lo + (hi - lo) * double(i) / double(nbins);
        return edges;
    }

    std::sort(clipped.begin(), clipped.end());
    auto quant = [&](double p) {
        if (clipped.empty()) return lo;
        double pos = p * (clipped.size() - 1);
        size_t i0 = static_cast<size_t>(std::floor(pos));
        size_t i1 = std::min(i0 + 1, clipped.size() - 1);
        double f = pos - double(i0);
        return clipped[i0] * (1.0 - f) + clipped[i1] * f;
    };

    for (int i = 1; i < nbins; ++i) edges[i] = quant(double(i) / double(nbins));

    // Ensure monotonicity; if quantiles collapse, fall back to uniform bins.
    bool ok = true;
    for (int i = 1; i <= nbins; ++i) {
        if (!(edges[i] > edges[i - 1])) { ok = false; break; }
    }
    if (!ok) {
        for (int i = 1; i < nbins; ++i) edges[i] = lo + (hi - lo) * double(i) / double(nbins);
    }
    return edges;
}

static int find_bin(const std::vector<double>& edges, double x, bool periodic_phi = false) {
    if (edges.size() < 2) return -1;
    if (!std::isfinite(x)) return -1;
    if (periodic_phi) {
        x = wrap_phi(x);
    }

    if (x < edges.front() || x > edges.back()) return -1;
    if (x == edges.back()) return static_cast<int>(edges.size()) - 2;

    auto it = std::upper_bound(edges.begin(), edges.end(), x);
    int idx = static_cast<int>(it - edges.begin()) - 1;
    if (idx < 0 || idx >= static_cast<int>(edges.size()) - 1) return -1;
    return idx;
}

static std::vector<std::string> split_csv_line(const std::string&) { return {}; }

struct LinearFitResult {
    bool ok = false;
    std::vector<double> p;
    std::vector<double> perr;
    TMatrixD cov;
    double chi2 = 0.0;
    double ndf = 0.0;
};

static std::vector<double> phi_basis_means(double phi1, double phi2) {
    const double d = phi2 - phi1;
    if (d <= 0.0) return {1.0, 0.0, 0.0, 0.0};
    double c1 = (std::sin(phi2) - std::sin(phi1)) / d;
    double c2 = (std::sin(2.0 * phi2) - std::sin(2.0 * phi1)) / (2.0 * d);
    double s1 = (-std::cos(phi2) + std::cos(phi1)) / d;
    return {1.0, c1, c2, s1};
}

static std::array<double, 4> sigma_model_basis(double phi, double epsilon, double gamma_flux, int helicity_sign) {
    const double inv2pi = 1.0 / (2.0 * TMath::Pi());
    const double eps = clamp(epsilon, 0.0, 1.0);
    (void)gamma_flux;
    // const double gamma = (std::isfinite(gamma_flux) && gamma_flux > 0.0) ? gamma_flux : 0.0;
    // Gamma is intentionally not multiplied here. With w_base = full_weight/sigcm,
    // the retained factor siglab/sigcm already contains SIMC's davejac*gtpr*fac.
    // Reapplying this macro's Gamma would double-count the flux-like part.
    const double k_tl = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 + eps)));
    const double k_tlp = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 - eps)));
    return {
        inv2pi,
        inv2pi * k_tl * std::cos(phi),
        inv2pi * eps * std::cos(2.0 * phi),
        inv2pi * static_cast<double>(helicity_sign) * k_tlp * std::sin(phi)
    };
}

static double sigma_model_value(double phi,
                                double epsilon,
                                double gamma_flux,
                                double sigma_u,
                                double sigma_tl,
                                double sigma_tt,
                                double sigma_tlp,
                                int helicity_sign) {
    const auto b = sigma_model_basis(phi, epsilon, gamma_flux, helicity_sign);
    return b[0] * sigma_u + b[1] * sigma_tl + b[2] * sigma_tt + b[3] * sigma_tlp;
}

static void style_current_pad(double left = 0.13, double right = 0.04, double bottom = 0.12, double top = 0.08, bool grid = false) {
    if (!gPad) return;
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetBottomMargin(bottom);
    gPad->SetTopMargin(top);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(grid, grid);
}

static void style_axis(TAxis* ax, double title_size = 0.045, double label_size = 0.04, double title_offset = 1.1) {
    if (!ax) return;
    ax->SetTitleSize(title_size);
    ax->SetLabelSize(label_size);
    ax->SetTitleOffset(title_offset);
}

static void style_hist_axes(TH1* h, double title_size = 0.045, double label_size = 0.04) {
    if (!h) return;
    style_axis(h->GetXaxis(), title_size, label_size, 1.05);
    style_axis(h->GetYaxis(), title_size, label_size, 1.25);
    h->SetTitleSize(title_size);
}

static void style_graph_axes(TGraphErrors* g, double title_size = 0.045, double label_size = 0.04) {
    if (!g) return;
    style_axis(g->GetXaxis(), title_size, label_size, 1.05);
    style_axis(g->GetYaxis(), title_size, label_size, 1.25);
}

static TLegend* make_compact_legend(double x1, double y1, double x2, double y2, double text_size = 0.035) {
    auto leg = new TLegend(x1, y1, x2, y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(text_size);
    return leg;
}

static void draw_pad_message(const std::string& title, const std::string& line1, const std::string& line2 = "") {
    if (!gPad) return;
    gPad->Clear();
    style_current_pad();
    TLatex note;
    note.SetNDC(true);
    note.SetTextFont(42);
    note.SetTextSize(0.052);
    note.DrawLatex(0.12, 0.78, title.c_str());
    note.SetTextSize(0.038);
    note.DrawLatex(0.12, 0.62, line1.c_str());
    if (!line2.empty()) note.DrawLatex(0.12, 0.52, line2.c_str());
}

static FourierFit weighted_linear_fit(const std::vector<std::vector<double>>& X,
                                           const std::vector<double>& y,
                                           const std::vector<double>& ey) {
    FourierFit r;
    if (X.empty() || y.size() != X.size() || ey.size() != y.size()) return r;
    const int npar = static_cast<int>(X.front().size());
    const int npts = static_cast<int>(X.size());
    if (npts < npar) return r;

    std::vector<std::vector<double>> M(npar, std::vector<double>(npar, 0.0));
    std::vector<double> b(npar, 0.0);
    int used = 0;
    for (int i = 0; i < npts; ++i) {
        if (!(ey[i] > 0.0) || !std::isfinite(y[i])) continue;
        const double w = 1.0 / (ey[i] * ey[i]);
        ++used;
        for (int a = 0; a < npar; ++a) {
            b[a] += w * X[i][a] * y[i];
            for (int c = 0; c < npar; ++c) M[a][c] += w * X[i][a] * X[i][c];
        }
    }
    if (used < npar) return r;

    // Gauss-Jordan inversion of the normal matrix (npar <= 4 in this analysis).
    std::vector<std::vector<double>> A(npar, std::vector<double>(2 * npar, 0.0));
    for (int i = 0; i < npar; ++i) {
        for (int j = 0; j < npar; ++j) A[i][j] = M[i][j];
        A[i][i + npar] = 1.0;
    }
    for (int col = 0; col < npar; ++col) {
        int piv = col;
        double best = std::fabs(A[col][col]);
        for (int row = col + 1; row < npar; ++row) {
            double v = std::fabs(A[row][col]);
            if (v > best) { best = v; piv = row; }
        }
        if (best == 0.0) return r;
        if (piv != col) std::swap(A[piv], A[col]);

        double diag = A[col][col];
        for (int j = 0; j < 2 * npar; ++j) A[col][j] /= diag;
        for (int row = 0; row < npar; ++row) {
            if (row == col) continue;
            double f = A[row][col];
            for (int j = 0; j < 2 * npar; ++j) A[row][j] -= f * A[col][j];
        }
    }

    std::vector<std::vector<double>> inv(npar, std::vector<double>(npar, 0.0));
    for (int i = 0; i < npar; ++i)
        for (int j = 0; j < npar; ++j)
            inv[i][j] = A[i][j + npar];

    std::vector<double> p(npar, 0.0);
    for (int i = 0; i < npar; ++i) {
        for (int j = 0; j < npar; ++j) p[i] += inv[i][j] * b[j];
    }

    r.ok = true;
    r.p = p;
    r.perr.resize(npar, 0.0);
    // Safe assignment to TMatrixD cov with dimension check
    r.cov.ResizeTo(npar, npar);
    bool dim_ok = (r.cov.GetNrows() == npar && r.cov.GetNcols() == npar);
    if (!dim_ok) {
        std::cerr << "[DEBUG] TMatrixD cov dimension mismatch: "
                  << "r.cov is " << r.cov.GetNrows() << "x" << r.cov.GetNcols()
                  << ", expected " << npar << "x" << npar << std::endl;
    }
    for (int i = 0; i < npar && dim_ok; ++i) {
        for (int j = 0; j < npar; ++j) r.cov(i, j) = inv[i][j];
        r.perr[i] = (inv[i][i] > 0.0) ? std::sqrt(inv[i][i]) : 0.0;
    }

    double chi2 = 0.0;
    int nused = 0;
    for (int i = 0; i < npts; ++i) {
        if (!(ey[i] > 0.0) || !std::isfinite(y[i])) continue;
        double yhat = 0.0;
        for (int a = 0; a < npar; ++a) yhat += r.p[a] * X[i][a];
        double pull = (y[i] - yhat) / ey[i];
        chi2 += pull * pull;
        ++nused;
    }
    r.chi2 = chi2;
    r.ndf = std::max(0, nused - npar);
    return r;
}

class ExclPi0XSecAnalysis {

public:
    ~ExclPi0XSecAnalysis() {
        // Reset histograms before closing files
        h_q2_data.reset();
        h_q2_sim.reset();
        h_xb_data.reset();
        h_xb_sim.reset();
        h_tprime_data.reset();
        h_tprime_sim.reset();
        h_phi_data.reset();
        h_phi_sim.reset();
        h_q2_xb_data.reset();
        h_q2_xb_sim.reset();
        h_tprime_phi_data.reset();
        h_tprime_phi_sim.reset();
        // Close and null files
        cleanup();
        t_sim = nullptr;
        t_data = nullptr;
    }
    explicit ExclPi0XSecAnalysis(const AnalysisConfig& c) : cfg(c) {}
    void Run();

private:
    AnalysisConfig cfg;
    CutFlow cutflow;

    TFile* f_sim = nullptr;
    TFile* f_data = nullptr;
    TTree* t_sim = nullptr;
    TTree* t_data = nullptr;
    TFile* fout = nullptr;

    bool has_helicity = false;
    bool has_sim_helicity = false;
    bool has_model_xsec = false;
    std::string model_xsec_branch;

    std::string combined_pdf_path;
    bool combined_pdf_open = false;

    std::vector<double> phi_edges, tprime_edges, q2_edges, xb_edges;
    std::vector<std::vector<double>> xb_edges_by_q2;
    std::vector<SliceResult> slices;

    // Global QA histograms
    std::unique_ptr<TH1D> h_q2_data, h_q2_sim, h_xb_data, h_xb_sim, h_tprime_data, h_tprime_sim, h_phi_data, h_phi_sim;
    std::unique_ptr<TH2D> h_q2_xb_data, h_q2_xb_sim, h_tprime_phi_data, h_tprime_phi_sim;

    int slice_index(int it, int iq, int ix) const {
        return (it * cfg.n_q2 + iq) * cfg.n_xb + ix;
    }

    SliceResult& slice(int it, int iq, int ix) {
        return slices[slice_index(it, iq, ix)];
    }
    const SliceResult& slice(int it, int iq, int ix) const {
        return slices[slice_index(it, iq, ix)];
    }

    void log(const std::string& s) const { if (cfg.verbose) std::cout << "[INFO] " << s << "\n"; }
    void warn(const std::string& s) const { std::cerr << "[WARN] " << s << "\n"; }
    [[noreturn]] void die(const std::string& s) const { throw std::runtime_error(s); }

    void load_input();
    void detect_optional_branches();
    void build_binning();
    void init_storage();
    void fill_from_trees();
    void compute_ratios_and_xsec();
    void fit_slices();
    void make_global_plots();
    void make_epsilon_plots();
    void make_slice_plots();
    void make_sigma_vs_tprime_plots();
    void write_results();
    void write_csv();
    void write_slice_csv();
    void cleanup();
    void init_combined_pdf();
    void append_to_combined_pdf(TCanvas* c);
    void close_combined_pdf();
    std::string format_slice_bin_label(int it, int iq, int ix) const;
    void draw_slice_bin_label(int it, int iq, int ix, double y_ndc = 0.92) const;

    void fill_data_event(double q2, double t, double tmin, double xb, double phi, double pi0_weight, float scale, double charge_uC, double total_charge_uC, int is_exclusive, int helicity, bool use_helicity, double W);
    void fill_sim_event(float q2,
                        float t,
                        float tmin,
                        float xb,
                        float phi,
                        float full_weight,
                        int is_exclusive,
                        float model_xsec,
                        int helicity,
                        bool use_helicity,
                        float W);
    void finalize_slice_means(SliceResult& s);

    bool slice_passes_kin(const float q2, const float xb, const double tprime) const {
        return (q2 >= cfg.q2_min && q2 <= cfg.q2_max &&
                xb >= cfg.xb_min && xb <= cfg.xb_max &&
                tprime >= cfg.tprime_min && tprime <= cfg.tprime_max);
    }

    double calc_tprime(float t, float tmin) const { return static_cast<double>(t) - static_cast<double>(tmin); }

    double calc_w(float q2, float xb) const {
        return std::sqrt(std::max(0.0, q2_xb_to_w2(q2, xb, cfg.mp)));
    }

    void accumulate_global_histograms(float q2, float xb, double tprime, double phi, double weight_data, double weight_sim);
    void write_canvas_pdf_png(TCanvas* c, const std::string& base);
};

void ExclPi0XSecAnalysis::load_input() {
    f_sim = TFile::Open(cfg.simc_file.c_str(), "READ");
    if (!f_sim || f_sim->IsZombie()) die("Cannot open SIMC input file.");

    f_data = TFile::Open(cfg.data_file.c_str(), "READ");
    if (!f_data || f_data->IsZombie()) die("Cannot open data input file.");

    t_sim = dynamic_cast<TTree*>(f_sim->Get(cfg.simc_tree.c_str()));
    t_data = dynamic_cast<TTree*>(f_data->Get(cfg.data_tree.c_str()));
    if (!t_sim) die("Cannot find SIMC tree.");
    if (!t_data) die("Cannot find data tree.");
}

void ExclPi0XSecAnalysis::detect_optional_branches() {
    has_helicity = (t_data->GetBranch("helicity") != nullptr);
    has_sim_helicity = (t_sim->GetBranch("helicity") != nullptr);

    for (const auto& b : cfg.model_xsec_candidates) {
        if (t_sim->GetBranch(b.c_str())) {
            has_model_xsec = true;
            model_xsec_branch = b;
            break;
        }
    }

    if (!has_model_xsec) {
        die("No SIMC CM model cross-section branch 'sigcm' found. This extraction requires de-modeling with full_weight/sigcm.");
    }
    if (!has_helicity) warn("No helicity branch found in data tree. TL' will not be optimized.");
    if (has_helicity && !has_sim_helicity) {
        warn("Data helicity exists but SIMC helicity branch is missing. TL' optimization is disabled.");
    }
}

void ExclPi0XSecAnalysis::build_binning() {
    // Adaptive binning: quantile Q2 bins, then xB quantiles built independently in each Q2 bin.
    // This keeps Q2:xB 2D cells closer in occupancy than a single global rectangular split.
    std::vector<double> tprime_vals;
    std::vector<std::pair<double, double>> q2_xb_all;
    std::vector<std::pair<double, double>> q2_xb_data;
    auto collect = [&](TTree* t, bool is_sim) {
        if (!is_sim) {
            double q2 = 0, tval = 0, tmin = 0, xb = 0, phi = 0, pi0_weight = 0;
            float scale = 0;
            int is_exclusive = 0, helicity = 0;
            t->SetBranchAddress("Q2", &q2);
            t->SetBranchAddress("t", &tval);
            t->SetBranchAddress("tmin", &tmin);
            t->SetBranchAddress("xB", &xb);
            t->SetBranchAddress("phi", &phi);
            t->SetBranchAddress("pi0_weight", &pi0_weight);
            t->SetBranchAddress("scale", &scale);
            t->SetBranchAddress("is_exclusive", &is_exclusive);
            if (has_helicity) t->SetBranchAddress("helicity", &helicity);
            const Long64_t n = t->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                if (!is_exclusive) continue;
                if (!(q2 >= cfg.q2_min && q2 <= cfg.q2_max)) continue;
                if (!(xb >= cfg.xb_min && xb <= cfg.xb_max)) continue;
                double tprime = calc_tprime(tval, tmin);
                if (!(tprime >= cfg.tprime_min && tprime <= cfg.tprime_max)) continue;
                if (!std::isfinite(phi)) continue;
                tprime_vals.push_back(tprime);
                q2_xb_all.emplace_back(q2, xb);
                q2_xb_data.emplace_back(q2, xb);
            }
        } else {
            float q2 = 0, tval = 0, tmin = 0, xb = 0, phi = 0, full_weight = 0, model_xsec = 0;
            int is_exclusive = 0, helicity = 0;
            t->SetBranchAddress("Q2", &q2);
            t->SetBranchAddress("t", &tval);
            t->SetBranchAddress("tmin", &tmin);
            t->SetBranchAddress("xB", &xb);
            t->SetBranchAddress("phi", &phi);
            t->SetBranchAddress("full_weight", &full_weight);
            t->SetBranchAddress("is_exclusive", &is_exclusive);
            t->SetBranchAddress(model_xsec_branch.c_str(), &model_xsec);
            if (has_sim_helicity) t->SetBranchAddress("helicity", &helicity);
            const Long64_t n = t->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                t->GetEntry(i);
                if (!is_exclusive) continue;
                if (!(q2 >= cfg.q2_min && q2 <= cfg.q2_max)) continue;
                if (!(xb >= cfg.xb_min && xb <= cfg.xb_max)) continue;
                if (!std::isfinite(full_weight) || !std::isfinite(model_xsec) || std::fabs(model_xsec) < 1e-20f) continue;
                double tprime = calc_tprime(tval, tmin);
                if (!(tprime >= cfg.tprime_min && tprime <= cfg.tprime_max)) continue;
                if (!std::isfinite(phi)) continue;
                tprime_vals.push_back(tprime);
                q2_xb_all.emplace_back(q2, xb);
            }
        }
    };
    collect(t_data, false);
    collect(t_sim, true);

    const std::vector<std::pair<double, double>>& q2_xb_source = !q2_xb_data.empty() ? q2_xb_data : q2_xb_all;
    std::vector<double> q2_source;
    std::vector<double> xb_source;
    q2_source.reserve(q2_xb_source.size());
    xb_source.reserve(q2_xb_source.size());
    for (const auto& p : q2_xb_source) {
        q2_source.push_back(p.first);
        xb_source.push_back(p.second);
    }

    phi_edges = std::vector<double>(cfg.n_phi + 1, 0.0);
    for (int i = 0; i <= cfg.n_phi; ++i) phi_edges[i] = cfg.phi_min + (cfg.phi_max - cfg.phi_min) * double(i) / double(cfg.n_phi);

    tprime_edges = quantile_edges(tprime_vals, cfg.n_tprime, cfg.tprime_min, cfg.tprime_max);
    q2_edges = quantile_edges(q2_source, cfg.n_q2, cfg.q2_min, cfg.q2_max);

    xb_edges_by_q2.assign(cfg.n_q2, std::vector<double>(cfg.n_xb + 1, 0.0));
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        const double qlo = q2_edges[iq];
        const double qhi = q2_edges[iq + 1];
        std::vector<double> xb_in_q2;
        xb_in_q2.reserve(q2_xb_source.size() / std::max(1, cfg.n_q2));

        for (const auto& p : q2_xb_source) {
            const bool in_q2 = (p.first >= qlo) && ((iq == cfg.n_q2 - 1) ? (p.first <= qhi) : (p.first < qhi));
            if (in_q2) xb_in_q2.push_back(p.second);
        }
        xb_edges_by_q2[iq] = quantile_edges(xb_in_q2, cfg.n_xb, cfg.xb_min, cfg.xb_max);
    }

    // Keep a global xb edge vector for backward-compatible metadata products.
    xb_edges = quantile_edges(xb_source, cfg.n_xb, cfg.xb_min, cfg.xb_max);

    std::ostringstream source_note;
    source_note << "Computed adaptive Q2:xB bins using " << (!q2_xb_data.empty() ? "DATA" : "DATA+SIM") << " source.";
    log(source_note.str());
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        std::ostringstream oss;
        oss << "Q2 bin " << iq << " [" << q2_edges[iq] << ", " << q2_edges[iq + 1] << "] xB edges:";
        for (double e : xb_edges_by_q2[iq]) oss << " " << e;
        log(oss.str());
    }

    std::vector<std::vector<size_t>> occ(cfg.n_q2, std::vector<size_t>(cfg.n_xb, 0));
    for (const auto& p : q2_xb_source) {
        int iq = find_bin(q2_edges, p.first, false);
        if (iq < 0 || iq >= cfg.n_q2) continue;
        int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], p.second, false);
        if (ix < 0 || ix >= cfg.n_xb) continue;
        occ[static_cast<size_t>(iq)][static_cast<size_t>(ix)]++;
    }
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        std::ostringstream oss;
        oss << "2D occupancy row iq=" << iq << ":";
        for (int ix = 0; ix < cfg.n_xb; ++ix) oss << " " << occ[static_cast<size_t>(iq)][static_cast<size_t>(ix)];
        log(oss.str());
    }

    log("t' quantile edges: " + std::to_string(tprime_edges.front()) + " ... " + std::to_string(tprime_edges.back()));
    log("phi uniform bins: n=" + std::to_string(cfg.n_phi));
    // Binning: Q2 quantile bins, xB quantile bins within each Q2 bin, t' quantile, phi uniform.
}

void ExclPi0XSecAnalysis::init_storage() {
    slices.assign(cfg.n_tprime * cfg.n_q2 * cfg.n_xb, SliceResult{});
    for (auto& s : slices) s.phi.resize(cfg.n_phi);

    h_q2_data.reset(new TH1D("h_q2_data", "Data;Q^{2} [GeV^{2}];Weighted counts", 100, cfg.q2_min, cfg.q2_max));
    h_q2_sim.reset(new TH1D("h_q2_sim", "SIMC;Q^{2} [GeV^{2}];Weighted counts", 100, cfg.q2_min, cfg.q2_max));
    h_xb_data.reset(new TH1D("h_xb_data", "Data;x_{B};Weighted counts", 100, cfg.xb_min, cfg.xb_max));
    h_xb_sim.reset(new TH1D("h_xb_sim", "SIMC;x_{B};Weighted counts", 100, cfg.xb_min, cfg.xb_max));
    h_tprime_data.reset(new TH1D("h_tprime_data", "Data;t' [GeV^{2}];Weighted counts", 120, cfg.tprime_min, cfg.tprime_max));
    h_tprime_sim.reset(new TH1D("h_tprime_sim", "SIMC;t' [GeV^{2}];Weighted counts", 120, cfg.tprime_min, cfg.tprime_max));
    h_phi_data.reset(new TH1D("h_phi_data", "Data;#phi [rad];Weighted counts", cfg.n_phi, cfg.phi_min, cfg.phi_max));
    h_phi_sim.reset(new TH1D("h_phi_sim", "SIMC;#phi [rad];Weighted counts", cfg.n_phi, cfg.phi_min, cfg.phi_max));
    h_q2_xb_data.reset(new TH2D("h_q2_xb_data", "Data;Q^{2} [GeV^{2}];x_{B}", 100, cfg.q2_min, cfg.q2_max, 100, cfg.xb_min, cfg.xb_max));
    h_q2_xb_sim.reset(new TH2D("h_q2_xb_sim", "SIMC;Q^{2} [GeV^{2}];x_{B}", 100, cfg.q2_min, cfg.q2_max, 100, cfg.xb_min, cfg.xb_max));
    h_tprime_phi_data.reset(new TH2D("h_tprime_phi_data", "Data;t' [GeV^{2}];#phi [rad]", 100, cfg.tprime_min, cfg.tprime_max, cfg.n_phi, cfg.phi_min, cfg.phi_max));
    h_tprime_phi_sim.reset(new TH2D("h_tprime_phi_sim", "SIMC;t' [GeV^{2}];#phi [rad]", 100, cfg.tprime_min, cfg.tprime_max, cfg.n_phi, cfg.phi_min, cfg.phi_max));
}

void ExclPi0XSecAnalysis::accumulate_global_histograms(float q2, float xb, double tprime, double phi, double weight_data, double weight_sim) {
    h_q2_data->Fill(q2, weight_data);
    h_xb_data->Fill(xb, weight_data);
    h_tprime_data->Fill(tprime, weight_data);
    h_phi_data->Fill(phi, weight_data);
    h_q2_xb_data->Fill(q2, xb, weight_data);
    h_tprime_phi_data->Fill(tprime, phi, weight_data);

    h_q2_sim->Fill(q2, weight_sim);
    h_xb_sim->Fill(xb, weight_sim);
    h_tprime_sim->Fill(tprime, weight_sim);
    h_phi_sim->Fill(phi, weight_sim);
    h_q2_xb_sim->Fill(q2, xb, weight_sim);
    h_tprime_phi_sim->Fill(tprime, phi, weight_sim);
}

void ExclPi0XSecAnalysis::fill_data_event(double q2, double t, double tmin, double xb, double phi, double pi0_weight, float scale, double charge_uC, double total_charge_uC, int is_exclusive, int helicity, bool use_helicity, double W) {
    cutflow.n_data_total++;
    if (!is_exclusive) return;
    double tprime = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tprime)) return;
    if (!std::isfinite(phi)) return;
    if (!std::isfinite(pi0_weight) || !std::isfinite(scale)) return;

    cutflow.n_data_pass++;

    int it = find_bin(tprime_edges, tprime, false);
    int iq = find_bin(q2_edges, q2, false);
    if (iq < 0 || iq >= cfg.n_q2) return;
    int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], xb, false);
    int ip = find_bin(phi_edges, phi, true);
    if (it < 0 || ix < 0 || ip < 0) return;
    cutflow.n_data_inrange++;

    // NOTE: scripts/combine_analysis_branches.py defines scale as:
    //   scale = float(ps_val) / (float(cput_val) * float(charge_mC))
    // Therefore, scale already includes the run charge normalization.
    // To properly normalize per event, use (charge_uC/total_charge_uC) here.
    double charge_fraction = 1.0;
    if (std::isfinite(charge_uC) && charge_uC > 0.0 && std::isfinite(total_charge_uC) && total_charge_uC > 0.0) {
        charge_fraction = charge_uC / total_charge_uC;
    }
    double w = is_exclusive * pi0_weight * static_cast<double>(scale) * charge_fraction;

    PhiBin& pb = slice(it, iq, ix).phi[ip];
    pb.data += w;
    pb.data_sumw2 += w * w;
    pb.n_data += 1;
    pb.mean_q2_data += q2 * w;
    pb.mean_xb_data += xb * w;
    pb.mean_tprime_data += tprime * w;
    if (use_helicity && helicity > 0) {
        pb.data_plus += w;
        pb.data_plus_sumw2 += w * w;
    } else if (use_helicity && helicity < 0) {
        pb.data_minus += w;
        pb.data_minus_sumw2 += w * w;
    }

    SliceResult& s = slice(it, iq, ix);
    s.sumw_data += w;
    s.sumw2_data += w * w;
    s.mean_q2_data += q2 * w;
    s.mean_xb_data += xb * w;
    s.mean_tprime_data += tprime * w;
    accumulate_global_histograms(q2, xb, tprime, wrap_phi(phi), w, 0.0);
}

void ExclPi0XSecAnalysis::fill_sim_event(float q2,
                                         float t,
                                         float tmin,
                                         float xb,
                                         float phi,
                                         float full_weight,
                                         int is_exclusive,
                                         float model_xsec,
                                         int helicity,
                                         bool use_helicity,
                                         float W) {
    (void)W;
    cutflow.n_sim_total++;
    if (!is_exclusive) return;
    double tprime = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tprime)) return;
    if (!std::isfinite(phi)) return;
    if (!std::isfinite(full_weight) || !std::isfinite(model_xsec) || std::fabs(model_xsec) < 1e-20f) return;

    cutflow.n_sim_pass++;

    int it = find_bin(tprime_edges, tprime, false);
    int iq = find_bin(q2_edges, q2, false);
    if (iq < 0 || iq >= cfg.n_q2) return;
    int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], xb, false);
    int ip = find_bin(phi_edges, phi, true);
    if (it < 0 || ix < 0 || ip < 0) return;
    cutflow.n_sim_inrange++;

    // No bin migration is applied in this macro. Each reconstructed SIMC event
    // contributes only to its reconstructed analysis bin, so the de-modeled
    // SIMC sum is the diagonal effective response for that bin.
    const double base_w = static_cast<double>(full_weight) / static_cast<double>(model_xsec);
    if (!std::isfinite(base_w)) return;

    const double phiw = wrap_phi(phi);
    const double c1 = std::cos(phiw);
    const double c2 = std::cos(2.0 * phiw);
    const double s1 = std::sin(phiw);

    PhiBin& pb = slice(it, iq, ix).phi[ip];
    pb.sim_base += base_w;
    pb.sim_base_cos1 += base_w * c1;
    pb.sim_base_cos2 += base_w * c2;
    if (use_helicity && helicity > 0) {
        pb.sim_base_plus += base_w;
        pb.sim_base_cos1_plus += base_w * c1;
        pb.sim_base_cos2_plus += base_w * c2;
        pb.sim_base_sin1_plus += base_w * s1;
    } else if (use_helicity && helicity < 0) {
        pb.sim_base_minus += base_w;
        pb.sim_base_cos1_minus += base_w * c1;
        pb.sim_base_cos2_minus += base_w * c2;
        pb.sim_base_sin1_minus += base_w * s1;
    }

    pb.n_sim += 1;
    pb.mean_q2_sim += q2 * base_w;
    pb.mean_xb_sim += xb * base_w;
    pb.mean_tprime_sim += tprime * base_w;

    SliceResult& s = slice(it, iq, ix);
    s.sumw_sim += base_w;
    s.sumw2_sim += base_w * base_w;
    s.mean_q2_sim += q2 * base_w;
    s.mean_xb_sim += xb * base_w;
    s.mean_tprime_sim += tprime * base_w;
    accumulate_global_histograms(q2, xb, tprime, phiw, 0.0, base_w);
}

void ExclPi0XSecAnalysis::fill_from_trees() {
    // Data loop
    double dq2 = 0, dt = 0, dtmin = 0, dxb = 0, dphi = 0, dpi0_weight = 0, dW = 0;
    float dcharge_uC = 0;
    float dscale = 0;
    int dis_exclusive = 0, dhelicity = 0;
    int drun_number = 0;
    t_data->SetBranchAddress("Q2", &dq2);
    t_data->SetBranchAddress("t", &dt);
    t_data->SetBranchAddress("tmin", &dtmin);
    t_data->SetBranchAddress("xB", &dxb);
    t_data->SetBranchAddress("phi", &dphi);
    t_data->SetBranchAddress("pi0_weight", &dpi0_weight);
    t_data->SetBranchAddress("scale", &dscale);
    t_data->SetBranchAddress("is_exclusive", &dis_exclusive);
    t_data->SetBranchAddress("charge_uC", &dcharge_uC);
    t_data->SetBranchAddress("run_number", &drun_number);
    t_data->SetBranchAddress("W", &dW);
    if (has_helicity) t_data->SetBranchAddress("helicity", &dhelicity);

    // First pass: sum charge_uC per run_number
    // std::map<int, double> run_charge_map;
    const Long64_t ndata = t_data->GetEntries();
    std::unordered_set<int> seen_runs;
    double total_charge_uC = 0.0;

    for (Long64_t i = 0; i < ndata; ++i) {
        t_data->GetEntry(i);
        if (seen_runs.insert(drun_number).second && std::isfinite(dcharge_uC) && dcharge_uC > 0.0f) {
            total_charge_uC += dcharge_uC;
        }
    }

    if (!(total_charge_uC > 0.0)) {
        warn("Total charge_uC from combined data is invalid; falling back to neutral per-run charge factor.");
    }

    // Second pass: fill events using total_charge_uC
    for (Long64_t i = 0; i < ndata; ++i) {
        t_data->GetEntry(i);
        fill_data_event(dq2, dt, dtmin, dxb, dphi, dpi0_weight, dscale, dcharge_uC, total_charge_uC, dis_exclusive, dhelicity, has_helicity, dW);
    }

    // SIMC loop
    float sim_q2 = 0, sim_t = 0, sim_tmin = 0, sim_xb = 0, sim_phi = 0, full_weight = 0, sim_model_xsec = 0, sim_W = 0;
    int sim_is_exclusive = 0, sim_helicity = 0;
    t_sim->SetBranchAddress("Q2", &sim_q2);
    t_sim->SetBranchAddress("t", &sim_t);
    t_sim->SetBranchAddress("tmin", &sim_tmin);
    t_sim->SetBranchAddress("xB", &sim_xb);
    t_sim->SetBranchAddress("phi", &sim_phi);
    t_sim->SetBranchAddress("full_weight", &full_weight);
    t_sim->SetBranchAddress("is_exclusive", &sim_is_exclusive);
    t_sim->SetBranchAddress("W", &sim_W);
    t_sim->SetBranchAddress(model_xsec_branch.c_str(), &sim_model_xsec);
    if (has_sim_helicity) t_sim->SetBranchAddress("helicity", &sim_helicity);

    const Long64_t nsim = t_sim->GetEntries();
    for (Long64_t i = 0; i < nsim; ++i) {
        t_sim->GetEntry(i);
        fill_sim_event(sim_q2,
                       sim_t,
                       sim_tmin,
                       sim_xb,
                       sim_phi,
                       full_weight,
                       sim_is_exclusive,
                       sim_model_xsec,
                       sim_helicity,
                       has_sim_helicity,
                       sim_W);
    }

    // finalize all weighted means and per-phi errors
    for (auto& s : slices) {
        if (s.sumw_data > 0.0) {
            s.mean_q2_data /= s.sumw_data;
            s.mean_xb_data /= s.sumw_data;
            s.mean_tprime_data /= s.sumw_data;
        }
        if (s.sumw_sim > 0.0) {
            s.mean_q2_sim /= s.sumw_sim;
            s.mean_xb_sim /= s.sumw_sim;
            s.mean_tprime_sim /= s.sumw_sim;
        }
        for (auto& pb : s.phi) {
            if (pb.data > 0.0) {
                pb.mean_q2_data /= pb.data;
                pb.mean_xb_data /= pb.data;
                pb.mean_tprime_data /= pb.data;
            }
            if (pb.sim_base > 0.0) {
                pb.mean_q2_sim /= pb.sim_base;
                pb.mean_xb_sim /= pb.sim_base;
                pb.mean_tprime_sim /= pb.sim_base;
            }
        }
    }
}

void ExclPi0XSecAnalysis::compute_ratios_and_xsec() {
    const double contam = (std::isfinite(cfg.tgt_contam) && cfg.tgt_contam > 0.0) ? cfg.tgt_contam : 1.0;

    // In this method, model-xsec branch is mandatory for de-modeling SIMC.
    for (auto& s : slices) s.has_model_xsec = true;

    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                SliceResult& s = slice(it, iq, ix);

                double q2_for_eps = s.mean_q2_sim;
                if (!(std::isfinite(q2_for_eps) && q2_for_eps > 0.0)) q2_for_eps = s.mean_q2_data;
                if (!(std::isfinite(q2_for_eps) && q2_for_eps > 0.0)) q2_for_eps = 0.5 * (q2_edges[iq] + q2_edges[iq + 1]);

                double xb_for_eps = s.mean_xb_sim;
                if (!(std::isfinite(xb_for_eps) && xb_for_eps > 0.0)) xb_for_eps = s.mean_xb_data;
                if (!(std::isfinite(xb_for_eps) && xb_for_eps > 0.0)) xb_for_eps = 0.5 * (xb_edges_by_q2[iq][ix] + xb_edges_by_q2[iq][ix + 1]);

                s.epsilon = epsilon_virtual(cfg.ebeam, std::max(1e-12, q2_for_eps), std::max(1e-12, xb_for_eps), cfg.mp);
                s.gamma_flux = virtual_photon_flux(cfg.ebeam,
                                                   std::max(1e-12, q2_for_eps),
                                                   std::max(1e-12, xb_for_eps),
                                                   cfg.mp,
                                                   s.epsilon);
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    PhiBin& pb = s.phi[ip];
                    pb.data /= contam;
                    pb.data_sumw2 /= (contam * contam);
                    pb.data_plus /= contam;
                    pb.data_minus /= contam;
                    pb.data_plus_sumw2 /= (contam * contam);
                    pb.data_minus_sumw2 /= (contam * contam);

                    // These are solved in fit_slices() from the diagonal
                    // response-weighted basis K_r F_n; no N_rv migration matrix
                    // is formed in this version.
                    pb.sim = 0.0;
                    pb.sim_sumw2 = 0.0;
                    pb.ratio = 0.0;
                    pb.ratio_err = 0.0;
                    pb.xsec = 0.0;
                    pb.xsec_err = 0.0;
                    pb.xsec_sys_tgt = 0.0;
                    pb.mean_q2_xsec = pb.mean_q2_sim;
                    pb.mean_xb_xsec = pb.mean_xb_sim;
                    pb.mean_tprime_xsec = pb.mean_tprime_sim;
                }
            }
        }
    }
}

void ExclPi0XSecAnalysis::finalize_slice_means(SliceResult& s) {
    if (s.sumw_data > 0.0) {
        s.mean_q2_data /= s.sumw_data;
        s.mean_xb_data /= s.sumw_data;
        s.mean_tprime_data /= s.sumw_data;
    }
    if (s.sumw_sim > 0.0) {
        s.mean_q2_sim /= s.sumw_sim;
        s.mean_xb_sim /= s.sumw_sim;
        s.mean_tprime_sim /= s.sumw_sim;
    }
}

void ExclPi0XSecAnalysis::fit_slices() {
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                SliceResult& s = slice(it, iq, ix);
                s.fit_xsec = FourierFit();
                s.fit_ratio = FourierFit();
                s.fit_asym = FourierFit();

                const double eps = clamp(s.epsilon, 0.0, 1.0);
                // const double gamma = (std::isfinite(s.gamma_flux) && s.gamma_flux > 0.0) ? s.gamma_flux : 0.0;
                // Do not apply this macro's Gamma in the yield matrix. The
                // full_weight/sigcm SIMC base keeps siglab/sigcm =
                // davejac*gtpr*fac, so the needed flux/Jacobian convention is
                // already inside the per-event acceptance sum.
                const double inv2pi = 1.0 / (2.0 * TMath::Pi());
                const double k_tl = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 + eps)));
                const double k_tlp = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 - eps)));
                const bool fit_tlp = (has_helicity && has_sim_helicity);
                const int npar = fit_tlp ? 4 : 3;

                auto as_row = [&](const std::array<double, 4>& a) {
                    std::vector<double> v(npar, 0.0);
                    for (int i = 0; i < npar; ++i) v[i] = a[static_cast<size_t>(i)];
                    return v;
                };

                auto row_norm = [&](const std::array<double, 4>& a) {
                    double sabs = 0.0;
                    for (int i = 0; i < npar; ++i) sabs += std::fabs(a[static_cast<size_t>(i)]);
                    return sabs;
                };

                std::vector<std::vector<double>> X;
                std::vector<double> y;
                std::vector<double> ey;
                X.reserve(cfg.n_phi * (fit_tlp ? 3 : 1));
                y.reserve(cfg.n_phi * (fit_tlp ? 3 : 1));
                ey.reserve(cfg.n_phi * (fit_tlp ? 3 : 1));

                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    const PhiBin& pb = s.phi[ip];
                    std::array<double, 4> row_total = {
                        inv2pi * pb.sim_base,
                        inv2pi * k_tl * pb.sim_base_cos1,
                        inv2pi * eps * pb.sim_base_cos2,
                        inv2pi * k_tlp * (pb.sim_base_sin1_plus - pb.sim_base_sin1_minus)
                    };

                    const double ey_data = std::sqrt(std::max(0.0, pb.data_sumw2));
                    if (ey_data > 0.0 && row_norm(row_total) > 0.0 && std::isfinite(pb.data)) {
                        X.push_back(as_row(row_total));
                        y.push_back(pb.data);
                        ey.push_back(ey_data);
                    }

                    if (fit_tlp) {
                        std::array<double, 4> row_plus = {
                            inv2pi * pb.sim_base_plus,
                            inv2pi * k_tl * pb.sim_base_cos1_plus,
                            inv2pi * eps * pb.sim_base_cos2_plus,
                            inv2pi * k_tlp * pb.sim_base_sin1_plus
                        };
                        const double ey_plus = std::sqrt(std::max(0.0, pb.data_plus_sumw2));
                        if (ey_plus > 0.0 && row_norm(row_plus) > 0.0 && std::isfinite(pb.data_plus)) {
                            X.push_back(as_row(row_plus));
                            y.push_back(pb.data_plus);
                            ey.push_back(ey_plus);
                        }

                        std::array<double, 4> row_minus = {
                            inv2pi * pb.sim_base_minus,
                            inv2pi * k_tl * pb.sim_base_cos1_minus,
                            inv2pi * eps * pb.sim_base_cos2_minus,
                            -inv2pi * k_tlp * pb.sim_base_sin1_minus
                        };
                        const double ey_minus = std::sqrt(std::max(0.0, pb.data_minus_sumw2));
                        if (ey_minus > 0.0 && row_norm(row_minus) > 0.0 && std::isfinite(pb.data_minus)) {
                            X.push_back(as_row(row_minus));
                            y.push_back(pb.data_minus);
                            ey.push_back(ey_minus);
                        }
                    }
                }

                if (X.size() >= static_cast<size_t>(npar)) {
                    s.fit_xsec = weighted_linear_fit(X, y, ey);
                }

                if (s.fit_xsec.ok && s.fit_xsec.p.size() >= static_cast<size_t>(npar)) {
                    s.fit_xsec.absolute_xsec_fit = true;
                    s.fit_xsec.sigmaU = s.fit_xsec.p[0];
                    s.fit_xsec.sigmaU_err = (s.fit_xsec.perr.size() > 0) ? s.fit_xsec.perr[0] : 0.0;
                    s.fit_xsec.sigmaTL = s.fit_xsec.p[1];
                    s.fit_xsec.sigmaTL_err = (s.fit_xsec.perr.size() > 1) ? s.fit_xsec.perr[1] : 0.0;
                    s.fit_xsec.sigmaTT = s.fit_xsec.p[2];
                    s.fit_xsec.sigmaTT_err = (s.fit_xsec.perr.size() > 2) ? s.fit_xsec.perr[2] : 0.0;
                    s.fit_xsec.sigmaTLp = (npar > 3) ? s.fit_xsec.p[3] : 0.0;
                    s.fit_xsec.sigmaTLp_err = (npar > 3 && s.fit_xsec.perr.size() > 3) ? s.fit_xsec.perr[3] : 0.0;

                    if (fit_tlp) {
                        s.fit_asym.ok = true;
                        s.fit_asym.helicity_asymmetry_fit = true;
                        s.fit_asym.chi2 = s.fit_xsec.chi2;
                        s.fit_asym.ndf = s.fit_xsec.ndf;
                        s.fit_asym.sigmaTLp = s.fit_xsec.sigmaTLp;
                        s.fit_asym.sigmaTLp_err = s.fit_xsec.sigmaTLp_err;
                        s.fit_asym.p = {s.fit_xsec.sigmaTLp};
                        s.fit_asym.perr = {s.fit_xsec.sigmaTLp_err};
                    }
                }

                auto dot_model = [&](const std::array<double, 4>& a) {
                    double v = 0.0;
                    if (!s.fit_xsec.ok) return v;
                    const int nuse = std::min(npar, static_cast<int>(s.fit_xsec.p.size()));
                    for (int i = 0; i < nuse; ++i) v += a[static_cast<size_t>(i)] * s.fit_xsec.p[static_cast<size_t>(i)];
                    return v;
                };

                auto cov_variance = [&](const std::array<double, 4>& a) {
                    if (!s.fit_xsec.ok) return 0.0;
                    const int nrows = s.fit_xsec.cov.GetNrows();
                    const int ncols = s.fit_xsec.cov.GetNcols();
                    const int nuse = std::min(npar, std::min(nrows, ncols));
                    double var = 0.0;
                    for (int i = 0; i < nuse; ++i) {
                        for (int j = 0; j < nuse; ++j) {
                            var += a[static_cast<size_t>(i)] * s.fit_xsec.cov(i, j) * a[static_cast<size_t>(j)];
                        }
                    }
                    return std::max(0.0, var);
                };

                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    PhiBin& pb = s.phi[ip];

                    const std::array<double, 4> row_total = {
                        inv2pi * pb.sim_base,
                        inv2pi * k_tl * pb.sim_base_cos1,
                        inv2pi * eps * pb.sim_base_cos2,
                        inv2pi * k_tlp * (pb.sim_base_sin1_plus - pb.sim_base_sin1_minus)
                    };

                    pb.sim = dot_model(row_total);
                    pb.sim_sumw2 = cov_variance(row_total);

                    if (pb.sim > 0.0) {
                        pb.ratio = pb.data / pb.sim;
                        pb.ratio_err = std::sqrt(std::max(0.0, pb.data_sumw2)) / pb.sim;
                    } else {
                        pb.ratio = 0.0;
                        pb.ratio_err = 0.0;
                    }

                    const double phi_c = 0.5 * (phi_edges[ip] + phi_edges[ip + 1]);
                    if (s.fit_xsec.ok) {
                        pb.xsec = sigma_model_value(phi_c,
                                                    eps,
                                                    s.gamma_flux,
                                                    s.fit_xsec.sigmaU,
                                                    s.fit_xsec.sigmaTL,
                                                    s.fit_xsec.sigmaTT,
                                                    s.fit_xsec.sigmaTLp,
                                                    0);
                        const std::array<double, 4> basis_xsec = sigma_model_basis(phi_c, eps, s.gamma_flux, 0);
                        pb.xsec_err = std::sqrt(cov_variance(basis_xsec));
                    } else {
                        pb.xsec = 0.0;
                        pb.xsec_err = 0.0;
                    }

                    pb.xsec_sys_tgt = pb.xsec * (cfg.tgt_contam_err / std::max(1e-12, cfg.tgt_contam));
                    pb.mean_q2_xsec = pb.mean_q2_sim;
                    pb.mean_xb_xsec = pb.mean_xb_sim;
                    pb.mean_tprime_xsec = pb.mean_tprime_sim;
                }

                // Keep the ratio-Fourier diagnostic for continuity with previous outputs.
                std::vector<std::vector<double>> Xratio;
                std::vector<double> yratio;
                std::vector<double> eyratio;
                Xratio.reserve(cfg.n_phi);
                yratio.reserve(cfg.n_phi);
                eyratio.reserve(cfg.n_phi);
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    const PhiBin& pb = s.phi[ip];
                    if (!(pb.ratio_err > 0.0) || !std::isfinite(pb.ratio)) continue;
                    const std::vector<double> basis = phi_basis_means(phi_edges[ip], phi_edges[ip + 1]);
                    Xratio.push_back({basis[0], basis[1], basis[2]});
                    yratio.push_back(pb.ratio);
                    eyratio.push_back(pb.ratio_err);
                }
                if (Xratio.size() >= 3) {
                    s.fit_ratio = weighted_linear_fit(Xratio, yratio, eyratio);
                    if (s.fit_ratio.ok) {
                        s.fit_ratio.absolute_xsec_fit = false;
                        s.fit_ratio.p.resize(3);
                        s.fit_ratio.perr.resize(3);
                    }
                }
            }
        }
    }
}

void ExclPi0XSecAnalysis::write_canvas_pdf_png(TCanvas* c, const std::string& base) {
    if (cfg.write_pdf) c->SaveAs((base + ".pdf").c_str());
    append_to_combined_pdf(c);
    if (cfg.write_png) c->SaveAs((base + ".png").c_str());
}

std::string ExclPi0XSecAnalysis::format_slice_bin_label(int it, int iq, int ix) const {
    if (it < 0 || it >= cfg.n_tprime || iq < 0 || iq >= cfg.n_q2 || ix < 0 || ix >= cfg.n_xb) {
        return "Invalid bin index";
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3)
        << "it=" << it << " iq=" << iq << " ix=" << ix
        << " | Q^{2}[" << q2_edges[iq] << "," << q2_edges[iq + 1] << "]"
        << " x_{B}[" << xb_edges_by_q2[iq][ix] << "," << xb_edges_by_q2[iq][ix + 1] << "]"
        << " t'[" << tprime_edges[it] << "," << tprime_edges[it + 1] << "]";
    return oss.str();
}

void ExclPi0XSecAnalysis::draw_slice_bin_label(int it, int iq, int ix, double y_ndc) const {
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.026);
    lat.DrawLatex(0.11, y_ndc, format_slice_bin_label(it, iq, ix).c_str());
}

void ExclPi0XSecAnalysis::init_combined_pdf() {
    fs::path p(cfg.out_all_plots_pdf);
    if (p.is_relative()) p = fs::path(cfg.out_dir) / p;
    combined_pdf_path = p.string();
    combined_pdf_open = false;
}

void ExclPi0XSecAnalysis::append_to_combined_pdf(TCanvas* c) {
    if (!cfg.write_pdf || !c) return;
    if (combined_pdf_path.empty()) init_combined_pdf();
    if (combined_pdf_path.empty()) return;

    if (!combined_pdf_open) {
        c->SaveAs((combined_pdf_path + "[").c_str());
        combined_pdf_open = true;
    }
    c->SaveAs(combined_pdf_path.c_str());
}

void ExclPi0XSecAnalysis::close_combined_pdf() {
    if (!cfg.write_pdf || !combined_pdf_open || combined_pdf_path.empty()) return;
    TCanvas cclose("cclose_allplots", "cclose_allplots", 10, 10);
    cclose.SaveAs((combined_pdf_path + "]").c_str());
    combined_pdf_open = false;
    log("Wrote combined plot PDF: " + combined_pdf_path);
}

void ExclPi0XSecAnalysis::make_global_plots() {
    fs::create_directories(fs::path(cfg.out_dir) / "global");

    // Global overlays should use the same target-contamination normalization
    // as the slice-level extraction to avoid apparent data/sim mismatches.
    const double contam_scale = (std::isfinite(cfg.tgt_contam) && cfg.tgt_contam > 0.0)
                                    ? (1.0 / cfg.tgt_contam)
                                    : 1.0;
    auto clone_scaled_data_hist = [&](const TH1D* src, const char* name) {
        TH1D* h = dynamic_cast<TH1D*>(src->Clone(name));
        if (h) {
            h->SetDirectory(nullptr);
            h->Scale(contam_scale);
        }
        return std::unique_ptr<TH1D>(h);
    };

    auto h_q2_data_corr = clone_scaled_data_hist(h_q2_data.get(), "h_q2_data_corr");
    auto h_xb_data_corr = clone_scaled_data_hist(h_xb_data.get(), "h_xb_data_corr");
    auto h_tprime_data_corr = clone_scaled_data_hist(h_tprime_data.get(), "h_tprime_data_corr");
    auto h_phi_data_corr = clone_scaled_data_hist(h_phi_data.get(), "h_phi_data_corr");

    auto clone_area_scaled_sim_hist = [](const TH1D* src, const char* name, const TH1D* data_ref) {
        TH1D* h = dynamic_cast<TH1D*>(src->Clone(name));
        if (h) {
            h->SetDirectory(nullptr);
            const double sim_int = h->Integral();
            const double data_int = data_ref ? data_ref->Integral() : 0.0;
            if (sim_int > 0.0 && data_int > 0.0) h->Scale(data_int / sim_int);
        }
        return std::unique_ptr<TH1D>(h);
    };

    auto h_q2_sim_shape = clone_area_scaled_sim_hist(h_q2_sim.get(), "h_q2_sim_shape", h_q2_data_corr.get());
    auto h_xb_sim_shape = clone_area_scaled_sim_hist(h_xb_sim.get(), "h_xb_sim_shape", h_xb_data_corr.get());
    auto h_tprime_sim_shape = clone_area_scaled_sim_hist(h_tprime_sim.get(), "h_tprime_sim_shape", h_tprime_data_corr.get());
    auto h_phi_sim_shape = clone_area_scaled_sim_hist(h_phi_sim.get(), "h_phi_sim_shape", h_phi_data_corr.get());

    auto draw_shape_overlay = [](TH1D* hdata, TH1D* hsim, const char* data_label, const char* sim_label) {
        hdata->SetLineWidth(2);
        hdata->SetLineColor(kBlack);
        hsim->SetLineWidth(2);
        hsim->SetLineColor(kRed);
        hdata->GetYaxis()->SetTitle("Weighted counts (SIMC area norm.)");
        const double ymax = std::max(hdata->GetMaximum(), hsim->GetMaximum());
        hdata->SetMinimum(0.0);
        hdata->SetMaximum((ymax > 0.0) ? 1.18 * ymax : 1.0);
        hdata->Draw("hist");
        hsim->Draw("hist same");
        hdata->Draw("hist same");
        auto leg = make_compact_legend(0.50, 0.72, 0.88, 0.88, 0.038);
        leg->AddEntry(hdata, data_label, "l");
        leg->AddEntry(hsim, sim_label, "l");
        leg->Draw();
    };

    TCanvas c1("c1", "global", 1200, 900);
    c1.Divide(2,2);

    c1.cd(1);
    draw_shape_overlay(h_q2_data_corr.get(), h_q2_sim_shape.get(), "Data (tgt corrected)", "SIMC (area normalized)");

    c1.cd(2);
    draw_shape_overlay(h_xb_data_corr.get(), h_xb_sim_shape.get(), "Data (tgt corrected)", "SIMC (area normalized)");

    c1.cd(3);
    draw_shape_overlay(h_tprime_data_corr.get(), h_tprime_sim_shape.get(), "Data (tgt corrected)", "SIMC (area normalized)");

    c1.cd(4);
    draw_shape_overlay(h_phi_data_corr.get(), h_phi_sim_shape.get(), "Data (tgt corrected)", "SIMC (area normalized)");

    c1.Update();
    write_canvas_pdf_png(&c1, (fs::path(cfg.out_dir) / "global" / "global_1d_distributions").string());

    // --- DEBUG: Q2:xB bin grid overlay ---
    TCanvas c2("c2", "Q2:xB binning", 800, 700);
    h_q2_xb_data->Draw("COLZ");
    // Draw Q2 and xB bin edges
    for (size_t i = 1; i < q2_edges.size() - 1; ++i) {
        TLine* l = new TLine(q2_edges[i], cfg.xb_min, q2_edges[i], cfg.xb_max);
        l->SetLineColor(kBlue+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
    }
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        const double qlo = q2_edges[iq];
        const double qhi = q2_edges[iq + 1];
        const auto& xrow = xb_edges_by_q2[iq];
        for (size_t i = 1; i < xrow.size() - 1; ++i) {
            TLine* l = new TLine(qlo, xrow[i], qhi, xrow[i]);
            l->SetLineColor(kRed+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
        }
    }
    c2.Update();
    write_canvas_pdf_png(&c2, (fs::path(cfg.out_dir) / "global" / "q2_xb_binning_debug").string());

    // --- DEBUG: t' quantile bin edges ---
    TCanvas c3("c3", "t' binning", 800, 700);
    h_tprime_data_corr->GetYaxis()->SetTitle("Tgt-corrected weighted counts");
    h_tprime_data_corr->SetMinimum(0.0);
    h_tprime_data_corr->SetMaximum(1.15 * std::max(1e-12, h_tprime_data_corr->GetMaximum()));
    h_tprime_data_corr->Draw("hist");
    for (size_t i = 1; i < tprime_edges.size() - 1; ++i) {
        TLine* l = new TLine(tprime_edges[i], 0, tprime_edges[i], h_tprime_data_corr->GetMaximum());
        l->SetLineColor(kGreen+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
    }
    c3.Update();
    write_canvas_pdf_png(&c3, (fs::path(cfg.out_dir) / "global" / "tprime_binning_debug").string());

    // --- DEBUG: phi bin edges ---
    TCanvas c4("c4", "phi binning", 800, 700);
    h_phi_data_corr->GetYaxis()->SetTitle("Tgt-corrected weighted counts");
    h_phi_data_corr->SetMinimum(0.0);
    h_phi_data_corr->SetMaximum(1.15 * std::max(1e-12, h_phi_data_corr->GetMaximum()));
    h_phi_data_corr->Draw("hist");
    for (size_t i = 1; i < phi_edges.size() - 1; ++i) {
        TLine* l = new TLine(phi_edges[i], 0, phi_edges[i], h_phi_data_corr->GetMaximum());
        l->SetLineColor(kMagenta+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
    }
    c4.Update();
    write_canvas_pdf_png(&c4, (fs::path(cfg.out_dir) / "global" / "phi_binning_debug").string());

    // --- Existing 2D plots ---
    TCanvas c5("c5", "2d", 1200, 500);
    c5.Divide(2,1);
    c5.cd(1); h_q2_xb_data->Draw("COLZ");
    c5.cd(2); h_tprime_phi_data->Draw("COLZ");
    c5.Update();
    write_canvas_pdf_png(&c5, (fs::path(cfg.out_dir) / "global" / "data_occupancy_2d").string());

    TCanvas c6("c6", "2d_sim", 1200, 500);
    c6.Divide(2,1);
    c6.cd(1); h_q2_xb_sim->Draw("COLZ");
    c6.cd(2); h_tprime_phi_sim->Draw("COLZ");
    c6.Update();
    write_canvas_pdf_png(&c6, (fs::path(cfg.out_dir) / "global" / "simc_occupancy_2d").string());
}

void ExclPi0XSecAnalysis::make_epsilon_plots() {
    fs::create_directories(fs::path(cfg.out_dir) / "global");

    TCanvas c_eps_t("c_eps_t", "epsilon_vs_tprime", 1200, 900);
    TH1D hframe("hframe_eps_t", "Virtual photon #epsilon vs t';t' [GeV^{2}];#epsilon", 100, cfg.tprime_min, cfg.tprime_max);
    hframe.SetMinimum(0.0);
    hframe.SetMaximum(1.05);
    hframe.Draw("AXIS");

    std::vector<std::unique_ptr<TGraphErrors>> graphs;
    TLegend leg(0.52, 0.58, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 2, kOrange + 7, kCyan + 2};

    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix < cfg.n_xb; ++ix) {
            std::vector<double> xt, yeps, ey;
            xt.reserve(cfg.n_tprime);
            yeps.reserve(cfg.n_tprime);
            ey.reserve(cfg.n_tprime);

            for (int it = 0; it < cfg.n_tprime; ++it) {
                const SliceResult& s = slice(it, iq, ix);
                if (!std::isfinite(s.epsilon)) continue;
                xt.push_back(0.5 * (tprime_edges[it] + tprime_edges[it + 1]));
                yeps.push_back(s.epsilon);
                ey.push_back(0.0);
            }

            if (xt.empty()) continue;

            auto g = std::make_unique<TGraphErrors>(static_cast<int>(xt.size()), xt.data(), yeps.data(), nullptr, ey.data());
            int series = iq * cfg.n_xb + ix;
            g->SetLineColor(colors[series % 6]);
            g->SetMarkerColor(colors[series % 6]);
            g->SetMarkerStyle(20 + (series % 10));
            g->SetLineWidth(2);
            g->Draw("PL SAME");

            std::ostringstream lbl;
            lbl << "Q^{2}[" << std::fixed << std::setprecision(2) << q2_edges[iq] << "," << q2_edges[iq + 1]
                << "], x_{B}[" << xb_edges_by_q2[iq][ix] << "," << xb_edges_by_q2[iq][ix + 1] << "]";
            leg.AddEntry(g.get(), lbl.str().c_str(), "lp");
            graphs.emplace_back(std::move(g));
        }
    }

    leg.Draw();
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.12, 0.93, Form("E_{beam} = %.3f GeV", cfg.ebeam));
    c_eps_t.Update();
    write_canvas_pdf_png(&c_eps_t, (fs::path(cfg.out_dir) / "global" / "epsilon_vs_tprime_by_q2_xb").string());

    const int nslices = cfg.n_tprime * cfg.n_q2 * cfg.n_xb;
    TH1D h_eps_slice("h_eps_slice", "#epsilon by (t',Q^{2},x_{B}) slice;slice index;#epsilon", nslices, 0.5, nslices + 0.5);
    int bin = 1;
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                h_eps_slice.SetBinContent(bin, s.epsilon);
                std::ostringstream bl;
                bl << std::fixed << std::setprecision(2)
                   << "t[" << tprime_edges[it] << "," << tprime_edges[it + 1] << "] "
                   << "Q2[" << q2_edges[iq] << "," << q2_edges[iq + 1] << "] "
                   << "xB[" << xb_edges_by_q2[iq][ix] << "," << xb_edges_by_q2[iq][ix + 1] << "]";
                h_eps_slice.GetXaxis()->SetBinLabel(bin, bl.str().c_str());
                ++bin;
            }
        }
    }

    TCanvas c_eps_idx("c_eps_idx", "epsilon_by_slice", 1300, 700);
    h_eps_slice.SetMinimum(0.0);
    h_eps_slice.SetMaximum(1.05);
    h_eps_slice.SetStats(0);
    h_eps_slice.SetMarkerStyle(20);
    h_eps_slice.SetLineWidth(2);
    h_eps_slice.GetXaxis()->SetLabelSize(0.018);
    h_eps_slice.GetXaxis()->LabelsOption("v");
    h_eps_slice.Draw("P HIST");
    c_eps_idx.SetBottomMargin(0.22);
    c_eps_idx.Update();
    write_canvas_pdf_png(&c_eps_idx, (fs::path(cfg.out_dir) / "global" / "epsilon_by_slice_index").string());
}

void ExclPi0XSecAnalysis::make_slice_plots() {
    fs::create_directories(fs::path(cfg.out_dir) / "slices");
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                std::ostringstream tag;
                tag << "t" << it << "_q" << iq << "_x" << ix;
                TCanvas c(("c_"+tag.str()).c_str(), tag.str().c_str(), 1600, 1100);
                c.Divide(2, 2, 0.01, 0.02);

                std::vector<double> phi_c(cfg.n_phi), ydata(cfg.n_phi), yerr(cfg.n_phi), ysim(cfg.n_phi), ysimerr(cfg.n_phi), ratio(cfg.n_phi), ratioerr(cfg.n_phi);
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    phi_c[ip] = 0.5 * (phi_edges[ip] + phi_edges[ip+1]);
                    const PhiBin& pb = s.phi[ip];
                    ydata[ip] = pb.data; yerr[ip] = std::sqrt(std::max(0.0, pb.data_sumw2));
                    ysim[ip]  = pb.sim;  ysimerr[ip] = std::sqrt(std::max(0.0, pb.sim_sumw2));
                    ratio[ip] = pb.ratio; ratioerr[ip] = pb.ratio_err;
                }

                c.cd(1);
                style_current_pad(0.13, 0.04, 0.12, 0.10);
                auto gdata = new TGraphErrors(cfg.n_phi, phi_c.data(), ydata.data(), nullptr, yerr.data());
                auto gsim  = new TGraphErrors(cfg.n_phi, phi_c.data(), ysim.data(), nullptr, ysimerr.data());
                gdata->SetTitle("Data/SIMC yields;#phi [rad];Weighted counts");
                gdata->SetMarkerStyle(20); gdata->SetLineColor(kBlack); gdata->SetMarkerColor(kBlack);
                gsim->SetMarkerStyle(24); gsim->SetLineColor(kRed); gsim->SetMarkerColor(kRed);
                gdata->Draw("AP");
                style_graph_axes(gdata, 0.046, 0.04);
                gsim->Draw("P SAME");
                auto l1 = make_compact_legend(0.68, 0.76, 0.90, 0.89, 0.035);
                l1->AddEntry(gdata,"Data","p"); l1->AddEntry(gsim,"SIMC","p"); l1->Draw();
                draw_slice_bin_label(it, iq, ix);

                c.cd(2);
                style_current_pad(0.13, 0.04, 0.12, 0.10);
                auto grat = new TGraphErrors(cfg.n_phi, phi_c.data(), ratio.data(), nullptr, ratioerr.data());
                grat->SetTitle("Ratio data/SIMC;#phi [rad];Ratio");
                grat->SetMarkerStyle(20); grat->Draw("AP");
                style_graph_axes(grat, 0.046, 0.04);
                if (s.fit_ratio.ok) {
                    auto f = new TF1(("fr_"+tag.str()).c_str(), "[0] + [1]*cos(x) + [2]*cos(2*x)", cfg.phi_min, cfg.phi_max);
                    f->SetParameters(s.fit_ratio.p[0], s.fit_ratio.p[1], s.fit_ratio.p[2]);
                    f->SetLineColor(kRed);
                    f->Draw("same");
                }
                draw_slice_bin_label(it, iq, ix);

                c.cd(3);
                style_current_pad(0.13, 0.04, 0.12, 0.10);
                if (s.fit_xsec.ok && s.fit_xsec.p.size() >= 3) {
                    const bool fit_tlp = (has_helicity && has_sim_helicity && s.fit_xsec.p.size() >= 4);
                    const double eps = clamp(s.epsilon, 0.0, 1.0);
                    // const double gamma = (std::isfinite(s.gamma_flux) && s.gamma_flux > 0.0) ? s.gamma_flux : 0.0;
                    // Gamma is recorded for diagnostics, but not used in this
                    // CM-response overlay because full_weight/sigcm already
                    // carries SIMC's davejac*gtpr*fac through siglab/sigcm.
                    const double inv2pi = 1.0 / (2.0 * TMath::Pi());
                    const double k_lt = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 + eps)));
                    const double k_tlp = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 - eps)));

                    std::vector<double> xsec(cfg.n_phi), xsecerr(cfg.n_phi);
                    for (int ip = 0; ip < cfg.n_phi; ++ip) {
                        const PhiBin& pb = s.phi[ip];
                        xsec[ip] = pb.xsec;
                        xsecerr[ip] = pb.xsec_err;
                    }

                    const int ncurve = 361;
                    std::vector<double> ph_curve(ncurve);
                    std::vector<double> y_tot(ncurve), y_u(ncurve), y_lt(ncurve), y_tt(ncurve), y_tlp(ncurve, 0.0);
                    for (int i = 0; i < ncurve; ++i) {
                        const double ph = cfg.phi_min + (cfg.phi_max - cfg.phi_min) * (static_cast<double>(i) / static_cast<double>(ncurve - 1));
                        ph_curve[i] = ph;
                        y_u[i] = inv2pi * s.fit_xsec.sigmaU;
                        y_lt[i] = inv2pi * k_lt * s.fit_xsec.sigmaTL * std::cos(ph);
                        y_tt[i] = inv2pi * eps * s.fit_xsec.sigmaTT * std::cos(2.0 * ph);
                        if (fit_tlp) y_tlp[i] = inv2pi * k_tlp * s.fit_xsec.sigmaTLp * std::sin(ph);
                        y_tot[i] = y_u[i] + y_lt[i] + y_tt[i] + (fit_tlp ? y_tlp[i] : 0.0);
                    }

                    double ymin = std::numeric_limits<double>::infinity();
                    double ymax = -std::numeric_limits<double>::infinity();
                    for (int ip = 0; ip < cfg.n_phi; ++ip) {
                        if (!std::isfinite(xsec[ip])) continue;
                        ymin = std::min(ymin, xsec[ip] - std::max(0.0, xsecerr[ip]));
                        ymax = std::max(ymax, xsec[ip] + std::max(0.0, xsecerr[ip]));
                    }
                    for (int i = 0; i < ncurve; ++i) {
                        ymin = std::min({ymin, y_tot[i], y_u[i], y_lt[i], y_tt[i]});
                        ymax = std::max({ymax, y_tot[i], y_u[i], y_lt[i], y_tt[i]});
                        if (fit_tlp) {
                            ymin = std::min(ymin, y_tlp[i]);
                            ymax = std::max(ymax, y_tlp[i]);
                        }
                    }
                    if (!std::isfinite(ymin) || !std::isfinite(ymax)) {
                        ymin = -1.0;
                        ymax = 1.0;
                    }
                    const double span = std::max(1e-12, ymax - ymin);

                    auto gx = new TGraphErrors(cfg.n_phi, phi_c.data(), xsec.data(), nullptr, xsecerr.data());
                    gx->SetTitle("d#sigma/d#phi decomposition in slice;#phi [rad];Cross section");
                    gx->SetMarkerStyle(20);
                    gx->SetMarkerColor(kBlack);
                    gx->SetLineColor(kBlack);
                    gx->SetMinimum(ymin - 0.14 * span);
                    gx->SetMaximum(ymax + 0.32 * span);
                    gx->Draw("AP");
                    style_graph_axes(gx, 0.046, 0.04);

                    auto g_tot = new TGraphErrors(ncurve, ph_curve.data(), y_tot.data(), nullptr, nullptr);
                    auto g_u = new TGraphErrors(ncurve, ph_curve.data(), y_u.data(), nullptr, nullptr);
                    auto g_lt = new TGraphErrors(ncurve, ph_curve.data(), y_lt.data(), nullptr, nullptr);
                    auto g_tt = new TGraphErrors(ncurve, ph_curve.data(), y_tt.data(), nullptr, nullptr);
                    TGraphErrors* g_tlp = nullptr;
                    if (fit_tlp) g_tlp = new TGraphErrors(ncurve, ph_curve.data(), y_tlp.data(), nullptr, nullptr);

                    g_tot->SetLineColor(kBlue + 1); g_tot->SetLineWidth(3);
                    g_u->SetLineColor(kBlack); g_u->SetLineStyle(2); g_u->SetLineWidth(2);
                    g_lt->SetLineColor(kRed + 1); g_lt->SetLineStyle(7); g_lt->SetLineWidth(2);
                    g_tt->SetLineColor(kGreen + 2); g_tt->SetLineStyle(9); g_tt->SetLineWidth(2);
                    if (g_tlp) {
                        g_tlp->SetLineColor(kMagenta + 2);
                        g_tlp->SetLineStyle(5);
                        g_tlp->SetLineWidth(2);
                    }

                    g_tot->Draw("L SAME");
                    g_u->Draw("L SAME");
                    g_lt->Draw("L SAME");
                    g_tt->Draw("L SAME");
                    if (g_tlp) g_tlp->Draw("L SAME");

                    auto l3 = make_compact_legend(0.63, 0.58, 0.91, 0.88, 0.030);
                    l3->AddEntry(gx, "Extracted d#sigma/d#phi", "p");
                    l3->AddEntry(g_tot, "Total model", "l");
                    l3->AddEntry(g_u, "#sigma_{U} term", "l");
                    l3->AddEntry(g_lt, "#sigma_{LT} cos#phi term", "l");
                    l3->AddEntry(g_tt, "#sigma_{TT} cos2#phi term", "l");
                    if (g_tlp) l3->AddEntry(g_tlp, "#sigma_{TL'} sin#phi term", "l");
                    l3->Draw();

                    auto comp = new TPaveText(0.13, 0.73, 0.57, 0.84, "NDC");
                    comp->SetBorderSize(0);
                    comp->SetFillColorAlpha(kWhite, 0.85);
                    comp->SetTextAlign(12);
                    comp->SetTextFont(42);
                    comp->SetTextSize(0.022);
                    if (fit_tlp) {
                        comp->AddText(Form("#sigma_{U}=%.3g#pm%.2g, #sigma_{LT}=%.3g#pm%.2g",
                                           s.fit_xsec.sigmaU, std::max(0.0, s.fit_xsec.sigmaU_err),
                                           s.fit_xsec.sigmaTL, std::max(0.0, s.fit_xsec.sigmaTL_err)));
                        comp->AddText(Form("#sigma_{TT}=%.3g#pm%.2g, #sigma_{TL'}=%.3g#pm%.2g",
                                           s.fit_xsec.sigmaTT, std::max(0.0, s.fit_xsec.sigmaTT_err),
                                           s.fit_xsec.sigmaTLp, std::max(0.0, s.fit_xsec.sigmaTLp_err)));
                    } else {
                        comp->AddText(Form("#sigma_{U}=%.3g#pm%.2g, #sigma_{LT}=%.3g#pm%.2g",
                                           s.fit_xsec.sigmaU, std::max(0.0, s.fit_xsec.sigmaU_err),
                                           s.fit_xsec.sigmaTL, std::max(0.0, s.fit_xsec.sigmaTL_err)));
                        comp->AddText(Form("#sigma_{TT}=%.3g#pm%.2g",
                                           s.fit_xsec.sigmaTT, std::max(0.0, s.fit_xsec.sigmaTT_err)));
                    }
                    comp->Draw();
                } else {
                    draw_pad_message("Fit unavailable", "Not enough populated #phi bins in this slice.");
                }
                draw_slice_bin_label(it, iq, ix);

                c.cd(4);
                style_current_pad(0.13, 0.04, 0.12, 0.10);
                if (has_helicity) {
                    std::vector<double> asym(cfg.n_phi), asymerr(cfg.n_phi);
                    for (int ip = 0; ip < cfg.n_phi; ++ip) {
                        const PhiBin& pb = s.phi[ip];
                        double sum = pb.data_plus + pb.data_minus;
                        asym[ip] = (sum > 0.0) ? (pb.data_plus - pb.data_minus) / sum : 0.0;
                        asymerr[ip] = (sum > 0.0) ? std::sqrt(pb.data_plus_sumw2 + pb.data_minus_sumw2) / sum : 0.0;
                    }
                    auto ga = new TGraphErrors(cfg.n_phi, phi_c.data(), asym.data(), nullptr, asymerr.data());
                    ga->SetTitle("Helicity asymmetry diagnostic;#phi [rad];A_{h}");
                    ga->SetMarkerStyle(20);
                    ga->SetMarkerColor(kBlack);
                    ga->SetLineColor(kBlack);
                    ga->Draw("AP");
                    style_graph_axes(ga, 0.046, 0.04);

                    if (has_sim_helicity && s.fit_xsec.ok && s.fit_xsec.p.size() >= 4) {
                        const double eps = clamp(s.epsilon, 0.0, 1.0);
                        // const double gamma = (std::isfinite(s.gamma_flux) && s.gamma_flux > 0.0) ? s.gamma_flux : 0.0;
                        // Same convention as the fit rows: no extra Gamma on
                        // top of the full_weight/sigcm event base.
                        const double inv2pi = 1.0 / (2.0 * TMath::Pi());
                        const double k_tl = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 + eps)));
                        const double k_tlp = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 - eps)));

                        std::vector<double> asym_model(cfg.n_phi, 0.0);
                        for (int ip = 0; ip < cfg.n_phi; ++ip) {
                            const PhiBin& pb = s.phi[ip];
                            const double yp = inv2pi * (
                                s.fit_xsec.sigmaU * pb.sim_base_plus +
                                k_tl * s.fit_xsec.sigmaTL * pb.sim_base_cos1_plus +
                                eps * s.fit_xsec.sigmaTT * pb.sim_base_cos2_plus +
                                k_tlp * s.fit_xsec.sigmaTLp * pb.sim_base_sin1_plus
                            );
                            const double ym = inv2pi * (
                                s.fit_xsec.sigmaU * pb.sim_base_minus +
                                k_tl * s.fit_xsec.sigmaTL * pb.sim_base_cos1_minus +
                                eps * s.fit_xsec.sigmaTT * pb.sim_base_cos2_minus -
                                k_tlp * s.fit_xsec.sigmaTLp * pb.sim_base_sin1_minus
                            );
                            const double denom = yp + ym;
                            asym_model[ip] = (denom > 0.0) ? ((yp - ym) / denom) : 0.0;
                        }

                        auto gmodel = new TGraphErrors(cfg.n_phi, phi_c.data(), asym_model.data(), nullptr, nullptr);
                        gmodel->SetLineColor(kRed + 1);
                        gmodel->SetLineWidth(2);
                        gmodel->SetMarkerStyle(24);
                        gmodel->SetMarkerColor(kRed + 1);
                        gmodel->Draw("PL SAME");

                        auto l4 = make_compact_legend(0.56, 0.76, 0.90, 0.89, 0.033);
                        l4->AddEntry(ga, "Data A_{h}", "p");
                        l4->AddEntry(gmodel, "Optimized model A_{h}", "lp");
                        l4->Draw();
                    } else if (!has_sim_helicity) {
                        TLatex note;
                        note.SetNDC(true);
                        note.SetTextFont(42);
                        note.SetTextSize(0.034);
                        note.DrawLatex(0.13, 0.83, "SIMC helicity branch missing");
                        note.SetTextSize(0.030);
                        note.DrawLatex(0.13, 0.75, "No TL' optimization drawn.");
                    }
                } else {
                    draw_pad_message("No helicity branch", "TL' is not extracted for this input.");
                }
                draw_slice_bin_label(it, iq, ix);

                c.Update();
                write_canvas_pdf_png(&c, (fs::path(cfg.out_dir) / "slices" / ("slice_" + tag.str())).string());
            }
        }
    }
}

void ExclPi0XSecAnalysis::make_sigma_vs_tprime_plots() {
    fs::create_directories(fs::path(cfg.out_dir) / "sigma_vs_tprime");

    const bool have_tlp = (has_helicity && has_sim_helicity);

    auto update_range = [](double y, double ey, double& ymin, double& ymax) {
        if (!std::isfinite(y)) return;
        const double err = (std::isfinite(ey) && ey > 0.0) ? ey : 0.0;
        ymin = std::min(ymin, y - err);
        ymax = std::max(ymax, y + err);
    };

    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix < cfg.n_xb; ++ix) {
            std::vector<double> x, ex;
            std::vector<double> sigma_u, sigma_u_err;
            std::vector<double> sigma_tl, sigma_tl_err;
            std::vector<double> sigma_tt, sigma_tt_err;
            std::vector<double> sigma_tlp, sigma_tlp_err;

            x.reserve(cfg.n_tprime);
            ex.reserve(cfg.n_tprime);
            sigma_u.reserve(cfg.n_tprime);
            sigma_u_err.reserve(cfg.n_tprime);
            sigma_tl.reserve(cfg.n_tprime);
            sigma_tl_err.reserve(cfg.n_tprime);
            sigma_tt.reserve(cfg.n_tprime);
            sigma_tt_err.reserve(cfg.n_tprime);
            sigma_tlp.reserve(cfg.n_tprime);
            sigma_tlp_err.reserve(cfg.n_tprime);

            for (int it = 0; it < cfg.n_tprime; ++it) {
                const SliceResult& s = slice(it, iq, ix);
                if (!s.fit_xsec.ok) continue;

                const double tcenter = 0.5 * (tprime_edges[it] + tprime_edges[it + 1]);
                const double thalf = 0.5 * (tprime_edges[it + 1] - tprime_edges[it]);

                x.push_back(-tcenter);
                ex.push_back(std::fabs(thalf));

                sigma_u.push_back(s.fit_xsec.sigmaU);
                sigma_u_err.push_back(std::max(0.0, s.fit_xsec.sigmaU_err));
                sigma_tl.push_back(s.fit_xsec.sigmaTL);
                sigma_tl_err.push_back(std::max(0.0, s.fit_xsec.sigmaTL_err));
                sigma_tt.push_back(s.fit_xsec.sigmaTT);
                sigma_tt_err.push_back(std::max(0.0, s.fit_xsec.sigmaTT_err));
                sigma_tlp.push_back(s.fit_xsec.sigmaTLp);
                sigma_tlp_err.push_back(std::max(0.0, s.fit_xsec.sigmaTLp_err));
            }

            if (x.empty()) continue;

            std::ostringstream tag;
            tag << "q" << iq << "_x" << ix;
            TCanvas c(("c_sigma_vs_t_" + tag.str()).c_str(), "sigma_vs_tprime", 1100, 820);
            c.cd();
            style_current_pad(0.16, 0.04, 0.13, 0.18, true);
            double xmin = std::numeric_limits<double>::infinity();
            double xmax = -std::numeric_limits<double>::infinity();
            double ymin = std::numeric_limits<double>::infinity();
            double ymax = -std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < x.size(); ++i) {
                xmin = std::min(xmin, x[i] - ex[i]);
                xmax = std::max(xmax, x[i] + ex[i]);
                update_range(sigma_u[i], sigma_u_err[i], ymin, ymax);
                update_range(sigma_tl[i], sigma_tl_err[i], ymin, ymax);
                update_range(sigma_tt[i], sigma_tt_err[i], ymin, ymax);
                if (have_tlp) update_range(sigma_tlp[i], sigma_tlp_err[i], ymin, ymax);
            }

            if (!(std::isfinite(xmin) && std::isfinite(xmax) && std::isfinite(ymin) && std::isfinite(ymax))) {
                continue;
            }

            const double xspan = std::max(1e-6, xmax - xmin);
            const double yspan = std::max(1e-9, ymax - ymin);
            TH1D h_coeff_frame(("h_coeff_frame_" + tag.str()).c_str(),
                               ";-t' [GeV^{2}];d#sigma/dt [b/GeV^{2}]",
                               100,
                               xmin - 0.08 * xspan,
                               xmax + 0.08 * xspan);
            h_coeff_frame.SetMinimum(ymin - 0.18 * yspan);
            h_coeff_frame.SetMaximum(ymax + 0.24 * yspan);
            const Int_t old_max_digits = TGaxis::GetMaxDigits();
            TGaxis::SetMaxDigits(2);
            h_coeff_frame.Draw();
            style_hist_axes(&h_coeff_frame, 0.045, 0.04);

            auto g_u = new TGraphErrors(static_cast<int>(x.size()), x.data(), sigma_u.data(), ex.data(), sigma_u_err.data());
            auto g_tl = new TGraphErrors(static_cast<int>(x.size()), x.data(), sigma_tl.data(), ex.data(), sigma_tl_err.data());
            auto g_tt = new TGraphErrors(static_cast<int>(x.size()), x.data(), sigma_tt.data(), ex.data(), sigma_tt_err.data());
            g_u->SetLineColor(kBlack); g_u->SetMarkerColor(kBlack); g_u->SetMarkerStyle(20); g_u->SetLineWidth(2);
            g_tl->SetLineColor(kRed + 1); g_tl->SetMarkerColor(kRed + 1); g_tl->SetMarkerStyle(21); g_tl->SetLineWidth(2);
            g_tt->SetLineColor(kBlue + 1); g_tt->SetMarkerColor(kBlue + 1); g_tt->SetMarkerStyle(22); g_tt->SetLineWidth(2);

            g_u->Draw("P SAME");
            g_tl->Draw("P SAME");
            g_tt->Draw("P SAME");

            TGraphErrors* g_tlp = nullptr;
            if (have_tlp) {
                g_tlp = new TGraphErrors(static_cast<int>(x.size()), x.data(), sigma_tlp.data(), ex.data(), sigma_tlp_err.data());
                g_tlp->SetLineColor(kGreen + 2);
                g_tlp->SetMarkerColor(kGreen + 2);
                g_tlp->SetMarkerStyle(23);
                g_tlp->SetLineWidth(2);
                g_tlp->Draw("P SAME");
            }

            auto leg_coeff = make_compact_legend(0.70, 0.72, 0.90, 0.88, 0.040);
            leg_coeff->AddEntry(g_u, "#sigma_{U}", "p");
            leg_coeff->AddEntry(g_tl, "#sigma_{LT}", "p");
            leg_coeff->AddEntry(g_tt, "#sigma_{TT}", "p");
            if (g_tlp) leg_coeff->AddEntry(g_tlp, "#sigma_{TL'}", "p");
            leg_coeff->Draw();

            TLatex lat1;
            lat1.SetNDC(true);
            lat1.SetTextFont(42);
            lat1.SetTextAlign(22);
            lat1.SetTextSize(0.040);
            lat1.DrawLatex(0.50, 0.965, "Fitted response terms vs -t'");
            lat1.SetTextSize(0.030);
            lat1.DrawLatex(0.50, 0.925,
                           Form("Q^{2} #in [%.3f, %.3f] GeV^{2}, x_{B} #in [%.3f, %.3f]",
                                q2_edges[iq], q2_edges[iq + 1], xb_edges_by_q2[iq][ix], xb_edges_by_q2[iq][ix + 1]));
            lat1.SetTextSize(0.026);
            lat1.DrawLatex(0.50, 0.885, "No bin-migration correction; SIMC supplies diagonal effective response");

            c.Update();
            write_canvas_pdf_png(&c, (fs::path(cfg.out_dir) / "sigma_vs_tprime" / ("sigma_terms_vs_minus_tprime_" + tag.str())).string());
            TGaxis::SetMaxDigits(old_max_digits);
        }
    }
}

void ExclPi0XSecAnalysis::write_results() {
    fout->cd();
    // Save bin edges
    TVectorD vphi(phi_edges.size()), vt(tprime_edges.size()), vq(q2_edges.size()), vx(xb_edges.size());
    TMatrixD vxb2d(cfg.n_q2, cfg.n_xb + 1);
    for (size_t i = 0; i < phi_edges.size(); ++i) vphi[i] = phi_edges[i];
    for (size_t i = 0; i < tprime_edges.size(); ++i) vt[i] = tprime_edges[i];
    for (size_t i = 0; i < q2_edges.size(); ++i) vq[i] = q2_edges[i];
    for (size_t i = 0; i < xb_edges.size(); ++i) vx[i] = xb_edges[i];
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix <= cfg.n_xb; ++ix) {
            vxb2d(iq, ix) = xb_edges_by_q2[iq][ix];
        }
    }
    vphi.Write("phi_edges");
    vt.Write("tprime_edges");
    vq.Write("q2_edges");
    vx.Write("xb_edges");
    vxb2d.Write("xb_edges_by_q2");

    h_q2_data->Write();
    h_q2_sim->Write();
    h_xb_data->Write();
    h_xb_sim->Write();
    h_tprime_data->Write();
    h_tprime_sim->Write();
    h_phi_data->Write();
    h_phi_sim->Write();
    h_q2_xb_data->Write();
    h_q2_xb_sim->Write();
    h_tprime_phi_data->Write();
    h_tprime_phi_sim->Write();

    TParameter<Long64_t>("n_data_total", cutflow.n_data_total).Write();
    TParameter<Long64_t>("n_data_pass", cutflow.n_data_pass).Write();
    TParameter<Long64_t>("n_data_inrange", cutflow.n_data_inrange).Write();
    TParameter<Long64_t>("n_sim_total", cutflow.n_sim_total).Write();
    TParameter<Long64_t>("n_sim_pass", cutflow.n_sim_pass).Write();
    TParameter<Long64_t>("n_sim_inrange", cutflow.n_sim_inrange).Write();

    // write slice fit summaries as TObjString for portability
    std::ostringstream meta;
    meta << std::setprecision(10);
    meta << "tprime_edges:";
    for (double x : tprime_edges) meta << " " << x;
    meta << "\nq2_edges:";
    for (double x : q2_edges) meta << " " << x;
    meta << "\nxb_edges:";
    for (double x : xb_edges) meta << " " << x;
    meta << "\nxb_edges_by_q2:";
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        meta << "\n  iq=" << iq << ":";
        for (double x : xb_edges_by_q2[iq]) meta << " " << x;
    }
    meta << "\nphi_edges:";
    for (double x : phi_edges) meta << " " << x;
    meta << "\nmodel_xsec_branch=" << (has_model_xsec ? model_xsec_branch : "NONE");
    meta << "\nresponse_mode=diagonal_no_bin_migration";
    meta << "\nfit_basis=(1/(2pi))*{sigmaU,sqrt(2eps(1+eps))*sigmaLT*cos(phi),eps*sigmaTT*cos(2phi),h*sqrt(2eps(1-eps))*sigmaTLp*sin(phi)}";
    meta << "\ngamma_flux_usage=calculated_and_written_only_not_multiplied_in_fit";
    meta << "\ngamma_flux_note=full_weight_over_sigcm_retains_siglab_over_sigcm_equals_davejac_times_gtpr_times_fac";
    meta << "\nsim_weight_mode=full_weight_over_sigcm";
    meta << "\nhelicity_data=" << (has_helicity ? "yes" : "no");
    meta << "\nhelicity_sim=" << (has_sim_helicity ? "yes" : "no");
    TObjString(meta.str().c_str()).Write("analysis_metadata");
}

void ExclPi0XSecAnalysis::write_csv() {
    std::ofstream out(cfg.out_csv);
    out << std::setprecision(10);
    out << "it,iq,ix,ip,phi_lo,phi_hi,phi_center,q2_mean_data,xb_mean_data,tprime_mean_data,q2_mean_sim,xb_mean_sim,tprime_mean_sim,"
        << "epsilon,gamma_flux,data,data_err,sim,sim_err,ratio,ratio_err,xsec,xsec_err,xsec_sys_tgt\n";
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    const PhiBin& pb = s.phi[ip];
                    double phi_lo = phi_edges[ip];
                    double phi_hi = phi_edges[ip+1];
                    double phi_c = 0.5 * (phi_lo + phi_hi);
                    out << it << "," << iq << "," << ix << "," << ip << ","
                        << phi_lo << "," << phi_hi << "," << phi_c << ","
                        << pb.mean_q2_data << "," << pb.mean_xb_data << "," << pb.mean_tprime_data << ","
                        << pb.mean_q2_sim << "," << pb.mean_xb_sim << "," << pb.mean_tprime_sim << ","
                        << s.epsilon << "," << s.gamma_flux << ","
                        << pb.data << "," << std::sqrt(std::max(0.0, pb.data_sumw2)) << ","
                        << pb.sim << "," << std::sqrt(std::max(0.0, pb.sim_sumw2)) << ","
                        << pb.ratio << "," << pb.ratio_err << ","
                        << pb.xsec << "," << pb.xsec_err << "," << pb.xsec_sys_tgt << "\n";
                }
            }
        }
    }
}

void ExclPi0XSecAnalysis::write_slice_csv() {
    std::ofstream out(cfg.out_slice_csv);
    out << std::setprecision(10);
    out << "it,iq,ix,tprime_lo,tprime_hi,tprime_center,q2_lo,q2_hi,xb_lo,xb_hi,"
        << "mean_q2_data,mean_xb_data,mean_tprime_data,mean_q2_sim,mean_xb_sim,mean_tprime_sim,"
        << "epsilon,gamma_flux,has_model_xsec,sumw_data,sumw_sim,"
        << "fit_ratio_ok,fit_ratio_chi2,fit_ratio_ndf,fit_ratio_A,fit_ratio_Aerr,fit_ratio_B,fit_ratio_Berr,fit_ratio_C,fit_ratio_Cerr,"
        << "fit_xsec_ok,fit_xsec_chi2,fit_xsec_ndf,fit_xsec_sigmaU,fit_xsec_sigmaUerr,fit_xsec_sigmaTL,fit_xsec_sigmaTLerr,fit_xsec_sigmaTT,fit_xsec_sigmaTTerr,"
        << "fit_asym_ok,fit_asym_chi2,fit_asym_ndf,fit_asym_sigmaTLp,fit_asym_sigmaTLperr\n";
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                double tlo = tprime_edges[it], thi = tprime_edges[it+1];
                double qlo = q2_edges[iq], qhi = q2_edges[iq+1];
                double xlo = xb_edges_by_q2[iq][ix], xhi = xb_edges_by_q2[iq][ix+1];
                double tcenter = 0.5 * (tlo + thi);
                out << it << "," << iq << "," << ix << ","
                    << tlo << "," << thi << "," << tcenter << ","
                    << qlo << "," << qhi << ","
                    << xlo << "," << xhi << ","
                    << s.mean_q2_data << "," << s.mean_xb_data << "," << s.mean_tprime_data << ","
                    << s.mean_q2_sim << "," << s.mean_xb_sim << "," << s.mean_tprime_sim << ","
                    << s.epsilon << "," << s.gamma_flux << "," << (s.has_model_xsec ? 1 : 0) << ","
                    << s.sumw_data << "," << s.sumw_sim << ","
                    << (s.fit_ratio.ok ? 1 : 0) << "," << s.fit_ratio.chi2 << "," << s.fit_ratio.ndf << ","
                    << (s.fit_ratio.p.size() > 0 ? s.fit_ratio.p[0] : 0.0) << "," << (s.fit_ratio.perr.size() > 0 ? s.fit_ratio.perr[0] : 0.0) << ","
                    << (s.fit_ratio.p.size() > 1 ? s.fit_ratio.p[1] : 0.0) << "," << (s.fit_ratio.perr.size() > 1 ? s.fit_ratio.perr[1] : 0.0) << ","
                    << (s.fit_ratio.p.size() > 2 ? s.fit_ratio.p[2] : 0.0) << "," << (s.fit_ratio.perr.size() > 2 ? s.fit_ratio.perr[2] : 0.0) << ","
                    << (s.fit_xsec.ok ? 1 : 0) << "," << s.fit_xsec.chi2 << "," << s.fit_xsec.ndf << ","
                    << s.fit_xsec.sigmaU << "," << s.fit_xsec.sigmaU_err << ","
                    << s.fit_xsec.sigmaTL << "," << s.fit_xsec.sigmaTL_err << ","
                    << s.fit_xsec.sigmaTT << "," << s.fit_xsec.sigmaTT_err << ","
                    << (s.fit_asym.ok ? 1 : 0) << "," << s.fit_asym.chi2 << "," << s.fit_asym.ndf << ","
                    << s.fit_asym.sigmaTLp << "," << s.fit_asym.sigmaTLp_err << "\n";
            }
        }
    }
}

void ExclPi0XSecAnalysis::cleanup() {

    close_combined_pdf();

    if (fout) { fout->Write(); fout->Close(); fout = nullptr; }
    if (f_sim) { f_sim->Close(); f_sim = nullptr; }
    if (f_data) { f_data->Close(); f_data = nullptr; }
}

void ExclPi0XSecAnalysis::Run() {
    fs::create_directories(cfg.out_dir);
    fs::create_directories(fs::path(cfg.out_dir) / "global");
    fs::create_directories(fs::path(cfg.out_dir) / "slices");
    init_combined_pdf();
    load_input();
    detect_optional_branches();
    build_binning();
    init_storage();
    fill_from_trees();
    compute_ratios_and_xsec();
    fit_slices();

    fout = TFile::Open(cfg.out_root.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) die("Cannot open output ROOT file.");

    make_global_plots();
    make_epsilon_plots();
    make_slice_plots();
    make_sigma_vs_tprime_plots();
    close_combined_pdf();
    write_results();
    write_csv();
    write_slice_csv();
    // cleanup() is now handled by the destructor; do not call here to avoid double-free.

    log("Analysis complete.");
}

static std::string getenv_or_default(const char* key, const std::string& fallback = "") {
    const char* v = std::getenv(key);
    if (!v || !v[0]) return fallback;
    return std::string(v);
}

static std::string trim_copy_cli(const std::string& s) {
    size_t first = 0;
    while (first < s.size() && std::isspace(static_cast<unsigned char>(s[first]))) ++first;
    size_t last = s.size();
    while (last > first && std::isspace(static_cast<unsigned char>(s[last - 1]))) --last;
    return s.substr(first, last - first);
}

static std::string sanitize_token_cli(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-') out.push_back(c);
        else out.push_back('_');
    }
    return out;
}

static std::string require_next_arg(int& i, int argc, char** argv, const std::string& opt_name) {
    if (i + 1 >= argc) throw std::runtime_error("Missing value for option " + opt_name);
    ++i;
    return std::string(argv[i]);
}

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << R"( [options]

Primary options:
  --kin <Kin_old>                 Kinematic setting (used to derive standard paths)
  --target <name>                 Combined target token (default: LH2)
  --output-base <path>            Output base directory (default: $NPS_OUTPUT_BASE or output)
  --root-dir <path>               Root directory for combined/sim files (default from kin)
  --data-file <path>              Combined data ROOT file
  --sim-file <path>               Simulation ROOT file
  --out-dir <path>                Output directory for xsec products
  --out-root <path>               Output ROOT file
  --out-csv <path>                Output CSV summary
  --out-slice-csv <path>          Output per-slice CSV summary
  --all-plots-pdf <path>          Combined plots PDF path

Binning/range options:
  --n-phi <int>                   Phi bins
  --n-tprime <int>                t' bins
  --n-q2 <int>                    Q2 quantile bins
  --n-xb <int>                    xB quantile bins
  --phi-min <float>               Phi min
  --phi-max <float>               Phi max
  --tprime-min <float>            t' min
  --tprime-max <float>            t' max
  --q2-min <float>                Q2 min
  --q2-max <float>                Q2 max
  --xb-min <float>                xB min
  --xb-max <float>                xB max

Behavior options:
  --ebeam <float>                 Beam energy
  --quiet                         Disable verbose logging
  --no-diagnostics                Disable diagnostic outputs
  --no-png                        Disable PNG output
  --no-pdf                        Disable PDF output
  --help                          Show this message

Compatibility:
  positional arg1                 out_root (legacy)
  positional arg2                 out_csv (legacy)
)";
}

int main(int argc, char** argv) {
    try {
        TH1::SetDefaultSumw2(kTRUE);
        gStyle->SetOptStat(0);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTitleSize(0.045, "XYZ");
        gStyle->SetLabelSize(0.04, "XYZ");
        gStyle->SetTitleOffset(1.15, "X");
        gStyle->SetTitleOffset(1.25, "Y");
        gStyle->SetLegendBorderSize(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        AnalysisConfig cfg;

        std::string cli_kin = getenv_or_default("NPS_XSEC_KIN", "");
        std::string cli_target = getenv_or_default("NPS_XSEC_TARGET", "LH2");
        std::string cli_output_base = getenv_or_default("NPS_OUTPUT_BASE", "output");
        std::string cli_root_dir = getenv_or_default("NPS_XSEC_ROOT_DIR", "");

        bool data_file_set = false;
        bool sim_file_set = false;
        bool out_dir_set = false;
        bool out_root_set = false;
        bool out_csv_set = false;
        bool out_slice_csv_set = false;
        bool out_all_plots_pdf_set = false;

        auto apply_env_path = [&](const char* key, std::string& dst, bool& set_flag) {
            std::string v = getenv_or_default(key, "");
            if (!v.empty()) {
                dst = v;
                set_flag = true;
            }
        };

        apply_env_path("NPS_XSEC_DATA_FILE", cfg.data_file, data_file_set);
        apply_env_path("NPS_XSEC_SIM_FILE", cfg.simc_file, sim_file_set);
        apply_env_path("NPS_XSEC_OUT_DIR", cfg.out_dir, out_dir_set);
        apply_env_path("NPS_XSEC_OUT_ROOT", cfg.out_root, out_root_set);
        apply_env_path("NPS_XSEC_OUT_CSV", cfg.out_csv, out_csv_set);
        apply_env_path("NPS_XSEC_OUT_SLICE_CSV", cfg.out_slice_csv, out_slice_csv_set);
        apply_env_path("NPS_XSEC_ALL_PLOTS_PDF", cfg.out_all_plots_pdf, out_all_plots_pdf_set);

        int positional_idx = 0;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") { print_usage(argv[0]); return 0; }
            if (arg == "--kin") { cli_kin = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--target") { cli_target = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--output-base") { cli_output_base = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--root-dir") { cli_root_dir = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--data-file") { cfg.data_file = require_next_arg(i, argc, argv, arg); data_file_set = true; continue; }
            if (arg == "--sim-file") { cfg.simc_file = require_next_arg(i, argc, argv, arg); sim_file_set = true; continue; }
            if (arg == "--out-dir") { cfg.out_dir = require_next_arg(i, argc, argv, arg); out_dir_set = true; continue; }
            if (arg == "--out-root") { cfg.out_root = require_next_arg(i, argc, argv, arg); out_root_set = true; continue; }
            if (arg == "--out-csv") { cfg.out_csv = require_next_arg(i, argc, argv, arg); out_csv_set = true; continue; }
            if (arg == "--out-slice-csv") { cfg.out_slice_csv = require_next_arg(i, argc, argv, arg); out_slice_csv_set = true; continue; }
            if (arg == "--all-plots-pdf") { cfg.out_all_plots_pdf = require_next_arg(i, argc, argv, arg); out_all_plots_pdf_set = true; continue; }
            if (arg == "--n-phi") { cfg.n_phi = std::stoi(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--n-tprime") { cfg.n_tprime = std::stoi(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--n-q2") { cfg.n_q2 = std::stoi(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--n-xb") { cfg.n_xb = std::stoi(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--phi-min") { cfg.phi_min = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--phi-max") { cfg.phi_max = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--tprime-min") { cfg.tprime_min = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--tprime-max") { cfg.tprime_max = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--q2-min") { cfg.q2_min = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--q2-max") { cfg.q2_max = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--xb-min") { cfg.xb_min = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--xb-max") { cfg.xb_max = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--ebeam") { cfg.ebeam = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--quiet") { cfg.verbose = false; continue; }
            if (arg == "--no-diagnostics") { cfg.diagnostics = false; continue; }
            if (arg == "--no-png") { cfg.write_png = false; continue; }
            if (arg == "--no-pdf") { cfg.write_pdf = false; continue; }

            if (!arg.empty() && arg[0] == '-') throw std::runtime_error("Unknown option: " + arg);
            if (positional_idx == 0) { cfg.out_root = arg; out_root_set = true; }
            else if (positional_idx == 1) { cfg.out_csv = arg; out_csv_set = true; }
            else throw std::runtime_error("Unexpected positional argument: " + arg);
            ++positional_idx;
        }

        cli_kin = trim_copy_cli(cli_kin);
        cli_target = trim_copy_cli(cli_target);
        cli_output_base = trim_copy_cli(cli_output_base);
        cli_root_dir = trim_copy_cli(cli_root_dir);
        if (cli_target.empty()) cli_target = "LH2";
        if (cli_output_base.empty()) cli_output_base = "output";

        if (!cli_kin.empty()) {
            const std::string kin_safe = sanitize_token_cli(cli_kin);
            if (kin_safe.empty()) throw std::runtime_error("Invalid --kin value after sanitization: '" + cli_kin + "'");
            if (cli_root_dir.empty()) cli_root_dir = (fs::path(cli_output_base) / kin_safe / "root").string();
            if (!out_dir_set) { cfg.out_dir = (fs::path(cli_output_base) / kin_safe / "xsec").string(); out_dir_set = true; }
        }

        if (!cli_root_dir.empty()) {
            if (!data_file_set) {
                cfg.data_file = (fs::path(cli_root_dir) / ("combined_branches_" + sanitize_token_cli(cli_target) + ".root")).string();
                data_file_set = true;
            }
            if (!sim_file_set) {
                const fs::path smeared = fs::path(cli_root_dir) / "simc_pi0_analysis_output_smeared.root";
                const fs::path unsmeared = fs::path(cli_root_dir) / "simc_pi0_analysis_output.root";
                cfg.simc_file = fs::exists(smeared) ? smeared.string() : unsmeared.string();
                sim_file_set = true;
            }
            if (!out_dir_set) { cfg.out_dir = (fs::path(cli_root_dir).parent_path() / "xsec").string(); out_dir_set = true; }
        }

        if (!out_dir_set || trim_copy_cli(cfg.out_dir).empty()) cfg.out_dir = "output_pi0_xsec_no_simc_model";
        if (!out_root_set || trim_copy_cli(cfg.out_root).empty()) cfg.out_root = (fs::path(cfg.out_dir) / "excl_xsec_pi0_analysis_no_simc_model_output.root").string();
        if (!out_csv_set || trim_copy_cli(cfg.out_csv).empty()) cfg.out_csv = (fs::path(cfg.out_dir) / "excl_xsec_pi0_analysis_no_simc_model_summary.csv").string();
        if (!out_slice_csv_set || trim_copy_cli(cfg.out_slice_csv).empty()) cfg.out_slice_csv = (fs::path(cfg.out_dir) / "excl_xsec_pi0_analysis_no_simc_model_slice_summary.csv").string();
        if (!out_all_plots_pdf_set || trim_copy_cli(cfg.out_all_plots_pdf).empty()) cfg.out_all_plots_pdf = (fs::path(cfg.out_dir) / "all_generated_plots_no_simc_model.pdf").string();

        if (trim_copy_cli(cfg.data_file).empty()) throw std::runtime_error("Data file path is empty. Provide --data-file or --kin/--root-dir.");
        if (trim_copy_cli(cfg.simc_file).empty()) throw std::runtime_error("Simulation file path is empty. Provide --sim-file or --kin/--root-dir.");
        if (!fs::exists(cfg.data_file)) throw std::runtime_error("Data ROOT file not found: " + cfg.data_file);
        if (!fs::exists(cfg.simc_file)) throw std::runtime_error("Simulation ROOT file not found: " + cfg.simc_file);
        if (cfg.n_phi <= 0 || cfg.n_tprime <= 0 || cfg.n_q2 <= 0 || cfg.n_xb <= 0) {
            throw std::runtime_error("Binning values must be positive (n_phi, n_tprime, n_q2, n_xb).");
        }
        if (!(cfg.phi_min < cfg.phi_max) || !(cfg.tprime_min < cfg.tprime_max) ||
            !(cfg.q2_min < cfg.q2_max) || !(cfg.xb_min < cfg.xb_max)) {
            throw std::runtime_error("Invalid analysis range: each min must be strictly smaller than max.");
        }

        std::cout << "XSec input configuration:" << std::endl
                  << "  data file: " << cfg.data_file << std::endl
                  << "  sim file:  " << cfg.simc_file << std::endl
                  << "  out dir:   " << cfg.out_dir << std::endl
                  << "  out root:  " << cfg.out_root << std::endl
                  << "  out csv:   " << cfg.out_csv << std::endl
                  << "  out slice: " << cfg.out_slice_csv << std::endl
                  << "  plots pdf: " << cfg.out_all_plots_pdf << std::endl;

        ExclPi0XSecAnalysis ana(cfg);
        ana.Run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[FATAL] " << e.what() << std::endl;
        return 1;
    }
}
