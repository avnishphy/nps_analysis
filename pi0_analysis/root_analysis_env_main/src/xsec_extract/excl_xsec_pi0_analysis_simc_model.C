// Iterative SIMC ratio extraction in physical-t bins.
// Ysim(p)=C*sum[(full_weight/sigcm)*model(vertex,p)].
// g++ -O2 -std=c++17 -o excl_xsec_simc_model excl_xsec_pi0_analysis_simc_model.C `root-config --cflags --libs`

#include <TFile.h>
#include <TDirectory.h>
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
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TDecompLU.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TMath.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TPad.h>
#include <TAxis.h>

#ifdef NPS_ENABLE_PARTONS
#include "partons_pi0_projection.h"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
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
#include "simc_pi0_model.h"
#include "simc_pi0_reweight.h"
#include "xsec_config.h"

namespace fs = std::filesystem;

static constexpr const char* kSigmaTLpStatus =
    "unavailable_missing_helicity_luminosities_and_beam_polarization";

struct PhiBin {
    double data = 0.0, data_sumw2 = 0.0;
    double sim  = 0.0, sim_sumw2  = 0.0;
    double model_wsum = 0.0, model_xsec_wsum = 0.0;
    double model_phi_center = 0.0; // fitted model at phi bin center
    double ratio_before = std::numeric_limits<double>::quiet_NaN();
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
    std::vector<double> p;
    std::vector<double> perr;
    TMatrixD cov;
    double chi2 = 0.0;
    double ndf  = 0.0;
    double sigmaU = 0.0, sigmaU_err = 0.0;
    double sigmaTL = 0.0, sigmaTL_err = 0.0;
    double sigmaTT = 0.0, sigmaTT_err = 0.0;

    // Custom copy constructor
    FourierFit(const FourierFit& other)
        : ok(other.ok), absolute_xsec_fit(other.absolute_xsec_fit),
          p(other.p), perr(other.perr), chi2(other.chi2), ndf(other.ndf),
          sigmaU(other.sigmaU), sigmaU_err(other.sigmaU_err),
          sigmaTL(other.sigmaTL), sigmaTL_err(other.sigmaTL_err),
          sigmaTT(other.sigmaTT), sigmaTT_err(other.sigmaTT_err)
    {
        cov.ResizeTo(other.cov.GetNrows(), other.cov.GetNcols());
        cov = other.cov;
    }

    // Custom assignment operator
    FourierFit& operator=(const FourierFit& other) {
        if (this != &other) {
            ok = other.ok;
            absolute_xsec_fit = other.absolute_xsec_fit;
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
            cov.ResizeTo(other.cov.GetNrows(), other.cov.GetNcols());
            cov = other.cov;
        }
        return *this;
    }

    FourierFit() = default;
};

struct SliceResult {
    nps_simc_pi0::Model reference_model;
    std::vector<PhiBin> phi;
    double sumw_data = 0.0, sumw2_data = 0.0;
    double sumw_sim  = 0.0, sumw2_sim  = 0.0;
    double mean_q2_data = 0.0, mean_xb_data = 0.0, mean_tprime_data = 0.0;
    double mean_q2_sim  = 0.0, mean_xb_sim  = 0.0, mean_tprime_sim  = 0.0;
    double mean_t_sim = 0.0;
    double ref_base_sum = 0.0, ref_w_sum = 0.0, ref_xb_sum = 0.0, ref_q2_sum = 0.0;
    double mean_q2_abs  = 0.0, mean_xb_abs  = 0.0, mean_tprime_abs  = 0.0;
    bool has_model_xsec = false, reference_supported = false;
    double q2_direct_base_mean = std::numeric_limits<double>::quiet_NaN();
    double epsilon = 0.0;
    FourierFit fit_ratio;
    FourierFit fit_xsec;
    bool partons_ok = false;
    double partons_epsilon = 0.0, partons_electron_flux_xbq2 = 0.0;
    double partons_sigmaU = 0.0, partons_sigmaLT = 0.0, partons_sigmaTT = 0.0;
};

struct CutFlow {
    long long n_data_total = 0, n_data_pass = 0, n_data_inrange = 0;
    long long n_sim_total = 0, n_sim_pass = 0, n_sim_inrange = 0;
    long long n_sim_bad_full_weight = 0, n_sim_bad_sigcm = 0;
    long long n_sim_bad_vertex = 0, n_sim_bad_model = 0;
    long long n_sim_epsilon_fallback = 0;
};

struct ReweightEvent {
    int islice = -1, iphi = -1;
    double q2 = 0, w = 0, t = 0, phi = 0, eps = 0, base = 0;
    double rec_q2 = 0, rec_xb = 0, rec_t = 0, rec_tprime = 0, rec_phi = 0;
};

struct MissingMassDiagnostic {
    std::unique_ptr<TH1D> data;
    std::unique_ptr<TH1D> exclusive;
    double exclusive_shape_scale = 0.0;
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

static std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out += c;
    }
    out += "'";
    return out;
}

static void apply_publication_style() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);
    gStyle->SetTitleSize(0.050, "XYZ");
    gStyle->SetLabelSize(0.042, "XYZ");
    gStyle->SetTitleOffset(1.10, "X");
    gStyle->SetTitleOffset(1.35, "Y");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    gStyle->SetLegendFont(42);
    gStyle->SetEndErrorSize(4);
}

static void set_pub_pad(double left = 0.13, double right = 0.04, double bottom = 0.13, double top = 0.08) {
    if (!gPad) return;
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetBottomMargin(bottom);
    gPad->SetTopMargin(top);
    gPad->SetTicks(1, 1);
}

static void style_axes(TAxis* x, TAxis* y, double label_size = 0.040, double title_size = 0.046) {
    if (x) {
        x->SetLabelFont(42);
        x->SetTitleFont(42);
        x->SetLabelSize(label_size);
        x->SetTitleSize(title_size);
        x->SetTitleOffset(1.05);
    }
    if (y) {
        y->SetLabelFont(42);
        y->SetTitleFont(42);
        y->SetLabelSize(label_size);
        y->SetTitleSize(title_size);
        y->SetTitleOffset(1.35);
    }
}

static void style_legend(TLegend* leg, double text_size = 0.036) {
    if (!leg) return;
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(text_size);
}

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
    explicit ExclPi0XSecAnalysis(const AnalysisConfig& c) : cfg(c), model_parameters(nps_pi0_reweight::choose(c.model_identifier).default_parameters()) {}
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
    bool has_model_xsec = false, has_vertex_epsilon = false;
    std::string model_xsec_branch;
    std::vector<ReweightEvent> reweight_events;
    std::vector<double> model_parameters, model_errors;
    std::vector<int> free_model_indices;
    TMatrixD model_covariance;
    std::string model_fit_status = "not_run";
    double model_chi2_before = 0, model_chi2_after = 0, model_fit_pvalue = 0;
    int model_fit_ndf = 0, model_fit_bins = 0;
    double default_sigcm_max_relative_difference = 0;
    long long default_sigcm_mismatch_count = 0;

    std::string combined_pdf_path;
    std::vector<std::string> generated_pdf_paths;

    std::vector<double> phi_edges, tprime_edges, t_edges, q2_edges, xb_edges;
    std::vector<std::vector<double>> xb_edges_by_q2;
    std::vector<SliceResult> slices;
    std::vector<MissingMassDiagnostic> mmiss_diagnostics;

    double yield_norm_data = 0.0, yield_norm_data_sumw2 = 0.0;
    double yield_norm_sim = 0.0, yield_norm_sim_sumw2 = 0.0;
    double yield_norm_scale = 1.0, yield_norm_scale_err = 0.0;

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
    void apply_simc_to_data_yield_normalization();
    double model_objective(const std::vector<double>& p, bool strict) const;
    void fit_model();
    void rebuild_simulation(const std::vector<double>& p);
    void compute_mmiss_shape_scales();
    void compute_ratios_and_xsec();
    void fit_slices();
    void compute_partons_projection();
    void make_global_plots();
    void make_mmiss_comparison_plots();
    void make_yield_diagnostic_plots();
    void make_epsilon_plots();
    void make_slice_plots();
    void make_sigma_vs_t_plots();
    void make_partons_projection_plots();
    void write_results();
    void write_csv();
    void write_slice_csv();
    void cleanup();
    void init_combined_pdf();
    void close_combined_pdf();
    std::string format_slice_bin_label(int it, int iq, int ix) const;
    void draw_slice_bin_label(int it, int iq, int ix, double y_ndc = 0.92, double x_ndc = 0.16) const;

    void fill_data_event(double q2, double t, double tmin, double xb, double phi, double mmiss_all, double pi0_weight, float scale, double charge_uC, double total_charge_uC, int helicity, bool use_helicity, double W);
    void fill_sim_event(float q2, float t, float tmin, float xb, float phi, float mmiss,
                        float full_weight, float model_xsec, int is_exclusive, float W,
                        float vq2, float vw, float vt, float vphi, float veps);
    void fill_mmiss_data_diagnostic(double q2, double t, double tmin, double xb,
                                    double phi, double mmiss, double weight);
    void fill_mmiss_sim_diagnostic(float q2, float t, float tmin, float xb,
                                   float phi, float mmiss, float full_weight,
                                   float model_xsec, int is_exclusive);
    void finalize_slice_means(SliceResult& s);

    bool slice_passes_kin(const float q2, const float xb, const double tprime, const double t) const {
        return (q2 >= cfg.q2_min && q2 <= cfg.q2_max &&
                xb >= cfg.xb_min && xb <= cfg.xb_max &&
                tprime >= cfg.tprime_min && tprime <= cfg.tprime_max &&
                t >= cfg.t_min && t <= cfg.t_max &&
                xsec_inside_diamond(cfg, xb, q2));
    }

    bool passes_mmiss_cut(double mmiss) const {
        return std::isfinite(mmiss) && mmiss > cfg.mmiss_lower_gev && mmiss < cfg.mmiss_upper_gev;
    }

    double calc_tprime(float t, float tmin) const { return static_cast<double>(t) - static_cast<double>(tmin); }

    double calc_w(float q2, float xb) const {
        return std::sqrt(std::max(0.0, q2_xb_to_w2(q2, xb, cfg.mp)));
    }

    double calc_sigma_model_slice(const PhiBin& pb) const {
        if (pb.sim <= 0.0 || pb.model_wsum <= 0.0) return 0.0;
        return pb.model_xsec_wsum / pb.model_wsum;
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
    const std::vector<std::string> data_required = {
        "Q2", "t", "tmin", "xB", "phi", "mmiss_all", "pi0_weight",
        "scale", "charge_uC", "run_number", "W"
    };
    const std::vector<std::string> sim_required = {
        "Q2", "t", "tmin", "xB", "phi", "mmiss", "full_weight",
        "sigcm", "is_exclusive", "W", "Q2i", "Wi", "ti", "phipqi"
    };
    for (const auto& name : data_required)
        if (!t_data->GetBranch(name.c_str())) die("Missing mandatory data branch: " + name);
    for (const auto& name : sim_required)
        if (!t_sim->GetBranch(name.c_str())) die("Missing mandatory SIMC branch: " + name);
}

void ExclPi0XSecAnalysis::detect_optional_branches() {
    has_helicity = (t_data->GetBranch("helicity") != nullptr);
    model_xsec_branch = cfg.model_xsec_mode;
    has_model_xsec = (t_sim->GetBranch(model_xsec_branch.c_str()) != nullptr);
    has_vertex_epsilon = t_sim->GetBranch("epsilon_i") != nullptr;
    if (!has_vertex_epsilon)
        warn("SIMC epsilon_i absent: vertex epsilon is reconstructed from Q2i, Wi and fixed --ebeam; exact event-level reweighting cannot be claimed.");

    if (has_helicity)
        warn("LT' is unavailable in the SIMC-model extraction: helicity-specific luminosities and beam polarization are not inputs.");
    if (!has_model_xsec) die("Mandatory SIMC model cross-section branch 'sigcm' not found.");
}

void ExclPi0XSecAnalysis::build_binning() {
    // The two methods read identical fixed Q2, xB, phi and tprime edges.
    // The ratio method uses its separately configured physical-t edges.
    phi_edges = cfg.phi_bin_edges;
    tprime_edges = cfg.tprime_bin_edges;
    t_edges = cfg.t_bin_edges;
    q2_edges = cfg.q2_bin_edges;
    xb_edges_by_q2 = cfg.xb_bin_edges_by_q2;
    // Retain the legacy global xB metadata field. Per-Q2 rows are authoritative.
    xb_edges = xb_edges_by_q2.front();

    auto log_edges = [&](const char* name, const std::vector<double>& edges) {
        std::ostringstream out;
        out << name << " edges:";
        for (double edge : edges) out << " " << edge;
        log(out.str());
    };
    log_edges("phi", phi_edges);
    log_edges("tprime", tprime_edges);
    log_edges("physical t", t_edges);
    log_edges("Q2", q2_edges);
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        std::string name = "xB Q2 bin " + std::to_string(iq);
        log_edges(name.c_str(), xb_edges_by_q2[iq]);
    }
}

void ExclPi0XSecAnalysis::init_storage() {
    slices.assign(cfg.n_tprime * cfg.n_q2 * cfg.n_xb, SliceResult{});
    for (auto& s : slices) s.phi.resize(cfg.n_phi);
    mmiss_diagnostics.resize(slices.size());
    for (size_t i = 0; i < mmiss_diagnostics.size(); ++i) {
        auto make_hist = [&](const std::string& channel) {
            auto h = std::make_unique<TH1D>(
                ("h_mmiss_" + channel + "_slice" + std::to_string(i)).c_str(),
                ";Reconstructed missing mass [GeV];Weighted yield / 25 MeV",
                100, 0.0, 2.5);
            h->SetDirectory(nullptr);
            h->Sumw2();
            return h;
        };
        mmiss_diagnostics[i].data = make_hist("data");
        mmiss_diagnostics[i].exclusive = make_hist("exclusive");
    }

    h_q2_data.reset(new TH1D("h_q2_data", "Data;Q^{2} [GeV^{2}];Weighted counts", 100, cfg.q2_min, cfg.q2_max));
    h_q2_sim.reset(new TH1D("h_q2_sim", "SIMC;Q^{2} [GeV^{2}];Weighted counts", 100, cfg.q2_min, cfg.q2_max));
    h_xb_data.reset(new TH1D("h_xb_data", "Data;x_{B};Weighted counts", 100, cfg.xb_min, cfg.xb_max));
    h_xb_sim.reset(new TH1D("h_xb_sim", "SIMC;x_{B};Weighted counts", 100, cfg.xb_min, cfg.xb_max));
    h_tprime_data.reset(new TH1D("h_tprime_data", "Data;t' [GeV^{2}];Weighted counts", 120, cfg.tprime_min, cfg.tprime_max));
    h_tprime_sim.reset(new TH1D("h_tprime_sim", "SIMC;t' [GeV^{2}];Weighted counts", 120, cfg.tprime_min, cfg.tprime_max));
    h_phi_data.reset(new TH1D("h_phi_data", "Data;#phi [rad];Weighted counts", cfg.n_phi, phi_edges.data()));
    h_phi_sim.reset(new TH1D("h_phi_sim", "SIMC;#phi [rad];Weighted counts", cfg.n_phi, phi_edges.data()));
    h_q2_xb_data.reset(new TH2D("h_q2_xb_data", "Data;Q^{2} [GeV^{2}];x_{B}", 100, cfg.q2_min, cfg.q2_max, 100, cfg.xb_min, cfg.xb_max));
    h_q2_xb_sim.reset(new TH2D("h_q2_xb_sim", "SIMC;Q^{2} [GeV^{2}];x_{B}", 100, cfg.q2_min, cfg.q2_max, 100, cfg.xb_min, cfg.xb_max));
    h_tprime_phi_data.reset(new TH2D("h_tprime_phi_data", "Data;t' [GeV^{2}];#phi [rad]", 100, cfg.tprime_min, cfg.tprime_max, cfg.n_phi, phi_edges.data()));
    h_tprime_phi_sim.reset(new TH2D("h_tprime_phi_sim", "SIMC;t' [GeV^{2}];#phi [rad]", 100, cfg.tprime_min, cfg.tprime_max, cfg.n_phi, phi_edges.data()));
}

void ExclPi0XSecAnalysis::fill_mmiss_data_diagnostic(double q2, double t, double tmin,
        double xb, double phi, double mmiss, double weight) {
    if (!cfg.diagnostics || !std::isfinite(mmiss) || mmiss <= 0.0 || mmiss >= 2.5 ||
        !std::isfinite(weight)) return;
    const double tp = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tp, t) || !std::isfinite(phi)) return;
    const int it = find_bin(t_edges, t, false);
    const int iq = find_bin(q2_edges, q2, false);
    if (it < 0 || iq < 0 || iq >= cfg.n_q2 || find_bin(phi_edges, phi, true) < 0) return;
    const int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], xb, false);
    if (ix < 0) return;
    mmiss_diagnostics[static_cast<size_t>(slice_index(it, iq, ix))].data->Fill(mmiss, weight);
}

void ExclPi0XSecAnalysis::fill_mmiss_sim_diagnostic(float q2, float t, float tmin,
        float xb, float phi, float mmiss, float full_weight, float model_xsec,
        int is_exclusive) {
    if (!cfg.diagnostics || !is_exclusive || !std::isfinite(mmiss) ||
        mmiss <= 0.0f || mmiss >= 2.5f || !std::isfinite(full_weight) ||
        full_weight < 0.0f || !std::isfinite(model_xsec) || model_xsec < 0.0f) return;
    const double tp = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tp, t) || !std::isfinite(phi)) return;
    const int it = find_bin(t_edges, t, false);
    const int iq = find_bin(q2_edges, q2, false);
    if (it < 0 || iq < 0 || iq >= cfg.n_q2 || find_bin(phi_edges, phi, true) < 0) return;
    const int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], xb, false);
    if (ix < 0) return;
    mmiss_diagnostics[static_cast<size_t>(slice_index(it, iq, ix))].exclusive->Fill(
        mmiss, static_cast<double>(full_weight) * cfg.simc_yield_scale);
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

void ExclPi0XSecAnalysis::fill_data_event(double q2, double t, double tmin, double xb, double phi, double mmiss_all, double pi0_weight, float scale, double charge_uC, double total_charge_uC, int helicity, bool use_helicity, double W) {
    (void)W;
    cutflow.n_data_total++;
    if (!passes_mmiss_cut(mmiss_all)) return;
    double tprime = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tprime, t)) return;
    if (!std::isfinite(phi)) return;
    if (!std::isfinite(pi0_weight) || !std::isfinite(scale)) return;

    cutflow.n_data_pass++;

    int it = find_bin(t_edges, t, false);
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
    double w = pi0_weight * static_cast<double>(scale) * charge_fraction;

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

void ExclPi0XSecAnalysis::fill_sim_event(float q2, float t, float tmin, float xb, float phi,
                                         float mmiss, float full_weight, float model_xsec,
                                         int is_exclusive, float W,
                                         float vq2, float vw, float vt, float vphi, float veps) {
    (void)W;
    cutflow.n_sim_total++;
    if (!is_exclusive) return;
    if (!passes_mmiss_cut(mmiss)) return;
    double tprime = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tprime, t)) return;
    if (!std::isfinite(full_weight) || full_weight < 0) {
        ++cutflow.n_sim_bad_full_weight; return;
    }
    if (!std::isfinite(model_xsec) || model_xsec <= 0) {
        ++cutflow.n_sim_bad_sigcm; return;
    }
    if (!(std::isfinite(phi) && std::isfinite(vq2) && std::isfinite(vw) &&
          std::isfinite(vt) && std::isfinite(vphi))) {
        ++cutflow.n_sim_bad_vertex; return;
    }

    cutflow.n_sim_pass++;

    int it = find_bin(t_edges, t, false);
    int iq = find_bin(q2_edges, q2, false);
    if (iq < 0 || iq >= cfg.n_q2) return;
    int ix = find_bin(xb_edges_by_q2[static_cast<size_t>(iq)], xb, false);
    int ip = find_bin(phi_edges, phi, true);
    if (it < 0 || ix < 0 || ip < 0) return;
    cutflow.n_sim_inrange++;

    const double base = static_cast<double>(full_weight) / model_xsec;
    if (!std::isfinite(base)) { ++cutflow.n_sim_bad_full_weight; return; }
    double event_eps = veps;
    if (!has_vertex_epsilon ||
        !(std::isfinite(event_eps) && event_eps >= 0 && event_eps <= 1)) {
        try { event_eps = nps_pi0_reweight::epsilon(vq2, vw, cfg.ebeam, cfg.mp); }
        catch (const std::exception&) { ++cutflow.n_sim_bad_vertex; return; }
        ++cutflow.n_sim_epsilon_fallback;
    }
    if (!(std::isfinite(event_eps) && event_eps >= 0 && event_eps <= 1)) {
        ++cutflow.n_sim_bad_vertex; return;
    }
    double original_model = 0;
    try {
        original_model = nps_pi0_reweight::choose(cfg.model_identifier).cross_section(
            vq2, vw, -vt, vphi, event_eps, cfg.mp, cfg.mpi0, model_parameters);
    } catch (const std::exception&) { ++cutflow.n_sim_bad_model; return; }
    if (!(original_model > 0)) { ++cutflow.n_sim_bad_model; return; }
    const double rel = std::abs(original_model/model_xsec-1);
    default_sigcm_max_relative_difference = std::max(default_sigcm_max_relative_difference, rel);
    if (rel > 0.01) ++default_sigcm_mismatch_count;
    ReweightEvent event;
    event.islice = slice_index(it,iq,ix); event.iphi=ip;
    event.q2=vq2; event.w=vw; event.t=-vt; event.phi=vphi;
    event.eps=event_eps; event.base=base;
    event.rec_q2=q2; event.rec_xb=xb; event.rec_t=t; event.rec_tprime=tprime; event.rec_phi=wrap_phi(phi);
    reweight_events.push_back(event);
    double w = static_cast<double>(full_weight) * cfg.simc_yield_scale;
    PhiBin& pb = slice(it, iq, ix).phi[ip];
    pb.sim += w;
    pb.sim_sumw2 += w * w;
    pb.model_wsum += w;
    // Retain the old event-weighted model mean as a diagnostic only.
    pb.model_xsec_wsum += w * static_cast<double>(model_xsec);
    pb.n_sim += 1;
    pb.mean_q2_sim += q2 * w;
    pb.mean_xb_sim += xb * w;
    pb.mean_tprime_sim += tprime * w;

    SliceResult& s = slice(it, iq, ix);
    s.sumw_sim += w;
    s.sumw2_sim += w * w;
    s.mean_q2_sim += q2 * w;
    s.mean_xb_sim += xb * w;
    s.mean_t_sim += t * w;
    s.mean_tprime_sim += tprime * w;
    s.ref_base_sum += base;
    s.ref_w_sum += base * vw;
    s.ref_xb_sum += base * (vq2 / (vw*vw - cfg.mp*cfg.mp + vq2));
    s.ref_q2_sum += base * vq2;
    accumulate_global_histograms(q2, xb, tprime, wrap_phi(phi), 0.0, w);
}

void ExclPi0XSecAnalysis::fill_from_trees() {
    // Data loop
    double dq2 = 0, dt = 0, dtmin = 0, dxb = 0, dphi = 0, dmmiss_all = 0, dpi0_weight = 0, dW = 0;
    float dcharge_uC = 0;
    float dscale = 0;
    int dhelicity = 0;
    int drun_number = 0;
    t_data->SetBranchAddress("Q2", &dq2);
    t_data->SetBranchAddress("t", &dt);
    t_data->SetBranchAddress("tmin", &dtmin);
    t_data->SetBranchAddress("xB", &dxb);
    t_data->SetBranchAddress("phi", &dphi);
    t_data->SetBranchAddress("mmiss_all", &dmmiss_all);
    t_data->SetBranchAddress("pi0_weight", &dpi0_weight);
    t_data->SetBranchAddress("scale", &dscale);
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
        double charge_fraction = 1.0;
        if (std::isfinite(dcharge_uC) && dcharge_uC > 0.0 && total_charge_uC > 0.0)
            charge_fraction = dcharge_uC / total_charge_uC;
        fill_mmiss_data_diagnostic(dq2, dt, dtmin, dxb, dphi, dmmiss_all,
            dpi0_weight * static_cast<double>(dscale) * charge_fraction);
        fill_data_event(dq2, dt, dtmin, dxb, dphi, dmmiss_all, dpi0_weight, dscale, dcharge_uC, total_charge_uC, dhelicity, has_helicity, dW);
    }

    // SIMC loop
    float sim_q2 = 0, sim_t = 0, sim_tmin = 0, sim_xb = 0, sim_phi = 0, sim_mmiss = 0, full_weight = 0, model_xsec = 0, sim_W = 0;
    float vq2=0, vw=0, vt=0, vphi=0, veps=0;
    int sim_is_exclusive = 0;
    t_sim->SetBranchAddress("Q2", &sim_q2);
    t_sim->SetBranchAddress("t", &sim_t);
    t_sim->SetBranchAddress("tmin", &sim_tmin);
    t_sim->SetBranchAddress("xB", &sim_xb);
    t_sim->SetBranchAddress("phi", &sim_phi);
    t_sim->SetBranchAddress("mmiss", &sim_mmiss);
    t_sim->SetBranchAddress("full_weight", &full_weight);
    t_sim->SetBranchAddress("is_exclusive", &sim_is_exclusive);
    t_sim->SetBranchAddress("W", &sim_W);
    t_sim->SetBranchAddress("Q2i", &vq2);
    t_sim->SetBranchAddress("Wi", &vw);
    t_sim->SetBranchAddress("ti", &vt);
    t_sim->SetBranchAddress("phipqi", &vphi);
    if (has_vertex_epsilon) t_sim->SetBranchAddress("epsilon_i", &veps);
    if (has_model_xsec) t_sim->SetBranchAddress(model_xsec_branch.c_str(), &model_xsec);

    const Long64_t nsim = t_sim->GetEntries();
    for (Long64_t i = 0; i < nsim; ++i) {
        t_sim->GetEntry(i);
        fill_mmiss_sim_diagnostic(sim_q2, sim_t, sim_tmin, sim_xb, sim_phi,
                                  sim_mmiss, full_weight, model_xsec, sim_is_exclusive);
        fill_sim_event(sim_q2, sim_t, sim_tmin, sim_xb, sim_phi, sim_mmiss,
                       full_weight, model_xsec, sim_is_exclusive, sim_W,
                       vq2, vw, vt, vphi, veps);
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
            s.mean_t_sim /= s.sumw_sim;
            s.mean_tprime_sim /= s.sumw_sim;
        }
        for (auto& pb : s.phi) {
            if (pb.data > 0.0) {
                pb.mean_q2_data /= pb.data;
                pb.mean_xb_data /= pb.data;
                pb.mean_tprime_data /= pb.data;
            }
            if (pb.sim > 0.0) {
                pb.mean_q2_sim /= pb.sim;
                pb.mean_xb_sim /= pb.sim;
                pb.mean_tprime_sim /= pb.sim;
            }
        }
    }
}

void ExclPi0XSecAnalysis::apply_simc_to_data_yield_normalization() {
    if (!cfg.normalize_mmiss) return;

    // Use exactly the yields entering the extraction: all selected events in
    // the configured (Q2,xB,t',phi,Mmiss) phase space. Data are not divided by
    // tgt_contam in this mode; parse_config forces that correction to unity.
    for (const auto& s : slices) {
        yield_norm_data += s.sumw_data;
        yield_norm_data_sumw2 += s.sumw2_data;
        yield_norm_sim += s.sumw_sim;
        yield_norm_sim_sumw2 += s.sumw2_sim;
    }
    if (!(std::isfinite(yield_norm_data) && yield_norm_data > 0.0 &&
          std::isfinite(yield_norm_sim) && yield_norm_sim > 0.0))
        die("--normalize_mmiss: nonpositive data or exclusive-SIMC yield in the selected phase space.");

    yield_norm_scale = yield_norm_data / yield_norm_sim;
    if (!(std::isfinite(yield_norm_scale) && yield_norm_scale > 0.0))
        die("--normalize_mmiss: invalid data/SIMC yield scale.");
    yield_norm_scale_err = yield_norm_scale * std::sqrt(
        yield_norm_data_sumw2 / (yield_norm_data * yield_norm_data) +
        yield_norm_sim_sumw2 / (yield_norm_sim * yield_norm_sim));

    const double k = yield_norm_scale;
    const double k2 = k * k;
    for (auto& s : slices) {
        // Kinematic means were finalized from the default or final event weights and must not move.
        s.sumw_sim *= k;
        s.sumw2_sim *= k2;
        for (auto& pb : s.phi) {
            pb.sim *= k;
            pb.sim_sumw2 *= k2;
            // Preserve <sigcm> = sum(w*sigcm)/sum(w).
            pb.model_wsum *= k;
            pb.model_xsec_wsum *= k;
        }
    }
    for (TH1D* h : {h_q2_sim.get(), h_xb_sim.get(), h_tprime_sim.get(), h_phi_sim.get()})
        if (h) h->Scale(k);
    for (TH2D* h : {h_q2_xb_sim.get(), h_tprime_phi_sim.get()})
        if (h) h->Scale(k);

    std::cout << std::setprecision(10)
              << "SIMC-to-data yield normalization [" << cfg.mmiss_lower_gev << ", " << cfg.mmiss_upper_gev
              << "] GeV: data=" << yield_norm_data << " +/- " << std::sqrt(yield_norm_data_sumw2)
              << ", exclusive SIMC(reweighted*simc_yield_scale)=" << yield_norm_sim << " +/- "
              << std::sqrt(yield_norm_sim_sumw2) << ", scale=" << yield_norm_scale
              << " +/- " << yield_norm_scale_err
              << " (independent-sum diagnostic; extraction errors are conditional on this scale)\n";
}

void ExclPi0XSecAnalysis::compute_mmiss_shape_scales() {
    if (!cfg.diagnostics) return;
    // Same display normalization as the no-SIMC-model extractor:
    //   scale = Integral(data) / Integral(exclusive SIMC)
    // using the configured missing-mass window. This scale is applied only
    // to a plot clone in make_mmiss_comparison_plots(); it never modifies
    // the event-level reweighted sums used by the cross-section extraction.
    for (auto& m : mmiss_diagnostics) {
        double data_integral = 0.0;
        double sim_integral = 0.0;
        for (int bin = 1; bin <= m.data->GetNbinsX(); ++bin) {
            const double center = m.data->GetBinCenter(bin);
            if (center >= cfg.mmiss_lower_gev && center <= cfg.mmiss_upper_gev) {
                data_integral += m.data->GetBinContent(bin);
                sim_integral += m.exclusive->GetBinContent(bin);
            }
        }
        if (data_integral > 0.0 && sim_integral > 0.0)
            m.exclusive_shape_scale = data_integral / sim_integral;
    }
}

double ExclPi0XSecAnalysis::model_objective(const std::vector<double>& p, bool strict) const {
    std::vector<double> sums(slices.size()*cfg.n_phi, 0.0);
    std::vector<double> sumw2(sums.size(), 0.0);
    for (const auto& e : reweight_events) {
        double sigma = 0;
        try {
            sigma = nps_pi0_reweight::choose(cfg.model_identifier).cross_section(
                e.q2,e.w,e.t,e.phi,e.eps,cfg.mp,cfg.mpi0,p);
        } catch (const std::exception&) {
            if (strict) throw;
            return 1.e30;
        }
        if (!(std::isfinite(sigma) && sigma > 0)) {
            if (strict) throw std::runtime_error("Nonpositive model prediction for selected SIMC event");
            return 1.e30;
        }
        const double weight = cfg.simc_yield_scale*e.base*sigma;
        if (!std::isfinite(weight)) {
            if (strict) throw std::runtime_error("Nonfinite reweighted SIMC event");
            return 1.e30;
        }
        const size_t index = size_t(e.islice)*cfg.n_phi+e.iphi;
        sums[index] += weight;
        sumw2[index] += weight*weight;
    }
    double chi2 = 0;
    for (size_t is=0;is<slices.size();++is)
        for (int ip=0;ip<cfg.n_phi;++ip) {
            const auto& pb=slices[is].phi[ip];
            const size_t index=is*cfg.n_phi+ip;
            if (pb.n_data==0 || pb.n_sim==0) continue;
            const double d=pb.data/cfg.tgt_contam;
            const double var=pb.data_sumw2/(cfg.tgt_contam*cfg.tgt_contam)+sumw2[index];
            if (!(std::isfinite(var) && var>0 && std::isfinite(sums[index]))) {
                if (strict) throw std::runtime_error("Invalid yield-space fit variance");
                return 1.e30;
            }
            chi2 += (d-sums[index])*(d-sums[index])/var;
        }
    return std::isfinite(chi2) ? chi2 : 1.e30;
}

void ExclPi0XSecAnalysis::rebuild_simulation(const std::vector<double>& p) {
    for (auto& slice : slices) {
        slice.sumw_sim=0; slice.sumw2_sim=0;
        slice.mean_q2_sim=0; slice.mean_xb_sim=0;
        slice.mean_t_sim=0; slice.mean_tprime_sim=0;
        for (auto& pb : slice.phi) {
            pb.sim=0; pb.sim_sumw2=0;
            pb.mean_q2_sim=0; pb.mean_xb_sim=0; pb.mean_tprime_sim=0;
        }
    }
    for (TH1D* h : {h_q2_sim.get(),h_xb_sim.get(),h_tprime_sim.get(),h_phi_sim.get()})
        if (h) h->Reset();
    for (TH2D* h : {h_q2_xb_sim.get(),h_tprime_phi_sim.get()})
        if (h) h->Reset();
    for (const auto& e : reweight_events) {
        const double sigma=nps_pi0_reweight::choose(cfg.model_identifier).cross_section(
            e.q2,e.w,e.t,e.phi,e.eps,cfg.mp,cfg.mpi0,p);
        if (!(std::isfinite(sigma) && sigma>0))
            die("Nonpositive final model prediction for selected SIMC event");
        const double w=cfg.simc_yield_scale*e.base*sigma;
        if (!std::isfinite(w)) die("Nonfinite final SIMC weight");
        auto& slice=slices[e.islice];
        auto& pb=slice.phi[e.iphi];
        pb.sim+=w; pb.sim_sumw2+=w*w;
        pb.mean_q2_sim+=w*e.rec_q2;
        pb.mean_xb_sim+=w*e.rec_xb;
        pb.mean_tprime_sim+=w*e.rec_tprime;
        slice.sumw_sim+=w; slice.sumw2_sim+=w*w;
        slice.mean_q2_sim+=w*e.rec_q2;
        slice.mean_xb_sim+=w*e.rec_xb;
        slice.mean_t_sim+=w*e.rec_t;
        slice.mean_tprime_sim+=w*e.rec_tprime;
        accumulate_global_histograms(e.rec_q2,e.rec_xb,e.rec_tprime,e.rec_phi,0,w);
    }
    for (auto& slice : slices) {
        if (slice.sumw_sim>0) {
            slice.mean_q2_sim/=slice.sumw_sim;
            slice.mean_xb_sim/=slice.sumw_sim;
            slice.mean_t_sim/=slice.sumw_sim;
            slice.mean_tprime_sim/=slice.sumw_sim;
        }
        for (auto& pb : slice.phi)
            if (pb.sim>0) {
                pb.mean_q2_sim/=pb.sim;
                pb.mean_xb_sim/=pb.sim;
                pb.mean_tprime_sim/=pb.sim;
            }
    }
}

void ExclPi0XSecAnalysis::fit_model() {
    if (reweight_events.empty()) die("No usable SIMC events for model reweighting");
    model_chi2_before=model_objective(model_parameters,true);
    rebuild_simulation(model_parameters);
    for (auto& slice : slices)
        for (auto& pb : slice.phi)
            if (pb.n_data && pb.sim>0)
                pb.ratio_before=(pb.data/cfg.tgt_contam)/pb.sim;
    const auto spec=nps_pi0_reweight::choose(cfg.model_identifier).parameter_spec();
    if (!cfg.fixed_default_model) {
        std::stringstream names(cfg.model_free_parameters);
        std::string token;
        while (std::getline(names,token,',')) {
            token.erase(0,token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t")+1);
            if (token.empty()) continue;
            auto it=std::find_if(spec.begin(),spec.end(),
                                 [&](const auto& x){return x.name==token;});
            if (it==spec.end()) die("Unknown model parameter: "+token);
            const int index=int(it-spec.begin());
            const int fortran_coefficient=index%17+1;
            if (fortran_coefficient==1 || fortran_coefficient==2 ||
                fortran_coefficient==3 || fortran_coefficient==4 ||
                fortran_coefficient==11 || fortran_coefficient==15 ||
                fortran_coefficient==17)
                die("Parameter "+token+" is inactive when pi0 fpifact=0");
            if (std::find(free_model_indices.begin(),free_model_indices.end(),index)!=free_model_indices.end())
                die("Repeated model parameter: "+token);
            free_model_indices.push_back(index);
        }
        if (free_model_indices.empty()) die("No free model parameters; use --fixed-default-model");
    }
    for (const auto& slice : slices)
        for (const auto& pb : slice.phi)
            if (pb.n_data && pb.n_sim) ++model_fit_bins;
    model_fit_ndf=model_fit_bins-int(free_model_indices.size());
    if (model_fit_ndf<=0 && !cfg.fixed_default_model)
        die("Insufficient supported bins for selected model parameters");
    model_errors.assign(spec.size(),std::numeric_limits<double>::quiet_NaN());
    if (cfg.fixed_default_model) {
        model_fit_status="fixed_default";
    } else {
        auto minimizer=std::unique_ptr<ROOT::Math::Minimizer>(
            ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
        if (!minimizer) die("Cannot create ROOT Minuit2 minimizer");
        minimizer->SetMaxIterations(cfg.model_max_iterations);
        minimizer->SetMaxFunctionCalls(cfg.model_max_evaluations);
        minimizer->SetTolerance(cfg.model_tolerance);
        minimizer->SetPrintLevel(cfg.verbose ? 0 : -1);
        const unsigned n=free_model_indices.size();
        auto objective=[&](const double* x) {
            auto trial=model_parameters;
            for (unsigned j=0;j<n;++j) trial[free_model_indices[j]]=x[j];
            return model_objective(trial,false);
        };
        ROOT::Math::Functor fn(objective,n);
        minimizer->SetFunction(fn);
        for (unsigned j=0;j<n;++j) {
            const auto& v=spec[free_model_indices[j]];
            minimizer->SetLimitedVariable(j,v.name.c_str(),v.value,
                std::max(1.e-5,0.01*std::abs(v.value)),v.lower,v.upper);
        }
        const bool ok=minimizer->Minimize();
        if (!ok || minimizer->Status()!=0) {
            std::ostringstream failure;
            failure << "SIMC model fit failed: Minuit2 status=" << minimizer->Status()
                    << " EDM=" << minimizer->Edm() << " objective=" << minimizer->MinValue()
                    << " calls=" << minimizer->NCalls();
            for (unsigned j=0;j<n;++j) failure << " " << spec[free_model_indices[j]].name
                                                << "=" << minimizer->X()[j];
            die(failure.str());
        }
        for (unsigned j=0;j<n;++j) {
            model_parameters[free_model_indices[j]]=minimizer->X()[j];
            model_errors[free_model_indices[j]]=minimizer->Errors()[j];
        }
        model_covariance.ResizeTo(n,n);
        for (unsigned i=0;i<n;++i)
            for (unsigned j=0;j<n;++j)
                model_covariance(i,j)=minimizer->CovMatrix(i,j);
        model_fit_status="converged";
    }
    model_chi2_after=model_objective(model_parameters,true);
    model_fit_pvalue=model_fit_ndf>0 ? TMath::Prob(model_chi2_after,model_fit_ndf)
                                       : std::numeric_limits<double>::quiet_NaN();
    if (model_fit_status=="converged" && std::isfinite(model_fit_pvalue) &&
        model_fit_pvalue < 0.01) {
        model_fit_status="converged_residual_mismatch";
        warn("Model optimizer converged, but yield residuals exceed statistical precision (chi2/ndf="+
             std::to_string(model_chi2_after)+"/"+std::to_string(model_fit_ndf)+").");
    }
    rebuild_simulation(model_parameters);
    if (has_vertex_epsilon && cutflow.n_sim_epsilon_fallback)
        warn("epsilon_i contains invalid or missing values for "+
             std::to_string(cutflow.n_sim_epsilon_fallback)+
             " selected events; those use fixed-beam epsilon approximation.");
    std::cout << "[MODEL_REJECT] bad_full_weight=" << cutflow.n_sim_bad_full_weight
              << " bad_sigcm=" << cutflow.n_sim_bad_sigcm
              << " bad_vertex=" << cutflow.n_sim_bad_vertex
              << " epsilon_fallback=" << cutflow.n_sim_epsilon_fallback
              << " unsupported_or_nonphysical_model=" << cutflow.n_sim_bad_model << "\n";
    std::cout << "[MODEL_FIT] status=" << model_fit_status
              << " chi2_before=" << model_chi2_before
              << " chi2_after=" << model_chi2_after
              << " ndf=" << model_fit_ndf << " pvalue=" << model_fit_pvalue
              << " bins=" << model_fit_bins
              << " cached_events=" << reweight_events.size()
              << " default_sigcm_mismatch_gt_1pct=" << default_sigcm_mismatch_count
              << " default_sigcm_max_rel=" << default_sigcm_max_relative_difference << "\n";
    for (int index : free_model_indices)
        std::cout << "[MODEL_PARAM] " << spec[index].name << "="
                  << model_parameters[index] << " +/- " << model_errors[index] << "\n";
}

void ExclPi0XSecAnalysis::compute_ratios_and_xsec() {
    // Propagate has_model_xsec to all slices before loop
    for (auto& s : slices) s.has_model_xsec = has_model_xsec;
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                SliceResult& s = slice(it, iq, ix);
                // Accepted SIMC base-weight means are fixed during fitting and
                // shared by every phi bin. Q2 is derived from independent W,xB.
                if (s.ref_base_sum > 0) {
                    const double wref=s.ref_w_sum/s.ref_base_sum;
                    const double xbref=s.ref_xb_sum/s.ref_base_sum;
                    const double q2ref=xbref*(wref*wref-cfg.mp*cfg.mp)/(1-xbref);
                    const double tcenter=0.5*(t_edges[it]+t_edges[it+1]);
                    s.q2_direct_base_mean=s.ref_q2_sum/s.ref_base_sum;
                    try {
                        const double eps=nps_pi0_reweight::epsilon(q2ref,wref,cfg.ebeam,cfg.mp);
                        auto eval=[&](double phi) {
                            return nps_pi0_reweight::choose(cfg.model_identifier).cross_section(q2ref,wref,tcenter,phi,eps,
                                cfg.mp,cfg.mpi0,model_parameters);
                        };
                        const double f0=eval(0), f90=eval(0.5*TMath::Pi()), f180=eval(TMath::Pi());
                        const double even=0.5*(f0+f180);
                        s.reference_model.angular[0]=0.5*(even+f90);
                        s.reference_model.angular[1]=0.5*(f0-f180);
                        s.reference_model.angular[2]=0.5*(even-f90);
                        s.reference_model.q2=q2ref;
                        s.reference_model.xb=xbref;
                        s.reference_model.w=wref;
                        s.reference_model.t=tcenter;
                        s.reference_model.epsilon=eps;
                        s.epsilon=eps;
                        const double wsq=wref*wref;
                        const double eg=(wsq-cfg.mp*cfg.mp-q2ref)/(2*wref);
                        const double pg=std::sqrt(eg*eg+q2ref);
                        const double epi=(wsq+cfg.mpi0*cfg.mpi0-cfg.mp*cfg.mp)/(2*wref);
                        const double pp=std::sqrt(epi*epi-cfg.mpi0*cfg.mpi0);
                        const double tmin=-q2ref+cfg.mpi0*cfg.mpi0-2*eg*epi+2*pg*pp;
                        s.reference_model.tprime=tcenter-tmin;
                        s.reference_supported=true;
                    } catch (const std::exception& e) {
                        warn("Unsupported model reporting point for slice "+
                            std::to_string(it)+"/"+std::to_string(iq)+"/"+
                            std::to_string(ix)+": "+e.what());
                    }
                }
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    PhiBin& pb = s.phi[ip];
                    if (s.reference_supported) {
                        const double phi_center=0.5*(phi_edges[ip]+phi_edges[ip+1]);
                        pb.model_phi_center=nps_pi0_reweight::choose(cfg.model_identifier).cross_section(
                            s.reference_model.q2,s.reference_model.w,s.reference_model.t,
                            phi_center,s.epsilon,cfg.mp,cfg.mpi0,model_parameters);
                        if (!(std::isfinite(pb.model_phi_center) && pb.model_phi_center>0))
                            pb.model_phi_center=std::numeric_limits<double>::quiet_NaN();
                    }
                    pb.data /= cfg.tgt_contam;
                    pb.data_sumw2 /= (cfg.tgt_contam * cfg.tgt_contam);
                    pb.data_plus /= cfg.tgt_contam;
                    pb.data_minus /= cfg.tgt_contam;
                    pb.data_plus_sumw2 /= (cfg.tgt_contam * cfg.tgt_contam);
                    pb.data_minus_sumw2 /= (cfg.tgt_contam * cfg.tgt_contam);

                    if (pb.n_data > 0 && pb.sim > 0.0) {
                        pb.ratio = pb.data / pb.sim;
                        pb.ratio_err = std::sqrt((pb.data_sumw2 / (pb.sim * pb.sim)) +
                                                 (pb.data * pb.data * pb.sim_sumw2 / (pb.sim * pb.sim * pb.sim * pb.sim)));
                    }

                    // Only compute pb.xsec if has_model_xsec is true
                    if (s.reference_supported && pb.n_data > 0 && pb.sim > 0.0 && pb.model_phi_center > 0.0) {
                        double sigma_model = pb.model_phi_center;
                        pb.xsec = pb.ratio * sigma_model;
                        pb.xsec_err = pb.ratio_err * sigma_model;
                        pb.mean_q2_xsec = s.reference_model.q2;
                        pb.mean_xb_xsec = s.reference_model.xb;
                        pb.mean_tprime_xsec = s.reference_model.tprime;
                    } else {
                        pb.ratio = std::numeric_limits<double>::quiet_NaN();
                        pb.ratio_err = std::numeric_limits<double>::quiet_NaN();
                        pb.xsec = std::numeric_limits<double>::quiet_NaN();
                        pb.xsec_err = std::numeric_limits<double>::quiet_NaN();
                    }
                    pb.xsec_sys_tgt = pb.xsec * (cfg.tgt_contam_err / cfg.tgt_contam);
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
                std::vector<std::vector<double>> Xxsec, Xratio;
                std::vector<double> yxsec, eyxsec, yratio, eyratio;
                Xxsec.reserve(cfg.n_phi);
                Xratio.reserve(cfg.n_phi);
                yxsec.reserve(cfg.n_phi);
                eyxsec.reserve(cfg.n_phi);
                yratio.reserve(cfg.n_phi);
                eyratio.reserve(cfg.n_phi);

                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    double phi1 = phi_edges[ip];
                    double phi2 = phi_edges[ip + 1];
                    std::vector<double> basis = phi_basis_means(phi1, phi2);
                    std::vector<double> basis3 = {basis[0], basis[1], basis[2]};
                    std::vector<double> basis4 = basis;
                    PhiBin& pb = s.phi[ip];
                    if (pb.ratio_err > 0.0 && std::isfinite(pb.ratio)) {
                        Xratio.push_back(basis3);
                        yratio.push_back(pb.ratio);
                        eyratio.push_back(pb.ratio_err);
                    }
                    if (s.has_model_xsec && pb.xsec_err > 0.0 && std::isfinite(pb.xsec)) {
                        Xxsec.push_back(basis3);
                        yxsec.push_back(pb.xsec);
                        eyxsec.push_back(pb.xsec_err);
                    }
                    (void)basis4;
                }

                s.fit_ratio = FourierFit();
                if (Xratio.size() >= 3) {
                    s.fit_ratio = weighted_linear_fit(Xratio, yratio, eyratio);
                    if (s.fit_ratio.ok) {
                        s.fit_ratio.absolute_xsec_fit = false;
                        s.fit_ratio.p.resize(3);
                        s.fit_ratio.perr.resize(3);
                    }
                }

                if (s.has_model_xsec && Xxsec.size() >= 3) {
                    s.fit_xsec = weighted_linear_fit(Xxsec, yxsec, eyxsec);
                    if (s.fit_xsec.ok) {
                        s.fit_xsec.absolute_xsec_fit = true;
                        double eps = std::max(1e-12, s.epsilon);
                        // sigcm is d2sigma/(dt dphi). Convert angular-fit
                        // coefficients to U/LT/TT in
                        // [U + k_LT LT cos(phi) + eps TT cos(2phi)]/(2pi).
                        const double two_pi = 2.0 * TMath::Pi();
                        s.fit_xsec.sigmaU = two_pi * s.fit_xsec.p[0];
                        s.fit_xsec.sigmaU_err = two_pi * s.fit_xsec.perr[0];
                        s.fit_xsec.sigmaTL  = two_pi * s.fit_xsec.p[1] / std::sqrt(std::max(1e-12, 2.0 * eps * (1.0 + eps)));
                        s.fit_xsec.sigmaTL_err = two_pi * s.fit_xsec.perr[1] / std::sqrt(std::max(1e-12, 2.0 * eps * (1.0 + eps)));
                        s.fit_xsec.sigmaTT  = two_pi * s.fit_xsec.p[2] / std::max(1e-12, eps);
                        s.fit_xsec.sigmaTT_err = two_pi * s.fit_xsec.perr[2] / std::max(1e-12, eps);
                    }
                }

                // LT' is intentionally not fitted in this method.  The input
                // contract lacks helicity-specific luminosities and beam
                // polarization, so a helicity-difference amplitude cannot be
                // converted into a physical response function.
            }
        }
    }
}

void ExclPi0XSecAnalysis::compute_partons_projection() {
    if (!cfg.partons_projection) return;
#ifdef NPS_ENABLE_PARTONS
    nps_partons_pi0::Model model;
    model.initialize(cfg.partons_executable_path, cfg.partons_warmups, cfg.partons_calls);
    int projected = 0;
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                SliceResult& s = slice(it, iq, ix);
                if (!s.fit_xsec.ok || !(s.sumw_sim > 0.0)) continue;
                try {
                    const auto p = model.predict(
                        s.reference_model.q2, s.reference_model.xb, s.reference_model.t, cfg.ebeam);
                    if (!p.valid) {
                        warn("PARTONS prediction unavailable for slice " +
                             std::to_string(it) + "/" + std::to_string(iq) + "/" +
                             std::to_string(ix));
                        continue;
                    }
                    s.partons_ok = true;
                    s.partons_epsilon = p.epsilon;
                    s.partons_electron_flux_xbq2 = p.electron_flux_xbq2;
                    s.partons_sigmaU = p.sigmaU;
                    s.partons_sigmaLT = p.sigmaLT;
                    s.partons_sigmaTT = p.sigmaTT;
                    ++projected;
                } catch (const std::exception& e) {
                    warn("PARTONS failed for slice " + std::to_string(it) + "/" +
                         std::to_string(iq) + "/" + std::to_string(ix) + ": " + e.what());
                }
            }
        }
    }
    if (projected == 0) die("PARTONS produced no valid pi0 slice predictions.");
    log("PARTONS GK06/GPDGK19 projected " + std::to_string(projected) +
        " common extraction reference points.");
#else
    die("--partons requires compilation with NPS_ENABLE_PARTONS and native PARTONS libraries.");
#endif
}

void ExclPi0XSecAnalysis::write_canvas_pdf_png(TCanvas* c, const std::string& base) {
    if (cfg.write_pdf) {
        const std::string pdf_path = base + ".pdf";
        c->SaveAs(pdf_path.c_str());
        generated_pdf_paths.push_back(pdf_path);
    }
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
        << " t[" << t_edges[it] << "," << t_edges[it + 1] << "]";
    return oss.str();
}

void ExclPi0XSecAnalysis::draw_slice_bin_label(int it, int iq, int ix, double y_ndc, double x_ndc) const {
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.026);
    lat.DrawLatex(x_ndc, y_ndc, format_slice_bin_label(it, iq, ix).c_str());
}

void ExclPi0XSecAnalysis::init_combined_pdf() {
    fs::path p(cfg.out_all_plots_pdf);
    if (p.is_relative()) p = fs::path(cfg.out_dir) / p;
    combined_pdf_path = p.string();
    generated_pdf_paths.clear();
}

void ExclPi0XSecAnalysis::close_combined_pdf() {
    if (!cfg.write_pdf || combined_pdf_path.empty() || generated_pdf_paths.empty()) return;

    fs::create_directories(fs::path(combined_pdf_path).parent_path());

    if (generated_pdf_paths.size() == 1) {
        fs::copy_file(generated_pdf_paths.front(), combined_pdf_path, fs::copy_options::overwrite_existing);
        log("Wrote combined plot PDF: " + combined_pdf_path);
        generated_pdf_paths.clear();
        return;
    }

    if (gSystem->Exec("command -v pdfunite >/dev/null 2>&1") != 0) {
        warn("pdfunite not found; individual plot PDFs were written, but combined PDF was not created.");
        generated_pdf_paths.clear();
        return;
    }

    std::ostringstream cmd;
    cmd << "pdfunite";
    for (const std::string& p : generated_pdf_paths) {
        cmd << " " << shell_quote(p);
    }
    cmd << " " << shell_quote(combined_pdf_path);

    const int status = gSystem->Exec(cmd.str().c_str());
    if (status == 0) {
        log("Wrote combined plot PDF: " + combined_pdf_path);
    } else {
        warn("pdfunite failed; individual plot PDFs were written, but combined PDF was not created.");
    }
    generated_pdf_paths.clear();
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

    // Match the no-SIMC-model diagnostic convention explicitly:
    //   SIMC plot scale = Integral(target-corrected data) / Integral(SIMC).
    // Scale detached clones only. The stored histograms and every per-bin
    // reweighted sum used by sigma_data=(Y_data/Y_SIMC)*sigma_model remains
    // at their physical analysis normalization.
    auto clone_area_normalized_sim_hist = [](const TH1D* src, const char* name,
                                             const TH1D* data_reference) {
        TH1D* h = dynamic_cast<TH1D*>(src->Clone(name));
        if (h) {
            h->SetDirectory(nullptr);
            const double sim_integral = h->Integral();
            const double data_integral = data_reference ? data_reference->Integral() : 0.0;
            if (sim_integral > 0.0 && data_integral > 0.0)
                h->Scale(data_integral / sim_integral);
        }
        return std::unique_ptr<TH1D>(h);
    };

    auto h_q2_sim_shape = clone_area_normalized_sim_hist(
        h_q2_sim.get(), "h_q2_sim_area_normalized", h_q2_data_corr.get());
    auto h_xb_sim_shape = clone_area_normalized_sim_hist(
        h_xb_sim.get(), "h_xb_sim_area_normalized", h_xb_data_corr.get());
    auto h_tprime_sim_shape = clone_area_normalized_sim_hist(
        h_tprime_sim.get(), "h_tprime_sim_area_normalized", h_tprime_data_corr.get());
    auto h_phi_sim_shape = clone_area_normalized_sim_hist(
        h_phi_sim.get(), "h_phi_sim_area_normalized", h_phi_data_corr.get());

    TCanvas c1("c1", "global", 1500, 1100);
    c1.Divide(2,2);

    c1.cd(1);
    set_pub_pad(0.13, 0.04, 0.13, 0.09);
    h_q2_data_corr->SetTitle(";Q^{2} [GeV^{2}];Weighted counts (SIMC area norm.)");
    style_axes(h_q2_data_corr->GetXaxis(), h_q2_data_corr->GetYaxis());
    h_q2_data_corr->SetLineWidth(2); h_q2_data_corr->SetLineColor(kBlack);
    h_q2_sim_shape->SetLineWidth(2); h_q2_sim_shape->SetLineColor(kRed);
    double max_q2 = std::max(h_q2_data_corr->GetMaximum(), h_q2_sim_shape->GetMaximum());
    h_q2_data_corr->SetMaximum(1.18 * max_q2);
    h_q2_data_corr->Draw("hist");
    h_q2_sim_shape->Draw("hist same");
    auto leg1 = new TLegend(0.20,0.75,0.58,0.90); style_legend(leg1, 0.035); leg1->AddEntry(h_q2_data_corr.get(),"Data (tgt corrected)","l"); leg1->AddEntry(h_q2_sim_shape.get(),"SIMC (area normalized)","l"); leg1->Draw();

    c1.cd(2);
    set_pub_pad(0.13, 0.04, 0.13, 0.09);
    h_xb_data_corr->SetTitle(";x_{B};Weighted counts (SIMC area norm.)");
    style_axes(h_xb_data_corr->GetXaxis(), h_xb_data_corr->GetYaxis());
    h_xb_data_corr->SetLineWidth(2); h_xb_data_corr->SetLineColor(kBlack);
    h_xb_sim_shape->SetLineWidth(2); h_xb_sim_shape->SetLineColor(kRed);
    double max_xb = std::max(h_xb_data_corr->GetMaximum(), h_xb_sim_shape->GetMaximum());
    h_xb_data_corr->SetMaximum(1.18 * max_xb);
    h_xb_data_corr->Draw("hist");
    h_xb_sim_shape->Draw("hist same");
    auto leg2 = new TLegend(0.20,0.75,0.58,0.90); style_legend(leg2, 0.035); leg2->AddEntry(h_xb_data_corr.get(),"Data (tgt corrected)","l"); leg2->AddEntry(h_xb_sim_shape.get(),"SIMC (area normalized)","l"); leg2->Draw();

    c1.cd(3);
    set_pub_pad(0.13, 0.04, 0.13, 0.09);
    h_tprime_data_corr->SetTitle(";t' [GeV^{2}];Weighted counts (SIMC area norm.)");
    style_axes(h_tprime_data_corr->GetXaxis(), h_tprime_data_corr->GetYaxis());
    h_tprime_data_corr->SetLineWidth(2); h_tprime_data_corr->SetLineColor(kBlack);
    h_tprime_sim_shape->SetLineWidth(2); h_tprime_sim_shape->SetLineColor(kRed);
    double max_tprime = std::max(h_tprime_data_corr->GetMaximum(), h_tprime_sim_shape->GetMaximum());
    h_tprime_data_corr->SetMaximum(1.18 * max_tprime);
    h_tprime_data_corr->Draw("hist");
    h_tprime_sim_shape->Draw("hist same");
    auto leg3 = new TLegend(0.20,0.75,0.58,0.90); style_legend(leg3, 0.035); leg3->AddEntry(h_tprime_data_corr.get(),"Data (tgt corrected)","l"); leg3->AddEntry(h_tprime_sim_shape.get(),"SIMC (area normalized)","l"); leg3->Draw();

    c1.cd(4);
    set_pub_pad(0.13, 0.04, 0.13, 0.09);
    h_phi_data_corr->SetTitle(";#phi [rad];Weighted counts (SIMC area norm.)");
    style_axes(h_phi_data_corr->GetXaxis(), h_phi_data_corr->GetYaxis());
    h_phi_data_corr->SetLineWidth(2); h_phi_data_corr->SetLineColor(kBlack);
    h_phi_sim_shape->SetLineWidth(2); h_phi_sim_shape->SetLineColor(kRed);
    double max_phi = std::max(h_phi_data_corr->GetMaximum(), h_phi_sim_shape->GetMaximum());
    h_phi_data_corr->SetMaximum(1.18 * max_phi);
    h_phi_data_corr->Draw("hist");
    h_phi_sim_shape->Draw("hist same");
    auto leg4 = new TLegend(0.20,0.75,0.58,0.90); style_legend(leg4, 0.035); leg4->AddEntry(h_phi_data_corr.get(),"Data (tgt corrected)","l"); leg4->AddEntry(h_phi_sim_shape.get(),"SIMC (area normalized)","l"); leg4->Draw();

    c1.Update();
    write_canvas_pdf_png(&c1, (fs::path(cfg.out_dir) / "global" / "global_1d_distributions").string());

    // --- DEBUG: Q2:xB bin grid overlay ---
    TCanvas c2("c2", "Q2:xB binning", 1000, 850);
    set_pub_pad(0.13, 0.15, 0.13, 0.08);
    h_q2_xb_data->SetTitle(";Q^{2} [GeV^{2}];x_{B}");
    style_axes(h_q2_xb_data->GetXaxis(), h_q2_xb_data->GetYaxis());
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

    // --- DEBUG: configured t' bin edges ---
    TCanvas c3("c3", "t' binning", 1000, 850);
    set_pub_pad(0.13, 0.04, 0.13, 0.08);
    h_tprime_data_corr->SetTitle(";t' [GeV^{2}];Weighted counts");
    style_axes(h_tprime_data_corr->GetXaxis(), h_tprime_data_corr->GetYaxis());
    h_tprime_data_corr->Draw("hist");
    for (size_t i = 1; i < tprime_edges.size() - 1; ++i) {
        TLine* l = new TLine(tprime_edges[i], 0, tprime_edges[i], h_tprime_data_corr->GetMaximum());
        l->SetLineColor(kGreen+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
    }
    c3.Update();
    write_canvas_pdf_png(&c3, (fs::path(cfg.out_dir) / "global" / "tprime_distribution_quantiles").string());

    // --- Diagnostic: phi bin edges ---
    TCanvas c4("c4", "phi binning", 1000, 850);
    set_pub_pad(0.13, 0.04, 0.13, 0.08);
    h_phi_data_corr->SetTitle(";#phi [rad];Weighted counts");
    style_axes(h_phi_data_corr->GetXaxis(), h_phi_data_corr->GetYaxis());
    h_phi_data_corr->Draw("hist");
    for (size_t i = 1; i < phi_edges.size() - 1; ++i) {
        TLine* l = new TLine(phi_edges[i], 0, phi_edges[i], h_phi_data_corr->GetMaximum());
        l->SetLineColor(kMagenta+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
    }
    c4.Update();
    write_canvas_pdf_png(&c4, (fs::path(cfg.out_dir) / "global" / "phi_binning_debug").string());

    // --- Existing 2D plots ---
    TCanvas c5("c5", "2d", 1500, 650);
    c5.Divide(2,1);
    c5.cd(1); set_pub_pad(0.12, 0.16, 0.14, 0.08); h_q2_xb_data->SetTitle(";Q^{2} [GeV^{2}];x_{B}"); style_axes(h_q2_xb_data->GetXaxis(), h_q2_xb_data->GetYaxis()); h_q2_xb_data->Draw("COLZ");
    c5.cd(2); set_pub_pad(0.12, 0.16, 0.14, 0.08); h_tprime_phi_data->SetTitle(";t' [GeV^{2}];#phi [rad]"); style_axes(h_tprime_phi_data->GetXaxis(), h_tprime_phi_data->GetYaxis()); h_tprime_phi_data->Draw("COLZ");
    c5.Update();
    write_canvas_pdf_png(&c5, (fs::path(cfg.out_dir) / "global" / "data_occupancy_2d").string());

    TCanvas c6("c6", "2d_sim", 1500, 650);
    c6.Divide(2,1);
    c6.cd(1); set_pub_pad(0.12, 0.16, 0.14, 0.08); h_q2_xb_sim->SetTitle(";Q^{2} [GeV^{2}];x_{B}"); style_axes(h_q2_xb_sim->GetXaxis(), h_q2_xb_sim->GetYaxis()); h_q2_xb_sim->Draw("COLZ");
    c6.cd(2); set_pub_pad(0.12, 0.16, 0.14, 0.08); h_tprime_phi_sim->SetTitle(";t' [GeV^{2}];#phi [rad]"); style_axes(h_tprime_phi_sim->GetXaxis(), h_tprime_phi_sim->GetYaxis()); h_tprime_phi_sim->Draw("COLZ");
    c6.Update();
    write_canvas_pdf_png(&c6, (fs::path(cfg.out_dir) / "global" / "simc_occupancy_2d").string());
}

void ExclPi0XSecAnalysis::make_mmiss_comparison_plots() {
    if (!cfg.diagnostics || (!cfg.write_pdf && !cfg.write_png)) return;
    const fs::path plot_dir = fs::path(cfg.out_dir) / "mmiss_comparison";
    fs::create_directories(plot_dir);
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                auto& m = mmiss_diagnostics[static_cast<size_t>(slice_index(it, iq, ix))];
                const std::string suffix = "t" + std::to_string(it) + "_q" +
                    std::to_string(iq) + "_x" + std::to_string(ix);
                auto sim = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(
                    m.exclusive->Clone(("plot_mmiss_sim_" + suffix).c_str())));
                sim->SetDirectory(nullptr);
                sim->Scale(m.exclusive_shape_scale);

                TCanvas canvas(("c_mmiss_" + suffix).c_str(), "Missing-mass comparison", 1100, 700);
                set_pub_pad(0.12, 0.03, 0.13, 0.10);
                m.data->SetMarkerStyle(20); m.data->SetMarkerSize(0.65);
                m.data->SetLineColor(kBlack);
                sim->SetLineColor(kBlue + 1); sim->SetLineWidth(2); sim->SetLineStyle(2);
                const double ymax = 1.35 * std::max({m.data->GetMaximum(), sim->GetMaximum(), 1e-6});
                m.data->SetMinimum(0.0); m.data->SetMaximum(ymax);
                m.data->SetTitle(("Full reconstructed missing mass: " +
                                  format_slice_bin_label(it, iq, ix)).c_str());
                style_axes(m.data->GetXaxis(), m.data->GetYaxis());
                m.data->Draw("E1"); sim->Draw("HIST SAME"); m.data->Draw("E1 SAME");
                TLine cut_lo(cfg.mmiss_lower_gev, 0.0, cfg.mmiss_lower_gev, ymax);
                TLine cut_hi(cfg.mmiss_upper_gev, 0.0, cfg.mmiss_upper_gev, ymax);
                for (TLine* cut : {&cut_lo, &cut_hi}) {
                    cut->SetLineColor(kMagenta + 2); cut->SetLineStyle(3);
                    cut->SetLineWidth(2); cut->Draw();
                }
                TLegend legend(0.57, 0.68, 0.96, 0.89);
                style_legend(&legend, 0.034);
                legend.AddEntry(m.data.get(), "Data: all candidates", "lep");
                legend.AddEntry(sim.get(), "Original exclusive SIMC: area normalized", "l");
                legend.AddEntry(&cut_hi, "Extraction window", "l");
                legend.Draw();
                TLatex note; note.SetNDC(true); note.SetTextFont(42); note.SetTextSize(0.027);
                note.DrawLatex(0.14, 0.86, "Shape diagnostic only; scale does not enter extraction");
                canvas.Update();
                write_canvas_pdf_png(&canvas, (plot_dir / ("mmiss_" + suffix)).string());
            }
        }
    }
}

void ExclPi0XSecAnalysis::make_yield_diagnostic_plots() {
    if (!cfg.diagnostics || (!cfg.write_pdf && !cfg.write_png)) return;
    // Deliberately do not area-normalize these maps: unlike the 1D shape and
    // missing-mass overlays, they show the actual per-bin yields entering
    // Y_data/Y_SIMC and final reweighted-SIMC effective statistics.
    const fs::path plot_dir = fs::path(cfg.out_dir) / "diagnostics";
    fs::create_directories(plot_dir);
    TDirectory* root_dir = gDirectory;
    TDirectory* output_dir = fout->GetDirectory("diagnostics");
    if (!output_dir) output_dir = fout->mkdir("diagnostics");

    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix < cfg.n_xb; ++ix) {
            const std::string suffix = "q" + std::to_string(iq) + "_x" + std::to_string(ix);
            auto make_map = [&](const std::string& name, const std::string& title) {
                auto h = std::make_unique<TH2D>(
                    (name + "_" + suffix).c_str(), title.c_str(),
                    cfg.n_tprime, t_edges.data(), cfg.n_phi, phi_edges.data());
                h->SetDirectory(nullptr);
                return h;
            };
            auto data = make_map("data_yield", "Data yield;t [GeV^{2}];#phi [rad]");
            auto sim = make_map("simc_yield", "SIMC yield;t [GeV^{2}];#phi [rad]");
            auto ratio = make_map("data_simc_ratio", "Data/SIMC yield ratio;t [GeV^{2}];#phi [rad]");
            auto neff = make_map("simc_effective_events", "SIMC effective events;t [GeV^{2}];#phi [rad]");
            for (int it = 0; it < cfg.n_tprime; ++it) {
                const auto& s = slice(it, iq, ix);
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    const auto& p = s.phi[ip];
                    const int bx = it + 1, by = ip + 1;
                    data->SetBinContent(bx, by, p.data);
                    sim->SetBinContent(bx, by, p.sim);
                    ratio->SetBinContent(bx, by, p.n_data > 0 && p.sim > 0.0 && s.reference_supported ? p.data / p.sim : std::numeric_limits<double>::quiet_NaN());
                    neff->SetBinContent(bx, by,
                        p.sim_sumw2 > 0.0 ? p.sim * p.sim / p.sim_sumw2 : 0.0);
                }
            }

            TCanvas canvas(("c_yield_diag_" + suffix).c_str(), "Yield diagnostics", 1500, 1050);
            canvas.Divide(2, 2);
            std::array<TH2D*, 4> maps{data.get(), sim.get(), ratio.get(), neff.get()};
            for (size_t i = 0; i < maps.size(); ++i) {
                canvas.cd(static_cast<int>(i + 1));
                set_pub_pad(0.12, 0.16, 0.13, 0.10);
                style_axes(maps[i]->GetXaxis(), maps[i]->GetYaxis(), 0.035, 0.042);
                maps[i]->Draw("COLZ TEXT");
            }
            canvas.Update();
            write_canvas_pdf_png(&canvas, (plot_dir / ("yield_support_" + suffix)).string());

            if (output_dir) {
                output_dir->cd();
                for (TH2D* h : maps) h->Write();
                if (root_dir) root_dir->cd();
            }
        }
    }
}

void ExclPi0XSecAnalysis::make_epsilon_plots() {
    fs::create_directories(fs::path(cfg.out_dir) / "global");

    TCanvas c_eps_t("c_eps_t", "epsilon_vs_tprime", 1200, 900);
    set_pub_pad(0.13, 0.05, 0.13, 0.08);
    TH1D hframe("hframe_eps_t", ";t [GeV^{2}];Virtual photon #epsilon", 100, cfg.t_min, cfg.t_max);
    hframe.SetMinimum(0.0);
    hframe.SetMaximum(1.05);
    style_axes(hframe.GetXaxis(), hframe.GetYaxis());
    hframe.Draw("AXIS");

    std::vector<std::unique_ptr<TGraphErrors>> graphs;
    TLegend leg(0.48, 0.60, 0.89, 0.88);
    style_legend(&leg, 0.030);

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
                xt.push_back(0.5 * (t_edges[it] + t_edges[it + 1]));
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
    lat.SetTextSize(0.034);
    lat.DrawLatex(0.16, 0.92, Form("E_{beam} = %.3f GeV", cfg.ebeam));
    c_eps_t.Update();
    write_canvas_pdf_png(&c_eps_t, (fs::path(cfg.out_dir) / "global" / "epsilon_vs_tprime_by_q2_xb").string());

    const int nslices = cfg.n_tprime * cfg.n_q2 * cfg.n_xb;
    TH1D h_eps_slice("h_eps_slice", "#epsilon by (t',Q^{2},x_{B}) slice;;#epsilon", nslices, 0.5, nslices + 0.5);
    int bin = 1;
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                h_eps_slice.SetBinContent(bin, s.epsilon);
                std::ostringstream bl;
                bl << std::fixed << std::setprecision(2)
                   << "t[" << t_edges[it] << "," << t_edges[it + 1] << "] "
                   << "Q2[" << q2_edges[iq] << "," << q2_edges[iq + 1] << "] "
                   << "xB[" << xb_edges_by_q2[iq][ix] << "," << xb_edges_by_q2[iq][ix + 1] << "]";
                h_eps_slice.GetXaxis()->SetBinLabel(bin, bl.str().c_str());
                ++bin;
            }
        }
    }

    TCanvas c_eps_idx("c_eps_idx", "epsilon_by_slice", 1500, 800);
    h_eps_slice.SetMinimum(0.0);
    h_eps_slice.SetMaximum(1.05);
    h_eps_slice.SetStats(0);
    h_eps_slice.SetMarkerStyle(20);
    h_eps_slice.SetMarkerSize(0.9);
    h_eps_slice.SetLineWidth(2);
    style_axes(h_eps_slice.GetXaxis(), h_eps_slice.GetYaxis(), 0.034, 0.044);
    h_eps_slice.GetXaxis()->SetLabelSize(0.013);
    h_eps_slice.GetXaxis()->LabelsOption("v");
    c_eps_idx.SetLeftMargin(0.11);
    c_eps_idx.SetRightMargin(0.03);
    c_eps_idx.SetBottomMargin(0.30);
    c_eps_idx.SetTopMargin(0.08);
    c_eps_idx.SetTicks(1, 1);
    h_eps_slice.Draw("P HIST");
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
                TCanvas c(("c_"+tag.str()).c_str(), tag.str().c_str(), 1700, 1150);
                c.Divide(2,2);

                std::vector<double> phi_c(cfg.n_phi), ydata(cfg.n_phi), yerr(cfg.n_phi), ysim(cfg.n_phi), ysimerr(cfg.n_phi), ratio(cfg.n_phi), ratioerr(cfg.n_phi), xsec(cfg.n_phi), xsecerr(cfg.n_phi);
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    phi_c[ip] = 0.5 * (phi_edges[ip] + phi_edges[ip+1]);
                    const PhiBin& pb = s.phi[ip];
                    ydata[ip] = pb.data; yerr[ip] = std::sqrt(std::max(0.0, pb.data_sumw2));
                    ysim[ip]  = pb.sim;  ysimerr[ip] = std::sqrt(std::max(0.0, pb.sim_sumw2));
                    ratio[ip] = pb.ratio; ratioerr[ip] = pb.ratio_err;
                    xsec[ip] = pb.xsec; xsecerr[ip] = pb.xsec_err;
                }

                c.cd(1);
                set_pub_pad(0.13, 0.04, 0.13, 0.11);
                auto gdata = new TGraphErrors(cfg.n_phi, phi_c.data(), ydata.data(), nullptr, yerr.data());
                auto gsim  = new TGraphErrors(cfg.n_phi, phi_c.data(), ysim.data(), nullptr, ysimerr.data());
                gdata->SetTitle("Data and SIMC yields;#phi [rad];Weighted counts");
                gdata->SetMarkerStyle(20); gdata->SetMarkerSize(0.85); gdata->SetLineColor(kBlack); gdata->SetMarkerColor(kBlack);
                gsim->SetMarkerStyle(24); gsim->SetMarkerSize(0.85); gsim->SetLineColor(kRed + 1); gsim->SetMarkerColor(kRed + 1);
                gdata->Draw("AP");
                style_axes(gdata->GetXaxis(), gdata->GetYaxis(), 0.036, 0.043);
                gsim->Draw("P SAME");
                auto l1 = new TLegend(0.68,0.76,0.91,0.90); style_legend(l1, 0.035); l1->AddEntry(gdata,"Data","p"); l1->AddEntry(gsim,"SIMC","p"); l1->Draw();
                draw_slice_bin_label(it, iq, ix);

                c.cd(2);
                set_pub_pad(0.13, 0.04, 0.13, 0.11);
                auto grat = new TGraphErrors(cfg.n_phi, phi_c.data(), ratio.data(), nullptr, ratioerr.data());
                grat->SetTitle("Ratio data/SIMC;#phi [rad];Ratio");
                grat->SetMarkerStyle(20); grat->SetMarkerSize(0.85); grat->SetLineColor(kBlack); grat->SetMarkerColor(kBlack); grat->Draw("AP");
                style_axes(grat->GetXaxis(), grat->GetYaxis(), 0.036, 0.043);
                if (s.fit_ratio.ok) {
                    auto f = new TF1(("fr_"+tag.str()).c_str(), "[0] + [1]*cos(x) + [2]*cos(2*x)", cfg.phi_min, cfg.phi_max);
                    f->SetParameters(s.fit_ratio.p[0], s.fit_ratio.p[1], s.fit_ratio.p[2]);
                    f->SetLineColor(kRed + 1);
                    f->SetLineWidth(3);
                    f->Draw("same");
                }
                draw_slice_bin_label(it, iq, ix);

                c.cd(3);
                set_pub_pad(0.13, 0.04, 0.13, 0.11);
                if (s.has_model_xsec) {
                    auto gx = new TGraphErrors(cfg.n_phi, phi_c.data(), xsec.data(), nullptr, xsecerr.data());
                    gx->SetTitle("Extracted #sigma with fit components;#phi [rad];Cross section");
                    gx->SetMarkerStyle(20);
                    gx->SetMarkerSize(0.85);
                    gx->SetMarkerColor(kBlack);
                    gx->SetLineColor(kBlack);

                    if (s.fit_xsec.ok && s.fit_xsec.p.size() >= 3) {
                        const double p0 = s.fit_xsec.p[0];
                        const double p1 = s.fit_xsec.p[1];
                        const double p2 = s.fit_xsec.p[2];
                        const int ncurve = 361;
                        std::vector<double> ph_curve(ncurve);
                        std::vector<double> y_tot(ncurve), y_u(ncurve), y_lt(ncurve), y_tt(ncurve);

                        double ymin = std::numeric_limits<double>::infinity();
                        double ymax = -std::numeric_limits<double>::infinity();

                        for (int ip = 0; ip < cfg.n_phi; ++ip) {
                            if (!std::isfinite(xsec[ip])) continue;
                            ymin = std::min(ymin, xsec[ip] - std::max(0.0, xsecerr[ip]));
                            ymax = std::max(ymax, xsec[ip] + std::max(0.0, xsecerr[ip]));
                        }

                        for (int isamp = 0; isamp < ncurve; ++isamp) {
                            const double ph = cfg.phi_min + (cfg.phi_max - cfg.phi_min) * (static_cast<double>(isamp) / static_cast<double>(ncurve - 1));
                            ph_curve[isamp] = ph;
                            y_u[isamp] = p0;
                            y_lt[isamp] = p1 * std::cos(ph);
                            y_tt[isamp] = p2 * std::cos(2.0 * ph);
                            y_tot[isamp] = y_u[isamp] + y_lt[isamp] + y_tt[isamp];
                            ymin = std::min({ymin, y_u[isamp], y_lt[isamp], y_tt[isamp], y_tot[isamp]});
                            ymax = std::max({ymax, y_u[isamp], y_lt[isamp], y_tt[isamp], y_tot[isamp]});
                        }

                        if (std::isfinite(ymin) && std::isfinite(ymax)) {
                            const double span = std::max(1e-12, ymax - ymin);
                            gx->SetMinimum(ymin - 0.12 * span);
                            gx->SetMaximum(ymax + 0.55 * span);
                        }

                        gx->Draw("AP");
                        style_axes(gx->GetXaxis(), gx->GetYaxis(), 0.036, 0.043);

                        auto g_tot = new TGraphErrors(ncurve, ph_curve.data(), y_tot.data(), nullptr, nullptr);
                        auto g_u = new TGraphErrors(ncurve, ph_curve.data(), y_u.data(), nullptr, nullptr);
                        auto g_lt = new TGraphErrors(ncurve, ph_curve.data(), y_lt.data(), nullptr, nullptr);
                        auto g_tt = new TGraphErrors(ncurve, ph_curve.data(), y_tt.data(), nullptr, nullptr);

                        g_tot->SetLineColor(kBlue + 1); g_tot->SetLineWidth(3);
                        g_u->SetLineColor(kBlack); g_u->SetLineStyle(2); g_u->SetLineWidth(2);
                        g_lt->SetLineColor(kRed + 1); g_lt->SetLineStyle(7); g_lt->SetLineWidth(2);
                        g_tt->SetLineColor(kGreen + 2); g_tt->SetLineStyle(9); g_tt->SetLineWidth(2);

                        g_tot->Draw("L SAME");
                        g_u->Draw("L SAME");
                        g_lt->Draw("L SAME");
                        g_tt->Draw("L SAME");

                        auto l3 = new TLegend(0.52,0.64,0.89,0.87);
                        style_legend(l3, 0.027);
                        l3->AddEntry(gx, "Extracted #sigma", "p");
                        l3->AddEntry(g_tot, "Total fit", "l");
                        l3->AddEntry(g_u, "Constant term", "l");
                        l3->AddEntry(g_lt, "cos#phi term", "l");
                        l3->AddEntry(g_tt, "cos2#phi term", "l");
                        l3->Draw();

                    } else {
                        gx->Draw("AP");
                        style_axes(gx->GetXaxis(), gx->GetYaxis(), 0.036, 0.043);
                    }
                } else {
                    auto g2 = new TGraphErrors(cfg.n_phi, phi_c.data(), ratio.data(), nullptr, ratioerr.data());
                    g2->SetTitle("Ratio only (no model xsec branch);#phi [rad];Ratio");
                    g2->SetMarkerStyle(20); g2->SetMarkerSize(0.85); g2->Draw("AP");
                    style_axes(g2->GetXaxis(), g2->GetYaxis(), 0.036, 0.043);
                }
                c.cd(4);
                set_pub_pad(0.13, 0.04, 0.13, 0.11);
                auto h = new TH1D(("hltp_"+tag.str()).c_str(), "LT' status;#phi [rad];", 10, cfg.phi_min, cfg.phi_max);
                h->SetMinimum(0.0);
                h->SetMaximum(1.0);
                h->SetStats(0);
                h->Draw("AXIS");
                style_axes(h->GetXaxis(), h->GetYaxis(), 0.036, 0.043);
                TLatex note;
                note.SetNDC(true);
                note.SetTextFont(42);
                note.SetTextSize(0.045);
                note.SetTextAlign(22);
                note.DrawLatex(0.52, 0.56, "LT' unavailable in SIMC-model extraction");
                note.SetTextSize(0.031);
                note.DrawLatex(0.52, 0.47, "Requires helicity luminosities and beam polarization");
                draw_slice_bin_label(it, iq, ix);

                c.Update();
                write_canvas_pdf_png(&c, (fs::path(cfg.out_dir) / "slices" / ("slice_" + tag.str())).string());
            }
        }
    }
}

void ExclPi0XSecAnalysis::make_sigma_vs_t_plots() {
    fs::create_directories(fs::path(cfg.out_dir) / "sigma_vs_t");
    constexpr double unit_scale = 1.0e9; // microbarn/MeV^2 -> nb/GeV^2
    const std::array<const char*, 3> labels{"#sigma_{U}", "#sigma_{LT}", "#sigma_{TT}"};
    const std::array<int, 3> colors{kBlack, kRed + 1, kBlue + 1};
    const std::array<int, 3> markers{20, 21, 22};

    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix < cfg.n_xb; ++ix) {
            std::vector<double> x, ex;
            std::array<std::vector<double>, 3> sigma;
            std::array<std::vector<double>, 3> sigma_err;

            x.reserve(cfg.n_tprime);
            ex.reserve(cfg.n_tprime);
            for (int term = 0; term < 3; ++term) {
                sigma[term].reserve(cfg.n_tprime);
                sigma_err[term].reserve(cfg.n_tprime);
            }

            for (int it = 0; it < cfg.n_tprime; ++it) {
                const SliceResult& s = slice(it, iq, ix);
                if (!s.fit_xsec.ok) continue;

                const double tcenter = 0.5 * (t_edges[it] + t_edges[it + 1]);
                const double thalf = 0.5 * (t_edges[it + 1] - t_edges[it]);

                x.push_back(-tcenter);
                ex.push_back(std::fabs(thalf));
                sigma[0].push_back(unit_scale * s.fit_xsec.sigmaU);
                sigma_err[0].push_back(unit_scale * std::max(0.0, s.fit_xsec.sigmaU_err));
                sigma[1].push_back(unit_scale * s.fit_xsec.sigmaTL);
                sigma_err[1].push_back(unit_scale * std::max(0.0, s.fit_xsec.sigmaTL_err));
                sigma[2].push_back(unit_scale * s.fit_xsec.sigmaTT);
                sigma_err[2].push_back(unit_scale * std::max(0.0, s.fit_xsec.sigmaTT_err));
            }

            if (x.empty()) continue;

            std::ostringstream tag;
            tag << "q" << iq << "_x" << ix;
            TCanvas c(("c_sigma_vs_t_" + tag.str()).c_str(), "sigma_vs_t", 1800, 600);
            c.Divide(3, 1);
            // Independent y ranges prevent sigma_U from visually compressing
            // the smaller, signed sigma_LT and sigma_TT coefficients.
            std::array<std::unique_ptr<TGraphErrors>, 3> graphs;
            std::array<std::unique_ptr<TLine>, 3> zero_lines;
            for (int term = 0; term < 3; ++term) {
                c.cd(term + 1);
                set_pub_pad(0.17, 0.04, 0.14, 0.12);

                double ymin = std::numeric_limits<double>::infinity();
                double ymax = -std::numeric_limits<double>::infinity();
                for (size_t i = 0; i < x.size(); ++i) {
                    const double err = std::isfinite(sigma_err[term][i]) ? sigma_err[term][i] : 0.0;
                    ymin = std::min(ymin, sigma[term][i] - err);
                    ymax = std::max(ymax, sigma[term][i] + err);
                }
                if (term > 0) {
                    ymin = std::min(ymin, 0.0);
                    ymax = std::max(ymax, 0.0);
                }
                const double yspan = std::max(1e-9, ymax - ymin);

                graphs[term] = std::make_unique<TGraphErrors>(
                    static_cast<int>(x.size()), x.data(), sigma[term].data(),
                    ex.data(), sigma_err[term].data());
                auto* graph = graphs[term].get();
                graph->SetTitle(Form("%s;-t [GeV^{2}];response [nb/GeV^{2}]", labels[term]));
                graph->SetMinimum(ymin - 0.18 * yspan);
                graph->SetMaximum(ymax + 0.24 * yspan);
                graph->SetLineColor(colors[term]);
                graph->SetMarkerColor(colors[term]);
                graph->SetMarkerStyle(markers[term]);
                graph->SetMarkerSize(1.0);
                graph->SetLineWidth(2);
                graph->Draw("AP");
                style_axes(graph->GetXaxis(), graph->GetYaxis(), 0.04, 0.045);

                if (term > 0) {
                    const double xmin = graph->GetXaxis()->GetXmin();
                    const double xmax = graph->GetXaxis()->GetXmax();
                    zero_lines[term] = std::make_unique<TLine>(xmin, 0.0, xmax, 0.0);
                    zero_lines[term]->SetLineColor(kGray + 2);
                    zero_lines[term]->SetLineStyle(2);
                    zero_lines[term]->Draw();
                    graph->Draw("P SAME");
                }
            }

            c.cd(1);
            TLatex lat1;
            lat1.SetNDC(true);
            lat1.SetTextFont(42);
            lat1.SetTextSize(0.030);
            lat1.DrawLatex(0.20, 0.91,
                           Form("Q^{2} #in [%.3f, %.3f] GeV^{2}", q2_edges[iq], q2_edges[iq + 1]));
            lat1.DrawLatex(0.20, 0.85,
                           Form("x_{B} #in [%.3f, %.3f]", xb_edges_by_q2[iq][ix], xb_edges_by_q2[iq][ix + 1]));

            c.Update();
            write_canvas_pdf_png(&c, (fs::path(cfg.out_dir) / "sigma_vs_t" / ("sigma_terms_vs_minus_t_" + tag.str())).string());
        }
    }
}

void ExclPi0XSecAnalysis::make_partons_projection_plots() {
    if (!cfg.partons_projection || (!cfg.write_pdf && !cfg.write_png)) return;
    fs::create_directories(fs::path(cfg.out_dir) / "partons_projection");
    constexpr double unit_scale = 1.0e9; // microbarn/MeV^2 -> nb/GeV^2
    const std::array<const char*, 3> labels{"#sigma_{U}", "#sigma_{LT}", "#sigma_{TT}"};

    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix < cfg.n_xb; ++ix) {
            std::ostringstream tag;
            tag << "q" << iq << "_x" << ix;
            TCanvas canvas(("c_partons_" + tag.str()).c_str(), "PARTONS comparison", 1800, 600);
            canvas.Divide(3, 1);
            // Pads retain pointers to drawn objects; keep every graph and
            // legend alive until the completed canvas is written.
            std::array<std::unique_ptr<TGraphErrors>, 3> data_graphs;
            std::array<std::unique_ptr<TGraphErrors>, 3> model_graphs;
            std::array<std::unique_ptr<TLegend>, 3> legends;
            for (int term = 0; term < 3; ++term) {
                std::vector<double> x, ex, measured, measured_err;
                std::vector<double> theory_x, theory, theory_err;
                for (int it = 0; it < cfg.n_tprime; ++it) {
                    const SliceResult& s = slice(it, iq, ix);
                    if (!s.fit_xsec.ok) continue;
                    x.push_back(-s.reference_model.t);
                    ex.push_back(0.0);
                    measured.push_back(unit_scale * (term == 0 ? s.fit_xsec.sigmaU :
                                                      term == 1 ? s.fit_xsec.sigmaTL : s.fit_xsec.sigmaTT));
                    measured_err.push_back(unit_scale * std::max(0.0,
                        term == 0 ? s.fit_xsec.sigmaU_err :
                        term == 1 ? s.fit_xsec.sigmaTL_err : s.fit_xsec.sigmaTT_err));
                    if (s.partons_ok) {
                        theory_x.push_back(-s.reference_model.t);
                        theory.push_back(unit_scale * (term == 0 ? s.partons_sigmaU :
                                                       term == 1 ? s.partons_sigmaLT : s.partons_sigmaTT));
                        theory_err.push_back(0.0);
                    }
                }
                if (x.empty()) continue;
                canvas.cd(term + 1);
                set_pub_pad(0.17, 0.04, 0.14, 0.12);
                data_graphs[term] = std::make_unique<TGraphErrors>(
                    static_cast<int>(x.size()), x.data(), measured.data(), ex.data(), measured_err.data());
                auto* data_graph = data_graphs[term].get();
                data_graph->SetTitle(Form("%s;-#LTt'_{SIMC}#GT [GeV^{2}];response [nb/GeV^{2}]", labels[term]));
                data_graph->SetMarkerStyle(20); data_graph->SetMarkerColor(kBlack); data_graph->SetLineColor(kBlack);
                data_graph->Draw("AP");
                style_axes(data_graph->GetXaxis(), data_graph->GetYaxis(), 0.04, 0.045);
                model_graphs[term] = std::make_unique<TGraphErrors>(
                    static_cast<int>(theory_x.size()), theory_x.data(), theory.data(), theory_err.data(), theory_err.data());
                auto* model_graph = model_graphs[term].get();
                model_graph->SetMarkerStyle(25); model_graph->SetMarkerSize(1.2);
                model_graph->SetMarkerColor(kMagenta + 2); model_graph->SetLineColor(kMagenta + 2);
                if (!theory_x.empty()) model_graph->Draw("P SAME");
                legends[term] = std::make_unique<TLegend>(0.48, 0.72, 0.92, 0.88);
                style_legend(legends[term].get(), 0.035);
                legends[term]->AddEntry(data_graph, "Data/SIMC extraction", "p");
                if (!theory_x.empty()) legends[term]->AddEntry(model_graph, "PARTONS GK06/GPDGK19", "p");
                legends[term]->Draw();
            }
            canvas.Update();
            write_canvas_pdf_png(&canvas,
                (fs::path(cfg.out_dir) / "partons_projection" /
                 ("partons_coefficients_vs_minus_t_" + tag.str())).string());
        }
    }
}

void ExclPi0XSecAnalysis::write_results() {
    fout->cd();
    // Save bin edges
    TVectorD vphi(phi_edges.size()), vt(tprime_edges.size()), vphys_t(t_edges.size()), vq(q2_edges.size()), vx(xb_edges.size());
    TMatrixD vxb2d(cfg.n_q2, cfg.n_xb + 1);
    for (size_t i = 0; i < phi_edges.size(); ++i) vphi[i] = phi_edges[i];
    for (size_t i = 0; i < tprime_edges.size(); ++i) vt[i] = tprime_edges[i];
    for (size_t i = 0; i < t_edges.size(); ++i) vphys_t[i] = t_edges[i];
    for (size_t i = 0; i < q2_edges.size(); ++i) vq[i] = q2_edges[i];
    for (size_t i = 0; i < xb_edges.size(); ++i) vx[i] = xb_edges[i];
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix <= cfg.n_xb; ++ix) {
            vxb2d(iq, ix) = xb_edges_by_q2[iq][ix];
        }
    }
    vphi.Write("phi_edges");
    vt.Write("tprime_edges");
    vphys_t.Write("t_edges");
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
    if (cfg.diagnostics) {
        TDirectory* saved = gDirectory;
        TDirectory* mass_dir = fout->GetDirectory("missing_mass_diagnostics");
        if (!mass_dir) mass_dir = fout->mkdir("missing_mass_diagnostics");
        if (mass_dir) {
            mass_dir->cd();
            for (size_t i = 0; i < mmiss_diagnostics.size(); ++i) {
                mmiss_diagnostics[i].data->Write();
                mmiss_diagnostics[i].exclusive->Write();
                TParameter<double>(("exclusive_shape_scale_slice" + std::to_string(i)).c_str(),
                                   mmiss_diagnostics[i].exclusive_shape_scale).Write();
            }
            if (saved) saved->cd();
        }
    }

    TParameter<Long64_t>("n_data_total", cutflow.n_data_total).Write();
    TParameter<Long64_t>("n_data_pass", cutflow.n_data_pass).Write();
    TParameter<Long64_t>("n_data_inrange", cutflow.n_data_inrange).Write();
    TParameter<Long64_t>("n_sim_total", cutflow.n_sim_total).Write();
    TParameter<Long64_t>("n_sim_pass", cutflow.n_sim_pass).Write();
    TParameter<Long64_t>("n_sim_inrange", cutflow.n_sim_inrange).Write();

    TParameter<Long64_t>("n_sim_bad_full_weight", cutflow.n_sim_bad_full_weight).Write();
    TParameter<Long64_t>("n_sim_bad_sigcm", cutflow.n_sim_bad_sigcm).Write();
    TParameter<Long64_t>("n_sim_bad_vertex", cutflow.n_sim_bad_vertex).Write();
    TParameter<Long64_t>("n_sim_epsilon_fallback", cutflow.n_sim_epsilon_fallback).Write();
    TParameter<Long64_t>("n_sim_bad_model", cutflow.n_sim_bad_model).Write();
    TParameter<Long64_t>("n_sim_cached", reweight_events.size()).Write();
    TParameter<Long64_t>("default_sigcm_mismatch_gt_1pct", default_sigcm_mismatch_count).Write();
    TParameter<double>("default_sigcm_max_relative_difference", default_sigcm_max_relative_difference).Write();
    TParameter<double>("model_chi2_before", model_chi2_before).Write();
    TParameter<double>("model_chi2_after", model_chi2_after).Write();
    TParameter<double>("model_fit_pvalue", model_fit_pvalue).Write();
    TParameter<int>("model_fit_ndf", model_fit_ndf).Write();
    TParameter<int>("model_fit_bins", model_fit_bins).Write();
    TObjString(model_fit_status.c_str()).Write("model_fit_status");
    if (model_covariance.GetNrows()) model_covariance.Write("model_fit_covariance");
    {
        const auto specs=nps_pi0_reweight::choose(cfg.model_identifier).parameter_spec();
        std::string name;
        double value=0,error=0,lower=0,upper=0;
        int is_free=0;
        TTree tree("model_parameters","Fortran pi0 coefficient vector and fit");
        tree.Branch("name",&name); tree.Branch("value",&value);
        tree.Branch("error",&error); tree.Branch("lower",&lower);
        tree.Branch("upper",&upper); tree.Branch("free",&is_free);
        for (size_t i=0;i<specs.size();++i) {
            name=specs[i].name; value=model_parameters[i]; error=model_errors[i];
            lower=specs[i].lower; upper=specs[i].upper;
            is_free=std::find(free_model_indices.begin(),free_model_indices.end(),int(i))!=free_model_indices.end();
            tree.Fill();
        }
        tree.Write();
    }

    TParameter<int>("sigmaTLp_available", 0).Write();
    TObjString(kSigmaTLpStatus).Write("sigmaTLp_status");
    TParameter<int>("normalize_mmiss", cfg.normalize_mmiss ? 1 : 0).Write();
    TParameter<double>("yield_norm_data", yield_norm_data).Write();
    TParameter<double>("yield_norm_data_sumw2", yield_norm_data_sumw2).Write();
    TParameter<double>("yield_norm_sim", yield_norm_sim).Write();
    TParameter<double>("yield_norm_sim_sumw2", yield_norm_sim_sumw2).Write();
    TParameter<double>("yield_norm_scale", yield_norm_scale).Write();
    TParameter<double>("yield_norm_scale_err_independent", yield_norm_scale_err).Write();

    {
        int it_ref = 0, iq_ref = 0, ix_ref = 0;
        double q2_ref = 0, q2_direct=0, xb_ref = 0, tp_ref = 0, t_ref = 0, w_ref = 0, eps_ref = 0;
        TTree references("model_reference", "Common SIMC model reference per slice");
        references.Branch("it", &it_ref); references.Branch("iq", &iq_ref); references.Branch("ix", &ix_ref);
        references.Branch("Q2", &q2_ref); references.Branch("Q2_direct_base_mean", &q2_direct); references.Branch("xB", &xb_ref);
        references.Branch("tprime", &tp_ref); references.Branch("t", &t_ref);
        references.Branch("W", &w_ref); references.Branch("epsilon", &eps_ref);
        for (it_ref = 0; it_ref < cfg.n_tprime; ++it_ref)
            for (iq_ref = 0; iq_ref < cfg.n_q2; ++iq_ref)
                for (ix_ref = 0; ix_ref < cfg.n_xb; ++ix_ref) {
                    const auto& s = slice(it_ref, iq_ref, ix_ref);
                    if (!s.reference_supported) continue;
                    const auto& m = s.reference_model;
                    q2_ref = m.q2; q2_direct=s.q2_direct_base_mean; xb_ref = m.xb; tp_ref = m.tprime;
                    t_ref = m.t; w_ref = m.w; eps_ref = m.epsilon;
                    references.Fill();
                }
        references.Write();
    }

    if (cfg.partons_projection) {
        int it_out = 0, iq_out = 0, ix_out = 0, ok = 0;
        double q2 = 0.0, xb = 0.0, t = 0.0, tp = 0.0, eps = 0.0, flux = 0.0;
        double sigma_u = 0.0, sigma_lt = 0.0, sigma_tt = 0.0;
        double t_ref = 0.0;
        TTree model_tree("partons_projection", "GK06/GPDGK19 at extraction reference points");
        model_tree.Branch("it", &it_out); model_tree.Branch("iq", &iq_out); model_tree.Branch("ix", &ix_out);
        model_tree.Branch("partons_ok", &ok);
        model_tree.Branch("mean_q2_sim", &q2); model_tree.Branch("mean_xb_sim", &xb);
        model_tree.Branch("mean_t_sim", &t); model_tree.Branch("mean_tprime_sim", &tp);
        model_tree.Branch("t_ref", &t_ref);
        model_tree.Branch("partons_epsilon", &eps);
        model_tree.Branch("partons_electron_flux_xbq2", &flux);
        model_tree.Branch("partons_sigmaU", &sigma_u);
        model_tree.Branch("partons_sigmaLT", &sigma_lt);
        model_tree.Branch("partons_sigmaTT", &sigma_tt);
        for (it_out = 0; it_out < cfg.n_tprime; ++it_out) {
            for (iq_out = 0; iq_out < cfg.n_q2; ++iq_out) {
                for (ix_out = 0; ix_out < cfg.n_xb; ++ix_out) {
                    const auto& s = slice(it_out, iq_out, ix_out);
                    ok = s.partons_ok ? 1 : 0;
                    q2 = s.mean_q2_sim; xb = s.mean_xb_sim; t = s.mean_t_sim; tp = s.mean_tprime_sim;
                    t_ref = s.reference_model.t;
                    eps = s.partons_epsilon; flux = s.partons_electron_flux_xbq2;
                    sigma_u = s.partons_sigmaU; sigma_lt = s.partons_sigmaLT; sigma_tt = s.partons_sigmaTT;
                    model_tree.Fill();
                }
            }
        }
        model_tree.Write();
    }

    // write slice fit summaries as TObjString for portability
    std::ostringstream meta;
    meta << std::setprecision(10);
    meta << "tprime_edges_diagnostic:";
    for (double x : tprime_edges) meta << " " << x;
    meta << "\nq2_edges:";
    for (double x : q2_edges) meta << " " << x;
    meta << "\nxb_edges:";
    for (double x : xb_edges) meta << " " << x;
    meta << "\ndiamond_xb_q2_vertices=" << xsec_diamond_vertices_text(cfg);
    meta << "\nxb_edges_by_q2:";
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        meta << "\n  iq=" << iq << ":";
        for (double x : xb_edges_by_q2[iq]) meta << " " << x;
    }
    meta << "\nphi_edges:";
    for (double x : phi_edges) meta << " " << x;
    meta << "\nmodel_xsec_branch=" << (has_model_xsec ? model_xsec_branch : "NONE");
    meta << "\nmodel_xsec_mode=" << cfg.model_xsec_mode;
    meta << "\nconfigured_kinematic=" << cfg.configured_kinematic;
    meta << "\nhelicity=" << (has_helicity ? "yes" : "no");
    meta << "\nsigmaTLp_available=no";
    meta << "\nsigmaTLp_status=" << kSigmaTLpStatus;
    meta << "\nmissing_mass_cut=data:mmiss_all,sim:mmiss," << cfg.mmiss_lower_gev << "<Mmiss<" << cfg.mmiss_upper_gev;
    meta << "\ntarget_contam_factor=" << cfg.tgt_contam;
    meta << "\ntarget_contam_factor_err=" << cfg.tgt_contam_err;
    meta << "\nsimc_yield_scale=" << cfg.simc_yield_scale;
    meta << "\nnormalize_mmiss=" << (cfg.normalize_mmiss ? "yes" : "no");
    if (cfg.normalize_mmiss) {
        meta << "\nyield_norm_scope=all_selected_inrange_Q2_xB_tprime_phi_events_in_mmiss_window";
        meta << "\nyield_norm_data=" << yield_norm_data;
        meta << "\nyield_norm_data_sumw2=" << yield_norm_data_sumw2;
        meta << "\nyield_norm_sim_reweighted_times_manual_scale=" << yield_norm_sim;
        meta << "\nyield_norm_sim_reweighted_times_manual_scale_sumw2=" << yield_norm_sim_sumw2;
        meta << "\nyield_norm_scale_data_over_sim=" << yield_norm_scale;
        meta << "\nyield_norm_scale_err_independent=" << yield_norm_scale_err;
        meta << "\nyield_norm_target_correction=replaced_with_unity";
        meta << "\nyield_norm_uncertainty=extraction_errors_conditional_on_scale_bootstrap_needed_for_correlation";
        meta << "\nyield_norm_absolute_caveat=data_area_match_removes_independent_absolute_SIMC_normalization";
    }
    meta << "\nsim_selection=is_exclusive_and_shared_reconstructed_mmiss_window";
    meta << "\nsim_weight=simc_yield_scale*(full_weight/sigcm)*model(vertex,p_best)";
    meta << "\nmodel_xsec_reporting=point_evaluation_at_physical_t_and_phi_bin_centers";
    meta << "\nmodel_reference=accepted_SIMC_immutable_base_weight_means_W_and_xB_shared_by_phi_bins";
    meta << "\nmodel_reference_t=physical_t_bin_center;Q2_ref=xB_ref*(W_ref^2-Mp^2)/(1-xB_ref)";
    meta << "\nmodel_evaluator=parameterized_SIMC_sig_param_2021_pi0";
    meta << "\nmodel_units=microbarn_per_MeV2_per_radian";
    meta << "\nmodel_identifier=" << nps_pi0_reweight::choose(cfg.model_identifier).id;
    meta << "\nmodel_provenance=" << nps_pi0_reweight::choose(cfg.model_identifier).source;
    meta << "\nmodel_domain=W_ge_2_GeV_MAID_mixture_not_implemented";
    meta << "\nmodel_fit_free_selection=" << cfg.model_free_parameters;
    meta << "\nmodel_fit_fixed_default=" << (cfg.fixed_default_model ? "yes" : "no");
    meta << "\nmodel_fit_max_iterations=" << cfg.model_max_iterations;
    meta << "\nmodel_fit_max_evaluations=" << cfg.model_max_evaluations;
    meta << "\nmodel_fit_tolerance=" << cfg.model_tolerance;
    meta << "\nreconstructed_binning=Q2_xB_physical_t_phi; tprime_is_phase_space_cut_only";
    meta << "\nmodel_fit_status=" << model_fit_status;
    meta << "\nmodel_fit_objective=sum_b_(Ydata_b-Ysim_b)^2/(data_sumw2_b+C^2*sum_event_weight_b^2)";
    meta << "\nmodel_fit_scope=one_parameter_vector_shared_by_all_selected_bins";
    meta << "\nmodel_fit_chi2_before=" << model_chi2_before;
    meta << "\nmodel_fit_chi2_after=" << model_chi2_after;
    meta << "\nmodel_fit_pvalue=" << model_fit_pvalue;
    meta << "\nmodel_fit_ndf=" << model_fit_ndf;
    meta << "\nmodel_vertex_epsilon="
         << (cutflow.n_sim_epsilon_fallback==0 ? "epsilon_i_all_selected_events" :
             "mixed_or_all_fixed_beam_approximation_from_Q2i_Wi");
    meta << "\nmodel_vertex_epsilon_fallback_count=" << cutflow.n_sim_epsilon_fallback;
    meta << "\nmodel_vertex_t=minus_ti;model_vertex_phi=phipqi_radians";
    meta << "\nmodel_rejected_bad_full_weight=" << cutflow.n_sim_bad_full_weight;
    meta << "\nmodel_rejected_bad_sigcm=" << cutflow.n_sim_bad_sigcm;
    meta << "\nmodel_rejected_bad_vertex=" << cutflow.n_sim_bad_vertex;
    meta << "\nmodel_rejected_bad_model=" << cutflow.n_sim_bad_model;
    meta << "\nmodel_conditional_uncertainties=data_sumw2_plus_final_MC_sumw2_no_independent_fit_parameter_error";
    meta << "\nphysical_t_edges:";
    for (double edge : t_edges) meta << " " << edge;
    meta << "\nmodel_sigcm_full_weight_mean_role=legacy_diagnostic_only";
    meta << "\nextraction=sigma_data=(corrected_data_yield/reweighted_simc_yield)*fitted_model_at_bin_centers";
    meta << "\nangular_convention=d2sigma_dtdphi=(U+sqrt(2epsilon(1+epsilon))*LT*cos(phi)+epsilon*TT*cos(2phi))/(2pi)";
    meta << "\npartons_projection=" << (cfg.partons_projection ? "yes" : "no");
    meta << "\ndiagnostics=" << (cfg.diagnostics ? "yes" : "no");
    meta << "\nglobal_1d_plot_normalization=SIMC_clone_scaled_to_target_corrected_data_integral_per_observable";
    meta << "\nglobal_1d_plot_normalization_role=display_only_not_used_in_extraction";
    meta << "\nmissing_mass_diagnostic=full_range_data_and_original_full_weight_exclusive_SIMC_shape_comparison_only";
    meta << "\nmissing_mass_diagnostic_scale=window_integral_data_over_SIMC_not_used_in_extraction";
    meta << "\nyield_support_diagnostic=data_yield_simc_yield_ratio_and_simc_effective_events_by_physical_t_phi";
    meta << "\nsigma_vs_t_plot=independent_panels_sigmaU_sigmaLT_sigmaTT_nb_per_GeV2";
    if (cfg.partons_projection) {
        meta << "\npartons_model=DVMPProcessGK06_DVMPCFFGK06_GPDGK19_LO";
        meta << "\npartons_mc_warmups=" << cfg.partons_warmups;
        meta << "\npartons_mc_calls=" << cfg.partons_calls;
        meta << "\npartons_kinematics=extraction_reference_Q2_xB_t";
        meta << "\npartons_role=theory_comparison_only_not_extraction_normalization";
    }
    TObjString(meta.str().c_str()).Write("analysis_metadata");
}

void ExclPi0XSecAnalysis::write_csv() {
    std::ofstream out(cfg.out_csv);
    out << std::setprecision(10);
    out << "it,iq,ix,ip,phi_lo,phi_hi,phi_center,q2_mean_data,xb_mean_data,tprime_mean_data,q2_mean_sim,xb_mean_sim,tprime_mean_sim,"
        << "epsilon,data,data_err,sim,sim_err,ratio,ratio_err,model_sigcm_full_weight_mean,xsec,xsec_err,xsec_sys_tgt,"
        << "model_xsec_phi_center,q2_ref,xb_ref,tprime_ref,t_ref,W_ref,t_lo,t_hi,t_center,Q2_direct_base_mean,ratio_before,bin_status,model_fit_status,xsec_units\n";
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
                        << s.epsilon << ","
                        << pb.data << "," << std::sqrt(std::max(0.0, pb.data_sumw2)) << ","
                        << pb.sim << "," << std::sqrt(std::max(0.0, pb.sim_sumw2)) << ","
                        << pb.ratio << "," << pb.ratio_err << ","
                        << calc_sigma_model_slice(pb) << ","
                        << pb.xsec << "," << pb.xsec_err << "," << pb.xsec_sys_tgt << ","
                        << pb.model_phi_center << "," << s.reference_model.q2 << ","
                        << s.reference_model.xb << "," << s.reference_model.tprime << ","
                        << s.reference_model.t << "," << s.reference_model.w << ","
                        << t_edges[it] << "," << t_edges[it+1] << ","
                        << 0.5*(t_edges[it]+t_edges[it+1]) << ","
                        << s.q2_direct_base_mean << "," << pb.ratio_before << ","
                        << (pb.n_data==0 ? "empty_data" :
                            pb.n_sim==0 || !(pb.sim>0) ? "empty_sim" :
                            !s.reference_supported ? "unsupported_reference" :
                            !std::isfinite(pb.model_phi_center) ?
                                "nonphysical_reporting_model" : "measured")
                        << "," << model_fit_status << ",ub/MeV2/rad\n";
                }
            }
        }
    }
}

void ExclPi0XSecAnalysis::write_slice_csv() {
    std::ofstream out(cfg.out_slice_csv);
    out << std::setprecision(10);
    out << "it,iq,ix,t_lo,t_hi,t_center,q2_lo,q2_hi,xb_lo,xb_hi,"
        << "mean_q2_data,mean_xb_data,mean_tprime_data,mean_q2_sim,mean_xb_sim,mean_t_sim,mean_tprime_sim,"
        << "epsilon,has_model_xsec,sumw_data,sumw_sim,"
        << "fit_ratio_ok,fit_ratio_chi2,fit_ratio_ndf,fit_ratio_A,fit_ratio_Aerr,fit_ratio_B,fit_ratio_Berr,fit_ratio_C,fit_ratio_Cerr,"
        << "fit_xsec_ok,fit_xsec_chi2,fit_xsec_ndf,fit_xsec_sigmaU,fit_xsec_sigmaUerr,fit_xsec_sigmaTL,fit_xsec_sigmaTLerr,fit_xsec_sigmaTT,fit_xsec_sigmaTTerr,"
        << "sigmaTLp_available,sigmaTLp_status,sigmaTLp,sigmaTLp_err,"
        << "partons_ok,partons_epsilon,partons_electron_flux_xbq2,partons_sigmaU,partons_sigmaLT,partons_sigmaTT\n";
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                double tlo = t_edges[it], thi = t_edges[it+1];
                double qlo = q2_edges[iq], qhi = q2_edges[iq+1];
                double xlo = xb_edges_by_q2[iq][ix], xhi = xb_edges_by_q2[iq][ix+1];
                double tcenter = 0.5 * (tlo + thi);
                out << it << "," << iq << "," << ix << ","
                    << tlo << "," << thi << "," << tcenter << ","
                    << qlo << "," << qhi << ","
                    << xlo << "," << xhi << ","
                    << s.mean_q2_data << "," << s.mean_xb_data << "," << s.mean_tprime_data << ","
                    << s.mean_q2_sim << "," << s.mean_xb_sim << "," << s.mean_t_sim << "," << s.mean_tprime_sim << ","
                    << s.epsilon << "," << (s.has_model_xsec ? 1 : 0) << ","
                    << s.sumw_data << "," << s.sumw_sim << ","
                    << (s.fit_ratio.ok ? 1 : 0) << "," << s.fit_ratio.chi2 << "," << s.fit_ratio.ndf << ","
                    << (s.fit_ratio.p.size() > 0 ? s.fit_ratio.p[0] : 0.0) << "," << (s.fit_ratio.perr.size() > 0 ? s.fit_ratio.perr[0] : 0.0) << ","
                    << (s.fit_ratio.p.size() > 1 ? s.fit_ratio.p[1] : 0.0) << "," << (s.fit_ratio.perr.size() > 1 ? s.fit_ratio.perr[1] : 0.0) << ","
                    << (s.fit_ratio.p.size() > 2 ? s.fit_ratio.p[2] : 0.0) << "," << (s.fit_ratio.perr.size() > 2 ? s.fit_ratio.perr[2] : 0.0) << ","
                    << (s.fit_xsec.ok ? 1 : 0) << "," << s.fit_xsec.chi2 << "," << s.fit_xsec.ndf << ","
                    << s.fit_xsec.sigmaU << "," << s.fit_xsec.sigmaU_err << ","
                    << s.fit_xsec.sigmaTL << "," << s.fit_xsec.sigmaTL_err << ","
                    << s.fit_xsec.sigmaTT << "," << s.fit_xsec.sigmaTT_err << ","
                    << 0 << "," << kSigmaTLpStatus << ","
                    << std::numeric_limits<double>::quiet_NaN() << ","
                    << std::numeric_limits<double>::quiet_NaN() << ","
                    << (s.partons_ok ? 1 : 0) << "," << s.partons_epsilon << ","
                    << s.partons_electron_flux_xbq2 << "," << s.partons_sigmaU << ","
                    << s.partons_sigmaLT << "," << s.partons_sigmaTT << "\n";
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
    validate_xsec_binning(cfg, true);
    apply_publication_style();

    auto resolve_relative_output = [this](std::string& name) {
        fs::path p(name);
        if (p.is_relative()) name = (fs::path(cfg.out_dir) / p).string();
    };

    resolve_relative_output(cfg.out_root);
    resolve_relative_output(cfg.out_csv);
    resolve_relative_output(cfg.out_slice_csv);

    fs::create_directories(cfg.out_dir);
    fs::create_directories(fs::path(cfg.out_dir) / "global");
    fs::create_directories(fs::path(cfg.out_dir) / "slices");
    fs::create_directories(fs::path(cfg.out_dir) / "sigma_vs_t");
    init_combined_pdf();
    load_input();
    detect_optional_branches();
    build_binning();
    init_storage();
    fill_from_trees();
    fit_model();
    apply_simc_to_data_yield_normalization();
    compute_mmiss_shape_scales();
    compute_ratios_and_xsec();
    fit_slices();
    compute_partons_projection();

    fout = TFile::Open(cfg.out_root.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) die("Cannot open output ROOT file.");

    make_global_plots();
    make_mmiss_comparison_plots();
    make_yield_diagnostic_plots();
    make_epsilon_plots();
    make_slice_plots();
    make_sigma_vs_t_plots();
    make_partons_projection_plots();
    close_combined_pdf();
    write_results();
    write_csv();
    write_slice_csv();
    // cleanup() is now handled by the destructor; do not call here to avoid double-free.

    log("Analysis complete.");
}

static std::string require_arg(int& i, int argc, char** argv, const std::string& option) {
    if (++i >= argc) throw std::runtime_error("Missing value for " + option);
    return argv[i];
}

static std::string safe_token(const std::string& value) {
    std::string out;
    for (char c : value)
        out += (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-') ? c : '_';
    return out;
}

static void print_usage(const char* program) {
    std::cout << "Usage: " << program << R"( [options]

Hall-C extraction in physical-t, Q2, xB, phi bins:
  Ysim(p) = simc_yield_scale * sum[(full_weight/sigcm)*model(vertex,p)]
  sigma_data = (Ydata/Ysim(p_best))*model(W_ref,xB_ref,t_center,phi_center,p_best)
  W_ref and xB_ref are base-weight means of accepted SIMC events per slice.
  Q2_ref is derived from W_ref and xB_ref; model units: ub/MeV^2/radian.
  For W<2 GeV the Fortran MAID branch is unsupported.

Inputs/outputs:
  --kin <name>                 Derive standard repository paths only; physics
                               defaults remain those in xsec_config.h
  --target <name>              Combined target token (default: LH2)
  --output-base <path>         Output base (default: output)
  --root-dir <path>            Directory containing input ROOT files
  --data-file <path>           Data ROOT file; tree "physics"
  --sim-file <path>            SIMC ROOT file; tree "simulation"
  --out-dir <path>             Output directory
  --out-root <path>            Output ROOT file
  --out-csv <path>             Per-phi CSV
  --out-slice-csv <path>       Fourier-coefficient CSV
  --all-plots-pdf <path>       Combined PDF

Selection/binning:
  --mmiss-lower <GeV>          Default: 0.80
  --mmiss-upper <GeV>          Default: 1.10
  --target-contam <factor>     Divide data yield by factor (default: 0.584)
  --target-contam-err <value>  Absolute factor uncertainty (default: 0.014)
  Bin edges are fixed in xsec_config.h (including physical t).
  --ebeam <GeV>

Behavior:
  --simc-yield-scale <factor>  Extra SIMC scale; default 1
  --normalize_mmiss            Historical yield rescaling; only with
                               --fixed-default-model. It erases absolute sensitivity.
  --normalize-simc-to-data     Clear-name alias for --normalize_mmiss
  --model <id>                 Model identifier (currently SIMC_sig_param_2021_pi0_W_ge_2)
  --fixed-default-model        Evaluate original SIMC model without fitting
  --model-free <names>         Comma-separated coefficient names (default:
                               plus.p5,plus.p7,plus.p9), shared across bins
  --model-max-iterations <n>   Minuit2 iteration limit
  --model-max-evaluations <n>  Minuit2 objective evaluation limit
  --model-tolerance <x>        Minuit2 convergence tolerance
  --partons                    Add native GK06/GPDGK19 comparison
  --partons-warmups <int>      PARTONS MC warm-up calls; default 10000
  --partons-calls <int>        PARTONS MC calls; default 100000
  --quiet --no-diagnostics --no-png --no-pdf --help

Reported responses:
  U, LT, and TT are extracted. LT' is unavailable because this method has no
  helicity-specific luminosities or beam-polarization input.

Required SIMC branches include full_weight, sigcm, and is_exclusive.
)";
}

static AnalysisConfig parse_config(int argc, char** argv) {
    AnalysisConfig cfg = make_simc_model_config();
    cfg.partons_executable_path = argv[0];
    std::string kin, target = "LH2", output_base = "output", root_dir;
    bool data_set = false, sim_set = false, out_dir_set = false;
    bool out_root_set = false, out_csv_set = false, out_slice_set = false, all_pdf_set = false;
    int positional = 0;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") { print_usage(argv[0]); std::exit(0); }
        if (arg == "--kin") { kin = require_arg(i, argc, argv, arg); continue; }
        if (arg == "--target") { target = require_arg(i, argc, argv, arg); continue; }
        if (arg == "--output-base") { output_base = require_arg(i, argc, argv, arg); continue; }
        if (arg == "--root-dir") { root_dir = require_arg(i, argc, argv, arg); continue; }
        if (arg == "--data-file") { cfg.data_file = require_arg(i, argc, argv, arg); data_set = true; continue; }
        if (arg == "--sim-file") { cfg.simc_file = require_arg(i, argc, argv, arg); sim_set = true; continue; }
        if (arg == "--out-dir") { cfg.out_dir = require_arg(i, argc, argv, arg); out_dir_set = true; continue; }
        if (arg == "--out-root") { cfg.out_root = require_arg(i, argc, argv, arg); out_root_set = true; continue; }
        if (arg == "--out-csv") { cfg.out_csv = require_arg(i, argc, argv, arg); out_csv_set = true; continue; }
        if (arg == "--out-slice-csv") { cfg.out_slice_csv = require_arg(i, argc, argv, arg); out_slice_set = true; continue; }
        if (arg == "--all-plots-pdf") { cfg.out_all_plots_pdf = require_arg(i, argc, argv, arg); all_pdf_set = true; continue; }
        if (arg == "--mmiss-lower") { cfg.mmiss_lower_gev = std::stod(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--mmiss-upper") { cfg.mmiss_upper_gev = std::stod(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--target-contam") { cfg.tgt_contam = std::stod(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--target-contam-err") { cfg.tgt_contam_err = std::stod(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--simc-yield-scale") { cfg.simc_yield_scale = std::stod(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--normalize_mmiss" || arg == "--normalize-mmiss" ||
            arg == "--normalize-simc-to-data") { cfg.normalize_mmiss = true; continue; }
        if (arg == "--partons") { cfg.partons_projection = true; continue; }
        if (arg == "--partons-warmups") { cfg.partons_warmups = std::stoi(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--partons-calls") { cfg.partons_calls = std::stoi(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--ebeam") { cfg.ebeam = std::stod(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--model") { cfg.model_identifier = require_arg(i, argc, argv, arg); continue; }
        if (arg == "--fixed-default-model") { cfg.fixed_default_model = true; continue; }
        if (arg == "--model-free") { cfg.model_free_parameters = require_arg(i, argc, argv, arg); continue; }
        if (arg == "--model-max-iterations") { cfg.model_max_iterations = std::stoi(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--model-max-evaluations") { cfg.model_max_evaluations = std::stoi(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--model-tolerance") { cfg.model_tolerance = std::stod(require_arg(i, argc, argv, arg)); continue; }
        if (arg == "--quiet") { cfg.verbose = false; continue; }
        if (arg == "--no-diagnostics") { cfg.diagnostics = false; continue; }
        if (arg == "--no-png") { cfg.write_png = false; continue; }
        if (arg == "--no-pdf") { cfg.write_pdf = false; continue; }
        if (!arg.empty() && arg[0] == '-') throw std::runtime_error("Unknown option: " + arg);
        if (positional == 0) { cfg.out_root = arg; out_root_set = true; }
        else if (positional == 1) { cfg.out_csv = arg; out_csv_set = true; }
        else throw std::runtime_error("Unexpected positional argument: " + arg);
        ++positional;
    }

    if (!kin.empty() && kin != cfg.configured_kinematic) {
        std::cerr << "[WARN] --kin " << kin << " changes paths only; physics defaults remain "
                  << cfg.configured_kinematic
                  << " from xsec_config.h unless explicitly overridden.\n";
    }
    if (!kin.empty() && root_dir.empty()) root_dir = (fs::path(output_base) / safe_token(kin) / "root").string();
    if (!root_dir.empty()) {
        if (!data_set) cfg.data_file = (fs::path(root_dir) / ("combined_branches_" + safe_token(target) + ".root")).string();
        if (!sim_set) cfg.simc_file = (fs::path(root_dir) / "simc_pi0_analysis_output_smeared.root").string();
        if (!out_dir_set) cfg.out_dir = (fs::path(root_dir).parent_path() / "xsec_simc_model").string();
    }
    if (!out_root_set) cfg.out_root = "excl_xsec_pi0_analysis_simc_model_output.root";
    if (!out_csv_set) cfg.out_csv = "excl_xsec_pi0_analysis_simc_model_summary.csv";
    if (!out_slice_set) cfg.out_slice_csv = "excl_xsec_pi0_analysis_simc_model_slice_summary.csv";
    if (!all_pdf_set) cfg.out_all_plots_pdf = "all_generated_plots_simc_model.pdf";

    validate_xsec_binning(cfg, true);
    if (cfg.model_max_iterations<=0 || cfg.model_max_evaluations<=0 ||
        !(std::isfinite(cfg.model_tolerance) && cfg.model_tolerance>0))
        throw std::runtime_error("Invalid model optimization controls.");
    if (cfg.normalize_mmiss && !cfg.fixed_default_model)
        throw std::runtime_error("--normalize_mmiss conflicts with iterative absolute model fitting; use display-only plot normalization.");
    if (!(cfg.mmiss_lower_gev < cfg.mmiss_upper_gev)) throw std::runtime_error("Invalid missing-mass window.");
    if (!(cfg.tgt_contam > 0.0) || cfg.tgt_contam_err < 0.0)
        throw std::runtime_error("Target factor must be positive and its uncertainty nonnegative.");
    if (!(cfg.simc_yield_scale > 0.0)) throw std::runtime_error("SIMC yield scale must be positive.");
    if (cfg.normalize_mmiss) {
        // Historical normalization mode: the data integral sets the global
        // SIMC scale, so applying the separate 0.584 divisor would double-count.
        cfg.tgt_contam = 1.0;
        cfg.tgt_contam_err = 0.0;
    }
    if (cfg.partons_warmups <= 0 || cfg.partons_calls <= 0)
        throw std::runtime_error("PARTONS warmups and calls must be positive.");
#ifndef NPS_ENABLE_PARTONS
    if (cfg.partons_projection)
        throw std::runtime_error("--partons requires a native PARTONS-linked executable.");
#endif
    return cfg;
}

int main(int argc, char** argv) {
    try {
        gROOT->SetBatch(kTRUE);
        TH1::SetDefaultSumw2(kTRUE);
        apply_publication_style();
        AnalysisConfig cfg = parse_config(argc, argv);
        ExclPi0XSecAnalysis ana(cfg);
        ana.Run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[FATAL] " << e.what() << "\n";
        return 1;
    }
}
