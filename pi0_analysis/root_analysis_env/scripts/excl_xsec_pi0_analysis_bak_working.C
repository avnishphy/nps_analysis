// g++ -O2 -std=c++17 -o excl_xsec scripts/excl_xsec_pi0_analysis.C `root-config --cflags --libs`

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
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TDecompLU.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TMath.h>

#include <algorithm>
#include <array>
#include <cmath>
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
#include <vector>

namespace fs = std::filesystem;

struct AnalysisConfig {
    std::string simc_file = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/simc_pi0_analysis_output_smeared.root";
    std::string data_file = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root";
    std::string simc_tree = "simulation";
    std::string data_tree = "physics";

    std::string out_root = "excl_xsec_pi0_analysis_output.root";
    std::string out_csv  = "excl_xsec_pi0_analysis_summary.csv";
    std::string out_slice_csv = "excl_xsec_pi0_analysis_slice_summary.csv";
    std::string out_dir  = "output_pi0_xsec";

    double ebeam = 10.604; // GeV, set from Hall C beam energy used in the run period
    double mp    = 0.9382720813;
    double mpi0  = 0.1349768;

    int n_phi    = 12;
    int n_tprime = 5;
    int n_q2     = 5;
    int n_xb     = 5;

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

    // If the SIMC file contains an event-level model cross-section branch,
    // the code will use it to convert data/MC ratios into absolute cross sections.
    // Candidate branch names are searched in this order.
    std::vector<std::string> model_xsec_candidates = {
        "sigcm"
    };
};

struct PhiBin {
    double data = 0.0, data_sumw2 = 0.0;
    double sim  = 0.0, sim_sumw2  = 0.0;
    double model_wsum = 0.0, model_xsec_wsum = 0.0;
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
    // ...existing code...
    // ...existing code...
    // Only one, public destructor should exist. If you see this comment, the destructor is now public and unique.
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
    bool has_model_xsec = false;
    std::string model_xsec_branch;

    std::vector<double> phi_edges, tprime_edges, q2_edges, xb_edges;
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
    void make_slice_plots();
    void write_results();
    void write_csv();
    void write_slice_csv();
    void cleanup();

    void fill_data_event(double q2, double t, double tmin, double xb, double phi, double pi0_weight, float scale, double charge_uC, int is_exclusive, int helicity, bool use_helicity);
    void fill_sim_event(float q2, float t, float tmin, float xb, float phi, float full_weight, int is_exclusive, float model_xsec);
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
}

void ExclPi0XSecAnalysis::detect_optional_branches() {
    has_helicity = (t_data->GetBranch("helicity") != nullptr);

    for (const auto& b : cfg.model_xsec_candidates) {
        if (t_sim->GetBranch(b.c_str())) {
            has_model_xsec = true;
            model_xsec_branch = b;
            break;
        }
    }

    if (!has_helicity) warn("No helicity branch found in data tree. TL' extraction will be reported only as a diagnostic asymmetry.");
    if (!has_model_xsec) warn("No model cross-section branch found in SIMC. Absolute cross sections will not be formed; ratio plots will still be produced.");
}

void ExclPi0XSecAnalysis::build_binning() {
    // First pass: collect accepted values from the combined data+SIMC samples
    // inside the requested kinematic windows, then use those to define quantile bins.
    std::vector<double> tprime_vals, q2_vals, xb_vals;
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
                q2_vals.push_back(q2);
                xb_vals.push_back(xb);
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
            if (has_model_xsec) t->SetBranchAddress(model_xsec_branch.c_str(), &model_xsec);
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
                q2_vals.push_back(q2);
                xb_vals.push_back(xb);
            }
        }
    };

    collect(t_data, false);
    collect(t_sim, true);

    phi_edges = std::vector<double>(cfg.n_phi + 1, 0.0);
    for (int i = 0; i <= cfg.n_phi; ++i) phi_edges[i] = cfg.phi_min + (cfg.phi_max - cfg.phi_min) * double(i) / double(cfg.n_phi);

    tprime_edges = quantile_edges(tprime_vals, cfg.n_tprime, cfg.tprime_min, cfg.tprime_max);
    q2_edges     = quantile_edges(q2_vals,     cfg.n_q2,     cfg.q2_min,     cfg.q2_max);
    xb_edges     = quantile_edges(xb_vals,     cfg.n_xb,     cfg.xb_min,     cfg.xb_max);

    log("Computed bin edges.");
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

void ExclPi0XSecAnalysis::fill_data_event(double q2, double t, double tmin, double xb, double phi, double pi0_weight, float scale, double charge_uC, int is_exclusive, int helicity, bool use_helicity) {
    cutflow.n_data_total++;
    if (!is_exclusive) return;
    double tprime = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tprime)) return;
    if (!std::isfinite(phi)) return;

    cutflow.n_data_pass++;

    int it = find_bin(tprime_edges, tprime, false);
    int iq = find_bin(q2_edges, q2, false);
    int ix = find_bin(xb_edges, xb, false);
    int ip = find_bin(phi_edges, phi, true);
    if (it < 0 || iq < 0 || ix < 0 || ip < 0) return;
    cutflow.n_data_inrange++;

    double w = 1000.0 * is_exclusive * pi0_weight * static_cast<double>(scale) / charge_uC;
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

void ExclPi0XSecAnalysis::fill_sim_event(float q2, float t, float tmin, float xb, float phi, float full_weight, int is_exclusive, float model_xsec) {
    cutflow.n_sim_total++;
    if (!is_exclusive) return;
    double tprime = calc_tprime(t, tmin);
    if (!slice_passes_kin(q2, xb, tprime)) return;
    if (!std::isfinite(phi)) return;

    cutflow.n_sim_pass++;

    int it = find_bin(tprime_edges, tprime, false);
    int iq = find_bin(q2_edges, q2, false);
    int ix = find_bin(xb_edges, xb, false);
    int ip = find_bin(phi_edges, phi, true);
    if (it < 0 || iq < 0 || ix < 0 || ip < 0) return;
    cutflow.n_sim_inrange++;

    double w = static_cast<double>(full_weight);
    PhiBin& pb = slice(it, iq, ix).phi[ip];
    pb.sim += w;
    pb.sim_sumw2 += w * w;
    pb.model_wsum += w;
    if (has_model_xsec) pb.model_xsec_wsum += w * static_cast<double>(model_xsec);
    pb.n_sim += 1;
    pb.mean_q2_sim += q2 * w;
    pb.mean_xb_sim += xb * w;
    pb.mean_tprime_sim += tprime * w;

    SliceResult& s = slice(it, iq, ix);
    s.sumw_sim += w;
    s.sumw2_sim += w * w;
    s.mean_q2_sim += q2 * w;
    s.mean_xb_sim += xb * w;
    s.mean_tprime_sim += tprime * w;
    accumulate_global_histograms(q2, xb, tprime, wrap_phi(phi), 0.0, w);
}

void ExclPi0XSecAnalysis::fill_from_trees() {
    // Data loop
    double dq2 = 0, dt = 0, dtmin = 0, dxb = 0, dphi = 0, dpi0_weight = 0;
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
    if (has_helicity) t_data->SetBranchAddress("helicity", &dhelicity);

    // First pass: sum charge_uC per run_number
    // std::map<int, double> run_charge_map;
    const Long64_t ndata = t_data->GetEntries();
    std::unordered_set<int> seen_runs;
    double total_charge_uC = 0.0;

    for (Long64_t i = 0; i < ndata; ++i) {
        t_data->GetEntry(i);
        if (seen_runs.insert(drun_number).second) {
            total_charge_uC += dcharge_uC;
        }
    }

    // Second pass: fill events using total_charge_uC
    for (Long64_t i = 0; i < ndata; ++i) {
        t_data->GetEntry(i);
        fill_data_event(dq2, dt, dtmin, dxb, dphi, dpi0_weight, dscale, total_charge_uC, dis_exclusive, dhelicity, has_helicity);
    }

    // SIMC loop
    float sim_q2 = 0, sim_t = 0, sim_tmin = 0, sim_xb = 0, sim_phi = 0, full_weight = 0, model_xsec = 0;
    int sim_is_exclusive = 0;
    t_sim->SetBranchAddress("Q2", &sim_q2);
    t_sim->SetBranchAddress("t", &sim_t);
    t_sim->SetBranchAddress("tmin", &sim_tmin);
    t_sim->SetBranchAddress("xB", &sim_xb);
    t_sim->SetBranchAddress("phi", &sim_phi);
    t_sim->SetBranchAddress("full_weight", &full_weight);
    t_sim->SetBranchAddress("is_exclusive", &sim_is_exclusive);
    if (has_model_xsec) t_sim->SetBranchAddress(model_xsec_branch.c_str(), &model_xsec);

    const Long64_t nsim = t_sim->GetEntries();
    for (Long64_t i = 0; i < nsim; ++i) {
        t_sim->GetEntry(i);
        fill_sim_event(sim_q2, sim_t, sim_tmin, sim_xb, sim_phi, full_weight, sim_is_exclusive, model_xsec);
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
            if (pb.sim > 0.0) {
                pb.mean_q2_sim /= pb.sim;
                pb.mean_xb_sim /= pb.sim;
                pb.mean_tprime_sim /= pb.sim;
            }
        }
    }
}

void ExclPi0XSecAnalysis::compute_ratios_and_xsec() {
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                SliceResult& s = slice(it, iq, ix);
                s.epsilon = epsilon_virtual(cfg.ebeam, std::max(1e-12, s.mean_q2_sim), std::max(1e-12, s.mean_xb_sim), cfg.mp);
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    PhiBin& pb = s.phi[ip];
                    pb.data /= cfg.tgt_contam;
                    pb.data_sumw2 /= (cfg.tgt_contam * cfg.tgt_contam);
                    pb.data_plus /= cfg.tgt_contam;
                    pb.data_minus /= cfg.tgt_contam;
                    pb.data_plus_sumw2 /= (cfg.tgt_contam * cfg.tgt_contam);
                    pb.data_minus_sumw2 /= (cfg.tgt_contam * cfg.tgt_contam);

                    if (pb.sim > 0.0) {
                        pb.ratio = pb.data / pb.sim;
                        pb.ratio_err = std::sqrt((pb.data_sumw2 / (pb.sim * pb.sim)) +
                                                 (pb.data * pb.data * pb.sim_sumw2 / (pb.sim * pb.sim * pb.sim * pb.sim)));
                    }

                    pb.xsec_sys_tgt = pb.xsec * (cfg.tgt_contam_err / cfg.tgt_contam);

                    if (s.has_model_xsec && pb.sim > 0.0 && pb.model_wsum > 0.0) {
                        double sigma_model = pb.model_xsec_wsum / pb.model_wsum;
                        pb.xsec = pb.ratio * sigma_model;
                        pb.xsec_err = pb.ratio_err * sigma_model;
                        pb.xsec_sys_tgt = pb.xsec * (cfg.tgt_contam_err / cfg.tgt_contam);
                        pb.mean_q2_xsec = pb.mean_q2_sim;
                        pb.mean_xb_xsec = pb.mean_xb_sim;
                        pb.mean_tprime_xsec = pb.mean_tprime_sim;
                    }
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
                        double eps = epsilon_virtual(cfg.ebeam, std::max(1e-12, s.mean_q2_sim), std::max(1e-12, s.mean_xb_sim), cfg.mp);
                        s.fit_xsec.sigmaU = s.fit_xsec.p[0];
                        s.fit_xsec.sigmaU_err = s.fit_xsec.perr[0];
                        s.fit_xsec.sigmaTL  = s.fit_xsec.p[1] / std::sqrt(std::max(1e-12, 2.0 * eps * (1.0 + eps)));
                        s.fit_xsec.sigmaTL_err = s.fit_xsec.perr[1] / std::sqrt(std::max(1e-12, 2.0 * eps * (1.0 + eps)));
                        s.fit_xsec.sigmaTT  = s.fit_xsec.p[2] / std::max(1e-12, eps);
                        s.fit_xsec.sigmaTT_err = s.fit_xsec.perr[2] / std::max(1e-12, eps);
                    }
                }

                // Helicity diagnostic: build asymmetry A(phi) = (Y+ - Y-) / (Y+ + Y-)
                // This is saved as a diagnostic if helicity is available.
                if (has_helicity && s.has_model_xsec) {
                    std::vector<std::vector<double>> Xdiff;
                    std::vector<double> ydiff, eydiff;
                    double eps = std::max(1e-12, s.epsilon);
                    double kfac = std::sqrt(std::max(1e-12, 2.0 * eps * (1.0 - eps)));
                    for (int ip = 0; ip < cfg.n_phi; ++ip) {
                        const PhiBin& pb = s.phi[ip];
                        if (pb.sim <= 0.0) continue;
                        double sigma_model = pb.model_xsec_wsum / pb.model_wsum;
                        double diff_y = pb.data_plus - pb.data_minus;
                        double diff_yerr = std::sqrt(std::max(0.0, pb.data_plus_sumw2 + pb.data_minus_sumw2));
                        double diff_xsec = (diff_y / pb.sim) * sigma_model;
                        double diff_xsec_err = (diff_yerr / pb.sim) * sigma_model;
                        double sinbar = phi_basis_means(phi_edges[ip], phi_edges[ip+1])[3];
                        Xdiff.push_back({sinbar});
                        ydiff.push_back(diff_xsec);
                        eydiff.push_back(std::max(diff_xsec_err, 1e-12));
                    }
                    if (Xdiff.size() >= 1) {
                        s.fit_asym = weighted_linear_fit(Xdiff, ydiff, eydiff);
                        if (s.fit_asym.ok) {
                            s.fit_asym.helicity_asymmetry_fit = true;
                            s.fit_asym.sigmaTLp = s.fit_asym.p[0] / std::max(1e-12, 2.0 * kfac);
                            s.fit_asym.sigmaTLp_err = s.fit_asym.perr[0] / std::max(1e-12, 2.0 * kfac);
                        }
                    }
                } else if (has_helicity && !s.has_model_xsec) {
                    std::vector<std::vector<double>> Xasym;
                    std::vector<double> yasym, eyasym;
                    for (int ip = 0; ip < cfg.n_phi; ++ip) {
                        const PhiBin& pb = s.phi[ip];
                        double sum = pb.data_plus + pb.data_minus;
                        if (sum <= 0.0) continue;
                        double diff = pb.data_plus - pb.data_minus;
                        double asym = diff / sum;
                        double err = std::sqrt((pb.data_plus_sumw2 + pb.data_minus_sumw2)) / std::max(1e-12, sum);
                        Xasym.push_back({phi_basis_means(phi_edges[ip], phi_edges[ip+1])[3]});
                        yasym.push_back(asym);
                        eyasym.push_back(std::max(err, 1e-12));
                    }
                    if (Xasym.size() >= 1) {
                        s.fit_asym = weighted_linear_fit(Xasym, yasym, eyasym);
                        if (s.fit_asym.ok) s.fit_asym.helicity_asymmetry_fit = true;
                    }
                }
            }
        }
    }
}

void ExclPi0XSecAnalysis::write_canvas_pdf_png(TCanvas* c, const std::string& base) {
    if (cfg.write_pdf) c->SaveAs((base + ".pdf").c_str());
    if (cfg.write_png) c->SaveAs((base + ".png").c_str());
}

void ExclPi0XSecAnalysis::make_global_plots() {
    fs::create_directories(fs::path(cfg.out_dir) / "global");
    TCanvas c1("c1", "global", 1200, 900);
    c1.Divide(2,2);

    c1.cd(1);
    h_q2_data->SetLineWidth(2); h_q2_data->SetLineColor(kBlack);
    h_q2_sim->SetLineWidth(2); h_q2_sim->SetLineColor(kRed);
    h_q2_data->Draw("hist");
    h_q2_sim->Draw("hist same");
    auto leg1 = new TLegend(0.62,0.75,0.88,0.88); leg1->AddEntry(h_q2_data.get(),"Data","l"); leg1->AddEntry(h_q2_sim.get(),"SIMC","l"); leg1->Draw();

    c1.cd(2);
    h_xb_data->SetLineWidth(2); h_xb_data->SetLineColor(kBlack);
    h_xb_sim->SetLineWidth(2); h_xb_sim->SetLineColor(kRed);
    h_xb_data->Draw("hist");
    h_xb_sim->Draw("hist same");
    auto leg2 = new TLegend(0.62,0.75,0.88,0.88); leg2->AddEntry(h_xb_data.get(),"Data","l"); leg2->AddEntry(h_xb_sim.get(),"SIMC","l"); leg2->Draw();

    c1.cd(3);
    h_tprime_data->SetLineWidth(2); h_tprime_data->SetLineColor(kBlack);
    h_tprime_sim->SetLineWidth(2); h_tprime_sim->SetLineColor(kRed);
    h_tprime_data->Draw("hist");
    h_tprime_sim->Draw("hist same");
    auto leg3 = new TLegend(0.62,0.75,0.88,0.88); leg3->AddEntry(h_tprime_data.get(),"Data","l"); leg3->AddEntry(h_tprime_sim.get(),"SIMC","l"); leg3->Draw();

    c1.cd(4);
    h_phi_data->SetLineWidth(2); h_phi_data->SetLineColor(kBlack);
    h_phi_sim->SetLineWidth(2); h_phi_sim->SetLineColor(kRed);
    h_phi_data->Draw("hist");
    h_phi_sim->Draw("hist same");
    auto leg4 = new TLegend(0.62,0.75,0.88,0.88); leg4->AddEntry(h_phi_data.get(),"Data","l"); leg4->AddEntry(h_phi_sim.get(),"SIMC","l"); leg4->Draw();

    c1.Update();
    write_canvas_pdf_png(&c1, (fs::path(cfg.out_dir) / "global" / "global_1d_distributions").string());

    TCanvas c2("c2", "2d", 1200, 500);
    c2.Divide(2,1);
    c2.cd(1); h_q2_xb_data->Draw("COLZ");
    c2.cd(2); h_tprime_phi_data->Draw("COLZ");
    c2.Update();
    write_canvas_pdf_png(&c2, (fs::path(cfg.out_dir) / "global" / "data_occupancy_2d").string());

    TCanvas c3("c3", "2d_sim", 1200, 500);
    c3.Divide(2,1);
    c3.cd(1); h_q2_xb_sim->Draw("COLZ");
    c3.cd(2); h_tprime_phi_sim->Draw("COLZ");
    c3.Update();
    write_canvas_pdf_png(&c3, (fs::path(cfg.out_dir) / "global" / "simc_occupancy_2d").string());
}

void ExclPi0XSecAnalysis::make_slice_plots() {
    fs::create_directories(fs::path(cfg.out_dir) / "slices");
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                std::ostringstream tag;
                tag << "t" << it << "_q" << iq << "_x" << ix;
                TCanvas c(("c_"+tag.str()).c_str(), tag.str().c_str(), 1400, 900);
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
                auto gdata = new TGraphErrors(cfg.n_phi, phi_c.data(), ydata.data(), nullptr, yerr.data());
                auto gsim  = new TGraphErrors(cfg.n_phi, phi_c.data(), ysim.data(), nullptr, ysimerr.data());
                gdata->SetTitle("Data/SIMC yields;#phi [rad];Weighted counts");
                gdata->SetMarkerStyle(20); gdata->SetLineColor(kBlack); gdata->SetMarkerColor(kBlack);
                gsim->SetMarkerStyle(24); gsim->SetLineColor(kRed); gsim->SetMarkerColor(kRed);
                gdata->Draw("AP");
                gsim->Draw("P SAME");
                auto l1 = new TLegend(0.63,0.74,0.88,0.88); l1->AddEntry(gdata,"Data","p"); l1->AddEntry(gsim,"SIMC","p"); l1->Draw();

                c.cd(2);
                auto grat = new TGraphErrors(cfg.n_phi, phi_c.data(), ratio.data(), nullptr, ratioerr.data());
                grat->SetTitle("Ratio data/SIMC;#phi [rad];Ratio");
                grat->SetMarkerStyle(20); grat->Draw("AP");
                if (s.fit_ratio.ok) {
                    auto f = new TF1(("fr_"+tag.str()).c_str(), "[0] + [1]*cos(x) + [2]*cos(2*x)", cfg.phi_min, cfg.phi_max);
                    f->SetParameters(s.fit_ratio.p[0], s.fit_ratio.p[1], s.fit_ratio.p[2]);
                    f->SetLineColor(kRed);
                    f->Draw("same");
                }

                c.cd(3);
                if (s.has_model_xsec) {
                    auto gx = new TGraphErrors(cfg.n_phi, phi_c.data(), xsec.data(), nullptr, xsecerr.data());
                    gx->SetTitle("Extracted #sigma;#phi [rad];Cross section");
                    gx->SetMarkerStyle(20); gx->Draw("AP");
                    if (s.fit_xsec.ok) {
                        auto f = new TF1(("fx_"+tag.str()).c_str(), "[0] + [1]*cos(x) + [2]*cos(2*x)", cfg.phi_min, cfg.phi_max);
                        f->SetParameters(s.fit_xsec.p[0], s.fit_xsec.p[1], s.fit_xsec.p[2]);
                        f->SetLineColor(kBlue);
                        f->Draw("same");
                    }
                } else {
                    auto g2 = new TGraphErrors(cfg.n_phi, phi_c.data(), ratio.data(), nullptr, ratioerr.data());
                    g2->SetTitle("Ratio only (no model xsec branch);#phi [rad];Ratio");
                    g2->SetMarkerStyle(20); g2->Draw("AP");
                }

                c.cd(4);
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
                    ga->SetMarkerStyle(20); ga->Draw("AP");
                } else {
                    auto h = new TH1D(("hnd_"+tag.str()).c_str(), "No helicity branch;TL' not extracted", 10, 0, 1);
                    h->Draw();
                }

                c.Update();
                write_canvas_pdf_png(&c, (fs::path(cfg.out_dir) / "slices" / ("slice_" + tag.str())).string());
            }
        }
    }
}

void ExclPi0XSecAnalysis::write_results() {
    fout->cd();
    // Save bin edges
    TVectorD vphi(phi_edges.size()), vt(tprime_edges.size()), vq(q2_edges.size()), vx(xb_edges.size());
    for (size_t i = 0; i < phi_edges.size(); ++i) vphi[i] = phi_edges[i];
    for (size_t i = 0; i < tprime_edges.size(); ++i) vt[i] = tprime_edges[i];
    for (size_t i = 0; i < q2_edges.size(); ++i) vq[i] = q2_edges[i];
    for (size_t i = 0; i < xb_edges.size(); ++i) vx[i] = xb_edges[i];
    vphi.Write("phi_edges");
    vt.Write("tprime_edges");
    vq.Write("q2_edges");
    vx.Write("xb_edges");

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
    meta << "\nphi_edges:";
    for (double x : phi_edges) meta << " " << x;
    meta << "\nmodel_xsec_branch=" << (has_model_xsec ? model_xsec_branch : "NONE");
    meta << "\nhelicity=" << (has_helicity ? "yes" : "no");
    TObjString(meta.str().c_str()).Write("analysis_metadata");
}

void ExclPi0XSecAnalysis::write_csv() {
    std::ofstream out(cfg.out_csv);
    out << std::setprecision(10);
    out << "it,iq,ix,ip,phi_lo,phi_hi,phi_center,q2_mean_data,xb_mean_data,tprime_mean_data,q2_mean_sim,xb_mean_sim,tprime_mean_sim,"
        << "data,data_err,sim,sim_err,ratio,ratio_err,xsec,xsec_err,xsec_sys_tgt\n";
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
        << "epsilon,has_model_xsec,sumw_data,sumw_sim,"
        << "fit_ratio_ok,fit_ratio_chi2,fit_ratio_ndf,fit_ratio_A,fit_ratio_Aerr,fit_ratio_B,fit_ratio_Berr,fit_ratio_C,fit_ratio_Cerr,"
        << "fit_xsec_ok,fit_xsec_chi2,fit_xsec_ndf,fit_xsec_sigmaU,fit_xsec_sigmaUerr,fit_xsec_sigmaTL,fit_xsec_sigmaTLerr,fit_xsec_sigmaTT,fit_xsec_sigmaTTerr,"
        << "fit_asym_ok,fit_asym_chi2,fit_asym_ndf,fit_asym_sigmaTLp,fit_asym_sigmaTLperr\n";
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                double tlo = tprime_edges[it], thi = tprime_edges[it+1];
                double qlo = q2_edges[iq], qhi = q2_edges[iq+1];
                double xlo = xb_edges[ix], xhi = xb_edges[ix+1];
                double tcenter = 0.5 * (tlo + thi);
                out << it << "," << iq << "," << ix << ","
                    << tlo << "," << thi << "," << tcenter << ","
                    << qlo << "," << qhi << ","
                    << xlo << "," << xhi << ","
                    << s.mean_q2_data << "," << s.mean_xb_data << "," << s.mean_tprime_data << ","
                    << s.mean_q2_sim << "," << s.mean_xb_sim << "," << s.mean_tprime_sim << ","
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
                    << (s.fit_asym.ok ? 1 : 0) << "," << s.fit_asym.chi2 << "," << s.fit_asym.ndf << ","
                    << s.fit_asym.sigmaTLp << "," << s.fit_asym.sigmaTLp_err << "\n";
            }
        }
    }
}

void ExclPi0XSecAnalysis::cleanup() {

    if (fout) { fout->Write(); fout->Close(); fout = nullptr; }
    if (f_sim) { f_sim->Close(); f_sim = nullptr; }
    if (f_data) { f_data->Close(); f_data = nullptr; }
}

void ExclPi0XSecAnalysis::Run() {
    fs::create_directories(cfg.out_dir);
    fs::create_directories(fs::path(cfg.out_dir) / "global");
    fs::create_directories(fs::path(cfg.out_dir) / "slices");
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
    make_slice_plots();
    write_results();
    write_csv();
    write_slice_csv();
    // cleanup() is now handled by the destructor; do not call here to avoid double-free.

    log("Analysis complete.");
}

int main(int argc, char** argv) {
    try {
        TH1::SetDefaultSumw2(kTRUE);
        gStyle->SetOptStat(0);
        AnalysisConfig cfg;
        if (argc > 1) cfg.out_root = argv[1];
        if (argc > 2) cfg.out_csv = argv[2];
        ExclPi0XSecAnalysis ana(cfg);
        ana.Run();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "[FATAL] " << e.what() << "\n";
        return 1;
    }
}
