#pragma once

// Owned result records. Reconstructed observations and generated-bin parameters are distinct concepts.
#include "xsec_config.h"

inline std::string xsec_csv_quote(const std::string& value) {
    std::string quoted = "\"";
    for (char c : value) { if (c == '"') quoted += '"'; quoted += c; }
    return quoted + '"';
}

// Fixed upstream pi0_weight and run corrections define these event weights.
struct XsecWeightMoments {
    long long n = 0;
    double sum = 0., sum2 = 0.;
    double min = std::numeric_limits<double>::infinity();
    double max = -std::numeric_limits<double>::infinity();
    double max2 = 0.;
    void add(double w) {
        ++n; sum += w; sum2 += w*w;
        min = std::min(min,w); max = std::max(max,w);
        max2 = std::max(max2,w*w);
    }
    void merge(const XsecWeightMoments& o) {
        n += o.n; sum += o.sum; sum2 += o.sum2;
        min = std::min(min,o.min); max = std::max(max,o.max);
        max2 = std::max(max2,o.max2);
    }
    void divide(double f) {
        sum /= f; sum2 /= f*f;
        if (n) { min /= f; max /= f; max2 /= f*f; }
    }
    double scale() const { return sum > 0. ? sum2/sum : std::numeric_limits<double>::quiet_NaN(); }
    double effective_n() const { return sum2 > 0. ? sum*sum/sum2 : 0.; }
};

struct XsecScaledRow {
    double s_observed = std::numeric_limits<double>::quiet_NaN();
    double s_used = std::numeric_limits<double>::quiet_NaN();
    double deviance = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    std::string scale_source, exclusion_reason;
    bool support = false, included = false;
};

struct PhiBin {
    double data = 0.0, data_sumw2 = 0.0;
    // After a successful fit, sim is the reconstructed forward prediction and
    // sim_sumw2 its conditional same-fit variance, not a raw MC sum(w^2).
    double sim  = 0.0, sim_sumw2  = 0.0;

    // Reconstructed response weight for QA means only. Fourier integrals and
    // finite-MC covariance live in the reconstructed-by-generated matrix.
    double sim_base = 0.0;

    double data_plus = 0.0, data_plus_sumw2 = 0.0;
    double data_minus = 0.0, data_minus_sumw2 = 0.0;
    double ratio = 0.0, ratio_err = 0.0;
    double xsec = 0.0, xsec_err = 0.0;
    double xsec_sys_tgt = 0.0;
    double mean_q2_data = 0.0, mean_xb_data = 0.0, mean_tprime_data = 0.0;
    double mean_q2_sim = 0.0, mean_xb_sim = 0.0, mean_tprime_sim = 0.0;
    double mean_q2_xsec = 0.0, mean_xb_xsec = 0.0, mean_tprime_xsec = 0.0;
    int n_data = 0, n_sim = 0;
    XsecWeightMoments weights;
};

// One Eq. 5.30-5.31 point for the corresponding reconstructed/truth phi cell.
// The fitted curve in PhiBin::xsec keeps its historical meaning.
struct ExperimentalPoint {
    int row = -1, truth_cell = -1;
    std::string status = "unavailable_fit";
    double contribution = std::numeric_limits<double>::quiet_NaN();
    double subtracted_yield = std::numeric_limits<double>::quiet_NaN();
    double correction = std::numeric_limits<double>::quiet_NaN();
    double sigma_reference = std::numeric_limits<double>::quiet_NaN();
    double sigma_exp = std::numeric_limits<double>::quiet_NaN();
    double sigma_exp_err = std::numeric_limits<double>::quiet_NaN();
    double epsilon_reference = std::numeric_limits<double>::quiet_NaN();
    double q2_reference = std::numeric_limits<double>::quiet_NaN();
    double xb_reference = std::numeric_limits<double>::quiet_NaN();
    double tprime_reference = std::numeric_limits<double>::quiet_NaN();
    double reco_prediction_q2 = std::numeric_limits<double>::quiet_NaN();
    double reco_prediction_xb = std::numeric_limits<double>::quiet_NaN();
    double reco_prediction_tprime = std::numeric_limits<double>::quiet_NaN();
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
    std::string fit_scope = "global_migration", fit_failure_reason;
    int fit_rank = 0, fit_mc_iterations = 0;
    double fit_condition = std::numeric_limits<double>::quiet_NaN();
    std::vector<PhiBin> phi;
    double sumw_data = 0.0, sumw2_data = 0.0;
    // Legacy names: these remain sums of full_weight/sigcm and its square.
    // They normalize reconstructed MC means; they are NOT predicted yields
    // and must not be compared directly with sumw_data or sum(phi[].sim).
    double sumw_sim  = 0.0, sumw2_sim  = 0.0;
    double mean_q2_data = 0.0, mean_xb_data = 0.0, mean_tprime_data = 0.0;
    double mean_q2_sim  = 0.0, mean_xb_sim  = 0.0, mean_tprime_sim  = 0.0;
    double mean_t_sim = 0.0; // Physical t, unlike t'=t-tmin, is PARTONS input.
    double mean_q2_vertex_sim = 0.0, mean_xb_vertex_sim = 0.0;
    // Generated-bin reference means, not reconstructed-bin mixture means.
    double mean_tprime_vertex_sim = 0.0, truth_response_sum = 0.0;
    double mean_q2_abs  = 0.0, mean_xb_abs  = 0.0, mean_tprime_abs  = 0.0;
    bool has_model_xsec = false;
    double epsilon = 0.0;
    double gamma_flux = 0.0;
    FourierFit fit_ratio;
    FourierFit fit_xsec;
    FourierFit fit_asym; // Legacy output record; LT' is not fitted in this pipeline.
    bool partons_ok = false;
    double partons_epsilon = 0.0, partons_electron_flux_xbq2 = 0.0;
    double partons_sigmaU = 0.0, partons_sigmaLT = 0.0, partons_sigmaTT = 0.0;
};

struct CutFlow {
    long long n_data_total = 0, n_data_pass = 0, n_data_inrange = 0;
    long long n_sim_total = 0, n_sim_pass = 0, n_sim_inrange = 0;
    long long n_sim_zero_sigcm = 0;
    long long n_sim_vertex_fallback = 0;
};

struct MissingMassSliceDiagnostic {
    std::unique_ptr<TH1D> data, exclusive;
    double exclusive_normalization_scale = 0.0;
};
