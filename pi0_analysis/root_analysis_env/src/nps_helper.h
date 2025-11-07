// ============================================================================
// File: nps_helper.h
// Author: Avnish Singh (helpers refined)
// Purpose: Helper and utility functions for NPS π0 and DVCS analysis
// ============================================================================

#ifndef NPS_HELPER_H
#define NPS_HELPER_H

#include <array>
#include <cmath>
#include <algorithm>

using Vec3 = std::array<double,3>;
using Vec4 = std::array<double,4>;

namespace nps {

// -----------------------------------------------------------------
// Physical constants (GeV units, PDG 2024 values)
// -----------------------------------------------------------------
constexpr double kElectronMass_GeV = 0.0005109989461; // GeV
constexpr double kProtonMass_GeV   = 0.9382720813;    // GeV
constexpr double kPi0Mass_GeV      = 0.1349768;       // GeV
constexpr double kSpeedOfLight_cm_ns = 29.9792458;    // cm/ns
constexpr double kDefaultZ_NPS_cm = 407.0;
constexpr double kDeg2Rad = M_PI / 180.0;
constexpr double kRad2Deg = 180.0 / M_PI;

// small utility
inline double sqr(double x) { return x*x; }
inline double safe_sqrt(double x) { return (x <= 0.0) ? 0.0 : std::sqrt(x); }

// --------------------
// 3-vector & 4-vector helpers
// --------------------
inline double dot3(const Vec3 &a, const Vec3 &b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline double mass2_4vec(const Vec4 &p4) {
    // metric (+,-,-,-)
    return p4[0]*p4[0] - (p4[1]*p4[1] + p4[2]*p4[2] + p4[3]*p4[3]);
}
inline Vec4 add4(const Vec4 &a, const Vec4 &b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3]};
}
inline Vec4 sub4(const Vec4 &a, const Vec4 &b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3]};
}

// Rotate a 3-vector about Y axis by theta_deg (right-hand rule).
inline Vec3 rotate_y(const Vec3 &in, double theta_deg) {
    const double theta = theta_deg * kDeg2Rad;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    return Vec3{ c*in[0] - s*in[2], in[1], s*in[0] + c*in[2] };
}

// --------------------
// Cuts
// --------------------
inline bool hms_electron_cuts(double edtm_tdc, double h_delta, double h_gtr_th,
                              double h_gtr_ph, double h_cer_npeSum, double h_cal_etotnorm,
                              double h_react_z) noexcept
{
    // fast, ordered checks (fail-fast)
    if (edtm_tdc > 0.1) return false;
    if (h_react_z < -4.0 || h_react_z > 4.0) return false;
    if (h_delta < -15.0 || h_delta > 15.0) return false;
    if (h_gtr_th < -0.1 || h_gtr_th > 0.1) return false;
    if (h_gtr_ph < -0.04 || h_gtr_ph > 0.04) return false;
    if (h_cer_npeSum < 1.5) return false;
    if (h_cal_etotnorm < 0.7 || h_cal_etotnorm > 1.2) return false;
    return true;
}

inline bool nps_spatial_energy_cuts(double clusE, double clusX, double clusY, double clusT, double time_diff) noexcept {
    return (clusE > 0.8 && clusX > -30.0 && clusX < 30.0 && clusY > -36.0 && clusY < 36.0 && clusT > (150 - time_diff) && clusT < (150 + time_diff));
}

// --------------------
// Merge clusters in-place (arrays E,X,Y,T length n).
// - space2_thresh_cm2: squared spatial threshold in cm^2 (compare dx^2+dy^2 to it)
// - time_thresh_ns: time window in ns
// Implementation is O(n^2) which is fine for small n (typical nclus << 100).
// If you expect large n, use spatial hashing / kd-tree.
// --------------------
inline void mergeClusters(double *E, double *X, double *Y, double *T, int n,
                          double space2_thresh_cm2 = 50.0, double time_thresh_ns = 2.0)
{
    if (n <= 1) return;
    bool merged_any = true;
    // repeat passes to allow cascade merges (small n -> acceptable)
    while (merged_any) {
        merged_any = false;
        for (int i = 0; i < n-1; ++i) {
            if (!(E[i] > 0.0)) continue;
            for (int j = i+1; j < n; ++j) {
                if (!(E[j] > 0.0)) continue;
                const double dx = X[i] - X[j];
                const double dy = Y[i] - Y[j];
                const double space2 = dx*dx + dy*dy;
                const double dt = std::fabs(T[i] - T[j]);
                if (space2 <= space2_thresh_cm2 && dt <= time_thresh_ns) {
                    const double Ei = E[i], Ej = E[j];
                    const double Etot = Ei + Ej;
                    if (Etot <= 0.0) {
                        E[i] = 0.0; E[j] = 0.0;
                    } else {
                        X[i] = (X[i]*Ei + X[j]*Ej) / Etot;
                        Y[i] = (Y[i]*Ei + Y[j]*Ej) / Etot;
                        T[i] = (T[i]*Ei + T[j]*Ej) / Etot;
                        E[i] = Etot;
                        E[j] = 0.0;
                    }
                    merged_any = true;
                }
            }
        }
    }
}

// Pack clusters: move E>0 entries to front, return new count
inline int packClusters(double *E, double *X, double *Y, double *T, int n) {
    int w = 0;
    for (int r = 0; r < n; ++r) {
        if (E[r] > 0.0) {
            if (w != r) {
                E[w] = E[r]; X[w] = X[r]; Y[w] = Y[r]; T[w] = T[r];
            }
            ++w;
        }
    }
    return w;
}

// --------------------
// Photon & invariant-mass utilities
// --------------------
inline Vec4 photon4vector(double E, double x, double y, double z_nps = kDefaultZ_NPS_cm) {
    Vec3 r = {x, y, z_nps};
    const double rnorm = std::sqrt(dot3(r,r));
    if (rnorm <= 0.0) return Vec4{E, 0.0, 0.0, E}; // fallback
    const double invnorm = 1.0 / rnorm;
    Vec3 u = { r[0]*invnorm, r[1]*invnorm, r[2]*invnorm };
    return Vec4{ E, E*u[0], E*u[1], E*u[2] };
}

inline double invariant_mass_pi0(double e1, double e2,
                                 double x1, double x2,
                                 double y1, double y2,
                                 double z_nps = kDefaultZ_NPS_cm)
{
    const Vec4 p1 = photon4vector(e1, x1, y1, z_nps);
    const Vec4 p2 = photon4vector(e2, x2, y2, z_nps);
    const Vec4 tot = add4(p1, p2);
    const double m2 = mass2_4vec(tot);
    return safe_sqrt(m2);
}

// --------------------
// Missing mass for ep -> e' + pi0 + X (proton missing mass) using:
// p_in = e_beam + proton_at_rest
// p_out = e_out + pi0 (where pion is two photons)
// Note: theta_e_deg is the electron scattering angle in degrees.
// e_out energy Ee must be provided (GeV) and electron momentum (if needed) is computed.
// Theta conventions must match the tree (check units: HgtrTh often in radians).
// --------------------
inline double missing_mass_proton_pi0(double e_beam,
                                      double Ee, double px_e,
                                      double py_e, double pz_e,
                                      double e1, double e2,
                                      double x1, double y1,
                                      double x2, double y2,
                                      double z_nps = kDefaultZ_NPS_cm,
                                      double theta_nps_deg = -17.51)
{
    // incident electron (assume beam along +z). Use relativistic energy ~ momentum for beam since m_e negligible.
    const double pbeam = std::sqrt(std::max(0.0, e_beam*e_beam - kElectronMass_GeV*kElectronMass_GeV));
    Vec4 p4_ein = { e_beam, 0.0, 0.0, pbeam };

    // target proton at rest
    Vec4 p4_pin = { kProtonMass_GeV, 0.0, 0.0, 0.0 };

    // outgoing electron: Ee and theta
    const double p_eout = std::sqrt(std::max(0.0, Ee*Ee - kElectronMass_GeV*kElectronMass_GeV));
    Vec4 p4_eout = { Ee, px_e, py_e, pz_e };

    // photons -> pion
    const Vec4 ph1 = photon4vector(e1, x1, y1, z_nps);
    const Vec4 ph2 = photon4vector(e2, x2, y2, z_nps);

    // rotate photon momenta from NPS -> hall coords
    Vec3 p1_nps = { ph1[1], ph1[2], ph1[3] };
    Vec3 p2_nps = { ph2[1], ph2[2], ph2[3] };
    Vec3 p1_hall = rotate_y(p1_nps, theta_nps_deg);
    Vec3 p2_hall = rotate_y(p2_nps, theta_nps_deg);

    Vec4 p4_ph1 = { ph1[0], p1_hall[0], p1_hall[1], p1_hall[2] };
    Vec4 p4_ph2 = { ph2[0], p2_hall[0], p2_hall[1], p2_hall[2] };

    Vec4 p4_pi = add4(p4_ph1, p4_ph2);

    Vec4 p4_in = add4(p4_ein, p4_pin);
    Vec4 p4_out = add4(p4_eout, p4_pi);

    Vec4 p4_miss = sub4(p4_in, p4_out);
    const double m2 = mass2_4vec(p4_miss);
    return safe_sqrt(m2);
}

// --------------------
// Missing mass for DVCS-like events (ep -> e' gamma p')
// Single photon case. Returns missing mass (GeV) of the undetected proton.
// --------------------
inline double missing_mass_dvcs(double e_beam,
                                double Ee, double px_e,
                                double py_e, double pz_e,
                                double e_photon, double x_ph, double y_ph,
                                double z_nps = kDefaultZ_NPS_cm,
                                double theta_nps_deg = -17.51)
{
    const double pbeam = std::sqrt(std::max(0.0, e_beam*e_beam - kElectronMass_GeV*kElectronMass_GeV));
    Vec4 p4_ein = { e_beam, 0.0, 0.0, pbeam };
    Vec4 p4_pin = { kProtonMass_GeV, 0.0, 0.0, 0.0 };

    const double p_eout = std::sqrt(std::max(0.0, Ee*Ee - kElectronMass_GeV*kElectronMass_GeV));
    Vec4 p4_eout = { Ee, px_e, py_e, pz_e };

    const Vec4 ph = photon4vector(e_photon, x_ph, y_ph, z_nps);
    Vec3 p_nps = { ph[1], ph[2], ph[3] };
    Vec3 p_hall = rotate_y(p_nps, theta_nps_deg);
    Vec4 p4_gamma = { ph[0], p_hall[0], p_hall[1], p_hall[2] };

    Vec4 p4_miss = sub4(add4(p4_ein, p4_pin), add4(p4_eout, p4_gamma));
    const double m2 = mass2_4vec(p4_miss);
    return safe_sqrt(m2);
}

// --------------------
// fit_pi0_peak: fit invariant mass hist in window [fit_lo, fit_hi] (GeV)
// - h: pointer to TH1D (mass histogram in GeV)
// - fit_lo/fit_hi: fit window (default 0.120 - 0.140 GeV)
// - draw: whether to draw an annotated canvas and save it to PNG/root file
// - outPlotDir: used for saving PNG and writing canvas into output file if not empty
// Returns FitResult struct with numbers and also prints a verbose report.
// --------------------

// Result struct returned by fit_pi0_peak
struct FitResult {
    double mean = 0.0;               // narrow gaussian mean [GeV]
    double mean_err = 0.0;
    double sigma = 0.0;              // narrow gaussian sigma [GeV]
    double sigma_err = 0.0;
    double amp = 0.0;                // narrow gaussian amplitude (height)
    double amp_err = 0.0;
    double chi2 = 0.0;
    int ndf = 0;

    // Raw histogram integral in ROI (0.12-0.14 GeV)
    double integral = 0.0;
    double integral_err = 0.0;

    // Fitted contributions (integrated over ROI)
    double gauss_integral = 0.0;        // narrow gaussian (signal)
    double gauss_integral_err = 0.0;
    double broad_gauss_integral = 0.0;  // broad gaussian (background)
    double broad_gauss_integral_err = 0.0;
    double pedestal_integral = 0.0;
    double pedestal_integral_err = 0.0;

    // Combined background (broad + pedestal)
    double bg_integral = 0.0;
    double bg_integral_err = 0.0;

    // Background-subtracted yield and error (raw - bg)
    double yield_bgsub = 0.0;
    double yield_bgsub_err = 0.0;

    double significance = 0.0;        // S / sqrt(S+B) approx

    std::string canvas_name;          // name of the created canvas (if draw==true)
};


// Fit pi0 peak using: narrow gaussian (signal) + broad gaussian (low-mass background) + pedestal.
// IMPORTANT: fit is performed in [fit_lo, fit_hi] (user-specified), but ALL integrals & returned yields
// are calculated in the FIXED ROI [yield_lo, yield_hi] = [0.120, 0.140] GeV.
//
// Parameters:
//  - h         : pointer to TH1D (invariant mass histogram in GeV, with reasonable binning)
//  - fit_lo    : lower bound of FIT window (GeV)  -> this is where the function is fit
//  - fit_hi    : upper bound of FIT window (GeV)
//  - draw      : if true draw annotated canvas and save PNG to outPlotDir (if provided)
//  - outPlotDir: directory to save PNG (if not empty)
//  - run       : run number used for naming the canvas/PNG
inline FitResult fit_pi0_peak(TH1D *h,
                              double fit_lo = 0.02, double fit_hi = 0.30,
                              bool draw = true, const TString &outPlotDir = "", int run = -1)
{
    FitResult R;
    if (!h) {
        std::cerr << "[fit_pi0_peak] ERROR: null histogram pointer\n";
        return R;
    }
    if (fit_lo >= fit_hi) {
        std::cerr << "[fit_pi0_peak] ERROR: fit_lo >= fit_hi\n";
        return R;
    }

    // Fixed yield ROI (always use this for integrals & yield)
    const double yield_lo = 0.120; // GeV
    const double yield_hi = 0.140; // GeV

    // Make sure histogram has valid bin errors
    h->Sumw2();

    // --- initial guesses from histogram (local peak approximations) ---
    const int bin_max = h->GetMaximumBin();
    const double x_max = h->GetBinCenter(bin_max);
    const double y_max = h->GetBinContent(bin_max);

    // sensible starting values
    const double A1_init = std::max(1.0, y_max);        // narrow peak height guess
    const double mean1_init = nps::kPi0Mass_GeV;       // ~0.13498
    const double sigma1_init = 0.004;                  // ~4 MeV

    const double A2_init = 0.3 * y_max;                // broad gaussian height
    const double mean2_init = std::max(fit_lo, (fit_lo + mean1_init) * 0.5);
    const double sigma2_init = 0.02;                   // ~20 MeV

    const double ped_init = 0.0;                       // pedestal initial

    // --- define fit function: gaus(0) + gaus(3) + pol0(6) over user fit range ---
    TF1 *f = new TF1("f_pi0", "gaus(0) + gaus(3) + pol0(6)", fit_lo, fit_hi);
    f->SetParameters(A1_init, mean1_init, sigma1_init,
                     A2_init, mean2_init, sigma2_init,
                     ped_init);

    // set safe parameter limits to keep fit stable
    f->SetParLimits(0, 0.0, y_max * 50.0);                     // A1 >= 0
    f->SetParLimits(1, fit_lo + 1e-4, fit_hi - 1e-4);         // mean1 inside fit window
    f->SetParLimits(2, 5e-4, 0.05);                           // sigma1 between 0.5 MeV and 50 MeV

    f->SetParLimits(3, 0.0, y_max * 50.0);                     // A2 >= 0
    f->SetParLimits(4, std::max(fit_lo - 0.02, 0.0), mean1_init - 0.002);
    f->SetParLimits(5, 0.005, 0.12);                           // broad sigma: 5 MeV - 120 MeV

    f->SetParLimits(6, -y_max, y_max);                         // pedestal can be +/- y_max

    // --- perform fit quietly (R = use fit range, Q = quiet, N = do not store) ---
    int fit_status = h->Fit(f, "RQN");

    // --- extract fitted parameters & their errors ---
    const double A1  = f->GetParameter(0);
    const double mean1 = f->GetParameter(1);
    const double sigma1 = f->GetParameter(2);

    const double A2  = f->GetParameter(3);
    const double mean2 = f->GetParameter(4);
    const double sigma2 = f->GetParameter(5);

    const double ped = f->GetParameter(6);

    const double errA1 = f->GetParError(0);
    const double errMean1 = f->GetParError(1);
    const double errSigma1 = f->GetParError(2);
    const double errA2 = f->GetParError(3);
    const double errMean2 = f->GetParError(4);
    const double errSigma2 = f->GetParError(5);
    const double errPed = f->GetParError(6);

    const double chi2 = f->GetChisquare();
    const int ndf = f->GetNDF();

    // --- compute raw histogram integral in fixed ROI [0.120,0.140] GeV (use IntegralAndError) ---
    // Use bin indices from FindBin (IntegralAndError expects bin numbers)
    const int bin_y1 = h->FindBin(yield_lo);
    const int bin_y2 = h->FindBin(yield_hi);

    double hist_int = 0.0;
    double hist_int_err = 0.0;
    // IntegralAndError returns the integral and sets hist_int_err by reference
    hist_int = h->IntegralAndError(bin_y1, bin_y2, hist_int_err);

    // --- Compute analytic integrals of fit components over the ROI ---
    const double sqrt2 = std::sqrt(2.0);
    const double sqrt2pi = std::sqrt(2.0 * M_PI);
    // Corrected binwise gaussian integral that returns model COUNTS (comparable to TH1::Integral)
    auto gauss_integral = [&](double Amp, double mu, double sig) -> double {
        if (sig <= 0.0) return 0.0;

        // bin indices covering ROI
        const int bin_start = h->FindBin(yield_lo);
        const int bin_end   = h->FindBin(yield_hi);

        static int __tmpGcount = 0;
        ++__tmpGcount;
        TString tmpname = TString::Format("tmp_gaus_binint_%d", __tmpGcount);

        // create TF1 using ROOT 'gaus' semantics (par0 = amplitude in same units as histogram bin content)
        TF1 *g = new TF1(tmpname, "gaus", yield_lo, yield_hi);
        g->SetParameters(Amp, mu, sig);

        double model_counts = 0.0;

        for (int b = bin_start; b <= bin_end; ++b) {
            const double xlow  = h->GetXaxis()->GetBinLowEdge(b);
            const double xhigh = h->GetXaxis()->GetBinUpEdge(b);
            const double lo = std::max(xlow, yield_lo);
            const double hi = std::min(xhigh, yield_hi);
            if (hi <= lo) continue;

            // TF1::Integral gives continuous integral ∫ f(x) dx (units: f * x)
            double cont_int = g->Integral(lo, hi);

            // **convert to counts**: since f(x) was fit to bin counts (counts per bin),
            // the relation is: counts_in_bin ≈ cont_int / bin_width.
            const double binw = (xhigh - xlow);
            if (binw <= 0) continue;
            double counts_in_bin = cont_int / binw;

            model_counts += counts_in_bin;
        }

        delete g;
        return model_counts;
    };



    // narrow gaussian (signal) integral in ROI
    const double gauss1_area = (sigma1 > 0.0) ? gauss_integral(A1, mean1, sigma1) : 0.0;
    // broad gaussian (background) integral in ROI
    const double gauss2_area = (sigma2 > 0.0) ? gauss_integral(A2, mean2, sigma2) : 0.0;
    // pedestal integral in ROI
    const double ped_area = ped * (yield_hi - yield_lo);

    // --- approximate errors on integrals (uncorrelated propagation) ---
    auto sqr = [](double x){ return x*x; };

    // For gaussian area A * sigma * sqrt(2*pi) * factor(erf difference) the partials:
    // partial w.r.t A = sigma * sqrt(2*pi) * factor
    // partial w.r.t sigma ≈ A * sqrt(2*pi) * factor  (ignores derivative of erf argument; good enough as approximation)
    double erf1_hi = (sigma1>0) ? std::erf((yield_hi - mean1) / (sqrt2 * sigma1)) : 0.0;
    double erf1_lo = (sigma1>0) ? std::erf((yield_lo - mean1) / (sqrt2 * sigma1)) : 0.0;
    double factor1 = 0.5 * (erf1_hi - erf1_lo);
    double partialA1 = sigma1 * sqrt2pi * factor1;
    double partialSigma1 = A1 * sqrt2pi * factor1;
    const double gauss1_err = std::sqrt( sqr(errA1 * partialA1) + sqr(errSigma1 * partialSigma1) );

    double erf2_hi = (sigma2>0) ? std::erf((yield_hi - mean2) / (sqrt2 * sigma2)) : 0.0;
    double erf2_lo = (sigma2>0) ? std::erf((yield_lo - mean2) / (sqrt2 * sigma2)) : 0.0;
    double factor2 = 0.5 * (erf2_hi - erf2_lo);
    double partialA2 = sigma2 * sqrt2pi * factor2;
    double partialSigma2 = A2 * sqrt2pi * factor2;
    const double gauss2_err = std::sqrt( sqr(errA2 * partialA2) + sqr(errSigma2 * partialSigma2) );

    const double ped_err = std::fabs(errPed) * (yield_hi - yield_lo);

    // background = broad gaussian + pedestal
    const double bg = gauss2_area + ped_area;
    const double bg_err = std::sqrt( sqr(gauss2_err) + sqr(ped_err) );

    // background-subtracted yield and propagated error
    const double yield_bgsub = hist_int - bg;
    const double yield_bgsub_err = std::sqrt( sqr(hist_int_err) + sqr(bg_err) );

    // approximate signal significance (narrow gaussian is S, broad+ped is B)
    double significance = 0.0;
    if ((gauss1_area + bg) > 0.0) significance = gauss1_area / std::sqrt(gauss1_area + bg);

    // --- fill result structure ---
    R.mean = mean1; R.mean_err = errMean1;
    R.sigma = sigma1; R.sigma_err = errSigma1;
    R.amp = A1; R.amp_err = errA1;
    R.chi2 = chi2; R.ndf = ndf;

    R.integral = hist_int; R.integral_err = hist_int_err;
    R.gauss_integral = gauss1_area; R.gauss_integral_err = gauss1_err;
    R.broad_gauss_integral = gauss2_area; R.broad_gauss_integral_err = gauss2_err;
    R.pedestal_integral = ped_area; R.pedestal_integral_err = ped_err;

    R.bg_integral = bg; R.bg_integral_err = bg_err;
    R.yield_bgsub = yield_bgsub; R.yield_bgsub_err = yield_bgsub_err;
    R.significance = significance;

    // --- Print a verbose, informative report (units: MeV printed for mean/sigma) ---
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\n========================================\n";
    std::cout << Form("π0 fit summary (fit range = [%.3f, %.3f] GeV, ROI = [%.3f, %.3f] GeV)\n", fit_lo, fit_hi, yield_lo, yield_hi);
    std::cout << "----------------------------------------\n";
    std::cout << Form("Fit status: %d    chi2 = %.3f   ndf = %d   chi2/ndf = %.3f\n",
                      fit_status, chi2, ndf, (ndf>0 ? chi2/ndf : -1.0));
    std::cout << Form("Narrow Gaussian (signal): mean = %.2f ± %.2f MeV, sigma = %.2f ± %.2f MeV, A = %.3f ± %.3f\n",
                      mean1*1000.0, errMean1*1000.0, sigma1*1000.0, errSigma1*1000.0, A1, errA1);
    std::cout << Form("Broad Gaussian (bg)    : mean = %.2f ± %.2f MeV, sigma = %.2f ± %.2f MeV, A = %.3f ± %.3f\n",
                      mean2*1000.0, errMean2*1000.0, sigma2*1000.0, errSigma2*1000.0, A2, errA2);
    std::cout << Form("Pedestal (constant)     : integral = %.3f ± %.3f (over ROI)\n", ped_area, ped_err);
    std::cout << "----------------------------------------\n";
    std::cout << Form("Raw histogram integral in ROI [%g,%g] GeV = %.3f ± %.3f\n", yield_lo, yield_hi, hist_int, hist_int_err);
    std::cout << Form("Fitted narrow-gauss integral (signal) = %.3f ± %.3f\n", gauss1_area, gauss1_err);
    std::cout << Form("Fitted broad-gauss integral (bg)      = %.3f ± %.3f\n", gauss2_area, gauss2_err);
    std::cout << Form("Fitted pedestal integral (bg)         = %.3f ± %.3f\n", ped_area, ped_err);
    std::cout << Form("Total fitted background (bg)          = %.3f ± %.3f\n", bg, bg_err);
    std::cout << Form("Background-subtracted yield (raw - bg)= %.3f ± %.3f\n", yield_bgsub, yield_bgsub_err);
    std::cout << Form("Approx. significance S/sqrt(S+B) = %.3f\n", significance);
    std::cout << "========================================\n\n";

    // --- Optionally draw annotated canvas showing histogram, fit, components, ROI lines ---
    if (draw) {
        TString cname = (run >= 0) ? TString::Format("c_pi0_fit_run%d", run) : TString("c_pi0_fit");
        R.canvas_name = std::string(cname.Data());
        TCanvas *c = new TCanvas(cname, "π0 fit", 900, 700);
        gPad->SetRightMargin(0.05);
        h->Draw("E"); // draw with errors

        // draw the total fit (clone to ensure independent style)
        TF1 *fdisplay = (TF1*)f->Clone(TString::Format("%s_display", cname.Data()));
        fdisplay->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        fdisplay->SetLineColor(kRed); fdisplay->SetLineWidth(2);
        fdisplay->Draw("same");

        // narrow gaussian component (blue)
        TF1 *g1 = new TF1(TString::Format("%s_g1", cname.Data()), "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        g1->SetParameters(A1, mean1, sigma1); g1->SetLineColor(kBlue); g1->SetLineWidth(2); g1->Draw("same");

        // broad gaussian component (green)
        TF1 *g2 = new TF1(TString::Format("%s_g2", cname.Data()), "gaus", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        g2->SetParameters(A2, mean2, sigma2); g2->SetLineColor(kGreen+2); g2->SetLineWidth(2); g2->Draw("same");

        // pedestal line
        TF1 *pf = new TF1(TString::Format("%s_ped", cname.Data()), Form("%g", ped), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        pf->SetLineColor(kOrange+1); pf->SetLineStyle(2); pf->SetLineWidth(2); pf->Draw("same");

        // ROI vertical lines (yield range)
        TLine *Llo = new TLine(yield_lo, h->GetMinimum(), yield_lo, h->GetMaximum());
        TLine *Lhi = new TLine(yield_hi, h->GetMinimum(), yield_hi, h->GetMaximum());
        Llo->SetLineStyle(2); Llo->SetLineColor(kViolet); Llo->SetLineWidth(2); Llo->Draw("same");
        Lhi->SetLineStyle(2); Lhi->SetLineColor(kViolet); Lhi->SetLineWidth(2); Lhi->Draw("same");

        // annotation text (top-left)
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03);
        tx.DrawLatex(0.12, 0.92, Form("Run: %d", run));
        tx.DrawLatex(0.12, 0.88, Form("Fit range: [%.3f, %.3f] GeV", fit_lo, fit_hi));
        tx.DrawLatex(0.12, 0.84, Form("ROI (yield): [%.3f, %.3f] GeV", yield_lo, yield_hi));
        tx.DrawLatex(0.12, 0.80, Form("Yield (raw-bg) = %.3f ± %.3f", yield_bgsub, yield_bgsub_err));
        tx.DrawLatex(0.12, 0.76, Form("Signal (gaus) = %.3f ± %.3f", gauss1_area, gauss1_err));
        tx.DrawLatex(0.12, 0.72, Form("BG (broad+ped) = %.3f ± %.3f", bg, bg_err));
        tx.DrawLatex(0.12, 0.68, Form("S/sqrt(S+B) = %.3f", significance));

        // legend with components
        TLegend *leg = new TLegend(0.60, 0.62, 0.92, 0.88);
        leg->SetBorderSize(0); leg->SetFillColor(0);
        leg->AddEntry(h, "Mass spectrum", "lep");
        leg->AddEntry(fdisplay, "Total fit", "l");
        leg->AddEntry(g1, "Narrow Gauss (signal)", "l");
        leg->AddEntry(g2, "Broad Gauss (bg)", "l");
        leg->AddEntry(pf, "Pedestal", "l");
        leg->Draw("same");

        // Save PNG if directory provided
        if (outPlotDir != "") {
            TString png = TString::Format("%s/pi0_fit_run%d.png", outPlotDir.Data(), run);
            c->SaveAs(png);
        }

        // Keep canvas in memory so caller can write it to output TFile if desired.
        // (Don't delete c/fdisplay/g1/g2/pf/Llo/Lhi/leg/tx here.)
    } // end draw

    // If we didn't draw, we can clean f to avoid leaks
    if (!draw) delete f;

    return R;
}



} // namespace nps

#endif // NPS_HELPER_H
