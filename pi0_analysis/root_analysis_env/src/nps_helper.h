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
    return (clusE > 0.55 && clusX > -30.0 && clusX < 30.0 && clusY > -36.0 && clusY < 36.0 && clusT > (150 - time_diff) && clusT < (150 + time_diff));
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



} // namespace nps

#endif // NPS_HELPER_H
