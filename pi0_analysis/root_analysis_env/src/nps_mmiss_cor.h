// ============================================================================
// File: nps_mmiss_fixed.h
// Purpose: Missing-mass correction factors and corrected Mx using Mazouz & Avnish
// ============================================================================

#ifndef NPS_MMISS_FIXED_H
#define NPS_MMISS_FIXED_H

#include "nps_helper.h"   // for default constants (optional)
#include <array>
#include <cmath>
#include <algorithm>
#include <cstdio>

namespace nps {

// Build π⁰ 4-vector from two photon detector positions (follows Python logic).
// Inputs:
//   e1,e2 : photon energies (GeV)
//   x1,y1,x2,y2 : cluster positions on NPS plane (cm)
//   z_nps : distance to NPS plane (cm) (positive)
//   theta_nps_deg : rotation angle in degrees (NPS -> hall), same sign convention as Python (-17.51)
// Output:
//   p4_pi_out : {E, px, py, pz}
// Returns true always (no external failure mode)
inline bool build_pi0_4vec_from_detector_python_logic(
    double e1, double x1, double y1,
    double e2, double x2, double y2,
    double z_nps,
    double theta_nps_deg,
    std::array<double,4> &p4_pi_out)
{
    // Norms (same as Python)
    double norm1 = std::sqrt(x1*x1 + y1*y1 + z_nps*z_nps);
    double norm2 = std::sqrt(x2*x2 + y2*y2 + z_nps*z_nps);

    // Avoid division by zero - if norm == 0 we'll set direction to (0,0,1)
    double u1x = 0.0, u1y = 0.0, u1z = 1.0;
    double u2x = 0.0, u2y = 0.0, u2z = 1.0;
    if (norm1 > 0.0) { u1x = x1 / norm1; u1y = y1 / norm1; u1z = z_nps / norm1; }
    if (norm2 > 0.0) { u2x = x2 / norm2; u2y = y2 / norm2; u2z = z_nps / norm2; }

    // Photon 3-momenta in NPS frame (E * unit vector)
    double ph1x = e1 * u1x;
    double ph1y = e1 * u1y;
    double ph1z = e1 * u1z;
    double ph2x = e2 * u2x;
    double ph2y = e2 * u2y;
    double ph2z = e2 * u2z;

    // Rotation NPS -> hall about Y by angle theta (degrees -> radians)
    const double theta = theta_nps_deg * M_PI / 180.0;
    const double cos_t = std::cos(theta);
    const double sin_t = std::sin(theta);

    // Apply same rotation as Python:
    // x' = cosθ * x - sinθ * z
    // y' = y
    // z' = sinθ * x + cosθ * z
    double p1x = cos_t * ph1x - sin_t * ph1z;
    double p1y = ph1y;
    double p1z = sin_t * ph1x + cos_t * ph1z;

    double p2x = cos_t * ph2x - sin_t * ph2z;
    double p2y = ph2y;
    double p2z = sin_t * ph2x + cos_t * ph2z;

    // π⁰ 4-vector (E, px, py, pz)
    p4_pi_out[0] = e1 + e2;
    p4_pi_out[1] = p1x + p2x;
    p4_pi_out[2] = p1y + p2y;
    p4_pi_out[3] = p1z + p2z;

    return true;
}

// ---------------------------------------------------------------------------
// Mazouz-style correction factor (mirror of Python invariant_missing_mass_correction_fac_mazouz).
// Inputs correspond to the Python function ordering (miss_mass_sq_pr kept for compatibility but not used).
// Returns corr_fac (double). If inputs invalid (m_inv==0) returns 0.
// ---------------------------------------------------------------------------
inline double invariant_missing_mass_correction_fac_mazouz_from_detector(
    double /*miss_mass_sq_pr*/,
    double e_beam,
    double Ee,
    double px_e,
    double py_e,
    double pz_e,
    double e1,
    double e2,
    double x1,
    double y1,
    double x2,
    double y2,
    double z_nps = kDefaultZ_NPS_cm,
    double theta_nps_deg = -17.51,
    double M_proton = kProtonMass_GeV,
    double M_pion0  = kPi0Mass_GeV,
    bool verbose = false)
{
    // Build pi0 4-vector exactly as in Python
    std::array<double,4> p4_pi0 = {0.0,0.0,0.0,0.0};
    build_pi0_4vec_from_detector_python_logic(e1,x1,y1, e2,x2,y2, z_nps, theta_nps_deg, p4_pi0);

    // q vector: k - k'  with k = (E_beam, 0,0, E_beam) and k' = (Ee, px_e, py_e, pz_e)
    double qx = -px_e;
    double qy = -py_e;
    double qz = e_beam - pz_e;
    double q_mag = std::sqrt(qx*qx + qy*qy + qz*qz);

    // pi0 momentum magnitude
    double pi_px = p4_pi0[1], pi_py = p4_pi0[2], pi_pz = p4_pi0[3];
    double pi_mag = std::sqrt(pi_px*pi_px + pi_py*pi_py + pi_pz*pi_pz);

    // angle theta between q and pi0: use clamped dot product
    double dot = 1.0;
    if (q_mag > 0.0 && pi_mag > 0.0) {
        double u1x = qx / q_mag, u1y = qy / q_mag, u1z = qz / q_mag;
        double u2x = pi_px / pi_mag, u2y = pi_py / pi_mag, u2z = pi_pz / pi_mag;
        dot = u1x*u2x + u1y*u2y + u1z*u2z;
        dot = std::max(-1.0, std::min(1.0, dot));
    }
    double theta = std::acos(dot);

    // Invariant mass of pi0 from its 4-vector
    double m_inv_sq = p4_pi0[0]*p4_pi0[0] - (pi_px*pi_px + pi_py*pi_py + pi_pz*pi_pz);
    double m_inv = (m_inv_sq > 0.0) ? std::sqrt(m_inv_sq) : 0.0;

    // energy transfer nu
    double nu = e_beam - Ee;

    if (m_inv == 0.0) {
        if (verbose) std::printf("[Mazouz-from-det] m_inv==0 -> returning 0\n");
        return 0.0;
    }

    // denom e1+e2 and pair factor (same as Python)
    double denom_e_sum = (e1 + e2);
    double e_pair_factor = 0.0;
    if (denom_e_sum != 0.0) e_pair_factor = (e1 * e2) / denom_e_sum;
    else if (verbose) std::printf("[Mazouz-from-det] e1+e2==0 -> e_pair_factor=0\n");

    // Eq. (7) in your Python: includes 2*sqrt(2) factor
    double corr_fac = (2.0 / m_inv) * (
        m_inv*m_inv - 2.0 * std::sqrt(2.0) * (nu + M_proton - q_mag * std::cos(theta)) * e_pair_factor
    );

    if (verbose) {
        std::printf("[Mazouz-from-det] q_mag=%.6g pi_mag=%.6g theta=%.6g nu=%.6g m_inv=%.6g\n",
                    q_mag, pi_mag, theta, nu, m_inv);
        std::printf("[Mazouz-from-det] e1=%.6g e2=%.6g e_pair_factor=%.6g corr_fac=%.12g\n",
                    e1, e2, e_pair_factor, corr_fac);
    }

    return corr_fac;
}

// ---------------------------------------------------------------------------
// Avnish-style correction factor (mirror of Python invariant_missing_mass_correction_fac_avnish).
// ---------------------------------------------------------------------------
inline double invariant_missing_mass_correction_fac_avnish_from_detector(
    double /*miss_mass_sq_pr*/,
    double e_beam,
    double Ee,
    double px_e,
    double py_e,
    double pz_e,
    double e1,
    double e2,
    double x1,
    double y1,
    double x2,
    double y2,
    double z_nps = kDefaultZ_NPS_cm,
    double theta_nps_deg = -17.51,
    double M_proton = kProtonMass_GeV,
    double M_pion0  = kPi0Mass_GeV,
    bool verbose = false)
{
    // Build pi0 4-vector exactly as in Python
    std::array<double,4> p4_pi0 = {0.0,0.0,0.0,0.0};
    build_pi0_4vec_from_detector_python_logic(e1,x1,y1, e2,x2,y2, z_nps, theta_nps_deg, p4_pi0);

    // q vector
    double qx = -px_e;
    double qy = -py_e;
    double qz = e_beam - pz_e;
    double q_mag = std::sqrt(qx*qx + qy*qy + qz*qz);

    // pi0 magnitude
    double pi_px = p4_pi0[1], pi_py = p4_pi0[2], pi_pz = p4_pi0[3];
    double pi_mag = std::sqrt(pi_px*pi_px + pi_py*pi_py + pi_pz*pi_pz);

    // angle theta between q and pi0
    double dot = 1.0;
    if (q_mag > 0.0 && pi_mag > 0.0) {
        double u1x = qx / q_mag, u1y = qy / q_mag, u1z = qz / q_mag;
        double u2x = pi_px / pi_mag, u2y = pi_py / pi_mag, u2z = pi_pz / pi_mag;
        dot = u1x*u2x + u1y*u2y + u1z*u2z;
        dot = std::max(-1.0, std::min(1.0, dot));
    }
    double theta = std::acos(dot);

    // invariant mass of pi0
    double m_inv_sq = p4_pi0[0]*p4_pi0[0] - (pi_px*pi_px + pi_py*pi_py + pi_pz*pi_pz);
    double m_inv = (m_inv_sq > 0.0) ? std::sqrt(m_inv_sq) : 0.0;

    double nu = e_beam - Ee;

    if (m_inv == 0.0) {
        if (verbose) std::printf("[Avnish-from-det] m_inv==0 -> returning 0\n");
        return 0.0;
    }

    double denom = (e1 + e2);
    double last_term = 0.0;
    if (denom != 0.0) {
        last_term = q_mag * std::cos(theta) * pi_mag / denom;
    } else if (verbose) {
        std::printf("[Avnish-from-det] e1+e2==0 -> last_term=0\n");
    }

    double corr_fac = (2.0 / m_inv) * ( m_inv*m_inv - (e1 + e2) * (nu + M_proton - last_term) );

    if (verbose) {
        std::printf("[Avnish-from-det] q_mag=%.6g pi_mag=%.6g theta=%.6g nu=%.6g m_inv=%.6g\n",
                    q_mag, pi_mag, theta, nu, m_inv);
        std::printf("[Avnish-from-det] e1=%.6g e2=%.6g denom=%.6g last_term=%.6g corr_fac=%.12g\n",
                    e1, e2, denom, last_term, corr_fac);
    }

    return corr_fac;
}

// ---------------------------------------------------------------------------
// Apply Avnish correction using detector inputs:
//   miss_mass_sq_temp = miss_mass_sq_pr - corr_fac*(m_inv - M_pion0)
//   return sqrt(max(0, miss_mass_sq_temp))
// This mirrors the Python apply path (calls Avnish corr function).
// ---------------------------------------------------------------------------
inline double invariant_missing_mass_corrected_avnish_from_detector(
    double miss_mass_sq_pr,
    double e_beam,
    double Ee,
    double px_e,
    double py_e,
    double pz_e,
    double e1,
    double e2,
    double x1,
    double y1,
    double x2,
    double y2,
    double z_nps = kDefaultZ_NPS_cm,
    double theta_nps_deg = -17.51,
    double M_proton = kProtonMass_GeV,
    double M_pion0  = kPi0Mass_GeV,
    bool verbose = false)
{
    // Build pi0 4-vector (we only need m_inv for application)
    std::array<double,4> p4_pi0 = {0.0,0.0,0.0,0.0};
    build_pi0_4vec_from_detector_python_logic(e1,x1,y1, e2,x2,y2, z_nps, theta_nps_deg, p4_pi0);

    // Compute Avnish correction factor (call the Avnish implementation)
    double corr_fac = invariant_missing_mass_correction_fac_avnish_from_detector(
        miss_mass_sq_pr, e_beam, Ee, px_e, py_e, pz_e,
        e1, e2, x1, y1, x2, y2,
        z_nps, theta_nps_deg, M_proton, M_pion0, verbose);

    // m_inv
    double pi_px = p4_pi0[1], pi_py = p4_pi0[2], pi_pz = p4_pi0[3];
    double m_inv_sq = p4_pi0[0]*p4_pi0[0] - (pi_px*pi_px + pi_py*pi_py + pi_pz*pi_pz);
    double m_inv = (m_inv_sq > 0.0) ? std::sqrt(m_inv_sq) : 0.0;

    double miss_mass_sq_temp = miss_mass_sq_pr - corr_fac * (m_inv - M_pion0);

    if (verbose) {
        std::printf("[Avnish-apply-from-det] corr_fac=%.12g m_inv=%.6g M_pion0=%.6g msq_pr=%.12g msq_temp=%.12g\n",
                    corr_fac, m_inv, M_pion0, miss_mass_sq_pr, miss_mass_sq_temp);
    }

    return (miss_mass_sq_temp > 0.0) ? std::sqrt(miss_mass_sq_temp) : 0.0;
}

} // namespace nps

#endif // NPS_MMISS_FIXED_H
