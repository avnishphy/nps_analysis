// ============================================================================
// File: nps_mmiss.h
// Purpose: Missing-mass correction factors and corrected Mx using Mazouz & Avnish
//          Uses detector inputs (beam, scattered-e, two photons) — builds p4_pi0
// Author: adapted for user's codebase
// ============================================================================

#ifndef NPS_MMISS_H
#define NPS_MMISS_H

// reuse helpers/constants from your project; do NOT re-declare shared symbols here
#include "nps_helper.h"

// ROOT drawing headers (make sure ROOT include paths available when compiling)
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TString.h>

#include <vector>
#include <cmath>
#include <cstdio>
#include <string>

namespace nps {

// ---------------------------------------------------------------------------
// Event struct (light-weight) used by the plotting helper
// ---------------------------------------------------------------------------
struct MMEvent {
    double miss_mass_sq_pr; // Mx^2 (GeV^2) measured/constructed (if unknown supply <=0 and we'll compute)
    double e_beam;          // beam energy (GeV)
    double Ee;              // scattered electron energy (GeV)
    double px_e, py_e, pz_e;// scattered electron momentum components (GeV)
    double e1, x1, y1;      // photon1 (energy GeV, hit x,y cm on NPS plane)
    double e2, x2, y2;      // photon2
};

// ---------------------------------------------------------------------------
// Helper: build pi0 4-vector (hall coordinates) from two photons measured in NPS
// Follows the same logic as missing_mass_proton_pi0 in your codebase.
// Returns true on success; p4_pi_out is (E,px,py,pz).
// ---------------------------------------------------------------------------
inline bool build_pi0_4vec_from_detector(double e1, double x1, double y1,
                                         double e2, double x2, double y2,
                                         double z_nps,
                                         double theta_nps_deg,
                                         std::array<double,4> &p4_pi_out)
{
    // photon four-vectors in NPS frame (photon4vector available in nps_helper.h)
    const std::array<double,4> ph1 = photon4vector(e1, x1, y1, z_nps);
    const std::array<double,4> ph2 = photon4vector(e2, x2, y2, z_nps);

    // convert 3-momenta to Vec3 and rotate to hall
    Vec3 p1_nps = { ph1[1], ph1[2], ph1[3] };
    Vec3 p2_nps = { ph2[1], ph2[2], ph2[3] };
    Vec3 p1_hall = rotate_y(p1_nps, theta_nps_deg);
    Vec3 p2_hall = rotate_y(p2_nps, theta_nps_deg);

    p4_pi_out[0] = ph1[0] + ph2[0];
    p4_pi_out[1] = p1_hall[0] + p2_hall[0];
    p4_pi_out[2] = p1_hall[1] + p2_hall[1];
    p4_pi_out[3] = p1_hall[2] + p2_hall[2];
    return true;
}

// ---------------------------------------------------------------------------
// Mazouz-style correction factor (detector-input wrapper).
// Uses the same detector inputs as missing_mass_proton_pi0 to build p4_pi0.
// Returns corr_fac (does NOT change Mx^2).
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
    // build pi0 4-vector
    std::array<double,4> p4_pi0 = {0,0,0,0};
    build_pi0_4vec_from_detector(e1,x1,y1, e2,x2,y2, z_nps, theta_nps_deg, p4_pi0);

    // q vector convention used in your Python snippet: q = ( -px_e, -py_e, E_beam - pz_e )
    double qx = -px_e, qy = -py_e, qz = e_beam - pz_e;
    double q_mag = std::sqrt(qx*qx + qy*qy + qz*qz);

    double pi_px = p4_pi0[1], pi_py = p4_pi0[2], pi_pz = p4_pi0[3];
    double pi_mag = std::sqrt(pi_px*pi_px + pi_py*pi_py + pi_pz*pi_pz);

    double dot = 1.0;
    if (q_mag > 0.0 && pi_mag > 0.0) {
        double u1x = qx/q_mag, u1y = qy/q_mag, u1z = qz/q_mag;
        double u2x = pi_px/pi_mag, u2y = pi_py/pi_mag, u2z = pi_pz/pi_mag;
        dot = u1x*u2x + u1y*u2y + u1z*u2z;
        dot = std::max(-1.0, std::min(1.0, dot));
    }
    double theta = std::acos(dot);

    double m_inv_sq = p4_pi0[0]*p4_pi0[0] - (pi_px*pi_px + pi_py*pi_py + pi_pz*pi_pz);
    double m_inv = (m_inv_sq > 0.0) ? std::sqrt(m_inv_sq) : 0.0;
    double nu = e_beam - Ee;

    if (m_inv == 0.0) {
        if (verbose) std::printf("[Mazouz-from-det] m_inv==0 -> returning 0\n");
        return 0.0;
    }

    double denom_e_sum = (e1 + e2);
    double e_pair_factor = 0.0;
    if (denom_e_sum != 0.0) e_pair_factor = (e1 * e2) / denom_e_sum;
    else if (verbose) std::printf("[Mazouz-from-det] e1+e2==0 -> e_pair_factor=0\n");

    // Eq. (7) from your snippet:
    double corr_fac = (2.0 / m_inv) *
        ( m_inv*m_inv - 2.0*std::sqrt(2.0) * (nu + M_proton - q_mag * std::cos(theta)) * e_pair_factor );

    if (verbose) {
        std::printf("[Mazouz-from-det] q_mag=%.6g pi_mag=%.6g theta=%.6g nu=%.6g m_inv=%.6g\n",
                    q_mag, pi_mag, theta, nu, m_inv);
        std::printf("[Mazouz-from-det] e1=%.6g e2=%.6g e_pair_factor=%.6g corr_fac=%.12g\n",
                    e1, e2, e_pair_factor, corr_fac);
    }

    return corr_fac;
}

// ---------------------------------------------------------------------------
// Avnish-style correction factor (detector-input wrapper).
// Returns corr_fac. Use invariant_missing_mass_corrected_avnish_from_detector(...) to apply.
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
    // build pi0
    std::array<double,4> p4_pi0 = {0,0,0,0};
    build_pi0_4vec_from_detector(e1,x1,y1, e2,x2,y2, z_nps, theta_nps_deg, p4_pi0);

    double qx = -px_e, qy = -py_e, qz = e_beam - pz_e;
    double q_mag = std::sqrt(qx*qx + qy*qy + qz*qz);

    double pi_px = p4_pi0[1], pi_py = p4_pi0[2], pi_pz = p4_pi0[3];
    double pi_mag = std::sqrt(pi_px*pi_px + pi_py*pi_py + pi_pz*pi_pz);

    double dot = 1.0;
    if (q_mag > 0.0 && pi_mag > 0.0) {
        double u1x = qx/q_mag, u1y = qy/q_mag, u1z = qz/q_mag;
        double u2x = pi_px/pi_mag, u2y = pi_py/pi_mag, u2z = pi_pz/pi_mag;
        dot = u1x*u2x + u1y*u2y + u1z*u2z;
        dot = std::max(-1.0, std::min(1.0, dot));
    }
    double theta = std::acos(dot);

    double m_inv_sq = p4_pi0[0]*p4_pi0[0] - (pi_px*pi_px + pi_py*pi_py + pi_pz*pi_pz);
    double m_inv = (m_inv_sq > 0.0) ? std::sqrt(m_inv_sq) : 0.0;
    double nu = e_beam - Ee;

    if (m_inv == 0.0) {
        if (verbose) std::printf("[Avnish-from-det] m_inv==0 -> returning 0\n");
        return 0.0;
    }

    double denom = (e1 + e2);
    double last_term = 0.0;
    if (denom != 0.0) last_term = q_mag * std::cos(theta) * pi_mag / denom;
    else if (verbose) std::printf("[Avnish-from-det] e1+e2==0 -> last_term=0\n");

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
    std::array<double,4> p4_pi0 = {0,0,0,0};
    build_pi0_4vec_from_detector(e1,x1,y1, e2,x2,y2, z_nps, theta_nps_deg, p4_pi0);

    double corr_fac = invariant_missing_mass_correction_fac_avnish_from_detector(
        miss_mass_sq_pr, e_beam, Ee, px_e, py_e, pz_e, e1, e2, x1, y1, x2, y2,
        z_nps, theta_nps_deg, M_proton, M_pion0, verbose);

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

#endif // NPS_MMISS_H
