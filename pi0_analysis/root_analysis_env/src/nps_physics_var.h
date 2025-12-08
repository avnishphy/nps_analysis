// ============================================================================
// File: nps_physics_var.h
// Purpose: compute physics variables (t, pt, Q2, W, nu, theta, phi) for events
// Author: ChatGPT, adapted to user's codebase (uses nps_helper.h)
// ============================================================================

#ifndef NPS_PHYSICS_VAR_H
#define NPS_PHYSICS_VAR_H

// reuse shared helpers / constants
#include "nps_helper.h"

#include <cmath>
#include <cstdio>
#include <array>

namespace nps {

// return-type: commonly-used physics variables
struct PhysicsVars {
    double t    = 0.0;    // Mandelstam t (GeV^2) -- often negative
    double pt   = 0.0;    // transverse momentum of pi0 relative to virtual photon (GeV/c)
    double Q2   = 0.0;    // photon virtuality (GeV^2) = -q^2
    double W    = 0.0;    // invariant mass of gamma* + proton (GeV)
    double nu   = 0.0;    // energy transfer (GeV)
    double theta= 0.0;    // polar angle between pi0 and virtual photon (radians)
    double phi  = 0.0;    // azimuthal angle (radians), lepton-hadron plane angle
};

// small Vec3 utilities (useful local lambdas)
inline std::array<double,3> cross3(const std::array<double,3> &a, const std::array<double,3> &b) {
    return { a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] };
}

inline double norm3(const std::array<double,3> &a) {
    return std::sqrt(std::max(0.0, dot3(a,a)));
}

// ---------------------------------------------------------------------------
// Build pi0 four-vector in HALL coordinates from detector photon inputs,
// re-using photon4vector(...) and rotate_y(...) from nps_helper.h.
// Returns a Vec4-like std::array<double,4>: {E, px, py, pz}.
// ---------------------------------------------------------------------------
inline std::array<double,4> build_pi0_4vec_from_detector(
        double e1, double x1, double y1,
        double e2, double x2, double y2,
        double z_nps = kDefaultZ_NPS_cm,
        double theta_nps_deg = -17.51)
{
    // photon4vector returns {E, px, py, pz} in NPS coords (assumed)
    auto ph1 = photon4vector(e1, x1, y1, z_nps);
    auto ph2 = photon4vector(e2, x2, y2, z_nps);

    // rotate each photon's momentum from NPS -> hall coords
    Vec3 p1_nps = { ph1[1], ph1[2], ph1[3] };
    Vec3 p2_nps = { ph2[1], ph2[2], ph2[3] };
    Vec3 p1_hall = rotate_y(p1_nps, theta_nps_deg);
    Vec3 p2_hall = rotate_y(p2_nps, theta_nps_deg);

    // build four-vectors in hall frame
    std::array<double,4> p4_ph1 = { ph1[0], p1_hall[0], p1_hall[1], p1_hall[2] };
    std::array<double,4> p4_ph2 = { ph2[0], p2_hall[0], p2_hall[1], p2_hall[2] };

    // pi0 = ph1 + ph2
    std::array<double,4> p4_pi = add4(p4_ph1, p4_ph2); // uses helper's add4

    return p4_pi;
}

// ---------------------------------------------------------------------------
// Build missing-proton four-vector (p_miss = p_in - p_out) using the same
// inputs as missing_mass_proton_pi0. Returns p4_miss = {E, px, py, pz}.
// ---------------------------------------------------------------------------
inline std::array<double,4> missing_p4_from_detector(
        double Ebeam,
        double Ee, double px_e, double py_e, double pz_e,
        double e1, double x1, double y1,
        double e2, double x2, double y2,
        double z_nps = kDefaultZ_NPS_cm,
        double theta_nps_deg = -17.51)
{
    // incident electron 4-vector (beam along +z)
    const double pbeam = std::sqrt(std::max(0.0, Ebeam*Ebeam - kElectronMass_GeV*kElectronMass_GeV));
    std::array<double,4> p4_ein = { Ebeam, 0.0, 0.0, pbeam };

    // proton at rest
    std::array<double,4> p4_pin = { kProtonMass_GeV, 0.0, 0.0, 0.0 };

    // scattered electron 4-vector (use provided components)
    std::array<double,4> p4_eout = { Ee, px_e, py_e, pz_e };

    // build pi0 4-vector in hall coords
    std::array<double,4> p4_pi = build_pi0_4vec_from_detector(e1, x1, y1, e2, x2, y2, z_nps, theta_nps_deg);

    // p4_in and p4_out
    std::array<double,4> p4_in = add4(p4_ein, p4_pin);
    std::array<double,4> p4_out = add4(p4_eout, p4_pi);

    // missing 4-vector
    std::array<double,4> p4_miss = sub4(p4_in, p4_out);

    return p4_miss;
}

// ---------------------------------------------------------------------------
// Compute physics variables from detector inputs.
// Inputs are identical to missing_p4_from_detector (beam, scattered-e, two photons).
// Returns PhysicsVars with t, pt, Q2, W, nu, theta, phi. If a quantity is not
// computable (zero norms etc.) it will be returned as 0 and a debug print may
// be produced if verbose==true.
// ---------------------------------------------------------------------------
inline PhysicsVars compute_physics_vars_from_detector(
        double Ebeam,
        double Ee, double px_e, double py_e, double pz_e,
        double e1, double x1, double y1,
        double e2, double x2, double y2,
        double z_nps = kDefaultZ_NPS_cm,
        double theta_nps_deg = -17.51,
        bool verbose = false)
{
    PhysicsVars out;

    // build four-vectors we need
    const double pbeam = std::sqrt(std::max(0.0, Ebeam*Ebeam - kElectronMass_GeV*kElectronMass_GeV));
    std::array<double,4> p4_ein = { Ebeam, 0.0, 0.0, pbeam };
    std::array<double,4> p4_eout = { Ee, px_e, py_e, pz_e };
    std::array<double,4> p4_pin = { kProtonMass_GeV, 0.0, 0.0, 0.0 };

    // virtual photon q = k - k'
    std::array<double,4> q4 = sub4(p4_ein, p4_eout);
    double q0 = q4[0];
    std::array<double,3> q3 = { q4[1], q4[2], q4[3] };
    double q_mag = norm3(q3);

    // Q2 = - q^2 = -(q0^2 - |q|^2)
    double q2 = q0*q0 - q_mag*q_mag;
    out.Q2 = - q2;
    if (out.Q2 < 0 && std::fabs(out.Q2) < 1e-12) out.Q2 = 0.0; // numerical guard

    // nu = energy transfer
    out.nu = q0;

    // W^2 = (p + q)^2 = M_p^2 + 2*M_p*nu - Q2
    double W2 = kProtonMass_GeV*kProtonMass_GeV + 2.0*kProtonMass_GeV * out.nu - out.Q2;
    out.W = (W2 > 0.0) ? std::sqrt(W2) : 0.0;

    // pi0 4-vector in hall coordinates
    std::array<double,4> p4_pi = build_pi0_4vec_from_detector(e1,x1,y1,e2,x2,y2,z_nps,theta_nps_deg);
    std::array<double,3> p3_pi = { p4_pi[1], p4_pi[2], p4_pi[3] };
    double p_pi_mag = norm3(p3_pi);

    // missing proton 4-vector and recoil (p_recoil = p4_miss)
    std::array<double,4> p4_miss = missing_p4_from_detector(Ebeam, Ee, px_e, py_e, pz_e,
                                                            e1, x1, y1, e2, x2, y2,
                                                            z_nps, theta_nps_deg);

    // t = (p_initial_proton - p_recoil)^2
    std::array<double,4> diff_pr = sub4(p4_pin, p4_miss);
    double t4 = mass2_4vec(diff_pr); // this returns (E^2 - p^2)
    out.t = t4;

    // theta: angle between p_pi and q (use three-vectors)
    if (q_mag > 0.0 && p_pi_mag > 0.0) {
        double cos_theta = dot3(q3, p3_pi) / (q_mag * p_pi_mag);
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
        out.theta = std::acos(cos_theta);
    } else {
        out.theta = 0.0;
        if (verbose) std::printf("[compute_physics_vars] Warning: zero q_mag or p_pi_mag -> theta set to 0\n");
    }

    // pt: transverse momentum of pi0 wrt q (|p_pi - (p_pi·q_hat) q_hat|)
    if (q_mag > 0.0) {
        std::array<double,3> q_hat = { q3[0]/q_mag, q3[1]/q_mag, q3[2]/q_mag };
        double p_par = dot3(p3_pi, q_hat);
        std::array<double,3> p_par_vec = { p_par * q_hat[0], p_par * q_hat[1], p_par * q_hat[2] };
        std::array<double,3> pt_vec = { p3_pi[0] - p_par_vec[0], p3_pi[1] - p_par_vec[1], p3_pi[2] - p_par_vec[2] };
        out.pt = norm3(pt_vec);
    } else {
        out.pt = 0.0;
        if (verbose) std::printf("[compute_physics_vars] Warning: q_mag == 0 -> pt set to 0\n");
    }

    // phi: azimuthal angle between lepton plane and hadron plane.
    // Use standard construction:
    //   n_lep = k x k'
    //   n_had = q x p_pi
    //   phi = atan2( ( q_hat · (n_lep x n_had) ), ( n_lep · n_had ) )
    std::array<double,3> k3 = { p4_ein[1], p4_ein[2], p4_ein[3] };
    std::array<double,3> k3p = { p4_eout[1], p4_eout[2], p4_eout[3] };
    std::array<double,3> n_lep = cross3(k3, k3p);
    std::array<double,3> n_had = cross3(q3, p3_pi);

    double nlep_norm = norm3(n_lep);
    double nhad_norm = norm3(n_had);

    if (nlep_norm > 0.0 && nhad_norm > 0.0 && q_mag > 0.0) {
        // normalized normals
        std::array<double,3> nlep_hat = { n_lep[0]/nlep_norm, n_lep[1]/nlep_norm, n_lep[2]/nlep_norm };
        std::array<double,3> nhad_hat = { n_had[0]/nhad_norm, n_had[1]/nhad_norm, n_had[2]/nhad_norm };

        std::array<double,3> nlep_x_nhad = cross3(nlep_hat, nhad_hat);
        // q_hat:
        std::array<double,3> q_hat = { 0.0, 0.0, 0.0 };
        if (q_mag > 0.0) { q_hat = { q3[0]/q_mag, q3[1]/q_mag, q3[2]/q_mag }; }

        double num = dot3(q_hat, nlep_x_nhad);
        double den = dot3(nlep_hat, nhad_hat);

        out.phi = std::atan2(num, den); // range (-pi, pi)
        if (out.phi < 0); // keep as [-pi,pi] (user can map to 0-2pi if desired)
    } else {
        out.phi = 0.0;
        if (verbose) std::printf("[compute_physics_vars] Warning: degenerate plane normals -> phi set to 0\n");
    }

    // debug print
    if (verbose) {
        std::printf("[PhysicsVars DEBUG] Q2=%.6g nu=%.6g W=%.6g t=%.6g pt=%.6g theta(deg)=%.3g phi(deg)=%.3g\n",
                    out.Q2, out.nu, out.W, out.t, out.pt, out.theta*180.0/M_PI, out.phi*180.0/M_PI);
    }

    return out;
}

} // namespace nps

#endif // NPS_PHYSICS_VAR_H
