// // ============================================================================
// // File: nps_physics_var.h
// // Purpose: compute physics variables (t, pt, Q2, W, nu, theta, phi) for events
// //          Extended: adds s, xB, tmin, z
// // Author: ChatGPT, adapted to user's codebase (uses nps_helper.h)
// // ============================================================================

// #ifndef NPS_PHYSICS_VAR_H
// #define NPS_PHYSICS_VAR_H

// // reuse shared helpers / constants
// #include "nps_helper.h"

// #include <cmath>
// #include <cstdio>
// #include <array>
// #include <algorithm>

// namespace nps {

// // return-type: commonly-used physics variables
// struct PhysicsVars {
//     double t     = 0.0;    // Mandelstam t (GeV^2) -- often negative
//     double pt    = 0.0;    // transverse momentum of pi0 relative to virtual photon (GeV/c)
//     double Q2    = 0.0;    // photon virtuality (GeV^2) = -q^2
//     double W     = 0.0;    // invariant mass of gamma* + proton (GeV)
//     double nu    = 0.0;    // energy transfer (GeV)
//     double theta = 0.0;    // polar angle between pi0 and virtual photon (radians)
//     double phi   = 0.0;    // azimuthal angle (radians), lepton-hadron plane angle

//     // NEW variables
//     double s     = 0.0;    // invariant s = (k + p)^2 (GeV^2)
//     double xB    = 0.0;    // Bjorken x
//     double tmin  = 0.0;    // minimal kinematic t (GeV^2)
//     double z     = 0.0;    // z = E_pi / nu (target rest frame)
// };

// // small Vec3 utilities (useful local lambdas)
// inline std::array<double,3> cross3(const std::array<double,3> &a, const std::array<double,3> &b) {
//     return { a[1]*b[2] - a[2]*b[1],
//              a[2]*b[0] - a[0]*b[2],
//              a[0]*b[1] - a[1]*b[0] };
// }

// inline double norm3(const std::array<double,3> &a) {
//     return std::sqrt(std::max(0.0, dot3(a,a)));
// }

// // ---------------------------------------------------------------------------
// // Build pi0 four-vector in HALL coordinates from detector photon inputs,
// // re-using photon4vector(...) and rotate_y(...) from nps_helper.h.
// // Returns a Vec4-like std::array<double,4>: {E, px, py, pz}.
// // ---------------------------------------------------------------------------
// inline std::array<double,4> build_pi0_4vec_from_detector(
//         double e1, double x1, double y1,
//         double e2, double x2, double y2,
//         double z_nps = kDefaultZ_NPS_cm,
//         double theta_nps_deg = -17.51)
// {
//     // photon4vector returns {E, px, py, pz} in NPS coords (assumed)
//     auto ph1 = photon4vector(e1, x1, y1, z_nps);
//     auto ph2 = photon4vector(e2, x2, y2, z_nps);

//     // rotate each photon's momentum from NPS -> hall coords
//     Vec3 p1_nps = { ph1[1], ph1[2], ph1[3] };
//     Vec3 p2_nps = { ph2[1], ph2[2], ph2[3] };
//     Vec3 p1_hall = rotate_y(p1_nps, theta_nps_deg);
//     Vec3 p2_hall = rotate_y(p2_nps, theta_nps_deg);

//     // build four-vectors in hall frame
//     std::array<double,4> p4_ph1 = { ph1[0], p1_hall[0], p1_hall[1], p1_hall[2] };
//     std::array<double,4> p4_ph2 = { ph2[0], p2_hall[0], p2_hall[1], p2_hall[2] };

//     // pi0 = ph1 + ph2
//     std::array<double,4> p4_pi = add4(p4_ph1, p4_ph2); // uses helper's add4

//     return p4_pi;
// }

// // ---------------------------------------------------------------------------
// // Build missing-proton four-vector (p_miss = p_in - p_out) using the same
// // inputs as missing_mass_proton_pi0. Returns p4_miss = {E, px, py, pz}.
// // ---------------------------------------------------------------------------
// inline std::array<double,4> missing_p4_from_detector(
//         double Ebeam,
//         double Ee, double px_e, double py_e, double pz_e,
//         double e1, double x1, double y1,
//         double e2, double x2, double y2,
//         double z_nps = kDefaultZ_NPS_cm,
//         double theta_nps_deg = -17.51)
// {
//     // incident electron 4-vector (beam along +z)
//     const double pbeam = std::sqrt(std::max(0.0, Ebeam*Ebeam - kElectronMass_GeV*kElectronMass_GeV));
//     std::array<double,4> p4_ein = { Ebeam, 0.0, 0.0, pbeam };

//     // proton at rest
//     std::array<double,4> p4_pin = { kProtonMass_GeV, 0.0, 0.0, 0.0 };

//     // scattered electron 4-vector (use provided components)
//     std::array<double,4> p4_eout = { Ee, px_e, py_e, pz_e };

//     // build pi0 4-vector in hall coords
//     std::array<double,4> p4_pi = build_pi0_4vec_from_detector(e1, x1, y1, e2, x2, y2, z_nps, theta_nps_deg);

//     // p4_in and p4_out
//     std::array<double,4> p4_in = add4(p4_ein, p4_pin);
//     std::array<double,4> p4_out = add4(p4_eout, p4_pi);

//     // missing 4-vector
//     std::array<double,4> p4_miss = sub4(p4_in, p4_out);

//     return p4_miss;
// }

// // Helper: Kallen lambda: sqrt((a-(b+c))^2 - 4 b c)
// inline double kallen_lambda(double a, double b, double c) {
//     double arg = (a - (b + c))*(a - (b + c)) - 4.0*b*c;
//     return (arg > 0.0) ? std::sqrt(arg) : 0.0;
// }

// // ---------------------------------------------------------------------------
// // Compute physics variables from detector inputs.
// // Inputs are identical to missing_p4_from_detector (beam, scattered-e, two photons).
// // Returns PhysicsVars with t, pt, Q2, W, nu, theta, phi, and new ones s, xB, tmin, z.
// // If a quantity is not computable (zero norms etc.) it will be returned as 0.
// // ---------------------------------------------------------------------------
// inline PhysicsVars compute_physics_vars_from_detector(
//         double Ebeam,
//         double Ee, double px_e, double py_e, double pz_e,
//         double e1, double x1, double y1,
//         double e2, double x2, double y2,
//         double z_nps = kDefaultZ_NPS_cm,
//         double theta_nps_deg = -17.51,
//         bool verbose = false)
// {
//     PhysicsVars out;

//     // build four-vectors we need
//     const double pbeam = std::sqrt(std::max(0.0, Ebeam*Ebeam - kElectronMass_GeV*kElectronMass_GeV));
//     std::array<double,4> p4_ein = { Ebeam, 0.0, 0.0, pbeam };
//     std::array<double,4> p4_eout = { Ee, px_e, py_e, pz_e };
//     std::array<double,4> p4_pin = { kProtonMass_GeV, 0.0, 0.0, 0.0 };

//     // s = (k + p)^2 (invariant)
//     std::array<double,4> p4_initial_sum = add4(p4_ein, p4_pin);
//     out.s = mass2_4vec(p4_initial_sum); // GeV^2

//     // virtual photon q = k - k'
//     std::array<double,4> q4 = sub4(p4_ein, p4_eout);
//     double q0 = q4[0];
//     std::array<double,3> q3 = { q4[1], q4[2], q4[3] };
//     double q_mag = norm3(q3);

//     // Q2 = - q^2 = -(q0^2 - |q|^2)
//     double q2 = q0*q0 - q_mag*q_mag;
//     out.Q2 = - q2;
//     if (out.Q2 < 0 && std::fabs(out.Q2) < 1e-12) out.Q2 = 0.0; // numerical guard

//     // nu = energy transfer
//     out.nu = q0;

//     // xB = Q2 / (2 M_p nu)
//     double denom_xb = 2.0 * kProtonMass_GeV * out.nu;
//     if (denom_xb != 0.0) out.xB = out.Q2 / denom_xb;
//     else out.xB = 0.0;

//     // W^2 = (p + q)^2 = M_p^2 + 2*M_p*nu - Q2
//     double W2 = kProtonMass_GeV*kProtonMass_GeV + 2.0*kProtonMass_GeV * out.nu - out.Q2;
//     out.W = (W2 > 0.0) ? std::sqrt(W2) : 0.0;

//     // pi0 4-vector in hall coordinates
//     std::array<double,4> p4_pi = build_pi0_4vec_from_detector(e1,x1,y1,e2,x2,y2,z_nps,theta_nps_deg);
//     std::array<double,3> p3_pi = { p4_pi[1], p4_pi[2], p4_pi[3] };
//     double p_pi_mag = norm3(p3_pi);

//     // z = E_pi_lab / nu (if target at rest, P·p_pi / P·q reduces to E_pi/nu)
//     if (out.nu != 0.0) out.z = p4_pi[0] / out.nu;
//     else out.z = 0.0;

//     // missing proton 4-vector and recoil (p_recoil = p4_miss)
//     std::array<double,4> p4_miss = missing_p4_from_detector(Ebeam, Ee, px_e, py_e, pz_e,
//                                                             e1, x1, y1, e2, x2, y2,
//                                                             z_nps, theta_nps_deg);

//     // t = (p_initial_proton - p_recoil)^2
//     std::array<double,4> diff_pr = sub4(p4_pin, p4_miss);
//     double t4 = mass2_4vec(diff_pr); // this returns (E^2 - p^2)
//     out.t = t4;

//     // theta: angle between p_pi and q (use three-vectors)
//     if (q_mag > 0.0 && p_pi_mag > 0.0) {
//         double cos_theta = dot3(q3, p3_pi) / (q_mag * p_pi_mag);
//         cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
//         out.theta = std::acos(cos_theta);
//     } else {
//         out.theta = 0.0;
//         if (verbose) std::printf("[compute_physics_vars] Warning: zero q_mag or p_pi_mag -> theta set to 0\n");
//     }

//     // pt: transverse momentum of pi0 wrt q (|p_pi - (p_pi·q_hat) q_hat|)
//     if (q_mag > 0.0) {
//         std::array<double,3> q_hat = { q3[0]/q_mag, q3[1]/q_mag, q3[2]/q_mag };
//         double p_par = dot3(p3_pi, q_hat);
//         std::array<double,3> p_par_vec = { p_par * q_hat[0], p_par * q_hat[1], p_par * q_hat[2] };
//         std::array<double,3> pt_vec = { p3_pi[0] - p_par_vec[0], p3_pi[1] - p_par_vec[1], p3_pi[2] - p_par_vec[2] };
//         out.pt = norm3(pt_vec);
//     } else {
//         out.pt = 0.0;
//         if (verbose) std::printf("[compute_physics_vars] Warning: q_mag == 0 -> pt set to 0\n");
//     }

//     // phi: azimuthal angle between lepton plane and hadron plane.
//     // Use standard construction:
//     //   n_lep = k x k'
//     //   n_had = q x p_pi
//     //   phi = atan2( ( q_hat · (n_lep x n_had) ), ( n_lep · n_had ) )
//     std::array<double,3> k3 = { p4_ein[1], p4_ein[2], p4_ein[3] };
//     std::array<double,3> k3p = { p4_eout[1], p4_eout[2], p4_eout[3] };
//     std::array<double,3> n_lep = cross3(k3, k3p);
//     std::array<double,3> n_had = cross3(q3, p3_pi);

//     double nlep_norm = norm3(n_lep);
//     double nhad_norm = norm3(n_had);

//     if (nlep_norm > 0.0 && nhad_norm > 0.0 && q_mag > 0.0) {
//         // normalized normals
//         std::array<double,3> nlep_hat = { n_lep[0]/nlep_norm, n_lep[1]/nlep_norm, n_lep[2]/nlep_norm };
//         std::array<double,3> nhad_hat = { n_had[0]/nhad_norm, n_had[1]/nhad_norm, n_had[2]/nhad_norm };

//         std::array<double,3> nlep_x_nhad = cross3(nlep_hat, nhad_hat);
//         // q_hat:
//         std::array<double,3> q_hat = { 0.0, 0.0, 0.0 };
//         if (q_mag > 0.0) { q_hat = { q3[0]/q_mag, q3[1]/q_mag, q3[2]/q_mag }; }

//         double num = dot3(q_hat, nlep_x_nhad);
//         double den = dot3(nlep_hat, nhad_hat);

//         out.phi = std::atan2(num, den); // range (-pi, pi)
//         if (out.phi < 0); // keep as [-pi,pi] (user can map to 0-2pi if desired)
//     } else {
//         out.phi = 0.0;
//         if (verbose) std::printf("[compute_physics_vars] Warning: degenerate plane normals -> phi set to 0\n");
//     }

//     // ---------- compute tmin for exclusive gamma* p -> p' pi0 ----------
//     // Use CM (gamma*-p) two-body kinematics:
//     //   p* = lambda(W^2, Mp^2, Mpi^2) / (2 W)
//     //   E_p' = (W^2 + Mp^2 - Mpi^2) / (2 W)
//     //   tmin = 2 * Mp * (Mp - E_p')   (<= 0)
//     if (out.W > 0.0) {
//         const double Wval = out.W;
//         const double Mp = kProtonMass_GeV;
//         const double Mpi = kPi0Mass_GeV;
//         double lambda = kallen_lambda(Wval*Wval, Mp*Mp, Mpi*Mpi);
//         // p_star not explicitly needed for tmin via E_p', but we keep lambda in case debug needed
//         double E_pprime = (Wval*Wval + Mp*Mp - Mpi*Mpi) / (2.0 * Wval);
//         out.tmin = 2.0 * Mp * (Mp - E_pprime);
//     } else {
//         out.tmin = 0.0;
//         if (verbose) std::printf("[compute_physics_vars] Warning: W <= 0 -> tmin set to 0\n");
//     }

//     // debug print
//     if (verbose) {
//         std::printf("[PhysicsVars DEBUG] s=%.6g Q2=%.6g xB=%.6g nu=%.6g W=%.6g t=%.6g tmin=%.6g pt=%.6g z=%.6g theta(deg)=%.3g phi(deg)=%.3g\n",
//                     out.s, out.Q2, out.xB, out.nu, out.W, out.t, out.tmin, out.pt, out.z, out.theta*180.0/M_PI, out.phi*180.0/M_PI);
//     }

//     return out;
// }

// } // namespace nps

// #endif // NPS_PHYSICS_VAR_H





// ===============================================================================================================================================================




// ============================================================================
// File: nps_physics_var.h
// Purpose: compute event-level physics variables for exclusive pi0
//          electroproduction:
//
//              e(k) + p(P) -> e'(k') + p'(P') + pi0(p_pi)
//
//          Definitions used here are the standard invariant ones:
//
//              q      = k - k'
//              Q^2    = -q^2
//              nu     = q^0
//              W^2    = (P + q)^2
//              s      = (k + P)^2
//              xB     = Q^2 / (2 P·q) = Q^2 / (2 M_p nu) in the target rest frame
//              t      = (q - p_pi)^2 = (P - P')^2   [canonical exclusive definition]
//              pt     = |p_pi,⊥| with respect to q
//              theta  = angle between p_pi and q
//              phi    = signed azimuthal angle between lepton and hadron planes
//              z      = (P·p_pi)/(P·q) = E_pi/nu in the target rest frame
//              tmin   = forward-limit t for gamma* p -> pi0 p' at fixed W,Q^2
//
// Notes:
//  1) The code assumes the beam travels along +z in the hall frame.
//  2) The NPS photon positions are converted to a pi0 four-vector in the hall
//     frame using the existing helper photon4vector(...) + rotate_y(...).
//  3) The azimuth phi is returned in (-pi, pi] using the standard atan2 form.
//     This is a convention; your analysis can remap it to [0, 2pi) if desired.
//  4) For exclusive kinematics, the canonical t = (q - p_pi)^2 is preferred.
//     It is exactly equal to (P - P')^2 only in the fully exclusive limit.
// ============================================================================

#ifndef NPS_PHYSICS_VAR_H
#define NPS_PHYSICS_VAR_H

#include "nps_helper.h"

#include <cmath>
#include <cstdio>
#include <array>
#include <algorithm>

namespace nps {

// ---------------------------------------------------------------------------
// Container for commonly used kinematic variables.
// All momenta are in GeV, all angles in radians, and all invariants in GeV^2.
// ---------------------------------------------------------------------------
struct PhysicsVars {
    double t     = 0.0;    // Mandelstam t = (q - p_pi)^2, typically negative
    double pt    = 0.0;    // pi0 transverse momentum wrt virtual photon q
    double Q2    = 0.0;    // Q^2 = -q^2
    double W     = 0.0;    // W = sqrt((P + q)^2)
    double nu    = 0.0;    // energy transfer = q^0 in lab
    double theta = 0.0;    // angle between q and p_pi
    double phi   = 0.0;    // signed azimuth between lepton and hadron planes

    // Additional useful scalars
    double s     = 0.0;    // s = (k + P)^2
    double xB    = 0.0;    // Bjorken x = Q^2 / (2 P·q)
    double tmin  = 0.0;    // forward-limit t at fixed W,Q^2 for gamma* p -> pi0 p'
    double z     = 0.0;    // z = (P·p_pi)/(P·q) = E_pi/nu in target rest frame
};

// ---------------------------------------------------------------------------
// Small local helpers.
// ---------------------------------------------------------------------------
inline std::array<double,3> cross3(const std::array<double,3> &a,
                                   const std::array<double,3> &b)
{
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

inline double norm3(const std::array<double,3> &a) {
    return std::sqrt(std::max(0.0, dot3(a, a)));
}

// Standard Källén function sqrt(lambda(a,b,c)) with
// lambda = a^2 + b^2 + c^2 - 2ab - 2ac - 2bc
// written in a numerically stable equivalent form.
inline double kallen_lambda(double a, double b, double c) {
    const double arg = (a - (b + c)) * (a - (b + c)) - 4.0 * b * c;
    return (arg > 0.0) ? std::sqrt(arg) : 0.0;
}

// ---------------------------------------------------------------------------
// Build the pi0 four-vector in the hall frame from the two photon clusters.
//
// Inputs:
//   e1, x1, y1  : photon 1 energy and cluster position in NPS coordinates
//   e2, x2, y2  : photon 2 energy and cluster position in NPS coordinates
//   z_nps       : NPS z-position in cm
//   theta_nps_deg: rotation angle that maps NPS momentum axes into hall axes
//
// Returns:
//   p4_pi = {E, px, py, pz} in the hall frame
//
// Important:
//   The helper photon4vector(...) is assumed to produce a massless photon
//   four-vector in the NPS coordinate system. We rotate only the momentum
//   three-vector into the hall frame; the photon energy is unchanged.
// ---------------------------------------------------------------------------
inline std::array<double,4> build_pi0_4vec_from_detector(
        double e1, double x1, double y1,
        double e2, double x2, double y2,
        double z_nps = kDefaultZ_NPS_cm,
        double theta_nps_deg = -17.51)
{
    // Photon four-vectors in NPS coordinates.
    auto ph1 = photon4vector(e1, x1, y1, z_nps);
    auto ph2 = photon4vector(e2, x2, y2, z_nps);

    // Rotate momentum components from NPS frame into hall frame.
    Vec3 p1_nps = { ph1[1], ph1[2], ph1[3] };
    Vec3 p2_nps = { ph2[1], ph2[2], ph2[3] };

    Vec3 p1_hall = rotate_y(p1_nps, theta_nps_deg);
    Vec3 p2_hall = rotate_y(p2_nps, theta_nps_deg);

    // Reassemble full four-vectors in the hall frame.
    std::array<double,4> p4_ph1 = { ph1[0], p1_hall[0], p1_hall[1], p1_hall[2] };
    std::array<double,4> p4_ph2 = { ph2[0], p2_hall[0], p2_hall[1], p2_hall[2] };

    // pi0 = gamma1 + gamma2
    return add4(p4_ph1, p4_ph2);
}

// ---------------------------------------------------------------------------
// Build the missing four-vector for the exclusive reaction.
//
// p_miss = (k + P) - (k' + p_pi)
//
// For a perfectly exclusive event, p_miss should equal the recoil proton
// four-vector P' up to resolution effects.
//
// This utility is useful for diagnostics, but t should be computed canonically
// from q and p_pi, not from the missing proton, when possible.
// ---------------------------------------------------------------------------
inline std::array<double,4> missing_p4_from_detector(
        double Ebeam,
        double Ee, double px_e, double py_e, double pz_e,
        double e1, double x1, double y1,
        double e2, double x2, double y2,
        double z_nps = kDefaultZ_NPS_cm,
        double theta_nps_deg = -17.51)
{
    // Incident electron four-vector in the hall frame.
    // Beam is taken along +z.
    const double pbeam = std::sqrt(std::max(0.0,
                              Ebeam*Ebeam - kElectronMass_GeV*kElectronMass_GeV));
    const std::array<double,4> p4_ein = { Ebeam, 0.0, 0.0, pbeam };

    // Proton at rest in the lab.
    const std::array<double,4> p4_pin = { kProtonMass_GeV, 0.0, 0.0, 0.0 };

    // Scattered electron four-vector, already expressed in hall coordinates.
    const std::array<double,4> p4_eout = { Ee, px_e, py_e, pz_e };

    // Reconstructed pi0 in hall coordinates.
    const std::array<double,4> p4_pi =
        build_pi0_4vec_from_detector(e1, x1, y1, e2, x2, y2, z_nps, theta_nps_deg);

    const std::array<double,4> p4_in  = add4(p4_ein,  p4_pin);
    const std::array<double,4> p4_out = add4(p4_eout, p4_pi);

    return sub4(p4_in, p4_out);
}

// ---------------------------------------------------------------------------
// Compute physics variables from detector inputs.
//
// The inputs are the same as for missing_p4_from_detector(...):
//   - beam energy
//   - scattered-electron energy and lab-frame momentum components
//   - photon cluster energies and positions
//
// Physics definitions used:
//
//   q = k - k'
//   Q^2 = -q^2 = |\vec q|^2 - q_0^2
//   nu  = q_0
//   s   = (k + P)^2
//   W^2 = (P + q)^2 = M_p^2 + 2 M_p nu - Q^2
//   xB  = Q^2 / (2 M_p nu)   [target rest frame]
//   z   = E_pi / nu          [target rest frame]
//   theta = angle between q and p_pi
//   pt    = component of p_pi transverse to q
//   phi   = signed azimuth between lepton plane and hadron plane
//   t     = (q - p_pi)^2
//
// tmin is computed in the gamma*-p center-of-mass frame at forward pi0
// production (cos(theta*) = 1), which is the standard exclusive-electro-
// production kinematic endpoint.
//
// Returns:
//   PhysicsVars with all variables filled. Non-computable quantities are
//   left as 0.
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

    // -----------------------------------------------------------------------
    // Construct the primary four-vectors in the hall frame.
    // -----------------------------------------------------------------------
    const double pbeam = std::sqrt(std::max(0.0,
                              Ebeam*Ebeam - kElectronMass_GeV*kElectronMass_GeV));

    const std::array<double,4> p4_ein  = { Ebeam, 0.0, 0.0, pbeam };
    const std::array<double,4> p4_eout  = { Ee, px_e, py_e, pz_e };
    const std::array<double,4> p4_pin   = { kProtonMass_GeV, 0.0, 0.0, 0.0 };

    // -----------------------------------------------------------------------
    // Invariant s = (k + P)^2.
    // This is the total invariant mass squared of the initial e-p system.
    // -----------------------------------------------------------------------
    const std::array<double,4> p4_initial_sum = add4(p4_ein, p4_pin);
    out.s = mass2_4vec(p4_initial_sum);

    // -----------------------------------------------------------------------
    // Virtual photon: q = k - k'
    // -----------------------------------------------------------------------
    const std::array<double,4> q4 = sub4(p4_ein, p4_eout);
    const double q0 = q4[0];
    const std::array<double,3> q3 = { q4[1], q4[2], q4[3] };
    const double q_mag = norm3(q3);

    // Q^2 = -q^2 = |\vec q|^2 - q0^2
    const double q2 = q0*q0 - q_mag*q_mag;
    out.Q2 = -q2;
    if (out.Q2 < 0.0 && std::fabs(out.Q2) < 1e-12) out.Q2 = 0.0;

    // Energy transfer.
    out.nu = q0;

    // Bjorken x in the target rest frame.
    const double denom_xb = 2.0 * kProtonMass_GeV * out.nu;
    out.xB = (denom_xb != 0.0) ? (out.Q2 / denom_xb) : 0.0;

    // W^2 = (P + q)^2 = M_p^2 + 2 M_p nu - Q^2
    const double W2 = kProtonMass_GeV*kProtonMass_GeV
                    + 2.0 * kProtonMass_GeV * out.nu
                    - out.Q2;
    out.W = (W2 > 0.0) ? std::sqrt(W2) : 0.0;

    // -----------------------------------------------------------------------
    // Reconstruct the pi0 four-vector in the hall frame.
    // -----------------------------------------------------------------------
    const std::array<double,4> p4_pi =
        build_pi0_4vec_from_detector(e1, x1, y1, e2, x2, y2, z_nps, theta_nps_deg);

    const std::array<double,3> p3_pi = { p4_pi[1], p4_pi[2], p4_pi[3] };
    const double p_pi_mag = norm3(p3_pi);

    // z = (P·p_pi)/(P·q) = E_pi/nu in the target rest frame.
    // This is not a fundamental exclusive observable, but it is useful for
    // diagnostics and for consistency checks against semi-inclusive variables.
    out.z = (out.nu != 0.0) ? (p4_pi[0] / out.nu) : 0.0;

    // -----------------------------------------------------------------------
    // Canonical exclusive Mandelstam t:
    //   t = (q - p_pi)^2
    //
    // This is the invariant momentum transfer to the hadronic final state.
    // In the exactly exclusive limit it equals (P - P')^2.
    // -----------------------------------------------------------------------
    out.t = mass2_4vec(sub4(q4, p4_pi));

    // -----------------------------------------------------------------------
    // theta = angle between the virtual photon q and the reconstructed pi0.
    // -----------------------------------------------------------------------
    if (q_mag > 0.0 && p_pi_mag > 0.0) {
        double cos_theta = dot3(q3, p3_pi) / (q_mag * p_pi_mag);
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
        out.theta = std::acos(cos_theta);
    } else {
        out.theta = 0.0;
        if (verbose) {
            std::printf("[compute_physics_vars] Warning: zero q_mag or p_pi_mag -> theta = 0\n");
        }
    }

    // -----------------------------------------------------------------------
    // pt = transverse momentum of pi0 with respect to q.
    //
    // Construction:
    //   p_pi = p_parallel + p_transverse
    //   p_parallel = (p_pi · q_hat) q_hat
    //   pt = |p_transverse|
    // -----------------------------------------------------------------------
    if (q_mag > 0.0) {
        const std::array<double,3> q_hat = { q3[0]/q_mag, q3[1]/q_mag, q3[2]/q_mag };
        const double p_par = dot3(p3_pi, q_hat);

        const std::array<double,3> p_par_vec = {
            p_par * q_hat[0],
            p_par * q_hat[1],
            p_par * q_hat[2]
        };

        const std::array<double,3> pt_vec = {
            p3_pi[0] - p_par_vec[0],
            p3_pi[1] - p_par_vec[1],
            p3_pi[2] - p_par_vec[2]
        };

        out.pt = norm3(pt_vec);
    } else {
        out.pt = 0.0;
        if (verbose) {
            std::printf("[compute_physics_vars] Warning: q_mag == 0 -> pt = 0\n");
        }
    }

    // -----------------------------------------------------------------------
    // phi: signed azimuthal angle between the lepton plane and hadron plane.
    //
    // Lepton plane normal:
    //   n_lep = k x k'
    //
    // Hadron plane normal:
    //   n_had = q x p_pi
    //
    // Signed angle:
    //   phi = atan2( q_hat · (n_lep x n_had), n_lep · n_had )
    //
    // This returns phi in (-pi, pi].  If your analysis prefers [0, 2pi),
    // remap it downstream with:
    //   if (phi < 0) phi += 2*pi;
    // -----------------------------------------------------------------------
    const std::array<double,3> k3  = { p4_ein[1],  p4_ein[2],  p4_ein[3]  };
    const std::array<double,3> kp3 = { p4_eout[1], p4_eout[2], p4_eout[3] };

    const std::array<double,3> n_lep = cross3(k3, kp3);
    const std::array<double,3> n_had = cross3(q3, p3_pi);

    const double nlep_norm = norm3(n_lep);
    const double nhad_norm = norm3(n_had);

    if (nlep_norm > 0.0 && nhad_norm > 0.0 && q_mag > 0.0) {
        const std::array<double,3> nlep_hat = {
            n_lep[0] / nlep_norm,
            n_lep[1] / nlep_norm,
            n_lep[2] / nlep_norm
        };
        const std::array<double,3> nhad_hat = {
            n_had[0] / nhad_norm,
            n_had[1] / nhad_norm,
            n_had[2] / nhad_norm
        };

        const std::array<double,3> nlep_x_nhad = cross3(nlep_hat, nhad_hat);

        const std::array<double,3> q_hat = {
            q3[0] / q_mag,
            q3[1] / q_mag,
            q3[2] / q_mag
        };

        const double num = dot3(q_hat, nlep_x_nhad);
        const double den = dot3(nlep_hat, nhad_hat);

        out.phi = std::atan2(num, den);

        // Keep the default output in (-pi, pi].  If you want [0, 2pi),
        // perform the remapping outside or add a dedicated flag.
    } else {
        out.phi = 0.0;
        if (verbose) {
            std::printf("[compute_physics_vars] Warning: degenerate plane normals -> phi = 0\n");
        }
    }

    // -----------------------------------------------------------------------
    // tmin for exclusive gamma* p -> pi0 p' at fixed W and Q^2.
    //
    // In the gamma*-p CM frame:
    //   q* = (q0*, 0, 0, |q*|)
    //   p_pi* = (E_pi*, p_pi* sin theta*, 0, p_pi* cos theta*)
    //
    // The forward-limit value corresponds to cos(theta*) = +1:
    //
    //   tmin = m_pi^2 - Q^2 - 2 ( q0* E_pi* - |q*| |p_pi*| )
    //
    // with
    //   q0*    = (W^2 - M^2 - Q^2) / (2W)
    //   |q*|   = sqrt((q0*)^2 + Q^2)
    //   E_pi*  = (W^2 + m_pi^2 - M^2) / (2W)
    //   |p_pi*|= lambda^{1/2}(W^2, M^2, m_pi^2) / (2W)
    //
    // This quantity is only meaningful above threshold:
    //   W > M_p + m_pi
    // -----------------------------------------------------------------------
    if (out.W > (kProtonMass_GeV + kPi0Mass_GeV)) {
        const double W   = out.W;
        const double Mp  = kProtonMass_GeV;
        const double Mpi = kPi0Mass_GeV;

        const double q0_star  = (W*W - Mp*Mp - out.Q2) / (2.0 * W);
        const double q_star   = std::sqrt(std::max(0.0, q0_star*q0_star + out.Q2));

        const double Epi_star  = (W*W + Mpi*Mpi - Mp*Mp) / (2.0 * W);
        const double ppi_star  = kallen_lambda(W*W, Mp*Mp, Mpi*Mpi) / (2.0 * W);

        out.tmin = Mpi*Mpi - out.Q2 - 2.0 * (q0_star * Epi_star - q_star * ppi_star);
    } else {
        out.tmin = 0.0;
        if (verbose) {
            std::printf("[compute_physics_vars] Warning: W below pi0 production threshold -> tmin = 0\n");
        }
    }

    // -----------------------------------------------------------------------
    // Optional debug printout.
    // -----------------------------------------------------------------------
    if (verbose) {
        std::printf(
            "[PhysicsVars DEBUG] s=%.6g Q2=%.6g xB=%.6g nu=%.6g W=%.6g "
            "t=%.6g tmin=%.6g pt=%.6g z=%.6g theta(deg)=%.3g phi(deg)=%.3g\n",
            out.s, out.Q2, out.xB, out.nu, out.W,
            out.t, out.tmin, out.pt, out.z,
            out.theta * 180.0 / M_PI, out.phi * 180.0 / M_PI
        );
    }

    return out;
}

} // namespace nps

#endif // NPS_PHYSICS_VAR_H