#pragma once

#include <cmath>
#include <stdexcept>

namespace nps_xsec {

// Virtual-photon polarization at the hard-scattering vertex of an exclusive
// hydrogen SIMC event. The original h10 "epsilon" branch is reconstructed
// epsilon; it is not the vertex epsilon used by SIMC's pion model and sigcm.
//
// Inputs from the matching original exclusive h10 entry:
//   q2i_gev2, wi_gev       : generated Q2i [GeV^2] and Wi [GeV]
//   hsxptari, hsyptari    : original HMS electron x/y slopes (dimensionless
//                           tangents, as stored; do not convert from mrad)
// Other inputs:
//   hms_theta_deg          : electron-arm central angle from that SIMC run's
//                           .hist file [degrees], not a reconstructed angle
//   proton_mass_gev       : proton mass used for the hydrogen W definition [GeV]
//
// For the HMS electron arm, SIMC physics_angles uses phi0=3*pi/2, giving
// cos(theta_e)=(cos(theta0)+hsyptari*sin(theta0)) /
//              sqrt(1+hsxptari^2+hsyptari^2).
// The hydrogen identity W_i^2=M_p^2+2*M_p*nu_i-Q2_i supplies vertex nu_i.
// These are the inputs to SIMC's main%epsilon; no nominal beam energy enters.
// Invalid inputs throw rather than silently producing a response coefficient.
inline double vertex_epsilon_from_exclusive_simc(double q2i_gev2,
                                                 double wi_gev,
                                                 double hsxptari,
                                                 double hsyptari,
                                                 double hms_theta_deg,
                                                 double proton_mass_gev) {
    constexpr double pi = 3.141592653589793238462643383279502884;
    if (!(std::isfinite(q2i_gev2) && q2i_gev2 > 0.0 &&
          std::isfinite(wi_gev) && wi_gev > 0.0 &&
          std::isfinite(hsxptari) && std::isfinite(hsyptari) &&
          std::isfinite(hms_theta_deg) && hms_theta_deg > 0.0 && hms_theta_deg < 180.0 &&
          std::isfinite(proton_mass_gev) && proton_mass_gev > 0.0))
        throw std::domain_error("Invalid exclusive SIMC vertex-epsilon input");

    const double theta0 = hms_theta_deg * pi / 180.0;
    const double cos_theta =
        (std::cos(theta0) + hsyptari * std::sin(theta0)) /
        std::sqrt(1.0 + hsxptari * hsxptari + hsyptari * hsyptari);
    const double nu_gev = (wi_gev * wi_gev + q2i_gev2 -
                           proton_mass_gev * proton_mass_gev) /
                          (2.0 * proton_mass_gev);
    if (!(std::isfinite(cos_theta) && cos_theta > -1.0 && cos_theta < 1.0 &&
          std::isfinite(nu_gev) && nu_gev > 0.0))
        throw std::domain_error("Unphysical exclusive SIMC vertex kinematics for epsilon");

    // tan^2(theta_e/2)=(1-cos(theta_e))/(1+cos(theta_e)).
    const double tan_half_theta_sq = (1.0 - cos_theta) / (1.0 + cos_theta);
    const double epsilon =
        1.0 / (1.0 + 2.0 * (1.0 + nu_gev * nu_gev / q2i_gev2) * tan_half_theta_sq);
    if (!(std::isfinite(epsilon) && epsilon > 0.0 && epsilon < 1.0))
        throw std::domain_error("Unphysical exclusive SIMC vertex epsilon");
    return epsilon;
}

} // namespace nps_xsec
