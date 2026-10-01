#ifndef NPS_PI0_RESPONSE_CONVENTIONS_H
#define NPS_PI0_RESPONSE_CONVENTIONS_H

// Pure arithmetic shared by the PARTONS adapter and its regression test.
// No ROOT or PARTONS dependency is needed to check the angular normalization.
#include <array>
#include <cmath>
#include <limits>

namespace nps_pi0_conventions {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double nb_per_gev2_to_ub_per_mev2 = 1.0e-9;

// Hand virtual-photon flux for the electron variables (xB,Q2):
//   Gamma = alpha (W2-M2)/(8*pi*E2*M2*xB2*(1-epsilon)).
// Q2 and W2 are in GeV2; E and M are in GeV.  Multiplying Gamma
// by d2sigma_gamma/(dt dphi) gives d4sigma_e/(dxB dQ2 dt dphi).
// PARTONS DVMPProcessGK06 uses Gamma/(2*pi) times a response bracket;
// its additional target-azimuth divisor is canceled by UUUMinus itself.
inline double hand_flux_xbq2(double q2, double xb, double ebeam,
                            double epsilon, double mass, double alpha) {
    if (!(std::isfinite(q2) && q2 > 0.0 && std::isfinite(xb) && xb > 0.0 && xb < 1.0 &&
          std::isfinite(ebeam) && ebeam > 0.0 && std::isfinite(epsilon) &&
          epsilon > 0.0 && epsilon < 1.0 && std::isfinite(mass) && mass > 0.0 &&
          std::isfinite(alpha) && alpha > 0.0))
        return std::numeric_limits<double>::quiet_NaN();
    const double w2_minus_m2 = q2 * (1.0 / xb - 1.0);
    return alpha * w2_minus_m2 /
        (8.0 * pi * ebeam * ebeam * mass * mass * xb * xb * (1.0 - epsilon));
}

struct Responses {
    // Phi-integrated U=T+epsilon*L and interference responses, in ub/MeV2.
    // A rejected projection retains NaNs so it cannot resemble a zero model.
    bool valid = false;
    double U = std::numeric_limits<double>::quiet_NaN();
    double LT = std::numeric_limits<double>::quiet_NaN();
    double TT = std::numeric_limits<double>::quiet_NaN();
};

// Input: native four-fold electron observables in nb, ordered phi=0,pi/2,pi.
// The extractor convention is
//   d2sigma_gamma/(dt dphi) = [U + sqrt(2*eps*(1+eps))*LT*cos(phi)
//                                 + eps*TT*cos(2*phi)]/(2*pi).
// First divide by the FULL Hand flux; then the 2*pi below converts Fourier
// coefficients into the phi-integrated responses. Dividing instead by the
// native PARTONS prefactor Gamma/(2*pi) would count 2*pi twice.
// For d/dt, nb/GeV2 -> ub/MeV2 contributes 10^-3 * 10^-6 = 10^-9.
inline Responses project_electron_observables(const std::array<double, 3>& electron_nb,
                                               double full_hand_flux, double epsilon) {
    Responses out;
    if (!(std::isfinite(full_hand_flux) && full_hand_flux > 0.0 &&
          std::isfinite(epsilon) && epsilon > 0.0 && epsilon < 1.0)) return out;
    for (double value : electron_nb) if (!std::isfinite(value)) return out;
    const double scale = nb_per_gev2_to_ub_per_mev2 / full_hand_flux;
    const double f0 = electron_nb[0] * scale;
    const double f90 = electron_nb[1] * scale;
    const double f180 = electron_nb[2] * scale;
    const double c0 = (f0 + f180 + 2.0 * f90) / 4.0;
    const double c1 = (f0 - f180) / 2.0;
    const double c2 = (f0 + f180 - 2.0 * f90) / 4.0;
    out.U = 2.0 * pi * c0;
    out.LT = 2.0 * pi * c1 / std::sqrt(2.0 * epsilon * (1.0 + epsilon));
    out.TT = 2.0 * pi * c2 / epsilon;
    out.valid = std::isfinite(out.U) && out.U > 0.0 &&
                std::isfinite(out.LT) && std::isfinite(out.TT);
    return out;
}

} // namespace nps_pi0_conventions
#endif
