#pragma once

#include <cmath>
#include <limits>

namespace nps_xsec {

struct YieldAsymmetry {
    bool valid = false;
    double value = std::numeric_limits<double>::quiet_NaN();
    double error = std::numeric_limits<double>::quiet_NaN();
};

// Diagnostic for two independent weighted yields. This is not a beam-spin
// asymmetry: helicity-specific luminosities and beam polarization are absent.
// Variances are sums of squared event weights, conditional on those weights.
inline YieldAsymmetry helicity_yield_asymmetry(double plus, double minus,
                                             double variance_plus, double variance_minus) {
    YieldAsymmetry result;
    const double sum = plus + minus;
    if (!(std::isfinite(plus) && std::isfinite(minus) && std::isfinite(sum) && sum > 0.0 &&
          std::isfinite(variance_plus) && variance_plus > 0.0 &&
          std::isfinite(variance_minus) && variance_minus > 0.0))
        return result;

    // A=(P-M)/(P+M): dA/dP=2M/(P+M)^2, dA/dM=-2P/(P+M)^2.
    // Propagating only the numerator would neglect its correlation with the
    // denominator. Require weighted entries of both signs; an absent helicity
    // sample must not appear as a measured A=+1 or -1. Signed subtraction
    // yields remain allowed, so the diagnostic is not clipped to [-1,1].
    const double derivative_plus = 2.0 * (minus / sum) / sum;
    const double derivative_minus = -2.0 * (plus / sum) / sum;
    result.value = (plus - minus) / sum;
    result.error = std::hypot(derivative_plus * std::sqrt(variance_plus),
                              derivative_minus * std::sqrt(variance_minus));
    result.valid = std::isfinite(result.value) && std::isfinite(result.error);
    return result;
}

} // namespace nps_xsec
