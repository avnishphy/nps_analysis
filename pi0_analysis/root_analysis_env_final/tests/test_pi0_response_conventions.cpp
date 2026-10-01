// Standalone angular-normalization regression; no ROOT/PARTONS installation.
// Build/run from this directory:
//   c++ -std=c++17 -O2 -Wall -Wextra -pedantic test_pi0_response_conventions.cpp -o /tmp/test_pi0_response_conventions
//   /tmp/test_pi0_response_conventions
#include "../src/xsec_extract/pi0_response_conventions.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>

namespace {
void require(bool condition, const char* label) {
    if (!condition) { std::cerr << "FAIL: " << label << '\n'; std::exit(1); }
}
void close(double measured, double expected, const char* label) {
    require(std::isfinite(measured) &&
            std::abs(measured - expected) <= 2.0e-13 * std::max(1.0e-15, std::abs(expected)), label);
}
}

int main() {
    using namespace nps_pi0_conventions;
    const double q2 = 4.2, xb = 0.60, ebeam = 10.6, mass = 0.9382720813;
    const double alpha = 1.0 / 137.035999084;
    // Construct electron samples from the independently written native
    // PARTONS prefactor, not from hand_flux_xbq2, to catch an extra/missing 2*pi.
    // U,LT,TT here have nb/GeV2 units; include both interference signs and
    // several polarization values so cos(phi), cos(2phi), and epsilon factors
    // cannot silently interchange. These responses remain positive at all phi.
    for (double epsilon : {0.15, 0.65, 0.92}) {
        for (double sign : {-1.0, 1.0}) {
            const double U = 31.0, LT = sign * 2.3, TT = -sign * 4.7;
            const double w2 = mass * mass + q2 * (1.0 / xb - 1.0);
            const double native_prefactor = alpha * (w2 - mass * mass) /
                (16.0 * pi * pi * ebeam * ebeam * mass * mass * q2 * (1.0 - epsilon)) *
                q2 / (xb * xb);
            const double klt = std::sqrt(2.0 * epsilon * (1.0 + epsilon));
            std::array<double, 3> samples{};
            const std::array<double, 3> angles{0.0, pi / 2.0, pi};
            for (std::size_t i = 0; i < samples.size(); ++i)
                samples[i] = native_prefactor *
                    (U + klt * LT * std::cos(angles[i]) + epsilon * TT * std::cos(2.0 * angles[i]));
            const double flux = hand_flux_xbq2(q2, xb, ebeam, epsilon, mass, alpha);
            close(flux, 2.0 * pi * native_prefactor, "Hand flux versus native prefactor");
            const auto responses = project_electron_observables(samples, flux, epsilon);
            require(responses.valid, "finite physical model accepted");
            close(responses.U, U * 1.0e-9, "U normalization and units");
            close(responses.LT, LT * 1.0e-9, "LT sign and normalization");
            close(responses.TT, TT * 1.0e-9, "TT sign and normalization");
            // Test the recovered angular curve at points NOT used to project.
            // Integrating it numerically must give U, with zero LT/TT average.
            double integral = 0.0;
            constexpr int nphi = 128;
            for (int i = 0; i < nphi; ++i) {
                const double phi = (i + 0.5) * 2.0 * pi / nphi;
                const double recovered = (responses.U + klt * responses.LT * std::cos(phi) +
                    epsilon * responses.TT * std::cos(2.0 * phi)) / (2.0 * pi);
                const double expected = native_prefactor / flux * 1.0e-9 *
                    (U + klt * LT * std::cos(phi) + epsilon * TT * std::cos(2.0 * phi));
                close(recovered, expected, "independent angular reconstruction");
                integral += recovered * 2.0 * pi / nphi;
            }
            close(integral, U * 1.0e-9, "full azimuth integral");
        }
    }
    // Invalid/model-domain inputs must not be advertised as usable predictions.
    require(!project_electron_observables({0.0, 0.0, 0.0}, 1.0, 0.5).valid, "GK zero rejected");
    require(!project_electron_observables({1.0, 1.0, 1.0}, 0.0, 0.5).valid, "zero flux rejected");
    require(!project_electron_observables({1.0, 1.0, 1.0}, 1.0, 0.0).valid, "singular LT/TT rejected");
    require(!project_electron_observables({NAN, 1.0, 1.0}, 1.0, 0.5).valid, "nonfinite sample rejected");
    require(!std::isfinite(hand_flux_xbq2(q2, 1.2, ebeam, 0.5, mass, alpha)), "invalid xB rejected");
    std::cout << "PASS: GK flux, Fourier convention, units, signs, angular reconstruction, invalid inputs\n";
}
