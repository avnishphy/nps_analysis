#ifndef NPS_SIMC_PI0_MODEL_H
#define NPS_SIMC_PI0_MODEL_H

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

// Port of simc_gfortran_updated/physics_pion.f: sig_param_2021/exclfit
// (inspected 2026-09-19). Only doing_pizero at W >= 2 GeV is supported;
// physics_pion.f uses a separate MAID path below that boundary.
// Inputs: GeV, GeV^2, radians. Output: microbarn/MeV^2/radian, like sigcm.
// Preserve the source's 0.938 mass and 3.1415928 normalization constants.
namespace nps_simc_pi0 {
struct Model {
    double q2 = 0, xb = 0, tprime = 0, t = 0, w = 0, epsilon = 0;
    // Coefficients of 1, cos(phi), cos(2phi), including 1/(2pi).
    std::array<double, 3> angular{};

    double phi_average(double lo, double hi) const {
        if (!std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo))
            throw std::runtime_error("SIMC model: invalid phi bin");
        const double width = hi - lo;
        return angular[0] + angular[1] * (std::sin(hi) - std::sin(lo)) / width
            + angular[2] * (std::sin(2*hi) - std::sin(2*lo)) / (2*width);
    }
};

inline std::array<double, 3> angular_coefficients(double q2, double wsq,
                                                 double t, double theta,
                                                 double eps) {
    if (!(std::isfinite(q2) && q2 > 0 && std::isfinite(wsq) && wsq >= 4 &&
          std::isfinite(t) && std::isfinite(theta) && std::isfinite(eps) &&
          eps >= 0 && eps <= 1))
        throw std::runtime_error("SIMC pi0 model requires physical kinematics and W >= 2 GeV");
    constexpr std::array<double, 17> pp{{
        1.60077, -0.01523, 37.08142, -4.11060, 23.26192, 0.00983,
        0.87073, -5.77115, -271.08678, 0.13766, -0.00855, 0.27885,
        -1.13212, -1.50415, -6.34766, 0.55769, -0.01709}};
    constexpr std::array<double, 17> pm{{
        1.75169, 0.11144, 47.35877, -4.69434, 1.60552, 0.00800,
        0.44194, -2.29188, -41.67194, 0.69475, 0.02527, -0.50178,
        -1.22825, -1.16878, 5.75825, -1.00355, 0.05055}};
    std::array<double, 3> result{};
    const double at = std::abs(t), st = std::sin(theta);
    const double norm = 8.539 / std::pow(wsq - 0.938*0.938, 2)
                        / 2.0 / 3.1415928 / 1.e6;
    for (const auto& p : {pp, pm}) {
        // fpifact=0 for pi0, hence sigL=0 in SIMC's present model.
        const double sigT = p[4]/q2 * std::exp(p[5]*q2*q2)
            / (std::pow(wsq, p[11]) + std::pow(std::sqrt(wsq), p[15]))
            * std::exp(p[13]*at);
        const double sigLT = p[6]/(1+p[9]*q2) * std::exp(p[7]*at)
            * st / std::pow(wsq, p[12]);
        const double sigTT = p[8]/(1+q2) * std::exp(-7*at) * st*st;
        result[0] += 0.5 * norm * sigT;
        result[1] += 0.5 * norm * std::sqrt(2*eps*(1+eps)) * sigLT;
        result[2] += 0.5 * norm * eps * sigTT;
    }
    for (double c : result)
        if (!std::isfinite(c)) throw std::runtime_error("Nonfinite SIMC pi0 model coefficient");
    return result;
}

inline Model at_reference(double q2, double xb, double tprime, double ebeam,
                          double mp, double mpi0) {
    if (!(std::isfinite(q2) && q2 > 0 && std::isfinite(xb) && xb > 0 && xb < 1 &&
          std::isfinite(tprime) && tprime <= 0 && std::isfinite(ebeam) && ebeam > 0 &&
          std::isfinite(mp) && mp > 0 && std::isfinite(mpi0) && mpi0 > 0))
        throw std::runtime_error("Invalid SIMC model reference kinematics");
    Model model;
    model.q2 = q2; model.xb = xb; model.tprime = tprime;
    const double wsq = mp*mp + q2*(1/xb - 1);
    model.w = std::sqrt(wsq);
    if (model.w < 2 || model.w <= mp + mpi0)
        throw std::runtime_error("SIMC pi0 reference W < 2 GeV: MAID evaluation is not implemented");
    const double nu = q2/(2*mp*xb), eprime = ebeam - nu;
    if (!(eprime > 0) || q2 > 4*ebeam*eprime)
        throw std::runtime_error("Unphysical electron kinematics at SIMC reference point");
    const double y = nu/ebeam, z = q2/(4*ebeam*ebeam);
    model.epsilon = (1-y-z)/(1-y+0.5*y*y+z);
    const double egamma = (wsq-mp*mp-q2)/(2*model.w);
    const double pgamma = std::sqrt(egamma*egamma+q2);
    const double epi = (wsq+mpi0*mpi0-mp*mp)/(2*model.w);
    const double ppi = std::sqrt(epi*epi-mpi0*mpi0);
    const double tmin = -q2+mpi0*mpi0-2*egamma*epi+2*pgamma*ppi;
    model.t = tmin+tprime;
    const double costheta = 1+tprime/(2*pgamma*ppi);
    if (costheta < -1-1e-12 || costheta > 1+1e-12)
        throw std::runtime_error("Unphysical tprime at SIMC reference point");
    const double theta = std::acos(std::clamp(costheta, -1.0, 1.0));
    model.angular = angular_coefficients(q2, wsq, model.t, theta, model.epsilon);
    return model;
}
} // namespace nps_simc_pi0
#endif
