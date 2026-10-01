#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// physics_pion.f sig_param_2021/exclfit, doing_pizero=.true., W>=2 GeV.
// Fortran uses the arithmetic mean of the pi+ and pi- fits with fpifact=0.
// q2,t are GeV^2, W is GeV, phi is radians; result is ub/MeV^2/radian.
namespace nps_pi0_reweight {
struct Parameter {
    std::string name;
    double value, lower, upper;
    bool free_by_default;
};
inline const char* identifier() { return "SIMC_sig_param_2021_pi0_W_ge_2"; }
inline const char* provenance() { return "simc_gfortran_updated/physics_pion.f:738-855"; }
inline std::vector<Parameter> parameters() {
    const std::array<double, 17> plus{{
        1.60077, -0.01523, 37.08142, -4.11060, 23.26192, 0.00983,
        0.87073, -5.77115, -271.08678, 0.13766, -0.00855, 0.27885,
        -1.13212, -1.50415, -6.34766, 0.55769, -0.01709}};
    const std::array<double, 17> minus{{
        1.75169, 0.11144, 47.35877, -4.69434, 1.60552, 0.00800,
        0.44194, -2.29188, -41.67194, 0.69475, 0.02527, -0.50178,
        -1.22825, -1.16878, 5.75825, -1.00355, 0.05055}};
    std::vector<Parameter> result;
    for (int charge=0; charge<2; ++charge)
        for (int j=0; j<17; ++j) {
            const double v = charge == 0 ? plus[j] : minus[j];
            const double span = std::max(1.0, 5.0*std::abs(v));
            // Default fit is a modest identifiable T, LT, TT subset of
            // the actual plus-fit coefficients, shared by all selected bins.
            const bool active = charge == 0 && (j == 4 || j == 6 || j == 8);
            result.push_back({std::string(charge == 0 ? "plus.p" : "minus.p")
                              + std::to_string(j+1), v,
                              j == 4 ? 0.0 : v-span, v+span, active});
        }
    return result;
}
inline std::vector<double> defaults() {
    std::vector<double> result;
    for (const auto& p : parameters()) result.push_back(p.value);
    return result;
}
inline double epsilon(double q2, double w, double ebeam, double mp) {
    const double nu = (w*w-mp*mp+q2)/(2*mp);
    const double eprime = ebeam-nu;
    if (!(std::isfinite(q2) && q2 > 0 && std::isfinite(w) && w >= 2 &&
          std::isfinite(ebeam) && ebeam > 0 && eprime > 0 &&
          q2 <= 4*ebeam*eprime))
        throw std::runtime_error("Unphysical electron kinematics for epsilon");
    const double y=nu/ebeam, z=q2/(4*ebeam*ebeam);
    return (1-y-z)/(1-y+0.5*y*y+z);
}
inline double evaluate(double q2, double w, double t, double phi, double eps,
                       double mp, double mpi0, const std::vector<double>& coeff) {
    if (coeff.size()!=34 || !(std::isfinite(q2) && q2>0 &&
          std::isfinite(w) && w>=2 && std::isfinite(t) &&
          std::isfinite(phi) && std::isfinite(eps) && eps>=0 && eps<=1))
        throw std::runtime_error("Unsupported SIMC pi0 event kinematics");
    const double wsq=w*w, eg=(wsq-mp*mp-q2)/(2*w);
    const double pg=std::sqrt(eg*eg+q2);
    const double epi=(wsq+mpi0*mpi0-mp*mp)/(2*w);
    const double pp2=epi*epi-mpi0*mpi0;
    if (!(pp2>0)) throw std::runtime_error("SIMC pi0 below threshold");
    const double ct=(t+q2-mpi0*mpi0+2*eg*epi)/(2*pg*std::sqrt(pp2));
    if (!std::isfinite(ct) || ct < -1-1e-5 || ct > 1+1e-5)
        throw std::runtime_error("SIMC pi0 unphysical CM angle");
    const double st=std::sqrt(std::max(0.0,1-std::clamp(ct,-1.0,1.0)*std::clamp(ct,-1.0,1.0)));
    const double at=std::abs(t);
    // Preserve exclfit's 0.938 mass, 3.1415928, and 1e6 unit conversion.
    const double norm=8.539/std::pow(wsq-0.938*0.938,2)/2/3.1415928/1.e6;
    double sigma=0;
    for (int k=0;k<2;++k) {
        const double* p=coeff.data()+17*k;
        // sigL=0 because fpifact=0. Other coefficients are untouched.
        const double sigT=p[4]/q2*std::exp(p[5]*q2*q2)/
            (std::pow(wsq,p[11])+std::pow(w,p[15]))*std::exp(p[13]*at);
        const double sigLT=p[6]/(1+p[9]*q2)*std::exp(p[7]*at)*st/
            std::pow(wsq,p[12]);
        const double sigTT=p[8]/(1+q2)*std::exp(-7*at)*st*st;
        sigma += 0.5*norm*(sigT+std::sqrt(2*eps*(1+eps))*sigLT*std::cos(phi)
                           +eps*sigTT*std::cos(2*phi));
    }
    if (!std::isfinite(sigma)) throw std::runtime_error("Nonfinite SIMC pi0 prediction");
    return sigma;
}
// A model registry keeps the extractor independent of the physics formula.
// Add another Definition here when a separately validated alternative exists.
struct Definition {
    const char* id;
    const char* source;
    std::vector<Parameter> (*parameter_spec)();
    std::vector<double> (*default_parameters)();
    double (*cross_section)(double,double,double,double,double,double,double,
                            const std::vector<double>&);
};
inline const Definition& choose(const std::string& requested) {
    static const Definition simc{identifier(),provenance(),parameters,defaults,evaluate};
    if (requested == simc.id) return simc;
    throw std::runtime_error("Unknown pi0 model identifier: "+requested);
}
} // namespace nps_pi0_reweight
