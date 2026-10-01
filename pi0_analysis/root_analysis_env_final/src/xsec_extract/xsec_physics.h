#pragma once

// Kinematic conventions and Fourier basis. Cross sections are microbarn/MeV^2; angles are radians.
#include "xsec_types.h"

// ROOT can otherwise report an incompatible/missing branch and leave the
// destination variable unchanged, silently creating repeated event values.
template<class T> inline void bind_branch(TTree* tree,const char* name,T* address) {
    if(tree->SetBranchAddress(name,address)<0)
        throw std::runtime_error(std::string("Cannot bind branch ")+name+" in "+tree->GetName());
}

static double wrap_phi(double x) {
    double y = std::fmod(x, 2.0 * TMath::Pi());
    if (y < 0) y += 2.0 * TMath::Pi();
    // map exact 2pi to 0
    if (y >= 2.0 * TMath::Pi()) y = 0.0;
    return y;
}

static double clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

static double kallen_lambda_sqrt(double a, double b, double c) {
    double arg = (a - (b + c)) * (a - (b + c)) - 4.0 * b * c;
    return (arg > 0.0) ? std::sqrt(arg) : 0.0;
}

static double q2_xb_to_w2(double q2, double xb, double mp) {
    // W^2 = M^2 + Q^2(1/xB - 1)
    return mp * mp + q2 * (1.0 / xb - 1.0);
}

static double epsilon_virtual(double ebeam, double q2, double xb, double mp) {
    // Using y = nu/E = Q2/(2 M xB E)
    // epsilon = [1 - y - Q2/(4E^2)] / [1 - y + y^2/2 + Q2/(4E^2)]
    // This is the standard electron-scattering form for negligible electron mass.
    if (ebeam <= 0.0 || q2 <= 0.0 || xb <= 0.0) return 0.0;
    double y = q2 / (2.0 * mp * xb * ebeam);
    if (y <= 0.0) return 0.0;
    double e2 = ebeam * ebeam;
    double term = q2 / (4.0 * e2);
    double num = 1.0 - y - term;
    double den = 1.0 - y + 0.5 * y * y + term;
    if (den <= 0.0) return 0.0;
    return clamp(num / den, 0.0, 1.0);
}

static double virtual_photon_flux(double ebeam, double q2, double xb, double mp, double epsilon) {
    // Gamma_{gamma*} = alpha/(8pi) * Q2/(M^2 E^2) * (1-xB)/xB^3 * 1/(1-epsilon).
    // This is a flux factor, not an additional SIMC generation Jacobian.
    constexpr double alpha = 1.0 / 137.035999084;
    if (ebeam <= 0.0 || q2 <= 0.0 || xb <= 0.0 || mp <= 0.0) return 0.0;
    const double eps = clamp(epsilon, 0.0, 1.0 - 1e-12);
    return (alpha / (8.0 * TMath::Pi())) *
           (q2 / (mp * mp * ebeam * ebeam)) *
           ((1.0 - xb) / (xb * xb * xb)) *
           (1.0 / (1.0 - eps));
}

static int find_bin(const std::vector<double>& edges, double x, bool periodic_phi = false) {
    if (edges.size() < 2) return -1;
    if (!std::isfinite(x)) return -1;
    if (periodic_phi) {
        x = wrap_phi(x);
    }

    if (x < edges.front() || x > edges.back()) return -1;
    if (x == edges.back()) return static_cast<int>(edges.size()) - 2;

    auto it = std::upper_bound(edges.begin(), edges.end(), x);
    int idx = static_cast<int>(it - edges.begin()) - 1;
    if (idx < 0 || idx >= static_cast<int>(edges.size()) - 1) return -1;
    return idx;
}

static std::vector<double> phi_basis_means(double phi1, double phi2) {
    const double d = phi2 - phi1;
    if (d <= 0.0) return {1.0, 0.0, 0.0, 0.0};
    double c1 = (std::sin(phi2) - std::sin(phi1)) / d;
    double c2 = (std::sin(2.0 * phi2) - std::sin(2.0 * phi1)) / (2.0 * d);
    double s1 = (-std::cos(phi2) + std::cos(phi1)) / d;
    return {1.0, c1, c2, s1};
}

static std::array<double, 4> sigma_model_basis(double phi, double epsilon, double gamma_flux, int helicity_sign) {
    const double inv2pi = 1.0 / (2.0 * TMath::Pi());
    const double eps = clamp(epsilon, 0.0, 1.0);
    (void)gamma_flux;
    // const double gamma = (std::isfinite(gamma_flux) && gamma_flux > 0.0) ? gamma_flux : 0.0;
    // Gamma is intentionally not multiplied here. With w_base = full_weight/sigcm,
    // the retained factor siglab/sigcm already contains SIMC's davejac*gtpr*fac.
    // Reapplying this macro's Gamma would double-count the flux-like part.
    const double k_tl = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 + eps)));
    const double k_tlp = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 - eps)));
    return {
        inv2pi,
        inv2pi * k_tl * std::cos(phi),
        inv2pi * eps * std::cos(2.0 * phi),
        inv2pi * static_cast<double>(helicity_sign) * k_tlp * std::sin(phi)
    };
}

static double sigma_model_value(double phi,
                                double epsilon,
                                double gamma_flux,
                                double sigma_u,
                                double sigma_tl,
                                double sigma_tt,
                                double sigma_tlp,
                                int helicity_sign) {
    const auto b = sigma_model_basis(phi, epsilon, gamma_flux, helicity_sign);
    return b[0] * sigma_u + b[1] * sigma_tl + b[2] * sigma_tt + b[3] * sigma_tlp;
}
