#pragma once

// Continuous nonnegative U_s + LT/TT angular response, with independent U_s
// and common LT/TT. Each constraint maps a setting/block to its three global
// parameter columns. As in solve_positive_response, the setting's event
// epsilon maximum suffices for the entire epsilon envelope and all phi.
#include "xsec_positive_solver.h"

namespace nps_xsec {

inline PositiveSolution solve_joint_positive_response(
    const std::vector<std::vector<double>>& X, const std::vector<double>& y,
    const std::vector<double>& variance, const std::vector<double>& epsilon_max,
    const std::vector<std::array<int,3>>& columns,
    double rank_tolerance = 1e-10) {
    using namespace positive_detail;
    PositiveSolution output;
    output.fit = solve_weighted_response(X, y, variance, rank_tolerance);
    const std::size_t n = output.fit.parameters.size(), blocks = epsilon_max.size();
    if (columns.size() != blocks) fail("joint epsilon/constraint count mismatch");
    std::vector<bool> is_u(n, false), is_interference(n, false);
    for (const auto& index : columns) {
        for (int col : index)
            if (col < 0 || static_cast<std::size_t>(col) >= n) fail("invalid joint constraint column");
        is_u[index[0]] = true;
        is_interference[index[1]] = is_interference[index[2]] = true;
    }
    for (std::size_t j = 0; j < n; ++j)
        if (is_u[j] == is_interference[j]) fail("invalid or missing joint coefficient role");
    for (double e : epsilon_max)
        if (!std::isfinite(e) || e <= 0 || e >= 1) fail("invalid event epsilon maximum");
    Vector p0(output.fit.parameters.begin(), output.fit.parameters.end()), p = p0;
    Matrix factor(n, Vector(n, 0));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            factor[i][j] = output.fit.inverse_information_factor[i*n+j];
    Real anchor_u = 0;
    for (std::size_t j = 0; j < n; ++j)
        anchor_u = std::max(anchor_u, std::abs(p0[j]) +
            std::sqrt(static_cast<Real>(output.fit.covariance[j*n+j])));
    if (!(anchor_u > 0) || !std::isfinite(anchor_u)) fail("invalid joint parameter scale");
    Vector anchor_delta(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        anchor_delta[i] = (is_u[i] ? anchor_u : 0) - p0[i];
    const Vector anchor = factor_coordinates(factor, anchor_delta);
    std::vector<Cut> cuts;
    const auto add_cut = [&](std::size_t b, Real z) {
        const Real e = epsilon_max[b];
        for (const auto& old : cuts)
            if (old.block == b && std::abs(old.cosine-z) < 1e-12L)
                fail("joint continuous constraint generation stalled");
        const std::array<Real,3> angular{{1,
            std::sqrt(2*e*(1+e))*z, e*(2*z*z-1)}};
        Cut cut;
        cut.block=b; cut.epsilon=e; cut.cosine=z; cut.normal.assign(n,0);
        for (int a=0; a<3; ++a) {
            const auto col = static_cast<std::size_t>(columns[b][a]);
            cut.offset += angular[a]*p0[col];
            for (std::size_t j=0; j<n; ++j)
                cut.normal[j] += angular[a]*factor[col][j];
        }
        const Real length=norm(cut.normal);
        if (!(length>0) || !std::isfinite(length)) fail("invalid joint constraint norm");
        cut.offset/=length;
        for (Real& value:cut.normal) value/=length;
        cuts.push_back(std::move(cut));
    };
    for (std::size_t b=0; b<blocks; ++b) {
        add_cut(b,-1);
        add_cut(b,0);
        add_cut(b,1);
    }
    bool solved=false;
    for (int exchange=0; exchange<256; ++exchange) {
        bool feasible=true;
        output.minima.clear(); output.min_cos.clear();
        output.boundary_tolerances.clear(); output.feasibility_tolerances.clear();
        output.boundary_active=false;
        for (std::size_t b=0; b<blocks; ++b) {
            const auto& index = columns[b];
            const auto m=minimum(p[index[0]],p[index[1]],p[index[2]],epsilon_max[b]);
            const Real e=epsilon_max[b];
            const Real scale=std::max({std::abs(p[index[0]]),
                std::sqrt(2*e*(1+e))*std::abs(p[index[1]]),e*std::abs(p[index[2]]),anchor_u*1e-6L});
            output.minima.push_back(static_cast<double>(m.first));
            output.min_cos.push_back(static_cast<double>(m.second));
            output.boundary_tolerances.push_back(static_cast<double>(2e-8L*scale));
            output.feasibility_tolerances.push_back(static_cast<double>(2e-10L*scale));
            if (m.first < -2e-10L*scale) {
                feasible=false;
                if (exchange>0) add_cut(b,m.second);
            }
            if (m.first <= 2e-8L*scale) output.boundary_active=true;
        }
        if (feasible) {solved=true;break;}
        const Vector z=project(cuts,anchor,output.iterations);
        p=p0;
        for (std::size_t i=0; i<n; ++i)
            for (std::size_t j=0; j<n; ++j) p[i]+=factor[i][j]*z[j];
    }
    if (!solved) fail("joint continuous epsilon/angle constraints did not converge");
    output.fit.parameters.assign(p.begin(),p.end());
    Real chi2=0;
    for (std::size_t row=0; row<X.size(); ++row) {
        Real prediction=0;
        for (std::size_t j=0; j<n; ++j) prediction+=X[row][j]*p[j];
        const Real residual=y[row]-prediction;
        chi2+=residual*residual/variance[row];
    }
    output.fit.chi2=static_cast<double>(chi2);
    if (!std::isfinite(output.fit.chi2)) fail("joint constrained chi-square is nonfinite");
    return output;
}
} // namespace nps_xsec
