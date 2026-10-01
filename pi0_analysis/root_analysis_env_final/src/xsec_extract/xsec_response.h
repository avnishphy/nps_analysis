#pragma once

// Integrated forward response: rows are reconstructed bins; columns are
// generated bins times U/LT/TT. Weights already include luminosity, generation
// volume, radiation and acceptance. Never probability-normalize the columns
// or divide by reconstructed-bin widths again.
// Ali thesis p.142, Eqs.(5.2)-(5.5): vertex component factors are integrated
// into reconstructed yields. README.md records both requested thesis sources
// and the separate provenance/limitations of our exterior nuisance model.
#include "xsec_physics.h"
#include "xsec_linear_solver.h"

namespace nps_xsec {
// Keep the full outer product of the three correlated Fourier terms. Events
// belong to one row/block. Finite MC variance uses the Poissonized convention;
// fixed-generated-count multinomial covariance is not included.
struct ResponseCell {
    std::array<double, 3> basis{};
    std::array<double, 9> covariance{};
    std::array<double, 3> reco_q2{}, reco_xb{}, reco_tprime{};
    // Cross moments with sum(base_w*epsilon) retain the event-level
    // longitudinal response for a joint T/L fit and propagate its covariance
    // with U/LT/TT. They also support truth-phi reference-epsilon errors.
    double epsilon_weight = 0., epsilon_weight_covariance = 0.;
    std::array<double, 3> basis_epsilon_covariance{};
    long long events = 0;
    void add(const std::array<double, 3> &b, double q2 = 0, double xb = 0,
             double tp = 0, double eps_weight = 0) {
        ++events;
        epsilon_weight += eps_weight;
        epsilon_weight_covariance += eps_weight * eps_weight;
        for (int i = 0; i < 3; ++i) {
            basis[i] += b[i];
            basis_epsilon_covariance[i] += b[i] * eps_weight;
            reco_q2[i] += b[i] * q2;
            reco_xb[i] += b[i] * xb;
            reco_tprime[i] += b[i] * tp;
            for (int j = 0; j < 3; ++j)
                covariance[3 * i + j] += b[i] * b[j];
        }
    }
};

// Response-weighted generated coordinates define display/GK reference points.
// They do not convert the fitted binwise constants to point cross sections.
struct TruthMoments {
    double weight = 0, q2 = 0, xb = 0, t = 0, tprime = 0, epsilon = 0;
    // Positivity must cover the event-level epsilon in the response, not just
    // its weighted mean. This envelope includes accepted events in all rows.
    double epsilon_max = 0;
    long long events = 0;
    void add(double w, double q, double x, double tt, double tp, double eps) {
        weight += w;
        q2 += w * q;
        xb += w * x;
        t += w * tt;
        tprime += w * tp;
        epsilon += w * eps;
        epsilon_max = std::max(epsilon_max, eps);
        ++events;
    }
};

// Exact two-body signed forward limit. SIMC ti stores positive -t.
inline double forward_t(double q2, double W, double mp, double mpi) {
    if (!(q2 > 0 && W > mp + mpi))
        throw std::runtime_error("Invalid vertex Q2/W for pi0 t_min");
    const double q0 = (W * W - mp * mp - q2) / (2 * W), epi = (W * W + mpi * mpi - mp * mp) / (2 * W);
    return mpi * mpi - q2 -
           2 * (q0 * epi - std::sqrt(q0 * q0 + q2) * std::sqrt(std::max(0.0, epi * epi - mpi * mpi)));
}

inline const char *guard_name(int face) {
    static const char *names[] = {"tprime_below", "tprime_above", "q2_below",
                                  "q2_above",     "xb_below",     "xb_above"};
    return names[face];
}

// Six disjoint exterior regions receive free U/LT/TT nuisance coefficients.
// Corners follow explicit priority t', Q2, xB; only populated blocks are fit.
// No discarded feed-in, fixed generator-model background, or edge clamping.
// Coarse exterior shapes still require guard-definition/closure variations.
// Small migration probability does not imply a small fitted contribution:
// unconstrained exterior coefficients may compensate other bins. Full rank
// establishes a numerical solution, not precise or physical nuisance values.
inline int truth_block(double q, double x, double tp, const std::vector<double> &qe,
                       const std::vector<std::vector<double>> &xe, const std::vector<double> &te) {
    const int nq = static_cast<int>(qe.size()) - 1, nx = static_cast<int>(xe.front().size()) - 1;
    const int count = (static_cast<int>(te.size()) - 1) * nq * nx;
    if (!(std::isfinite(q) && std::isfinite(x) && std::isfinite(tp)))
        throw std::runtime_error("Nonfinite generated coordinate");
    if (tp < te.front())
        return count;
    if (tp > te.back())
        return count + 1;
    if (q < qe.front())
        return count + 2;
    if (q > qe.back())
        return count + 3;
    const int iq = find_bin(qe, q);
    if (x < xe[iq].front())
        return count + 4;
    if (x > xe[iq].back())
        return count + 5;
    return (find_bin(te, tp) * nq + iq) * nx + find_bin(xe[iq], x);
}

// Conditional response variance for one reconstructed prediction. The stored
// event outer products are positive semidefinite, although summation in double
// precision can leave a roundoff-sized negative quadratic form for a nearly
// null direction. Sum in long double with compensation and compare a negative
// result to the absolute term sum; never hide NaN, overflow, or a material
// covariance failure behind an unconditional max(0, variance).
inline double mc_prediction_variance(const std::vector<ResponseCell> &row, const std::vector<int> &blocks,
                                     const std::vector<double> &p) {
    if (p.size() != 3 * blocks.size())
        throw std::runtime_error("MC prediction variance: parameter/block size mismatch");
    long double variance = 0.L;
    for (size_t b = 0; b < blocks.size(); ++b) {
        if (blocks[b] < 0 || static_cast<size_t>(blocks[b]) >= row.size())
            throw std::runtime_error("MC prediction variance: invalid truth block");
        for (int i = 0; i < 3; ++i)
            if (!std::isfinite(p[3 * b + i]))
                throw std::runtime_error("MC prediction variance: nonfinite coefficient");
        long double cell_variance = 0.L, correction = 0.L, absolute_terms = 0.L;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                const double cov = row[blocks[b]].covariance[3 * i + j];
                if (!std::isfinite(cov))
                    throw std::runtime_error("MC prediction variance: nonfinite event covariance");
                const long double term = static_cast<long double>(p[3 * b + i]) * cov * p[3 * b + j];
                if (!std::isfinite(term))
                    throw std::runtime_error("MC prediction variance: quadratic term overflow");
                absolute_terms += std::abs(term);
                const long double increment = term - correction;
                const long double next = cell_variance + increment;
                correction = (next - cell_variance) - increment;
                cell_variance = next;
            }
        const long double roundoff = 128.L * std::numeric_limits<double>::epsilon() * absolute_terms;
        if (!std::isfinite(cell_variance) || !std::isfinite(absolute_terms) || cell_variance < -roundoff)
            throw std::runtime_error(
                "MC prediction variance: invalid negative/overflowed covariance quadratic form");
        variance += std::max(0.L, cell_variance);
    }
    const double result = static_cast<double>(variance);
    if (!std::isfinite(result))
        throw std::runtime_error("MC prediction variance exceeds finite double range");
    return result;
}

// Minimum response over phi, used only to flag nonphysical solutions; never
// modifies coefficients or conceals a negative response with clipping.
inline double minimum_response(double u, double lt, double tt, double eps) {
    // Put z=cos(phi): the response bracket is a*z^2+b*z+c on [-1,1].
    // Check both endpoints and any interior minimum. The omitted positive
    // 1/(2pi) does not affect the sign; this returns the bracket, not d2sigma.
    const double a = 2 * eps * tt, b = std::sqrt(2 * eps * (1 + eps)) * lt, c = u - eps * tt;
    double answer = std::min(a + b + c, a - b + c);
    if (a > 0 && std::abs(b / (2 * a)) <= 1)
        answer = std::min(answer, c - b * b / (4 * a));
    return answer;
}
} // namespace nps_xsec
