#ifndef NPS_XSEC_POSITIVE_SOLVER_H
#define NPS_XSEC_POSITIVE_SOLVER_H

// Optional physical constraint on the FULL angular cross section. LT and TT
// remain signed. This header has no event selection, binning, or MC iteration;
// the caller repeats this solve with updated finite-MC variances as for GLS.
#include "xsec_linear_solver.h"

#include <algorithm>
#include <array>
#include <utility>

namespace nps_xsec {

struct PositiveSolution {
  LinearSolution fit;
  int iterations = 0;  // Total finite-QP working-set iterations.
  std::vector<double> minima, min_cos; // Bracket minima, with no 1/(2*pi).
  // Per-block native bracket tolerances. Export these with minima so plots,
  // CSV flags, and the global boundary flag use precisely the same test.
  std::vector<double> boundary_tolerances, feasibility_tolerances;
  bool boundary_active = false;
};

namespace positive_detail {
using Real = long double;
using Vector = std::vector<Real>;
using Matrix = std::vector<Vector>;

inline Real dot(const Vector& a, const Vector& b) {
  Real sum = 0;
  for (std::size_t i = 0; i < a.size(); ++i) sum += a[i] * b[i];
  return sum;
}
inline Real norm(const Vector& a) { return std::sqrt(dot(a, a)); }
inline void fail(const std::string& why) {
  throw std::runtime_error("Positive response fit: " + why);
}

// Solve the dense SVD-factor system F z = delta without forming F F^T.
// Normalize each equation first: physical-unit changes may span hundreds of
// orders of magnitude, while only the scaled response decides identifiability.
// Long-double partial-pivot elimination does not square the condition number.
inline Vector factor_coordinates(Matrix factor, Vector delta) {
  const std::size_t n = delta.size();
  if (factor.size() != n) fail("factor/coordinate dimension mismatch");
  for (std::size_t i = 0; i < n; ++i) {
    if (factor[i].size() != n) fail("non-square inverse-information factor");
    Real scale = 0;
    for (Real v : factor[i]) {
      if (!std::isfinite(v)) fail("nonfinite inverse-information factor");
      scale = std::max(scale, std::abs(v));
    }
    if (!(scale > 0) || !std::isfinite(delta[i])) fail("invalid factor coordinate equation");
    for (Real& v : factor[i]) v /= scale;
    delta[i] /= scale;
  }
  for (std::size_t j = 0; j < n; ++j) {
    std::size_t pivot = j;
    for (std::size_t i = j + 1; i < n; ++i)
      if (std::abs(factor[i][j]) > std::abs(factor[pivot][j])) pivot = i;
    if (!(std::abs(factor[pivot][j]) > 0)) fail("singular inverse-information factor");
    std::swap(factor[j], factor[pivot]);
    std::swap(delta[j], delta[pivot]);
    for (std::size_t i = j + 1; i < n; ++i) {
      const Real multiplier = factor[i][j] / factor[j][j];
      for (std::size_t k = j + 1; k < n; ++k) factor[i][k] -= multiplier * factor[j][k];
      delta[i] -= multiplier * delta[j];
    }
  }
  Vector z(n, 0);
  for (std::size_t ii = n; ii > 0; --ii) {
    const std::size_t i = ii - 1;
    Real rhs = delta[i];
    for (std::size_t j = i + 1; j < n; ++j) rhs -= factor[i][j] * z[j];
    z[i] = rhs / factor[i][i];
    if (!std::isfinite(z[i])) fail("nonfinite factor coordinate");
  }
  return z;
}

// Exact all-angle minimization, z=cos(phi) in [-1,1]. No angular mesh is used.
inline std::pair<Real, Real> minimum(Real u, Real lt, Real tt, Real e) {
  const Real a = 2 * e * tt, b = std::sqrt(2 * e * (1 + e)) * lt;
  const auto value = [&](Real z) { return u + b * z + e * tt * (2 * z * z - 1); };
  Real z = b > 0 ? -1 : 1;
  Real result = value(z);
  if (a > 0) {
    const Real vertex = -b / (2 * a);
    if (vertex > -1 && vertex < 1 && value(vertex) < result) {
      z = vertex;
      result = value(z);
    }
  }
  return {result, z};
}

struct Cut {
  Vector normal; // Unit normal in covariance-whitened coordinates.
  Real offset = 0; // offset + normal.dot(z) >= 0.
  std::size_t block = 0;
  Real cosine = 0;
  Real epsilon = 0; // Used by the four-component joint T/L constraint.
};

// QR of the active normals, reorthogonalized twice. Working directly with QR
// avoids squaring the condition number of nearly parallel angular cuts. No
// pseudoinverse/rank truncation is used to disguise a degenerate active set.
inline bool active_qr(const std::vector<Cut>& cuts, const std::vector<std::size_t>& active,
                      std::size_t n, Matrix& q, Matrix& r) {
  const std::size_t k = active.size();
  q.clear();
  r.assign(k, Vector(k, 0));
  if (k > n) return false;
  for (std::size_t j = 0; j < k; ++j) {
    Vector v = cuts[active[j]].normal;
    for (int pass = 0; pass < 2; ++pass)
      for (std::size_t i = 0; i < j; ++i) {
        const Real projection = dot(q[i], v);
        r[i][j] += projection;
        for (std::size_t l = 0; l < n; ++l) v[l] -= projection * q[i][l];
      }
    r[j][j] = norm(v);
    if (!(r[j][j] > 1e-13L)) return false;
    for (Real& x : v) x /= r[j][j];
    q.push_back(std::move(v));
  }
  return true;
}

// Strictly convex finite QP: min ||z||^2/2 subject to all supplied cuts.
// A feasible interior anchor initializes a primal active-set method. Every
// accepted step stays feasible; negative KKT multipliers release constraints.
inline Vector project(const std::vector<Cut>& cuts, const Vector& anchor, int& iterations) {
  const std::size_t n = anchor.size();
  Vector z = anchor;
  std::vector<std::size_t> active;
  Matrix q, r;
  for (const auto& cut : cuts)
    if (!(cut.offset + dot(cut.normal, z) > 0)) fail("interior anchor is not feasible");
  for (int step = 0; step < 20000; ++step) {
    ++iterations;
    if (!active_qr(cuts, active, n, q, r)) fail("dependent QP working constraints");
    const std::size_t k = active.size();
    Vector coordinates(k, 0), candidate(n, 0), lambda(k, 0);
    // R^T coordinates = -offset, then candidate = Q coordinates.
    for (std::size_t j = 0; j < k; ++j) {
      Real rhs = -cuts[active[j]].offset;
      for (std::size_t i = 0; i < j; ++i) rhs -= r[i][j] * coordinates[i];
      coordinates[j] = rhs / r[j][j];
      for (std::size_t l = 0; l < n; ++l) candidate[l] += coordinates[j] * q[j][l];
    }
    Vector direction(n, 0);
    for (std::size_t j = 0; j < n; ++j) direction[j] = candidate[j] - z[j];
    const Real tolerance = 2e-12L * (1 + norm(z));
    if (norm(direction) <= tolerance) {
      // At the active-face minimum, z = A_active^T lambda. Nonnegative
      // lambda plus primal feasibility certifies the convex QP optimum.
      for (std::size_t jj = k; jj > 0; --jj) {
        const std::size_t j = jj - 1;
        Real rhs = coordinates[j];
        for (std::size_t i = j + 1; i < k; ++i) rhs -= r[j][i] * lambda[i];
        lambda[j] = rhs / r[j][j];
      }
      const auto worst = std::min_element(lambda.begin(), lambda.end());
      if (worst == lambda.end() || *worst >= -tolerance) {
        for (const auto& cut : cuts)
          if (cut.offset + dot(cut.normal, candidate) < -8 * tolerance)
            fail("QP failed final primal feasibility check");
        return candidate;
      }
      active.erase(active.begin() + std::distance(lambda.begin(), worst));
      continue;
    }
    Real fraction = 1;
    std::size_t blocker = cuts.size();
    for (std::size_t i = 0; i < cuts.size(); ++i) {
      if (std::find(active.begin(), active.end(), i) != active.end()) continue;
      const Real slope = dot(cuts[i].normal, direction);
      if (slope >= -2e-14L * (1 + norm(direction))) continue;
      const Real gap = cuts[i].offset + dot(cuts[i].normal, z);
      if (gap < -8 * tolerance) fail("QP lost primal feasibility");
      const Real limit = std::max(0.L, gap) / -slope;
      if (limit < fraction) { fraction = limit; blocker = i; }
    }
    for (std::size_t j = 0; j < n; ++j) z[j] += fraction * direction[j];
    if (blocker != cuts.size()) active.push_back(blocker);
  }
  fail("finite QP did not converge within 20000 working-set iterations");
  return {};
}
} // namespace positive_detail

// For fixed epsilon e, the bracket minimum is U+e*TT-s*|LT| for TT<=0
// or s*|LT|>=4e*TT, otherwise U-e*TT-(1+e)*LT^2/(4TT), with
// s=sqrt(2e(1+e)). The lower bound on U increases monotonically with e.
// Therefore positivity at each block's largest event epsilon is EXACTLY
// sufficient for all phi and every 0<=epsilon<=epsilon_max in that block.
// This applies equally to all fitted nuisance/guard triplets. The caller must
// supply EVENT maxima, not weighted mean epsilon or reconstructed epsilon.
//
// IMPORTANT: fit.covariance remains (X' W X)^-1 from the unconstrained SVD.
// It is a curvature diagnostic, NOT the sampling covariance of a constrained
// estimate on a boundary. Boundary-aware intervals need repeated full fits or
// suitable profile calibration; callers must not label these errors as such.
inline PositiveSolution solve_positive_response(
    const std::vector<std::vector<double>>& X, const std::vector<double>& y,
    const std::vector<double>& variance, const std::vector<double>& epsilon_max,
    double rank_tolerance = 1e-10) {
  using namespace positive_detail;
  PositiveSolution output;
  // Always retain the original full-rank identifiability test. Positivity is
  // not permission to recover an unsupported column or singular truth bin.
  output.fit = solve_weighted_response(X, y, variance, rank_tolerance);
  const std::size_t n = output.fit.parameters.size(), blocks = epsilon_max.size();
  if (n != 3 * blocks) fail("epsilon/block count does not match coefficient triplets");
  for (double e : epsilon_max)
    if (!std::isfinite(e) || e < 0 || e > 1) fail("event epsilon maximum must lie in [0,1]");
  Vector p0(output.fit.parameters.begin(), output.fit.parameters.end()), p = p0;
  Matrix factor(n, Vector(n, 0));
  // p=p0+F*z makes the objective ||z||^2/2. F is dense, not triangular;
  // use the SVD factor itself, never Cholesky of the rounded covariance.
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      factor[i][j] = output.fit.inverse_information_factor[i * n + j];
  std::vector<Cut> cuts;
  const auto add_cut = [&](std::size_t block, Real cosine) {
    for (const auto& old : cuts)
      if (old.block == block && std::abs(old.cosine - cosine) < 1e-12L)
        fail("angular constraint generation stalled at a repeated violating cut");
    const Real e = epsilon_max[block];
    const std::array<Real, 3> angular{{1, std::sqrt(2 * e * (1 + e)) * cosine,
                                      e * (2 * cosine * cosine - 1)}};
    Cut cut;
    cut.block = block; cut.cosine = cosine; cut.normal.assign(n, 0);
    for (int a = 0; a < 3; ++a) {
      cut.offset += angular[a] * p0[3 * block + a];
      for (std::size_t j = 0; j < n; ++j)
        cut.normal[j] += angular[a] * factor[3 * block + a][j];
    }
    const Real length = norm(cut.normal);
    if (!(length > 0) || !std::isfinite(length)) fail("invalid angular constraint norm");
    cut.offset /= length;
    for (Real& v : cut.normal) v /= length;
    cuts.push_back(std::move(cut));
  };
  Real anchor_u = 0;
  for (std::size_t j = 0; j < n; ++j)
    anchor_u = std::max(anchor_u, std::abs(p0[j]) + std::sqrt(static_cast<Real>(output.fit.covariance[j*n+j])));
  if (!(anchor_u > 0) || !std::isfinite(anchor_u)) fail("invalid physical coefficient scale");
  Vector anchor_delta(n, 0);
  for (std::size_t i = 0; i < n; ++i)
    anchor_delta[i] = (i % 3 == 0 ? anchor_u : 0) - p0[i];
  const Vector anchor = factor_coordinates(factor, anchor_delta);
  for (std::size_t b = 0; b < blocks; ++b) {
    add_cut(b, -1);
    if (epsilon_max[b] > 0) { add_cut(b, 0); add_cut(b, 1); }
  }
  bool solved = false;
  for (int exchange = 0; exchange < 128; ++exchange) {
    bool feasible = true;
    output.minima.clear(); output.min_cos.clear();
    output.boundary_tolerances.clear(); output.feasibility_tolerances.clear();
    output.boundary_active = false;
    for (std::size_t b = 0; b < blocks; ++b) {
      const Real e = epsilon_max[b], s = std::sqrt(2 * e * (1 + e));
      const auto m = minimum(p[3*b], p[3*b+1], p[3*b+2], e);
      output.minima.push_back(static_cast<double>(m.first));
      output.min_cos.push_back(static_cast<double>(m.second));
      const Real scale = std::max({std::abs(p[3*b]), s*std::abs(p[3*b+1]),
                                  e*std::abs(p[3*b+2]), anchor_u * 1e-6L});
      output.boundary_tolerances.push_back(static_cast<double>(2e-8L * scale));
      output.feasibility_tolerances.push_back(static_cast<double>(2e-10L * scale));
      if (m.first < -2e-10L * scale) {
        feasible = false;
        // Initial endpoint cuts already exist; they are solved first. Later
        // exchanges add only the exact newly discovered angular minimum.
        if (exchange > 0) add_cut(b, m.second);
      }
      if (m.first <= 2e-8L * scale) output.boundary_active = true;
    }
    if (feasible) { solved = true; break; }
    const Vector z = project(cuts, anchor, output.iterations);
    p = p0;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j) p[i] += factor[i][j] * z[j];
  }
  if (!solved) fail("continuous-angle constraint generation did not converge");
  output.fit.parameters.assign(p.begin(), p.end());
  Real chi2 = 0;
  for (std::size_t row = 0; row < X.size(); ++row) {
    Real prediction = 0;
    for (std::size_t j = 0; j < n; ++j) prediction += X[row][j] * p[j];
    const Real residual = y[row] - prediction;
    chi2 += residual * residual / variance[row];
  }
  output.fit.chi2 = static_cast<double>(chi2);
  if (!std::isfinite(output.fit.chi2)) fail("constrained chi-square is not finite");
  // ndf/rank/condition continue to describe the original response geometry;
  // active inequality boundaries do not have a simple fixed chi2 DOF rule.
  return output;
}
} // namespace nps_xsec
#endif // NPS_XSEC_POSITIVE_SOLVER_H
