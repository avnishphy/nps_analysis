// Standalone deterministic checks for the migration-response numerical solver.
// Compile with ROOT's flags and libraries; no production events are required.
#include "../src/xsec_extract/xsec_linear_solver.h"
#include "../src/xsec_extract/xsec_joint_positive_solver.h"
#include "../src/xsec_extract/xsec_response.h"
#include <algorithm>
#include <iostream>

namespace {
void require(bool passed, const char* explanation) {
  if (!passed) throw std::runtime_error(explanation);
}
bool close(double observed, double expected, double tolerance = 1e-10) {
  return std::abs(observed - expected) <= tolerance * std::max(1., std::abs(expected));
}
std::vector<double> predict(const std::vector<std::vector<double>>& X,
                            const std::vector<double>& parameters) {
  std::vector<double> y(X.size(), 0.);
  for (std::size_t r = 0; r < X.size(); ++r)
    for (std::size_t j = 0; j < parameters.size(); ++j) y[r] += X[r][j] * parameters[j];
  return y;
}
template <typename Action> void rejects(Action action, const char* explanation) {
  bool rejected = false;
  try { action(); } catch (const std::runtime_error&) { rejected = true; }
  require(rejected, explanation);
}
void check_factor(const nps_xsec::LinearSolution& fit,
                  const std::vector<std::vector<double>>& X,
                  const std::vector<double>& variance, double tolerance = 1e-10) {
  const std::size_t n = fit.parameters.size();
  require(fit.inverse_information_factor.size() == n * n, "SVD factor dimensions failed");
  std::vector<std::vector<long double>> orthogonal(X.size(), std::vector<long double>(n, 0));
  for (std::size_t r = 0; r < X.size(); ++r)
    for (std::size_t k = 0; k < n; ++k) {
      for (std::size_t j = 0; j < n; ++j)
        orthogonal[r][k] += static_cast<long double>(X[r][j]) * fit.inverse_information_factor[j*n+k];
      orthogonal[r][k] /= std::sqrt(static_cast<long double>(variance[r]));
    }
  for (std::size_t j = 0; j < n; ++j)
    for (std::size_t k = 0; k < n; ++k) {
      long double gram = 0, covariance = 0;
      for (std::size_t r = 0; r < X.size(); ++r) gram += orthogonal[r][j] * orthogonal[r][k];
      for (std::size_t l = 0; l < n; ++l)
        covariance += static_cast<long double>(fit.inverse_information_factor[j*n+l]) *
                      fit.inverse_information_factor[k*n+l];
      require(close(static_cast<double>(gram), j == k ? 1. : 0., tolerance),
              "SVD factor does not whiten the response curvature");
      require(close(static_cast<double>(covariance), fit.covariance[j*n+k], tolerance),
              "SVD factor/covariance mismatch");
    }
}
}

int main() {
  // Two independent copies of identity provide six parameters and positive
  // degrees of freedom; their known unit-variance covariance is I/2.
  std::vector<std::vector<double>> identity(12, std::vector<double>(6, 0.));
  for (std::size_t r = 0; r < 12; ++r) identity[r][r % 6] = 1.;
  const std::vector<double> truth{3., -0.4, 0.2, 7., 0.8, -0.3};
  const std::vector<double> unit_variance(12, 1.);
  const auto exact = nps_xsec::solve_weighted_response(identity, predict(identity, truth), unit_variance);
  for (std::size_t j = 0; j < 6; ++j) {
    require(close(exact.parameters[j], truth[j]), "identity recovery failed");
    for (std::size_t k = 0; k < 6; ++k)
      require(close(exact.covariance[j * 6 + k], j == k ? 0.5 : 0.), "identity covariance failed");
  }
  require(exact.rank == 6 && exact.ndf == 6 && exact.chi2 < 1e-20, "identity diagnostics failed");
  check_factor(exact, identity, unit_variance);

  // Each truth bin has U, LT, TT coefficients. Reconstructed angular yields
  // mix both truth bins with a deliberately asymmetric 2x2 migration matrix.
  // Negative interference coefficients remain admissible in the recovered fit.
  std::vector<std::vector<double>> migration(12, std::vector<double>(6, 0.));
  const double response[2][2] = {{0.75, 0.35}, {0.25, 0.65}};
  const double pi = std::acos(-1.);
  std::vector<double> variance(12);
  for (std::size_t b = 0; b < 2; ++b)
    for (std::size_t p = 0; p < 6; ++p) {
      const std::size_t r = b * 6 + p;
      const double phi = (p + 0.5) * 2. * pi / 6.;
      const double angular[3] = {1., std::cos(phi), std::cos(2. * phi)};
      variance[r] = 1. + 0.1 * r;
      for (std::size_t t = 0; t < 2; ++t)
        for (std::size_t k = 0; k < 3; ++k) migration[r][3 * t + k] = response[b][t] * angular[k];
    }
  const auto y = predict(migration, truth);
  const auto fit = nps_xsec::solve_weighted_response(migration, y, variance);
  for (std::size_t j = 0; j < 6; ++j)
    require(close(fit.parameters[j], truth[j]), "migrated recovery failed");
  require(std::abs(fit.covariance[3]) > 0.01, "cross-bin covariance missing");
  check_factor(fit, migration, variance);
  // Independently verify C (X^T W X) = I; this checks the complete covariance,
  // including cross-bin and cross-structure-function correlations.
  for (std::size_t j = 0; j < 6; ++j)
    for (std::size_t k = 0; k < 6; ++k) {
      double product = 0.;
      for (std::size_t l = 0; l < 6; ++l)
        for (std::size_t r = 0; r < 12; ++r)
          product += fit.covariance[j * 6 + l] * migration[r][l] * migration[r][k] / variance[r];
      require(close(product, j == k ? 1. : 0.), "migration covariance inverse check failed");
    }

  // Unit changes span 240 orders of magnitude. Predictions, rank, condition,
  // and covariance transformed back to the original units must remain stable.
  auto rescaled = migration;
  const double units[6] = {1e120, 1e-120, 1e70, 1e-70, 1e20, 1e-20};
  for (auto& row : rescaled)
    for (std::size_t j = 0; j < 6; ++j) row[j] *= units[j];
  const auto scaled_fit = nps_xsec::solve_weighted_response(rescaled, y, variance);
  for (std::size_t j = 0; j < 6; ++j) {
    require(close(scaled_fit.parameters[j] * units[j], fit.parameters[j]), "unit-scaled parameter failed");
    for (std::size_t k = 0; k < 6; ++k)
      require(close(scaled_fit.covariance[j * 6 + k] * units[j] * units[k], fit.covariance[j * 6 + k]),
              "unit-scaled covariance failed");
  }
  require(close(scaled_fit.condition, fit.condition), "unit-dependent rank conditioning");
  check_factor(scaled_fit, rescaled, variance);

  // These full-rank responses have condition 2e8--6.7e9. Their rounded
  // covariances lose the small eigenvalue (~0.125) beside entries 1/(4*d^2).
  // Cholesky of those covariances failed for d=1e-8, 5e-10 and 3e-10 with
  // ROOT 6.30.04. No response singular direction may be dropped to fix that.
  for (const double delta : {1e-8, 1e-9, 5e-10, 3e-10}) {
    const std::vector<std::vector<double>> ill_conditioned{
        {1., 1.+delta, 0.}, {1., 1.-delta, 0.}, {0., 0., 1.},
        {1., 1.+delta, 0.}, {1., 1.-delta, 0.}, {0., 0., 1.}};
    const std::vector<double> ill_y = predict(ill_conditioned, {-1., 2., .2});
    const std::vector<double> ill_variance(6, 1.);
    const auto ill_linear = nps_xsec::solve_weighted_response(ill_conditioned, ill_y, ill_variance);
    require(ill_linear.rank == 3 && ill_linear.condition > 1e8, "ill-conditioned regression fixture failed");
    check_factor(ill_linear, ill_conditioned, ill_variance, 2e-6);
    const auto zero_epsilon = nps_xsec::solve_positive_response(ill_conditioned, ill_y, ill_variance, {0.});
    require(zero_epsilon.boundary_active && close(zero_epsilon.fit.parameters[0], 0., 1e-7) &&
            close(zero_epsilon.fit.parameters[1], 1., 1e-6) &&
            close(zero_epsilon.fit.parameters[2], .2, 1e-6),
            "ill-conditioned U>=0 boundary projection failed");
    // At epsilon=1/2 the analytic optimum is the phi=pi boundary. The
    // well-measured sum U+LT=1 and TT=.2 remain unchanged to O(delta^2).
    const double expected_lt = 1.1 / (1. + std::sqrt(1.5));
    const auto positive = nps_xsec::solve_positive_response(ill_conditioned, ill_y, ill_variance, {.5});
    const auto joint = nps_xsec::solve_joint_positive_response(ill_conditioned, ill_y, ill_variance,
                                                              {.5}, {{{0, 1, 2}}});
    for (const auto* solution : {&positive, &joint}) {
      require(solution->boundary_active && close(solution->fit.parameters[0], 1.-expected_lt, 2e-6) &&
              close(solution->fit.parameters[1], expected_lt, 2e-6) &&
              close(solution->fit.parameters[2], .2, 2e-6),
              "ill-conditioned all-angle boundary projection failed");
      require(solution->minima[0] >= -solution->feasibility_tolerances[0],
              "ill-conditioned boundary is not feasible");
      require(solution->fit.chi2 < 1e-12, "ill-conditioned boundary has wrong objective");
    }
    rejects([&] { nps_xsec::solve_positive_response(ill_conditioned, ill_y, ill_variance, {.5}, 1e-8); },
            "positivity bypassed the relative rank threshold");
  }

  auto degenerate = migration;
  for (auto& row : degenerate) row[3] = row[0];
  rejects([&] { nps_xsec::solve_weighted_response(degenerate, y, variance); }, "duplicate column accepted");
  // A nearly dependent column must fail the configured relative threshold,
  // while the same finite direction is retained when explicitly resolvable.
  for (std::size_t r = 0; r < degenerate.size(); ++r)
    degenerate[r][3] = migration[r][0] + 1e-7 * migration[r][3];
  rejects([&] { nps_xsec::solve_weighted_response(degenerate, y, variance, 1e-6); },
          "configured relative rank threshold ignored");
  const auto near_fit = nps_xsec::solve_weighted_response(degenerate, predict(degenerate, truth), variance, 1e-10);
  require(near_fit.rank == 6 && near_fit.condition > 1e7, "resolvable near-dependent direction discarded");
  for (auto& row : degenerate) row[3] = 0.;
  rejects([&] { nps_xsec::solve_weighted_response(degenerate, y, variance); }, "unsupported column accepted");
  auto invalid_variance = variance;
  invalid_variance[0] = 0.;
  rejects([&] { nps_xsec::solve_weighted_response(migration, y, invalid_variance); }, "zero variance accepted");
  auto invalid_response = migration;
  invalid_response[0][0] = std::numeric_limits<double>::quiet_NaN();
  rejects([&] { nps_xsec::solve_weighted_response(invalid_response, y, variance); }, "NaN response accepted");
  auto square = identity;
  square.resize(6);
  rejects([&] { nps_xsec::solve_weighted_response(square, truth, std::vector<double>(6, 1.)); }, "zero ndf accepted");
  // Finite-MC variance is the sum of squared event predictions. Validate its
  // independent outer-product implementation and ensure invalid covariance
  // cannot turn into a spurious zero uncertainty through max(0, NaN).
  std::vector<nps_xsec::ResponseCell> cells(1);
  cells[0].add({1., 2., -3.});
  cells[0].add({2., -1., 0.5});
  const std::vector<double> parameters{3., -0.4, 0.2};
  const double expected_mc=std::pow(3.-0.8-0.6,2)+std::pow(6.+0.4+0.1,2);
  require(close(nps_xsec::mc_prediction_variance(cells,{0},parameters),expected_mc),
          "MC covariance does not match squared event predictions");
  cells[0].covariance[0]=std::numeric_limits<double>::quiet_NaN();
  rejects([&] { nps_xsec::mc_prediction_variance(cells,{0},parameters); }, "NaN MC covariance hidden");
  cells[0].covariance.fill(0.);
  cells[0].covariance[0]=-1.;
  rejects([&] { nps_xsec::mc_prediction_variance(cells,{0},parameters); }, "negative MC covariance hidden");
  cells[0].covariance[0]=1e200;
  rejects([&] { nps_xsec::mc_prediction_variance(cells,{0},{1e200,0.,0.}); }, "MC covariance overflow hidden");
  std::cout << "PASS: identity and two-bin migration recovery; full covariance/SVD factor; 240-order unit rescaling; "
               "ill-conditioned single/joint positivity; "
               "rank/support/input rejection; finite-MC covariance checks. condition=" << fit.condition << " chi2=" << fit.chi2 << '\n';
  return 0;
}
