#ifndef NPS_XSEC_LINEAR_SOLVER_H
#define NPS_XSEC_LINEAR_SOLVER_H

// Weighted forward-response fit, independent of event selection and binning.
// The caller supplies X[reconstructed row][truth coefficient], the measured
// reconstructed yields y, and their variances. Consequently migrations are
// fitted simultaneously; this routine never divides a yield by a diagonal
// acceptance or assumes reconstructed and generated bins coincide.
#include <TDecompSVD.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace nps_xsec {

struct LinearSolution {
  std::vector<double> parameters;       // Original physical coefficient units.
  std::vector<double> covariance;       // Row-major n_parameters squared.
  std::vector<double> singular_values;  // Whitened, unit-column-norm response.
  std::size_t rank = 0;
  double condition = 0.;               // Largest/smallest scaled singular value.
  double chi2 = 0.;
  int ndf = 0;                         // Included rows minus fitted parameters.
  // Row-major dense F = D^-1 V S^-1, so covariance = F F^T and
  // X_whitened F = U. Retained directly from SVD: factoring the rounded
  // covariance squares the condition number and can lose valid directions.
  std::vector<double> inverse_information_factor;
};

// Solve min sum_r (y_r - sum_j X_rj a_j)^2 / variance_r.
// Variances are treated as supplied, known measurement variances: covariance
// is (X^T W X)^-1, WITHOUT multiplication by chi2/ndf. This supports background-
// subtracted, potentially negative yields and signed interference coefficients.
// Response-MC statistics and correlations between measured rows are not included
// here; they require a separate uncertainty treatment by the caller.
//
// Numerical conditioning is assessed after whitening and column normalization,
// so a mere change of cross-section units cannot decide which physics parameter
// is identifiable. All singular directions must pass the relative rank test;
// no truncation, regularization, nonnegativity constraint, or silent row removal
// is allowed. Unsupported or degenerate truth coefficients abort the fit.
inline LinearSolution solve_weighted_response(
    const std::vector<std::vector<double>>& X,
    const std::vector<double>& y,
    const std::vector<double>& variance,
    double rank_tolerance = 1e-10) {
  const auto fail = [](const std::string& reason) {
    throw std::runtime_error("Weighted response fit: " + reason);
  };
  if (!std::isfinite(rank_tolerance) || rank_tolerance <= 0. || rank_tolerance >= 1.)
    fail("rank_tolerance must be finite and strictly between zero and one");
  if (X.empty() || X.front().empty()) fail("empty response matrix");
  const std::size_t nr = X.size(), np = X.front().size();
  if (y.size() != nr || variance.size() != nr) fail("yield/variance row-count mismatch");
  if (nr <= np) fail("number of rows must exceed number of fitted parameters");
  if (nr > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      np > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      np > std::numeric_limits<std::size_t>::max() / np)
    fail("matrix dimensions exceed supported index range");

  // Use long-double arithmetic during whitening/scaling. In particular, an
  // otherwise valid response with very large/small physical units should not
  // overflow while computing squared column norms in double precision.
  std::vector<long double> sigma(nr), column_norm(np, 0.L);
  for (std::size_t r = 0; r < nr; ++r) {
    if (X[r].size() != np) fail("ragged response matrix at row " + std::to_string(r));
    if (!std::isfinite(y[r])) fail("nonfinite yield at row " + std::to_string(r));
    if (!std::isfinite(variance[r]) || variance[r] <= 0.)
      fail("variance must be finite and positive at row " + std::to_string(r));
    sigma[r] = std::sqrt(static_cast<long double>(variance[r]));
    for (std::size_t j = 0; j < np; ++j) {
      if (!std::isfinite(X[r][j]))
        fail("nonfinite response at row " + std::to_string(r) + ", column " + std::to_string(j));
      column_norm[j] = std::hypot(column_norm[j], static_cast<long double>(X[r][j]) / sigma[r]);
    }
  }
  for (std::size_t j = 0; j < np; ++j)
    if (!std::isfinite(column_norm[j]) || column_norm[j] <= 0.)
      fail("unsupported truth coefficient (zero/invalid response column) " + std::to_string(j));

  TMatrixD scaled(static_cast<int>(nr), static_cast<int>(np));
  for (std::size_t r = 0; r < nr; ++r)
    for (std::size_t j = 0; j < np; ++j)
      scaled(static_cast<int>(r), static_cast<int>(j)) =
          static_cast<double>((static_cast<long double>(X[r][j]) / sigma[r]) / column_norm[j]);

  TDecompSVD svd(scaled);
  if (!svd.Decompose()) fail("ROOT SVD decomposition failed");
  const TVectorD& singular = svd.GetSig();
  const TMatrixD& U = svd.GetU();
  const TMatrixD& V = svd.GetV();
  LinearSolution result;
  result.singular_values.resize(np);
  double largest = 0., smallest = std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k < np; ++k) {
    const double s = singular(static_cast<int>(k));
    if (!std::isfinite(s) || s < 0.) fail("invalid SVD singular value");
    result.singular_values[k] = s;
    if (s > largest) largest = s;
    if (s < smallest) smallest = s;
  }
  for (double s : result.singular_values)
    if (s > rank_tolerance * largest) ++result.rank;
  if (result.rank != np) {
    std::ostringstream threshold;
    threshold << std::scientific << rank_tolerance;
    fail("rank-deficient migration response: rank " + std::to_string(result.rank) +
         " of " + std::to_string(np) + "; relative SVD threshold " +
         threshold.str() +
         ". Merge unsupported truth bins or improve response statistics; no fit was returned.");
  }
  result.condition = largest / smallest;

  // X_whitened = U S V^T D, with D the physical column norms.
  // a = D^-1 V S^-1 U^T y_whitened. Implement explicitly so ROOT's default
  // pseudoinverse tolerance cannot silently discard a singular direction.
  std::vector<long double> projected(np, 0.L);
  for (std::size_t k = 0; k < np; ++k) {
    for (std::size_t r = 0; r < nr; ++r)
      projected[k] += static_cast<long double>(U(static_cast<int>(r), static_cast<int>(k))) * y[r] / sigma[r];
    projected[k] /= result.singular_values[k];
  }
  result.parameters.resize(np);
  result.covariance.resize(np * np);
  result.inverse_information_factor.resize(np * np);
  for (std::size_t j = 0; j < np; ++j) {
    long double value = 0.L;
    for (std::size_t k = 0; k < np; ++k)
      value += static_cast<long double>(V(static_cast<int>(j), static_cast<int>(k))) * projected[k];
    result.parameters[j] = static_cast<double>(value / column_norm[j]);
    if (!std::isfinite(result.parameters[j])) fail("fitted parameter exceeds finite double range");
    for (std::size_t k = 0; k < np; ++k) {
      const double factor = static_cast<double>(
          (static_cast<long double>(V(static_cast<int>(j), static_cast<int>(k))) /
           result.singular_values[k]) / column_norm[j]);
      if (!std::isfinite(factor)) fail("inverse-information factor exceeds finite double range");
      result.inverse_information_factor[j * np + k] = factor;
    }
    for (std::size_t l = 0; l <= j; ++l) {
      long double cov = 0.L;
      for (std::size_t k = 0; k < np; ++k) {
        const long double vk_j = V(static_cast<int>(j), static_cast<int>(k));
        const long double vk_l = V(static_cast<int>(l), static_cast<int>(k));
        const long double s = result.singular_values[k];
        cov += (vk_j / s) * (vk_l / s);
      }
      const double physical_cov = static_cast<double>((cov / column_norm[j]) / column_norm[l]);
      if (!std::isfinite(physical_cov) || (j == l && physical_cov <= 0.))
        fail("parameter covariance exceeds representable double range");
      result.covariance[j * np + l] = physical_cov;
      result.covariance[l * np + j] = physical_cov;
    }
  }
  long double chi2 = 0.L;
  for (std::size_t r = 0; r < nr; ++r) {
    long double prediction = 0.L;
    for (std::size_t j = 0; j < np; ++j)
      prediction += static_cast<long double>(X[r][j]) * result.parameters[j];
    const long double residual = (static_cast<long double>(y[r]) - prediction) / sigma[r];
    chi2 += residual * residual;
  }
  result.chi2 = static_cast<double>(chi2);
  if (!std::isfinite(result.chi2)) fail("chi-square exceeds finite double range");
  result.ndf = static_cast<int>(nr - np);
  return result;
}

}  // namespace nps_xsec
#endif  // NPS_XSEC_LINEAR_SOLVER_H
