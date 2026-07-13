#pragma once

#include <cmath>
#include <limits>

#include "eff_math_utils.h"
// this header file calculates the livetime using EDTM 
namespace effstuff {

struct NewGenEdtmAccumulator {
  long long numerator_events = 0;
  double denominator_scaler_edtm = 0.0;
};

inline double NewGen_EDTM_livetime(const NewGenEdtmAccumulator& acc,
                                   double ps_factor) {
  return safe_div(static_cast<double>(acc.numerator_events) * ps_factor,
                  acc.denominator_scaler_edtm) ;
}

inline double NewGen_EDTM_livetime_err(const NewGenEdtmAccumulator& acc,
                                       double ps_factor) {
  const double denominator = acc.denominator_scaler_edtm;
  if (!std::isfinite(ps_factor) || !std::isfinite(denominator) ||
      denominator <= 0.0 || acc.numerator_events < 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double numerator = static_cast<double>(acc.numerator_events);
  const double sigma_n = std::sqrt(std::max(0.0, numerator));
  const double sigma_d = std::sqrt(denominator);

  const double term_n = std::pow((ps_factor / denominator) * sigma_n, 2.0);
  const double term_d = std::pow((ps_factor * numerator / (denominator * denominator)) * sigma_d, 2.0);
  return std::sqrt(std::max(0.0, term_n + term_d));
}

}  // namespace effstuff
