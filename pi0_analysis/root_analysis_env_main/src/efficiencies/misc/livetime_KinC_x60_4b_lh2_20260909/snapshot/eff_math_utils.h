#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>

namespace effstuff {

inline double safe_div(double numerator, double denominator) {
  if (!std::isfinite(denominator) || denominator == 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return numerator / denominator;
}

inline std::string trim_copy(const std::string& input) {
  size_t first = 0;
  while (first < input.size() && std::isspace(static_cast<unsigned char>(input[first])) != 0) {
    ++first;
  }

  if (first == input.size()) {
    return "";
  }

  size_t last = input.size() - 1;
  while (last > first && std::isspace(static_cast<unsigned char>(input[last])) != 0) {
    --last;
  }

  return input.substr(first, last - first + 1);
}

inline std::string to_lower_copy(const std::string& input) {
  std::string out = input;
  std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return out;
}

inline bool iequals(const std::string& a, const std::string& b) {
  return to_lower_copy(a) == to_lower_copy(b);
}

inline double binomial_efficiency_error(long long numerator,
                                        long long denominator) {
  if (denominator <= 0 || numerator < 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double p = static_cast<double>(numerator) /
                   static_cast<double>(denominator);
  if (!std::isfinite(p)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double variance = std::max(0.0, p * (1.0 - p)) /
                          static_cast<double>(denominator);
  return std::sqrt(variance);
}

}  // namespace effstuff
