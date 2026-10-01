#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Scaled-Poisson deviance for a weighted compound-Poisson approximation.
// Y is a weighted sum, never an integer Poisson observation. At mu=Y,
// D=(mu-Y)^2/(s*Y)+higher orders and s*Y=sum(w^2).
namespace nps_xsec {
inline double scaled_poisson_deviance(double y, double mu, double s) {
    if (!(std::isfinite(y) && y >= 0. && std::isfinite(mu) && mu > 0. &&
          std::isfinite(s) && s > 0.))
        throw std::runtime_error("scaled-Poisson deviance requires Y>=0, mu>0, s>0");
    if (y == 0.) return 2.*mu/s;
    const double delta=(mu-y)/y;
    const double bracket=(delta > -0.5 && delta < 0.5)
        ? y*(delta-std::log1p(delta))
        : mu-y+y*std::log(y/mu);
    return 2.*std::max(0.,bracket)/s;
}
}
