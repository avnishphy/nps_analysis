#include "../xsec_scaled_poisson_stat.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

int main() {
    using nps_xsec::scaled_poisson_deviance;
    auto near=[](double a,double b,double tol=1e-12) {
        return std::abs(a-b)<tol*std::max(1.,std::abs(b));
    };
    assert(near(scaled_poisson_deviance(0.,0.1,1.),0.2));
    assert(near(scaled_poisson_deviance(0.,1.,1.),2.));
    assert(near(scaled_poisson_deviance(3.,2.,0.5),
                4.*(2.-3.+3.*std::log(3./2.))));
    assert(near(scaled_poisson_deviance(5.,5.,0.3),0.));
    const double y=100.,mu=101.,s=0.2;
    const double gaussian=(y-mu)*(y-mu)/(s*y);
    assert(std::abs(scaled_poisson_deviance(y,mu,s)/gaussian-1.)<0.01);
    bool invalid=false;
    try { (void)scaled_poisson_deviance(1.,-1.,1.); }
    catch (const std::runtime_error&) { invalid=true; }
    assert(invalid);
}
