#include "xsec_sigparam2021_pi0_model.h"
#include <iostream>
#include <iomanip>
int main() {
    double q,w,t,tp,eps;
    std::cout<<std::setprecision(17);
    while(std::cin>>q>>w>>t>>tp>>eps) {
        const auto x=nps_xsec::sigparam2021::kinematics(q,w,t,tp,eps,.9382720813,.1349768);
        const auto f=nps_xsec::sigparam2021::baseline(x);
        std::cout<<f.sigma_U<<' '<<f.sigma_LT<<' '<<f.sigma_TT<<'\n';
    }
}
