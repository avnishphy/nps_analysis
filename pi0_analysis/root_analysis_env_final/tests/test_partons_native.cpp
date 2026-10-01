// Optional integration smoke test: requires native PARTONS/ElementaryUtils/
// NumA++ and the same libraries/properties used by run_xsec_pipeline.sh.
// Low MC statistics keep this test short; coefficients are NOT compared to
// fixed values because the GK convolution is stochastic. The algebra-only
// regression in test_pi0_response_conventions.cpp supplies exact expectations.
#include "../src/xsec_extract/partons_pi0_projection.h"
#include <algorithm>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    (void)argc;
    try {
        nps_partons_pi0::Model model;
        model.initialize(argv[0], 1000, 10000);
        const auto p=model.predict(5.8, .58, -1.0, 10.538);
        if(!p.valid) throw std::runtime_error("Native PARTONS point unexpectedly invalid");
        const double pi=nps_pi0_conventions::pi;
        const double klt=std::sqrt(2*p.epsilon*(1+p.epsilon));
        const std::array<double,3> phi{0,pi/2,pi};
        double max_relative_error=0;
        for(size_t i=0;i<phi.size();++i) {
            // Restore nb/GeV2 and the electron prefactor from the extracted
            // phi-integrated responses. FULL Hand Gamma/(2*pi) is required.
            const double reconstructed=p.electron_flux_xbq2/(2*pi)*1e9*
                (p.sigmaU+klt*p.sigmaLT*std::cos(phi[i])+p.epsilon*p.sigmaTT*std::cos(2*phi[i]));
            const double error=std::abs(reconstructed-p.electron_nb[i])/std::max(1e-30,std::abs(p.electron_nb[i]));
            max_relative_error=std::max(max_relative_error,error);
            if(error>1e-12) throw std::runtime_error("Native electron-observable round-trip failed");
        }
        std::cout<<std::setprecision(17)<<"PASS native GK smoke: epsilon="<<p.epsilon
                 <<", full_Hand_flux="<<p.electron_flux_xbq2<<", U_LT_TT_nb_per_GeV2="
                 <<p.sigmaU*1e9<<","<<p.sigmaLT*1e9<<","<<p.sigmaTT*1e9
                 <<", max_roundtrip_relative_error="<<max_relative_error<<'\n';
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"FAIL native GK smoke: "<<error.what()<<'\n';return 1;
    }
}
