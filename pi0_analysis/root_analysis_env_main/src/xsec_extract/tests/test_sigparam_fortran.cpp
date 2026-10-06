#include "xsec_sigparam2021_pi0_model.h"
#include <iostream>
#include <iomanip>
#include <fstream>
// SIMC's -ff2c plus -fdefault-real-8 promotes the function return ABI to
// REAL(16); explicit REAL*8 arguments remain eight-byte references.
extern "C" __float128 sig_param_2021_(double*,double*,double*,double*,double*,double*,int*,int*) asm("sig_param_2021__");
int main(int argc,char** argv) {
    double worst=0,worst_pi=0;int count=0;
    // Direct LT sign isolation: U and TT cancel between phi=0 and phi=pi.
    for(double tp:{-.65,-.325,-.065})for(int charge:{0,1}) {
        double q=4.,s=8.,t=-.16+tp,th=.3,eps=.515,zero=0.,pi=std::acos(-1.);int neutral=0;
        const double f0=sig_param_2021_(&th,&zero,&t,&q,&s,&eps,&charge,&neutral);
        const double fpi=sig_param_2021_(&th,&pi,&t,&q,&s,&eps,&charge,&neutral);
        const nps_xsec::XSecModelKinematics x{q,s,t,tp,-tp,th,eps};
        const auto f=nps_xsec::sigparam2021::components(x,charge==0?nps_xsec::sigparam2021::pp:nps_xsec::sigparam2021::pm,true);
        const double isolated=(f0-fpi)*nps_xsec::sigparam2021::source_pi/std::sqrt(2*eps*(1+eps));
        if(std::abs(isolated-f.sigma_LT)>1e-12*std::abs(f.sigma_LT))return 7;
        std::cout<<std::setprecision(17)<<"LT_SIGN charge="<<charge<<" t="<<t<<" Fortran_LT="<<isolated<<" CPP_LT="<<f.sigma_LT<<'\n';
    }
    auto compare=[&](double q,double s,double t,double theta,double phi,double eps,int charge) {
        int neutral=0;double th=theta,ph=phi,tt=t,qq=q,ss=s,ee=eps;
        const double original=sig_param_2021_(&th,&ph,&tt,&qq,&ss,&ee,&charge,&neutral);
        nps_xsec::XSecModelKinematics x{q,s,t,0,0,theta,eps};
        const auto f=nps_xsec::sigparam2021::components(x,charge==0?nps_xsec::sigparam2021::pp:nps_xsec::sigparam2021::pm,true);
        const double raw=f.sigma_U+std::sqrt(2*eps*(1+eps))*std::cos(phi)*f.sigma_LT+eps*std::cos(2*phi)*f.sigma_TT;
        const double port=raw/(2*nps_xsec::sigparam2021::source_pi);
        const double scale=std::max(std::abs(original),1e-20);
        worst=std::max(worst,std::abs(port-original)/scale);
        worst_pi=std::max(worst_pi,std::abs(raw/(2*std::acos(-1.))-original)/scale);++count;
    };
    for(double q:{1.2,3.3,4.,4.7,7.})for(double s:{4.1,7.5,11.})
    for(double t:{-.01,-.15,-.5,-1.1,.5})for(double theta:{0.,.2,.7})
    for(double phi:{0.,.6,1.8,3.1,5.8})for(double eps:{.25,.55,.8})for(int charge:{0,1})
        compare(q,s,t,theta,phi,eps,charge);
    if(argc>1) {
        std::ifstream input(argv[1]);if(!input)throw std::runtime_error("Missing generated-kinematics input");
        double q,s,t,theta,phi,eps;
        while(input>>q>>s>>t>>theta>>phi>>eps)for(int charge:{0,1})compare(q,s,t,theta,phi,eps,charge);
        if(!input.eof())throw std::runtime_error("Malformed generated-kinematics input");
    }
    std::cout<<std::setprecision(17)<<"CHARGED_FORTRAN cases="<<count<<" max_relative="<<worst
             <<" extractor_exact_pi_max_relative="<<worst_pi<<'\n';
    if(worst>1e-10 || worst_pi>5e-8)return 1;
    nps_xsec::XSecModelKinematics x{4.,8.,0.,0.,0.,0.,.5};
    const auto neutral=nps_xsec::sigparam2021::baseline(x);
    if(!(neutral.sigma_L>0 && std::isfinite(neutral.sigma_L)))return 2;
    for(auto invalid:{0.,-1.}) {x.Q2=invalid;try{nps_xsec::sigparam2021::baseline(x);return 3;}catch(const std::domain_error&) {}}
    x.Q2=4.;x.W2=3.99;try{nps_xsec::sigparam2021::baseline(x);return 4;}catch(const std::domain_error&) {}
    const auto& pp=nps_xsec::sigparam2021::pp;
    const double pole=(-pp[0]-std::sqrt(pp[0]*pp[0]-4*pp[1]))/(2*pp[1]);
    try{nps_xsec::sigparam2021::damping(pole,pp);return 5;}catch(const std::domain_error&) {}
    x.Q2=1e4;x.W2=8.;try{nps_xsec::sigparam2021::baseline(x);return 6;}catch(const std::domain_error&) {}
    std::cout<<"PASS: charged normalization, sign, W2, theta/phi; smooth nonzero neutral L; Q2 guards\n";
}
