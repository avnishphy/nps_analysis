#pragma once
#include <array>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <string>

namespace nps_xsec {
struct XSecModelKinematics {
    double Q2=0, W2=0, t=0, tprime=0, tau=0, theta_cm=0, epsilon=0;
};
struct StructureFunctions {
    double sigma_U=0, sigma_LT=0, sigma_TT=0, sigma_T=0, sigma_L=0;
    std::array<double,3> unseparated() const {return {sigma_U,sigma_LT,sigma_TT};}
};
namespace sigparam2021 {
// Authoritative source: /u/group/nps/singhav/simc_gfortran_updated/physics_pion.f
// sig_param_2021/exclfit, lines 738-855; coefficients fixed, never fitted here.
inline constexpr std::array<double,17> pp{{
    1.60077,-0.01523,37.08142,-4.11060,23.26192,0.00983,
    0.87073,-5.77115,-271.08678,0.13766,-0.00855,0.27885,
    -1.13212,-1.50415,-6.34766,0.55769,-0.01709}};
inline constexpr std::array<double,17> pm{{
    1.75169,0.11144,47.35877,-4.69434,1.60552,0.00800,
    0.44194,-2.29188,-41.67194,0.69475,0.02527,-0.50178,
    -1.22825,-1.16878,5.75825,-1.00355,0.05055}};
inline constexpr double source_mp=0.938, source_pi=3.1415928, Q0sq=1.;
inline double denominator(double d,double scale,const char* name) {
    if(!std::isfinite(d) || std::abs(d)<=1e-12*std::max(1.,scale))
        throw std::domain_error(std::string("SigParam2021 singular denominator: ")+name);
    return d;
}
inline double damping(double q,const std::array<double,17>& p) {
    return 1./denominator(1+p[0]*q+p[1]*q*q,1+std::abs(p[0]*q)+std::abs(p[1]*q*q),"G(Q2)");
}
inline double kernel(double q,double at,const std::array<double,17>& p,bool charged) {
    const double g=damping(q,p);
    return q*g*g*(charged?at/((at+.02)*(at+.02)):1./Q0sq);
}
inline StructureFunctions components(const XSecModelKinematics& x,
                                    const std::array<double,17>& p,bool charged=false) {
    for(double v:{x.Q2,x.W2,x.t,x.theta_cm,x.epsilon})
        if(!std::isfinite(v))throw std::domain_error("Nonfinite SigParam2021 kinematics");
    if(!(x.Q2>0 && x.W2>0 && x.epsilon>=0 && x.epsilon<=1 && x.theta_cm>=0 && x.theta_cm<=std::acos(-1.)))
        throw std::domain_error("Invalid SigParam2021 Q2/W2/epsilon/theta");
    const double q=x.Q2,s=x.W2,w=std::sqrt(s),at=std::abs(x.t),st=std::sin(x.theta_cm);
    const double dl=denominator(std::pow(s,p[10])+std::pow(w,p[16]),1.,"L W2"),
                 dt=denominator(std::pow(s,p[11])+std::pow(w,p[15]),1.,"T W2"),
                 dlt=denominator((1+p[9]*q)*std::pow(s,p[12]),1.,"LT Q2/W2"),
                 rw=denominator(s-source_mp*source_mp,s,"W2-Mp2");
    // SigParam2021 was fitted to charged-pion electroproduction. Its L term
    // contains |t|/(|t|+0.02)^2 * Q2*G(Q2)^2. For this provisional pi0 proxy,
    // remove that pole shape without forcing L=0: use (Q2/Q0sq)*G(Q2)^2.
    // G is empirical charged-fit damping, NOT a pi0 electromagnetic or
    // transition form factor. Q0sq=1 GeV^2 is only a dimensional reference.
    // The charged=true path exists solely for original-Fortran validation.
    StructureFunctions f;
    f.sigma_L=(p[2]+p[14]/q)*kernel(q,at,p,charged)*std::exp(p[3]*at)/dl;
    f.sigma_T=p[4]/q*std::exp(p[5]*q*q)*std::exp(p[13]*at)/dt;
    f.sigma_LT=p[6]*std::exp(p[7]*at)*st/dlt;
    f.sigma_TT=p[8]/(1+q)*std::exp(-7*at)*st*st;
    // Components in microbarn/MeV^2. No angular 1/(2*pi) here: that factor
    // already appears exactly once in the detector response event basis.
    const double norm=8.539/(rw*rw)/1e6;
    f.sigma_T*=norm;f.sigma_L*=norm;f.sigma_LT*=norm;f.sigma_TT*=norm;
    f.sigma_U=f.sigma_T+x.epsilon*f.sigma_L;
    for(double v:{f.sigma_U,f.sigma_T,f.sigma_L,f.sigma_LT,f.sigma_TT})
        if(!std::isfinite(v))throw std::domain_error("Nonfinite SigParam2021 exponential/component");
    return f;
}
inline StructureFunctions baseline(const XSecModelKinematics& x) {
    if(x.W2<4.)throw std::domain_error("SigParam2021 pi0 proxy requires W >= 2 GeV; MAID region unsupported");
    const auto a=components(x,pp),b=components(x,pm);
    return { .5*(a.sigma_U+b.sigma_U), .5*(a.sigma_LT+b.sigma_LT),
             .5*(a.sigma_TT+b.sigma_TT), .5*(a.sigma_T+b.sigma_T), .5*(a.sigma_L+b.sigma_L)};
}
inline XSecModelKinematics kinematics(double q,double w,double t,double tp,double eps,
                                      double mp,double mpi) {
    if(!(q>0 && w>mp+mpi && tp<=1e-10))throw std::domain_error("Invalid generated exclusive model kinematics");
    const double eg=(w*w-mp*mp-q)/(2*w),pg=std::sqrt(eg*eg+q),
                 epi=(w*w+mpi*mpi-mp*mp)/(2*w),p2=epi*epi-mpi*mpi;
    if(!(p2>0))throw std::domain_error("Invalid generated two-body threshold");
    const double ct=(t+q-mpi*mpi+2*eg*epi)/(2*pg*std::sqrt(p2));
    if(!std::isfinite(ct) || ct < -1-1e-5 || ct > 1+1e-5)
        throw std::domain_error("Invalid generated theta_cm from exclusive two-body invariants");
    return {q,w*w,t,tp,-tp,std::acos(std::clamp(ct,-1.,1.)),eps};
}
}

inline constexpr size_t kModelParameters=4;
struct ModelEvaluation {
    std::array<double,3> value{};
    std::array<std::array<double,kModelParameters>,3> jacobian{};
};
// One production model: a pivoted U slope and U/LT/TT normalizations.
// The baseline is fixed; the detector solver also fits migration nuisances.
struct XSecModelDefinition {
    std::string id,label;
    std::array<std::string,kModelParameters> names;
    std::array<double,kModelParameters> initial,lower,upper;
    StructureFunctions (*baseline)(const XSecModelKinematics&);
};
inline const XSecModelDefinition& xsec_model() {
    static const double inf=std::numeric_limits<double>::infinity();
    static const XSecModelDefinition model{"sigparam2021_pi0","SigParam2021-inspired pi0", 
        {"N_U","DeltaB_U","N_LT","N_TT"},{1.,0.,1.,1.},
        {1e-12,-20.,-inf,-inf},{inf,20.,inf,inf},sigparam2021::baseline};
    return model;
}
inline ModelEvaluation evaluate_cached_model(const std::array<double,3>& baseline,const double* p,
                                             double tau,double tau0) {
    ModelEvaluation f;
    const double dt=tau-tau0,shape=std::exp(-p[1]*dt);
    f.value={p[0]*shape*baseline[0],p[2]*baseline[1],p[3]*baseline[2]};
    f.jacobian[0][0]=shape*baseline[0];f.jacobian[0][1]=-dt*f.value[0];
    f.jacobian[1][2]=baseline[1];f.jacobian[2][3]=baseline[2];
    for(int k=0;k<3;++k) {
        if(!std::isfinite(f.value[k]) || !std::isfinite(baseline[k]))throw std::domain_error("Nonfinite normalized model");
    }
    return f;
}
struct ModelEvent {
    int row=-1,block=-1,truth_phi=-1;
    double weight=0,phi=0;
    XSecModelKinematics kinematics;
    StructureFunctions baseline;
    std::array<double,3> basis{};
};
}
