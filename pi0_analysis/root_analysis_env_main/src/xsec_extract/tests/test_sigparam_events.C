#include "xsec_proxy_solver.h"
#include <iostream>
#include <iomanip>

int main() {
    using namespace nps_xsec;
    ProxyProblem q;
    q.blocks={0,1,5,7};q.physical={true,true,false,false};q.fixed={false,false,false,true};q.tprime={-.1,-.3,0,0};q.epsilon={.7,.7,.7,.7};
    q.design.assign(72,std::vector<double>(12,0));q.events.resize(72);q.reporting.resize(4);
    q.variance.assign(72,1e-5);q.scale.assign(72,.01);
    std::vector<ModelEvent> cache;cache.reserve(72*12);
    const std::vector<double> truth{1.2,.4,-.8,.5,4e-8,-2e-9,1e-9};
    for(int r=0;r<72;++r)for(int n=0;n<12;++n) {
        ModelEvent e;e.row=r;const int b=n==11?2:n==10?3:n%2;e.block=q.blocks[b];
        e.phi=2*std::acos(-1.)*(r%12+.5)/12.+.01*n;
        e.weight=1e6*(b==2?.3*(1+.3*std::sin(r)):.8+.03*n);
        const double tau=.035+.07*(r/12)+.002*n;
        e.kinematics={3.7+.04*n,7.8+.05*n,-.13-tau,-tau,tau,.1+.02*n,.6+.003*n};
        if(b!=2)e.baseline=xsec_model().baseline(e.kinematics);
        const double eps=e.kinematics.epsilon;
        e.basis={e.weight/(2*std::acos(-1.)),e.weight*std::sqrt(2*eps*(1+eps))*std::cos(e.phi)/(2*std::acos(-1.)),e.weight*eps*std::cos(2*e.phi)/(2*std::acos(-1.))};
        for(int k=0;k<3;++k)q.design[r][3*b+k]+=e.basis[k];
        cache.push_back(e);q.events[r].push_back(&cache.back());q.reporting[b].push_back(&cache.back());
    }
    q.set_pivot();if(q.parameter_count()!=7)return 6;q.y=q.prediction(truth);
    // Zero slope reproduces all three old components exactly at every event.
    for(const auto& e:cache) {
        auto zero=truth;zero[1]=0.;
        const auto f=evaluate_cached_model(e.baseline.unseparated(),zero.data(),e.kinematics.tau,q.tau0);
        if(f.value[0]!=zero[0]*e.baseline.sigma_U || f.value[1]!=zero[2]*e.baseline.sigma_LT || f.value[2]!=zero[3]*e.baseline.sigma_TT)return 4;
        // Moving a fixed pivot changes only the definition of N_U.
        auto shifted=truth;const double delta=.1;shifted[0]*=std::exp(-truth[1]*delta);
        const auto a=evaluate_cached_model(e.baseline.unseparated(),truth.data(),e.kinematics.tau,q.tau0);
        const auto b=evaluate_cached_model(e.baseline.unseparated(),shifted.data(),e.kinematics.tau,q.tau0+delta);
        if(std::abs(a.value[0]-b.value[0])>1e-14*std::max(1e-30,std::abs(a.value[0])))return 5;
    }
    double derivative_error=0,mc_error=0;
    for(size_t r=0;r<q.events.size();++r) {
        const auto v=q.event_row(r,truth);double square=0;
        for(const auto* e:q.events[r]) {
            const auto f=q.event_value(*e,truth);double x=0;
            for(int k=0;k<3;++k)x+=e->basis[k]*f.value[k];square+=x*x;
        }
        mc_error=std::max(mc_error,std::abs(square-v.mc_variance));
        for(size_t j=0;j<truth.size();++j) {
            auto a=truth,b=truth;double h=1e-5*std::max(std::abs(truth[j]),j<kModelParameters?.1:1e-9);
            a[j]+=h;b[j]-=h;
            const double numeric=(q.event_row(r,a).prediction-q.event_row(r,b).prediction)/(2*h);
            derivative_error=std::max(derivative_error,std::abs(numeric-v.jacobian[j])/std::max(1e-8,std::abs(v.jacobian[j])));
        }
    }
    if(derivative_error>2e-5 || mc_error>1e-15)return 1;
    ProxyOptions o;o.tolerance=1e-10;
    auto seed=proxy_seed(q,o);double closure=0;
    for(bool poisson:{false,true}) {
        q.poisson=poisson;
        const auto fit=minimize_proxy(q,o,seed);
        for(size_t j=0;j<truth.size();++j)closure=std::max(closure,std::abs(fit.parameters[j]-truth[j])/std::max(std::abs(truth[j]),1e-10));
        std::cout<<"EVENT_CLOSURE poisson="<<poisson<<" status="<<fit.status<<" cov="<<fit.covariance_status<<" objective="<<fit.objective<<" EDM="<<fit.edm<<'\n';
        if(!fit.converged || fit.status!=0 || fit.covariance_status!=3 || fit.objective>1e-5)return 2;
        seed=fit.parameters;
    }
    std::cout<<std::setprecision(17)<<"EVENT_PASS max_derivative_relative="<<derivative_error<<" mc_square_error="<<mc_error<<" max_parameter_relative="<<closure<<'\n';
    return closure<.01?0:3;
}
