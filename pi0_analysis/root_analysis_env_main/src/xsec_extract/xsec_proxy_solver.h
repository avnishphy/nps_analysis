#pragma once
#include "xsec_sigparam2021_pi0_model.h"
#include "xsec_response.h"
#include "xsec_scaled_poisson_stat.h"
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <memory>
#include <numeric>

namespace nps_xsec {
struct ProxyOptions {
    std::string fit_strategy="joint_minuit", staged_trace;
    std::vector<size_t> staged_order; // Empty = natural order; diagnostic permutations.
    std::array<double,kModelParameters> initial{{1.,0.,1.,1.}};
    bool fix_u_slope=false; // Validation route: reproduce three normalizations.
    int starts=6, max_calls=250000, max_iterations=50000;
    // Minuit EDM target is 0.001*tolerance*ErrorDef. Avoid demanding EDM
    // below numerical noise at a positivity boundary.
    double tolerance=1e-6;
};
struct ProxyProblem {
    double tau0=0.; // Computed once from positive physical response weights.
    void set_pivot() {
        double sum=0.,weighted=0.;
        for(size_t b=0;b<blocks.size();++b)if(physical[b])for(const auto* e:reporting[b]) {
            if(!(std::isfinite(e->weight) && e->weight>=0 && std::isfinite(e->kinematics.tau)))
                throw std::runtime_error("Invalid physical response weight/kinematics for pivot");
            sum+=e->weight;weighted+=e->weight*e->kinematics.tau;
        }
        if(!(sum>0))throw std::runtime_error("No positive physical response for pivot");
        tau0=weighted/sum;
    }
    // Row-local pointers into an immutable accepted-event cache. Reporting
    // averages use all accepted events, independently of objective row masks.
    std::vector<std::vector<const ModelEvent*>> events, reporting;
    bool is_physics(size_t b) const { return physical.at(b); }
    bool is_fixed(size_t b) const { return fixed.empty()?false:fixed.at(b); }
    bool is_nuisance(size_t b) const { return !is_physics(b) && !is_fixed(b); }
    size_t nuisance_index(size_t b) const {
        size_t index=kModelParameters;
        for(size_t i=0;i<b;++i)if(is_nuisance(i))index+=3;
        return index;
    }
    ModelEvaluation event_value(const ModelEvent& e,const std::vector<double>& p) const {
        const auto it=std::find(blocks.begin(),blocks.end(),e.block);
        if(it==blocks.end())throw std::logic_error("Inactive model event");
        const size_t b=it-blocks.begin();
        if(is_physics(b))return evaluate_cached_model(e.baseline.unseparated(),p.data(),e.kinematics.tau,tau0);
        if(is_fixed(b)) {
            ModelEvaluation f;f.value=e.baseline.unseparated();return f;
        }
        ModelEvaluation f;const size_t j=nuisance_index(b);
        for(int k=0;k<3;++k)f.value[k]=p[j+k];
        return f;
    }
    struct RowEvaluation {
        double prediction=0,mc_variance=0;
        std::array<double,3> components{};
        std::vector<double> jacobian;
    };
    RowEvaluation event_row(size_t r,const std::vector<double>& p) const {
        RowEvaluation out;out.jacobian.assign(parameter_count(),0.);
        for(const auto* e:events.at(r)) {
            const size_t b=std::find(blocks.begin(),blocks.end(),e->block)-blocks.begin();
            const auto f=event_value(*e,p);double v=0.;
            for(int k=0;k<3;++k) {
                const double c=e->basis[k]*f.value[k];v+=c;out.components[k]+=c;
                if(is_physics(b))for(size_t j=0;j<kModelParameters;++j)out.jacobian[j]+=e->basis[k]*f.jacobian[k][j];
                else if(is_nuisance(b))out.jacobian[nuisance_index(b)+k]+=e->basis[k];
            }
            out.prediction+=v;out.mc_variance+=v*v;
        }
        return out;
    }
    std::vector<std::vector<double>> design;
    std::vector<double> y, variance, scale, tprime, epsilon;
    // A block has exactly one treatment: fitted physics, fitted low-tprime
    // feed-in nuisance, or fixed event-level model feed-in.
    std::vector<bool> physical, fixed;
    std::vector<int> blocks;
    bool poisson=false, positive=false;
    size_t parameter_count() const {
        size_t count=kModelParameters;
        for(size_t b=0;b<blocks.size();++b)if(is_nuisance(b))count+=3;
        return count;
    }
    std::vector<std::string> names() const {
        const auto& names=xsec_model().names;
        std::vector<std::string> n(names.begin(),names.end());
        for(size_t b=0;b<blocks.size();++b) if(is_nuisance(b))
            for(const char* k:{"U","LT","TT"}) n.push_back("feedin_tprime_below_"+std::string(k));
        return n;
    }
    // Reference: independent generated-bin coefficients. Model: theta supplies
    // physical U/LT/TT; exterior (or excluded-group) blocks stay independent.
    std::vector<double> coefficients(const std::vector<double>& p,
                                     std::vector<std::vector<double>>* jac=nullptr) const {
        const size_t np=parameter_count();
        if(p.size()!=np) throw std::runtime_error("Proxy parameter size mismatch");
        std::vector<double> c(3*blocks.size());
        if(jac) jac->assign(c.size(),std::vector<double>(np,0.));
        size_t nuisance=kModelParameters;
        for(size_t b=0;b<blocks.size();++b) {
            if(is_physics(b)) {
                    double weight=0.;
                    for(const auto* e:reporting.at(b)) {
                        const auto f=event_value(*e,p);weight+=e->weight;
                        for(int k=0;k<3;++k) {
                            c[3*b+k]+=e->weight*f.value[k];
                            if(jac)for(size_t j=0;j<kModelParameters;++j)(*jac)[3*b+k][j]+=e->weight*f.jacobian[k][j];
                        }
                    }
                    if(!(weight>0))throw std::runtime_error("Empty physical event reporting bin");
                    for(int k=0;k<3;++k) {
                        c[3*b+k]/=weight;
                        if(jac)for(size_t j=0;j<kModelParameters;++j)(*jac)[3*b+k][j]/=weight;
                    }
            } else if(is_fixed(b)) {
                    double weight=0.;
                    for(const auto* e:reporting.at(b)) {
                        weight+=e->weight;
                        for(int k=0;k<3;++k)c[3*b+k]+=e->weight*e->baseline.unseparated()[k];
                    }
                    if(!(weight>0))throw std::runtime_error("Empty fixed model feed-in block");
                    for(int k=0;k<3;++k)c[3*b+k]/=weight;
            } else for(int k=0;k<3;++k) {
                c[3*b+k]=p[nuisance];
                if(jac) (*jac)[3*b+k][nuisance]=1.;
                ++nuisance;
            }
        }
        return c;
    }
    std::vector<double> prediction(const std::vector<double>& p) const {
        std::vector<double> mu(events.size());
        for(size_t r=0;r<mu.size();++r)mu[r]=event_row(r,p).prediction;
        return mu;
    }
    double objective(const std::vector<double>& p) const {
        const auto& model=xsec_model();
        for(size_t j=0;j<kModelParameters;++j)if(!std::isfinite(p[j]) || p[j]<model.lower[j] || p[j]>model.upper[j])return 1e100;
            if(positive || poisson)for(size_t b=0;b<blocks.size();++b)if(!is_fixed(b))for(const auto* e:reporting[b]) {
                const auto f=event_value(*e,p);
                if(minimum_response(f.value[0],f.value[1],f.value[2],is_physics(b)?e->kinematics.epsilon:epsilon[b])<0.)return 1e100;
            }
            double total=0.;
            for(size_t r=0;r<events.size();++r) {
                const double mu=event_row(r,p).prediction;
                if(!std::isfinite(mu) || (poisson && mu<=0.))return 1e100;
                total+=poisson?scaled_poisson_deviance(y[r],mu,scale[r]):(y[r]-mu)*(y[r]-mu)/variance[r];
            }
            return std::isfinite(total)?total:1e100;
    }
    std::vector<double> gradient(const std::vector<double>& p) const {
        std::vector<double> g(parameter_count(),0.);
        if(objective(p)>=1e99)return g;
            for(size_t r=0;r<events.size();++r) {
                const auto v=event_row(r,p);
                const double factor=poisson?2*(1-y[r]/v.prediction)/scale[r]:2*(v.prediction-y[r])/variance[r];
                for(size_t j=0;j<g.size();++j)g[j]+=factor*v.jacobian[j];
            }
            return g;
    }
};
struct ProxyResult {
    std::vector<double> parameters,covariance,errors,initial;
    int status=-1,covariance_status=-1,calls=0;
    bool converged=false,boundary=false;
    double edm=0,objective=1e100;
};
inline std::vector<double> proxy_seed(const ProxyProblem& q,const ProxyOptions& opts) {
    // Infer an amplitude scale from measured yield/unit-U response. This is
    // initialization only; independent extracted points never enter the objective.
    double y=std::accumulate(q.y.begin(),q.y.end(),0.), response=0.;
    for(const auto& row:q.design) for(size_t j=0;j<row.size();j+=3) response+=row[j];
    if(!(response>0.)) throw std::runtime_error("No unit-U response for proxy seed");
    const double amplitude=std::max(1e-20,y/response);
    std::vector<double> p(q.parameter_count(),0.);
    std::copy(opts.initial.begin(),opts.initial.end(),p.begin());
    if(opts.fix_u_slope)p[1]=0.;
    for(size_t j=kModelParameters;j<p.size();j+=3) p[j]=amplitude;
    return p;
}
inline ProxyResult minimize_proxy(const ProxyProblem& q,const ProxyOptions& opts,
                                  const std::vector<double>& seed) {
    const size_t np=q.parameter_count();
    if(q.y.size()<=np-size_t(opts.fix_u_slope)) throw std::runtime_error("Proxy fit has no positive nominal DOF");
    std::vector<double> unit(np,1.);
    for(size_t j=0;j<kModelParameters;++j)
        unit[j]=std::max({std::abs(seed[j]),std::abs(seed[0])*0.1,1e-20});
    unit[1]=1.; // GeV^-2 slope scale, independent of the zero initial value.
    for(size_t j=kModelParameters;j<np;++j)
        unit[j]=std::max({std::abs(seed[j]),std::abs(seed[kModelParameters+(j-kModelParameters)/3*3])*0.1,1e-20});
    auto physical=[&](const double* z) {
        std::vector<double> p(np);for(size_t j=0;j<np;++j)p[j]=z[j]*unit[j];return p;
    };
    auto objective=[&](const double* z) {
        try { return q.objective(physical(z)); } catch(const std::exception&) {return 1e100;}
    };
    std::unique_ptr<ROOT::Math::Minimizer> m(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
    if(!m) throw std::runtime_error("Minuit2 unavailable");
    auto gradient=[&](const double* z,double* out) {
        try {const auto g=q.gradient(physical(z));for(size_t j=0;j<np;++j)out[j]=g[j]*unit[j];}
        catch(const std::exception&){std::fill(out,out+np,0.);}
    };
    ROOT::Math::GradFunctor fn(objective,np,gradient);
    ROOT::Math::Functor constrained_fn(objective,np);
    // The rejection boundary is nonsmooth. Numerical derivatives let Minuit
    // probe that boundary; the smooth unconstrained Gaussian uses exact gradients.
    if(q.positive || q.poisson)m->SetFunction(constrained_fn);else m->SetFunction(fn);
    m->SetErrorDef(1.);
    m->SetMaxFunctionCalls(opts.max_calls);m->SetMaxIterations(opts.max_iterations);
    m->SetTolerance(opts.tolerance);m->SetPrintLevel(-1);m->SetStrategy(q.positive||q.poisson?1:2);
    const auto names=q.names();
    for(size_t j=0;j<np;++j) {
        const double start=seed[j]/unit[j];
        const auto& definition=xsec_model();
        if(j==1 && opts.fix_u_slope)m->SetFixedVariable(j,names[j],0.);
        else if(j<kModelParameters && std::isfinite(definition.lower[j]) && std::isfinite(definition.upper[j]))
            m->SetLimitedVariable(j,names[j],start,0.03,definition.lower[j]/unit[j],definition.upper[j]/unit[j]);
        else if(j<kModelParameters && std::isfinite(definition.lower[j]))
            m->SetLowerLimitedVariable(j,names[j],start,0.03,definition.lower[j]/unit[j]);
        else m->SetVariable(j,names[j],start,0.03);
    }
    ProxyResult out;out.initial=seed;
    out.converged=m->Minimize();
    m->Hesse();out.status=m->Status();out.covariance_status=m->CovMatrixStatus();
    out.edm=m->Edm();out.calls=m->NCalls();out.parameters=physical(m->X());
    out.objective=q.objective(out.parameters);out.errors.resize(np);out.covariance.resize(np*np);
    for(size_t i=0;i<np;++i) {
        out.errors[i]=m->Errors()[i]*unit[i];
        for(size_t j=0;j<np;++j) out.covariance[i*np+j]=m->CovMatrix(i,j)*unit[i]*unit[j];
        if(i<kModelParameters) {
            const auto& definition=xsec_model();
            out.boundary=out.boundary || (out.parameters[i]-definition.lower[i])/unit[i]<1e-6 ||
                (definition.upper[i]-out.parameters[i])/unit[i]<1e-6;
        }
    }
    if(q.positive || q.poisson) {
        for(size_t b=0;b<q.blocks.size();++b)if(!q.is_fixed(b))for(const auto* e:q.reporting[b]) {
            const auto f=q.event_value(*e,out.parameters);
            if(minimum_response(f.value[0],f.value[1],f.value[2],q.is_physics(b)?e->kinematics.epsilon:q.epsilon[b])<=
               2e-8*std::max({std::abs(f.value[0]),std::abs(f.value[1]),std::abs(f.value[2]),1e-30}))out.boundary=true;
        }
        const auto c=q.coefficients(out.parameters);
        for(size_t b=0;b<q.blocks.size();++b)
            if(q.is_nuisance(b))
            if(minimum_response(c[3*b],c[3*b+1],c[3*b+2],q.epsilon[b])<=
               2e-8*std::max({std::abs(c[3*b]),std::abs(c[3*b+1]),std::abs(c[3*b+2]),1e-30})) out.boundary=true;
    }
    return out;
}
inline std::vector<double> proxy_covariance(const std::vector<std::vector<double>>& jac,
                                           const std::vector<double>& covariance) {
    const size_t nc=jac.size(),np=jac.front().size();
    std::vector<double> result(nc*nc,0.);
    for(size_t a=0;a<nc;++a) for(size_t b=0;b<nc;++b)
        for(size_t i=0;i<np;++i) for(size_t j=0;j<np;++j)
            result[a*nc+b]+=jac[a][i]*covariance[i*np+j]*jac[b][j];
    return result;
}
}
#include "xsec_proxy_staged.h"
