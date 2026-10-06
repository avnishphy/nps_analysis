// Candidate objective screening; production integration only after toy gates.
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <cstdlib>
#include <cstdio>

namespace {
struct Cell {
    std::vector<double> n,a;
    double lo=-1,hi=6,D=0,V=0;bool onoff=false;
    // Profile independent Poisson means subject to sum a_k lambda_k=mu.
    // Concave dual: max_t sum n_k log(1+t*a_k)-t*mu.
    std::array<double,2> poisson(double mu) const {
        if(onoff && a[0]==1. && a[1]<0.) {
            const double tau=-1/a[1],c=n[0],s=n[1],u=c+s-(1+tau)*mu;
            const double disc=std::sqrt(u*u+4*(1+tau)*s*mu);
            const double beta=u>=0?(u+disc)/(2*(1+tau)):(2*s*mu)/(disc-u);
            const double lc=mu+beta,ls=tau*beta;
            double val=lc+ls-c-s;
            if(c>0)val+=c*std::log(c/lc);
            if(s>0)val+=s*std::log(s/ls);
            const double deriv=lc>0?1-c/lc:1.;
            return {std::max(0.,val),-deriv};
        }
        auto score=[&](double t) {
            double g=0;
            for(size_t k=0;k<n.size();++k)if(n[k]>0) {
                const double d=1+t*a[k];
                if(d<=0)return a[k]>0?std::numeric_limits<double>::infinity():-std::numeric_limits<double>::infinity();
                g+=a[k]*n[k]/d;
            }
            return g;
        };
        double t=0,L=lo,H=hi;
        if(score(lo)<=mu)t=lo;
        else if(score(hi)>=mu)t=hi;
        else {
            t=std::clamp((D-mu)/std::max(1.,V),lo*.99,hi*.99);
            for(int it=0;it<70;++it) {
                double g=-mu,d=0;
                for(size_t k=0;k<n.size();++k)if(n[k]>0) {
                    const double w=a[k]/(1+t*a[k]);g+=n[k]*w;d-=n[k]*w*w;
                }
                if(std::abs(g)<1e-11*(1+std::abs(mu)))break;
                if(g>0)L=t;else H=t;
                double next=t-g/d;
                if(!(next>L && next<H))next=(L+H)/2;
                if(next==t)break;
                t=next;
            }
        }
        double value=-t*mu;
        for(size_t k=0;k<n.size();++k)if(n[k]>0)value+=n[k]*std::log1p(t*a[k]);
        return {std::max(0.,value),t};
    }
};
struct Problem {
    std::vector<Cell> cells;
    std::vector<double> x;
    int mode=0;double beta=0,variance_factor=0;
    std::array<double,2> profile(double turn,double width,double upper,double seed)const {
        std::vector<double> f(x.size());for(size_t i=0;i<x.size();++i)f[i]=1/(1+std::exp((x[i]-turn)/width));
        auto score=[&](double A) {
            double q=0,d=0;
            for(size_t i=0;i<x.size();++i) {
                const auto hz=cells[i].poisson(A*f[i]);const double z=hz[1];q+=f[i]*z;
                if(z>cells[i].lo && z<cells[i].hi) {
                    double v=0;for(size_t k=0;k<cells[i].n.size();++k)if(cells[i].n[k]>0) {
                        const double w=cells[i].a[k]/(1+z*cells[i].a[k]);v+=cells[i].n[k]*w*w;
                    }
                    if(v>0)d-=f[i]*f[i]/v;
                }
            }
            return std::array<double,2>{q,d};
        };
        double A=0;
        if(score(0)[0]>0) {
            double lo=0,hi=upper;A=std::clamp(seed,.001,upper*.99);
            for(int it=0;it<60;++it) {
                const auto q=score(A);
                if(std::abs(q[0])<1e-9)break;
                if(q[0]>0)lo=A;else hi=A;
                double next=A-q[0]/q[1];
                if(!(next>lo&&next<hi))next=(lo+hi)/2;
                if(std::abs(next-A)<1e-10*(1+A)){A=next;break;}A=next;
            }
        }
        double p[3]={A,turn,width};return {A,(*this)(p)};
    }
    double operator()(const double* p)const {
        double value=0;
        for(size_t i=0;i<x.size();++i) {
            const double f=1/(1+std::exp((x[i]-p[1])/p[2])),mu=p[0]*f;
            if(mode==0)value+=2*cells[i].poisson(mu)[0];
            else if(mode==1) {
                const double v=mu+variance_factor*beta;
                if(!(v>0))return 1e100;
                value+=(cells[i].D-mu)*(cells[i].D-mu)/v+std::log(v);
            } else value+=(cells[i].D-mu)*(cells[i].D-mu)/std::max(cells[i].V,1e-300);
        }
        return value;
    }
};
}

// Diagnostic objective evaluation for independent optimizer/profile checks.
extern "C" double pi0_primitive_deviance(int nk,const double* counts,
                                        const double* coefficients,const double* p) {
    Problem problem;
    const double pos=*std::max_element(coefficients,coefficients+nk),neg=*std::min_element(coefficients,coefficients+nk);
    for(int i=0;i<200;++i) {
        const double x=(i+.5)*.002;
        if(!((x>=.01&&x<=.11)||(x>=.15&&x<=.4)))continue;
        Cell c;c.lo=-1/pos;c.hi=-1/neg;c.onoff=nk==2;
        for(int k=0;k<nk;++k) {
            const double n=counts[200*k+i],a=coefficients[k];
            if(nk==2||n>0){c.n.push_back(n);c.a.push_back(a);}
            c.D+=a*n;c.V+=a*a*n;
        }
        problem.x.push_back(x);problem.cells.push_back(c);
    }
    return problem(p);
}

// Counts layout category-major, 200 fixed production mass bins. params:
// [A,turn,width,objective,zero_branch,status,evaluations,total_background].
extern "C" int pi0_fit_primitive(int nk,const double* counts,const double* coefficients,
                                 int mode,int fixed_shape,double* result) {
    if(nk<2 || mode<0 || mode>2)return 10;
    Problem problem;problem.mode=mode;
    const double pos=*std::max_element(coefficients,coefficients+nk),neg=*std::min_element(coefficients,coefficients+nk);
    if(!(pos>0 && neg<0))return 11;
    double coinmax=0,side_total=0;
    for(int i=0;i<200;++i) {
        const double x=(i+.5)*.002;
        if(!((x>=.01&&x<=.11)||(x>=.15&&x<=.4)))continue;
        Cell c;c.lo=-1/pos;c.hi=-1/neg;c.onoff=nk==2;
        for(int k=0;k<nk;++k) {
            const double n=counts[200*k+i],a=coefficients[k];
            if(!(std::isfinite(n)&&n>=0&&std::isfinite(a)))return 12;
            // Keep on/off zeros for its closed-form profile. General case
            // zero cells enter via the full-domain extrema c.lo/c.hi.
            if(nk==2 || n>0) {c.n.push_back(n);c.a.push_back(a);}
            c.D+=a*n;c.V+=a*a*n;
        }
        coinmax=std::max(coinmax,c.D);
        side_total+=counts[200+i];
        problem.x.push_back(x);problem.cells.push_back(c);
    }
    if(mode==1 && !(nk==2 && coefficients[0]==1. && coefficients[1]<0.))return 13;
    problem.beta=-coefficients[1]*side_total/problem.x.size();problem.variance_factor=1-coefficients[1];
    if(mode==2) {
        double normalization=0,normvar=0;
        for(int i=0;i<200;++i)for(int k=1;k<nk;++k) {normalization-=coefficients[k]*counts[200*k+i];normvar+=coefficients[k]*coefficients[k]*counts[200*k+i];}
        for(size_t i=0;i<problem.x.size();++i) {
            const int bin=int(problem.x[i]/.002);const double bg=counts[bin]-problem.cells[i].D;
            if(normalization>0)problem.cells[i].V+=bg*bg*normvar/(normalization*normalization);
            if(!(problem.cells[i].V>0))problem.cells[i].V=problem.cells[i].D>0?problem.cells[i].D:1.;
        }
    }
    const double upper=std::max(10.,20*std::max(1.,coinmax));
    // With fixed theta the exact Poisson profile is convex in A. Solve its
    // monotone score including the genuine constrained endpoint.
    if(fixed_shape && mode==0) {
        auto score=[&](double A) {double q=0;for(size_t i=0;i<problem.x.size();++i) {
            double f=1/(1+std::exp((problem.x[i]-.145)/.025));q+=f*problem.cells[i].poisson(A*f)[1];}return q;};
        double A=0;
        if(score(0)>0) {
            double lo=0,hi=upper;
            if(score(hi)>0)return 14;
            for(int i=0;i<60;++i) {double mid=(lo+hi)/2;if(score(mid)>0)lo=mid;else hi=mid;}A=(lo+hi)/2;
        }
        double p[3]={A,.145,.025};result[0]=A;result[1]=.145;result[2]=.025;result[3]=problem(p);result[4]=A==0;result[5]=0;result[6]=0;
        result[7]=0;for(int i=0;i<200;++i)result[7]+=A/(1+std::exp(((i+.5)*.002-.145)/.025));return 0;
    }
    double best=std::numeric_limits<double>::infinity();std::array<double,3> pars{};int status=1,calls=0;
    std::vector<std::array<double,3>> extra;
    if(mode==0&&!fixed_shape) {
        // Resolve narrow-width local modes at the unchanged mass-bin scale.
        // Profile A exactly on the grid; retain distinct local minima as seeds.
        std::vector<std::array<double,3>> grid;
        double seed=1;
        for(int j=0;j<=55;++j) {
            const double turn=.11+.002*j;const auto r=problem.profile(turn,.001,upper,seed);
            grid.push_back({r[1],r[0],turn});seed=r[0];
        }
        std::vector<std::array<double,3>> minima;
        for(size_t j=0;j<grid.size();++j)if((j==0||grid[j][0]<=grid[j-1][0])&&(j+1==grid.size()||grid[j][0]<=grid[j+1][0]))minima.push_back(grid[j]);
        std::sort(minima.begin(),minima.end());
        for(size_t j=0;j<std::min<size_t>(5,minima.size());++j)extra.push_back({minima[j][1],minima[j][2],.001});
    }
    for(int start=0;start<(fixed_shape?1:5+int(extra.size()));++start) {
        std::unique_ptr<ROOT::Math::Minimizer> minimizer(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
        minimizer->SetPrintLevel(-1);minimizer->SetMaxFunctionCalls(4000);minimizer->SetTolerance(1e-5);
        ROOT::Math::Functor fun([&](const double* p){return problem(p);},3);minimizer->SetFunction(fun);
        double amplitude=std::max(.02,problem.cells.empty()?1.:problem.cells[0].D);
        if(start>=5)amplitude=extra[start-5][0];
        minimizer->SetLimitedVariable(0,"A",std::min(amplitude,upper/2),.05,0,upper);
        if(fixed_shape) {minimizer->SetFixedVariable(1,"turn",.145);minimizer->SetFixedVariable(2,"width",.025);}
        else {
            static const double turn[]={.145,.11,.22,.145,.18},width[]={.025,.001,.1,.1,.001};
            minimizer->SetLimitedVariable(1,"turn",start<5?turn[start]:extra[start-5][1],.005,.11,.22);
            minimizer->SetLimitedVariable(2,"width",start<5?width[start]:extra[start-5][2],.002,.001,.10);
        }
        minimizer->Minimize();calls+=minimizer->NCalls();
        if(std::isfinite(minimizer->MinValue()) && minimizer->MinValue()<best) {
            best=minimizer->MinValue();std::copy(minimizer->X(),minimizer->X()+3,pars.begin());status=minimizer->Status();
        }
    }
    // Deterministic bounded pattern polish for a flat/active shape direction.
    // Verify local objective convergence rather than requiring a shape Hessian.
    // All 26 neighboring directions are searched, including coupled moves.
    if(status!=0) {
        std::array<double,3> step={.02*std::max(1.,pars[0]),fixed_shape?0:.002,fixed_shape?0:.002};
        const std::array<double,3> low={0,.11,.001},high={upper,.22,.1};
        int halvings=0;
        for(int it=0;it<2000 && halvings<22;++it) {
            auto next=pars;double value=best;
            for(int da=-1;da<=1;++da)for(int dt=-1;dt<=1;++dt)for(int dw=-1;dw<=1;++dw) {
                if((da==0&&dt==0&&dw==0)||(fixed_shape&&(dt!=0||dw!=0)))continue;
                const int dir[]={da,dt,dw};auto trial=pars;
                for(int k=0;k<3;++k)trial[k]=std::clamp(pars[k]+dir[k]*step[k],low[k],high[k]);
                const double v=problem(trial.data());++calls;
                if(v<value) {value=v;next=trial;}
            }
            if(value<best-1e-12) {best=value;pars=next;}
            else {for(auto& s:step)s*=.5;++halvings;}
        }
        if(halvings==22)status=0;
    }
    // Explicit zero is always an admissible nested model. Numerical equality
    // uses objective roundoff, not an amplitude threshold.
    double zero[3]={0,.145,.025};const double objective0=problem(zero);
    const double roundoff=64*std::numeric_limits<double>::epsilon()*problem.x.size()*std::max(1.,std::abs(objective0));
    if(objective0<=best+roundoff) {pars={0,.145,.025};best=objective0;status=0;}
    for(int i=0;i<3;++i)result[i]=pars[i];result[3]=best;result[4]=pars[0]==0;result[5]=status;result[6]=calls;
    result[7]=0;for(int i=0;i<200;++i)result[7]+=pars[0]/(1+std::exp(((i+.5)*.002-pars[1])/pars[2]));
    return status==0?0:1;
}
