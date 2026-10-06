#pragma once
// Opt-in Gaussian central-fit solver. Included after ProxyProblem/ProxyResult.
// The angular cone is unchanged; no confidence covariance is constructed.
#include <TDecompSVD.h>
#include <TMatrixDSymEigen.h>
#include <fstream>
#include <iomanip>

namespace nps_xsec { namespace staged_detail {
using V=std::vector<double>;
using Matrix=std::vector<V>;
inline double dot(const V&a,const V&b){return std::inner_product(a.begin(),a.end(),b.begin(),0.);}
inline V solve(const Matrix&a,const V&b){
    TMatrixD m(b.size(),b.size());TVectorD v(b.size());
    for(size_t i=0;i<b.size();++i){v[i]=b[i];for(size_t j=0;j<b.size();++j)m(i,j)=a[i][j];}
    Bool_t ok;TDecompSVD d(m);auto x=d.Solve(v,ok);
    if(!ok)throw std::runtime_error("staged_feasible singular block");
    return V(x.GetMatrixArray(),x.GetMatrixArray()+b.size());
}
inline double polynomial(const V&a,double z){double f=0;for(auto i=a.rbegin();i!=a.rend();++i)f=f*z+*i;return f;}
inline V roots(V a){
    double norm=0;for(double x:a)norm=std::max(norm,std::abs(x));if(norm==0)return {};
    for(double&x:a)x/=norm;
    while(a.size()>1 && std::abs(a.back())<2e-14)a.pop_back();
    if(a.size()==1)return {};
    if(a.size()==2){double z=-a[0]/a[1];return z>=-1&&z<=1?V{z}:V{};}
    V deriv;for(size_t i=1;i<a.size();++i)deriv.push_back(i*a[i]);
    V cuts=roots(deriv);cuts.push_back(-1);cuts.push_back(1);std::sort(cuts.begin(),cuts.end());V out;
    for(double z:cuts)if(std::abs(polynomial(a,z))<2e-13)out.push_back(z);
    for(size_t i=1;i<cuts.size();++i){
        double lo=cuts[i-1],hi=cuts[i],fl=polynomial(a,lo),fh=polynomial(a,hi);
        if(fl*fh>=0)continue;
        for(int n=0;n<60;++n){double mid=(lo+hi)/2,fm=polynomial(a,mid);if(fl*fm<=0)hi=mid;else{lo=mid;fl=fm;}}
        out.push_back((lo+hi)/2);
    }
    return out;
}
inline double closed_u(double l,double t,double e){
    double u=-minimum_response(0,l,t,e);
    for(int i=0;i<8;++i){if(minimum_response(u,l,t,e)>=0)return u;u=std::nextafter(u,INFINITY);}
    throw std::runtime_error("staged_feasible boundary representation failed");
}
inline V angular_gradient(const V&p,size_t j,double e){
    double aa=2*e*p[j+2],bb=std::sqrt(2*e*(1+e))*p[j+1];V zs{-1,1};
    if(aa>0 && std::abs(bb/(2*aa))<=1)zs.push_back(-bb/(2*aa));
    double z=zs[0];for(double x:zs)if(aa*x*x+bb*x<aa*z*z+bb*z)z=x;
    return {1,std::sqrt(2*e*(1+e))*z,e*(2*z*z-1)};
}
// Exact convex quadratic over the complete closed three-coefficient cone.
inline V cone_minimum(const Matrix&original_h,const V&original_b,double e,int* evaluations=nullptr){
    if(!(e>0 && e<1))throw std::runtime_error("staged_feasible needs 0 < epsilon < 1");
    const double scale=1/std::sqrt(original_h[0][0]),B=std::sqrt(2*e*(1+e));
    Matrix h=original_h;V b=original_b;for(size_t i=0;i<3;++i){b[i]*=scale;for(double&v:h[i])v*=scale*scale;}
    auto product=[&](const V&v){V x(3);for(size_t i=0;i<3;++i)x[i]=dot(h[i],v);return x;};
    auto value=[&](const V&v){if(evaluations)++*evaluations;return dot(b,v)+dot(v,product(v));};
    std::vector<V> candidates{V(3,0.)};V rhs=b;for(double&x:rhs)x*=-.5;V uncon=solve(h,rhs);
    // If feasible, the strictly convex unconstrained minimum is already the
    // answer. Comparing nearly equal absolute quadratics can spuriously pick
    // a boundary point in a weak direction (loss of objective significance).
    if(minimum_response(uncon[0],uncon[1],uncon[2],e)>=0){
        for(double&x:uncon)x*=scale;
        uncon[0]=std::max(uncon[0],closed_u(uncon[1],uncon[2],e));return uncon;
    }
    auto ray=[&](V w){double t=std::max(0.,-dot(b,w)/(2*dot(w,product(w))));for(double&x:w)x*=t;candidates.push_back(w);};
    for(double sign:{-1.,1.}){
        std::vector<V> w{{3*e,sign*4*e/B,1},{e,0,-1}};Matrix hh(2,V(2));V bb(2);
        for(size_t i=0;i<2;++i){bb[i]=-.5*dot(w[i],b);for(size_t j=0;j<2;++j)hh[i][j]=dot(w[i],product(w[j]));ray(w[i]);}
        V c=solve(hh,bb);if(c[0]>=0 && c[1]>=0){V x(3);for(size_t i=0;i<3;++i)x[i]=c[0]*w[0][i]+c[1]*w[1][i];candidates.push_back(x);}
    }
    Matrix w{{e,0,2*e},{0,-4*e/B,0},{1,0,0}};V aa(3,0),bb(5,0),stationary(6,0);
    for(size_t i=0;i<3;++i)for(size_t n=0;n<3;++n){aa[n]+=b[i]*w[i][n];
        for(size_t j=0;j<3;++j)for(size_t m=0;m<3;++m)bb[n+m]+=h[i][j]*w[i][n]*w[j][m];}
    // d(A^2/B)/dz=0, excluding A=0 (the apex is already a candidate).
    for(size_t i=1;i<aa.size();++i)for(size_t j=0;j<bb.size();++j)stationary[i-1+j]+=2*i*aa[i]*bb[j];
    for(size_t i=0;i<aa.size();++i)for(size_t j=1;j<bb.size();++j)stationary[i+j-1]-=aa[i]*j*bb[j];
    V zs=roots(stationary);zs.push_back(-1);zs.push_back(1);
    for(double z:zs)ray({e*(1+2*z*z),-4*e*z/B,1});
    V best=*std::min_element(candidates.begin(),candidates.end(),[&](const V&x,const V&y){return value(x)<value(y);});
    for(double&x:best)x*=scale;
    best[0]=std::max(best[0],closed_u(best[1],best[2],e));return best;
}
struct Evaluation {double f=0;V g,unit;Matrix h;};
struct Solver {
    const ProxyProblem&q;const ProxyOptions&o;size_t np,nr;int calls=0;double last_slope=INFINITY;
    struct Term{size_t row;double u,dt;};std::vector<Term> terms;
    std::vector<const ModelEvent*> physical_events;V U,DU,D2U,L,T,fixed_prediction;Matrix nuisance;
    std::vector<std::pair<size_t,double>> blocks;
    Solver(const ProxyProblem&problem,const ProxyOptions&options):q(problem),o(options),np(q.parameter_count()),nr(q.y.size()),U(nr),DU(nr),D2U(nr),L(nr),T(nr),fixed_prediction(nr),nuisance(nr,V(np,0)){
        for(size_t b=0;b<q.blocks.size();++b)if(q.is_physics(b))physical_events.insert(physical_events.end(),q.reporting[b].begin(),q.reporting[b].end());else if(q.is_nuisance(b))blocks.push_back({q.nuisance_index(b),q.epsilon[b]});
        if(!o.staged_order.empty()){
            auto order=o.staged_order;std::sort(order.begin(),order.end());
            if(order.size()!=blocks.size())throw std::runtime_error("Invalid staged block permutation");
            for(size_t j=0;j<order.size();++j)if(order[j]!=j)throw std::runtime_error("Invalid staged block permutation");
            auto original=blocks;for(size_t j=0;j<blocks.size();++j)blocks[j]=original[o.staged_order[j]];
        }
        for(size_t r=0;r<nr;++r)for(const auto*e:q.events[r]){
            size_t b=std::find(q.blocks.begin(),q.blocks.end(),e->block)-q.blocks.begin();
            if(q.is_physics(b)){terms.push_back({r,e->basis[0]*e->baseline.sigma_U,e->kinematics.tau-q.tau0});L[r]+=e->basis[1]*e->baseline.sigma_LT;T[r]+=e->basis[2]*e->baseline.sigma_TT;}
            else if(q.is_nuisance(b))for(size_t k=0;k<3;++k)nuisance[r][q.nuisance_index(b)+k]+=e->basis[k];
            else for(size_t k=0;k<3;++k)fixed_prediction[r]+=e->basis[k]*e->baseline.unseparated()[k];
        }
    }
    double margin(const V&p,bool physics_only=false)const{
        double m=INFINITY;for(const auto*e:physical_events){auto f=evaluate_cached_model(e->baseline.unseparated(),p.data(),e->kinematics.tau,q.tau0);m=std::min(m,minimum_response(f.value[0],f.value[1],f.value[2],e->kinematics.epsilon));}
        if(!physics_only)for(auto b:blocks)m=std::min(m,minimum_response(p[b.first],p[b.first+1],p[b.first+2],b.second));return m;
    }
    bool feasible(const V&p,bool physics_only=false)const{
        for(size_t j=0;j<np;++j)if(!std::isfinite(p[j]))return false;
        const auto&d=xsec_model();for(size_t j=0;j<4;++j)if(p[j]<d.lower[j]||p[j]>d.upper[j])return false;
        return margin(p,physics_only)>=0;
    }
    Evaluation evaluate(const V&p){
        if(++calls>o.max_calls)throw std::runtime_error("staged_feasible objective call limit");
        if(last_slope!=p[1]){std::fill(U.begin(),U.end(),0);std::fill(DU.begin(),DU.end(),0);std::fill(D2U.begin(),D2U.end(),0);
            for(auto t:terms){double v=t.u*std::exp(-p[1]*t.dt);U[t.row]+=v;DU[t.row]-=t.dt*v;D2U[t.row]+=t.dt*t.dt*v;}last_slope=p[1];}
        Evaluation out;out.g.assign(np,0);out.h.assign(np,V(np,0));out.unit.resize(np);V j(np);
        for(size_t r=0;r<nr;++r){j=nuisance[r];j[0]=U[r];j[1]=p[0]*DU[r];j[2]=L[r];j[3]=T[r];
            double mu=fixed_prediction[r]+p[0]*U[r]+p[2]*L[r]+p[3]*T[r]+dot(nuisance[r],p),res=mu-q.y[r],weight=1/q.variance[r];out.f+=res*res*weight;
            for(size_t i=0;i<np;++i){out.g[i]+=2*res*weight*j[i];for(size_t k=0;k<np;++k)out.h[i][k]+=2*weight*j[i]*j[k];}
        }
        for(size_t j=0;j<np;++j)out.unit[j]=1/std::sqrt(out.h[j][j]/2);return out;
    }
    bool physics(V&p){
        for(int iteration=0;iteration<100;++iteration){
            auto ev=evaluate(p);Matrix h(4,V(4));V g(4),scale(4);double maxg=0;
            for(size_t i=0;i<4;++i){scale[i]=ev.unit[i]/std::sqrt(2.);g[i]=ev.g[i]*scale[i];if(i!=1||!o.fix_u_slope)maxg=std::max(maxg,std::abs(g[i]));}
            // Absolute precision near an exact closure; residual-relative
            // precision at nonzero chi-square. This never loosens Minuit.
            if(maxg<1e-12+1e-8*std::min(1.,std::sqrt(ev.f)))return true;
            // Exact residual curvature for the pivoted U exponential.
            for(size_t r=0;r<nr;++r){double res=(fixed_prediction[r]+p[0]*U[r]+p[2]*L[r]+p[3]*T[r]+dot(nuisance[r],p)-q.y[r])/q.variance[r];ev.h[0][1]+=2*res*DU[r];ev.h[1][0]+=2*res*DU[r];ev.h[1][1]+=2*res*p[0]*D2U[r];}
            TMatrixDSym m(4);for(size_t i=0;i<4;++i)for(size_t j=0;j<4;++j)m(i,j)=ev.h[i][j]*scale[i]*scale[j];
            TMatrixDSymEigen eig(m);double minval=eig.GetEigenValues().Min(),damping=std::max(0.,1e-6-minval);
            for(size_t i=0;i<4;++i)for(size_t j=0;j<4;++j)h[i][j]=m(i,j)+(i==j?damping:0.);
            if(o.fix_u_slope){g[1]=0;for(size_t i=0;i<4;++i)h[1][i]=h[i][1]=0;h[1][1]=1;}
            for(double&v:g)v=-v;V step=solve(h,g);for(size_t i=0;i<4;++i)step[i]*=scale[i];
            double gd=0;for(size_t i=0;i<4;++i)gd+=ev.g[i]*step[i];bool accepted=false;
            for(int ls=0;ls<60;++ls){double alpha=std::ldexp(1.,-ls);V v=p;for(size_t i=0;i<4;++i)v[i]+=alpha*step[i];
                if(!feasible(v,true))continue;
                if(evaluate(v).f<=ev.f+1e-4*alpha*gd+1e-13){p=v;accepted=true;break;}}
            if(!accepted)return false;
        }return false;
    }
    bool block(V&p,size_t j,double e){
        auto ev=evaluate(p);Matrix h(3,V(3));V b(3);
        for(size_t i=0;i<3;++i){b[i]=ev.g[j+i];for(size_t k=0;k<3;++k){h[i][k]=ev.h[j+i][j+k]/2;b[i]-=2*h[i][k]*p[j+k];}}
        V v=p,best=cone_minimum(h,b,e,&calls);for(size_t i=0;i<3;++i)v[j+i]=best[i];
        if(evaluate(v).f>ev.f+1e-11)return false;p=v;return true;
    }
    double kkt(const V&p,Evaluation ev)const{
        V g=ev.g;for(size_t j=0;j<np;++j)g[j]*=ev.unit[j];if(o.fix_u_slope)g[1]=0;
        const double denom=std::max(1.,std::sqrt(dot(g,g)));
        for(auto b:blocks){size_t j=b.first;double e=b.second,scale=std::max({std::abs(p[j]),std::abs(p[j+1]),std::abs(p[j+2]),1e-30});
            if(scale==1e-30){
                // Cone apex: the gradient must belong to the dual cone.
                const double B=std::sqrt(2*e*(1+e));V poly{e*ev.g[j]+ev.g[j+2],-4*e*ev.g[j+1]/B,2*e*ev.g[j]};V z{-1,1};if(poly[2]>0 && std::abs(poly[1]/(2*poly[2]))<=1)z.push_back(-poly[1]/(2*poly[2]));
                double mn=e*ev.g[j]-ev.g[j+2];for(double x:z)mn=std::min(mn,polynomial(poly,x));
                if(mn>=0)g[j]=g[j+1]=g[j+2]=0;continue;
            }
            if(minimum_response(p[j],p[j+1],p[j+2],e)>2e-8*scale)continue;
            V zs{-1,1};double aa=2*e*p[j+2],bb=std::sqrt(2*e*(1+e))*p[j+1];
            if(aa>0 && std::abs(bb/(2*aa))<=1)zs.push_back(-bb/(2*aa));
            std::vector<V> normals;
            for(double z:zs)if(std::abs(p[j]+bb*z+e*p[j+2]*(2*z*z-1))<=2e-8*scale){
                V c{1,std::sqrt(2*e*(1+e))*z,e*(2*z*z-1)};for(size_t k=0;k<3;++k)c[k]*=ev.unit[j+k];normals.push_back(c);
            }
            V localg{g[j],g[j+1],g[j+2]},best=localg;
            for(auto c:normals){double lambda=std::max(0.,dot(localg,c)/dot(c,c));V r=localg;for(size_t k=0;k<3;++k)r[k]-=lambda*c[k];if(dot(r,r)<dot(best,best))best=r;}
            // Both endpoint normals are needed at LT=0, TT<0.
            if(normals.size()==2){auto c=normals[0],d=normals[1];double cc=dot(c,c),dd=dot(d,d),cd=dot(c,d),det=cc*dd-cd*cd;
                if(det>1e-14*cc*dd){double gc=dot(localg,c),gd=dot(localg,d),x=(gc*dd-gd*cd)/det,y=(gd*cc-gc*cd)/det;
                    if(x>=0&&y>=0){V r=localg;for(size_t k=0;k<3;++k)r[k]-=x*c[k]+y*d[k];if(dot(r,r)<dot(best,best))best=r;}}
            }
            for(size_t k=0;k<3;++k)g[j+k]=best[k];
        }
        return std::sqrt(dot(g,g))/denom;
    }
};
} // namespace staged_detail

inline ProxyResult minimize_staged_proxy(const ProxyProblem&q,const ProxyOptions&o,const std::vector<double>&seed){
    using namespace staged_detail;
    if(q.poisson || !q.positive)throw std::runtime_error("staged_feasible is a Gaussian positive-xsec diagnostic; use joint_minuit for other objectives");
    if(q.y.size()<=q.parameter_count()-size_t(o.fix_u_slope))throw std::runtime_error("Proxy fit has no positive nominal DOF");
    for(double v:q.variance)if(!(std::isfinite(v)&&v>0))throw std::runtime_error("staged_feasible requires positive finite fitted-row variances");
    Solver solver(q,o);ProxyResult result;result.initial=seed;V p=seed;
    result.errors.assign(p.size(),NAN);result.covariance.assign(p.size()*p.size(),NAN);result.covariance_status=0;result.boundary=true;
    if(!solver.feasible(p)){result.parameters=p;result.status=10;return result;}
    static int solve_id=0;const int id=solve_id++;
    std::ofstream history;
    if(!o.staged_trace.empty()){
        bool exists=std::ifstream(o.staged_trace).good();history.open(o.staged_trace,std::ios::app);history<<std::setprecision(17);
        if(!exists){history<<"solve,cycle,block,objective_before,objective_after,scaled_change,positivity_margin,block_success,KKT_relative,penalty_calls";for(auto&name:q.names())history<<','<<name;history<<'\n';}
    }
    bool physics_ok=solver.physics(p);int stable=0;double kk=INFINITY;
    for(int cycle=0;cycle<o.max_iterations;++cycle){
        auto before=solver.evaluate(p);V previous=p;bool allok=physics_ok;
        for(size_t k=0;k<=solver.blocks.size();++k){
            V prior=p;double f=solver.evaluate(p).f;bool ok;
            if(k==solver.blocks.size())ok=solver.physics(p);
            else ok=solver.block(p,solver.blocks[k].first,solver.blocks[k].second);
            auto ev=solver.evaluate(p);allok=allok&&ok;double change=0;for(size_t j=0;j<p.size();++j)change=std::max(change,std::abs(p[j]-prior[j])/ev.unit[j]);
            if(history){history<<id<<','<<cycle<<','<<(k==solver.blocks.size()?-1:int(solver.blocks[k].first))<<','<<f<<','<<ev.f<<','<<change<<','<<solver.margin(p)<<','<<ok<<','<<solver.kkt(p,ev)<<",0";for(double v:p)history<<','<<v;history<<'\n';}
            if(ev.f>f+1e-9)throw std::runtime_error("staged_feasible nonmonotone update");
        }
        auto ev=solver.evaluate(p);kk=solver.kkt(p,ev);double change=0;for(size_t j=0;j<p.size();++j)change=std::max(change,std::abs(p[j]-previous[j])/ev.unit[j]);
        const double residual_scale=std::min(1.,std::sqrt(ev.f));
        stable=std::abs(ev.f-before.f)<1e-13+1e-9*residual_scale && change<1e-10+1e-6*residual_scale && allok && kk<3e-11+3e-7*residual_scale && solver.feasible(p)?stable+1:0;
        if(stable>=2){result.converged=true;break;}
        physics_ok=true;
    }
    result.parameters=p;result.objective=q.objective(p);result.calls=solver.calls+1;result.edm=NAN;result.status=result.converged?0:11;
    if(!(result.objective<1e99) || !std::isfinite(result.objective)){result.converged=false;result.status=12;}
    std::cout<<"[STAGED_FEASIBLE] status="<<result.status<<" objective="<<result.objective<<" KKT="<<kk<<" calls="<<result.calls<<" penalty_calls="<<(result.objective>=1e99)<<" covariance=unavailable\n";
    return result;
}
} // namespace nps_xsec
