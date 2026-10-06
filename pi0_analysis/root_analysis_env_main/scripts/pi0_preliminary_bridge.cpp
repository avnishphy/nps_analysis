// Preliminary ensemble adapter: reuse the production constrained M0 solver.
#include "../src/xsec_extract/xsec_proxy_solver.h"
#include <iostream>
#include <streambuf>
namespace {
struct NullBuffer:std::streambuf{int overflow(int c) override{return c;}};
std::vector<nps_xsec::ModelEvent> nominal;
std::vector<int> blocks;
std::vector<int> roles; // 0=physics, 1=fitted tprime feed-in, 2=fixed model feed-in
std::vector<int> included_rows;
int physical_blocks=0;
double pivot=0;
}
extern "C" int prelim_cone(const double*h,const double*b,double epsilon,double*out){
    try{nps_xsec::staged_detail::Matrix H(3,std::vector<double>(3));
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)H[i][j]=h[3*i+j];
        auto x=nps_xsec::staged_detail::cone_minimum(H,std::vector<double>(b,b+3),epsilon);
        std::copy(x.begin(),x.end(),out);return 0;
    }catch(...){return 1;}
}
extern "C" int prelim_init(int n,const double* a,double tau0,int nblocks,
        const int* block_values,const int* block_roles,int nphysical,int nrows,const int* included){
    nominal.resize(n);pivot=tau0;
    blocks.assign(block_values,block_values+nblocks);roles.assign(block_roles,block_roles+nblocks);physical_blocks=nphysical;
    included_rows.assign(included,included+nrows);
    for(int i=0;i<n;++i){auto&e=nominal[i];auto x=a+12*i;
        e.row=int(x[0]);e.block=int(x[1]);e.weight=x[2];e.kinematics.tau=x[3];
        e.kinematics.epsilon=x[4];e.baseline.sigma_U=x[5];
        e.baseline.sigma_LT=x[6];e.baseline.sigma_TT=x[7];
        for(int k=0;k<3;++k)e.basis[k]=x[8+k];
        e.kinematics.tprime=x[11];
    }return 0;
}
extern "C" int prelim_fit(const double*y,const double*v,const int*mult,
        const double*start,double*p,double*pub,double*info,double*prediction,double*variance,double*jac){
    NullBuffer null;auto prior=std::cout.rdbuf();if(!std::getenv("PRELIM_TRACE"))std::cout.rdbuf(&null);
    try{
        using namespace nps_xsec;
        ProxyProblem q;q.positive=true;q.blocks=blocks;q.tau0=pivot;
        q.reporting.resize(blocks.size());q.epsilon.assign(blocks.size(),0);q.tprime.assign(blocks.size(),0);
        q.physical.resize(blocks.size(),false);q.fixed.resize(blocks.size(),false);
        for(size_t b=0;b<blocks.size();++b){q.physical[b]=roles[b]==0;q.fixed[b]=roles[b]==2;}
        std::vector<int> rowmap(included_rows.size(),-1);int j=0;
        for(size_t r=0;r<included_rows.size();++r)if(included_rows[r]){
            rowmap[r]=j++;q.y.push_back(y[r]);q.variance.push_back(v[r]);}
        q.events.resize(j);
        // Repeated pointers represent repeated physical observations. This
        // preserves Poisson multiplicity m (not m^2) in finite-MC Sumw2.
        for(size_t i=0;i<nominal.size();++i){auto&e=nominal[i];
            int b=std::find(blocks.begin(),blocks.end(),e.block)-blocks.begin();
            for(int m=0;m<mult[i];++m){q.reporting[b].push_back(&e);
                if(rowmap[e.row]>=0)q.events[rowmap[e.row]].push_back(&e);}
            if(mult[i])q.epsilon[b]=std::max(q.epsilon[b],e.kinematics.epsilon);
        }
        // Keep the declared nominal pivot fixed across ensembles: this is a
        // coordinate convention, so physics parameters share one definition.
        auto dv=q.variance;ProxyOptions o;o.fit_strategy="staged_feasible";o.max_calls=12000;
        if(std::getenv("PRELIM_TRACE")){std::cout<<std::unitbuf;
            if(const char* trace=std::getenv("PRELIM_HISTORY"))o.staged_trace=trace;}
        const size_t nparameters=q.parameter_count();
        std::vector<double> seed(start,start+nparameters);
        auto fit=minimize_staged_proxy(q,o,seed);int iteration=0;bool settled=false;
        while(fit.converged&&fit.status==0&&iteration<100){
            ++iteration;auto nv=dv;
            for(size_t r=0;r<nv.size();++r)nv[r]+=q.event_row(r,fit.parameters).mc_variance;
            auto pv=q.variance;q.variance=nv;
            auto next=minimize_staged_proxy(q,o,fit.parameters);double change=0;
            for(size_t k=0;k<nparameters;++k)change=std::max(change,std::abs(next.parameters[k]-fit.parameters[k])/std::max(std::abs(next.parameters[k]),1e-30));
            for(size_t r=0;r<nv.size();++r)change=std::max(change,std::abs(nv[r]-pv[r])/nv[r]);
            fit=next;if(change<1e-6){settled=true;break;}
        }
        info[0]=fit.objective;info[1]=iteration;info[2]=settled;info[3]=fit.status;
        if(!fit.converged||!settled){std::cout.rdbuf(prior);return 1;}
        std::copy(fit.parameters.begin(),fit.parameters.end(),p);
        auto c=q.coefficients(fit.parameters);
        for(int k=0;k<3;++k)for(int b=0;b<physical_blocks;++b)
            pub[k*physical_blocks+b]=c[3*b+k]*1e9;
        staged_detail::Solver solver(q,o);info[4]=solver.margin(fit.parameters,true);
        info[5]=solver.margin(fit.parameters);info[6]=q.tau0;
        for(size_t r=0;r<rowmap.size();++r)if(rowmap[r]>=0){auto z=q.event_row(rowmap[r],fit.parameters);
            prediction[r]=z.prediction;variance[r]=q.variance[rowmap[r]];
            for(size_t k=0;k<nparameters;++k)jac[r*nparameters+k]=z.jacobian[k];}
        std::cout.rdbuf(prior);return 0;
    }catch(const std::exception&e){std::cout.rdbuf(prior);std::cerr<<"preliminary bridge: "<<e.what()<<'\n';return 2;}
}
