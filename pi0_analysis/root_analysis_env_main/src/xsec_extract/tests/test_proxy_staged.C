// Exact convex nuisance block tests, including weak and nonsmooth boundaries.
#include "xsec_proxy_solver.h"
#include <cassert>
#include <random>
#include <iostream>
int main(){
    using namespace nps_xsec;using namespace nps_xsec::staged_detail;
    std::mt19937_64 rng(20261004);std::uniform_real_distribution<double> unit(-1,1);
    double worst=0;
    for(int i=0;i<2000;++i){
        double e=.1+.8*(unit(rng)+1)/2;V truth{0,unit(rng),unit(rng)};
        if(i%7==0){truth[1]=0;truth[2]=-1;}
        truth[0]=closed_u(truth[1],truth[2],e);
        if(i%4==0)truth[0]+=.3;
        if(i%4==1)truth[0]+=1e-6;
        Matrix h{{1,.03,.0001},{.03,.2,.00001},{.0001,.00001,i%3==0?1e-6:.01}};
        V b(3);for(int j=0;j<3;++j)b[j]=-2*dot(h[j],truth);
        // Known constrained optimum with a strictly positive multiplier.
        if(i%4>=2){auto c=angular_gradient(truth,0,e);for(int j=0;j<3;++j)b[j]+=.2*c[j];}
        V p=cone_minimum(h,b,e);assert(minimum_response(p[0],p[1],p[2],e)>=0);
        double err=0;for(int j=0;j<3;++j)err=std::max(err,std::abs(p[j]-truth[j]));worst=std::max(worst,err);
        if(err>2e-7){std::cerr<<"cone closure failure "<<i<<" error="<<err<<'\n';return 1;}
    }
    std::cout<<"PASS: 2000 exact nuisance-cone closures; maximum coefficient error="<<worst<<'\n';
}
