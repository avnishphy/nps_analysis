#include "../src/xsec_extract/simc_pi0_model.h"
#include <TFile.h>
#include <TTree.h>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    assert(argc == 2);
    const std::filesystem::path out(argv[1]);
    std::filesystem::create_directories(out);
    const double pi = std::acos(-1.0), mp = 0.9382720813, mpi = 0.1349768;
    auto model = nps_simc_pi0::at_reference(5.7, 0.55, -0.15, 10.538, mp, mpi);
    assert(std::abs(model.phi_average(0, 2*pi)/model.angular[0]-1) < 1e-14);
    // Independent midpoint quadrature verifies the analytic angular integral.
    for (int bin = 0; bin < 12; ++bin) {
        double lo=bin*pi/6, hi=(bin+1)*pi/6, integral=0;
        for (int i=0; i<10000; ++i) {
            double phi=lo+(i+0.5)*(hi-lo)/10000;
            integral += model.angular[0]+model.angular[1]*std::cos(phi)
                +model.angular[2]*std::cos(2*phi);
        }
        assert(std::abs(model.phi_average(lo,hi)/(integral/10000)-1)<1e-9);
    }
    bool rejected=false;
    try { nps_simc_pi0::at_reference(1, 0.8, -0.1, 10.538, mp, mpi); }
    catch (const std::runtime_error&) { rejected=true; }
    assert(rejected);
    // Inputs/results for comparison against the original Fortran functions.
    std::ofstream cases(out/"model_cases.txt"), cpp(out/"model_cpp.txt");
    cases << std::setprecision(17); cpp << std::setprecision(17);
    for (double q2 : {2.0, 5.7, 8.0}) for (double wsq : {4.1, 6.0, 9.0})
        for (double t : {-0.2, -0.6}) for (double th : {0.0, 0.3, 0.7})
            for (double eps : {0.2, 0.8}) for (double phi : {0.0, 1.0, 3.0}) {
                auto a=nps_simc_pi0::angular_coefficients(q2,wsq,t,th,eps);
                cases << th << ' ' << phi << ' ' << t << ' ' << q2 << ' ' << wsq << ' ' << eps << '\n';
                cpp << a[0]+a[1]*std::cos(phi)+a[2]*std::cos(2*phi) << '\n';
            }

    // Reconstructed kinematics vary with phi. Data are exactly 1.7 times MC;
    // extraction must still use one reference point for the entire slice.
    TFile sim((out/"sim.root").c_str(),"RECREATE");
    TTree st("simulation","simulation");
    float q2,t,tmin,xb,phi,mmiss=0.95,weight,sigcm,W;
    float Q2i,Wi,ti,phipqi,epsilon_i;
    int exclusive=1;
    st.Branch("Q2",&q2); st.Branch("t",&t); st.Branch("tmin",&tmin);
    st.Branch("xB",&xb); st.Branch("phi",&phi); st.Branch("mmiss",&mmiss);
    st.Branch("full_weight",&weight); st.Branch("sigcm",&sigcm);
    st.Branch("is_exclusive",&exclusive); st.Branch("W",&W);
    st.Branch("Q2i",&Q2i); st.Branch("Wi",&Wi);
    st.Branch("ti",&ti); st.Branch("phipqi",&phipqi);
    st.Branch("epsilon_i",&epsilon_i);
    TFile data((out/"data.root").c_str(),"RECREATE");
    TTree dt("physics","physics");
    double dq,dtv,dtmin,dx,dp,dm=0.95,dw,dW;
    float scale=1,charge=1000;
    int run=1;
    dt.Branch("Q2",&dq); dt.Branch("t",&dtv); dt.Branch("tmin",&dtmin);
    dt.Branch("xB",&dx); dt.Branch("phi",&dp); dt.Branch("mmiss_all",&dm);
    dt.Branch("pi0_weight",&dw); dt.Branch("scale",&scale);
    dt.Branch("charge_uC",&charge); dt.Branch("run_number",&run); dt.Branch("W",&dW);
    for(int ip=0;ip<12;++ip) for(int k=0;k<4;++k) {
        q2=5.3+0.2*k+0.005*ip; xb=0.53+0.005*k; phi=(ip+0.5)*pi/6;
        auto m=nps_simc_pi0::at_reference(q2,xb,-0.08-0.01*k,10.538,mp,mpi);
        t=m.t; tmin=m.t-m.tprime; W=m.w; weight=1+0.1*k;
        Q2i=q2; Wi=W; ti=-t; phipqi=phi;
        epsilon_i=m.epsilon;
        sigcm=m.angular[0]+m.angular[1]*std::cos(phi)+m.angular[2]*std::cos(2*phi);
        dq=q2; dtv=t; dtmin=tmin; dx=xb; dp=phi; dW=W; dw=1.7*weight;
        st.Fill(); dt.Fill();
    }
    sim.cd(); st.Write(); data.cd(); dt.Write();
    std::cout << "PASS angular integration, domain guard, and fixture generation\n";
}
