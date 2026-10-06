// Statistical study only: reuse production kinematics and configured binning.
#include "xsec_response.h"
#include "xsec_vertex_epsilon.h"
extern "C" int pi0_mc_rows() {
    const AnalysisConfig c;return c.n_tprime*c.n_q2*c.n_xb*c.n_phi;
}
extern "C" int pi0_mc_terms(int n,const double* input,int* rows,int* blocks,double* basis) {
    try {
        const AnalysisConfig c;
        if(c.mmiss_select!="window")throw std::runtime_error("MC study requires explicit window selection");
        validate_xsec_binning(c);
        for(int i=0;i<n;++i) {
            // Reco Q2,t,tmin,xB,phi,mmiss,full_weight,sigcm;
            // matched raw Q2i,Wi,ti,phipqi,hsxptari,hsyptari.
            const double* a=input+14*i;rows[i]=blocks[i]=-1;
            const double tp=a[1]-a[2];
            if(!(a[5]>=c.mmiss_lower_gev && a[5]<=c.mmiss_upper_gev &&
                 a[0]>=c.q2_min && a[0]<=c.q2_max && a[3]>=c.xb_min && a[3]<=c.xb_max &&
                 tp>=c.tprime_min && tp<=c.tprime_max && xsec_inside_diamond(c,a[3],a[0])))continue;
            const int it=find_bin(c.tprime_bin_edges,tp),iq=find_bin(c.q2_bin_edges,a[0]);
            const int ix=find_bin(c.xb_bin_edges_by_q2.at(iq),a[3]),ip=find_bin(c.phi_bin_edges,a[4],true);
            if(it<0 || ix<0 || ip<0)continue;
            if(!std::isfinite(a[7]) || std::abs(a[7])<1e-20)throw std::runtime_error("Absent generator support");
            const double w=a[6]/a[7],q=a[8],W=a[9],x=q/(W*W-c.mp*c.mp+q);
            const double eps=nps_xsec::vertex_epsilon_from_exclusive_simc(q,W,a[12],a[13],c.hms_theta_deg,c.mp);
            const double truth_tp=-a[10]-nps_xsec::forward_t(q,W,c.mp,c.mpi0),phi=wrap_phi(a[11]);
            if(!(std::isfinite(w)&&w>=0&&x>0&&x<1))throw std::runtime_error("Invalid event basis");
            rows[i]=((it*c.n_q2+iq)*c.n_xb+ix)*c.n_phi+ip;
            blocks[i]=nps_xsec::truth_block(q,x,truth_tp,c.q2_bin_edges,c.xb_bin_edges_by_q2,c.tprime_bin_edges);
            basis[3*i]=w/(2*TMath::Pi());
            basis[3*i+1]=w*std::sqrt(2*eps*(1+eps))*std::cos(phi)/(2*TMath::Pi());
            basis[3*i+2]=w*eps*std::cos(2*phi)/(2*TMath::Pi());
        }
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
