#include "xsec_config.h"
#include "xsec_proxy_matching.h"
#include "xsec_physics.h"
int main(int argc,char** argv) {
    try {
        if(argc!=4)throw std::runtime_error("Usage: validate_vertex_epsilon SMEARED RAW OUTDIR");
        AnalysisConfig cfg;fs::create_directories(argv[3]);
        TFile smeared(argv[1],"READ"),raw(argv[2],"READ");
        auto* sim=dynamic_cast<TTree*>(smeared.Get("simulation"));auto* vertex=dynamic_cast<TTree*>(raw.Get("h10"));
        if(!sim || !vertex)throw std::runtime_error("Missing simulation/h10 tree");
        ULong64_t id=0;int exclusive=0;float mass=0,sig=0,recon=0,q=0,w=0,x=0,y=0,rawsig=0;
        bind_branch(sim,"event_id",&id);bind_branch(sim,"is_exclusive",&exclusive);bind_branch(sim,"mmiss",&mass);
        bind_branch(sim,"sigcm",&sig);bind_branch(sim,"epsilon_i",&recon);
        bind_branch(vertex,"Q2i",&q);bind_branch(vertex,"Wi",&w);bind_branch(vertex,"hsxptari",&x);bind_branch(vertex,"hsyptari",&y);bind_branch(vertex,"sigcm",&rawsig);
        ProxyMatchAudit audit(true,argv[3]);
        for(Long64_t i=0;i<sim->GetEntries();++i) {
            sim->GetEntry(i);if(!exclusive || mass<cfg.mmiss_lower_gev || mass>cfg.mmiss_upper_gev)continue;
            if(!audit.seen.insert(id).second){++audit.duplicate;continue;}
            if(id>=static_cast<ULong64_t>(vertex->GetEntries()) || vertex->GetEntry(id)<=0){++audit.unmatched;continue;}
            if(!std::isfinite(sig) || !std::isfinite(rawsig) || std::abs(double(sig)-rawsig)>1e-6*std::max(std::abs(double(sig)),std::abs(double(rawsig)))+1e-20){++audit.mismatch;continue;}
            try{audit.compare(recon,nps_xsec::vertex_epsilon_from_exclusive_simc(q,w,x,y,cfg.hms_theta_deg,cfg.mp));}
            catch(const std::exception&){++audit.invalid;}
        }
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
