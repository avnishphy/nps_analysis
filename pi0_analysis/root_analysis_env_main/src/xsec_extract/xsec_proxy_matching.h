#pragma once
#include "xsec_vertex_epsilon.h"
#include <set>

// Model-only audit; reference failure behavior is unchanged. All counters
// cover mass-selected exclusive smeared events before rectangular bin cuts.
struct ProxyMatchAudit {
    bool enabled;std::string directory;
    long long matched=0,unmatched=0,duplicate=0,mismatch=0,invalid=0,compared=0,missing_recon=0;
    double sum=0,sum2=0,maximum=0;
    std::set<ULong64_t> seen;
    std::unique_ptr<TH1D> difference;
    std::ofstream failures;
    ProxyMatchAudit(bool on,const std::string& dir):enabled(on),directory(dir) {
        if(on){difference=std::make_unique<TH1D>("epsilon_recon_minus_vertex",";epsilon_recon - epsilon_vertex;matched events",200,-1,1);difference->SetDirectory(nullptr);
            failures.open(fs::path(dir)/"model_vertex_rejections.csv");failures<<"smeared_entry,event_id,reason\n";}
    }
    void reject(Long64_t entry,ULong64_t id,const char* reason){failures<<entry<<','<<id<<','<<reason<<'\n';}
    void compare(double recon,double vertex) {
        ++matched;
        if(!std::isfinite(recon) || recon<0 || recon>1){++missing_recon;return;}
        const double d=recon-vertex;++compared;sum+=d;sum2+=d*d;maximum=std::max(maximum,std::abs(d));difference->Fill(d);
    }
    ~ProxyMatchAudit() {
        if(!enabled)return;
        std::ofstream f(fs::path(directory)/"model_vertex_matching.csv");
        f<<std::setprecision(17)<<"matched,unmatched,duplicate,sigcm_mismatch,invalid_vertex,compared,missing_reconstructed_epsilon,mean_difference,rms_difference,max_abs_difference\n"
         <<matched<<','<<unmatched<<','<<duplicate<<','<<mismatch<<','<<invalid<<','<<compared<<','<<missing_recon<<','
         <<(compared?sum/compared:0)<<','<<(compared?std::sqrt(sum2/compared):0)<<','<<maximum<<'\n';
        std::cout<<"[VERTEX_AUDIT] matched="<<matched<<" unmatched="<<unmatched<<" duplicate="<<duplicate
                 <<" sigcm_mismatch="<<mismatch<<" invalid="<<invalid<<" compared="<<compared
                 <<" mean(recon-vertex)="<<(compared?sum/compared:0)<<" RMS="<<(compared?std::sqrt(sum2/compared):0)
                 <<" max_abs="<<maximum<<'\n';
        TFile output((fs::path(directory)/"model_vertex_epsilon.root").c_str(),"RECREATE");difference->Write();output.Close();
        if(unmatched+duplicate+mismatch+invalid)std::cerr<<"[WARN] Model response rejected raw/smeared matching or vertex failures; inspect model_vertex_matching.csv\n";
    }
};
