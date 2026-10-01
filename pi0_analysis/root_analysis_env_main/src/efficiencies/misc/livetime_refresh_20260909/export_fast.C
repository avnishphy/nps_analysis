// Independent scalar export pilot: direct active-branch reads, no production edits.
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
void export_fast(int run=4259,int seg=0) {
 std::string stem="columns/run"+std::to_string(run)+"_seg"+std::to_string(seg);
 std::ifstream names(stem+"_T_columns.txt");std::vector<std::string> cols;std::string s;
 while(std::getline(names,s))cols.push_back(s);
 TFile f(("/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_"+std::to_string(run)+"_"+std::to_string(seg)+"_1_-1.root").c_str(),"READ");
 auto*t=f.Get<TTree>("T");if(!t||cols.empty())throw std::runtime_error("Missing input");
 t->SetBranchStatus("*",0);std::vector<double> values(cols.size());std::vector<TBranch*>bs;
 t->SetCacheSize(32*1024*1024);
 for(size_t j=0;j<cols.size();++j){t->SetBranchStatus(cols[j].c_str(),1);t->SetBranchAddress(cols[j].c_str(),&values[j]);bs.push_back(t->GetBranch(cols[j].c_str()));t->AddBranchToCache(cols[j].c_str(),false);}
 t->StopCacheLearningPhase();gSystem->mkdir("fast_check",true);
 std::ofstream out("fast_check/run"+std::to_string(run)+"_seg"+std::to_string(seg)+"_T.bin",std::ios::binary);
 for(Long64_t i=0;i<t->GetEntries();++i){for(auto*b:bs)if(b->GetEntry(i)<0)throw std::runtime_error("Read failed");out.write((char*)values.data(),values.size()*sizeof(double));}
 if(!out)throw std::runtime_error("Write failed");std::cout<<"FAST_EXPORTED "<<t->GetEntries()<<" "<<cols.size()<<std::endl;
}
