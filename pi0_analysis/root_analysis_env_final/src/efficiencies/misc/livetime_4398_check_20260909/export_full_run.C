#include "THaRunBase.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <regex>
void export_full_run() {
  const std::string output="/tmp/nps4398_resume/full_run";
  gSystem->mkdir(output.c_str(),true);
  for(int seg=0;seg<6;++seg) {
    const auto input=std::string("/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_4398_")+std::to_string(seg)+"_1_-1.root";
    TFile f(input.c_str(),"READ");
    if(f.IsZombie()) {std::cerr<<"Bad input "<<input<<std::endl;return;}
    if(auto *r=dynamic_cast<THaRunBase*>(f.Get("Run_Data"))) r->Print();
    for(const auto name : {"TSH","TSHelH","T"}) {
      auto *t=dynamic_cast<TTree*>(f.Get(name));
      if(!t) {std::cerr<<"Missing "<<name<<std::endl;return;}
      const std::regex re("g\\.(evnum|evtyp|trigbits|evtime)|T\\.hms\\.(hEDTM|hTRIG4)_(tdcTimeRaw|tdcTime|tdcMultiplicity)|H\\.(BCM4A\\.scalerCurrent|1MHz\\.scalerTime|EDTM\\.scaler)");
      std::vector<std::string> cols;
      auto *bs=t->GetListOfBranches();
      for(int i=0;i<bs->GetEntries();++i) {
        auto*b=static_cast<TBranch*>(bs->At(i));std::string bn=b->GetName();
        if(std::string(name)=="T"&&!std::regex_match(bn,re))continue;
        auto*l=b->GetLeaf(bn.c_str());
        if(l&&std::string(l->GetTypeName())=="Double_t"&&!l->GetLeafCount()&&l->GetLenStatic()==1)cols.push_back(bn);
      }
      std::vector<double> values(cols.size());t->SetBranchStatus("*",0);
      std::string stem=output+"/seg"+std::to_string(seg)+"_"+name;
      std::ofstream names(stem+"_columns.txt");
      for(size_t j=0;j<cols.size();++j) {names<<cols[j]<<'\n';t->SetBranchStatus(cols[j].c_str(),1);t->SetBranchAddress(cols[j].c_str(),&values[j]);}
      std::ofstream out(stem+".bin",std::ios::binary);
      for(Long64_t i=0;i<t->GetEntries();++i) {
        if(t->GetEntry(i)<0) {std::cerr<<"Read failure "<<seg<<" "<<name<<" "<<i<<std::endl;return;}
        out.write(reinterpret_cast<const char*>(values.data()),values.size()*sizeof(double));
      }
      t->ResetBranchAddresses();
      std::cout<<"EXPORTED segment="<<seg<<" tree="<<name<<" rows="<<t->GetEntries()<<" columns="<<cols.size()<<std::endl;
    }
  }
}
