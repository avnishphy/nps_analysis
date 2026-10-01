#include "snapshot/good_event_selection_helper.h"
#include "THaRunBase.h"
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TSystem.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <regex>
#include <vector>
#include <string>
#include <stdexcept>
// Run from this script's directory. Only reads paths frozen in manifest.tsv.
void export_cached(int onlyrun=0,int batch=0,int batches=1) {
  std::ifstream manifest("manifest.tsv"); int run,seg;std::string input;
  gSystem->mkdir("columns",true);
  const std::regex tre("g\\.(evnum|evtyp|trigbits|evtime)|T\\.hms\\.(hEDTM|hTRIG3|hTRIG4|hTRIG6)_(tdcTimeRaw|tdcTime|tdcMultiplicity)|H\\.(BCM4A\\.scalerCurrent|1MHz\\.scalerTime|EDTM\\.scaler)");
  const std::regex sre("evNumber|evcount|H\\.(1MHz\\.scalerTime|BCM4A\\.(scalerCurrent|scalerCharge)|EDTM\\.scaler|[hp](TRIG[1-6]|L1ACCP|EDTM_CP|PRE(40|100|150|200))\\.scaler|S1X\\.scalerRate)");
  while(manifest>>run>>seg>>input) {
    if((onlyrun&&run!=onlyrun)||(!onlyrun&&run%batches!=batch))continue;
    const std::string stem="columns/run"+std::to_string(run)+"_seg"+std::to_string(seg);
    if(!gSystem->AccessPathName((stem+".done").c_str())) continue;
    std::cout<<"START "<<run<<" "<<seg<<" "<<input<<std::endl;
    try {
      auto pick=effstuff::build_good_selection_summary(input,effstuff::make_default_selection_settings());
      std::ofstream sel(stem+"_selection.tsv");sel<<std::setprecision(17);
      sel<<"ok\t"<<pick.ok<<"\nmessage\t"<<pick.message<<"\nlow\t"<<pick.current_min_uA<<"\nhigh\t"<<pick.current_max_uA<<"\ni0\t"<<pick.i0_used_uA<<"\nmean_current\t"<<pick.mean_current_uA<<"\nhel_charge_before\t"<<pick.hel_charge_before_cut_uC<<"\nhel_charge_after\t"<<pick.hel_charge_after_cut_uC<<"\ngevnum_ranges\t"<<effstuff::ranges_to_string(pick.accepted_gevnum_ranges)<<"\nevcount_ranges\t"<<effstuff::ranges_to_string(pick.accepted_evcount_ranges)<<"\n";
      TFile f(input.c_str(),"READ");if(f.IsZombie())throw std::runtime_error("Cannot open input");
      if(auto*r=dynamic_cast<THaRunBase*>(f.Get("Run_Data")))r->Print();
      for(const std::string name:{"TSH","T"}) {
        auto*t=dynamic_cast<TTree*>(f.Get(name.c_str()));if(!t)throw std::runtime_error("Missing "+name);
        std::vector<std::string> cols;auto*bs=t->GetListOfBranches();
        for(int i=0;i<bs->GetEntries();++i) {
          auto*b=static_cast<TBranch*>(bs->At(i));std::string bn=b->GetName();
          if(!std::regex_match(bn,name=="T"?tre:sre))continue;
          auto*l=b->GetLeaf(bn.c_str());
          if(l&&std::string(l->GetTypeName())=="Double_t"&&!l->GetLeafCount()&&l->GetLenStatic()==1)cols.push_back(bn);
        }
        if(cols.empty())throw std::runtime_error("No columns "+name);
        std::vector<double> v(cols.size());std::vector<TBranch*> active;t->SetBranchStatus("*",0);
        std::ofstream names(stem+"_"+name+"_columns.txt");
        t->SetCacheSize(32*1024*1024);
        for(size_t j=0;j<cols.size();++j){names<<cols[j]<<'\n';t->SetBranchStatus(cols[j].c_str(),1);t->SetBranchAddress(cols[j].c_str(),&v[j]);active.push_back(t->GetBranch(cols[j].c_str()));t->AddBranchToCache(cols[j].c_str(),false);}
        t->StopCacheLearningPhase();
        std::ofstream out(stem+"_"+name+".bin",std::ios::binary);
        for(Long64_t i=0;i<t->GetEntries();++i){for(auto*b:active)if(b->GetEntry(i)<0)throw std::runtime_error("Read failure");out.write(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(double));}
        if(!out)throw std::runtime_error("Write failure");t->ResetBranchAddresses();
        std::cout<<"EXPORTED "<<run<<" "<<seg<<" "<<name<<" "<<t->GetEntries()<<" "<<cols.size()<<std::endl;
      }
      std::ofstream done(stem+".done");done<<"complete\n";
    } catch(const std::exception& e) {std::ofstream err(stem+".error");err<<e.what()<<'\n';std::cerr<<"FAILED "<<run<<" "<<seg<<" "<<e.what()<<std::endl;}
  }
}
