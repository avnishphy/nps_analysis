#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <regex>
void export_columns() {
  TFile f("/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_4398_0_1_-1.root", "READ");
  for (const auto name : {"TSH", "TSHelH", "T"}) {
    auto* t = dynamic_cast<TTree*>(f.Get(name));
    if (!t) continue;
    std::vector<std::string> cols;
    std::regex event_re("g\\..*|T\\.hms\\.(hEDTM|[hp]PRE[0-9]+|hTRIG[0-9]|npsTRIG[0-9])_(tdcTimeRaw|tdcTime|tdcMultiplicity)|H\\.(BCM4A|1MHz|EDTM|hTRIG4|hL1ACCP)\\.[^.]+");
    auto *branches = t->GetListOfBranches();
    for(int i=0;i<branches->GetEntries();++i) {
      auto *b=static_cast<TBranch*>(branches->At(i));
      std::string bn=b->GetName();
      if(std::string(name)=="T" && !std::regex_match(bn,event_re)) continue;
      auto *l=b->GetLeaf(bn.c_str());
      if(l && std::string(l->GetTypeName())=="Double_t" && !l->GetLeafCount() && l->GetLenStatic()==1) cols.push_back(bn);
    }
    std::vector<double> values(cols.size());
    t->SetBranchStatus("*",0);
    const std::string stem=std::string("/tmp/nps4398_resume/")+name;
    std::ofstream colout(stem+"_columns.txt");
    for(size_t j=0;j<cols.size();++j) {
      colout << cols[j] << '\n';
      t->SetBranchStatus(cols[j].c_str(),1);
      t->SetBranchAddress(cols[j].c_str(),&values[j]);
    }
    std::ofstream bin(stem+".bin",std::ios::binary);
    for(Long64_t i=0;i<t->GetEntries();++i) {
      if(t->GetEntry(i)<0) {std::cerr<<"Read failure "<<name<<" "<<i<<std::endl;return;}
      bin.write(reinterpret_cast<const char*>(values.data()),values.size()*sizeof(double));
    }
    std::cout << name << " rows=" << t->GetEntries() << " columns=" << cols.size() << std::endl;
    t->ResetBranchAddresses();
  }
}
