#include "THaRunBase.h"
#include "THaRunParameters.h"
#include <TFile.h>
#include <fstream>
#include <iostream>
void inspect_run() {
  TFile f("/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_4398_0_1_-1.root", "READ");
  auto *r = dynamic_cast<THaRunBase*>(f.Get("Run_Data"));
  if (!r) { std::cerr << "Missing Run_Data" << std::endl; return; }
  r->Print();
  if (r->GetParameters()) r->GetParameters()->Print();
  std::cout << "NCONFIG " << r->GetNConfig() << std::endl;
  std::ofstream out("/tmp/nps4398_resume/daq_config.txt");
  for (size_t i=0; i<r->GetNConfig(); ++i)
    out << "\nCONFIG " << i << "\n" << r->GetDAQConfig(i) << "\n";
}
