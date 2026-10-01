#define main nps_smearing_main_audit
#include "../src/simulation_smearing/nps_sim_smearing_new.C"
#undef main

int main() {
  TH1::AddDirectory(false);
  int mismatches=0, budget_failures=0;
  for (int which=0;which<3;++which) {
    const int n = which == 0 ? Config::MGGAMMA_NBINS
                            : which == 1 ? Config::MMISS_NBINS : Config::MPGG2_NBINS;
    const double lo = which == 0 ? Config::MGGAMMA_MIN
                                : which == 1 ? Config::MMISS_MIN : Config::MPGG2_MIN;
    const double hi = which == 0 ? Config::MGGAMMA_MAX
                                : which == 1 ? Config::MMISS_MAX : Config::MPGG2_MAX;
    TH1D h(Form("axis%d",which),"",n,lo,hi);
    FastHistogram1D f(h);
    int axis_bad=0;
    for (int i=1;i<=n+1;++i) {
      const double edge=h.GetXaxis()->GetBinLowEdge(i);
      for(double x:{std::nextafter(edge,-INFINITY),edge,std::nextafter(edge,INFINITY)}) {
        int root=h.GetXaxis()->FindFixBin(x)-1;
        if(root<0 || root>=n)root=-1;
        int fast=f.binIndex(x);
        if(root!=fast){++mismatches;++axis_bad;
          if(axis_bad<=3)std::cout<<std::setprecision(17)<<"axis="<<which<<" edge="<<i<<" x="<<x<<" root="<<root<<" fast="<<fast<<"\n";}
      }
    }
    std::cout<<"axis="<<which<<" boundary_mismatches="<<axis_bad<<"\n";
  }
  struct BudgetCase {const char* raw; size_t expected;};
  const size_t fallback=Config::FAST_PULL_CACHE_TOTAL_BUDGET_BYTES;
  const BudgetCase cases[]={{"-1",fallback},{" -1",fallback},{"\t-0",fallback},
    {" \n -12",fallback},{"0",0},{" 0",0},{"+0",0},{"4096",4096},
    {" \t4096",4096},{"",fallback},{" ",fallback},{"1x",fallback},
    {"184467440737095516160",fallback}};
  for(const auto& c:cases){setenv("NPS_FAST_PULL_CACHE_BUDGET_BYTES",c.raw,1);
    if(Config::fast_pull_cache_total_budget_bytes()!=c.expected)++budget_failures;}
  unsetenv("NPS_FAST_PULL_CACHE_BUDGET_BYTES");
  std::cout<<"boundary_mismatches_total="<<mismatches<<" budget_cases="<<sizeof(cases)/sizeof(cases[0])<<" budget_failures="<<budget_failures<<"\n";
  return mismatches || budget_failures ? 1:0;
}
