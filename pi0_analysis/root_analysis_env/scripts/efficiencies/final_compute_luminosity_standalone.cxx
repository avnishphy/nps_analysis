
// ============================================================================
// final_compute_luminosity_standalone.cxx
// Standalone improved luminosity and event selection analyzer for NPS
//
// Refactored for strict scaler delta logic, robust run parsing, and physics consistency
//
// Author: [Avnish Singh (LLM ofc!!)]
// Date: [Today's Date]
// ============================================================================



#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <limits>
#include <regex>

// For portable exit codes
#include <cstdlib>

using namespace std;
namespace fs = std::filesystem;

// --- Helper Structs ---
struct RunResult {
  int run = 0;
  int n_files = 0;
  double expected_current_uA = std::numeric_limits<double>::quiet_NaN();
  double avg_current_uA = std::numeric_limits<double>::quiet_NaN();
  double charge_uC = std::numeric_limits<double>::quiet_NaN();
  double hel0_charge = std::numeric_limits<double>::quiet_NaN();
  double pos_hel_charge = std::numeric_limits<double>::quiet_NaN();
  double neg_hel_charge = std::numeric_limits<double>::quiet_NaN();
  double total_hel_charge = std::numeric_limits<double>::quiet_NaN();
  double scaler_htrig_total = std::numeric_limits<double>::quiet_NaN();
  double scaler_edtm_total = std::numeric_limits<double>::quiet_NaN();
  double trig_accp_total = std::numeric_limits<double>::quiet_NaN();
  double scaler_edtm_no_cut = std::numeric_limits<double>::quiet_NaN();
  double trig_accp_no_cut = std::numeric_limits<double>::quiet_NaN();
  double beam_time = std::numeric_limits<double>::quiet_NaN();
  double beam_on_frac = std::numeric_limits<double>::quiet_NaN();
  long long tdc_edtm_accepted = 0;
  double NewGen_EDTM_livetime = std::numeric_limits<double>::quiet_NaN();
  std::string prescale_token;
  double ps_factor = 1.0;
  std::string which_TRIG;
  int negative_dt_intervals = 0;
  int non_monotonic_scaler_steps = 0;
  std::vector<std::string> missing_branches;
};


// --- Utility Functions ---


static int extract_prescale_r(const std::string &token_in) {
  std::string token = token_in;
  const size_t first = token.find_first_not_of(" \t\n\r");
  if (first == std::string::npos) return std::numeric_limits<int>::min();
  token.erase(0, first);
  const size_t last = token.find_last_not_of(" \t\n\r");
  token.erase(last + 1);
  const std::regex re("ps\\d*\\s*=\\s*(\\d+)", std::regex_constants::icase);
  std::smatch m;
  if (std::regex_search(token, m, re) && m.size() > 1) {
    try { return std::stoi(m[1].str()); } catch (...) { return std::numeric_limits<int>::min(); }
  }
  const size_t eq = token.find('=');
  if (eq != std::string::npos) {
    try {
      std::string rhs = token.substr(eq + 1);
      const size_t rhs_first = rhs.find_first_not_of(" \t\n\r");
      if (rhs_first == std::string::npos) return std::numeric_limits<int>::min();
      rhs.erase(0, rhs_first);
      const size_t rhs_last = rhs.find_last_not_of(" \t\n\r");
      rhs.erase(rhs_last + 1);
      return std::stoi(rhs);
    } catch (...) { return std::numeric_limits<int>::min(); }
  }
  return std::numeric_limits<int>::min();
}

static int extract_trig_number(const std::string &token_in) {
  std::string token = token_in;
  const size_t first = token.find_first_not_of(" \t\n\r");
  if (first == std::string::npos) return -1;
  token.erase(0, first);
  const size_t last = token.find_last_not_of(" \t\n\r");
  token.erase(last + 1);
  std::regex re("ps(\\d+)", std::regex_constants::icase);
  std::smatch m;
  if (std::regex_search(token, m, re) && m.size() > 1) {
    try { return std::stoi(m[1].str()); } catch (...) { return -1; }
  }
  return -1;
}

static double prescale_from_token(const std::string &token) {
  const int r = extract_prescale_r(token);
  if (r == std::numeric_limits<int>::min() || r <= 0) return 1.0;
  return std::pow(2.0, static_cast<double>(r - 1)) + 1.0;
}

// --- Main Analysis Logic ---
// [Implementation will follow the requirements and use the above helpers.]
// ...


// --- Helper: Parse prescale CSV ---
static std::map<int, std::string> parse_prescale_csv(const std::string &csv_path) {
  std::map<int, std::string> out;
  std::ifstream in(csv_path);
  if (!in) return out;
  std::string header_line;
  if (!std::getline(in, header_line)) return out;
  std::vector<std::string> headers;
  std::istringstream hss(header_line);
  std::string h;
  while (std::getline(hss, h, ',')) headers.push_back(h);
  int run_col = -1, ps_col = -1;
  for (size_t i = 0; i < headers.size(); ++i) {
    std::string lower = headers[i];
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower.find("run_number") != std::string::npos) run_col = static_cast<int>(i);
    if (lower.find("prescale") != std::string::npos) ps_col = static_cast<int>(i);
  }
  if (run_col == -1 || ps_col == -1) return out;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream lss(line);
    std::vector<std::string> fields;
    std::string f;
    while (std::getline(lss, f, ',')) fields.push_back(f);
    if ((int)fields.size() <= std::max(run_col, ps_col)) continue;
    int run = 0;
    try { run = std::stoi(fields[run_col]); } catch (...) { continue; }
    std::string token = fields[ps_col];
    out[run] = token;
  }
  return out;
}

// --- Helper: Index ROOT files by run ---
static std::map<int, std::vector<std::string>> index_run_files(const std::string &root_dir, const std::set<int> &target_runs) {
  std::map<int, std::vector<std::string>> indexed;
  if (!fs::exists(root_dir)) return indexed;
  for (const auto &entry : fs::directory_iterator(root_dir)) {
    if (!entry.is_regular_file()) continue;
    const std::string name = entry.path().filename().string();
    int run = 0;
    bool parsed = false;
    // Accept both legacy and skim naming
    if (name.find("nps_hms_coin_") == 0 && name.find(".root") != std::string::npos) {
      size_t p1 = std::string("nps_hms_coin_").size();
      size_t p2 = name.find('_', p1);
      if (p2 != std::string::npos) {
        try { run = std::stoi(name.substr(p1, p2 - p1)); parsed = true; } catch (...) { parsed = false; }
      }
    } else if (name.find("skim_run") == 0 && name.find(".root") != std::string::npos) {
      size_t p1 = std::string("skim_run").size();
      size_t p2 = name.find(".root");
      if (p2 > p1) {
        try { run = std::stoi(name.substr(p1, p2 - p1)); parsed = true; } catch (...) { parsed = false; }
      }
    }
    if (!parsed) continue;
    if (target_runs.find(run) == target_runs.end()) continue;
    indexed[run].push_back(entry.path().string());
  }
  for (auto &kv : indexed) std::sort(kv.second.begin(), kv.second.end());
  return indexed;
}

// --- Helper: Parse run list tokens ---
static std::vector<int> parse_run_list(const std::vector<std::string> &tokens) {
  std::set<int> unique;
  for (const auto &token_raw : tokens) {
    std::string token = token_raw;
    // Split comma-separated tokens
    std::replace(token.begin(), token.end(), ',', ' ');
    std::istringstream iss(token);
    std::string part;
    while (iss >> part) {
      auto dash = part.find('-');
      if (dash != std::string::npos) {
        int lo, hi;
        try {
          lo = std::stoi(part.substr(0, dash));
          hi = std::stoi(part.substr(dash + 1));
        } catch (...) { throw std::runtime_error("Invalid run range: " + part); }
        if (hi < lo) std::swap(lo, hi);
        for (int run = lo; run <= hi; ++run) unique.insert(run);
      } else {
        try { unique.insert(std::stoi(part)); }
        catch (...) { throw std::runtime_error("Invalid run token: " + part); }
      }
    }
  }
  return std::vector<int>(unique.begin(), unique.end());
}


// --- Helper: Compute median expected current from TSH branch (single pass, no repeated open) ---
static double compute_median_expected_current(const std::vector<double>& all_currents) {
  if (all_currents.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::vector<double> currents = all_currents;
  const size_t mid = currents.size() / 2;
  std::nth_element(currents.begin(), currents.begin() + mid, currents.end());
  double med_hi = currents[mid];
  if (currents.size() % 2 == 1) return med_hi;
  std::nth_element(currents.begin(), currents.begin() + (mid - 1), currents.end());
  const double med_lo = currents[mid - 1];
  return 0.5 * (med_lo + med_hi);
}

// --- Helper: Accumulate all scaler and T-tree quantities for a file in one pass ---
struct FileAccum {
  // For median expected current
  std::vector<double> scaler_currents;
  // For scaler delta-based accumulation
  double scaler_htrig_total = 0.0, scaler_edtm_total = 0.0, trig_accp_total = 0.0;
  double scaler_edtm_no_cut = 0.0, trig_accp_no_cut = 0.0;
  double charge_sum = 0.0, current_sum = 0.0, beam_time = 0.0, total_valid_time = 0.0;
  int negative_dt_intervals = 0, non_monotonic_scaler_steps = 0;
  double hel0_charge = 0.0, pos_hel_charge = 0.0, neg_hel_charge = 0.0;
  bool helicity_available = false;
  // For T-tree EDTM
  std::vector<double> edtm_tdc_values; // For histogramming
  std::set<std::string> missing_branches;
};


// --- Helper: Accumulate scaler deltas for a file (TSH tree) ---
void accumulate_scaler_deltas_for_file(const std::string& file, int trig_num, double expected_current, double window_frac, bool apply_current_window, FileAccum& accum) {
  TFile f(file.c_str(), "READ");
  if (f.IsZombie()) return;
  auto *tsh = dynamic_cast<TTree *>(f.Get("TSH"));
  if (!tsh) { accum.missing_branches.insert("TSH"); return; }
  double htrig = 0.0, edtm = 0.0, accp = 0.0, bcmCurrent = 0.0, scalerTime = 0.0, helicity = 0.0;
  bool has_bcm = false, has_hel = false, has_scalerTime = false;
  std::string htrig_name = std::string("H.hTRIG") + std::to_string(trig_num) + ".scaler";
  if (tsh->GetBranch(htrig_name.c_str())) { tsh->SetBranchStatus(htrig_name.c_str(), 1); tsh->SetBranchAddress(htrig_name.c_str(), &htrig); }
  else accum.missing_branches.insert(htrig_name);
  if (tsh->GetBranch("H.EDTM.scaler")) { tsh->SetBranchStatus("H.EDTM.scaler", 1); tsh->SetBranchAddress("H.EDTM.scaler", &edtm); }
  else accum.missing_branches.insert("H.EDTM.scaler");
  if (tsh->GetBranch("H.hL1ACCP.scaler")) { tsh->SetBranchStatus("H.hL1ACCP.scaler", 1); tsh->SetBranchAddress("H.hL1ACCP.scaler", &accp); }
  else accum.missing_branches.insert("H.hL1ACCP.scaler");
  if (tsh->GetBranch("H.BCM4A.scalerCurrent")) { has_bcm = true; tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1); tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent); }
  else accum.missing_branches.insert("H.BCM4A.scalerCurrent");
  if (tsh->GetBranch("H.1MHz.scalerTime")) { has_scalerTime = true; tsh->SetBranchStatus("H.1MHz.scalerTime", 1); tsh->SetBranchAddress("H.1MHz.scalerTime", &scalerTime); }
  else accum.missing_branches.insert("H.1MHz.scalerTime");
  if (tsh->GetBranch("helicity")) { has_hel = true; tsh->SetBranchStatus("helicity", 1); tsh->SetBranchAddress("helicity", &helicity); accum.helicity_available = true; }

  // HARD ERROR for missing critical branches
  if (!has_bcm || !has_scalerTime) {
    std::cerr << "[ERROR] File " << file << ": Missing critical branch (" << (!has_bcm ? "BCM" : "") << (!has_bcm && !has_scalerTime ? ", " : "") << (!has_scalerTime ? "scalerTime" : "") << "). Skipping file.\n";
    return;
  }

  const Long64_t n = tsh->GetEntries();
  if (n < 2) return;
  tsh->GetEntry(0);
  double prev_htrig = htrig;
  double prev_edtm = edtm;
  double prev_accp = accp;
  double prev_bcm = bcmCurrent;
  double prev_time = scalerTime;
  double prev_hel = helicity;
  // For median current, collect all scalerCurrent values
  if (std::isfinite(prev_bcm) && std::abs(prev_bcm) > 1.0) accum.scaler_currents.push_back(prev_bcm);
  for (Long64_t i = 1; i < n; ++i) {
    tsh->GetEntry(i);
    if (std::isfinite(bcmCurrent) && std::abs(bcmCurrent) > 1.0) accum.scaler_currents.push_back(bcmCurrent);
    double dt = scalerTime - prev_time;
    if (!std::isfinite(dt) || dt <= 0.0) {
      accum.negative_dt_intervals++;
      prev_htrig = htrig;
      prev_edtm = edtm;
      prev_accp = accp;
      prev_bcm = bcmCurrent;
      prev_time = scalerTime;
      prev_hel = helicity;
      continue;
    }
    double dhtrig = htrig - prev_htrig;
    double dedtm = edtm - prev_edtm;
    double daccp = accp - prev_accp;
    double ibcm = 0.5 * (bcmCurrent + prev_bcm);
    if (dhtrig < 0 || dedtm < 0 || daccp < 0) accum.non_monotonic_scaler_steps++;
    accum.total_valid_time += dt;
    // Robust current window logic
    double Imin = (1.0 - window_frac) * expected_current;
    double Imax = (1.0 + window_frac) * expected_current;
    bool in_window = apply_current_window ? (ibcm >= Imin && ibcm <= Imax) : true;
    if (in_window) {
      accum.scaler_htrig_total += std::max(0.0, dhtrig);
      accum.scaler_edtm_total += std::max(0.0, dedtm);
      accum.trig_accp_total += std::max(0.0, daccp);
      accum.charge_sum += std::max(0.0, ibcm) * dt;
      accum.beam_time += dt;
      accum.current_sum += std::max(0.0, ibcm) * dt;
      if (has_hel) {
        if (std::abs(prev_hel) < 1e-3) accum.hel0_charge += std::max(0.0, ibcm) * dt;
        else if (prev_hel > 0) accum.pos_hel_charge += std::max(0.0, ibcm) * dt;
        else if (prev_hel < 0) accum.neg_hel_charge += std::max(0.0, ibcm) * dt;
      }
    }
    accum.scaler_edtm_no_cut += std::max(0.0, dedtm);
    accum.trig_accp_no_cut += std::max(0.0, daccp);
    prev_htrig = htrig;
    prev_edtm = edtm;
    prev_accp = accp;
    prev_bcm = bcmCurrent;
    prev_time = scalerTime;
    prev_hel = helicity;
  }
}

// --- Helper: Accumulate T-tree EDTM TDC values for a file (T tree) ---
void accumulate_ttree_edtm_for_file(const std::string& file, double expected_current, double window_frac, FileAccum& accum) {
  TFile f(file.c_str(), "READ");
  if (f.IsZombie()) return;
  auto *ttree = dynamic_cast<TTree *>(f.Get("T"));
  if (!ttree) { accum.missing_branches.insert("T"); return; }
  // Use H.BCM4A.scalerCurrent for event-by-event current, matching scaler logic
  bool has_bcm = false, has_edtm_tdc = false;
  if (ttree->GetBranch("H.BCM4A.scalerCurrent")) has_bcm = true;
  else accum.missing_branches.insert("H.BCM4A.scalerCurrent");
  if (ttree->GetBranch("T.hms.hEDTM_tdcTimeRaw")) has_edtm_tdc = true;
  else accum.missing_branches.insert("T.hms.hEDTM_tdcTimeRaw");
  if (!(has_bcm && has_edtm_tdc)) return;
  double bcmCurrentT = 0.0, edtm_tdc = 0.0;
  ttree->SetBranchStatus("*", 0);
  ttree->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
  ttree->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
  ttree->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrentT);
  ttree->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtm_tdc);
  double Imin = (1.0 - window_frac) * expected_current;
  double Imax = (1.0 + window_frac) * expected_current;
  const Long64_t nT = ttree->GetEntries();
  for (Long64_t i = 0; i < nT; ++i) {
    ttree->GetEntry(i);
    if (bcmCurrentT >= Imin && bcmCurrentT <= Imax && edtm_tdc > 0.0)
      accum.edtm_tdc_values.push_back(edtm_tdc);
  }
}

// --- Helper: Find main peak and count events in ±500 window ---
long long find_edtm_peak_and_count(const std::vector<double>& tdc_values, int& main_peak_bin_out) {
  if (tdc_values.empty()) { main_peak_bin_out = 0; return 0; }
  // Use 10 ns/bin, auto-range, as in compute_luminosity_scaler.cxx
  const double bin_width = 10.0;
  std::map<int, int> tdc_hist;
  for (double val : tdc_values) {
    int bin = static_cast<int>(val / bin_width);
    tdc_hist[bin]++;
  }
  // Find main peak (bin with max count)
  int main_peak_bin = 0, max_count = 0;
  for (const auto& kv : tdc_hist) {
    if (kv.second > max_count) { max_count = kv.second; main_peak_bin = kv.first; }
  }
  double edtm_peak = (main_peak_bin + 0.5) * bin_width;
  // Count events within ±500 ns of main peak
  long long accepted = 0;
  for (double val : tdc_values) {
    if (std::abs(val - edtm_peak) < 500.0) accepted++;
  }
  main_peak_bin_out = main_peak_bin;
  return accepted;
}





// --- Main ---
int main(int argc, char **argv) {
  // --- Parse args ---
  std::string root_dir = ".";
  std::string out_csv = "final_luminosity_results.csv";
  std::vector<std::string> run_tokens;
  double window_frac = 0.15;
  double default_ps = 1.0;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--root-dir" && i + 1 < argc) root_dir = argv[++i];
    else if (arg == "--out-csv" && i + 1 < argc) out_csv = argv[++i];
    else if (arg == "--default-ps" && i + 1 < argc) default_ps = std::stod(argv[++i]);
    else if (arg == "--window-frac" && i + 1 < argc) window_frac = std::stod(argv[++i]);
    else if (arg == "--runs") {
      ++i;
      while (i < argc && std::string(argv[i]).rfind("--", 0) != 0) run_tokens.push_back(argv[i++]);
      --i;
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: " << argv[0] << " --runs <list...> [options]\n";
      return 0;
    }
  }
  if (run_tokens.empty()) {
    std::cerr << "Error: --runs is required\n";
    return 1;
  }
  std::vector<int> runs;
  try {
    runs = parse_run_list(run_tokens);
  } catch (const std::exception& e) {
    std::cerr << "Run list parse error: " << e.what() << std::endl;
    return 1;
  }
  if (runs.empty()) {
    std::cerr << "No valid runs parsed from --runs input\n";
    return 1;
  }
  const std::set<int> run_set(runs.begin(), runs.end());
  const auto indexed_files = index_run_files(root_dir, run_set);
  // --- Load prescale config CSV ---
  const std::string prescale_csv = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/config/nps_dvcs_all_kins_main.csv";
  auto prescale_lookup = parse_prescale_csv(prescale_csv);

  std::ofstream out(out_csv);
  if (!out) {
    std::cerr << "Failed to open output CSV: " << out_csv << "\n";
    return 2;
  }
  // CSV header: only real, computed fields
  out << "run,n_files,expected_current_uA,avg_current_uA,charge_uC,hel0_charge,pos_hel_charge,neg_hel_charge,total_hel_charge,scaler_htrig_total,scaler_edtm_total,trig_accp_total,scaler_edtm_no_cut,trig_accp_no_cut,beam_time,beam_on_frac,tdc_edtm_accepted,NewGen_EDTM_livetime,prescale_token,ps_factor,which_TRIG,negative_dt_intervals,non_monotonic_scaler_steps,missing_branches\n";
  out << std::fixed << std::setprecision(6);

  for (int run : runs) {
    auto it_files = indexed_files.find(run);
    if (it_files == indexed_files.end() || it_files->second.empty()) {
      std::cerr << "[WARN] No files for run " << run << ". Skipping.\n";
      continue;
    }
    const auto &files = it_files->second;
    std::string prescale_token = "";
    double ps = default_ps;
    int trig_num = -1;
    auto it_ps = prescale_lookup.find(run);
    if (it_ps != prescale_lookup.end()) {
      prescale_token = it_ps->second;
      ps = prescale_from_token(prescale_token);
      trig_num = extract_trig_number(prescale_token);
      if (trig_num < 0) {
        std::cerr << "[ERROR] Could not parse a valid trigger number from prescale token for run " << run << ". Aborting.\n";
        return 2;
      }
    } else {
      std::cerr << "[ERROR] No prescale info found for run " << run << ". Cannot determine trigger branch. Aborting.\n";
      return 2;
    }


    // --- Accumulate all file quantities in one pass (median current estimation) ---
    FileAccum accum;
    for (const auto &file : files) {
      accumulate_scaler_deltas_for_file(file, trig_num, 0.0, 0.0, true, accum);
    }
    // Improved expected current estimation
    std::vector<double> valid_currents;
    for (double val : accum.scaler_currents)
      if (std::isfinite(val) && val > 1.0) valid_currents.push_back(val);
    double expected_current = std::numeric_limits<double>::quiet_NaN();
    if (!valid_currents.empty()) {
      expected_current = compute_median_expected_current(valid_currents);
    } else {
      double sum = 0.0;
      for (double val : accum.scaler_currents)
        if (std::isfinite(val) && val > 1.0) sum += val;
      if (!accum.scaler_currents.empty())
        expected_current = sum / accum.scaler_currents.size();
    }
    if (!std::isfinite(expected_current)) {
      std::cerr << "[ERROR] Run " << run << ": Could not determine expected_current. Using 10 uA as fallback.\n";
      expected_current = 10.0;
    }

    // Robust current window handling
    bool apply_current_window = true;
    if (!std::isfinite(expected_current) || expected_current < 1.0) {
      std::cerr << "[WARN] Run " << run << ": expected_current=" << expected_current
                << " uA is invalid. Disabling current window cut for this run.\n";
      apply_current_window = false;
    }
    // (no longer needed: apply_current_window_global)

    // Now re-accumulate all physics with correct expected_current and window_frac
    FileAccum accum2;
    for (const auto &file : files) {
      accumulate_scaler_deltas_for_file(file, trig_num, expected_current, window_frac, apply_current_window, accum2);
      accumulate_ttree_edtm_for_file(file, expected_current, window_frac, accum2);
    }

    // DEBUG instrumentation
    std::cerr << "[DEBUG] Run " << run
              << ": expected_current=" << expected_current
              << " uA, scaler entries=" << accum.scaler_currents.size()
              << ", valid dt intervals=" << (accum2.total_valid_time > 0 ? accum2.total_valid_time : 0)
              << ", beam_time=" << accum2.beam_time
              << ", charge_sum=" << accum2.charge_sum
              << ", fraction in window=" << (accum2.total_valid_time > 0 ? accum2.beam_time / accum2.total_valid_time : 0)
              << "\n";

    // Fallback if beam_time == 0
    if (accum2.beam_time == 0) {
      std::cerr << "[WARN] Run " << run << ": beam_time == 0 after current window cut. Recomputing without current window.\n";
      apply_current_window = false;
      FileAccum fallback_accum;
      for (const auto &file : files)
        accumulate_scaler_deltas_for_file(file, trig_num, expected_current, window_frac, false, fallback_accum);
      accum2 = fallback_accum;
      std::cerr << "[DEBUG] Run " << run << " [fallback]: beam_time=" << accum2.beam_time << ", charge_sum=" << accum2.charge_sum << "\n";
    }

    // --- Find main peak and count events in ±500 window ---
    int main_peak_bin = 0;
    long long tdc_edtm_accepted = find_edtm_peak_and_count(accum2.edtm_tdc_values, main_peak_bin);

    // --- Physics-corrected accumulation and output ---
    RunResult rr;
    rr.run = run;
    rr.n_files = static_cast<int>(files.size());
    rr.expected_current_uA = expected_current;
    rr.ps_factor = ps;
    rr.prescale_token = prescale_token;
    char htrig_branch[64];
    snprintf(htrig_branch, sizeof(htrig_branch), "H.hTRIG%d.scaler", trig_num);
    rr.which_TRIG = std::string(htrig_branch);

    rr.charge_uC = (accum2.beam_time > 0) ? accum2.charge_sum : std::numeric_limits<double>::quiet_NaN();
    rr.hel0_charge = accum2.helicity_available ? accum2.hel0_charge : std::numeric_limits<double>::quiet_NaN();
    rr.pos_hel_charge = accum2.helicity_available ? accum2.pos_hel_charge : std::numeric_limits<double>::quiet_NaN();
    rr.neg_hel_charge = accum2.helicity_available ? accum2.neg_hel_charge : std::numeric_limits<double>::quiet_NaN();
    rr.total_hel_charge = accum2.helicity_available ? (accum2.hel0_charge + accum2.pos_hel_charge + accum2.neg_hel_charge) : std::numeric_limits<double>::quiet_NaN();
    rr.avg_current_uA = (accum2.beam_time > 0) ? (accum2.current_sum / accum2.beam_time) : std::numeric_limits<double>::quiet_NaN();
    rr.scaler_htrig_total = accum2.scaler_htrig_total;
    rr.scaler_edtm_total = accum2.scaler_edtm_total;
    rr.trig_accp_total = accum2.trig_accp_total;
    rr.scaler_edtm_no_cut = accum2.scaler_edtm_no_cut;
    rr.trig_accp_no_cut = accum2.trig_accp_no_cut;
    rr.beam_time = accum2.beam_time;
    rr.beam_on_frac = (accum2.total_valid_time > 0) ? (accum2.beam_time / accum2.total_valid_time) : std::numeric_limits<double>::quiet_NaN();
    rr.tdc_edtm_accepted = tdc_edtm_accepted;
    if (accum2.scaler_edtm_total == 0) {
      std::cerr << "[WARN] Run " << run << ": scaler_edtm_total == 0, cannot compute livetime.\n";
      rr.NewGen_EDTM_livetime = std::numeric_limits<double>::quiet_NaN();
    } else {
      rr.NewGen_EDTM_livetime = (static_cast<double>(tdc_edtm_accepted) / accum2.scaler_edtm_total * ps);
    }
    rr.negative_dt_intervals = accum2.negative_dt_intervals;
    rr.non_monotonic_scaler_steps = accum2.non_monotonic_scaler_steps;
    rr.missing_branches = std::vector<std::string>(accum2.missing_branches.begin(), accum2.missing_branches.end());

    // --- Output ---
    out << rr.run << ","
        << rr.n_files << ","
        << rr.expected_current_uA << ","
        << rr.avg_current_uA << ","
        << rr.charge_uC << ","
        << rr.hel0_charge << ","
        << rr.pos_hel_charge << ","
        << rr.neg_hel_charge << ","
        << rr.total_hel_charge << ","
        << rr.scaler_htrig_total << ","
        << rr.scaler_edtm_total << ","
        << rr.trig_accp_total << ","
        << rr.scaler_edtm_no_cut << ","
        << rr.trig_accp_no_cut << ","
        << rr.beam_time << ","
        << rr.beam_on_frac << ","
        << rr.tdc_edtm_accepted << ","
        << rr.NewGen_EDTM_livetime << ","
        << '"' << rr.prescale_token << '"' << ","
        << rr.ps_factor << ","
        << '"' << rr.which_TRIG << '"' << ","
        << rr.negative_dt_intervals << ","
        << rr.non_monotonic_scaler_steps << ","
        << '"';
    size_t j = 0;
    for (const auto &branch : rr.missing_branches) {
      if (j++) out << ";";
      std::string safe_branch = branch;
      std::replace(safe_branch.begin(), safe_branch.end(), '"', '\'');
      out << safe_branch;
    }
    out << '"' << "\n";
  }
  std::cout << "Wrote " << out_csv << "\n";
  return 0;
}
