#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <cstdio>


namespace fs = std::filesystem;

/*
 * Build and run (tcsh):
 *
 *   g++ -O3 -std=c++17 -Wall -Wextra -pedantic `root-config --cflags` compute_luminosity_scaler.cxx -o compute_luminosity_scaler `root-config --glibs`
 *
 *   ./compute_luminosity_scaler --runs 4076 4077 4080-4083 \
 *     --root-dir /cache/hallc/c-nps/analysis/pass2/replays/production/ \
 *     --out-csv livetime_results_recomputed.csv
 *     --out-csv livetime_results_recomputed.csv
 */



struct RunResult {
    // NewGEN NPS coin livetime (ps6 only)
    double NewGEN_nps_coin_livetime = std::numeric_limits<double>::quiet_NaN();
    std::string NewGEN_nps_coin_livetime_str = "N.A.";
    double NewGEN_nps_coin_livetime_trig6 = std::numeric_limits<double>::quiet_NaN();
    std::string NewGEN_nps_coin_livetime_trig6_str = "N.A.";
  int run = 0;
  int n_files = 0;
  double expected_current_uA = std::numeric_limits<double>::quiet_NaN();
  double avg_current_uA = std::numeric_limits<double>::quiet_NaN(); 
  double charge_uC = 0.0;
  double scaler_htrig_total = 0.0;
  double scaler_edtm_total = 0.0;
  double trig_accp_total = 0.0;
  double scaler_edtm_no_cut = 0.0;
  double trig_accp_no_cut = 0.0;
  double beam_on_percent_edtm = std::numeric_limits<double>::quiet_NaN();
  double beam_on_percent_trig_accp = std::numeric_limits<double>::quiet_NaN();
  long long t_edtm_accepted = 0;
  long long t_htrig4_all_accepted = 0;
  long long t_htrig_phy_accepted = 0;
  std::string which_TRIG;
  double clta_livetime_tsh = std::numeric_limits<double>::quiet_NaN();
  double clta_livetime_tdc = std::numeric_limits<double>::quiet_NaN();
  double cltp_livetime_tdc = std::numeric_limits<double>::quiet_NaN();
  double tlt_livetime_edtm = std::numeric_limits<double>::quiet_NaN();
  int negative_dt_intervals = 0;
  int non_monotonic_scaler_steps = 0;
  int missing_branches = 0;
  // OG pretrig livetime 3/4
  double OG_pretrig_livetime_3_4 = std::numeric_limits<double>::quiet_NaN();
  double OG_pretrig_livetime_3_4_err = std::numeric_limits<double>::quiet_NaN();
  std::string prescale_token; // from config CSV
  double ps_factor = std::numeric_limits<double>::quiet_NaN(); // calculated
  double beam_time = 0.0; // total beam time within current cut

  // NewGen EDTM livetime
  double NewGen_EDTM_livetime = std::numeric_limits<double>::quiet_NaN();
  long long newgen_edtm_numerator = 0;
};

struct Config {
  std::string root_dir = "/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/root_analysis_env_skim/x60_4b";
  std::string out_csv = "livetime_results_recomputed.csv";
  std::vector<std::string> run_tokens;
  double default_ps = 1.0;
  double current_window = 1.5;
};

static double safe_div(double n, double d) {
  if (d == 0.0 || !std::isfinite(d)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return n / d;
}

static double clamp01(double x) {
  if (!std::isfinite(x)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

static bool starts_with(const std::string &s, const std::string &prefix) {
  return s.rfind(prefix, 0) == 0;
}

static bool ends_with(const std::string &s, const std::string &suffix) {
  if (s.size() < suffix.size()) return false;
  return s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static int extract_prescale_r(const std::string &token_in) {
  std::string token = token_in;
  // Match Python logic: token.strip()
  const size_t first = token.find_first_not_of(" \t\n\r");
  if (first == std::string::npos) {
    return std::numeric_limits<int>::min();
  }
  token.erase(0, first);
  const size_t last = token.find_last_not_of(" \t\n\r");
  token.erase(last + 1);

  const std::regex re("ps\\d*\\s*=\\s*(\\d+)", std::regex_constants::icase);
  std::smatch m;
  if (std::regex_search(token, m, re) && m.size() > 1) {
    try {
      return std::stoi(m[1].str());
    } catch (...) {
      return std::numeric_limits<int>::min();
    }
  }

  const size_t eq = token.find('=');
  if (eq != std::string::npos) {
    try {
      std::string rhs = token.substr(eq + 1);
      const size_t rhs_first = rhs.find_first_not_of(" \t\n\r");
      if (rhs_first == std::string::npos) {
        return std::numeric_limits<int>::min();
      }
      rhs.erase(0, rhs_first);
      const size_t rhs_last = rhs.find_last_not_of(" \t\n\r");
      rhs.erase(rhs_last + 1);
      return std::stoi(rhs);
    } catch (...) {
      return std::numeric_limits<int>::min();
    }
  }

  return std::numeric_limits<int>::min();
}

// Extracts the TRIG number from a prescale token, e.g., "ps4=..." -> 4
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
    try {
      return std::stoi(m[1].str());
    } catch (...) {
      return -1;
    }
  }
  return -1;
}

static double prescale_from_token(const std::string &token) {
  // Exact rule from combine_analysis_branches.py:
  // if parse fails or r <= 0 -> 1; else 2**(r-1) + 1
  const int r = extract_prescale_r(token);
  if (r == std::numeric_limits<int>::min() || r <= 0) {
    return 1.0;
  }
  return std::pow(2.0, static_cast<double>(r - 1)) + 1.0;
}

// Reads a CSV file and returns a map from run number to prescale token (as string)
static std::map<int, std::string> parse_prescale_csv(const std::string &csv_path) {
  std::map<int, std::string> out;
  std::ifstream in(csv_path);
  if (!in) {
    std::cerr << "[WARN] Could not open prescale config CSV: " << csv_path << "\n";
    return out;
  }

  std::string header_line;
  if (!std::getline(in, header_line)) {
    std::cerr << "[WARN] Empty prescale config CSV: " << csv_path << "\n";
    return out;
  }

  std::vector<std::string> headers;
  std::istringstream hss(header_line);
  std::string h;
  while (std::getline(hss, h, ',')) {
    headers.push_back(h);
  }

  int run_col = -1;
  int ps_col = -1;
  for (size_t i = 0; i < headers.size(); ++i) {
    std::string lower = headers[i];
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower.find("run_number") != std::string::npos) run_col = static_cast<int>(i);
    if (lower.find("prescale") != std::string::npos) ps_col = static_cast<int>(i);
  }
  if (run_col == -1) {
    std::cerr << "[ERROR] Config CSV missing 'run_number' column.\n";
    return out;
  }
  if (ps_col == -1) {
    std::cerr << "[ERROR] Config CSV missing 'prescale' column.\n";
    return out;
  }
  std::cerr << "[INFO] Using prescale column '" << headers[ps_col] << "' from config CSV\n";

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream lss(line);
    std::vector<std::string> fields;
    std::string f;
    while (std::getline(lss, f, ',')) fields.push_back(f);
    if ((int)fields.size() <= std::max(run_col, ps_col)) continue;
    int run = 0;
    try {
      run = std::stoi(fields[run_col]);
    } catch (...) {
      continue;
    }
    std::string token = fields[ps_col];
    out[run] = token;
  }
  return out;
}

static std::vector<int> parse_run_list(const std::vector<std::string> &tokens) {
  std::set<int> unique;
  for (const auto &token_raw : tokens) {
    std::string token = token_raw;
    const size_t first = token.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) continue;
    token.erase(0, first);
    const size_t last = token.find_last_not_of(" \t\n\r");
    token.erase(last + 1);
    if (token.empty()) continue;

    auto dash = token.find('-');
    if (dash != std::string::npos) {
      int lo = std::stoi(token.substr(0, dash));
      int hi = std::stoi(token.substr(dash + 1));
      if (hi < lo) std::swap(lo, hi);
      for (int run = lo; run <= hi; ++run) unique.insert(run);
    } else {
      unique.insert(std::stoi(token));
    }
  }
  return std::vector<int>(unique.begin(), unique.end());
}



static std::map<int, std::vector<std::string>> index_run_files(
    const std::string &root_dir,
    const std::set<int> &target_runs) {
  std::map<int, std::vector<std::string>> indexed;
  if (!fs::exists(root_dir)) return indexed;

  const std::string legacy_prefix = "nps_hms_coin_";
  const std::string legacy_suffix = "_1_-1.root";
  const std::string skim_prefix = "skim_run";
  const std::string skim_suffix = ".root";

  for (const auto &entry : fs::directory_iterator(root_dir)) {
    if (!entry.is_regular_file()) continue;
    const std::string name = entry.path().filename().string();

    int run = 0;
    bool parsed = false;

    if (starts_with(name, legacy_prefix) && ends_with(name, legacy_suffix)) {
      const size_t run_start = legacy_prefix.size();
      const size_t run_end = name.find('_', run_start);
      if (run_end != std::string::npos && run_end > run_start) {
        try {
          run = std::stoi(name.substr(run_start, run_end - run_start));
          parsed = true;
        } catch (...) {
          parsed = false;
        }
      }
    } else if (starts_with(name, skim_prefix) && ends_with(name, skim_suffix)) {
      const size_t run_start = skim_prefix.size();
      const size_t run_end = name.size() - skim_suffix.size();
      if (run_end > run_start) {
        try {
          run = std::stoi(name.substr(run_start, run_end - run_start));
          parsed = true;
        } catch (...) {
          parsed = false;
        }
      }
    }

    if (!parsed) continue;

    if (target_runs.find(run) == target_runs.end()) continue;
    indexed[run].push_back(entry.path().string());
  }

  for (auto &kv : indexed) {
    std::sort(kv.second.begin(), kv.second.end());
  }
  return indexed;
}

static bool has_branch(TTree *tree, const char *name) {
  return tree && tree->GetBranch(name) != nullptr;
}

static double compute_median_expected_current(
  const std::vector<std::string> &files) {
  std::vector<double> currents;

  for (const auto &path : files) {
    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) continue;

    auto *tsh = dynamic_cast<TTree *>(file.Get("TSH"));
    if (!tsh) continue;
    if (!has_branch(tsh, "H.BCM4A.scalerCurrent")) continue;

    double bcmCurrent = 0.0;
    tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent);
    tsh->SetBranchStatus("*", 0);
    tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1);

    const Long64_t n = tsh->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
      tsh->GetEntry(i);
      // Only keep beam-on values: finite and above threshold (e.g., 1.0 uA)
      if (std::isfinite(bcmCurrent) && std::abs(bcmCurrent) > 1.0)
        currents.push_back(bcmCurrent);
    }
  }

  if (currents.empty()) return std::numeric_limits<double>::quiet_NaN();

  const size_t mid = currents.size() / 2;
  std::nth_element(currents.begin(), currents.begin() + mid, currents.end());
  double med_hi = currents[mid];

  if (currents.size() % 2 == 1) return med_hi;

  std::nth_element(currents.begin(), currents.begin() + (mid - 1), currents.end());
  const double med_lo = currents[mid - 1];
  return 0.5 * (med_lo + med_hi);
}

static RunResult analyze_run(

  int run,
  const std::vector<std::string> &files,
  double expected_current,
  double ps_factor,
  double current_window,
  int trig_num) {

  // NewGEN NPS coin livetime counters
  long long newgen_num = 0;
  long long newgen_den = 0;
  long long newgen_trig6_num = 0;
  long long newgen_trig6_den = 0;

  RunResult r;
  r.run = run;
  r.n_files = static_cast<int>(files.size());
  r.expected_current_uA = expected_current;

  double charge_uC = 0.0;
  double current_time_weighted = 0.0;
  double time_total_cut = 0.0;

  // Build branch names from trig_num
  char htrig_branch[64];
  snprintf(htrig_branch, sizeof(htrig_branch), "H.hTRIG%d.scaler", trig_num);
  r.which_TRIG = std::string(htrig_branch);

  // OG pretrig livetime 3/4 counters
  long long og3_numerator = 0;
  long long og3_denominator = 0;

  // --- First pass: collect all edtm_tdc values in current window for peak finding ---
  std::vector<double> edtm_tdc_all;
  for (const auto &path : files) {
    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) continue;
    auto *tevt = dynamic_cast<TTree *>(file.Get("T"));
    if (!tevt) continue;
    double edtm_tdc = 0.0;
    double bcmCurrent_evt = 0.0;
    tevt->SetBranchStatus("*", 0);
    tevt->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    tevt->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
    tevt->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtm_tdc);
    tevt->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent_evt);
    const Long64_t nevt = tevt->GetEntries();
    for (Long64_t i = 0; i < nevt; ++i) {
      tevt->GetEntry(i);
      bool in_current_window = std::abs(bcmCurrent_evt - expected_current) < current_window;
      if (in_current_window && edtm_tdc > 0.0)
        edtm_tdc_all.push_back(edtm_tdc);
    }
  }

  // Find the main peak position (mode) of edtm_tdc_all
  double edtm_peak = 0.0;
  if (!edtm_tdc_all.empty()) {
    // Use histogram binning to find the mode
    const double bin_width = 10.0;
    std::map<int, int> hist;
    for (double v : edtm_tdc_all) {
      int bin = static_cast<int>(v / bin_width);
      hist[bin]++;
    }
    int max_bin = 0, max_count = 0;
    for (const auto& kv : hist) {
      if (kv.second > max_count) {
        max_count = kv.second;
        max_bin = kv.first;
      }
    }
    edtm_peak = (max_bin + 0.5) * bin_width;
  } else {
    edtm_peak = 1000.0; // fallback default
  }

  // --- Second pass: main analysis as before, but use dynamic edtm_peak ---
  for (const auto &path : files) {
    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) {
      r.missing_branches++;
      continue;
    }

    auto *tsh = dynamic_cast<TTree *>(file.Get("TSH"));
    auto *tevt = dynamic_cast<TTree *>(file.Get("T"));
    if (!tsh || !tevt) {
      r.missing_branches++;
      continue;
    }

    // Build branch list dynamically
    std::vector<std::string> scaler_branches = {
      "H.1MHz.scalerTime",
      "H.BCM4A.scalerCurrent",
      std::string(htrig_branch),
      "H.EDTM.scaler",
      "H.hL1ACCP.scaler"
    };
    bool missing_scaler = false;
    for (const auto &b : scaler_branches) {
      if (!has_branch(tsh, b.c_str())) {
        missing_scaler = true;
        break;
      }
    }
    if (missing_scaler) {
      r.missing_branches++;
      continue;
    }

    double tMHz = 0.0;
    double bcmCurrent = 0.0;
    double htrig = 0.0;
    double edtm = 0.0;
    double l1accp = 0.0;

    tsh->SetBranchAddress("H.1MHz.scalerTime", &tMHz);
    tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent);
    tsh->SetBranchAddress(htrig_branch, &htrig);
    tsh->SetBranchAddress("H.EDTM.scaler", &edtm);
    tsh->SetBranchAddress("H.hL1ACCP.scaler", &l1accp);

    // Restrict IO to required branches only for faster reading.
    tsh->SetBranchStatus("*", 0);
    tsh->SetBranchStatus("H.1MHz.scalerTime", 1);
    tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
    tsh->SetBranchStatus(htrig_branch, 1);
    tsh->SetBranchStatus("H.EDTM.scaler", 1);
    tsh->SetBranchStatus("H.hL1ACCP.scaler", 1);

    const Long64_t n = tsh->GetEntries();
    if (n < 2) continue;

    tsh->GetEntry(0);
    double prev_t = tMHz;
    double prev_htrig = htrig;
    double prev_edtm = edtm;
    double prev_l1accp = l1accp;

    for (Long64_t i = 1; i < n; ++i) {
      tsh->GetEntry(i);

      double dt = tMHz - prev_t;
      if (!std::isfinite(dt) || dt < 0.0) {
        r.negative_dt_intervals++;
        prev_t = tMHz;
        prev_htrig = htrig;
        prev_edtm = edtm;
        prev_l1accp = l1accp;
        continue;
      }

      double d_htrig = htrig - prev_htrig;
      double d_edtm = edtm - prev_edtm;
      double d_accp = l1accp - prev_l1accp;

      if (d_htrig < 0.0 || d_edtm < 0.0 || d_accp < 0.0) {
        r.non_monotonic_scaler_steps++;
      }

      d_htrig = std::max(0.0, d_htrig);
      d_edtm = std::max(0.0, d_edtm);
      d_accp = std::max(0.0, d_accp);

      r.scaler_edtm_no_cut += d_edtm;
      r.trig_accp_no_cut += d_accp;

      bool in_window = std::abs(bcmCurrent - expected_current) < current_window;
      if (in_window) {
        charge_uC += bcmCurrent * dt;
        current_time_weighted += bcmCurrent * dt;
        time_total_cut += dt;
        r.scaler_htrig_total += d_htrig;
        r.scaler_edtm_total += d_edtm;
        r.trig_accp_total += d_accp;
      }

      prev_t = tMHz;
      prev_htrig = htrig;
      prev_edtm = edtm;
      prev_l1accp = l1accp;
    }

    // Event branch selection
    char trig_tdc_branch[64];
    snprintf(trig_tdc_branch, sizeof(trig_tdc_branch), "T.hms.hTRIG%d_tdcTimeRaw", trig_num);
    const char *event_branches[] = {
      "T.hms.hEDTM_tdcTimeRaw",
      trig_tdc_branch,
      "T.hms.hTRIG3_tdcTimeRaw",
      "T.hms.hPRE40_tdcTimeRaw",
      "T.hms.npsTRIG1_tdcTimeRaw",
      "T.hms.hTRIG4_tdcTimeRaw"
    };
    bool missing_evt = false;
    for (const auto *b : event_branches) {
      if (!has_branch(tevt, b)) {
        missing_evt = true;
        break;
      }
    }
    if (missing_evt) {
      r.missing_branches++;
      continue;
    }

    // Restrict IO to required branches only for faster reading.
    tevt->SetBranchStatus("*", 0);
    tevt->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    tevt->SetBranchStatus(trig_tdc_branch, 1);
    tevt->SetBranchStatus("T.hms.hTRIG3_tdcTimeRaw", 1);
    tevt->SetBranchStatus("T.hms.hPRE40_tdcTimeRaw", 1);
    tevt->SetBranchStatus("T.hms.npsTRIG1_tdcTimeRaw", 1);
    tevt->SetBranchStatus("T.hms.npsTRIG6_tdcTimeRaw", 1);
    tevt->SetBranchStatus("T.hms.hTRIG4_tdcTimeRaw", 1);
    tevt->SetBranchStatus("H.BCM4A.scalerCurrent", 1);

    double edtm_tdc = 0.0;
    double trig_tdc = 0.0;
    double trig3_tdc = 0.0;
    double pre40_tdc = 0.0;
    double npsTRIG1_tdc = 0.0;
    double npsTRIG6_tdc = 0.0;
    double hTRIG4_tdc = 0.0;
    double bcmCurrent_evt = 0.0;
    tevt->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtm_tdc);
    tevt->SetBranchAddress(trig_tdc_branch, &trig_tdc);
    tevt->SetBranchAddress("T.hms.hTRIG3_tdcTimeRaw", &trig3_tdc);
    tevt->SetBranchAddress("T.hms.hPRE40_tdcTimeRaw", &pre40_tdc);
    tevt->SetBranchAddress("T.hms.npsTRIG1_tdcTimeRaw", &npsTRIG1_tdc);
    tevt->SetBranchAddress("T.hms.npsTRIG6_tdcTimeRaw", &npsTRIG6_tdc);
    tevt->SetBranchAddress("T.hms.hTRIG4_tdcTimeRaw", &hTRIG4_tdc);
    tevt->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent_evt);
    // NewGEN_nps_coin_livetime (ps6 only)
    // Numerator: npsTRIG1_tdc > 0 && hTRIG4_tdc > 0 && edtm_tdc == 0
    // Denominator: npsTRIG1_tdc > 0 && hTRIG4_tdc > 0
    const Long64_t nevt = tevt->GetEntries();
    int debug_prints = 0;
    for (Long64_t i = 0; i < nevt; ++i) {
      tevt->GetEntry(i);
      if (trig_num == 6 && debug_prints < 10) {
        std::cout << "[DEBUG] Event " << i << ": hTRIG4_tdc=" << hTRIG4_tdc << ", npsTRIG1_tdc=" << npsTRIG1_tdc << ", edtm_tdc=" << edtm_tdc << ", bcmCurrent_evt=" << bcmCurrent_evt << std::endl;
        debug_prints++;
      }
      bool in_current_window = std::abs(bcmCurrent_evt - expected_current) < current_window;
      // Existing NewGEN (TRIG4 && NPS TRIG1)
      if (in_current_window && hTRIG4_tdc > 0.0 && (npsTRIG1_tdc > 0.0 || npsTRIG6_tdc > 0.0)) {
        newgen_den++;
        if (edtm_tdc == 0.0) newgen_num++;
      }
      // NewGEN using TRIG6 directly (hardware)
      if (in_current_window && trig_tdc > 0.0) {
        newgen_trig6_den++;
        if (edtm_tdc == 0.0) newgen_trig6_num++;
      }

      // NewGen_EDTM_livetime numerator: count events with hEDTM_tdcTimeRaw in current window and in main peak (±500)
      if (in_current_window && edtm_tdc > 0.0 && std::abs(edtm_tdc - edtm_peak) < 500.0) {
        r.newgen_edtm_numerator++;
      }
      bool edtm_acc = (edtm_tdc != 0.0);
      bool trig_acc = (trig_tdc != 0.0);
      // Physics triggers exclude EDTM: only count triggers not coincident with EDTM
      bool trig_phy = trig_acc && !edtm_acc;

      if (edtm_acc) r.t_edtm_accepted++;
      if (trig_acc) r.t_htrig4_all_accepted++;
      if (trig_phy) r.t_htrig_phy_accepted++;

      // OG_pretrig_livetime_3/4 calculation (denominator: PRE40 > 0)
      if (pre40_tdc > 0) {
        og3_denominator++;
        if (trig3_tdc > 0 && edtm_tdc == 0) og3_numerator++;
      }
    }
  }

  // Now, after all files processed, set NewGEN_nps_coin_livetime only for ps6 (trig_num == 6), else N.A.
  if (trig_num == 6) {
    r.NewGEN_nps_coin_livetime = safe_div(static_cast<double>(newgen_num), static_cast<double>(newgen_den));
    if (newgen_den > 0) {
      char buf[64];
      snprintf(buf, sizeof(buf), "%0.6f", r.NewGEN_nps_coin_livetime);
      r.NewGEN_nps_coin_livetime_str = buf;
    } else {
      r.NewGEN_nps_coin_livetime_str = "NaN";
    }
    r.NewGEN_nps_coin_livetime_trig6 = safe_div(static_cast<double>(newgen_trig6_num), static_cast<double>(newgen_trig6_den));
    if (newgen_trig6_den > 0) {
      char buf2[64];
      snprintf(buf2, sizeof(buf2), "%0.6f", r.NewGEN_nps_coin_livetime_trig6);
      r.NewGEN_nps_coin_livetime_trig6_str = buf2;
    } else {
      r.NewGEN_nps_coin_livetime_trig6_str = "NaN";
    }
    // NewGen_EDTM_livetime: numerator = r.newgen_edtm_numerator, denominator = r.scaler_edtm_total, apply ps_factor
    r.NewGen_EDTM_livetime = safe_div(static_cast<double>(r.newgen_edtm_numerator), r.scaler_edtm_total) * ps_factor;
  } else {
    r.NewGEN_nps_coin_livetime = std::numeric_limits<double>::quiet_NaN();
    r.NewGEN_nps_coin_livetime_str = "N.A.";
    r.NewGEN_nps_coin_livetime_trig6 = std::numeric_limits<double>::quiet_NaN();
    r.NewGEN_nps_coin_livetime_trig6_str = "N.A.";
    r.NewGen_EDTM_livetime = std::numeric_limits<double>::quiet_NaN();
    r.newgen_edtm_numerator = 0;
  }

  r.charge_uC = charge_uC;
  r.avg_current_uA = safe_div(current_time_weighted, time_total_cut);
  r.beam_on_percent_edtm = clamp01(safe_div(r.scaler_edtm_total, r.scaler_edtm_no_cut));
  r.beam_on_percent_trig_accp = clamp01(safe_div(r.trig_accp_total, r.trig_accp_no_cut));

  r.clta_livetime_tsh = safe_div(r.trig_accp_total, r.scaler_htrig_total) * ps_factor;
  r.clta_livetime_tdc = safe_div(r.t_htrig4_all_accepted, r.scaler_htrig_total)
                      * ps_factor * r.beam_on_percent_trig_accp;

  // Physics livetime: exclude EDTM triggers from both numerator and denominator
  const double denom_phy = r.scaler_htrig_total - r.scaler_edtm_total;
  r.cltp_livetime_tdc = safe_div(r.t_htrig_phy_accepted, denom_phy)
                      * ps_factor * r.beam_on_percent_trig_accp;

  const double edtm_denom = safe_div(r.scaler_edtm_total, r.beam_on_percent_edtm);
  r.tlt_livetime_edtm = safe_div(r.t_edtm_accepted, edtm_denom) * ps_factor;

  // OG_pretrig_livetime_3_4 calculation and error
  r.OG_pretrig_livetime_3_4 = safe_div(static_cast<double>(og3_numerator), static_cast<double>(og3_denominator));
  if (og3_denominator > 0) {
    double p = r.OG_pretrig_livetime_3_4;
    r.OG_pretrig_livetime_3_4_err = std::sqrt(p * (1.0 - p) / og3_denominator);
  } else {
    r.OG_pretrig_livetime_3_4_err = std::numeric_limits<double>::quiet_NaN();
  }

  r.beam_time = time_total_cut;

  return r;
}

static void print_usage(const char *prog) {
  std::cerr
      << "Usage: " << prog << " --runs <list...> [options]\n"
      << "Options:\n"
      << "  --root-dir <path>            Default: /cache/hallc/c-nps/analysis/pass2/replays/production/\n"
      << "  --out-csv <path>             Default: livetime_results_recomputed.csv\n"
      << "  --default-ps <value>         Default: 1.0\n"
      << "  --current-window <uA>        Default: 1.5\n"
      // removed current-correction option
      << "  --runs <tokens...>           Required. Example: 4076 4077 4080-4083\n";
}

static bool parse_args(int argc, char **argv, Config &cfg) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--root-dir" && i + 1 < argc) {
      cfg.root_dir = argv[++i];
    } else if (arg == "--out-csv" && i + 1 < argc) {
      cfg.out_csv = argv[++i];
    } else if (arg == "--default-ps" && i + 1 < argc) {
      cfg.default_ps = std::stod(argv[++i]);
    } else if (arg == "--current-window" && i + 1 < argc) {
      cfg.current_window = std::stod(argv[++i]);
    // removed current-correction option
    } else if (arg == "--runs") {
      ++i;
      while (i < argc && std::string(argv[i]).rfind("--", 0) != 0) {
        cfg.run_tokens.push_back(argv[i]);
        ++i;
      }
      --i;
    } else if (arg == "--help" || arg == "-h") {
      print_usage(argv[0]);
      return false;
    } else {
      std::cerr << "Unknown or incomplete argument: " << arg << "\n";
      print_usage(argv[0]);
      return false;
    }
  }

  if (cfg.run_tokens.empty()) {
    std::cerr << "Error: --runs is required\n";
    print_usage(argv[0]);
    return false;
  }

  return true;
}

int main(int argc, char **argv) {
  Config cfg;
  if (!parse_args(argc, argv, cfg)) {
    return 1;
  }


  std::vector<int> runs = parse_run_list(cfg.run_tokens);
  if (runs.empty()) {
    std::cerr << "No valid runs parsed from --runs input\n";
    return 1;
  }
  const std::set<int> run_set(runs.begin(), runs.end());
  const auto indexed_files = index_run_files(cfg.root_dir, run_set);

  // Load prescale config CSV
  const std::string prescale_csv = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/config/nps_dvcs_all_kins_main.csv";
  auto prescale_lookup = parse_prescale_csv(prescale_csv);

  std::vector<RunResult> results;
  results.reserve(runs.size());
  const size_t total_runs = runs.size();
  for (size_t i = 0; i < total_runs; ++i) {
    const int run = runs[i];
    std::cout << "[" << (i + 1) << '/' << total_runs << "] run " << run << ": processing..."
              << std::flush;

    auto it_files = indexed_files.find(run);
    if (it_files == indexed_files.end()) {
      std::cout << " skipped (no matching files)\n";
      continue;
    }
    const auto &files = it_files->second;
    if (files.empty()) {
      std::cout << " skipped (empty file list)\n";
      continue;
    }

    double expected_current = compute_median_expected_current(files);
    std::string prescale_token = "";
    double ps = cfg.default_ps;
    int trig_num = -1;

    // Require prescale info for every run
    auto it_ps = prescale_lookup.find(run);
    if (it_ps == prescale_lookup.end()) {
      std::cerr << "[WARN] No prescale info for run " << run << ". Skipping this run.\n";
      continue;
    }
    prescale_token = it_ps->second;
    ps = prescale_from_token(prescale_token);
    trig_num = extract_trig_number(prescale_token);
    if (trig_num < 0) {
      std::cerr << "[WARN] Could not extract trigger number from prescale token for run " << run << ". Skipping this run.\n";
      continue;
    }
    if (!std::isfinite(expected_current)) {
      expected_current = 10.0;
    }

    RunResult rr = analyze_run(
      run,
      files,
      expected_current,
      ps,
      cfg.current_window,
      trig_num);
    rr.prescale_token = prescale_token;
    rr.ps_factor = ps;
    char htrig_branch[64];
    snprintf(htrig_branch, sizeof(htrig_branch), "H.hTRIG%d.scaler", trig_num);
    rr.which_TRIG = std::string(htrig_branch);
    results.push_back(rr);

    std::cout << " processed (" << files.size() << " files)\n";
  }

  std::ofstream out(cfg.out_csv);
  if (!out) {
    std::cerr << "Failed to open output CSV: " << cfg.out_csv << "\n";
    return 2;
  }


    out << "run,n_files,expected_current_uA,avg_current_uA,charge_uC," 
      << "scaler_htrig_total,scaler_edtm_total,trig_accp_total," 
      << "scaler_edtm_no_cut,trig_accp_no_cut,beam_on_percent_edtm," 
      << "beam_on_percent_trig_accp,t_edtm_accepted,t_htrig4_all_accepted," 
      << "t_htrig_phy_accepted,clta_livetime_tsh,clta_livetime_tdc," 
      << "cltp_livetime_tdc,tlt_livetime_edtm,NewGEN_nps_coin_livetime,NewGEN_nps_coin_livetime_trig6,NewGen_EDTM_livetime,OG_pretrig_livetime_3_4,OG_pretrig_livetime_3_4_err,beam_time,negative_dt_intervals," 
      << "non_monotonic_scaler_steps,missing_branches,prescale_token,ps_factor,which_TRIG\n";

  out << std::fixed << std::setprecision(6);

  for (const auto &r : results) {
    out << r.run << ','
      << r.n_files << ','
      << r.expected_current_uA << ','
      << r.avg_current_uA << ','
      << r.charge_uC << ','
      << r.scaler_htrig_total << ','
      << r.scaler_edtm_total << ','
      << r.trig_accp_total << ','
      << r.scaler_edtm_no_cut << ','
      << r.trig_accp_no_cut << ','
      << r.beam_on_percent_edtm << ','
      << r.beam_on_percent_trig_accp << ','
      << r.t_edtm_accepted << ','
      << r.t_htrig4_all_accepted << ','
      << r.t_htrig_phy_accepted << ','
      << r.clta_livetime_tsh << ','
      << r.clta_livetime_tdc << ','
      << r.cltp_livetime_tdc << ','
      << r.tlt_livetime_edtm << ','
      << r.NewGEN_nps_coin_livetime_str << ','
      << r.NewGEN_nps_coin_livetime_trig6_str << ','
      << r.NewGen_EDTM_livetime << ','
      << r.OG_pretrig_livetime_3_4 << ','
      << r.OG_pretrig_livetime_3_4_err << ','
      << r.beam_time << ','
      << r.negative_dt_intervals << ','
      << r.non_monotonic_scaler_steps << ','
      << r.missing_branches << ','
      << '"' << r.prescale_token << '"' << ','
      << r.ps_factor << ','
      << '"' << r.which_TRIG << '"' << '\n';
  }

  std::cout << "Wrote " << cfg.out_csv << " with " << results.size() << " runs\n";
  return 0;
}
