#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

// compute_luminosity_scaler.cxx
//
// Purpose:
//   Recompute run-level and segment-level charge/current/livetime quantities from
//   replay ROOT files (TSH/T trees), with optional raw BCM4A calibration mode.
//
// Compile (tcsh/bash; ROOT environment must be loaded):
//   g++ -O2 -std=c++17 compute_luminosity_scaler.cxx -o compute_luminosity_scaler `root-config --cflags --libs`
//
// Required argument:
//   --runs <tokens...>
//     Run list tokens, e.g. "4077 4078 4080-4083".
//
// Optional arguments (with defaults):
//   --root-dir <path>
//     Default: /cache/hallc/c-nps/analysis/pass2/replays/production/
//     Directory containing nps_hms_coin_<run>_<seg>_1_-1.root files.
//
//   --db <path>
//     Default: /group/nps/jpcrafts/Pi_0/DataBase_production_runs_newBCMOffset.txt
//     DB file used for PS_factor fallback and Charge_yaopeng reference column.
//
//   --out-csv <path>
//     Default: livetime_results_recomputed.tsv
//     Run-level TSV output.
//
//   --out-segment-csv <path>
//     Default: segment_charge_results_recomputed.tsv
//     Segment-level TSV output.
//
//   --default-ps <value>
//     Default: 1.0
//     Used only when run PS is missing/invalid in DB.
//
//   --current-window <uA>
//     Default: 1.5
//     Beam window uses center ± width around expected run current.
//
//   --current-window-frac <f>
//     Default: 0.10
//     Beam window uses fractional half-width around expected current.
//
//   --use-absolute-window
//     Use --current-window (uA) instead of fractional window.
//
//   --current-correction <uA>
//     Default: 0.0
//     Additive correction applied to the selected current definition.
//
//   --beam-on-threshold <uA>
//     Default: 2.5
//     Defines charge_above_threshold_uC using current > threshold.
//
//   --use-raw-bcm4a
//     Default: disabled
//     Use raw BCM4A scaler branch with current reconstructed as:
//       I = (raw_rate - offset)/gain,
//     where raw_rate is computed interval-by-interval as dScaler/dt.
//
//   --use-rolling-stability
//     Default: disabled
//     Apply rolling stability window on TSHelH scaler charge/read.
//
//   --stability-window <N>
//     Default: 30
//     Rolling window size for charge stability filter.
//
//   --stability-frac-range <f>
//     Default: 0.08
//     Maximum fractional charge range (max-min)/mean in rolling window.
//
//   --bcm4a-gain <value>
//     Default: 9597.0
//     Gain used only with --use-raw-bcm4a.
//
//   --bcm4a-offset <value>
//     Default: 177.3
//     Offset used only with --use-raw-bcm4a.
//
//   --no-progress
//     Disable progress printing.
//
// Quick examples:
//   1) Minimal (uses all defaults):
//      ./compute_luminosity_scaler --runs 4077 4078 4079
//
//   2) Raw BCM4A calibrated mode (Christine constants):
//      ./compute_luminosity_scaler --runs 4077 --use-raw-bcm4a --bcm4a-gain 9597 --bcm4a-offset 177.3
//
//   3) Full custom output names and clean logs:
//      ./compute_luminosity_scaler --runs 4077-4087 --out-csv run.tsv --out-segment-csv seg.tsv --no-progress
namespace fs = std::filesystem;

struct DbRow {
  double ps = std::numeric_limits<double>::quiet_NaN();
  double ave_current = std::numeric_limits<double>::quiet_NaN();
  double charge_tot = std::numeric_limits<double>::quiet_NaN();
  int which_trigger = -1;
  std::string prescale_token;
};

struct RunResult {
  int run = 0;
  int n_files = 0;
  double expected_current_uA = std::numeric_limits<double>::quiet_NaN();
  double avg_current_uA = std::numeric_limits<double>::quiet_NaN();
  double Charge_yaopeng = std::numeric_limits<double>::quiet_NaN();
  double hel_charge_hm_uC = 0.0;
  double hel_charge_hp_uC = 0.0;
  double hel_charge_total_uC = 0.0;
  double hel0_charge_uC = 0.0;
  double charge_tshelh_uC = 0.0;
  double charge_uC = 0.0;
  double charge_above_threshold_uC = 0.0;
  double charge_above_3p0_uC = 0.0;
  double scaler_htrig4_total = 0.0;
  double scaler_edtm_total = 0.0;
  double trig_accp_total = 0.0;
  double scaler_edtm_no_cut = 0.0;
  double trig_accp_no_cut = 0.0;
  double beam_on_percent_edtm = std::numeric_limits<double>::quiet_NaN();
  double beam_on_percent_trig_accp = std::numeric_limits<double>::quiet_NaN();
  long long t_edtm_accepted = 0;
  long long t_htrig4_all_accepted = 0;
  long long t_htrig4_phy_accepted = 0;
  double clta_livetime_tsh = std::numeric_limits<double>::quiet_NaN();
  double clta_livetime_tdc = std::numeric_limits<double>::quiet_NaN();
  double cltp_livetime_tdc = std::numeric_limits<double>::quiet_NaN();
  double tlt_livetime_edtm = std::numeric_limits<double>::quiet_NaN();
  double tlt_livetime_edtm_sym_beamcut = std::numeric_limits<double>::quiet_NaN();
  double tlt_livetime_edtm_nocut_den = std::numeric_limits<double>::quiet_NaN();
  double tlt_livetime_edtm_ratio_sym = std::numeric_limits<double>::quiet_NaN();
  double tlt_livetime_edtm_ratio_nocut = std::numeric_limits<double>::quiet_NaN();
  int negative_dt_intervals = 0;
  int non_monotonic_scaler_steps = 0;
  int suspicious_scaler_jump_steps = 0;
  int missing_branches = 0;
  int hel_missing_branches = 0;
  double NewGEN_nps_coin_livetime = std::numeric_limits<double>::quiet_NaN();
  double NewGEN_nps_coin_livetime_trig6 = std::numeric_limits<double>::quiet_NaN();
  double NewGEN_electronic_livetime = std::numeric_limits<double>::quiet_NaN();
  double NewGEN_electronic_livetime_evt = std::numeric_limits<double>::quiet_NaN();
  double NewGEN_electronic_livetime_evt_noedtmden = std::numeric_limits<double>::quiet_NaN();
  long long HMS_pid_den_events = 0;
  double HMS_pid_eff = std::numeric_limits<double>::quiet_NaN();
  double HMS_cer_eff = std::numeric_limits<double>::quiet_NaN();
  double HMS_cal_eff = std::numeric_limits<double>::quiet_NaN();
  double HMS_cer_eff_tag_cal = std::numeric_limits<double>::quiet_NaN();
  double HMS_cal_eff_tag_cer = std::numeric_limits<double>::quiet_NaN();
  long long HMS_tracking_should_events = 0;
  long long HMS_tracking_did_events = 0;
  double HMS_tracking_eff = std::numeric_limits<double>::quiet_NaN();
  double HMS_tracking_eff_err = std::numeric_limits<double>::quiet_NaN();
  long long HMS_hodo_3of4_den_events = 0;
  double HMS_hodo_3of4_eff = std::numeric_limits<double>::quiet_NaN();
  long long HMS_hodo_3of4_nhits_den_events = 0;
  double HMS_hodo_3of4_nhits_eff = std::numeric_limits<double>::quiet_NaN();
  double OG_pretrig_livetime_3_4 = std::numeric_limits<double>::quiet_NaN();
  double OG_pretrig_livetime_3_4_err = std::numeric_limits<double>::quiet_NaN();
  double beam_time = 0.0;
  std::string prescale_token;
  double ps_factor = std::numeric_limits<double>::quiet_NaN();
  std::string which_TRIG;
};

struct SegmentResult {
  int run = 0;
  int segment = -1;
  std::string file;
  double expected_current_uA = std::numeric_limits<double>::quiet_NaN();
  double avg_current_uA = std::numeric_limits<double>::quiet_NaN();
  double charge_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_charge_hm_uC = 0.0;
  double hel_charge_hp_uC = 0.0;
  double hel_charge_total_uC = 0.0;
  double hel0_charge_uC = 0.0;
  double charge_tshelh_uC = 0.0;
  std::string good_evnum_ranges;
  std::string hel0_evnum_ranges;
  int negative_dt_intervals = 0;
  int non_monotonic_scaler_steps = 0;
  int missing_branches = 0;
  int hel_missing_branches = 0;
};

struct Config {
  std::string root_dir = "/cache/hallc/c-nps/analysis/pass2/replays/production/";
  std::string db_path = "/group/nps/jpcrafts/Pi_0/DataBase_production_runs_newBCMOffset.txt";
  std::string trigger_overrides_path = "/group/nps/jpcrafts/Pi_0/run_trigger_prescale_overrides.tsv";
  std::string report_dir = "/work/hallc/nps/nps-ana/REPORT_OUTPUT_pass2/COIN_old";
  bool use_report_efficiencies = false;
  std::string out_csv = "livetime_results_recomputed.tsv";
  std::string out_segment_csv = "segment_charge_results_recomputed.tsv";
  std::vector<std::string> run_tokens;
  double default_ps = 1.0;
  double current_window = 1.5; // uA
  double current_window_frac = 0.10; // fractional half-width
  bool use_fractional_window = true;
  double current_correction = 0.0;
  double beam_on_threshold_uA = 2.5;
  bool use_raw_bcm4a = false;
  double bcm4a_gain = 9597.0;
  double bcm4a_offset = 177.3;
  bool use_rolling_stability = false;
  int stability_window = 30;
  double stability_frac_range = 0.08;
  double pid_cer_npe_min = 1.0;
  double pid_cal_etotnorm_min = 0.6;
  double pid_cal_etotnorm_max = std::numeric_limits<double>::infinity();
  double hodo_track_max_chisq = 10.0;
  double hodo_track_slop_cm = 2.0;
  bool pid_require_track = true;
  bool progress = true;
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

static void compute_window_limits(double expected_current,
                                  double current_window,
                                  double current_window_frac,
                                  bool use_fractional_window,
                                  double &window_lo,
                                  double &window_hi) {
  if (use_fractional_window && std::isfinite(expected_current) && expected_current > 0.0) {
    window_lo = (1.0 - current_window_frac) * expected_current;
    window_hi = (1.0 + current_window_frac) * expected_current;
    return;
  }
  if (std::isfinite(expected_current)) {
    window_lo = expected_current - current_window;
    window_hi = expected_current + current_window;
    return;
  }
  window_lo = std::numeric_limits<double>::quiet_NaN();
  window_hi = std::numeric_limits<double>::quiet_NaN();
}

static std::vector<int> parse_run_list(const std::vector<std::string> &tokens) {
  std::set<int> unique;
  for (const auto &token_raw : tokens) {
    std::string token = token_raw;
    token.erase(0, token.find_first_not_of(" \t\n\r"));
    token.erase(token.find_last_not_of(" \t\n\r") + 1);
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

static std::map<int, DbRow> parse_db(const std::string &db_path) {
  std::map<int, DbRow> out;
  std::ifstream in(db_path);
  if (!in) return out;

  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    std::vector<std::string> parts;
    std::string p;
    while (iss >> p) parts.push_back(p);
    if (parts.size() < 11) continue;
    if (parts[0].empty() || !std::isdigit(parts[0][0])) continue;

    try {
      int run = std::stoi(parts[0]);
      DbRow row;
      row.which_trigger = std::stoi(parts[2]);
      row.ps = std::stod(parts[3]);
      row.ave_current = std::stod(parts[8]);
      row.charge_tot = std::stod(parts[10]);
      if (std::isfinite(row.ps) && row.ps > 0.0 && row.which_trigger > 0) {
        const long long ps_i = std::llround(row.ps);
        row.prescale_token = "ps" + std::to_string(row.which_trigger) + "=" + std::to_string(ps_i);
      }
      out[run] = row;
    } catch (...) {
      continue;
    }
  }
  return out;
}

static std::map<int, DbRow> parse_trigger_overrides(const std::string &path) {
  std::map<int, DbRow> out;
  if (path.empty()) return out;

  std::ifstream in(path);
  if (!in) return out;

  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    std::vector<std::string> parts;
    std::string p;
    while (iss >> p) parts.push_back(p);
    if (parts.size() < 3) continue;
    if (parts[0].empty() || !std::isdigit(parts[0][0])) continue;

    try {
      int run = std::stoi(parts[0]);
      DbRow row;
      row.which_trigger = std::stoi(parts[1]);
      row.ps = std::stod(parts[2]);
      if (parts.size() >= 4) {
        row.prescale_token = parts[3];
      }
      if (row.prescale_token.empty() && std::isfinite(row.ps) && row.ps > 0.0 && row.which_trigger > 0) {
        const long long ps_i = std::llround(row.ps);
        row.prescale_token = "ps" + std::to_string(row.which_trigger) + "=" + std::to_string(ps_i);
      }
      out[run] = row;
    } catch (...) {
      continue;
    }
  }

  return out;
}

static std::vector<std::string> find_run_files(const std::string &root_dir, int run) {
  std::vector<std::string> files;
  if (!fs::exists(root_dir)) return files;

  const std::string prefix = "nps_hms_coin_" + std::to_string(run) + "_";
  const std::string suffix = "_1_-1.root";

  for (const auto &entry : fs::directory_iterator(root_dir)) {
    if (!entry.is_regular_file()) continue;
    const std::string name = entry.path().filename().string();
    if (starts_with(name, prefix) && ends_with(name, suffix)) {
      files.push_back(entry.path().string());
    }
  }

  std::sort(files.begin(), files.end());
  return files;
}

static bool has_branch(TTree *tree, const char *name) {
  return tree && tree->GetBranch(name) != nullptr;
}

static int parse_segment_from_filename(const std::string &path, int run) {
  const std::string name = fs::path(path).filename().string();
  const std::string prefix = "nps_hms_coin_" + std::to_string(run) + "_";
  const std::string suffix = "_1_-1.root";
  if (!starts_with(name, prefix) || !ends_with(name, suffix)) {
    return -1;
  }
  const size_t lo = prefix.size();
  const size_t hi = name.size() - suffix.size();
  const std::string seg_str = name.substr(lo, hi - lo);
  try {
    return std::stoi(seg_str);
  } catch (...) {
    return -1;
  }
}

static std::vector<int> stable_evcounts_from_charge_window(
    const std::vector<int> &evcount_arr,
    const std::vector<double> &charge_arr,
    int window_n,
    double frac_range_max) {
  std::vector<int> out;
  if (evcount_arr.size() != charge_arr.size()) return out;
  if (evcount_arr.size() < static_cast<size_t>(window_n)) return out;

  size_t run_start = 0;
  for (size_t i = 1; i <= evcount_arr.size(); ++i) {
    const bool end_run = (i == evcount_arr.size()) || (evcount_arr[i] != evcount_arr[i - 1] + 1);
    if (!end_run) continue;

    const size_t run_end = i - 1;
    if (run_end + 1 - run_start >= static_cast<size_t>(window_n)) {
      double sum = 0.0;
      std::deque<size_t> dq_min;
      std::deque<size_t> dq_max;

      for (size_t j = run_start; j <= run_end; ++j) {
        sum += charge_arr[j];

        while (!dq_min.empty() && charge_arr[dq_min.back()] >= charge_arr[j]) dq_min.pop_back();
        dq_min.push_back(j);

        while (!dq_max.empty() && charge_arr[dq_max.back()] <= charge_arr[j]) dq_max.pop_back();
        dq_max.push_back(j);

        const long long old = static_cast<long long>(j) - window_n;
        if (old >= static_cast<long long>(run_start)) {
          sum -= charge_arr[static_cast<size_t>(old)];
          if (!dq_min.empty() && dq_min.front() == static_cast<size_t>(old)) dq_min.pop_front();
          if (!dq_max.empty() && dq_max.front() == static_cast<size_t>(old)) dq_max.pop_front();
        }

        if (j >= run_start + static_cast<size_t>(window_n - 1)) {
          const double mean = sum / static_cast<double>(window_n);
          if (mean > 0.0 && !dq_min.empty() && !dq_max.empty()) {
            const double qmin = charge_arr[dq_min.front()];
            const double qmax = charge_arr[dq_max.front()];
            const double frac_range = (qmax - qmin) / mean;
            if (frac_range <= frac_range_max) {
              out.push_back(evcount_arr[j]);
            }
          }
        }
      }
    }

    run_start = i;
  }

  return out;
}

static std::string ranges_from_evcount_evnum(const std::vector<int> &evcount,
                                             const std::vector<long long> &evnum) {
  if (evcount.empty() || evnum.empty() || evcount.size() != evnum.size()) return "";

  std::ostringstream oss;
  size_t start_idx = 0;
  bool first = true;

  for (size_t i = 1; i <= evcount.size(); ++i) {
    const bool end_run = (i == evcount.size()) || (evcount[i] != evcount[i - 1] + 1);
    if (!end_run) continue;

    const long long start_evnum = evnum[start_idx];
    const long long end_evnum = evnum[i - 1];
    if (!first) oss << ";";
    oss << start_evnum << "-" << end_evnum;
    first = false;
    start_idx = i;
  }

  return oss.str();
}

static std::vector<std::pair<long long, long long>> parse_evnum_ranges(const std::string &ranges) {
  std::vector<std::pair<long long, long long>> out;
  if (ranges.empty()) return out;

  std::stringstream ss(ranges);
  std::string token;
  while (std::getline(ss, token, ';')) {
    if (token.empty()) continue;
    const size_t dash = token.find('-');
    if (dash == std::string::npos) continue;
    try {
      long long lo = std::stoll(token.substr(0, dash));
      long long hi = std::stoll(token.substr(dash + 1));
      if (lo > hi) std::swap(lo, hi);
      out.emplace_back(lo, hi);
    } catch (...) {
      continue;
    }
  }

  std::sort(out.begin(), out.end());
  return out;
}

static bool evnum_in_ranges(long long evnum,
                            const std::vector<std::pair<long long, long long>> &ranges) {
  if (ranges.empty()) return true;
  auto it = std::upper_bound(
      ranges.begin(), ranges.end(), evnum,
      [](long long v, const std::pair<long long, long long> &r) { return v < r.first; });
  if (it == ranges.begin()) return false;
  --it;
  return evnum >= it->first && evnum <= it->second;
}

struct HelicitySelectionResult {
  double charge_hm_uC = 0.0;
  double charge_hp_uC = 0.0;
  double charge_total_uC = 0.0;
  double hel0_charge_uC = 0.0;
  std::string good_evnum_ranges;
  std::string hel0_evnum_ranges;
  int missing_branches = 0;
};

struct ReportSegmentMetrics {
  bool have_hdide = false;
  bool have_hscinshoulde = false;
  bool have_hodo3of4 = false;
  bool have_htrig = false;
  double hdide = 0.0;
  double hscinshoulde = 0.0;
  double hodo3of4 = std::numeric_limits<double>::quiet_NaN();
  double htrig = 0.0;
};

static ReportSegmentMetrics parse_report_segment_metrics(const std::string &report_path) {
  ReportSegmentMetrics m;
  std::ifstream in(report_path);
  if (!in) return m;

  const std::regex colon_number(R"(:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?))");
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("hdide", 0) == 0) {
      std::smatch match;
      if (std::regex_search(line, match, colon_number)) {
        m.hdide = std::stod(match[1].str());
        m.have_hdide = true;
      }
      continue;
    }
    if (line.rfind("hscinshoulde", 0) == 0) {
      std::smatch match;
      if (std::regex_search(line, match, colon_number)) {
        m.hscinshoulde = std::stod(match[1].str());
        m.have_hscinshoulde = true;
      }
      continue;
    }
    if (line.find("3_of_4 EFF") != std::string::npos) {
      std::smatch match;
      if (std::regex_search(line, match, colon_number)) {
        m.hodo3of4 = std::stod(match[1].str());
        m.have_hodo3of4 = true;
      }
      continue;
    }
    if (line.rfind("htrig", 0) == 0) {
      std::smatch match;
      if (std::regex_search(line, match, colon_number)) {
        m.htrig = std::stod(match[1].str());
        m.have_htrig = true;
      }
      continue;
    }
  }

  return m;
}

static HelicitySelectionResult analyze_segment_helicity(
    const std::string &path,
    double expected_current,
    const Config &cfg) {
  HelicitySelectionResult out;

  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) {
    out.missing_branches++;
    return out;
  }

  auto *tshelh = dynamic_cast<TTree *>(file.Get("TSHelH"));
  if (!tshelh) {
    out.missing_branches++;
    return out;
  }

  const char *hel_branches[] = {
    "H.BCM4A_Hel.scalerCurrent",
    "H.BCM4A_Hel.scalerCharge",
    "actualHelicity",
    "evcount"
  };
  for (const auto *b : hel_branches) {
    if (!has_branch(tshelh, b)) {
      out.missing_branches++;
      return out;
    }
  }

  tshelh->SetBranchStatus("*", 0);
  tshelh->SetBranchStatus("H.BCM4A_Hel.scalerCurrent", 1);
  tshelh->SetBranchStatus("H.BCM4A_Hel.scalerCharge", 1);
  tshelh->SetBranchStatus("actualHelicity", 1);
  tshelh->SetBranchStatus("evcount", 1);

  double hel_current = 0.0;
  double hel_charge = 0.0;
  double hel_state = 0.0;
  double evcount_raw = 0.0;

  tshelh->SetBranchAddress("H.BCM4A_Hel.scalerCurrent", &hel_current);
  tshelh->SetBranchAddress("H.BCM4A_Hel.scalerCharge", &hel_charge);
  tshelh->SetBranchAddress("actualHelicity", &hel_state);
  tshelh->SetBranchAddress("evcount", &evcount_raw);

  std::vector<long long> t_gevnum_sorted;
  {
    auto *t = dynamic_cast<TTree *>(file.Get("T"));
    if (!t || !has_branch(t, "g.evnum")) {
      out.missing_branches++;
    } else {
      t->SetBranchStatus("*", 0);
      t->SetBranchStatus("g.evnum", 1);
      double gevnum_raw = 0.0;
      t->SetBranchAddress("g.evnum", &gevnum_raw);
      const Long64_t nt = t->GetEntries();
      t_gevnum_sorted.reserve(static_cast<size_t>(nt));
      for (Long64_t i = 0; i < nt; ++i) {
        t->GetEntry(i);
        t_gevnum_sorted.push_back(static_cast<long long>(std::llround(gevnum_raw)));
      }
      std::sort(t_gevnum_sorted.begin(), t_gevnum_sorted.end());
    }
  }

  const Long64_t n = tshelh->GetEntries();
  if (n < 1) {
    return out;
  }

  double window_lo = std::numeric_limits<double>::quiet_NaN();
  double window_hi = std::numeric_limits<double>::quiet_NaN();
  compute_window_limits(expected_current,
                        cfg.current_window,
                        cfg.current_window_frac,
                        cfg.use_fractional_window,
                        window_lo,
                        window_hi);

  std::vector<int> evcount_arr;
  std::vector<double> charge_arr;
  std::vector<int> helicity_arr;

  std::vector<int> evcount_h0;
  std::vector<double> charge_h0;

  evcount_arr.reserve(static_cast<size_t>(n));
  charge_arr.reserve(static_cast<size_t>(n));
  helicity_arr.reserve(static_cast<size_t>(n));

  for (Long64_t i = 0; i < n; ++i) {
    tshelh->GetEntry(i);
    if (!std::isfinite(hel_current)) continue;

    const double calibrated_current = hel_current;
    const double corrected_current = calibrated_current + cfg.current_correction;
    (void)corrected_current;
    const bool in_window = (calibrated_current >= window_lo && calibrated_current <= window_hi);
    if (!in_window) continue;

    const int hel = static_cast<int>(std::llround(hel_state));
    const int evcount = static_cast<int>(std::llround(evcount_raw));
    if (hel == 0) {
      evcount_h0.push_back(evcount);
      charge_h0.push_back(hel_charge);
    } else {
      evcount_arr.push_back(evcount);
      charge_arr.push_back(hel_charge);
      helicity_arr.push_back(hel);
    }
  }

  if (evcount_arr.empty() && evcount_h0.empty()) {
    return out;
  }

  if (!evcount_arr.empty()) {
    std::vector<size_t> idx(evcount_arr.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
      return evcount_arr[a] < evcount_arr[b];
    });

    std::vector<int> evcount_sorted;
    std::vector<double> charge_sorted;
    std::vector<int> helicity_sorted;
    evcount_sorted.reserve(evcount_arr.size());
    charge_sorted.reserve(charge_arr.size());
    helicity_sorted.reserve(helicity_arr.size());

    for (size_t k : idx) {
      evcount_sorted.push_back(evcount_arr[k]);
      charge_sorted.push_back(charge_arr[k]);
      helicity_sorted.push_back(helicity_arr[k]);
    }

    evcount_arr.swap(evcount_sorted);
    charge_arr.swap(charge_sorted);
    helicity_arr.swap(helicity_sorted);
  }

  if (!evcount_h0.empty()) {
    std::vector<size_t> idx0(evcount_h0.size());
    std::iota(idx0.begin(), idx0.end(), 0);
    std::sort(idx0.begin(), idx0.end(), [&](size_t a, size_t b) {
      return evcount_h0[a] < evcount_h0[b];
    });

    std::vector<int> evcount_sorted;
    std::vector<double> charge_sorted;
    evcount_sorted.reserve(evcount_h0.size());
    charge_sorted.reserve(charge_h0.size());

    for (size_t k : idx0) {
      evcount_sorted.push_back(evcount_h0[k]);
      charge_sorted.push_back(charge_h0[k]);
    }

    evcount_h0.swap(evcount_sorted);
    charge_h0.swap(charge_sorted);
  }

  std::unordered_set<int> stable_set;
  std::unordered_set<int> stable_set_h0;
  if (cfg.use_rolling_stability) {
    const auto stable = stable_evcounts_from_charge_window(
        evcount_arr,
        charge_arr,
        cfg.stability_window,
        cfg.stability_frac_range);
    stable_set.insert(stable.begin(), stable.end());

    const auto stable_h0 = stable_evcounts_from_charge_window(
        evcount_h0,
        charge_h0,
        cfg.stability_window,
        cfg.stability_frac_range);
    stable_set_h0.insert(stable_h0.begin(), stable_h0.end());
  }

  std::vector<int> good_evcount;
  std::vector<long long> good_evnum;
  good_evcount.reserve(evcount_arr.size());
  good_evnum.reserve(evcount_arr.size());
  for (size_t i = 0; i < evcount_arr.size(); ++i) {
    if (cfg.use_rolling_stability) {
      if (stable_set.find(evcount_arr[i]) == stable_set.end()) continue;
    }
    const int hel = helicity_arr[i];
    const double q = charge_arr[i];
    if (hel < 0) out.charge_hm_uC += q;
    if (hel > 0) out.charge_hp_uC += q;
    out.charge_total_uC += q;
    const int evcount = evcount_arr[i];
    if (evcount >= 1 && static_cast<size_t>(evcount) <= t_gevnum_sorted.size()) {
      good_evcount.push_back(evcount);
      good_evnum.push_back(t_gevnum_sorted[static_cast<size_t>(evcount - 1)]);
    }
  }

  out.good_evnum_ranges = ranges_from_evcount_evnum(good_evcount, good_evnum);

  std::vector<int> hel0_evcount;
  std::vector<long long> hel0_evnum;
  hel0_evcount.reserve(evcount_h0.size());
  hel0_evnum.reserve(evcount_h0.size());
  for (size_t i = 0; i < evcount_h0.size(); ++i) {
    if (cfg.use_rolling_stability) {
      if (stable_set_h0.find(evcount_h0[i]) == stable_set_h0.end()) continue;
    }
    out.hel0_charge_uC += charge_h0[i];
    const int evcount = evcount_h0[i];
    if (evcount >= 1 && static_cast<size_t>(evcount) <= t_gevnum_sorted.size()) {
      hel0_evcount.push_back(evcount);
      hel0_evnum.push_back(t_gevnum_sorted[static_cast<size_t>(evcount - 1)]);
    }
  }

  out.hel0_evnum_ranges = ranges_from_evcount_evnum(hel0_evcount, hel0_evnum);

  return out;
}

static SegmentResult analyze_segment_charge(
    int run,
    const std::string &path,
    double expected_current,
    double current_window,
  double current_window_frac,
  bool use_fractional_window,
  double current_correction,
  bool use_raw_bcm4a,
  double bcm4a_gain,
  double bcm4a_offset) {

  SegmentResult s;
  s.run = run;
  s.segment = parse_segment_from_filename(path, run);
  s.file = fs::path(path).filename().string();
  s.expected_current_uA = expected_current;

  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) {
    s.missing_branches++;
    return s;
  }

  auto *tsh = dynamic_cast<TTree *>(file.Get("TSH"));
  if (!tsh) {
    s.missing_branches++;
    return s;
  }

  const char *scaler_branches[] = {
    "H.1MHz.scalerTime",
    "H.BCM4A.scalerCurrent",
    "H.BCM4A.scaler",
    "H.hTRIG4.scaler",
    "H.EDTM.scaler",
    "H.hL1ACCP.scaler"
  };
  for (const auto *b : scaler_branches) {
    if (!has_branch(tsh, b)) {
      s.missing_branches++;
      return s;
    }
  }

  tsh->SetBranchStatus("*", 0);
  tsh->SetBranchStatus("H.1MHz.scalerTime", 1);
  tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
  tsh->SetBranchStatus("H.BCM4A.scaler", 1);

  double tMHz = 0.0;
  double bcmCurrent = 0.0;
  double bcm4aScaler = 0.0;
  tsh->SetBranchAddress("H.1MHz.scalerTime", &tMHz);
  tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent);
  tsh->SetBranchAddress("H.BCM4A.scaler", &bcm4aScaler);

  const Long64_t n = tsh->GetEntries();
  if (n < 2) {
    return s;
  }

  double charge = 0.0;
  double current_time_weighted = 0.0;
  double time_total_cut = 0.0;

  double window_lo = std::numeric_limits<double>::quiet_NaN();
  double window_hi = std::numeric_limits<double>::quiet_NaN();
  compute_window_limits(expected_current,
                        current_window,
                        current_window_frac,
                        use_fractional_window,
                        window_lo,
                        window_hi);

  tsh->GetEntry(0);
  double prev_t = tMHz;
  double prev_bcm4a_scaler = bcm4aScaler;

  for (Long64_t i = 1; i < n; ++i) {
    tsh->GetEntry(i);
    const double dt = tMHz - prev_t;
    if (!std::isfinite(dt) || dt < 0.0) {
      s.negative_dt_intervals++;
      prev_t = tMHz;
      continue;
    }

    const double effective_gain = (std::isfinite(bcm4a_gain) && std::abs(bcm4a_gain) > 1e-12)
        ? bcm4a_gain
        : std::numeric_limits<double>::quiet_NaN();
    const double d_bcm4a = bcm4aScaler - prev_bcm4a_scaler;
    const double raw_rate = d_bcm4a / dt;
    const double calibrated_current = use_raw_bcm4a
      ? ((raw_rate - bcm4a_offset) / effective_gain)
      : bcmCurrent;
    const double corrected_current = calibrated_current + current_correction;
    const bool in_window = (calibrated_current >= window_lo && calibrated_current <= window_hi);
    if (in_window) {
      charge += corrected_current * dt;
      current_time_weighted += corrected_current * dt;
      time_total_cut += dt;
    }

    prev_t = tMHz;
    prev_bcm4a_scaler = bcm4aScaler;
  }

  s.charge_uC = charge;
  s.avg_current_uA = safe_div(current_time_weighted, time_total_cut);
  return s;
}

static bool compute_run_median_current(
    const std::vector<std::string> &files,
    double &median_current_uA,
  bool use_raw_bcm4a,
  double bcm4a_gain,
  double bcm4a_offset,
    bool show_progress) {
  std::vector<double> currents;

  for (size_t file_idx = 0; file_idx < files.size(); ++file_idx) {
    const auto &path = files[file_idx];
    if (show_progress) {
      std::cout << "[median] file " << (file_idx + 1) << "/" << files.size()
                << ": " << fs::path(path).filename().string() << "\n";
      std::cout.flush();
    }

    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) {
      continue;
    }

    auto *tsh = dynamic_cast<TTree *>(file.Get("TSH"));
    if (!tsh || !has_branch(tsh, "H.BCM4A.scalerCurrent")) {
      continue;
    }
    if (use_raw_bcm4a && (!has_branch(tsh, "H.BCM4A.scaler") || !has_branch(tsh, "H.1MHz.scalerTime"))) {
      continue;
    }

    tsh->SetBranchStatus("*", 0);
    tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
    if (use_raw_bcm4a) {
      tsh->SetBranchStatus("H.BCM4A.scaler", 1);
      tsh->SetBranchStatus("H.1MHz.scalerTime", 1);
    }

    double bcmCurrent = 0.0;
    double tMHz = 0.0;
    double bcm4aScaler = 0.0;
    tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent);
    if (use_raw_bcm4a) {
      tsh->SetBranchAddress("H.BCM4A.scaler", &bcm4aScaler);
      tsh->SetBranchAddress("H.1MHz.scalerTime", &tMHz);
    }

    const Long64_t n = tsh->GetEntries();
    if (n < 1) {
      continue;
    }

    tsh->GetEntry(0);
    double prev_t = tMHz;
    double prev_bcm4a_scaler = bcm4aScaler;

    for (Long64_t i = 0; i < n; ++i) {
      tsh->GetEntry(i);
      const double effective_gain = (std::isfinite(bcm4a_gain) && std::abs(bcm4a_gain) > 1e-12)
          ? bcm4a_gain
          : std::numeric_limits<double>::quiet_NaN();
      const double calibrated_current = use_raw_bcm4a
          ? ([&]() {
              const double dt = tMHz - prev_t;
              const double d_bcm4a = bcm4aScaler - prev_bcm4a_scaler;
              if (!std::isfinite(dt) || dt <= 0.0 || !std::isfinite(d_bcm4a)) {
                return std::numeric_limits<double>::quiet_NaN();
              }
              const double raw_rate = d_bcm4a / dt;
              return (raw_rate - bcm4a_offset) / effective_gain;
            })()
          : bcmCurrent;
      prev_t = tMHz;
      prev_bcm4a_scaler = bcm4aScaler;
      if (!std::isfinite(calibrated_current)) {
        continue;
      }
      if (calibrated_current <= 2.0) {
        continue;
      }
      currents.push_back(calibrated_current);
    }
  }

  if (currents.empty()) {
    return false;
  }

  const size_t n = currents.size();
  const size_t mid = n / 2;
  std::nth_element(currents.begin(), currents.begin() + mid, currents.end());
  double m = currents[mid];

  if (n % 2 == 0) {
    std::nth_element(currents.begin(), currents.begin() + mid - 1, currents.end());
    m = 0.5 * (m + currents[mid - 1]);
  }

  median_current_uA = m;
  return true;
}

static RunResult analyze_run(
    int run,
    const std::vector<std::string> &files,
  const std::string &report_dir,
  bool use_report_efficiencies,
    double expected_current,
    double ps_factor,
    double current_window,
  double current_window_frac,
  bool use_fractional_window,
  double current_correction,
  double beam_on_threshold_uA,
  bool use_raw_bcm4a,
  double bcm4a_gain,
  double bcm4a_offset,
  double pid_cer_npe_min,
  double pid_cal_etotnorm_min,
  double pid_cal_etotnorm_max,
  double hodo_track_max_chisq,
  double hodo_track_slop_cm,
  bool pid_require_track,
  bool use_rolling_stability,
  int stability_window,
  double stability_frac_range,
  int trig_number,
  const std::string &prescale_token,
  bool show_progress) {

  RunResult r;
  r.run = run;
  r.n_files = static_cast<int>(files.size());
  r.expected_current_uA = expected_current;
  r.ps_factor = ps_factor;
  r.prescale_token = prescale_token;

  char htrig_scaler_branch[64];
  char htrig_tdc_branch[64];
  std::snprintf(htrig_scaler_branch, sizeof(htrig_scaler_branch), "H.hTRIG%d.scaler", trig_number);
  std::snprintf(htrig_tdc_branch, sizeof(htrig_tdc_branch), "T.hms.hTRIG%d_tdcTimeRaw", trig_number);
  r.which_TRIG = htrig_scaler_branch;

  double charge_uC = 0.0;
  double charge_above_threshold_uC = 0.0;
  double charge_above_3p0_uC = 0.0;
  double current_time_weighted = 0.0;
  double time_total_cut = 0.0;
  long long newgen_num = 0;
  long long newgen_den = 0;
  long long newgen_trig6_num = 0;
  long long newgen_trig6_den = 0;
  long long newgen_elt_evt_num = 0;
  long long newgen_elt_evt_den = 0;
  long long newgen_elt_evt_noedtm_den = 0;
  long long pid_den = 0;
  long long pid_cer_num = 0;
  long long pid_cal_num = 0;
  long long pid_both_num = 0;
  long long pid_cer_tag_den = 0;
  long long pid_cal_tag_den = 0;
  long long hms_tracking_should = 0;
  long long hms_tracking_did = 0;
  long long hms_hodo_3of4_den = 0;
  long long hms_hodo_3of4_num = 0;
  long long hms_hodo_3of4_nhits_den = 0;
  long long hms_hodo_3of4_nhits_num = 0;
  long long hms_hodo_legacy_track_den = 0;
  long long hms_hodo_should_1x = 0;
  long long hms_hodo_should_1y = 0;
  long long hms_hodo_should_2x = 0;
  long long hms_hodo_should_2y = 0;
  long long hms_hodo_did_1x = 0;
  long long hms_hodo_did_1y = 0;
  long long hms_hodo_did_2x = 0;
  long long hms_hodo_did_2y = 0;
  double report_hdide_sum = 0.0;
  double report_hscinshoulde_sum = 0.0;
  double report_hodo_weight_sum = 0.0;
  double report_hodo_num_sum = 0.0;
  long long og3_numerator = 0;
  long long og3_denominator = 0;
  double newgen_elt_scaler_num = 0.0;
  double newgen_elt_scaler_den = 0.0;
  bool have_newgen_elt_scaler_num = false;
  bool have_newgen_elt_scaler_den = false;

  for (size_t file_idx = 0; file_idx < files.size(); ++file_idx) {
    const auto &path = files[file_idx];
    if (show_progress) {
      std::cout << "[run " << run << "] file " << (file_idx + 1) << "/" << files.size()
                << ": " << fs::path(path).filename().string() << "\n";
      std::cout.flush();
    }

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

    if (use_report_efficiencies) {
      const int segment = parse_segment_from_filename(path, run);
      if (segment >= 0) {
        const std::string report_path = report_dir + "/coin_NPS_HMS_report_" +
                                        std::to_string(run) + "_" +
                                        std::to_string(segment) + "_1_-1.report";
        const auto rep = parse_report_segment_metrics(report_path);
        if (rep.have_hdide && rep.have_hscinshoulde) {
          report_hdide_sum += rep.hdide;
          report_hscinshoulde_sum += rep.hscinshoulde;
        }
        if (rep.have_hodo3of4) {
          const double w = rep.have_htrig && std::isfinite(rep.htrig) && rep.htrig > 0.0 ? rep.htrig : 1.0;
          report_hodo_weight_sum += w;
          report_hodo_num_sum += rep.hodo3of4 * w;
        }
      }
    }

    std::vector<std::pair<long long, long long>> stable_evnum_ranges;
    if (use_rolling_stability) {
      Config stability_cfg;
      stability_cfg.current_window = current_window;
      stability_cfg.current_window_frac = current_window_frac;
      stability_cfg.use_fractional_window = use_fractional_window;
      stability_cfg.current_correction = current_correction;
      stability_cfg.use_rolling_stability = true;
      stability_cfg.stability_window = stability_window;
      stability_cfg.stability_frac_range = stability_frac_range;

      const HelicitySelectionResult hel = analyze_segment_helicity(path, expected_current, stability_cfg);
      stable_evnum_ranges = parse_evnum_ranges(hel.good_evnum_ranges);
    }

    const bool has_hpre100_scaler_cut = has_branch(tsh, "H.hPRE100.scalerCut");
    const bool has_hpre50_scaler_cut = has_branch(tsh, "H.hPRE50.scalerCut");
    const bool has_hpre40_scaler_cut = has_branch(tsh, "H.hPRE40.scalerCut");
    const char *hpre_ref_scaler_cut_branch = has_hpre50_scaler_cut
      ? "H.hPRE50.scalerCut"
      : (has_hpre40_scaler_cut ? "H.hPRE40.scalerCut" : nullptr);

    const char *scaler_branches[] = {
      "H.1MHz.scalerTime",
      "H.BCM4A.scalerCurrent",
      "H.BCM4A.scaler",
      htrig_scaler_branch,
      "H.EDTM.scaler",
      "H.hL1ACCP.scaler"
    };

    bool missing_scaler = false;
    for (const auto *b : scaler_branches) {
      if (!has_branch(tsh, b)) {
        missing_scaler = true;
        break;
      }
    }
    if (missing_scaler) {
      r.missing_branches++;
      continue;
    }

    tsh->SetBranchStatus("*", 0);
    tsh->SetBranchStatus("H.1MHz.scalerTime", 1);
    tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
    tsh->SetBranchStatus("H.BCM4A.scaler", 1);
    tsh->SetBranchStatus(htrig_scaler_branch, 1);
    tsh->SetBranchStatus("H.EDTM.scaler", 1);
    tsh->SetBranchStatus("H.hL1ACCP.scaler", 1);
    if (has_hpre100_scaler_cut) {
      tsh->SetBranchStatus("H.hPRE100.scalerCut", 1);
    }
    if (hpre_ref_scaler_cut_branch != nullptr) {
      tsh->SetBranchStatus(hpre_ref_scaler_cut_branch, 1);
    }
    tsh->SetCacheSize(20 * 1024 * 1024);
    tsh->AddBranchToCache("H.1MHz.scalerTime", true);
    tsh->AddBranchToCache("H.BCM4A.scalerCurrent", true);
    tsh->AddBranchToCache("H.BCM4A.scaler", true);
    tsh->AddBranchToCache(htrig_scaler_branch, true);
    tsh->AddBranchToCache("H.EDTM.scaler", true);
    tsh->AddBranchToCache("H.hL1ACCP.scaler", true);
    if (has_hpre100_scaler_cut) {
      tsh->AddBranchToCache("H.hPRE100.scalerCut", true);
    }
    if (hpre_ref_scaler_cut_branch != nullptr) {
      tsh->AddBranchToCache(hpre_ref_scaler_cut_branch, true);
    }

    double tMHz = 0.0;
    double bcmCurrent = 0.0;
    double bcm4aScaler = 0.0;
    double htrig = 0.0;
    double edtm = 0.0;
    double l1accp = 0.0;
    double hpre100_scaler_cut = 0.0;
    double hpre_ref_scaler_cut = 0.0;

    tsh->SetBranchAddress("H.1MHz.scalerTime", &tMHz);
    tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent);
    tsh->SetBranchAddress("H.BCM4A.scaler", &bcm4aScaler);
    tsh->SetBranchAddress(htrig_scaler_branch, &htrig);
    tsh->SetBranchAddress("H.EDTM.scaler", &edtm);
    tsh->SetBranchAddress("H.hL1ACCP.scaler", &l1accp);
    if (has_hpre100_scaler_cut) {
      tsh->SetBranchAddress("H.hPRE100.scalerCut", &hpre100_scaler_cut);
    }
    if (hpre_ref_scaler_cut_branch != nullptr) {
      tsh->SetBranchAddress(hpre_ref_scaler_cut_branch, &hpre_ref_scaler_cut);
    }

    const Long64_t n = tsh->GetEntries();
    if (n < 2) continue;

    const Long64_t scaler_progress_step = std::max<Long64_t>(1, n / 10);

    tsh->GetEntry(0);
    double prev_t = tMHz;
    double prev_bcm4a_scaler = bcm4aScaler;
    double prev_htrig = htrig;
    double prev_edtm = edtm;
    double prev_l1accp = l1accp;
    double prev_hpre100_scaler_cut = hpre100_scaler_cut;
    double prev_hpre_ref_scaler_cut = hpre_ref_scaler_cut;

    double window_lo = std::numeric_limits<double>::quiet_NaN();
    double window_hi = std::numeric_limits<double>::quiet_NaN();
    compute_window_limits(expected_current,
                current_window,
                current_window_frac,
                use_fractional_window,
                window_lo,
                window_hi);

    for (Long64_t i = 1; i < n; ++i) {
      if (show_progress && (i % scaler_progress_step == 0 || i == n - 1)) {
        const int pct = static_cast<int>((100.0 * i) / std::max<Long64_t>(1, n - 1));
        std::cout << "  [run " << run << "] scaler " << pct << "% ("
                  << i << "/" << (n - 1) << ")\r";
        std::cout.flush();
      }

      tsh->GetEntry(i);

      double dt = tMHz - prev_t;
      if (!std::isfinite(dt) || dt < 0.0) {
        r.negative_dt_intervals++;
        prev_t = tMHz;
        prev_bcm4a_scaler = bcm4aScaler;
        prev_htrig = htrig;
        prev_edtm = edtm;
        prev_l1accp = l1accp;
        prev_hpre100_scaler_cut = hpre100_scaler_cut;
        prev_hpre_ref_scaler_cut = hpre_ref_scaler_cut;
        continue;
      }

      double d_htrig = htrig - prev_htrig;
      double d_edtm = edtm - prev_edtm;
      double d_accp = l1accp - prev_l1accp;
        double d_hpre100_cut =
          has_hpre100_scaler_cut ? (hpre100_scaler_cut - prev_hpre100_scaler_cut) : 0.0;
        double d_hpre_ref_cut = (hpre_ref_scaler_cut_branch != nullptr)
          ? (hpre_ref_scaler_cut - prev_hpre_ref_scaler_cut)
          : 0.0;

        if (d_htrig < 0.0 || d_edtm < 0.0 || d_accp < 0.0 || d_hpre100_cut < 0.0 || d_hpre_ref_cut < 0.0) {
        r.non_monotonic_scaler_steps++;
      }

      d_htrig = std::max(0.0, d_htrig);
      d_edtm = std::max(0.0, d_edtm);
      d_accp = std::max(0.0, d_accp);
      d_hpre100_cut = std::max(0.0, d_hpre100_cut);
      d_hpre_ref_cut = std::max(0.0, d_hpre_ref_cut);

      // Reject pathological scaler-jump intervals (e.g. single-step ~2^32 spikes)
      // that would otherwise dominate livetime denominators.
      const double rate_htrig = d_htrig / dt;
      const double rate_edtm = d_edtm / dt;
      const double rate_accp = d_accp / dt;
      const bool suspicious_scaler_jump =
          (d_htrig > 1.0e8) || (d_edtm > 1.0e8) || (d_accp > 1.0e8) ||
          (rate_htrig > 1.0e7) || (rate_edtm > 1.0e7) || (rate_accp > 1.0e7);
      if (suspicious_scaler_jump) {
        r.suspicious_scaler_jump_steps++;
        prev_t = tMHz;
        prev_bcm4a_scaler = bcm4aScaler;
        prev_htrig = htrig;
        prev_edtm = edtm;
        prev_l1accp = l1accp;
        prev_hpre100_scaler_cut = hpre100_scaler_cut;
        prev_hpre_ref_scaler_cut = hpre_ref_scaler_cut;
        continue;
      }

      // OG-style electronic livetime uses PRE gate scalerCut ratios.
      if (has_hpre100_scaler_cut) {
        newgen_elt_scaler_num += d_hpre100_cut;
        have_newgen_elt_scaler_num = true;
      }
      if (hpre_ref_scaler_cut_branch != nullptr) {
        newgen_elt_scaler_den += d_hpre_ref_cut;
        have_newgen_elt_scaler_den = true;
      }

      r.scaler_edtm_no_cut += d_edtm;
      r.trig_accp_no_cut += d_accp;

      const double effective_gain = (std::isfinite(bcm4a_gain) && std::abs(bcm4a_gain) > 1e-12)
          ? bcm4a_gain
          : std::numeric_limits<double>::quiet_NaN();
      const double d_bcm4a = bcm4aScaler - prev_bcm4a_scaler;
      const double raw_rate = d_bcm4a / dt;
      const double calibrated_current = use_raw_bcm4a
          ? ((raw_rate - bcm4a_offset) / effective_gain)
          : bcmCurrent;
      double corrected_current = calibrated_current + current_correction;
      const bool in_window = (calibrated_current >= window_lo && calibrated_current <= window_hi);
      const bool above_threshold = (calibrated_current > beam_on_threshold_uA);
      if (in_window) {
        charge_uC += corrected_current * dt;
        current_time_weighted += corrected_current * dt;
        time_total_cut += dt;
        r.scaler_htrig4_total += d_htrig;
        r.scaler_edtm_total += d_edtm;
        r.trig_accp_total += d_accp;
      }
      if (above_threshold) {
        charge_above_threshold_uC += corrected_current * dt;
      }
      if (calibrated_current > 3.0) {
        charge_above_3p0_uC += corrected_current * dt;
      }

      prev_t = tMHz;
      prev_bcm4a_scaler = bcm4aScaler;
      prev_htrig = htrig;
      prev_edtm = edtm;
      prev_l1accp = l1accp;
      prev_hpre100_scaler_cut = hpre100_scaler_cut;
      prev_hpre_ref_scaler_cut = hpre_ref_scaler_cut;
    }
    if (show_progress) {
      std::cout << "\n";
      std::cout.flush();
    }

    const char *trig6_tdc_branch = has_branch(tevt, "T.hms.npsTRIG6_tdcTimeRaw")
      ? "T.hms.npsTRIG6_tdcTimeRaw"
      : (has_branch(tevt, "T.hms.hTRIG6_tdcTimeRaw") ? "T.hms.hTRIG6_tdcTimeRaw" : nullptr);
    const char *selected_htrig_tdc_branch =
      (trig_number == 6 && trig6_tdc_branch != nullptr) ? trig6_tdc_branch : htrig_tdc_branch;

    bool missing_evt =
        !has_branch(tevt, "T.hms.hEDTM_tdcTimeRaw") ||
      !has_branch(tevt, selected_htrig_tdc_branch) ||
      (trig6_tdc_branch == nullptr) ||
        !has_branch(tevt, "T.hms.hTRIG3_tdcTimeRaw") ||
        !has_branch(tevt, "T.hms.hPRE40_tdcTimeRaw") ||
        !has_branch(tevt, "T.hms.npsTRIG1_tdcTimeRaw") ||
        !has_branch(tevt, "T.hms.hTRIG4_tdcTimeRaw") ||
        !has_branch(tevt, "H.BCM4A.scalerCurrent") ||
        (use_rolling_stability && !has_branch(tevt, "g.evnum"));
    if (missing_evt) {
      r.missing_branches++;
      continue;
    }

    const bool have_pid_branches =
        has_branch(tevt, "H.cer.npeSum") &&
      has_branch(tevt, "H.cal.etotnorm") &&
      has_branch(tevt, "H.dc.ntrack") &&
      has_branch(tevt, "H.gtr.dp") &&
      has_branch(tevt, "H.gtr.th") &&
      has_branch(tevt, "H.gtr.ph");

    const bool have_tracking_hodo_branches =
      have_pid_branches &&
      has_branch(tevt, "H.hod.goodscinhit") &&
      has_branch(tevt, "H.hod.betanotrack");

    const bool have_hodo_3of4_branches =
      has_branch(tevt, "H.hod.1x.nhits") &&
      has_branch(tevt, "H.hod.1y.nhits") &&
      has_branch(tevt, "H.hod.2x.nhits") &&
      has_branch(tevt, "H.hod.2y.nhits") &&
      has_branch(tevt, "H.dc.ntrack");

    const bool have_hodo_legacy_branches =
      have_hodo_3of4_branches &&
      has_branch(tevt, "Ndata.H.dc.track_chisq") &&
      has_branch(tevt, "H.dc.track_chisq") &&
      has_branch(tevt, "H.hod.1x.DiffDisTrack") &&
      has_branch(tevt, "H.hod.1y.DiffDisTrack") &&
      has_branch(tevt, "H.hod.2x.DiffDisTrack") &&
      has_branch(tevt, "H.hod.2y.DiffDisTrack");

    tevt->SetBranchStatus("*", 0);
    tevt->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
    tevt->SetBranchStatus(selected_htrig_tdc_branch, 1);
    tevt->SetBranchStatus(trig6_tdc_branch, 1);
    tevt->SetBranchStatus("T.hms.hTRIG3_tdcTimeRaw", 1);
    tevt->SetBranchStatus("T.hms.hPRE40_tdcTimeRaw", 1);
    tevt->SetBranchStatus("T.hms.npsTRIG1_tdcTimeRaw", 1);
    tevt->SetBranchStatus("T.hms.hTRIG4_tdcTimeRaw", 1);
    tevt->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
    if (use_rolling_stability) {
      tevt->SetBranchStatus("g.evnum", 1);
    }
    if (have_pid_branches) {
      tevt->SetBranchStatus("H.cer.npeSum", 1);
      tevt->SetBranchStatus("H.cal.etotnorm", 1);
      tevt->SetBranchStatus("H.dc.ntrack", 1);
      tevt->SetBranchStatus("H.gtr.dp", 1);
      tevt->SetBranchStatus("H.gtr.th", 1);
      tevt->SetBranchStatus("H.gtr.ph", 1);
    } else if (have_hodo_3of4_branches) {
      tevt->SetBranchStatus("H.dc.ntrack", 1);
    }
    if (have_tracking_hodo_branches) {
      tevt->SetBranchStatus("H.hod.goodscinhit", 1);
      tevt->SetBranchStatus("H.hod.betanotrack", 1);
    }
    if (have_hodo_3of4_branches) {
      tevt->SetBranchStatus("H.hod.1x.nhits", 1);
      tevt->SetBranchStatus("H.hod.1y.nhits", 1);
      tevt->SetBranchStatus("H.hod.2x.nhits", 1);
      tevt->SetBranchStatus("H.hod.2y.nhits", 1);
    }
    if (have_hodo_legacy_branches) {
      tevt->SetBranchStatus("Ndata.H.dc.track_chisq", 1);
      tevt->SetBranchStatus("H.dc.track_chisq", 1);
      tevt->SetBranchStatus("H.hod.1x.DiffDisTrack", 1);
      tevt->SetBranchStatus("H.hod.1y.DiffDisTrack", 1);
      tevt->SetBranchStatus("H.hod.2x.DiffDisTrack", 1);
      tevt->SetBranchStatus("H.hod.2y.DiffDisTrack", 1);
    }
    tevt->SetCacheSize(40 * 1024 * 1024);
    tevt->AddBranchToCache("T.hms.hEDTM_tdcTimeRaw", true);
    tevt->AddBranchToCache(selected_htrig_tdc_branch, true);
    tevt->AddBranchToCache(trig6_tdc_branch, true);
    tevt->AddBranchToCache("T.hms.hTRIG3_tdcTimeRaw", true);
    tevt->AddBranchToCache("T.hms.hPRE40_tdcTimeRaw", true);
    tevt->AddBranchToCache("T.hms.npsTRIG1_tdcTimeRaw", true);
    tevt->AddBranchToCache("T.hms.hTRIG4_tdcTimeRaw", true);
    tevt->AddBranchToCache("H.BCM4A.scalerCurrent", true);
    if (use_rolling_stability) {
      tevt->AddBranchToCache("g.evnum", true);
    }
    if (have_pid_branches) {
      tevt->AddBranchToCache("H.cer.npeSum", true);
      tevt->AddBranchToCache("H.cal.etotnorm", true);
      tevt->AddBranchToCache("H.dc.ntrack", true);
      tevt->AddBranchToCache("H.gtr.dp", true);
      tevt->AddBranchToCache("H.gtr.th", true);
      tevt->AddBranchToCache("H.gtr.ph", true);
    } else if (have_hodo_3of4_branches) {
      tevt->AddBranchToCache("H.dc.ntrack", true);
    }
    if (have_tracking_hodo_branches) {
      tevt->AddBranchToCache("H.hod.goodscinhit", true);
      tevt->AddBranchToCache("H.hod.betanotrack", true);
    }
    if (have_hodo_3of4_branches) {
      tevt->AddBranchToCache("H.hod.1x.nhits", true);
      tevt->AddBranchToCache("H.hod.1y.nhits", true);
      tevt->AddBranchToCache("H.hod.2x.nhits", true);
      tevt->AddBranchToCache("H.hod.2y.nhits", true);
    }
    if (have_hodo_legacy_branches) {
      tevt->AddBranchToCache("Ndata.H.dc.track_chisq", true);
      tevt->AddBranchToCache("H.dc.track_chisq", true);
      tevt->AddBranchToCache("H.hod.1x.DiffDisTrack", true);
      tevt->AddBranchToCache("H.hod.1y.DiffDisTrack", true);
      tevt->AddBranchToCache("H.hod.2x.DiffDisTrack", true);
      tevt->AddBranchToCache("H.hod.2y.DiffDisTrack", true);
    }

    double edtm_tdc = 0.0;
    double trig_tdc = 0.0;
    double trig6_tdc = 0.0;
    double trig3_tdc = 0.0;
    double pre40_tdc = 0.0;
    double npsTRIG1_tdc = 0.0;
    double hTRIG4_tdc = 0.0;
    double bcmCurrent_evt = 0.0;
    double g_evnum = 0.0;
    double cer_npe_sum = std::numeric_limits<double>::quiet_NaN();
    double cal_etotnorm = std::numeric_limits<double>::quiet_NaN();
    double hdc_ntrack = 0.0;
    double hms_dp = std::numeric_limits<double>::quiet_NaN();
    double hms_th = std::numeric_limits<double>::quiet_NaN();
    double hms_ph = std::numeric_limits<double>::quiet_NaN();
    double hod_goodscinhit = std::numeric_limits<double>::quiet_NaN();
    double hod_betanotrack = std::numeric_limits<double>::quiet_NaN();
    double hod_1x_nhits = std::numeric_limits<double>::quiet_NaN();
    double hod_1y_nhits = std::numeric_limits<double>::quiet_NaN();
    double hod_2x_nhits = std::numeric_limits<double>::quiet_NaN();
    double hod_2y_nhits = std::numeric_limits<double>::quiet_NaN();
    int ndata_hdc_track_chisq = 0;
    double hdc_track_chisq[100] = {0.0};
    double hod_1x_diffdis = std::numeric_limits<double>::quiet_NaN();
    double hod_1y_diffdis = std::numeric_limits<double>::quiet_NaN();
    double hod_2x_diffdis = std::numeric_limits<double>::quiet_NaN();
    double hod_2y_diffdis = std::numeric_limits<double>::quiet_NaN();
    const bool trig6_aliases_selected =
        std::string(selected_htrig_tdc_branch) == std::string(trig6_tdc_branch);
    const bool trig3_aliases_selected =
      std::string(selected_htrig_tdc_branch) == std::string("T.hms.hTRIG3_tdcTimeRaw");
    const bool htrig4_aliases_selected =
        std::string(selected_htrig_tdc_branch) == std::string("T.hms.hTRIG4_tdcTimeRaw");
    tevt->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtm_tdc);
    tevt->SetBranchAddress(selected_htrig_tdc_branch, &trig_tdc);
    if (!trig6_aliases_selected) {
      tevt->SetBranchAddress(trig6_tdc_branch, &trig6_tdc);
    }
    if (!trig3_aliases_selected) {
      tevt->SetBranchAddress("T.hms.hTRIG3_tdcTimeRaw", &trig3_tdc);
    }
    tevt->SetBranchAddress("T.hms.hPRE40_tdcTimeRaw", &pre40_tdc);
    tevt->SetBranchAddress("T.hms.npsTRIG1_tdcTimeRaw", &npsTRIG1_tdc);
    if (!htrig4_aliases_selected) {
      tevt->SetBranchAddress("T.hms.hTRIG4_tdcTimeRaw", &hTRIG4_tdc);
    }
    tevt->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent_evt);
    if (use_rolling_stability) {
      tevt->SetBranchAddress("g.evnum", &g_evnum);
    }
    if (have_pid_branches) {
      tevt->SetBranchAddress("H.cer.npeSum", &cer_npe_sum);
      tevt->SetBranchAddress("H.cal.etotnorm", &cal_etotnorm);
      tevt->SetBranchAddress("H.dc.ntrack", &hdc_ntrack);
      tevt->SetBranchAddress("H.gtr.dp", &hms_dp);
      tevt->SetBranchAddress("H.gtr.th", &hms_th);
      tevt->SetBranchAddress("H.gtr.ph", &hms_ph);
    } else if (have_hodo_3of4_branches) {
      tevt->SetBranchAddress("H.dc.ntrack", &hdc_ntrack);
    }
    if (have_tracking_hodo_branches) {
      tevt->SetBranchAddress("H.hod.goodscinhit", &hod_goodscinhit);
      tevt->SetBranchAddress("H.hod.betanotrack", &hod_betanotrack);
    }
    if (have_hodo_3of4_branches) {
      tevt->SetBranchAddress("H.hod.1x.nhits", &hod_1x_nhits);
      tevt->SetBranchAddress("H.hod.1y.nhits", &hod_1y_nhits);
      tevt->SetBranchAddress("H.hod.2x.nhits", &hod_2x_nhits);
      tevt->SetBranchAddress("H.hod.2y.nhits", &hod_2y_nhits);
    }
    if (have_hodo_legacy_branches) {
      tevt->SetBranchAddress("Ndata.H.dc.track_chisq", &ndata_hdc_track_chisq);
      tevt->SetBranchAddress("H.dc.track_chisq", hdc_track_chisq);
      tevt->SetBranchAddress("H.hod.1x.DiffDisTrack", &hod_1x_diffdis);
      tevt->SetBranchAddress("H.hod.1y.DiffDisTrack", &hod_1y_diffdis);
      tevt->SetBranchAddress("H.hod.2x.DiffDisTrack", &hod_2x_diffdis);
      tevt->SetBranchAddress("H.hod.2y.DiffDisTrack", &hod_2y_diffdis);
    }

    const Long64_t nevt = tevt->GetEntries();
    const Long64_t event_progress_step = std::max<Long64_t>(1, nevt / 10);
    for (Long64_t i = 0; i < nevt; ++i) {
      if (show_progress && (i % event_progress_step == 0 || i == nevt - 1)) {
        const int pct = static_cast<int>((100.0 * (i + 1)) / std::max<Long64_t>(1, nevt));
        std::cout << "  [run " << run << "] events " << pct << "% ("
                  << (i + 1) << "/" << nevt << ")\r";
        std::cout.flush();
      }

      tevt->GetEntry(i);
      if (trig6_aliases_selected) {
        trig6_tdc = trig_tdc;
      }
      if (trig3_aliases_selected) {
        trig3_tdc = trig_tdc;
      }
      if (htrig4_aliases_selected) {
        hTRIG4_tdc = trig_tdc;
      }
      const bool in_current_window_evt =
          std::isfinite(bcmCurrent_evt) && (bcmCurrent_evt >= window_lo) && (bcmCurrent_evt <= window_hi);
      // Keep event-level NewGEN counters on the same current-window gate as scaler
      // denominators. TSHelH rolling-stability masks are computed on helicity windows
      // and are too sparse to apply directly event-by-event here.
      const bool in_beam_good_evt = in_current_window_evt;

      if (in_beam_good_evt && hTRIG4_tdc > 0.0 && npsTRIG1_tdc > 0.0) {
        newgen_den++;
        if (edtm_tdc == 0.0) newgen_num++;
      }
      if (in_beam_good_evt && trig_tdc > 0.0) {
        newgen_trig6_den++;
        if (edtm_tdc == 0.0) newgen_trig6_num++;
      }

      const bool trig6_or_nps1 = (trig6_tdc > 0.0) || (npsTRIG1_tdc > 0.0);
      if (in_beam_good_evt && pre40_tdc > 0.0) {
        newgen_elt_evt_den++;
        if (trig6_or_nps1 && edtm_tdc == 0.0) newgen_elt_evt_num++;
        if (edtm_tdc == 0.0) newgen_elt_evt_noedtm_den++;
      }

      if (have_pid_branches) {
        const bool edt_ok = std::isfinite(edtm_tdc) && (edtm_tdc < 0.1);
        const bool track_ok = (!pid_require_track) || (std::isfinite(hdc_ntrack) && hdc_ntrack > 0.5);
        const bool kinematic_ok = std::isfinite(hms_dp) && std::isfinite(hms_th) && std::isfinite(hms_ph) &&
                     (std::fabs(hms_dp) <= 8.5) && (std::fabs(hms_th) <= 0.09) && (std::fabs(hms_ph) <= 0.09);
        const bool pid_den_sel =
          in_beam_good_evt && edt_ok && track_ok && kinematic_ok;
        const bool cer_pass = std::isfinite(cer_npe_sum) && (cer_npe_sum > pid_cer_npe_min);
        const bool cal_pass = std::isfinite(cal_etotnorm) &&
                    (cal_etotnorm > pid_cal_etotnorm_min) &&
                    (cal_etotnorm <= pid_cal_etotnorm_max);

        if (pid_den_sel) {
          pid_den++;
          if (cer_pass) {
            pid_cer_num++;
            pid_cal_tag_den++;
          }
          if (cal_pass) {
            pid_cal_num++;
            pid_cer_tag_den++;
          }
          if (cer_pass && cal_pass) {
            pid_both_num++;
          }
        }
      }

      if (have_tracking_hodo_branches) {
        const bool good_scin = std::isfinite(hod_goodscinhit) && (std::llround(hod_goodscinhit) == 1);
        const bool electron_pid_should =
            std::isfinite(cal_etotnorm) && (cal_etotnorm > 0.6) && (cal_etotnorm < 2.0) &&
            std::isfinite(cer_npe_sum) && (cer_npe_sum > 0.5);
        const bool track_should =
            in_beam_good_evt && good_scin && electron_pid_should;
        if (track_should) {
          hms_tracking_should++;
          if (std::isfinite(hdc_ntrack) && hdc_ntrack > 0.5) {
            hms_tracking_did++;
          }
        }
      }

      if (have_hodo_3of4_branches) {
        const bool edt_ok_hodo = std::isfinite(edtm_tdc) && (edtm_tdc < 0.1);
        const bool hodo_den_nhits =
            in_beam_good_evt && edt_ok_hodo &&
            std::isfinite(hdc_ntrack) && (hdc_ntrack > 0.5);
        const bool good1x = std::isfinite(hod_1x_nhits) && (hod_1x_nhits > 0.0) && (hod_1x_nhits < 3.0);
        const bool good1y = std::isfinite(hod_1y_nhits) && (hod_1y_nhits > 0.0) && (hod_1y_nhits < 3.0);
        const bool good2x = std::isfinite(hod_2x_nhits) && (hod_2x_nhits > 0.0) && (hod_2x_nhits < 3.0);
        const bool good2y = std::isfinite(hod_2y_nhits) && (hod_2y_nhits > 0.0) && (hod_2y_nhits < 3.0);

        if (hodo_den_nhits) {
          hms_hodo_3of4_nhits_den++;
          const int good_planes = static_cast<int>(good1x) + static_cast<int>(good1y) +
            static_cast<int>(good2x) + static_cast<int>(good2y);
          if (good_planes >= 3) {
            hms_hodo_3of4_nhits_num++;
          }
        }

        if (have_hodo_legacy_branches && hodo_den_nhits) {
          double best_track_chisq = std::numeric_limits<double>::infinity();
          const int n_track_chisq = std::max(0, std::min(ndata_hdc_track_chisq, 100));
          for (int it = 0; it < n_track_chisq; ++it) {
            const double chi = hdc_track_chisq[it];
            if (std::isfinite(chi) && chi >= 0.0 && chi < best_track_chisq) {
              best_track_chisq = chi;
            }
          }
          const bool chisq_ok = std::isfinite(best_track_chisq) && (best_track_chisq < hodo_track_max_chisq);
          if (chisq_ok) {
            hms_hodo_legacy_track_den++;
            const bool should1x = std::isfinite(hod_1x_diffdis) && (std::abs(hod_1x_diffdis) <= hodo_track_slop_cm);
            const bool should1y = std::isfinite(hod_1y_diffdis) && (std::abs(hod_1y_diffdis) <= hodo_track_slop_cm);
            const bool should2x = std::isfinite(hod_2x_diffdis) && (std::abs(hod_2x_diffdis) <= hodo_track_slop_cm);
            const bool should2y = std::isfinite(hod_2y_diffdis) && (std::abs(hod_2y_diffdis) <= hodo_track_slop_cm);

            if (should1x) {
              hms_hodo_should_1x++;
              if (good1x) hms_hodo_did_1x++;
            }
            if (should1y) {
              hms_hodo_should_1y++;
              if (good1y) hms_hodo_did_1y++;
            }
            if (should2x) {
              hms_hodo_should_2x++;
              if (good2x) hms_hodo_did_2x++;
            }
            if (should2y) {
              hms_hodo_should_2y++;
              if (good2y) hms_hodo_did_2y++;
            }
          }
        }
      }

      bool edtm_acc = (edtm_tdc != 0.0);
      bool trig4_acc = (trig_tdc != 0.0);
      bool trig4_phy = trig4_acc && !edtm_acc;

      if (in_beam_good_evt) {
        if (edtm_acc) r.t_edtm_accepted++;
        if (trig4_acc) r.t_htrig4_all_accepted++;
        if (trig4_phy) r.t_htrig4_phy_accepted++;
      }

      if (in_beam_good_evt && pre40_tdc > 0.0) {
        og3_denominator++;
        if (trig3_tdc > 0.0 && edtm_tdc == 0.0) og3_numerator++;
      }
    }
    if (show_progress) {
      std::cout << "\n";
      std::cout.flush();
    }
  }

  if (trig_number == 6) {
    r.NewGEN_nps_coin_livetime = safe_div(static_cast<double>(newgen_num), static_cast<double>(newgen_den));
    r.NewGEN_nps_coin_livetime_trig6 =
        safe_div(static_cast<double>(newgen_trig6_num), static_cast<double>(newgen_trig6_den));
  } else {
    r.NewGEN_nps_coin_livetime = std::numeric_limits<double>::quiet_NaN();
    r.NewGEN_nps_coin_livetime_trig6 = std::numeric_limits<double>::quiet_NaN();
  }
  if (have_newgen_elt_scaler_num && have_newgen_elt_scaler_den) {
    // OG-style electronic livetime from pretrigger gate ratios, using PRE50 when
    // available and PRE40 otherwise, with PRE100 as the comparison gate.
    r.NewGEN_electronic_livetime = safe_div(newgen_elt_scaler_num, newgen_elt_scaler_den);
  } else {
    r.NewGEN_electronic_livetime = std::numeric_limits<double>::quiet_NaN();
  }
  r.NewGEN_electronic_livetime_evt =
      safe_div(static_cast<double>(newgen_elt_evt_num), static_cast<double>(newgen_elt_evt_den));
  r.NewGEN_electronic_livetime_evt_noedtmden =
      safe_div(static_cast<double>(newgen_elt_evt_num), static_cast<double>(newgen_elt_evt_noedtm_den));
  r.HMS_pid_den_events = pid_den;
  r.HMS_pid_eff = safe_div(static_cast<double>(pid_both_num), static_cast<double>(pid_den));
  r.HMS_cer_eff = safe_div(static_cast<double>(pid_cer_num), static_cast<double>(pid_den));
  r.HMS_cal_eff = safe_div(static_cast<double>(pid_cal_num), static_cast<double>(pid_den));
  r.HMS_cer_eff_tag_cal = safe_div(static_cast<double>(pid_both_num), static_cast<double>(pid_cer_tag_den));
  r.HMS_cal_eff_tag_cer = safe_div(static_cast<double>(pid_both_num), static_cast<double>(pid_cal_tag_den));
  if (report_hscinshoulde_sum > 0.0 && report_hdide_sum >= 0.0) {
    r.HMS_tracking_should_events = static_cast<long long>(std::llround(report_hscinshoulde_sum));
    r.HMS_tracking_did_events = static_cast<long long>(std::llround(report_hdide_sum));
    r.HMS_tracking_eff = safe_div(report_hdide_sum, report_hscinshoulde_sum);
  } else {
    r.HMS_tracking_should_events = hms_tracking_should;
    r.HMS_tracking_did_events = hms_tracking_did;
    r.HMS_tracking_eff = safe_div(static_cast<double>(hms_tracking_did), static_cast<double>(hms_tracking_should));
  }
  if (r.HMS_tracking_should_events > 0) {
    const double p = r.HMS_tracking_eff;
    r.HMS_tracking_eff_err = std::sqrt(std::max(0.0, p * (1.0 - p)) /
                                       static_cast<double>(r.HMS_tracking_should_events));
  } else {
    r.HMS_tracking_eff_err = std::numeric_limits<double>::quiet_NaN();
  }
  if (report_hodo_weight_sum > 0.0) {
    r.HMS_hodo_3of4_den_events = static_cast<long long>(std::llround(report_hodo_weight_sum));
    r.HMS_hodo_3of4_eff = safe_div(report_hodo_num_sum, report_hodo_weight_sum);
  } else {
    const double e1 = safe_div(static_cast<double>(hms_hodo_did_1x), static_cast<double>(hms_hodo_should_1x));
    const double e2 = safe_div(static_cast<double>(hms_hodo_did_1y), static_cast<double>(hms_hodo_should_1y));
    const double e3 = safe_div(static_cast<double>(hms_hodo_did_2x), static_cast<double>(hms_hodo_should_2x));
    const double e4 = safe_div(static_cast<double>(hms_hodo_did_2y), static_cast<double>(hms_hodo_should_2y));
    if (std::isfinite(e1) && std::isfinite(e2) && std::isfinite(e3) && std::isfinite(e4)) {
      const double p3of4 =
        e1 * e2 * e3 +
        e1 * e2 * e4 +
        e1 * e3 * e4 +
        e2 * e3 * e4 -
        3.0 * e1 * e2 * e3 * e4;
      r.HMS_hodo_3of4_den_events = hms_hodo_legacy_track_den;
      r.HMS_hodo_3of4_eff = p3of4;
    } else {
      r.HMS_hodo_3of4_den_events = hms_hodo_3of4_den;
      r.HMS_hodo_3of4_eff = safe_div(static_cast<double>(hms_hodo_3of4_num), static_cast<double>(hms_hodo_3of4_den));
    }
  }
  r.HMS_hodo_3of4_nhits_den_events = hms_hodo_3of4_nhits_den;
  r.HMS_hodo_3of4_nhits_eff = safe_div(static_cast<double>(hms_hodo_3of4_nhits_num),
                                       static_cast<double>(hms_hodo_3of4_nhits_den));
  r.OG_pretrig_livetime_3_4 = safe_div(static_cast<double>(og3_numerator), static_cast<double>(og3_denominator));
  if (og3_denominator > 0) {
    const double p = r.OG_pretrig_livetime_3_4;
    r.OG_pretrig_livetime_3_4_err = std::sqrt(p * (1.0 - p) / static_cast<double>(og3_denominator));
  } else {
    r.OG_pretrig_livetime_3_4_err = std::numeric_limits<double>::quiet_NaN();
  }

  r.charge_uC = charge_uC;
  r.charge_above_threshold_uC = charge_above_threshold_uC;
  r.charge_above_3p0_uC = charge_above_3p0_uC;
  r.avg_current_uA = safe_div(current_time_weighted, time_total_cut);
  r.beam_time = time_total_cut;
  r.beam_on_percent_edtm = clamp01(safe_div(r.scaler_edtm_total, r.scaler_edtm_no_cut));
  r.beam_on_percent_trig_accp = clamp01(safe_div(r.trig_accp_total, r.trig_accp_no_cut));

  r.clta_livetime_tsh = safe_div(r.trig_accp_total, r.scaler_htrig4_total) * ps_factor;
  r.clta_livetime_tdc = safe_div(r.t_htrig4_all_accepted, r.scaler_htrig4_total)
                      * ps_factor * r.beam_on_percent_trig_accp;

  const double denom_phy = r.scaler_htrig4_total - r.scaler_edtm_total;
  r.cltp_livetime_tdc = safe_div(r.t_htrig4_phy_accepted, denom_phy)
                      * ps_factor * r.beam_on_percent_trig_accp;

  const double edtm_denom = safe_div(r.scaler_edtm_total, r.beam_on_percent_edtm);
  r.tlt_livetime_edtm = safe_div(r.t_edtm_accepted, edtm_denom) * ps_factor;
    r.tlt_livetime_edtm_sym_beamcut =
      safe_div(static_cast<double>(r.t_edtm_accepted), r.scaler_edtm_total) * ps_factor;
    r.tlt_livetime_edtm_nocut_den =
      safe_div(static_cast<double>(r.t_edtm_accepted), r.scaler_edtm_no_cut) * ps_factor;
    r.tlt_livetime_edtm_ratio_sym =
      safe_div(r.tlt_livetime_edtm, r.tlt_livetime_edtm_sym_beamcut);
    r.tlt_livetime_edtm_ratio_nocut =
      safe_div(r.tlt_livetime_edtm, r.tlt_livetime_edtm_nocut_den);

  return r;
}

static void print_usage(const char *prog) {
  std::cerr
      << "Usage: " << prog << " --runs <list...> [options]\n"
      << "Options:\n"
      << "  --root-dir <path>            Default: /cache/hallc/c-nps/analysis/pass2/replays/production/\n"
      << "  --db <path>                  Default: /group/nps/jpcrafts/Pi_0/DataBase_production_runs_newBCMOffset.txt\n"
      << "  --trigger-overrides <path>   Default: /group/nps/jpcrafts/Pi_0/run_trigger_prescale_overrides.tsv\n"
      << "  --use-report-efficiencies    Optional: source tracking/hodo efficiencies from report files\n"
      << "  --report-dir <path>          Default: /work/hallc/nps/nps-ana/REPORT_OUTPUT_pass2/COIN_old (used only with --use-report-efficiencies)\n"
      << "  --out-csv <path>             Default: livetime_results_recomputed.tsv\n"
      << "  --out-segment-csv <path>     Default: segment_charge_results_recomputed.tsv\n"
      << "  --default-ps <value>         Default: 1.0\n"
      << "  --current-window <uA>        Uses center ± width (Default: 1.5, i.e. ±1.5 uA)\n"
      << "  --current-window-frac <f>    Uses fractional half-width (Default: 0.10, i.e. ±10%)\n"
      << "  --use-absolute-window        Use --current-window instead of fractional window\n"
      << "  --current-correction <uA>    Default: 0.0\n"
      << "  --beam-on-threshold <uA>     Additional charge metric uses current > threshold (Default: 2.5)\n"
      << "  --use-raw-bcm4a              Use H.BCM4A.scaler with I=(R-offset)/gain for charge/current cuts\n"
      << "  --bcm4a-gain <value>         Raw BCM4A gain for --use-raw-bcm4a (Default: 9597.0)\n"
      << "  --bcm4a-offset <value>       Raw BCM4A offset for --use-raw-bcm4a (Default: 177.3)\n"
      << "  --use-rolling-stability      Enable rolling stability window on TSHelH charge/read\n"
      << "  --stability-window <N>       Rolling window size (Default: 30)\n"
      << "  --stability-frac-range <f>   Max fractional charge range (Default: 0.08)\n"
      << "  --pid-cer-npe-min <value>    CER NPE threshold for PID efficiency (Default: 1.0; uses >)\n"
      << "  --pid-cal-etottracknorm-min <value>  CAL etotnorm lower cut for PID efficiency (Default: 0.6; uses >)\n"
      << "  --pid-cal-etottracknorm-max <value>  Optional CAL etotnorm upper cut (Default: disabled)\n"
      << "  --hodo-track-max-chisq <value>  Legacy-like hodo efficiency track chisq cut (Default: 10.0)\n"
      << "  --hodo-track-slop-cm <value>    Legacy-like hodo efficiency track slop in cm (Default: 2.0)\n"
      << "  --pid-no-track-required      Do not require H.dc.ntrack > 0 for PID denominator\n"
      << "  --no-progress                Disable progress output\n"
      << "  --runs <tokens...>           Required. Example: 4076 4077 4080-4083\n";
}

static bool parse_args(int argc, char **argv, Config &cfg) {
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--root-dir" && i + 1 < argc) {
      cfg.root_dir = argv[++i];
    } else if (arg == "--db" && i + 1 < argc) {
      cfg.db_path = argv[++i];
    } else if (arg == "--trigger-overrides" && i + 1 < argc) {
      cfg.trigger_overrides_path = argv[++i];
    } else if (arg == "--use-report-efficiencies") {
      cfg.use_report_efficiencies = true;
    } else if (arg == "--report-dir" && i + 1 < argc) {
      cfg.report_dir = argv[++i];
    } else if (arg == "--out-csv" && i + 1 < argc) {
      cfg.out_csv = argv[++i];
    } else if (arg == "--out-segment-csv" && i + 1 < argc) {
      cfg.out_segment_csv = argv[++i];
    } else if (arg == "--default-ps" && i + 1 < argc) {
      cfg.default_ps = std::stod(argv[++i]);
    } else if (arg == "--current-window" && i + 1 < argc) {
      cfg.current_window = std::stod(argv[++i]);
    } else if (arg == "--current-window-frac" && i + 1 < argc) {
      cfg.current_window_frac = std::stod(argv[++i]);
      cfg.use_fractional_window = true;
    } else if (arg == "--use-absolute-window") {
      cfg.use_fractional_window = false;
    } else if (arg == "--current-correction" && i + 1 < argc) {
      cfg.current_correction = std::stod(argv[++i]);
    } else if (arg == "--beam-on-threshold" && i + 1 < argc) {
      cfg.beam_on_threshold_uA = std::stod(argv[++i]);
    } else if (arg == "--use-raw-bcm4a") {
      cfg.use_raw_bcm4a = true;
    } else if (arg == "--bcm4a-gain" && i + 1 < argc) {
      cfg.bcm4a_gain = std::stod(argv[++i]);
    } else if (arg == "--bcm4a-offset" && i + 1 < argc) {
      cfg.bcm4a_offset = std::stod(argv[++i]);
    } else if (arg == "--use-rolling-stability") {
      cfg.use_rolling_stability = true;
    } else if (arg == "--stability-window" && i + 1 < argc) {
      cfg.stability_window = std::stoi(argv[++i]);
    } else if (arg == "--stability-frac-range" && i + 1 < argc) {
      cfg.stability_frac_range = std::stod(argv[++i]);
    } else if (arg == "--pid-cer-npe-min" && i + 1 < argc) {
      cfg.pid_cer_npe_min = std::stod(argv[++i]);
    } else if (arg == "--pid-cal-etottracknorm-min" && i + 1 < argc) {
      cfg.pid_cal_etotnorm_min = std::stod(argv[++i]);
    } else if (arg == "--pid-cal-etottracknorm-max" && i + 1 < argc) {
      cfg.pid_cal_etotnorm_max = std::stod(argv[++i]);
    } else if (arg == "--hodo-track-max-chisq" && i + 1 < argc) {
      cfg.hodo_track_max_chisq = std::stod(argv[++i]);
    } else if (arg == "--hodo-track-slop-cm" && i + 1 < argc) {
      cfg.hodo_track_slop_cm = std::stod(argv[++i]);
    } else if (arg == "--pid-no-track-required") {
      cfg.pid_require_track = false;
    } else if (arg == "--no-progress") {
      cfg.progress = false;
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
  auto db = parse_db(cfg.db_path);
  auto trigger_overrides = parse_trigger_overrides(cfg.trigger_overrides_path);

  std::vector<RunResult> results;
  std::vector<SegmentResult> segment_results;
  int found_runs = 0;
  for (int run : runs) {
    auto files = find_run_files(cfg.root_dir, run);
    if (files.empty()) {
      if (cfg.progress) {
        std::cout << "[run " << run << "] no files found in " << cfg.root_dir << "\n";
      }
      continue;
    }

    ++found_runs;
    if (cfg.progress) {
      std::cout << "[run " << run << "] starting with " << files.size() << " file(s)\n";
      std::cout.flush();
    }

    double expected_current = 10.0;
    double ps = cfg.default_ps;
    int trig_number = 4;
    std::string prescale_token;

    const bool have_median = compute_run_median_current(
      files,
      expected_current,
      cfg.use_raw_bcm4a,
      cfg.bcm4a_gain,
      cfg.bcm4a_offset,
      cfg.progress);

    auto it = db.find(run);
    if (!have_median && it != db.end()) {
      if (std::isfinite(it->second.ave_current)) expected_current = it->second.ave_current;
    }
    if (it != db.end()) {
      if (std::isfinite(it->second.ps) && it->second.ps > 0.0) ps = it->second.ps;
      if (it->second.which_trigger > 0) trig_number = it->second.which_trigger;
      prescale_token = it->second.prescale_token;
    }
    const auto ov_it = trigger_overrides.find(run);
    if (ov_it != trigger_overrides.end()) {
      if (std::isfinite(ov_it->second.ps) && ov_it->second.ps > 0.0) ps = ov_it->second.ps;
      if (ov_it->second.which_trigger > 0) trig_number = ov_it->second.which_trigger;
      if (!ov_it->second.prescale_token.empty()) prescale_token = ov_it->second.prescale_token;
    }

    if (cfg.progress) {
      double window_lo = std::numeric_limits<double>::quiet_NaN();
      double window_hi = std::numeric_limits<double>::quiet_NaN();
      compute_window_limits(expected_current,
                            cfg.current_window,
                            cfg.current_window_frac,
                            cfg.use_fractional_window,
                            window_lo,
                            window_hi);
      std::cout << "[run " << run << "] current-window center = "
                << std::fixed << std::setprecision(3) << expected_current
                << " uA (" << (have_median ? "median" : "fallback") << ")\n";
      if (ov_it != trigger_overrides.end()) {
        std::cout << "[run " << run << "] using trigger override: trig="
                  << trig_number << ", ps=" << ps
                  << ", token=" << prescale_token << "\n";
      }
      if (cfg.use_fractional_window) {
        std::cout << "[run " << run << "] fractional window = ±"
                  << std::fixed << std::setprecision(3) << (cfg.current_window_frac * 100.0)
                  << "% => [" << window_lo << ", " << window_hi << "] uA\n";
      } else {
        std::cout << "[run " << run << "] absolute window = ±"
                  << std::fixed << std::setprecision(3) << cfg.current_window
                  << " uA => [" << window_lo << ", " << window_hi << "] uA\n";
      }
      std::cout.flush();
    }

    double run_hel_hm = 0.0;
    double run_hel_hp = 0.0;
    double run_hel0 = 0.0;
    int run_hel_missing = 0;

    for (const auto &path : files) {
      SegmentResult seg = analyze_segment_charge(
          run,
          path,
          expected_current,
          cfg.current_window,
          cfg.current_window_frac,
          cfg.use_fractional_window,
          cfg.current_correction,
          cfg.use_raw_bcm4a,
          cfg.bcm4a_gain,
          cfg.bcm4a_offset);

      HelicitySelectionResult hel = analyze_segment_helicity(path, expected_current, cfg);
      seg.hel_charge_hm_uC = hel.charge_hm_uC;
      seg.hel_charge_hp_uC = hel.charge_hp_uC;
      seg.hel_charge_total_uC = hel.charge_total_uC;
      seg.hel0_charge_uC = hel.hel0_charge_uC;
      seg.charge_tshelh_uC = hel.charge_total_uC + hel.hel0_charge_uC;
      seg.good_evnum_ranges = hel.good_evnum_ranges;
      seg.hel0_evnum_ranges = hel.hel0_evnum_ranges;
      seg.hel_missing_branches = hel.missing_branches;

      run_hel_hm += hel.charge_hm_uC;
      run_hel_hp += hel.charge_hp_uC;
      run_hel0 += hel.hel0_charge_uC;
      run_hel_missing += hel.missing_branches;

      segment_results.push_back(seg);
    }

    results.push_back(analyze_run(
      run,
      files,
      cfg.report_dir,
      cfg.use_report_efficiencies,
      expected_current,
      ps,
      cfg.current_window,
      cfg.current_window_frac,
      cfg.use_fractional_window,
      cfg.current_correction,
      cfg.beam_on_threshold_uA,
      cfg.use_raw_bcm4a,
      cfg.bcm4a_gain,
      cfg.bcm4a_offset,
      cfg.pid_cer_npe_min,
      cfg.pid_cal_etotnorm_min,
      cfg.pid_cal_etotnorm_max,
      cfg.hodo_track_max_chisq,
      cfg.hodo_track_slop_cm,
      cfg.pid_require_track,
      cfg.use_rolling_stability,
      cfg.stability_window,
      cfg.stability_frac_range,
      trig_number,
      prescale_token,
      cfg.progress));
    results.back().hel_charge_hm_uC = run_hel_hm;
    results.back().hel_charge_hp_uC = run_hel_hp;
    results.back().hel_charge_total_uC = run_hel_hm + run_hel_hp;
    results.back().hel0_charge_uC = run_hel0;
    results.back().charge_tshelh_uC = results.back().hel_charge_total_uC + results.back().hel0_charge_uC;
    results.back().hel_missing_branches = run_hel_missing;
    if (it != db.end() && std::isfinite(it->second.charge_tot)) {
      results.back().Charge_yaopeng = it->second.charge_tot;
    }

    if (cfg.progress) {
      std::cout << "[run " << run << "] finished\n";
      std::cout.flush();
    }
  }

  if (cfg.progress) {
    std::cout << "Processed " << found_runs << " run(s) with available files. Writing TSV...\n";
    std::cout.flush();
  }

  std::ofstream out(cfg.out_csv);
  if (!out) {
    std::cerr << "Failed to open output TSV: " << cfg.out_csv << "\n";
    return 2;
  }

          out << "run\tn_files\texpected_current_uA\tavg_current_uA\tCharge_yaopeng\t"
            << "hel_charge_hm_uC\thel_charge_hp_uC\thel_charge_total_uC\t"
            << "hel0_charge_uC\tcharge_tshelh_uC\tcharge_uC\tcharge_above_threshold_uC\tcharge_above_3p0_uC\t"
      << "scaler_htrig4_total\tscaler_edtm_total\ttrig_accp_total\t"
      << "scaler_edtm_no_cut\ttrig_accp_no_cut\tbeam_on_percent_edtm\t"
      << "beam_on_percent_trig_accp\tt_edtm_accepted\tt_htrig4_all_accepted\t"
      << "t_htrig4_phy_accepted\tclta_livetime_tsh\tclta_livetime_tdc\t"
      << "cltp_livetime_tdc\ttlt_livetime_edtm\t"
      << "tlt_livetime_edtm_sym_beamcut\t"
      << "tlt_livetime_edtm_nocut_den\t"
      << "tlt_livetime_edtm_ratio_sym\t"
      << "tlt_livetime_edtm_ratio_nocut\t"
      << "NewGEN_nps_coin_livetime\tNewGEN_nps_coin_livetime_trig6\t"
      << "NewGEN_electronic_livetime\tNewGEN_electronic_livetime_evt\t"
      << "NewGEN_electronic_livetime_evt_noedtmden\t"
      << "HMS_pid_den_events\tHMS_pid_eff\tHMS_cer_eff\tHMS_cal_eff\t"
      << "HMS_cer_eff_tag_cal\tHMS_cal_eff_tag_cer\t"
      << "HMS_tracking_should_events\tHMS_tracking_did_events\t"
      << "HMS_tracking_eff\tHMS_tracking_eff_err\t"
      << "HMS_hodo_3of4_den_events\tHMS_hodo_3of4_eff\t"
      << "HMS_hodo_3of4_nhits_den_events\tHMS_hodo_3of4_nhits_eff\t"
      << "OG_pretrig_livetime_3_4\tOG_pretrig_livetime_3_4_err\t"
      << "beam_time\tprescale_token\tps_factor\twhich_TRIG\t"
      << "negative_dt_intervals\tnon_monotonic_scaler_steps\t"
      << "suspicious_scaler_jump_steps\tmissing_branches\thel_missing_branches\n";

  out << std::fixed << std::setprecision(6);
  for (const auto &r : results) {
    out << r.run << '\t'
      << r.n_files << '\t'
      << r.expected_current_uA << '\t'
      << r.avg_current_uA << '\t'
      << r.Charge_yaopeng << '\t'
      << r.hel_charge_hm_uC << '\t'
      << r.hel_charge_hp_uC << '\t'
      << r.hel_charge_total_uC << '\t'
      << r.hel0_charge_uC << '\t'
      << r.charge_tshelh_uC << '\t'
      << r.charge_uC << '\t'
      << r.charge_above_threshold_uC << '\t'
      << r.charge_above_3p0_uC << '\t'
      << r.scaler_htrig4_total << '\t'
      << r.scaler_edtm_total << '\t'
      << r.trig_accp_total << '\t'
      << r.scaler_edtm_no_cut << '\t'
      << r.trig_accp_no_cut << '\t'
      << r.beam_on_percent_edtm << '\t'
      << r.beam_on_percent_trig_accp << '\t'
      << r.t_edtm_accepted << '\t'
      << r.t_htrig4_all_accepted << '\t'
      << r.t_htrig4_phy_accepted << '\t'
      << r.clta_livetime_tsh << '\t'
      << r.clta_livetime_tdc << '\t'
      << r.cltp_livetime_tdc << '\t'
      << r.tlt_livetime_edtm << '\t'
      << r.tlt_livetime_edtm_sym_beamcut << '\t'
      << r.tlt_livetime_edtm_nocut_den << '\t'
      << r.tlt_livetime_edtm_ratio_sym << '\t'
      << r.tlt_livetime_edtm_ratio_nocut << '\t'
        << r.NewGEN_nps_coin_livetime << '\t'
        << r.NewGEN_nps_coin_livetime_trig6 << '\t'
      << r.NewGEN_electronic_livetime << '\t'
      << r.NewGEN_electronic_livetime_evt << '\t'
      << r.NewGEN_electronic_livetime_evt_noedtmden << '\t'
      << r.HMS_pid_den_events << '\t'
      << r.HMS_pid_eff << '\t'
      << r.HMS_cer_eff << '\t'
      << r.HMS_cal_eff << '\t'
      << r.HMS_cer_eff_tag_cal << '\t'
      << r.HMS_cal_eff_tag_cer << '\t'
      << r.HMS_tracking_should_events << '\t'
      << r.HMS_tracking_did_events << '\t'
      << r.HMS_tracking_eff << '\t'
      << r.HMS_tracking_eff_err << '\t'
      << r.HMS_hodo_3of4_den_events << '\t'
      << r.HMS_hodo_3of4_eff << '\t'
      << r.HMS_hodo_3of4_nhits_den_events << '\t'
      << r.HMS_hodo_3of4_nhits_eff << '\t'
        << r.OG_pretrig_livetime_3_4 << '\t'
        << r.OG_pretrig_livetime_3_4_err << '\t'
        << r.beam_time << '\t'
        << r.prescale_token << '\t'
        << r.ps_factor << '\t'
        << r.which_TRIG << '\t'
      << r.negative_dt_intervals << '\t'
      << r.non_monotonic_scaler_steps << '\t'
      << r.suspicious_scaler_jump_steps << '\t'
      << r.missing_branches << '\t'
      << r.hel_missing_branches << '\n';
  }

  std::cout << "Wrote " << cfg.out_csv << " with " << results.size() << " runs\n";

  std::ofstream out_seg(cfg.out_segment_csv);
  if (!out_seg) {
        std::cerr << "Failed to open segment output TSV: " << cfg.out_segment_csv << "\n";
    return 3;
  }

      out_seg << "run\tsegment\tfile\texpected_current_uA\tavg_current_uA\tcharge_uC\t"
        << "hel_charge_hm_uC\thel_charge_hp_uC\thel_charge_total_uC\t"
        << "hel0_charge_uC\tcharge_tshelh_uC\t"
        << "good_evnum_ranges\thel0_evnum_ranges\t"
        << "negative_dt_intervals\tnon_monotonic_scaler_steps\tmissing_branches\thel_missing_branches\n";
  out_seg << std::fixed << std::setprecision(6);
  for (const auto &s : segment_results) {
        out_seg << s.run << '\t'
          << s.segment << '\t'
          << s.file << '\t'
          << s.expected_current_uA << '\t'
          << s.avg_current_uA << '\t'
          << s.charge_uC << '\t'
          << s.hel_charge_hm_uC << '\t'
          << s.hel_charge_hp_uC << '\t'
          << s.hel_charge_total_uC << '\t'
          << s.hel0_charge_uC << '\t'
          << s.charge_tshelh_uC << '\t'
          << s.good_evnum_ranges << '\t'
          << s.hel0_evnum_ranges << '\t'
          << s.negative_dt_intervals << '\t'
          << s.non_monotonic_scaler_steps << '\t'
          << s.missing_branches << '\t'
          << s.hel_missing_branches << '\n';
  }

  std::cout << "Wrote " << cfg.out_segment_csv << " with " << segment_results.size() << " segments\n";
  return 0;
}
