// ============================================================================
// final_compute_luminosity_standalone.cxx
//
// Standalone Hall C NPS efficiency/luminosity calculator.
// Produces paired outputs:
//   1) good-event selected calculation
//   2) current-window-only calculation without NPS_good_events selection
//
// Build after loading the Hall C/NPS ROOT environment:
//   g++ -O2 -std=c++17 final_compute_luminosity_standalone.cxx \
//      -o final_compute_luminosity_standalone `root-config --cflags --libs`
// ============================================================================

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
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
#include <vector>

#include "NPS_selection_Christine_helper.h"

namespace fs = std::filesystem;

struct Config {
  std::string root_dir = ".";
  std::string prescale_csv =
      "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/config/nps_dvcs_all_kins_main.csv";
  std::string out_csv = "final_luminosity_results_good_events.csv";
  std::string out_csv_no_good;
  std::vector<std::string> run_tokens;

  double default_ps = 1.0;
  double window_frac = 0.15;
  double current_window_uA = 1.5;
  bool use_absolute_window = false;
  bool progress = true;

  std::string good_helicity = "quartet_pm";
  bool good_use_current_window = true;
  double good_i0_fixed_uA = -1.0;
  int stable_window_N = 30;
  double stability_frac_range = 0.08;

  double pid_cer_npe_min = 1.0;
  double pid_cal_etotnorm_min = 0.6;
  double pid_cal_etotnorm_max = std::numeric_limits<double>::infinity();
  bool pid_require_track = true;
  bool show_help = false;
};

struct HmsPidAccum {
  long long den = 0;
  long long both = 0;
  long long cer = 0;
  long long cal = 0;
  long long cer_tag_den = 0;
  long long cal_tag_den = 0;
};

struct TrackingAccum {
  long long should = 0;
  long long did = 0;
};

struct HodoAccum {
  long long den = 0;
  long long num = 0;
};

struct ScalerAccum {
  double charge_uC = 0.0;
  double current_weighted = 0.0;
  double beam_time = 0.0;
  double total_valid_time = 0.0;

  double scaler_htrig_total = 0.0;
  double scaler_edtm_total = 0.0;
  double trig_accp_total = 0.0;
  double scaler_edtm_no_cut = 0.0;
  double trig_accp_no_cut = 0.0;

  double s1x_scaler_total = 0.0;
  double s1x_rate_sample_sum = 0.0;
  double s1x_rate_sample_sum2 = 0.0;
  long long s1x_rate_samples = 0;

  int negative_dt_intervals = 0;
  int non_monotonic_scaler_steps = 0;
  int suspicious_scaler_jump_steps = 0;
  int missing_branch_segments = 0;
  int scaler_evcount_gate_segments = 0;
  std::set<std::string> missing_branches;
};

struct HelChargeAccum {
  double hel0_before = 0.0;
  double hel_neg_before = 0.0;
  double hel_pos_before = 0.0;
  double hel0_after = 0.0;
  double hel_neg_after = 0.0;
  double hel_pos_after = 0.0;
  int missing_branch_segments = 0;
};

struct EventAccum {
  std::map<int, long long> edtm_hist;
  long long tdc_edtm_accepted = 0;
  long long tdc_htrig_all_accepted = 0;
  long long tdc_htrig_phy_accepted = 0;
  long long selected_events = 0;
  long long current_window_events = 0;
  long long missing_branch_segments = 0;
  HmsPidAccum pid;
  TrackingAccum tracking;
  HodoAccum hodo;
};

struct SegmentSelection {
  bool ok = false;
  bool apply_good_events = false;
  std::string message;
  std::vector<RangeI> evcount_ranges;
  std::vector<RangeI> gevnum_ranges;
  double i0_used_uA = std::numeric_limits<double>::quiet_NaN();
  double current_min_uA = std::numeric_limits<double>::quiet_NaN();
  double current_max_uA = std::numeric_limits<double>::quiet_NaN();
  double mean_current_uA = std::numeric_limits<double>::quiet_NaN();
  int n_scaler_reads_pre = 0;
  int n_scaler_reads_post = 0;
};

struct RunResult {
  int run = 0;
  int n_files = 0;
  int n_segments_used = 0;
  int selection_failed_segments = 0;
  bool good_event_selection = false;

  double expected_current_uA = std::numeric_limits<double>::quiet_NaN();
  double selection_i0_uA = std::numeric_limits<double>::quiet_NaN();
  double selection_current_min_uA = std::numeric_limits<double>::quiet_NaN();
  double selection_current_max_uA = std::numeric_limits<double>::quiet_NaN();
  int selection_scaler_reads_pre = 0;
  int selection_scaler_reads_post = 0;

  double avg_current_uA = std::numeric_limits<double>::quiet_NaN();
  double charge_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_neg_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_pos_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel0_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_neg_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_pos_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel0_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();

  double scaler_htrig_total = 0.0;
  double scaler_edtm_total = 0.0;
  double trig_accp_total = 0.0;
  double scaler_edtm_no_cut = 0.0;
  double trig_accp_no_cut = 0.0;
  double beam_time = 0.0;
  double beam_on_frac = std::numeric_limits<double>::quiet_NaN();
  double beam_on_percent_edtm = std::numeric_limits<double>::quiet_NaN();
  double beam_on_percent_trig_accp = std::numeric_limits<double>::quiet_NaN();

  double htrig_rate_hz = std::numeric_limits<double>::quiet_NaN();
  double edtm_rate_hz = std::numeric_limits<double>::quiet_NaN();
  double trig_accp_rate_hz = std::numeric_limits<double>::quiet_NaN();
  double s1x_rate_hz_from_scaler = std::numeric_limits<double>::quiet_NaN();
  double s1x_scalerRate_mean_hz = std::numeric_limits<double>::quiet_NaN();
  double s1x_scalerRate_rms_hz = std::numeric_limits<double>::quiet_NaN();

  long long tdc_edtm_accepted = 0;
  long long tdc_htrig_all_accepted = 0;
  long long tdc_htrig_phy_accepted = 0;
  long long selected_events = 0;
  long long current_window_events = 0;
  double tdc_edtm_rate_hz = std::numeric_limits<double>::quiet_NaN();
  double tdc_htrig_all_rate_hz = std::numeric_limits<double>::quiet_NaN();
  double tdc_htrig_phy_rate_hz = std::numeric_limits<double>::quiet_NaN();

  double NewGen_EDTM_livetime = std::numeric_limits<double>::quiet_NaN();
  double clta_livetime_tsh = std::numeric_limits<double>::quiet_NaN();
  double clta_livetime_tdc = std::numeric_limits<double>::quiet_NaN();
  double cltp_livetime_tdc = std::numeric_limits<double>::quiet_NaN();

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

  std::string prescale_token;
  double ps_factor = 1.0;
  std::string which_TRIG;
  int negative_dt_intervals = 0;
  int non_monotonic_scaler_steps = 0;
  int suspicious_scaler_jump_steps = 0;
  int missing_branch_segments = 0;
  int hel_missing_branch_segments = 0;
  int scaler_evcount_gate_segments = 0;
  std::string missing_branches;
};

static std::string trim_copy(std::string s) {
  const size_t first = s.find_first_not_of(" \t\n\r");
  if (first == std::string::npos) return "";
  const size_t last = s.find_last_not_of(" \t\n\r");
  return s.substr(first, last - first + 1);
}

static std::string csv_quote(const std::string &s) {
  std::string out = "\"";
  for (char c : s) {
    if (c == '"') out += "\"\"";
    else out += c;
  }
  out += "\"";
  return out;
}

static double safe_div(double n, double d) {
  if (!std::isfinite(n) || !std::isfinite(d) || d == 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return n / d;
}

static double clamp01(double x) {
  if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
  return std::max(0.0, std::min(1.0, x));
}

static double binomial_err(long long num, long long den) {
  if (den <= 0) return std::numeric_limits<double>::quiet_NaN();
  const double p = safe_div(static_cast<double>(num), static_cast<double>(den));
  return std::sqrt(std::max(0.0, p * (1.0 - p)) / static_cast<double>(den));
}

static bool has_branch(TTree *tree, const char *name) {
  return tree != nullptr && tree->GetBranch(name) != nullptr;
}

static bool in_ranges(long long value, const std::vector<RangeI> &ranges) {
  if (ranges.empty()) return true;
  for (const auto &r : ranges) {
    if (value >= r.lo && value <= r.hi) return true;
  }
  return false;
}

static std::string ranges_to_string(const std::vector<RangeI> &ranges) {
  std::ostringstream out;
  bool first = true;
  for (const auto &r : ranges) {
    if (!first) out << ";";
    first = false;
    out << r.lo << "-" << r.hi;
  }
  return out.str();
}

static std::string with_suffix_before_extension(const std::string &path,
                                                const std::string &suffix) {
  fs::path p(path);
  const std::string stem = p.stem().string();
  const std::string ext = p.extension().string();
  fs::path out = p.parent_path() / (stem + suffix + ext);
  return out.string();
}

static int extract_prescale_r(const std::string &token_in) {
  const std::string token = trim_copy(token_in);
  const std::regex re("ps\\d*\\s*=\\s*(\\d+)", std::regex_constants::icase);
  std::smatch m;
  if (std::regex_search(token, m, re) && m.size() > 1) {
    try { return std::stoi(m[1].str()); } catch (...) { return std::numeric_limits<int>::min(); }
  }
  const size_t eq = token.find('=');
  if (eq != std::string::npos) {
    try { return std::stoi(trim_copy(token.substr(eq + 1))); }
    catch (...) { return std::numeric_limits<int>::min(); }
  }
  return std::numeric_limits<int>::min();
}

static int extract_trig_number(const std::string &token_in) {
  const std::string token = trim_copy(token_in);
  const std::regex re("ps(\\d+)", std::regex_constants::icase);
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

static std::map<int, std::string> parse_prescale_csv(const std::string &csv_path) {
  std::map<int, std::string> out;
  std::ifstream in(csv_path);
  if (!in) return out;

  std::string header_line;
  if (!std::getline(in, header_line)) return out;
  std::vector<std::string> headers;
  std::istringstream hss(header_line);
  std::string h;
  while (std::getline(hss, h, ',')) headers.push_back(trim_copy(h));

  int run_col = -1, ps_col = -1;
  for (size_t i = 0; i < headers.size(); ++i) {
    std::string lower = headers[i];
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower.find("run_number") != std::string::npos || lower == "run") run_col = static_cast<int>(i);
    if (lower.find("prescale") != std::string::npos) ps_col = static_cast<int>(i);
  }
  if (run_col < 0 || ps_col < 0) return out;

  std::string line;
  while (std::getline(in, line)) {
    if (trim_copy(line).empty()) continue;
    std::istringstream lss(line);
    std::vector<std::string> fields;
    std::string f;
    while (std::getline(lss, f, ',')) fields.push_back(trim_copy(f));
    if (static_cast<int>(fields.size()) <= std::max(run_col, ps_col)) continue;
    try {
      out[std::stoi(fields[run_col])] = fields[ps_col];
    } catch (...) {
      continue;
    }
  }
  return out;
}

static std::vector<int> parse_run_list(const std::vector<std::string> &tokens) {
  std::set<int> unique;
  for (const auto &token_raw : tokens) {
    std::string token = token_raw;
    std::replace(token.begin(), token.end(), ',', ' ');
    std::istringstream iss(token);
    std::string part;
    while (iss >> part) {
      const size_t dash = part.find('-');
      if (dash != std::string::npos) {
        int lo = std::stoi(part.substr(0, dash));
        int hi = std::stoi(part.substr(dash + 1));
        if (hi < lo) std::swap(lo, hi);
        for (int r = lo; r <= hi; ++r) unique.insert(r);
      } else {
        unique.insert(std::stoi(part));
      }
    }
  }
  return std::vector<int>(unique.begin(), unique.end());
}

static std::map<int, std::vector<std::string>> index_run_files(
    const std::string &root_dir,
    const std::set<int> &target_runs) {
  std::map<int, std::vector<std::string>> indexed;
  if (!fs::exists(root_dir)) return indexed;

  for (const auto &entry : fs::directory_iterator(root_dir)) {
    if (!entry.is_regular_file()) continue;
    const std::string name = entry.path().filename().string();
    int run = 0;
    bool parsed = false;

    if (name.rfind("nps_hms_coin_", 0) == 0 && name.find(".root") != std::string::npos) {
      const size_t p1 = std::string("nps_hms_coin_").size();
      const size_t p2 = name.find('_', p1);
      if (p2 != std::string::npos) {
        try { run = std::stoi(name.substr(p1, p2 - p1)); parsed = true; } catch (...) {}
      }
    } else if (name.rfind("skim_run", 0) == 0 && name.find(".root") != std::string::npos) {
      const size_t p1 = std::string("skim_run").size();
      const size_t p2 = name.find(".root");
      if (p2 > p1) {
        try { run = std::stoi(name.substr(p1, p2 - p1)); parsed = true; } catch (...) {}
      }
    }

    if (parsed && target_runs.count(run)) indexed[run].push_back(entry.path().string());
  }
  for (auto &kv : indexed) std::sort(kv.second.begin(), kv.second.end());
  return indexed;
}

static double median(std::vector<double> values) {
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  const size_t mid = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + mid, values.end());
  const double hi = values[mid];
  if (values.size() % 2 == 1) return hi;
  std::nth_element(values.begin(), values.begin() + mid - 1, values.end());
  return 0.5 * (values[mid - 1] + hi);
}

static double compute_expected_current(const std::vector<std::string> &files) {
  std::vector<double> currents;
  for (const auto &path : files) {
    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) continue;
    auto *tsh = dynamic_cast<TTree *>(file.Get("TSH"));
    if (!has_branch(tsh, "H.BCM4A.scalerCurrent")) continue;

    double current = 0.0;
    tsh->SetBranchStatus("*", 0);
    tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
    tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &current);
    const Long64_t n = tsh->GetEntries();
    currents.reserve(currents.size() + static_cast<size_t>(std::min<Long64_t>(n, 100000)));
    for (Long64_t i = 0; i < n; ++i) {
      tsh->GetEntry(i);
      if (std::isfinite(current) && current > 1.0) currents.push_back(current);
    }
  }
  return median(std::move(currents));
}

static HelicityMode parse_helicity_mode_or_die(const std::string &mode) {
  if (mode == "quartet_pm") return HelicityMode::QuartetPM;
  if (mode == "quartet_all") return HelicityMode::QuartetAll;
  if (mode == "ignore") return HelicityMode::IgnoreHelicity;
  throw std::runtime_error("Unknown --good-helicity mode: " + mode);
}

static SelectionSettings make_selection_settings(const Config &cfg,
                                                 double expected_current) {
  SelectionSettings sel;
  sel.helicity_mode = parse_helicity_mode_or_die(cfg.good_helicity);
  sel.junk_floor_uA = 2.5;
  sel.mean_current_calc_min = 2.5;
  sel.mean_current_min = 2.75;
  sel.use_current_window = cfg.good_use_current_window;
  if (cfg.use_absolute_window) {
    const double center = (cfg.good_i0_fixed_uA > 0.0) ? cfg.good_i0_fixed_uA : expected_current;
    sel.i0_mode = I0Mode::Fixed;
    sel.i0_fixed_uA = center;
    sel.window_frac = (std::isfinite(center) && center > 0.0)
        ? cfg.current_window_uA / center
        : cfg.window_frac;
  } else {
    sel.i0_mode = (cfg.good_i0_fixed_uA > 0.0) ? I0Mode::Fixed : I0Mode::Peak;
    sel.i0_fixed_uA = cfg.good_i0_fixed_uA;
    sel.window_frac = cfg.window_frac;
  }
  sel.stable_window_N = cfg.stable_window_N;
  sel.max_charge_frac_spread = cfg.stability_frac_range;
  sel.branch_scaler_current = "H.BCM4A_Hel.scalerCurrent";
  sel.branch_scaler_charge = "H.BCM4A_Hel.scalerCharge";
  sel.branch_helicity = "actualHelicity";
  sel.branch_evcount = "evcount";
  sel.branch_evNumber = "evNumber";
  sel.branch_evnum_T = "g.evnum";
  return sel;
}

static SegmentSelection build_segment_selection(const std::string &path,
                                                const Config &cfg,
                                                bool apply_good_events,
                                                double expected_current) {
  SegmentSelection sel;
  sel.apply_good_events = apply_good_events;
  if (!apply_good_events) {
    sel.ok = true;
    sel.message = "good_event_selection_disabled";
    sel.i0_used_uA = expected_current;
    if (cfg.use_absolute_window) {
      sel.current_min_uA = expected_current - cfg.current_window_uA;
      sel.current_max_uA = expected_current + cfg.current_window_uA;
    } else {
      sel.current_min_uA = (1.0 - cfg.window_frac) * expected_current;
      sel.current_max_uA = (1.0 + cfg.window_frac) * expected_current;
    }
    return sel;
  }

  try {
    SelectionResult pick = build_event_selection(path, make_selection_settings(cfg, expected_current));
    sel.ok = pick.ok;
    sel.message = pick.message;
    sel.evcount_ranges = pick.evcount_ranges;
    sel.gevnum_ranges = pick.gevnum_ranges;
    sel.i0_used_uA = pick.i0_used_uA;
    sel.current_min_uA = pick.current_min_uA;
    sel.current_max_uA = pick.current_max_uA;
    sel.mean_current_uA = pick.mean_current_uA;
    sel.n_scaler_reads_pre = pick.n_scaler_reads_pre;
    sel.n_scaler_reads_post = pick.n_scaler_reads_post;
  } catch (const std::exception &ex) {
    sel.ok = false;
    sel.message = ex.what();
  }
  return sel;
}

static bool current_in_mode_window(double current,
                                   const SegmentSelection &selection,
                                   double expected_current,
                                   double window_frac) {
  double lo = selection.current_min_uA;
  double hi = selection.current_max_uA;
  if (!std::isfinite(lo) || !std::isfinite(hi)) {
    lo = (1.0 - window_frac) * expected_current;
    hi = (1.0 + window_frac) * expected_current;
  }
  if (!std::isfinite(lo) || !std::isfinite(hi)) return true;
  return current >= lo && current <= hi;
}

static void accumulate_hel_charge(const std::string &path,
                                  const SegmentSelection &selection,
                                  HelChargeAccum &acc) {
  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) {
    ++acc.missing_branch_segments;
    return;
  }
  auto *tsh = dynamic_cast<TTree *>(file.Get("TSHelH"));
  if (!has_branch(tsh, "actualHelicity") ||
      !has_branch(tsh, "H.BCM4A_Hel.scalerCharge") ||
      !has_branch(tsh, "evcount")) {
    ++acc.missing_branch_segments;
    return;
  }

  double helicity = 0.0, charge = 0.0, evcount = 0.0;
  tsh->SetBranchStatus("*", 0);
  tsh->SetBranchStatus("actualHelicity", 1);
  tsh->SetBranchStatus("H.BCM4A_Hel.scalerCharge", 1);
  tsh->SetBranchStatus("evcount", 1);
  tsh->SetBranchAddress("actualHelicity", &helicity);
  tsh->SetBranchAddress("H.BCM4A_Hel.scalerCharge", &charge);
  tsh->SetBranchAddress("evcount", &evcount);

  const Long64_t n = tsh->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    tsh->GetEntry(i);
    const long long hel = std::llround(helicity);
    if (hel == 0) acc.hel0_before += charge;
    else if (hel < 0) acc.hel_neg_before += charge;
    else acc.hel_pos_before += charge;

    if (selection.apply_good_events &&
        !in_ranges(std::llround(evcount), selection.evcount_ranges)) {
      continue;
    }
    if (hel == 0) acc.hel0_after += charge;
    else if (hel < 0) acc.hel_neg_after += charge;
    else acc.hel_pos_after += charge;
  }
}

static void accumulate_scalers(const std::string &path,
                               int trig_num,
                               double expected_current,
                               const Config &cfg,
                               const SegmentSelection &selection,
                               ScalerAccum &acc) {
  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) {
    ++acc.missing_branch_segments;
    return;
  }
  auto *tsh = dynamic_cast<TTree *>(file.Get("TSH"));
  if (!tsh) {
    ++acc.missing_branch_segments;
    acc.missing_branches.insert("TSH");
    return;
  }

  const std::string htrig_name = "H.hTRIG" + std::to_string(trig_num) + ".scaler";
  const bool has_core =
      has_branch(tsh, "H.1MHz.scalerTime") &&
      has_branch(tsh, "H.BCM4A.scalerCurrent") &&
      has_branch(tsh, htrig_name.c_str()) &&
      has_branch(tsh, "H.EDTM.scaler") &&
      has_branch(tsh, "H.hL1ACCP.scaler");
  if (!has_core) {
    ++acc.missing_branch_segments;
    for (const char *b : {"H.1MHz.scalerTime", "H.BCM4A.scalerCurrent", "H.EDTM.scaler", "H.hL1ACCP.scaler"}) {
      if (!has_branch(tsh, b)) acc.missing_branches.insert(b);
    }
    if (!has_branch(tsh, htrig_name.c_str())) acc.missing_branches.insert(htrig_name);
    return;
  }

  const bool has_evcount = has_branch(tsh, "evcount");
  const bool has_s1x_scaler = has_branch(tsh, "H.S1X.scaler");
  const bool has_s1x_rate = has_branch(tsh, "H.S1X.scalerRate");
  if (!has_s1x_scaler) acc.missing_branches.insert("H.S1X.scaler");
  if (!has_s1x_rate) acc.missing_branches.insert("H.S1X.scalerRate");

  double time = 0.0, current = 0.0, htrig = 0.0, edtm = 0.0, accp = 0.0;
  double evcount = 0.0, s1x_scaler = 0.0, s1x_rate = 0.0;

  tsh->SetBranchStatus("*", 0);
  tsh->SetBranchStatus("H.1MHz.scalerTime", 1);
  tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
  tsh->SetBranchStatus(htrig_name.c_str(), 1);
  tsh->SetBranchStatus("H.EDTM.scaler", 1);
  tsh->SetBranchStatus("H.hL1ACCP.scaler", 1);
  if (has_evcount) tsh->SetBranchStatus("evcount", 1);
  if (has_s1x_scaler) tsh->SetBranchStatus("H.S1X.scaler", 1);
  if (has_s1x_rate) tsh->SetBranchStatus("H.S1X.scalerRate", 1);

  tsh->SetBranchAddress("H.1MHz.scalerTime", &time);
  tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &current);
  tsh->SetBranchAddress(htrig_name.c_str(), &htrig);
  tsh->SetBranchAddress("H.EDTM.scaler", &edtm);
  tsh->SetBranchAddress("H.hL1ACCP.scaler", &accp);
  if (has_evcount) tsh->SetBranchAddress("evcount", &evcount);
  if (has_s1x_scaler) tsh->SetBranchAddress("H.S1X.scaler", &s1x_scaler);
  if (has_s1x_rate) tsh->SetBranchAddress("H.S1X.scalerRate", &s1x_rate);

  const Long64_t n = tsh->GetEntries();
  if (n < 2) return;
  tsh->GetEntry(0);

  double prev_time = time;
  double prev_current = current;
  double prev_htrig = htrig;
  double prev_edtm = edtm;
  double prev_accp = accp;
  double prev_s1x_scaler = s1x_scaler;

  if (selection.apply_good_events && has_evcount) ++acc.scaler_evcount_gate_segments;

  for (Long64_t i = 1; i < n; ++i) {
    tsh->GetEntry(i);
    const double dt = time - prev_time;
    if (!std::isfinite(dt) || dt <= 0.0) {
      ++acc.negative_dt_intervals;
      prev_time = time;
      prev_current = current;
      prev_htrig = htrig;
      prev_edtm = edtm;
      prev_accp = accp;
      prev_s1x_scaler = s1x_scaler;
      continue;
    }

    double dhtrig = htrig - prev_htrig;
    double dedtm = edtm - prev_edtm;
    double daccp = accp - prev_accp;
    double ds1x = has_s1x_scaler ? (s1x_scaler - prev_s1x_scaler) : 0.0;
    if (dhtrig < 0.0 || dedtm < 0.0 || daccp < 0.0 || ds1x < 0.0) {
      ++acc.non_monotonic_scaler_steps;
    }
    dhtrig = std::max(0.0, dhtrig);
    dedtm = std::max(0.0, dedtm);
    daccp = std::max(0.0, daccp);
    ds1x = std::max(0.0, ds1x);

    const double rate_htrig = dhtrig / dt;
    const double rate_edtm = dedtm / dt;
    const double rate_accp = daccp / dt;
    const bool suspicious =
        (dhtrig > 1.0e8) || (dedtm > 1.0e8) || (daccp > 1.0e8) ||
        (rate_htrig > 1.0e7) || (rate_edtm > 1.0e7) || (rate_accp > 1.0e7);
    if (suspicious) {
      ++acc.suspicious_scaler_jump_steps;
      prev_time = time;
      prev_current = current;
      prev_htrig = htrig;
      prev_edtm = edtm;
      prev_accp = accp;
      prev_s1x_scaler = s1x_scaler;
      continue;
    }

    acc.total_valid_time += dt;
    acc.scaler_edtm_no_cut += dedtm;
    acc.trig_accp_no_cut += daccp;

    const double avg_current = 0.5 * (prev_current + current);
    bool accepted = current_in_mode_window(avg_current, selection, expected_current, cfg.window_frac);
    if (selection.apply_good_events && has_evcount) {
      accepted = accepted && in_ranges(std::llround(evcount), selection.evcount_ranges);
    }

    if (accepted) {
      acc.charge_uC += std::max(0.0, avg_current) * dt;
      acc.current_weighted += std::max(0.0, avg_current) * dt;
      acc.beam_time += dt;
      acc.scaler_htrig_total += dhtrig;
      acc.scaler_edtm_total += dedtm;
      acc.trig_accp_total += daccp;
      acc.s1x_scaler_total += ds1x;
      if (has_s1x_rate && std::isfinite(s1x_rate) && s1x_rate >= 0.0) {
        acc.s1x_rate_sample_sum += s1x_rate;
        acc.s1x_rate_sample_sum2 += s1x_rate * s1x_rate;
        ++acc.s1x_rate_samples;
      }
    }

    prev_time = time;
    prev_current = current;
    prev_htrig = htrig;
    prev_edtm = edtm;
    prev_accp = accp;
    prev_s1x_scaler = s1x_scaler;
  }
}

static bool hms_pid_kinematics_pass(double dp, double th, double ph) {
  return std::isfinite(dp) && std::isfinite(th) && std::isfinite(ph) &&
         std::fabs(dp) <= 8.5 && std::fabs(th) <= 0.09 && std::fabs(ph) <= 0.09;
}

static void accumulate_events_first_pass(const std::string &path,
                                         int trig_num,
                                         double expected_current,
                                         const Config &cfg,
                                         const SegmentSelection &selection,
                                         EventAccum &acc) {
  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) {
    ++acc.missing_branch_segments;
    return;
  }
  auto *t = dynamic_cast<TTree *>(file.Get("T"));
  if (!t) {
    ++acc.missing_branch_segments;
    return;
  }

  const std::string htrig_tdc_branch = "T.hms.hTRIG" + std::to_string(trig_num) + "_tdcTimeRaw";
  const char *trig6_branch = has_branch(t, "T.hms.npsTRIG6_tdcTimeRaw")
      ? "T.hms.npsTRIG6_tdcTimeRaw"
      : (has_branch(t, "T.hms.hTRIG6_tdcTimeRaw") ? "T.hms.hTRIG6_tdcTimeRaw" : nullptr);
  const std::string selected_htrig_branch =
      (trig_num == 6 && trig6_branch != nullptr) ? std::string(trig6_branch) : htrig_tdc_branch;

  const bool has_core =
      has_branch(t, "g.evnum") &&
      has_branch(t, "fEvtHdr.fEvtType") &&
      has_branch(t, "H.BCM4A.scalerCurrent") &&
      has_branch(t, "T.hms.hEDTM_tdcTimeRaw");
  if (!has_core) {
    ++acc.missing_branch_segments;
    return;
  }

  const bool has_htrig_tdc = has_branch(t, selected_htrig_branch.c_str());
  const bool has_pid =
      has_branch(t, "H.cer.npeSum") &&
      has_branch(t, "H.cal.etotnorm") &&
      has_branch(t, "H.dc.ntrack") &&
      has_branch(t, "H.gtr.dp") &&
      has_branch(t, "H.gtr.th") &&
      has_branch(t, "H.gtr.ph");
  const bool has_tracking = has_pid && has_branch(t, "H.hod.goodscinhit");
  const bool has_hodo =
      has_branch(t, "H.dc.ntrack") &&
      has_branch(t, "H.hod.1x.nhits") &&
      has_branch(t, "H.hod.1y.nhits") &&
      has_branch(t, "H.hod.2x.nhits") &&
      has_branch(t, "H.hod.2y.nhits");

  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("g.evnum", 1);
  t->SetBranchStatus("fEvtHdr.fEvtType", 1);
  t->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
  t->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
  if (has_htrig_tdc) t->SetBranchStatus(selected_htrig_branch.c_str(), 1);
  if (has_pid) {
    t->SetBranchStatus("H.cer.npeSum", 1);
    t->SetBranchStatus("H.cal.etotnorm", 1);
    t->SetBranchStatus("H.dc.ntrack", 1);
    t->SetBranchStatus("H.gtr.dp", 1);
    t->SetBranchStatus("H.gtr.th", 1);
    t->SetBranchStatus("H.gtr.ph", 1);
  } else if (has_hodo) {
    t->SetBranchStatus("H.dc.ntrack", 1);
  }
  if (has_tracking) t->SetBranchStatus("H.hod.goodscinhit", 1);
  if (has_hodo) {
    t->SetBranchStatus("H.hod.1x.nhits", 1);
    t->SetBranchStatus("H.hod.1y.nhits", 1);
    t->SetBranchStatus("H.hod.2x.nhits", 1);
    t->SetBranchStatus("H.hod.2y.nhits", 1);
  }

  double gevnum = 0.0, current = 0.0, edtm_tdc = 0.0, htrig_tdc = 0.0;
  unsigned int evt_type = 0;
  double cer = std::numeric_limits<double>::quiet_NaN();
  double cal = std::numeric_limits<double>::quiet_NaN();
  double ntrack = std::numeric_limits<double>::quiet_NaN();
  double dp = std::numeric_limits<double>::quiet_NaN();
  double th = std::numeric_limits<double>::quiet_NaN();
  double ph = std::numeric_limits<double>::quiet_NaN();
  double good_scin = std::numeric_limits<double>::quiet_NaN();
  double h1x = std::numeric_limits<double>::quiet_NaN();
  double h1y = std::numeric_limits<double>::quiet_NaN();
  double h2x = std::numeric_limits<double>::quiet_NaN();
  double h2y = std::numeric_limits<double>::quiet_NaN();

  t->SetBranchAddress("g.evnum", &gevnum);
  t->SetBranchAddress("fEvtHdr.fEvtType", &evt_type);
  t->SetBranchAddress("H.BCM4A.scalerCurrent", &current);
  t->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtm_tdc);
  if (has_htrig_tdc) t->SetBranchAddress(selected_htrig_branch.c_str(), &htrig_tdc);
  if (has_pid) {
    t->SetBranchAddress("H.cer.npeSum", &cer);
    t->SetBranchAddress("H.cal.etotnorm", &cal);
    t->SetBranchAddress("H.dc.ntrack", &ntrack);
    t->SetBranchAddress("H.gtr.dp", &dp);
    t->SetBranchAddress("H.gtr.th", &th);
    t->SetBranchAddress("H.gtr.ph", &ph);
  } else if (has_hodo) {
    t->SetBranchAddress("H.dc.ntrack", &ntrack);
  }
  if (has_tracking) t->SetBranchAddress("H.hod.goodscinhit", &good_scin);
  if (has_hodo) {
    t->SetBranchAddress("H.hod.1x.nhits", &h1x);
    t->SetBranchAddress("H.hod.1y.nhits", &h1y);
    t->SetBranchAddress("H.hod.2x.nhits", &h2x);
    t->SetBranchAddress("H.hod.2y.nhits", &h2y);
  }

  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    const bool physics_event = evt_type == 1;
    const bool good_event = !selection.apply_good_events ||
        in_ranges(std::llround(gevnum), selection.gevnum_ranges);
    const bool selected = physics_event && good_event;
    const bool current_ok = current_in_mode_window(current, selection, expected_current, cfg.window_frac);
    if (current_ok) ++acc.current_window_events;
    if (selected && current_ok) ++acc.selected_events;

    if (selected && current_ok && std::isfinite(edtm_tdc) && edtm_tdc > 0.0) {
      ++acc.edtm_hist[static_cast<int>(edtm_tdc / 10.0)];
    }
    if (selected && current_ok && has_htrig_tdc) {
      const bool edtm_acc = std::isfinite(edtm_tdc) && edtm_tdc != 0.0;
      const bool htrig_acc = std::isfinite(htrig_tdc) && htrig_tdc != 0.0;
      if (edtm_acc) ++acc.tdc_edtm_accepted;
      if (htrig_acc) ++acc.tdc_htrig_all_accepted;
      if (htrig_acc && !edtm_acc) ++acc.tdc_htrig_phy_accepted;
    }

    if (has_pid) {
      const bool edt_ok = std::isfinite(edtm_tdc) && edtm_tdc < 0.1;
      const bool track_ok = !cfg.pid_require_track || (std::isfinite(ntrack) && ntrack > 0.5);
      const bool den = selected && current_ok && edt_ok && track_ok && hms_pid_kinematics_pass(dp, th, ph);
      const bool cer_pass = std::isfinite(cer) && cer > cfg.pid_cer_npe_min;
      const bool cal_pass = std::isfinite(cal) &&
          cal > cfg.pid_cal_etotnorm_min &&
          cal <= cfg.pid_cal_etotnorm_max;
      if (den) {
        ++acc.pid.den;
        if (cer_pass) {
          ++acc.pid.cer;
          ++acc.pid.cal_tag_den;
        }
        if (cal_pass) {
          ++acc.pid.cal;
          ++acc.pid.cer_tag_den;
        }
        if (cer_pass && cal_pass) ++acc.pid.both;
      }
    }

    if (has_tracking) {
      const bool good_scin_hit = std::isfinite(good_scin) && std::llround(good_scin) == 1;
      const bool electron_pid_should =
          std::isfinite(cal) && cal > 0.6 && cal < 2.0 &&
          std::isfinite(cer) && cer > 0.5;
      if (selected && current_ok && good_scin_hit && electron_pid_should) {
        ++acc.tracking.should;
        if (std::isfinite(ntrack) && ntrack > 0.5) ++acc.tracking.did;
      }
    }

    if (has_hodo) {
      const bool edt_ok = std::isfinite(edtm_tdc) && edtm_tdc < 0.1;
      if (selected && current_ok && edt_ok && std::isfinite(ntrack) && ntrack > 0.5) {
        ++acc.hodo.den;
        const bool good1x = std::isfinite(h1x) && h1x > 0.0 && h1x < 3.0;
        const bool good1y = std::isfinite(h1y) && h1y > 0.0 && h1y < 3.0;
        const bool good2x = std::isfinite(h2x) && h2x > 0.0 && h2x < 3.0;
        const bool good2y = std::isfinite(h2y) && h2y > 0.0 && h2y < 3.0;
        const int good_planes = static_cast<int>(good1x) + static_cast<int>(good1y) +
                                static_cast<int>(good2x) + static_cast<int>(good2y);
        if (good_planes >= 3) ++acc.hodo.num;
      }
    }
  }
}

static double compute_edtm_peak(const std::map<int, long long> &hist) {
  if (hist.empty()) return std::numeric_limits<double>::quiet_NaN();
  int best_bin = 0;
  long long best = -1;
  for (const auto &kv : hist) {
    if (kv.second > best) {
      best = kv.second;
      best_bin = kv.first;
    }
  }
  if (best <= 0) return std::numeric_limits<double>::quiet_NaN();
  return (static_cast<double>(best_bin) + 0.5) * 10.0;
}

static long long count_edtm_peak_events(const std::string &path,
                                        double expected_current,
                                        const Config &cfg,
                                        const SegmentSelection &selection,
                                        double peak,
                                        EventAccum &acc) {
  if (!std::isfinite(peak)) return 0;
  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) {
    ++acc.missing_branch_segments;
    return 0;
  }
  auto *t = dynamic_cast<TTree *>(file.Get("T"));
  if (!has_branch(t, "g.evnum") ||
      !has_branch(t, "fEvtHdr.fEvtType") ||
      !has_branch(t, "H.BCM4A.scalerCurrent") ||
      !has_branch(t, "T.hms.hEDTM_tdcTimeRaw")) {
    ++acc.missing_branch_segments;
    return 0;
  }

  double gevnum = 0.0, current = 0.0, edtm_tdc = 0.0;
  unsigned int evt_type = 0;
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("g.evnum", 1);
  t->SetBranchStatus("fEvtHdr.fEvtType", 1);
  t->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
  t->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
  t->SetBranchAddress("g.evnum", &gevnum);
  t->SetBranchAddress("fEvtHdr.fEvtType", &evt_type);
  t->SetBranchAddress("H.BCM4A.scalerCurrent", &current);
  t->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtm_tdc);

  long long count = 0;
  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    const bool selected = (evt_type == 1) &&
        (!selection.apply_good_events || in_ranges(std::llround(gevnum), selection.gevnum_ranges));
    const bool current_ok = current_in_mode_window(current, selection, expected_current, cfg.window_frac);
    if (selected && current_ok && std::isfinite(edtm_tdc) && edtm_tdc > 0.0 &&
        std::fabs(edtm_tdc - peak) <= 500.0) {
      ++count;
    }
  }
  return count;
}

static RunResult analyze_run_mode(int run,
                                  const std::vector<std::string> &files,
                                  const Config &cfg,
                                  bool apply_good_events,
                                  double expected_current,
                                  int trig_num,
                                  double ps,
                                  const std::string &prescale_token) {
  RunResult rr;
  rr.run = run;
  rr.n_files = static_cast<int>(files.size());
  rr.good_event_selection = apply_good_events;
  rr.expected_current_uA = expected_current;
  rr.ps_factor = ps;
  rr.prescale_token = prescale_token;
  rr.which_TRIG = "H.hTRIG" + std::to_string(trig_num) + ".scaler";

  ScalerAccum scaler;
  HelChargeAccum hel;
  EventAccum events;
  std::vector<std::pair<std::string, SegmentSelection>> usable_segments;

  for (const auto &path : files) {
    SegmentSelection selection = build_segment_selection(path, cfg, apply_good_events, expected_current);
    if (!selection.ok) {
      ++rr.selection_failed_segments;
      if (cfg.progress) {
        std::cerr << "[run " << run << "] selection failed for "
                  << fs::path(path).filename().string() << ": "
                  << selection.message << "\n";
      }
      continue;
    }

    if (std::isfinite(selection.i0_used_uA)) {
      if (!std::isfinite(rr.selection_i0_uA)) rr.selection_i0_uA = 0.0;
      rr.selection_i0_uA += selection.i0_used_uA;
    }
    if (std::isfinite(selection.current_min_uA)) {
      if (!std::isfinite(rr.selection_current_min_uA)) rr.selection_current_min_uA = 0.0;
      rr.selection_current_min_uA += selection.current_min_uA;
    }
    if (std::isfinite(selection.current_max_uA)) {
      if (!std::isfinite(rr.selection_current_max_uA)) rr.selection_current_max_uA = 0.0;
      rr.selection_current_max_uA += selection.current_max_uA;
    }
    rr.selection_scaler_reads_pre += selection.n_scaler_reads_pre;
    rr.selection_scaler_reads_post += selection.n_scaler_reads_post;

    accumulate_hel_charge(path, selection, hel);
    accumulate_scalers(path, trig_num, expected_current, cfg, selection, scaler);
    accumulate_events_first_pass(path, trig_num, expected_current, cfg, selection, events);
    usable_segments.push_back({path, selection});
  }

  rr.n_segments_used = static_cast<int>(usable_segments.size());
  if (rr.n_segments_used > 0) {
    if (std::isfinite(rr.selection_i0_uA)) rr.selection_i0_uA /= rr.n_segments_used;
    if (std::isfinite(rr.selection_current_min_uA)) rr.selection_current_min_uA /= rr.n_segments_used;
    if (std::isfinite(rr.selection_current_max_uA)) rr.selection_current_max_uA /= rr.n_segments_used;
  }

  const double edtm_peak = compute_edtm_peak(events.edtm_hist);
  long long edtm_peak_count = 0;
  for (const auto &seg : usable_segments) {
    edtm_peak_count += count_edtm_peak_events(
        seg.first, expected_current, cfg, seg.second, edtm_peak, events);
  }

  rr.avg_current_uA = safe_div(scaler.current_weighted, scaler.beam_time);
  rr.charge_uC = scaler.beam_time > 0.0 ? scaler.charge_uC : std::numeric_limits<double>::quiet_NaN();
  rr.hel0_charge_before_cut_uC = hel.hel0_before;
  rr.hel_neg_charge_before_cut_uC = hel.hel_neg_before;
  rr.hel_pos_charge_before_cut_uC = hel.hel_pos_before;
  rr.hel_charge_before_cut_uC = hel.hel0_before + hel.hel_neg_before + hel.hel_pos_before;
  rr.hel0_charge_after_cut_uC = hel.hel0_after;
  rr.hel_neg_charge_after_cut_uC = hel.hel_neg_after;
  rr.hel_pos_charge_after_cut_uC = hel.hel_pos_after;
  rr.hel_charge_after_cut_uC = hel.hel0_after + hel.hel_neg_after + hel.hel_pos_after;

  rr.scaler_htrig_total = scaler.scaler_htrig_total;
  rr.scaler_edtm_total = scaler.scaler_edtm_total;
  rr.trig_accp_total = scaler.trig_accp_total;
  rr.scaler_edtm_no_cut = scaler.scaler_edtm_no_cut;
  rr.trig_accp_no_cut = scaler.trig_accp_no_cut;
  rr.beam_time = scaler.beam_time;
  rr.beam_on_frac = safe_div(scaler.beam_time, scaler.total_valid_time);
  rr.beam_on_percent_edtm = clamp01(safe_div(scaler.scaler_edtm_total, scaler.scaler_edtm_no_cut));
  rr.beam_on_percent_trig_accp = clamp01(safe_div(scaler.trig_accp_total, scaler.trig_accp_no_cut));

  rr.htrig_rate_hz = safe_div(scaler.scaler_htrig_total, scaler.beam_time);
  rr.edtm_rate_hz = safe_div(scaler.scaler_edtm_total, scaler.beam_time);
  rr.trig_accp_rate_hz = safe_div(scaler.trig_accp_total, scaler.beam_time);
  rr.s1x_rate_hz_from_scaler = safe_div(scaler.s1x_scaler_total, scaler.beam_time);
  if (scaler.s1x_rate_samples > 0) {
    rr.s1x_scalerRate_mean_hz =
        scaler.s1x_rate_sample_sum / static_cast<double>(scaler.s1x_rate_samples);
    const double mean2 = scaler.s1x_rate_sample_sum2 / static_cast<double>(scaler.s1x_rate_samples);
    rr.s1x_scalerRate_rms_hz =
        std::sqrt(std::max(0.0, mean2 - rr.s1x_scalerRate_mean_hz * rr.s1x_scalerRate_mean_hz));
  }

  rr.tdc_edtm_accepted = edtm_peak_count;
  rr.tdc_htrig_all_accepted = events.tdc_htrig_all_accepted;
  rr.tdc_htrig_phy_accepted = events.tdc_htrig_phy_accepted;
  rr.selected_events = events.selected_events;
  rr.current_window_events = events.current_window_events;
  rr.tdc_edtm_rate_hz = safe_div(static_cast<double>(rr.tdc_edtm_accepted), scaler.beam_time);
  rr.tdc_htrig_all_rate_hz = safe_div(static_cast<double>(rr.tdc_htrig_all_accepted), scaler.beam_time);
  rr.tdc_htrig_phy_rate_hz = safe_div(static_cast<double>(rr.tdc_htrig_phy_accepted), scaler.beam_time);

  rr.NewGen_EDTM_livetime =
      safe_div(static_cast<double>(rr.tdc_edtm_accepted), scaler.scaler_edtm_total) * ps;
  rr.clta_livetime_tsh = safe_div(scaler.trig_accp_total, scaler.scaler_htrig_total) * ps;
  rr.clta_livetime_tdc =
      safe_div(static_cast<double>(events.tdc_htrig_all_accepted), scaler.scaler_htrig_total) *
      ps * rr.beam_on_percent_trig_accp;
  rr.cltp_livetime_tdc =
      safe_div(static_cast<double>(events.tdc_htrig_phy_accepted),
               scaler.scaler_htrig_total - scaler.scaler_edtm_total) *
      ps * rr.beam_on_percent_trig_accp;

  rr.HMS_pid_den_events = events.pid.den;
  rr.HMS_pid_eff = safe_div(static_cast<double>(events.pid.both), static_cast<double>(events.pid.den));
  rr.HMS_cer_eff = safe_div(static_cast<double>(events.pid.cer), static_cast<double>(events.pid.den));
  rr.HMS_cal_eff = safe_div(static_cast<double>(events.pid.cal), static_cast<double>(events.pid.den));
  rr.HMS_cer_eff_tag_cal =
      safe_div(static_cast<double>(events.pid.both), static_cast<double>(events.pid.cer_tag_den));
  rr.HMS_cal_eff_tag_cer =
      safe_div(static_cast<double>(events.pid.both), static_cast<double>(events.pid.cal_tag_den));
  rr.HMS_tracking_should_events = events.tracking.should;
  rr.HMS_tracking_did_events = events.tracking.did;
  rr.HMS_tracking_eff =
      safe_div(static_cast<double>(events.tracking.did), static_cast<double>(events.tracking.should));
  rr.HMS_tracking_eff_err = binomial_err(events.tracking.did, events.tracking.should);
  rr.HMS_hodo_3of4_den_events = events.hodo.den;
  rr.HMS_hodo_3of4_eff =
      safe_div(static_cast<double>(events.hodo.num), static_cast<double>(events.hodo.den));

  rr.negative_dt_intervals = scaler.negative_dt_intervals;
  rr.non_monotonic_scaler_steps = scaler.non_monotonic_scaler_steps;
  rr.suspicious_scaler_jump_steps = scaler.suspicious_scaler_jump_steps;
  rr.missing_branch_segments = scaler.missing_branch_segments + static_cast<int>(events.missing_branch_segments);
  rr.hel_missing_branch_segments = hel.missing_branch_segments;
  rr.scaler_evcount_gate_segments = scaler.scaler_evcount_gate_segments;
  std::ostringstream missing;
  bool first = true;
  for (const auto &b : scaler.missing_branches) {
    if (!first) missing << ";";
    first = false;
    missing << b;
  }
  rr.missing_branches = missing.str();
  return rr;
}

static void write_csv(const std::string &path, const std::vector<RunResult> &rows) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Failed to open output CSV: " + path);

  out << "run,n_files,n_segments_used,selection_failed_segments,good_event_selection,"
      << "expected_current_uA,selection_i0_uA,selection_current_min_uA,selection_current_max_uA,"
      << "selection_scaler_reads_pre,selection_scaler_reads_post,"
      << "avg_current_uA,charge_uC,"
      << "HEL_charge_before_cut_uC,HEL_neg_charge_before_cut_uC,HEL_pos_charge_before_cut_uC,HEL0_charge_before_cut_uC,"
      << "HEL_charge_after_cut_uC,HEL_neg_charge_after_cut_uC,HEL_pos_charge_after_cut_uC,HEL0_charge_after_cut_uC,"
      << "scaler_htrig_total,scaler_edtm_total,trig_accp_total,scaler_edtm_no_cut,trig_accp_no_cut,"
      << "beam_time,beam_on_frac,beam_on_percent_edtm,beam_on_percent_trig_accp,"
      << "htrig_rate_hz,edtm_rate_hz,trig_accp_rate_hz,"
      << "s1x_rate_hz_from_scaler,s1x_scalerRate_mean_hz,s1x_scalerRate_rms_hz,"
      << "tdc_edtm_accepted,tdc_htrig_all_accepted,tdc_htrig_phy_accepted,selected_events,current_window_events,"
      << "tdc_edtm_rate_hz,tdc_htrig_all_rate_hz,tdc_htrig_phy_rate_hz,"
      << "NewGen_EDTM_livetime,clta_livetime_tsh,clta_livetime_tdc,cltp_livetime_tdc,"
      << "HMS_pid_den_events,HMS_pid_eff,HMS_cer_eff,HMS_cal_eff,HMS_cer_eff_tag_cal,HMS_cal_eff_tag_cer,"
      << "HMS_tracking_should_events,HMS_tracking_did_events,HMS_tracking_eff,HMS_tracking_eff_err,"
      << "HMS_hodo_3of4_den_events,HMS_hodo_3of4_eff,"
      << "prescale_token,ps_factor,which_TRIG,"
      << "negative_dt_intervals,non_monotonic_scaler_steps,suspicious_scaler_jump_steps,"
      << "missing_branch_segments,hel_missing_branch_segments,scaler_evcount_gate_segments,missing_branches\n";

  out << std::fixed << std::setprecision(8);
  for (const auto &r : rows) {
    out << r.run << ','
        << r.n_files << ','
        << r.n_segments_used << ','
        << r.selection_failed_segments << ','
        << (r.good_event_selection ? 1 : 0) << ','
        << r.expected_current_uA << ','
        << r.selection_i0_uA << ','
        << r.selection_current_min_uA << ','
        << r.selection_current_max_uA << ','
        << r.selection_scaler_reads_pre << ','
        << r.selection_scaler_reads_post << ','
        << r.avg_current_uA << ','
        << r.charge_uC << ','
        << r.hel_charge_before_cut_uC << ','
        << r.hel_neg_charge_before_cut_uC << ','
        << r.hel_pos_charge_before_cut_uC << ','
        << r.hel0_charge_before_cut_uC << ','
        << r.hel_charge_after_cut_uC << ','
        << r.hel_neg_charge_after_cut_uC << ','
        << r.hel_pos_charge_after_cut_uC << ','
        << r.hel0_charge_after_cut_uC << ','
        << r.scaler_htrig_total << ','
        << r.scaler_edtm_total << ','
        << r.trig_accp_total << ','
        << r.scaler_edtm_no_cut << ','
        << r.trig_accp_no_cut << ','
        << r.beam_time << ','
        << r.beam_on_frac << ','
        << r.beam_on_percent_edtm << ','
        << r.beam_on_percent_trig_accp << ','
        << r.htrig_rate_hz << ','
        << r.edtm_rate_hz << ','
        << r.trig_accp_rate_hz << ','
        << r.s1x_rate_hz_from_scaler << ','
        << r.s1x_scalerRate_mean_hz << ','
        << r.s1x_scalerRate_rms_hz << ','
        << r.tdc_edtm_accepted << ','
        << r.tdc_htrig_all_accepted << ','
        << r.tdc_htrig_phy_accepted << ','
        << r.selected_events << ','
        << r.current_window_events << ','
        << r.tdc_edtm_rate_hz << ','
        << r.tdc_htrig_all_rate_hz << ','
        << r.tdc_htrig_phy_rate_hz << ','
        << r.NewGen_EDTM_livetime << ','
        << r.clta_livetime_tsh << ','
        << r.clta_livetime_tdc << ','
        << r.cltp_livetime_tdc << ','
        << r.HMS_pid_den_events << ','
        << r.HMS_pid_eff << ','
        << r.HMS_cer_eff << ','
        << r.HMS_cal_eff << ','
        << r.HMS_cer_eff_tag_cal << ','
        << r.HMS_cal_eff_tag_cer << ','
        << r.HMS_tracking_should_events << ','
        << r.HMS_tracking_did_events << ','
        << r.HMS_tracking_eff << ','
        << r.HMS_tracking_eff_err << ','
        << r.HMS_hodo_3of4_den_events << ','
        << r.HMS_hodo_3of4_eff << ','
        << csv_quote(r.prescale_token) << ','
        << r.ps_factor << ','
        << csv_quote(r.which_TRIG) << ','
        << r.negative_dt_intervals << ','
        << r.non_monotonic_scaler_steps << ','
        << r.suspicious_scaler_jump_steps << ','
        << r.missing_branch_segments << ','
        << r.hel_missing_branch_segments << ','
        << r.scaler_evcount_gate_segments << ','
        << csv_quote(r.missing_branches) << '\n';
  }
}

static void print_usage(const char *prog) {
  std::cerr
      << "Usage: " << prog << " --runs <list...> [options]\n"
      << "Options:\n"
      << "  --root-dir <path>             ROOT replay directory\n"
      << "  --prescale-csv <path>         CSV containing run_number and prescale columns\n"
      << "  --out-csv <path>              Good-event output CSV\n"
      << "  --out-csv-no-good <path>      No-good-event-selection output CSV\n"
      << "  --default-ps <value>          Default prescale factor\n"
      << "  --window-frac <f>             Current-window half-width fraction (default 0.15)\n"
      << "  --current-window <uA>         Absolute half-width around expected current\n"
      << "  --use-absolute-window         Use --current-window instead of --window-frac\n"
      << "  --good-helicity <mode>        quartet_pm|quartet_all|ignore (default quartet_pm)\n"
      << "  --good-no-current-window      Disable current window inside good-event helper\n"
      << "  --good-I0 <uA>                Fixed I0 for good-event helper\n"
      << "  --stable-window-N <N>         Rolling stability window (default 30)\n"
      << "  --stability-frac-range <f>    Rolling charge spread limit (default 0.08)\n"
      << "  --pid-no-track-required       Do not require track in PID denominator\n"
      << "  --no-progress                 Reduce console output\n"
      << "  --runs <tokens...>            Required, e.g. 4076 4077 4080-4083\n";
}

static bool parse_args(int argc, char **argv, Config &cfg) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto need = [&](const std::string &key) -> std::string {
      if (i + 1 >= argc) throw std::runtime_error("Missing value after " + key);
      return argv[++i];
    };

    try {
      if (arg == "--root-dir") cfg.root_dir = need(arg);
      else if (arg == "--prescale-csv" || arg == "--db") cfg.prescale_csv = need(arg);
      else if (arg == "--out-csv") cfg.out_csv = need(arg);
      else if (arg == "--out-csv-no-good") cfg.out_csv_no_good = need(arg);
      else if (arg == "--default-ps") cfg.default_ps = std::stod(need(arg));
      else if (arg == "--window-frac" || arg == "--current-window-frac") cfg.window_frac = std::stod(need(arg));
      else if (arg == "--current-window") {
        cfg.current_window_uA = std::stod(need(arg));
        cfg.use_absolute_window = true;
      }
      else if (arg == "--use-absolute-window") cfg.use_absolute_window = true;
      else if (arg == "--good-helicity") cfg.good_helicity = need(arg);
      else if (arg == "--good-no-current-window") cfg.good_use_current_window = false;
      else if (arg == "--good-I0") cfg.good_i0_fixed_uA = std::stod(need(arg));
      else if (arg == "--stable-window-N" || arg == "--stability-window") cfg.stable_window_N = std::stoi(need(arg));
      else if (arg == "--stability-frac-range") cfg.stability_frac_range = std::stod(need(arg));
      else if (arg == "--pid-cer-npe-min") cfg.pid_cer_npe_min = std::stod(need(arg));
      else if (arg == "--pid-cal-etottracknorm-min") cfg.pid_cal_etotnorm_min = std::stod(need(arg));
      else if (arg == "--pid-cal-etottracknorm-max") cfg.pid_cal_etotnorm_max = std::stod(need(arg));
      else if (arg == "--pid-no-track-required") cfg.pid_require_track = false;
      else if (arg == "--no-progress") cfg.progress = false;
      else if (arg == "--runs") {
        ++i;
        while (i < argc && std::string(argv[i]).rfind("--", 0) != 0) {
          cfg.run_tokens.push_back(argv[i++]);
        }
        --i;
      } else if (arg == "--help" || arg == "-h") {
        cfg.show_help = true;
        print_usage(argv[0]);
        return false;
      } else {
        throw std::runtime_error("Unknown or incomplete argument: " + arg);
      }
    } catch (const std::exception &ex) {
      std::cerr << "Argument error: " << ex.what() << "\n";
      print_usage(argv[0]);
      return false;
    }
  }

  if (cfg.run_tokens.empty()) {
    std::cerr << "Error: --runs is required\n";
    print_usage(argv[0]);
    return false;
  }
  if (cfg.out_csv_no_good.empty()) {
    cfg.out_csv_no_good = with_suffix_before_extension(cfg.out_csv, "_no_good_events");
  }
  return true;
}

int main(int argc, char **argv) {
  Config cfg;
  if (!parse_args(argc, argv, cfg)) return cfg.show_help ? 0 : 1;

  std::vector<int> runs;
  try {
    runs = parse_run_list(cfg.run_tokens);
    (void)parse_helicity_mode_or_die(cfg.good_helicity);
  } catch (const std::exception &ex) {
    std::cerr << "Configuration error: " << ex.what() << "\n";
    return 1;
  }

  const std::set<int> run_set(runs.begin(), runs.end());
  const auto indexed_files = index_run_files(cfg.root_dir, run_set);
  const auto prescale_lookup = parse_prescale_csv(cfg.prescale_csv);
  if (prescale_lookup.empty()) {
    std::cerr << "[WARN] No prescale rows loaded from " << cfg.prescale_csv << "\n";
  }

  std::vector<RunResult> good_rows;
  std::vector<RunResult> no_good_rows;

  for (int run : runs) {
    const auto fit = indexed_files.find(run);
    if (fit == indexed_files.end() || fit->second.empty()) {
      std::cerr << "[WARN] No files for run " << run << " in " << cfg.root_dir << "\n";
      continue;
    }
    const auto &files = fit->second;
    if (cfg.progress) {
      std::cout << "[run " << run << "] files=" << files.size() << "\n";
    }

    std::string prescale_token;
    double ps = cfg.default_ps;
    int trig_num = 4;
    const auto pit = prescale_lookup.find(run);
    if (pit != prescale_lookup.end()) {
      prescale_token = pit->second;
      ps = prescale_from_token(prescale_token);
      trig_num = extract_trig_number(prescale_token);
      if (trig_num <= 0) {
        std::cerr << "[ERROR] Could not parse trigger number from prescale token for run "
                  << run << ": " << prescale_token << "\n";
        return 2;
      }
    } else {
      std::cerr << "[WARN] No prescale token for run " << run
                << "; using default ps=" << cfg.default_ps
                << " and H.hTRIG4\n";
    }

    double expected_current = compute_expected_current(files);
    if (!std::isfinite(expected_current) || expected_current <= 0.0) {
      expected_current = 10.0;
      std::cerr << "[WARN] Run " << run << ": using fallback expected current 10 uA\n";
    }

    good_rows.push_back(analyze_run_mode(
        run, files, cfg, true, expected_current, trig_num, ps, prescale_token));
    no_good_rows.push_back(analyze_run_mode(
        run, files, cfg, false, expected_current, trig_num, ps, prescale_token));
  }

  try {
    write_csv(cfg.out_csv, good_rows);
    write_csv(cfg.out_csv_no_good, no_good_rows);
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\n";
    return 3;
  }

  std::cout << "Wrote good-event CSV: " << cfg.out_csv
            << " (" << good_rows.size() << " runs)\n";
  std::cout << "Wrote no-good-event CSV: " << cfg.out_csv_no_good
            << " (" << no_good_rows.size() << " runs)\n";
  return 0;
}
