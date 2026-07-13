#ifndef NPS_HEL_SELECT_H
#define NPS_HEL_SELECT_H

/*
The code is taken from Christine's original code; https://github.com/cploen/CHARGE_ACCT/tree/main
Some modifications, though not necessary, might be made to fit into our framework.
===============================================================================
nps_hel_select.h
-------------------------------------------------------------------------------
Purpose:
  Self-contained reusable helper for beam-stability / helicity-aware event
  selection using Hall C NPS replay ROOT files containing trees "TSHelH" and "T".

What it returns:
  - accepted evcount ranges
  - accepted evNumber ranges
  - accepted g.evnum ranges
  - ready-to-use cut strings for TSHelH and T

Usage:
  #include "nps_hel_select.h"

  SelectionSettings sel;
  sel.helicity_mode = HelicityMode::QuartetPM;   // or QuartetAll / IgnoreHelicity
  sel.stable_window_N = 30;                      // N=1 effectively disables rolling stability
  sel.max_charge_frac_spread = 0.08;
  sel.window_frac = 0.15;

  auto pick = build_event_selection(run, seg, rootfile, sel);
  if (!pick.ok) {
    std::cerr << pick.message << std::endl;
    return;
  }
  
  ROOT::RDataFrame dT("T", rootfile);
  auto d_good_T = dT.Filter(pick.gevnum_cut);
  
  ROOT::RDataFrame dH("TSHelH", rootfile);
  auto d_good_H = dH.Filter(pick.evcount_cut);
    
Notes:
  - QuartetPM    : keep only helicity != 0 and quartet-snap accepted ranges
  - QuartetAll   : keep all helicity states and quartet-snap accepted ranges
  - IgnoreHelicity: ignore helicity entirely and do NOT quartet-snap
  - Rolling charge stability is always part of the helper.
    Set stable_window_N = 1 to make it effectively off.
  - This helper preserves the original evcount -> event-number indexing convention.
===============================================================================
*/

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using std::string;
using std::vector;

// =================================================================================
// Reusable event-selection helper types
// =================================================================================
enum class HelicityMode {
  QuartetPM,      // keep only helicity ±1, then quartet-snap accepted regions
  QuartetAll,     // keep helicity -1,0,+1, then quartet-snap accepted regions
  IgnoreHelicity  // do not cut on helicity, do not quartet-snap
};

enum class I0Mode {
  Peak,   // use mode of scalerCurrent histogram above junk floor
  Fixed   // use user-provided fixed I0
};

struct SelectionSettings {
  // --- tunables ---
  HelicityMode helicity_mode = HelicityMode::QuartetPM;

  double junk_floor_uA = 2.5;         // reject junk current reads under this threshold
  double mean_current_calc_min = 2.5; // floor for computing mean current
  double mean_current_min = 2.75;     // production runs were as low as 3 uA

  bool use_current_window = true;
  I0Mode i0_mode = I0Mode::Peak;
  double i0_fixed_uA = -1.0;          // used only when i0_mode == Fixed
  double window_frac = 0.15;          // keep |I-I0|/I0 <= window_frac

  int stable_window_N = 30;              // N=1 effectively disables the rolling requirement
  double max_charge_frac_spread = 0.08;  // (qmax-qmin)/mean within rolling window

  // branch names
  string branch_scaler_current = "H.BCM4A_Hel.scalerCurrent";
  string branch_scaler_charge  = "H.BCM4A_Hel.scalerCharge";
  string branch_helicity       = "actualHelicity";
  string branch_evcount        = "evcount";
  string branch_evNumber       = "evNumber";
  string branch_evnum_T        = "g.evnum";
};

struct RangeI {
  int lo = 0;
  int hi = 0;
};

struct SelectionResult {
  bool ok = false;
  string rootfile;
  string message;

  // bookkeeping / diagnostics
  HelicityMode helicity_mode = HelicityMode::QuartetPM;
  bool quartet_snap_applied = true;

  double mean_current_uA = 0.0;
  double i0_used_uA = -1.0;
  double current_min_uA = -1.0;
  double current_max_uA = -1.0;

  int n_scaler_reads_pre = 0;
  int n_scaler_reads_post = 0;

  // accepted regions
  vector<RangeI> evcount_ranges;   // accepted evcount ranges
  vector<RangeI> evNumber_ranges;  // mapped from TSHelH evNumber
  vector<RangeI> gevnum_ranges;    // mapped from T g.evnum

  // ready-to-use filter strings
  string evcount_cut;
  string evNumber_cut;
  string gevnum_cut;
};

// =================================================================================
// Small helper utilities
// =================================================================================
static inline string helicity_mode_to_string(HelicityMode mode) {
  switch (mode) {
    case HelicityMode::QuartetPM:      return "QuartetPM";
    case HelicityMode::QuartetAll:     return "QuartetAll";
    case HelicityMode::IgnoreHelicity: return "IgnoreHelicity";
  }
  return "Unknown";
}

static inline bool helicity_requires_quartet_snap(HelicityMode mode) {
  return (mode == HelicityMode::QuartetPM || mode == HelicityMode::QuartetAll);
}

static inline string build_helicity_filter_expr(const SelectionSettings& sel) {
  const string& h = sel.branch_helicity;

  switch (sel.helicity_mode) {
    case HelicityMode::QuartetPM:
      return Form("(int)round(%s) != 0", h.c_str());

    case HelicityMode::QuartetAll:
      return ""; // keep all helicity states

    case HelicityMode::IgnoreHelicity:
      return ""; // completely ignore helicity
  }
  return "";
}

static inline vector<RangeI> toRangeI(const vector<vector<int>>& vv) {
  vector<RangeI> out;
  out.reserve(vv.size());
  for (const auto& v : vv) {
    if (v.size() >= 2) out.push_back({v[0], v[1]});
  }
  return out;
}

static inline string build_evcount_cut_string(const vector<RangeI>& ranges,
                                              const string& branch_evcount) {
  string out;
  for (const auto& r : ranges) {
    string term = Form("(%d<=(int)round(%s)&&(int)round(%s)<=%d)",
                       r.lo, branch_evcount.c_str(), branch_evcount.c_str(), r.hi);
    out = out.empty() ? term : out + "||" + term;
  }
  return out;
}

static inline string build_event_cut_string(const vector<RangeI>& ranges,
                                            const string& branch_event) {
  string out;
  for (const auto& r : ranges) {
    string term = Form("(%d<=(int)round(%s)&&(int)round(%s)<=%d)",
                       r.lo, branch_event.c_str(), branch_event.c_str(), r.hi);
    out = out.empty() ? term : out + "||" + term;
  }
  return out;
}

static inline bool tree_has_named_branch(TTree* tree, const string& branch_name) {
  return tree != nullptr && tree->GetBranch(branch_name.c_str()) != nullptr;
}

// =================================================================================
// Quartet alignment helpers
// =================================================================================
template<typename RDF>
static inline int alignStartToQuartetEnd(int evnum, const RDF& /*dfHelper*/) {
  return evnum - (evnum % 4);
}

template<typename RDF>
static inline int alignEndToQuartetEnd(int evnum, const RDF& /*dfHelper*/) {
  return evnum + (3 - (evnum % 4));
}

template<typename RDF>
static inline vector<vector<int>> buildQuartetAlignedRanges(const vector<int>& evC_arr,
                                                            const RDF& dfHelper) {
  vector<vector<int>> evtCRanges;
  if (evC_arr.empty()) return evtCRanges;

  size_t run_start_idx = 0;

  for (size_t i = 1; i <= evC_arr.size(); ++i) {
    bool end_of_run = (i == evC_arr.size()) || (evC_arr[i] != evC_arr[i - 1] + 1);
    if (!end_of_run) continue;

    int run_start_ev = evC_arr[run_start_idx];
    int run_end_ev   = evC_arr[i - 1];
    int run_len      = run_end_ev - run_start_ev + 1;

    if (run_len > 0) {
      vector<int> range;
      range.push_back(alignStartToQuartetEnd(run_start_ev, dfHelper));
      range.push_back(alignEndToQuartetEnd(run_end_ev, dfHelper));
      evtCRanges.push_back(range);
    }

    run_start_idx = i;
  }

  return evtCRanges;
}

static inline vector<vector<int>> mergeOverlappingRanges(const vector<vector<int>>& ranges) {
  vector<vector<int>> merged;
  if (ranges.empty()) return merged;

  merged.push_back(ranges[0]);

  for (size_t i = 1; i < ranges.size(); ++i) {
    int start = ranges[i][0];
    int end   = ranges[i][1];

    if (start <= merged.back()[1] + 1) {
      merged.back()[1] = std::max(merged.back()[1], end);
    } else {
      merged.push_back(ranges[i]);
    }
  }

  return merged;
}

// =================================================================================
// Rolling charge-stability helper
// =================================================================================
static inline vector<int> stableEvcountsFromChargeWindow(const vector<int>& evC_arr,
                                                         const vector<double>& Q_arr,
                                                         int N,
                                                         double fracRangeMax)
{
  vector<int> out;

  if ((int)evC_arr.size() != (int)Q_arr.size()) return out;
  if (N < 1) return out;                  // invalid
  if ((int)evC_arr.size() < N) return out;

  int runStart = 0;
  for (int i = 1; i <= (int)evC_arr.size(); i++) {
    bool endRun = (i == (int)evC_arr.size()) || (evC_arr[i] != evC_arr[i - 1] + 1);
    if (!endRun) continue;

    int runEnd = i - 1;
    int runLen = runEnd - runStart + 1;
    if (runLen >= N) {
      double sum = 0.0;
      std::deque<int> dqMin, dqMax;

      for (int j = runStart; j <= runEnd; j++) {
        sum += Q_arr[j];

        while (!dqMin.empty() && Q_arr[dqMin.back()] >= Q_arr[j]) dqMin.pop_back();
        dqMin.push_back(j);

        while (!dqMax.empty() && Q_arr[dqMax.back()] <= Q_arr[j]) dqMax.pop_back();
        dqMax.push_back(j);

        int old = j - N;
        if (old >= runStart) {
          sum -= Q_arr[old];
          if (!dqMin.empty() && dqMin.front() == old) dqMin.pop_front();
          if (!dqMax.empty() && dqMax.front() == old) dqMax.pop_front();
        }

        if (j >= runStart + (N - 1)) {
          double mean = sum / (double)N;
          if (mean > 0) {
            double qmin = Q_arr[dqMin.front()];
            double qmax = Q_arr[dqMax.front()];
            double fracRange = (qmax - qmin) / mean;

            if (fracRange <= fracRangeMax) {
              out.push_back(evC_arr[j]);
            }
          }
        }
      }
    }

    runStart = i;
  }

  return out;
}

// =================================================================================
// Reusable selection engine
// =================================================================================
static inline SelectionResult build_event_selection(const string& rootfile, const SelectionSettings& sel)
{
  SelectionResult res;
  res.rootfile = rootfile;
  res.helicity_mode = sel.helicity_mode;
  res.quartet_snap_applied = helicity_requires_quartet_snap(sel.helicity_mode);

  // ---------------------------------------------------------
  // Stable window sanity check
  // ---------------------------------------------------------
  
  if (sel.stable_window_N < 1) {
    res.ok = false;
    res.message = "Invalid stable_window_N: must be >= 1.";
    return res;
  }

  TFile file(rootfile.c_str(), "READ");
  if (file.IsZombie()) {
    res.ok = false;
    res.message = Form("Could not open ROOT file: %s", rootfile.c_str());
    return res;
  }

  auto* d_shelp = dynamic_cast<TTree*>(file.Get("TSHelH"));
  auto* d_T = dynamic_cast<TTree*>(file.Get("T"));
  if (d_shelp == nullptr || d_T == nullptr) {
    res.ok = false;
    res.message = "Missing required tree(s): TSHelH and/or T.";
    return res;
  }

  const bool needs_helicity = (sel.helicity_mode == HelicityMode::QuartetPM);

  if (!tree_has_named_branch(d_shelp, sel.branch_scaler_current) ||
      !tree_has_named_branch(d_shelp, sel.branch_scaler_charge) ||
      !tree_has_named_branch(d_shelp, sel.branch_evcount) ||
      !tree_has_named_branch(d_shelp, sel.branch_evNumber) ||
      (needs_helicity && !tree_has_named_branch(d_shelp, sel.branch_helicity))) {
    res.ok = false;
    res.message = "Missing one or more required TSHelH branches for selection.";
    return res;
  }

  if (!tree_has_named_branch(d_T, sel.branch_evnum_T)) {
    res.ok = false;
    res.message = "Missing required T branch for g.evnum mapping.";
    return res;
  }

  d_shelp->SetBranchStatus("*", 0);
  d_shelp->SetBranchStatus(sel.branch_scaler_current.c_str(), 1);
  d_shelp->SetBranchStatus(sel.branch_scaler_charge.c_str(), 1);
  d_shelp->SetBranchStatus(sel.branch_evcount.c_str(), 1);
  d_shelp->SetBranchStatus(sel.branch_evNumber.c_str(), 1);
  if (needs_helicity) {
    d_shelp->SetBranchStatus(sel.branch_helicity.c_str(), 1);
  }

  double scaler_current = 0.0;
  double scaler_charge = 0.0;
  double helicity_raw = 0.0;
  double evcount_raw = 0.0;
  double evnumber_raw = 0.0;

  d_shelp->SetBranchAddress(sel.branch_scaler_current.c_str(), &scaler_current);
  d_shelp->SetBranchAddress(sel.branch_scaler_charge.c_str(), &scaler_charge);
  d_shelp->SetBranchAddress(sel.branch_evcount.c_str(), &evcount_raw);
  d_shelp->SetBranchAddress(sel.branch_evNumber.c_str(), &evnumber_raw);
  if (needs_helicity) {
    d_shelp->SetBranchAddress(sel.branch_helicity.c_str(), &helicity_raw);
  }

  struct ScalerRow {
    double current = 0.0;
    double charge = 0.0;
    int helicity_int = 0;
    int evcount_int = 0;
    int evnumber_int = 0;
  };

  vector<ScalerRow> scaler_rows;
  vector<int> evN_arr;
  const Long64_t n_shelp = d_shelp->GetEntries();
  if (n_shelp > 0) {
    scaler_rows.reserve(static_cast<size_t>(n_shelp));
    evN_arr.reserve(static_cast<size_t>(n_shelp));
  }

  constexpr int kI0Bins = 120;
  constexpr double kI0Min = 0.0;
  constexpr double kI0Max = 60.0;
  const double i0_bin_width = (kI0Max - kI0Min) / static_cast<double>(kI0Bins);
  int i0_hist[kI0Bins] = {0};

  double mean_current_sum = 0.0;
  long long mean_current_count = 0;

  for (Long64_t i = 0; i < n_shelp; ++i) {
    d_shelp->GetEntry(i);

    const int hel_int = needs_helicity ? static_cast<int>(std::llround(helicity_raw)) : 1;
    const int evcount_int = static_cast<int>(evcount_raw);
    const int evnumber_int = static_cast<int>(evnumber_raw);

    scaler_rows.push_back({scaler_current, scaler_charge, hel_int, evcount_int, evnumber_int});
    evN_arr.push_back(evnumber_int);

    if (scaler_current > sel.mean_current_calc_min) {
      mean_current_sum += scaler_current;
      ++mean_current_count;
    }

    if (sel.use_current_window && sel.i0_mode == I0Mode::Peak && scaler_current > sel.junk_floor_uA) {
      const int bin = static_cast<int>((scaler_current - kI0Min) / i0_bin_width);
      if (bin >= 0 && bin < kI0Bins) {
        ++i0_hist[bin];
      }
    }
  }

  // ---------------------------------------------------------
  // Mean current sanity check
  // ---------------------------------------------------------
  if (mean_current_count > 0) {
    res.mean_current_uA = mean_current_sum / static_cast<double>(mean_current_count);
  } else {
    res.mean_current_uA = std::numeric_limits<double>::quiet_NaN();
  }

  if (res.mean_current_uA < sel.mean_current_min) {
    res.ok = false;
    res.message = Form("Mean current too low (%g uA < %g uA)",
                       res.mean_current_uA, sel.mean_current_min);
    return res;
  }

  // ---------------------------------------------------------
  // Determine I0 if current window is requested
  // ---------------------------------------------------------
  if (sel.use_current_window) {
    if (sel.i0_mode == I0Mode::Fixed) {
      res.i0_used_uA = sel.i0_fixed_uA;
    } else {
      int max_bin = 0;
      int max_count = -1;
      for (int b = 0; b < kI0Bins; ++b) {
        if (i0_hist[b] > max_count) {
          max_count = i0_hist[b];
          max_bin = b;
        }
      }
      res.i0_used_uA = kI0Min + (static_cast<double>(max_bin) + 0.5) * i0_bin_width;
    }

    if (!(res.i0_used_uA > 0.0)) {
      res.ok = false;
      res.message = Form("Invalid I0 value (%g uA)", res.i0_used_uA);
      return res;
    }

    res.current_min_uA = (1.0 - sel.window_frac) * res.i0_used_uA;
    res.current_max_uA = (1.0 + sel.window_frac) * res.i0_used_uA;
  }

  // ---------------------------------------------------------
  // Preselection on TSHelH
  // ---------------------------------------------------------
  vector<int> evC_arr;
  vector<double> Q_arr;
  evC_arr.reserve(scaler_rows.size());
  Q_arr.reserve(scaler_rows.size());

  for (const auto& row : scaler_rows) {
    if (!(row.current > sel.junk_floor_uA)) {
      continue;
    }
    if (needs_helicity && row.helicity_int == 0) {
      continue;
    }
    if (sel.use_current_window &&
        (row.current < res.current_min_uA || row.current > res.current_max_uA)) {
      continue;
    }
    evC_arr.push_back(row.evcount_int);
    Q_arr.push_back(row.charge);
  }

  res.n_scaler_reads_pre = static_cast<int>(evC_arr.size());

  vector<int> idx(evC_arr.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b){ return evC_arr[a] < evC_arr[b]; });

  vector<int> evC_sorted;
  vector<double> Q_sorted;
  evC_sorted.reserve(evC_arr.size());
  Q_sorted.reserve(Q_arr.size());

  for (int k : idx) {
    evC_sorted.push_back(evC_arr[k]);
    Q_sorted.push_back(Q_arr[k]);
  }

  evC_arr.swap(evC_sorted);
  Q_arr.swap(Q_sorted);

  // ---------------------------------------------------------
  // Rolling charge stability (always part of helper)
  // N=1 effectively disables the rolling requirement
  // ---------------------------------------------------------
  evC_arr = stableEvcountsFromChargeWindow(
    evC_arr,
    Q_arr,
    sel.stable_window_N,
    sel.max_charge_frac_spread
  );

  res.n_scaler_reads_post = (int)evC_arr.size();

  if (evC_arr.empty()) {
    res.ok = false;
    res.message = "No scaler reads survived selection.";
    return res;
  }

  // ---------------------------------------------------------
  // Build accepted evcount ranges
  // ---------------------------------------------------------
  vector<vector<int>> evtCRanges_raw;

  if (helicity_requires_quartet_snap(sel.helicity_mode)) {
    evtCRanges_raw = buildQuartetAlignedRanges(evC_arr, 0);
  } else {
    size_t run_start_idx = 0;
    for (size_t i = 1; i <= evC_arr.size(); ++i) {
      bool end_of_run = (i == evC_arr.size()) || (evC_arr[i] != evC_arr[i - 1] + 1);
      if (!end_of_run) continue;

      evtCRanges_raw.push_back({evC_arr[run_start_idx], evC_arr[i - 1]});
      run_start_idx = i;
    }
  }

  auto evtCRanges_merged_vv = mergeOverlappingRanges(evtCRanges_raw);
  res.evcount_ranges = toRangeI(evtCRanges_merged_vv);

  // ---------------------------------------------------------
  // Build mapping arrays
  // ---------------------------------------------------------
  std::sort(evN_arr.begin(), evN_arr.end());

  d_T->SetBranchStatus("*", 0);
  d_T->SetBranchStatus(sel.branch_evnum_T.c_str(), 1);

  double gevnum_raw = 0.0;
  d_T->SetBranchAddress(sel.branch_evnum_T.c_str(), &gevnum_raw);

  vector<int> gEvnumArr;
  const Long64_t n_t = d_T->GetEntries();
  if (n_t > 0) {
    gEvnumArr.reserve(static_cast<size_t>(n_t));
  }

  for (Long64_t i = 0; i < n_t; ++i) {
    d_T->GetEntry(i);
    gEvnumArr.push_back(static_cast<int>(gevnum_raw));
  }
  std::sort(gEvnumArr.begin(), gEvnumArr.end());

  // ---------------------------------------------------------
  // Map evcount ranges -> evNumber and g.evnum ranges
  // NOTE: preserves original indexing convention
  // ---------------------------------------------------------
  for (const auto& r : res.evcount_ranges) {
    if (r.lo < 1 || r.hi < 1) continue;
    if (r.lo > (int)evN_arr.size() || r.hi > (int)evN_arr.size()) continue;
    if (r.lo > (int)gEvnumArr.size() || r.hi > (int)gEvnumArr.size()) continue;

    res.evNumber_ranges.push_back({evN_arr[r.lo - 1], evN_arr[r.hi - 1]});
    res.gevnum_ranges.push_back({gEvnumArr[r.lo - 1], gEvnumArr[r.hi - 1]});
  }

  // ---------------------------------------------------------
  // Build cut strings
  // ---------------------------------------------------------
  res.evcount_cut = build_evcount_cut_string(res.evcount_ranges, sel.branch_evcount);
  res.evNumber_cut = build_event_cut_string(res.evNumber_ranges, sel.branch_evNumber);
  res.gevnum_cut   = build_event_cut_string(res.gevnum_ranges, sel.branch_evnum_T);

  res.ok = true;
  res.message = "Selection built successfully.";
  return res;
}

#endif
