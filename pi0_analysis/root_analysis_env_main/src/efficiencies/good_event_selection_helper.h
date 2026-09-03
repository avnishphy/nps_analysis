#pragma once

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "efficiency_types.h"
#include "NPS_selection_helper.h"

namespace effstuff {

inline SelectionSettings make_default_selection_settings() {
  SelectionSettings settings;
  settings.helicity_mode = HelicityMode::QuartetPM;

  settings.junk_floor_uA = 2.5;
  settings.mean_current_calc_min = 2.5;
  settings.mean_current_min = 2.75;

  settings.use_current_window = true;
  settings.i0_mode = I0Mode::Peak;
  settings.i0_fixed_uA = -1.0;
  settings.window_frac = 0.15;

  settings.stable_window_N = 30;
  settings.max_charge_frac_spread = 0.08;

  settings.branch_scaler_current = "H.BCM4A_Hel.scalerCurrent";
  settings.branch_scaler_charge = "H.BCM4A_Hel.scalerCharge";
  settings.branch_helicity = "actualHelicity";
  settings.branch_evcount = "evcount";
  settings.branch_evNumber = "evNumber";
  settings.branch_evnum_T = "g.evnum";

  return settings;
}

inline bool event_value_in_ranges(long long value,
                                  const std::vector<EventRange>& ranges) {
  if (ranges.empty()) {
    return true;
  }

  for (const auto& range : ranges) {
    if (value >= range.lo && value < range.hi) {
      return true;
    }
  }

  return false;
}

inline bool evcount_value_in_ranges(long long value,
                                    const std::vector<EventRange>& ranges) {
  if (ranges.empty()) {
    return true;
  }

  for (const auto& range : ranges) {
    if (value > range.lo && value <= range.hi) {
      return true;
    }
  }

  return false;
}

inline std::string ranges_to_string(const std::vector<EventRange>& ranges) {
  std::ostringstream out;
  bool first = true;
  for (const auto& range : ranges) {
    if (!first) {
      out << ";";
    }
    first = false;
    out << range.lo << "-" << range.hi;
  }
  return out.str();
}

inline GoodSelectionSummary build_good_selection_summary(
    const std::string& rootfile,
    const SelectionSettings& settings) {
  GoodSelectionSummary summary;

  const SelectionResult pick = build_event_selection(rootfile, settings);
  if (!pick.ok) {
    summary.ok = false;
    summary.message = pick.message;
    return summary;
  }

  summary.ok = true;
  summary.message = pick.message;
  summary.helicity_mode = helicity_mode_to_string(pick.helicity_mode);
  summary.quartet_snap_applied = pick.quartet_snap_applied;
  summary.evcount_cut = pick.evcount_cut;
  summary.evnumber_cut = pick.evNumber_cut;
  summary.gevnum_cut = pick.gevnum_cut;
  summary.mean_current_uA = pick.mean_current_uA;
  summary.i0_used_uA = pick.i0_used_uA;
  summary.n_scaler_reads_pre = pick.n_scaler_reads_pre;
  summary.n_scaler_reads_post = pick.n_scaler_reads_post;
  summary.current_min_uA = pick.current_min_uA;
  summary.current_max_uA = pick.current_max_uA;
  summary.current_window_enabled = settings.use_current_window;

  for (const auto& range : pick.evcount_ranges) {
    summary.accepted_evcount_ranges.push_back({range.lo, range.hi});
  }
  for (const auto& range : pick.evNumber_ranges) {
    summary.accepted_evnumber_ranges.push_back({range.lo, range.hi});
  }
  for (const auto& range : pick.gevnum_ranges) {
    summary.accepted_gevnum_ranges.push_back({range.lo, range.hi});
  }

  try {
    TFile file(rootfile.c_str(), "READ");
    if (file.IsZombie()) {
      summary.ok = false;
      summary.message = "Failed to open ROOT file for HEL charge summary.";
      return summary;
    }

    auto* tsh = dynamic_cast<TTree*>(file.Get("TSHelH"));
    if (tsh == nullptr) {
      summary.ok = false;
      summary.message = "TSHelH tree not found for HEL charge summary.";
      return summary;
    }

    const bool has_required =
        tsh->GetBranch(settings.branch_helicity.c_str()) != nullptr &&
        tsh->GetBranch(settings.branch_scaler_charge.c_str()) != nullptr &&
        tsh->GetBranch(settings.branch_evcount.c_str()) != nullptr;
    if (!has_required) {
      summary.ok = false;
      summary.message = "Missing TSHelH branches needed for HEL charge summary.";
      return summary;
    }

    tsh->SetBranchStatus("*", 0);
    tsh->SetBranchStatus(settings.branch_helicity.c_str(), 1);
    tsh->SetBranchStatus(settings.branch_scaler_charge.c_str(), 1);
    tsh->SetBranchStatus(settings.branch_evcount.c_str(), 1);

    double helicity_raw = 0.0;
    double charge_raw = 0.0;
    double evcount_raw = 0.0;

    tsh->SetBranchAddress(settings.branch_helicity.c_str(), &helicity_raw);
    tsh->SetBranchAddress(settings.branch_scaler_charge.c_str(), &charge_raw);
    tsh->SetBranchAddress(settings.branch_evcount.c_str(), &evcount_raw);

    double hel0_before = 0.0;
    double hel_neg_before = 0.0;
    double hel_pos_before = 0.0;
    double hel0_after = 0.0;
    double hel_neg_after = 0.0;
    double hel_pos_after = 0.0;

    const bool apply_after_cut = !summary.accepted_evcount_ranges.empty();

    const Long64_t n = tsh->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
      tsh->GetEntry(i);

      const long long hel = std::llround(helicity_raw);
      if (hel == 0) {
        hel0_before += charge_raw;
      } else if (hel < 0) {
        hel_neg_before += charge_raw;
      } else {
        hel_pos_before += charge_raw;
      }

      if (!apply_after_cut) {
        continue;
      }

      const long long evcount = std::llround(evcount_raw);
      if (!evcount_value_in_ranges(evcount, summary.accepted_evcount_ranges)) {
        continue;
      }

      if (hel == 0) {
        hel0_after += charge_raw;
      } else if (hel < 0) {
        hel_neg_after += charge_raw;
      } else {
        hel_pos_after += charge_raw;
      }
    }

    summary.hel0_charge_before_cut_uC = hel0_before;
    summary.hel_neg_charge_before_cut_uC = hel_neg_before;
    summary.hel_pos_charge_before_cut_uC = hel_pos_before;
    summary.hel_charge_before_cut_uC = hel0_before + hel_neg_before + hel_pos_before;

    if (apply_after_cut) {
      summary.hel0_charge_after_cut_uC = hel0_after;
      summary.hel_neg_charge_after_cut_uC = hel_neg_after;
      summary.hel_pos_charge_after_cut_uC = hel_pos_after;
      summary.hel_charge_after_cut_uC = hel0_after + hel_neg_after + hel_pos_after;
    } else {
      summary.hel_charge_after_cut_uC = 0.0;
      summary.hel0_charge_after_cut_uC = 0.0;
      summary.hel_neg_charge_after_cut_uC = 0.0;
      summary.hel_pos_charge_after_cut_uC = 0.0;
    }
  } catch (const std::exception& ex) {
    summary.ok = false;
    summary.message = std::string("Failed to compute HEL charge summary: ") + ex.what();
  }

  return summary;
}

inline bool current_in_selection_window(double current_uA,
                                        const GoodSelectionSummary& selection) {
  if (!selection.current_window_enabled) {
    return true;
  }
  if (!std::isfinite(selection.current_min_uA) || !std::isfinite(selection.current_max_uA)) {
    return true;
  }
  return current_uA >= selection.current_min_uA && current_uA <= selection.current_max_uA;
}

}  // namespace effstuff
