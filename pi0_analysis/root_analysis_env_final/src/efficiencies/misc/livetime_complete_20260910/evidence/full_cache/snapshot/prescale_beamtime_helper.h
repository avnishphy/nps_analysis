#pragma once

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <limits>
#include <regex>
#include <string>
#include <vector>

#include "eff_math_utils.h"
#include "good_event_selection_helper.h"

namespace effstuff {

struct PrescaleInfo {
  std::string prescale_token;
  double ps_factor = 1.0;
  int trig_number = 4;
  std::string which_TRIG = "H.hTRIG4.scaler";
  bool valid = true;
  bool multiple_enabled = false;
  std::string message;
};

struct PrescaleAssignment {
  int trig_number = -1;
  int setting = -1;
};

struct BeamTimeAccumulator {
  double beam_time = 0.0;
  double scaler_edtm_total = 0.0;
  double scaler_edtm_total_evcount_gated = 0.0;
  double s1x_rate_time_sum = 0.0;
  double s1x_rate2_time_sum = 0.0;
  double s1x_sample_time = 0.0;
  int negative_dt_intervals = 0;
  int non_monotonic_scaler_steps = 0;
  int missing_branch_segments = 0;
  int missing_s1x_segments = 0;
  bool used_evcount_range_gate = false;
};

inline bool tree_has_branch(TTree* tree, const char* branch_name) {
  return tree != nullptr && tree->GetBranch(branch_name) != nullptr;
}

inline std::vector<PrescaleAssignment> parse_prescale_assignments(
    const std::string& token) {
  const std::regex re("ps\\s*(\\d+)\\s*=\\s*(-?\\d+)",
                      std::regex_constants::icase);
  std::vector<PrescaleAssignment> assignments;
  for (std::sregex_iterator it(token.begin(), token.end(), re), end;
       it != end; ++it) {
    assignments.push_back({std::stoi((*it)[1].str()),
                           std::stoi((*it)[2].str())});
  }
  return assignments;
}

inline double prescale_from_setting(int setting) {
  if (setting <= 0) {
    return 1.0;
  }
  return std::pow(2.0, static_cast<double>(setting - 1)) + 1.0;
}

inline PrescaleInfo build_prescale_info(const std::string& prescale_token,
                                        double default_ps_factor = 1.0,
                                        int default_trig_number = 4) {
  PrescaleInfo info;
  info.prescale_token = trim_copy(prescale_token);

  if (info.prescale_token.empty()) {
    info.ps_factor = default_ps_factor;
    info.trig_number = default_trig_number;
    info.message = "empty token: using configured default trigger and factor";
  } else {
    const auto assignments = parse_prescale_assignments(info.prescale_token);
    std::vector<PrescaleAssignment> enabled;
    for (const auto& assignment : assignments) {
      // Hall C convention: -1 disables a trigger; settings >= 0 enable it.
      if (assignment.setting >= 0) {
        enabled.push_back(assignment);
      }
    }

    if (!enabled.empty()) {
      info.trig_number = enabled.front().trig_number;
      info.ps_factor = prescale_from_setting(enabled.front().setting);
      info.multiple_enabled = enabled.size() > 1;
      info.message = info.multiple_enabled
          ? "multiple enabled assignments: provisionally selected first listed ps" +
                std::to_string(info.trig_number) + "=" +
                std::to_string(enabled.front().setting)
          : "selected sole enabled assignment ps" +
                std::to_string(info.trig_number) + "=" +
                std::to_string(enabled.front().setting);
    } else {
      info.valid = false;
      info.ps_factor = std::numeric_limits<double>::quiet_NaN();
      info.trig_number = -1;
      info.which_TRIG.clear();
      info.message = "no enabled psN=setting assignment found";
      return info;
    }
  }

  info.which_TRIG = "H.hTRIG" + std::to_string(info.trig_number) + ".scaler";
  return info;
}

inline void accumulate_beam_time_and_scaler(const std::string& filepath,
                                            const GoodSelectionSummary& selection,
                                            BeamTimeAccumulator& acc) {
  TFile file(filepath.c_str(), "READ");
  if (file.IsZombie()) {
    ++acc.missing_branch_segments;
    return;
  }

  auto* tsh = dynamic_cast<TTree*>(file.Get("TSH"));
  if (tsh == nullptr) {
    ++acc.missing_branch_segments;
    return;
  }

  if (!tree_has_branch(tsh, "H.1MHz.scalerTime") ||
      !tree_has_branch(tsh, "H.BCM4A.scalerCurrent") ||
      !tree_has_branch(tsh, "H.EDTM.scaler")) {
    ++acc.missing_branch_segments;
    return;
  }

  const bool use_evcount_gate =
      !selection.accepted_evcount_ranges.empty() && tree_has_branch(tsh, "evcount");
  const bool have_s1x_rate = tree_has_branch(tsh, "H.S1X.scalerRate");
  if (!have_s1x_rate) {
    ++acc.missing_s1x_segments;
  }

  tsh->SetBranchStatus("*", 0);
  tsh->SetBranchStatus("H.1MHz.scalerTime", 1);
  tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
  tsh->SetBranchStatus("H.EDTM.scaler", 1);
  if (use_evcount_gate) {
    tsh->SetBranchStatus("evcount", 1);
    acc.used_evcount_range_gate = true;
  }
  if (have_s1x_rate) {
    tsh->SetBranchStatus("H.S1X.scalerRate", 1);
  }

  double tMHz = 0.0;
  double bcmCurrent = 0.0;
  double edtmScaler = 0.0;
  double evcountRaw = 0.0;
  double s1xRate = std::numeric_limits<double>::quiet_NaN();

  tsh->SetBranchAddress("H.1MHz.scalerTime", &tMHz);
  tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent);
  tsh->SetBranchAddress("H.EDTM.scaler", &edtmScaler);
  if (use_evcount_gate) {
    tsh->SetBranchAddress("evcount", &evcountRaw);
  }
  if (have_s1x_rate) {
    tsh->SetBranchAddress("H.S1X.scalerRate", &s1xRate);
  }

  const Long64_t n = tsh->GetEntries();
  if (n < 2) {
    return;
  }

  tsh->GetEntry(0);
  double prev_t = tMHz;
  double prev_edtm = edtmScaler;

  for (Long64_t i = 1; i < n; ++i) {
    tsh->GetEntry(i);

    const double dt = tMHz - prev_t;
    if (!std::isfinite(dt) || dt < 0.0) {
      ++acc.negative_dt_intervals;
      prev_t = tMHz;
      prev_edtm = edtmScaler;
      continue;
    }

    double d_edtm = edtmScaler - prev_edtm;
    if (d_edtm < 0.0) {
      ++acc.non_monotonic_scaler_steps;
    }
    d_edtm = std::max(0.0, d_edtm);

    // The production EDTM denominator is pulser/current based. Beam time keeps
    // the analysis good-scaler gate, while the evcount-gated EDTM total is
    // retained only as a diagnostic because EDTM pulses are artificial livetime
    // probes, not selected physics events.
    const bool in_current_window = current_in_selection_window(bcmCurrent, selection);
    const bool in_evcount_range =
        !use_evcount_gate ||
        evcount_value_in_ranges(static_cast<long long>(std::llround(evcountRaw)),
                                selection.accepted_evcount_ranges);
    if (in_current_window) {
      acc.scaler_edtm_total += d_edtm;
      if (in_evcount_range) {
        acc.beam_time += dt;
        acc.scaler_edtm_total_evcount_gated += d_edtm;
        if (have_s1x_rate && std::isfinite(s1xRate) && s1xRate >= 0.0 && dt > 0.0) {
          acc.s1x_rate_time_sum += s1xRate * dt;
          acc.s1x_rate2_time_sum += s1xRate * s1xRate * dt;
          acc.s1x_sample_time += dt;
        }
      }
    }

    prev_t = tMHz;
    prev_edtm = edtmScaler;
  }
}

}  // namespace effstuff
