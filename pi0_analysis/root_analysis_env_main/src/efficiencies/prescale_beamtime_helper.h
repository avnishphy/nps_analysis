#pragma once

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <limits>
#include <regex>
#include <string>

#include "good_event_selection_helper.h"

namespace effstuff {

struct PrescaleInfo {
  std::string prescale_token;
  double ps_factor = 1.0;
  int trig_number = 4;
  std::string which_TRIG = "H.hTRIG4.scaler";
};

struct BeamTimeAccumulator {
  double beam_time = 0.0;
  double scaler_edtm_total = 0.0;
  double scaler_edtm_total_evcount_gated = 0.0;
  int negative_dt_intervals = 0;
  int non_monotonic_scaler_steps = 0;
  int missing_branch_segments = 0;
  bool used_evcount_range_gate = false;
};

inline bool tree_has_branch(TTree* tree, const char* branch_name) {
  return tree != nullptr && tree->GetBranch(branch_name) != nullptr;
}

inline int extract_prescale_r(const std::string& token_in) {
  std::string token = trim_copy(token_in);
  if (token.empty()) {
    return std::numeric_limits<int>::min();
  }

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
    const std::string rhs = trim_copy(token.substr(eq + 1));
    if (!rhs.empty()) {
      try {
        return std::stoi(rhs);
      } catch (...) {
        return std::numeric_limits<int>::min();
      }
    }
  }

  return std::numeric_limits<int>::min();
}

inline int extract_trig_number(const std::string& token_in) {
  std::string token = trim_copy(token_in);
  if (token.empty()) {
    return -1;
  }

  const std::regex re("ps(\\d+)", std::regex_constants::icase);
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

inline double prescale_from_token(const std::string& token) {
  const int r = extract_prescale_r(token);
  if (r == std::numeric_limits<int>::min() || r <= 0) {
    return 1.0;
  }
  return std::pow(2.0, static_cast<double>(r - 1)) + 1.0;
}

inline PrescaleInfo build_prescale_info(const std::string& prescale_token,
                                        double default_ps_factor = 1.0,
                                        int default_trig_number = 4) {
  PrescaleInfo info;
  info.prescale_token = trim_copy(prescale_token);

  if (info.prescale_token.empty()) {
    info.ps_factor = default_ps_factor;
    info.trig_number = default_trig_number;
  } else {
    info.ps_factor = prescale_from_token(info.prescale_token);
    info.trig_number = extract_trig_number(info.prescale_token);
    if (info.trig_number <= 0) {
      info.trig_number = default_trig_number;
    }
    if (!std::isfinite(info.ps_factor) || info.ps_factor <= 0.0) {
      info.ps_factor = default_ps_factor;
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

  tsh->SetBranchStatus("*", 0);
  tsh->SetBranchStatus("H.1MHz.scalerTime", 1);
  tsh->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
  tsh->SetBranchStatus("H.EDTM.scaler", 1);
  if (use_evcount_gate) {
    tsh->SetBranchStatus("evcount", 1);
    acc.used_evcount_range_gate = true;
  }

  double tMHz = 0.0;
  double bcmCurrent = 0.0;
  double edtmScaler = 0.0;
  double evcountRaw = 0.0;

  tsh->SetBranchAddress("H.1MHz.scalerTime", &tMHz);
  tsh->SetBranchAddress("H.BCM4A.scalerCurrent", &bcmCurrent);
  tsh->SetBranchAddress("H.EDTM.scaler", &edtmScaler);
  if (use_evcount_gate) {
    tsh->SetBranchAddress("evcount", &evcountRaw);
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
        value_in_ranges(static_cast<long long>(std::llround(evcountRaw)),
                        selection.accepted_evcount_ranges);
    if (in_current_window) {
      acc.scaler_edtm_total += d_edtm;
      if (in_evcount_range) {
        acc.beam_time += dt;
        acc.scaler_edtm_total_evcount_gated += d_edtm;
      }
    }

    prev_t = tMHz;
    prev_edtm = edtmScaler;
  }
}

}  // namespace effstuff
