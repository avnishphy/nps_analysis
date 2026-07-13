#pragma once

#include <cmath>

#include "eff_math_utils.h"
#include "hms_pid_common.h"

namespace effstuff {

struct HMSHodo3of4EffAccumulator {
  long long denominator_events = 0;
  long long numerator_events = 0;
};

// Hodoscope 3/4 efficiency uses clean selected, non-EDTM HMS events with a DC
// track as the denominator. The numerator requires at least N good hodo planes,
// where all hodo thresholds live in HMSPidCuts.
inline void accumulate_HMS_hodo_3of4_eff(HMSHodo3of4EffAccumulator& acc,
                                         const HMSEventContext& evt,
                                         const HMSPidCuts& cuts) {
  const bool den_sel = evt.selected && hms_edtm_pass(evt, cuts) &&
                       hms_track_pass(evt, cuts);
  if (!den_sel) {
    return;
  }

  ++acc.denominator_events;

  const bool good1x = hms_hodo_plane_pass(evt.hod_1x_nhits, cuts);
  const bool good1y = hms_hodo_plane_pass(evt.hod_1y_nhits, cuts);
  const bool good2x = hms_hodo_plane_pass(evt.hod_2x_nhits, cuts);
  const bool good2y = hms_hodo_plane_pass(evt.hod_2y_nhits, cuts);

  const int good_planes = static_cast<int>(good1x) + static_cast<int>(good1y) +
                          static_cast<int>(good2x) + static_cast<int>(good2y);
  if (good_planes >= cuts.hod_min_good_planes) {
    ++acc.numerator_events;
  }
}

inline double HMS_hodo_3of4_eff(const HMSHodo3of4EffAccumulator& acc) {
  return safe_div(static_cast<double>(acc.numerator_events),
                  static_cast<double>(acc.denominator_events));
}

}  // namespace effstuff
