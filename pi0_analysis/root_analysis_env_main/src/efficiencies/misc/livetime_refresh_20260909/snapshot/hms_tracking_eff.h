#pragma once

#include <cmath>
#include <limits>

#include "eff_math_utils.h"
#include "hms_pid_common.h"

namespace effstuff {

struct HMSTrackingEffAccumulator {
  long long should_events = 0;
  long long did_events = 0;
};

// Tracking efficiency is measured with a no-track hodoscope beta tag. The
// denominator asks whether the HMS should have reconstructed a DC track without
// using H.dc.ntrack; the numerator then checks whether a track was found.
inline void accumulate_HMS_tracking_eff(HMSTrackingEffAccumulator& acc,
                                        const HMSEventContext& evt,
                                        const HMSPidCuts& cuts) {
  const bool electron_pid_should =
      hms_cal_pass(evt, cuts) &&
      hms_cer_pass(evt, cuts) &&
      hms_good_scin_pass(evt, cuts) &&
      hms_hod_notrack_pass(evt, cuts);

  if (!(evt.selected && hms_edtm_pass(evt, cuts) && electron_pid_should)) {
    return;
  }

  ++acc.should_events;
  if (hms_track_pass(evt, cuts)) {
    ++acc.did_events;
  }
}

inline double HMS_tracking_eff(const HMSTrackingEffAccumulator& acc) {
  return safe_div(static_cast<double>(acc.did_events),
                  static_cast<double>(acc.should_events));
}

inline double HMS_tracking_eff_err(const HMSTrackingEffAccumulator& acc) {
  if (acc.should_events <= 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double p = HMS_tracking_eff(acc);
  return std::sqrt(std::max(0.0, p * (1.0 - p)) /
                   static_cast<double>(acc.should_events));
}

}  // namespace effstuff
