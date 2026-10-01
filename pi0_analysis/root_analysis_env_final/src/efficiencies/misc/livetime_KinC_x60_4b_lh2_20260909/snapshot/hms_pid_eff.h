#pragma once

#include "eff_math_utils.h"
#include "hms_pid_common.h"
namespace effstuff {

struct HMSPidEffAccumulator {
  long long denominator_events = 0;
  long long numerator_events = 0;
};

inline void accumulate_HMS_pid_eff(HMSPidEffAccumulator& acc,
                                   const HMSEventContext& evt,
                                   const HMSPidCuts& cuts) {
  // Inclusive HMS PID efficiency: among accepted, non-EDTM, tracked HMS events,
  // count the fraction passing both CER and CAL electron-identification cuts.
  if (!hms_pid_denominator_pass(evt, cuts)) {
    return;
  }

  ++acc.denominator_events;
  if (hms_cer_pass(evt, cuts) && hms_cal_pass(evt, cuts)) {
    ++acc.numerator_events;
  }
}

inline double HMS_pid_eff(const HMSPidEffAccumulator& acc) {
  return safe_div(static_cast<double>(acc.numerator_events),
                  static_cast<double>(acc.denominator_events));
}

}  // namespace effstuff
