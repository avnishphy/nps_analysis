#pragma once

#include "eff_math_utils.h"
#include "hms_pid_common.h"

namespace effstuff {

struct HMSCerEffTagCalAccumulator {
  long long denominator_events = 0;
  long long numerator_events = 0;
};

inline void accumulate_HMS_cer_eff_tag_cal(HMSCerEffTagCalAccumulator& acc,
                                           const HMSEventContext& evt,
                                           const HMSPidCuts& cuts) {
  // CER tag/probe efficiency: use CAL-tagged electron candidates as the
  // denominator, then ask what fraction also pass the CER cut.
  if (!hms_pid_denominator_pass(evt, cuts)) {
    return;
  }

  if (!hms_cal_pass(evt, cuts)) {
    return;
  }

  ++acc.denominator_events;
  if (hms_cer_pass(evt, cuts)) {
    ++acc.numerator_events;
  }
}

inline double HMS_cer_eff_tag_cal(const HMSCerEffTagCalAccumulator& acc) {
  return safe_div(static_cast<double>(acc.numerator_events),
                  static_cast<double>(acc.denominator_events));
}

}  // namespace effstuff
