#pragma once

#include "eff_math_utils.h"
#include "hms_pid_common.h"

namespace effstuff {

struct HMSCalEffTagCerAccumulator {
  long long denominator_events = 0;
  long long numerator_events = 0;
};

inline void accumulate_HMS_cal_eff_tag_cer(HMSCalEffTagCerAccumulator& acc,
                                           const HMSEventContext& evt,
                                           const HMSPidCuts& cuts) {
  // CAL tag/probe efficiency: use CER-tagged electron candidates as the
  // denominator, then ask what fraction also pass the CAL cut.
  if (!hms_pid_denominator_pass(evt, cuts)) {
    return;
  }

  if (!hms_cer_pass(evt, cuts)) {
    return;
  }

  ++acc.denominator_events;
  if (hms_cal_pass(evt, cuts)) {
    ++acc.numerator_events;
  }
}

inline double HMS_cal_eff_tag_cer(const HMSCalEffTagCerAccumulator& acc) {
  return safe_div(static_cast<double>(acc.numerator_events),
                  static_cast<double>(acc.denominator_events));
}

}  // namespace effstuff
