#pragma once

#include <cmath>
#include <limits>

#include "hms_event_context.h"

namespace effstuff {

// Centralized HMS efficiency cut values. Despite the historical HMSPidCuts
// name, this struct owns every HMS detector threshold used by these efficiency
// headers so the PID, tracking, and hodoscope definitions cannot drift apart.
struct HMSPidCuts {
  // Electron PID tag/probe cuts.
  double cer_npe_min = 4.0;
  double cal_etotnorm_min = 0.8;
  double cal_etotnorm_max = 1.2;

  // Shared event-quality cuts used to reject EDTM pulser events and require HMS
  // tracking where a metric denominator intentionally calls for a track.
  double edtm_tdc_max = 0.1;
  double hdc_ntrack_min = 0.5;

  // HMS acceptance cuts for PID-style denominators.
  double hms_dp_abs_max = 15.0;
  double hms_th_abs_max = 0.1;
  double hms_ph_abs_max = 0.04;

  // Tracking-efficiency should-track tag. H.hod.betanotrack is deliberately
  // independent of the drift-chamber track count being measured.
  double hod_notrack_min = 0.8;
  double hod_notrack_max = 1.2;
  long long hod_goodscinhit_value = 1;

  // Hodoscope 3-of-4 plane-hit definition.
  double hod_plane_nhits_min = 0.0;
  double hod_plane_nhits_max = 3.0;
  int hod_min_good_planes = 3;

  // PID efficiency denominators normally require an HMS track. This switch is
  // left here for controlled studies, but production defaults to true.
  bool require_track = true;
};

inline bool hms_edtm_pass(const HMSEventContext& evt,
                          const HMSPidCuts& cuts) {
  return std::isfinite(evt.edtm_tdc) && (evt.edtm_tdc < cuts.edtm_tdc_max);
}

inline bool hms_track_pass(const HMSEventContext& evt,
                           const HMSPidCuts& cuts) {
  return std::isfinite(evt.hdc_ntrack) && (evt.hdc_ntrack > cuts.hdc_ntrack_min);
}

inline bool hms_cer_pass(const HMSEventContext& evt,
                         const HMSPidCuts& cuts) {
  return std::isfinite(evt.cer_npe_sum) && (evt.cer_npe_sum > cuts.cer_npe_min);
}

inline bool hms_cal_pass(const HMSEventContext& evt,
                         const HMSPidCuts& cuts) {
  return std::isfinite(evt.cal_etotnorm) &&
         (evt.cal_etotnorm > cuts.cal_etotnorm_min) &&
         (evt.cal_etotnorm <= cuts.cal_etotnorm_max);
}

inline bool hms_pid_kinematics_pass(const HMSEventContext& evt,
                                    const HMSPidCuts& cuts) {
  return std::isfinite(evt.hms_dp) && std::isfinite(evt.hms_th) && std::isfinite(evt.hms_ph) &&
         (std::fabs(evt.hms_dp) <= cuts.hms_dp_abs_max) &&
         (std::fabs(evt.hms_th) <= cuts.hms_th_abs_max) &&
         (std::fabs(evt.hms_ph) <= cuts.hms_ph_abs_max);
}

inline bool hms_good_scin_pass(const HMSEventContext& evt,
                               const HMSPidCuts& cuts) {
  return std::isfinite(evt.hod_goodscinhit) &&
         (std::llround(evt.hod_goodscinhit) == cuts.hod_goodscinhit_value);
}

inline bool hms_hod_notrack_pass(const HMSEventContext& evt,
                                 const HMSPidCuts& cuts) {
  return std::isfinite(evt.hod_notrack) &&
         (evt.hod_notrack > cuts.hod_notrack_min) &&
         (evt.hod_notrack < cuts.hod_notrack_max);
}

inline bool hms_hodo_plane_pass(double nhits,
                                const HMSPidCuts& cuts) {
  return std::isfinite(nhits) &&
         (nhits > cuts.hod_plane_nhits_min) &&
         (nhits < cuts.hod_plane_nhits_max);
}

inline bool hms_pid_denominator_pass(const HMSEventContext& evt,
                                     const HMSPidCuts& cuts) {
  // Common denominator for HMS PID and tag/probe PID efficiencies:
  // selected physics event, non-EDTM, optional HMS track, and HMS acceptance.
  const bool track_ok = (!cuts.require_track) || hms_track_pass(evt, cuts);
  return evt.selected && hms_edtm_pass(evt, cuts) && track_ok &&
         hms_pid_kinematics_pass(evt, cuts);
}

}  // namespace effstuff
