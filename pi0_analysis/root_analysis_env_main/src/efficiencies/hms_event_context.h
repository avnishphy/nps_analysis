#pragma once

#include <limits>

namespace effstuff {

struct HMSEventContext {
  bool selected = false;

  double edtm_tdc = std::numeric_limits<double>::quiet_NaN();
  double cer_npe_sum = std::numeric_limits<double>::quiet_NaN();
  double cal_etotnorm = std::numeric_limits<double>::quiet_NaN();
  double hdc_ntrack = std::numeric_limits<double>::quiet_NaN();

  double hms_dp = std::numeric_limits<double>::quiet_NaN();
  double hms_th = std::numeric_limits<double>::quiet_NaN();
  double hms_ph = std::numeric_limits<double>::quiet_NaN();

  double hod_goodscinhit = std::numeric_limits<double>::quiet_NaN();
  double hod_notrack = std::numeric_limits<double>::quiet_NaN();

  double hod_1x_nhits = std::numeric_limits<double>::quiet_NaN();
  double hod_1y_nhits = std::numeric_limits<double>::quiet_NaN();
  double hod_2x_nhits = std::numeric_limits<double>::quiet_NaN();
  double hod_2y_nhits = std::numeric_limits<double>::quiet_NaN();
};

}  // namespace effstuff
