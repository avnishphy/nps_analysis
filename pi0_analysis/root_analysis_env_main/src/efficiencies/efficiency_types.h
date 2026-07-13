#pragma once

#include <limits>
#include <string>
#include <vector>

namespace effstuff {

struct EventRange {
  long long lo = 0;
  long long hi = 0;
};

struct GoodSelectionSummary {
  bool ok = false;
  std::string message;

  std::string helicity_mode;
  bool quartet_snap_applied = false;

  std::string evcount_cut;
  std::string evnumber_cut;
  std::string gevnum_cut;

  std::vector<EventRange> accepted_evcount_ranges;
  std::vector<EventRange> accepted_evnumber_ranges;
  std::vector<EventRange> accepted_gevnum_ranges;

  double hel_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel0_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_neg_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_pos_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel0_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_neg_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_pos_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();

  double mean_current_uA = std::numeric_limits<double>::quiet_NaN();
  double i0_used_uA = std::numeric_limits<double>::quiet_NaN();
  int n_scaler_reads_pre = 0;
  int n_scaler_reads_post = 0;

  double current_min_uA = std::numeric_limits<double>::quiet_NaN();
  double current_max_uA = std::numeric_limits<double>::quiet_NaN();
  bool current_window_enabled = false;
};

struct RunConfigRow {
  int run_number = 0;
  std::string kin_old;
  std::string run_type;
  std::string prescale_token;
};

struct LocatedRunFiles {
  std::vector<std::string> files;
  std::string source;
  int segment_count = 0;
};

struct RunProcessingRow {
  int run_number = 0;
  std::string kinematic_setting;
  std::string run_type;

  double hel_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_neg_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_pos_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_neg_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_pos_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();

  double HMS_pid_eff = std::numeric_limits<double>::quiet_NaN();
  double HMS_pid_eff_err = std::numeric_limits<double>::quiet_NaN();
  double HMS_cal_eff_tag_cer = std::numeric_limits<double>::quiet_NaN();
  double HMS_cal_eff_tag_cer_err = std::numeric_limits<double>::quiet_NaN();
  double HMS_cer_eff_tag_cal = std::numeric_limits<double>::quiet_NaN();
  double HMS_cer_eff_tag_cal_err = std::numeric_limits<double>::quiet_NaN();
  double HMS_tracking_eff = std::numeric_limits<double>::quiet_NaN();
  double HMS_tracking_eff_err = std::numeric_limits<double>::quiet_NaN();
  double HMS_hodo_3of4_eff = std::numeric_limits<double>::quiet_NaN();
  double HMS_hodo_3of4_eff_err = std::numeric_limits<double>::quiet_NaN();

  double NewGen_EDTM_livetime = std::numeric_limits<double>::quiet_NaN();
  double NewGen_EDTM_livetime_err = std::numeric_limits<double>::quiet_NaN();
  long long NewGen_EDTM_num = 0;
  double NewGen_EDTM_den = std::numeric_limits<double>::quiet_NaN();
  double NewGen_EDTM_peak = std::numeric_limits<double>::quiet_NaN();
  long long NewGen_EDTM_num_good_event_gated = 0;
  double NewGen_EDTM_den_evcount_gated_value = std::numeric_limits<double>::quiet_NaN();
  bool NewGen_EDTM_den_evcount_gated = false;

  std::string prescale_token;
  double ps_factor = std::numeric_limits<double>::quiet_NaN();
  std::string which_TRIG;
  double beam_time = std::numeric_limits<double>::quiet_NaN();

  std::string file_source_used;
  int segment_count_found = 0;
  int n_segments = 0;
  std::string run_processing_status;

  int missing_branch_segments = 0;
  int selection_failed_segments = 0;
  int negative_dt_intervals = 0;
  int non_monotonic_scaler_steps = 0;
};

struct SelectionReportRow {
  int run_number = 0;
  std::string kinematic_setting;
  std::string run_type;
  int segment_number = -1;
  std::string segment_file;
  std::string file_source_used;
  bool selection_ok = false;
  std::string selection_message;

  std::string helicity_mode;
  bool quartet_snap_applied = false;
  bool current_window_enabled = false;
  double mean_current_uA = std::numeric_limits<double>::quiet_NaN();
  double i0_used_uA = std::numeric_limits<double>::quiet_NaN();
  double current_min_uA = std::numeric_limits<double>::quiet_NaN();
  double current_max_uA = std::numeric_limits<double>::quiet_NaN();
  int n_scaler_reads_pre = 0;
  int n_scaler_reads_post = 0;

  double hel0_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_neg_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_pos_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_charge_before_cut_uC = std::numeric_limits<double>::quiet_NaN();

  double hel0_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_neg_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_pos_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();
  double hel_charge_after_cut_uC = std::numeric_limits<double>::quiet_NaN();

  std::string evcount_cut;
  std::string evnumber_cut;
  std::string gevnum_cut;

  std::string accepted_evcount_ranges;
  std::string accepted_evnumber_ranges;
  std::string accepted_gevnum_ranges;
};

struct KinematicProcessingSummary {
  std::string kinematic_setting;
  int total_runs_listed = 0;
  int production_runs = 0;
  int runs_selected_by_type = 0;
  int runs_processed = 0;
  int runs_not_found = 0;
  int malformed_rows_skipped = 0;
};

}  // namespace effstuff
