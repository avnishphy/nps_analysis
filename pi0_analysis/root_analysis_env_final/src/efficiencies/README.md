# NPS HMS Efficiency and Livetime Pipeline

This directory contains a modular C++/ROOT pipeline to compute run-level HMS efficiencies and livetime observables for NPS HMS coincidence replay files. The implementation is centered on `compute_efficiencies_stuff.cxx` and is designed to be maintainable, auditable, and easy to extend.

## 1. What This Pipeline Computes

For each selected run (grouped by kinematic setting `Kin_old`), the code computes:

- HEL charge before good-event cuts
- HEL charge after good-event cuts
- Accepted event ranges used in all downstream event/scaler filtering
- HMS_pid_eff
- HMS_cal_eff_tag_cer
- HMS_cer_eff_tag_cal
- HMS_tracking_eff
- HMS_hodo_3of4_eff
- NewGen_EDTM_livetime
- Uncertainties for all efficiency observables and NewGen EDTM livetime
- prescale_token
- ps_factor
- which_TRIG
- beam_time
- file source (`updated` preferred over `production`)
- segment count found
- run-processing status

Outputs:

- One CSV per processed kinematic setting
- One global summary CSV aggregating counts by kinematic setting

## 2. Data Inputs and Paths

Default paths in `compute_efficiencies_stuff.cxx`:

- Config CSV:
  - `/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/config/nps_dvcs_all_kins_main.csv`
- Updated replay directory:
  - `/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated`
- Production replay directory:
  - `/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/production`
- Output directory:
  - `/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/efficiency_stuff`

Replay file naming convention:

- `nps_hms_coin_<run>_<segment>_1_-1.root`

If a run exists in both `updated` and `production`, the code always uses `updated`.

## 3. Directory Contents (Core Files)

### Main driver

- `compute_efficiencies_stuff.cxx`
  - Program entrypoint
  - CLI parsing and kinematic selection (interactive or explicit CLI)
  - Per-run orchestration across all segments
  - CSV writing per kinematic setting and global summary writing

### Selection and filtering

- `NPS_selection_helper.h`
  - Canonical hel selection engine used by this pipeline
- `good_event_selection_helper.h`
  - Wraps `NPS_selection_helper.h`
  - Produces `GoodSelectionSummary` with:
    - cut strings
    - accepted ranges
    - HEL charge before/after cut
    - current-window limits

### Config parsing

- `config_csv_helper.h`
  - CSV parser with quoted-field support
  - Required headers:
    - `run_number`, `Kin_old`, `Type`, `prescale`
  - Row validation and malformed-row accounting
  - Type filtering helpers

### File discovery

- `root_file_discovery.h`
  - Finds all segments for a run
  - Sorts by parsed segment number
  - Implements updated-over-production preference

### Prescale and beam-time

- `prescale_beamtime_helper.h`
  - Parses prescale token and trigger number
  - Computes `ps_factor`
  - Accumulates beam time and scaler EDTM denominator from `TSH`

### HMS efficiencies (modular)

- `hms_event_context.h`
  - Shared event-level context object
- `hms_pid_common.h`
  - Shared PID cuts and denominator logic
- `hms_pid_eff.h`
- `hms_cal_eff_tag_cer.h`
- `hms_cer_eff_tag_cal.h`
- `hms_tracking_eff.h`
- `hms_hodo_3of4_eff.h`

### Livetime module

- `newgen_edtm_livetime.h`
  - `NewGen_EDTM_livetime = numerator / denominator * ps_factor`
  - propagated uncertainty helper for `NewGen_EDTM_livetime_err`

### Shared utility and data types

- `eff_math_utils.h`
  - `safe_div`, trim/lower helpers, binomial uncertainty helper
- `efficiency_types.h`
  - Shared structs for run rows, summaries, ranges, etc.

### Parallel launcher

- `run_efficiencies_stuff_parallel.sh`
  - Interactive kinematic prompt (if not specified via CLI)
  - Per-run parallel dispatch using `xargs -P`
  - Worker cap: `min(--jobs, total_jobs)`
  - Real-time progress monitor with spinner + started/running/done counts
  - Safe merge of per-run partial outputs into final per-kinematic CSVs

### Reference legacy files (for comparison/traceability)

- `compute_luminosity_scaler.cxx`
- `compute_luminosity_scaler_josh.cxx`
- `nps_hel_select.h` (retained legacy implementation; not used by the active pipeline)
- `nps_hel_good_events.C` (standalone diagnostic/snapshot wrapper using `NPS_selection_helper.h`)

## 4. End-to-End Algorithm

For each selected `Kin_old`:

1. Read config CSV and filter rows by kinematic setting.
2. Apply run type filter (default: `production`, `Production`; configurable).
3. For each selected run:
   1. Locate all segments in `updated` first, else `production`.
   2. For each segment:
      1. Build good-event selection using the canonical hel logic from `NPS_selection_helper.h`.
      2. Store accepted `evcount` and `g.evnum` ranges.
      3. Accumulate HEL charge before/after cuts.
      4. Accumulate beam time and EDTM scaler diagnostics from `TSH`, gated by:
         - current window from selection
         - accepted `evcount` ranges from the same good-event selection
      5. Accumulate the production EDTM scaler denominator from current-window `TSH` intervals.
      6. Collect current-window EDTM TDC values to estimate the EDTM peak, plus a good-event-gated diagnostic subset.
      7. Accumulate HMS efficiency numerators/denominators from selected events.
   3. Determine EDTM peak from current-window EDTM TDC values via histogram-mode binning (`bin width = 10`).
   4. Count cached current-window EDTM values within the peak window for the NewGen numerator.
   5. Compute final run-level observables and status.
4. Write per-kinematic CSV.
5. Write global summary CSV.

## 5. Good-Event Selection: Exact Defaults and Impact

Selection settings from `make_default_selection_settings()`:

- `helicity_mode = QuartetPM`
- `junk_floor_uA = 2.5`
- `mean_current_calc_min = 2.5`
- `mean_current_min = 2.75`
- `use_current_window = true`
- `i0_mode = Peak`
- `window_frac = 0.15`
- `stable_window_N = 30`
- `max_charge_frac_spread = 0.08`

Branch mapping defaults:

- `H.BCM4A_Hel.scalerCurrent`
- `H.BCM4A_Hel.scalerCharge`
- `actualHelicity`
- `evcount`
- `evNumber`
- `g.evnum`

The accepted ranges generated by this selection are the core gate for downstream event calculations. This is the key mechanism ensuring consistency from good-event determination to efficiency/livetime outputs.

## 6. Event Selection Applied in Efficiency Loops

An event is considered selected iff all are true:

- `round(fEvtHdr.fEvtType) == 1`
- `g.evnum` lies in accepted `g.evnum` ranges

This selected-flag then feeds each metric module.
NewGen EDTM uses the event-by-event current window only for the production livetime. The selected-event and accepted-`evcount` EDTM quantities are written as diagnostics.

## 7. Exact Metric Definitions

## 7.1 PID common denominator preconditions

From `hms_pid_common.h`:

- Selected event (as above)
- `edtm_tdc < 0.1`
- Optional track requirement (`H.dc.ntrack > 0.5`, enabled by default)
- Reconstructed-track acceptance cuts:
  - `|H.gtr.dp| <= 10.0`
  - `|H.gtr.th| <= 0.1`
  - `|H.gtr.ph| <= 0.04`
  - `|H.react.z| <= 8.0`

PID thresholds:

- CER pass: `H.cer.npeSum > 4.0`
- CAL pass: `0.8 < H.cal.etotnorm <= 1.2`

## 7.2 HMS_pid_eff

- Denominator: PID-common denominator events
- Numerator: denominator events passing both CER and CAL

Formula:

- `HMS_pid_eff = N(CER && CAL) / N(PID denominator)`
- `HMS_pid_eff_err = sqrt( p * (1 - p) / N_den )`

## 7.3 HMS_cal_eff_tag_cer

- Denominator: PID-common denominator events that pass CER
- Numerator: denominator events that also pass CAL

Formula:

- `HMS_cal_eff_tag_cer = N(CAL && CER) / N(CER-tag denominator)`
- `HMS_cal_eff_tag_cer_err = sqrt( p * (1 - p) / N_den )`

## 7.4 HMS_cer_eff_tag_cal

- Denominator: PID-common denominator events that pass CAL
- Numerator: denominator events that also pass CER

Formula:

- `HMS_cer_eff_tag_cal = N(CER && CAL) / N(CAL-tag denominator)`
- `HMS_cer_eff_tag_cal_err = sqrt( p * (1 - p) / N_den )`

## 7.5 HMS_tracking_eff

Tracking “should” condition:

- selected event
- `edtm_tdc < 0.1`
- `round(H.hod.goodscinhit) == 1`
- `0.8 < H.cal.etotnorm <= 1.2`
- `H.cer.npeSum > 4.0`
- `0.8 < H.hod.betanotrack < 1.2`

Tracking “did” condition:

- should-condition and `H.dc.ntrack > 0.5`

Formula:

- `HMS_tracking_eff = N(did) / N(should)`
- `HMS_tracking_eff_err = sqrt( p * (1 - p) / N_should )`

## 7.6 HMS_hodo_3of4_eff

Denominator condition:

- selected event
- `edtm_tdc < 0.1`
- `H.dc.ntrack > 0.5`
- reconstructed-track acceptance cuts listed in Section 7.1, including
  `|H.gtr.dp| <= 10.0` and `|H.react.z| <= 8.0`

Per-plane good-hit condition:

- each of `H.hod.{1x,1y,2x,2y}.nhits` in `(0, 3)`

Numerator:

- denominator event with at least 3 good planes

Formula:

- `HMS_hodo_3of4_eff = N(good_planes >= 3) / N(denominator)`
- `HMS_hodo_3of4_eff_err = sqrt( p * (1 - p) / N_den )`

## 7.7 NewGen_EDTM_livetime

Numerator:

- current-window events with finite `T.hms.hEDTM_tdcTimeRaw > 1`
- if EDTM peak is known: require `|edtm_tdc - peak| <= 500`
- if EDTM peak unavailable: numerator is 0

Denominator:

- summed positive scaler increments of `H.EDTM.scaler` in `TSH`
- gated by the current window

Final formula:

- `NewGen_EDTM_livetime = (N_edtm_current_window / Scaler_EDTM_current_window) * ps_factor`

Diagnostic columns:

- `NewGen_EDTM_num_good_event_gated`: current-window EDTM numerator after the event-level good selection
- `NewGen_EDTM_den_evcount_gated_value`: current-window EDTM scaler denominator after the `TSH evcount` gate
- these gated quantities are not used in the production livetime because EDTM pulses are artificial livetime probes rather than selected physics events

Uncertainty model used in code (first-order propagation, Poisson-like counts):

- Let $R = ps \cdot N / D$
- $\sigma_N = \sqrt{N}$
- $\sigma_D = \sqrt{D}$
- $\sigma_R = \sqrt{(\partial R/\partial N)^2 \sigma_N^2 + (\partial R/\partial D)^2 \sigma_D^2}$
- Implemented as:
  - `term_n = ((ps / D) * sigma_n)^2`
  - `term_d = ((ps * N / D^2) * sigma_d)^2`
  - `NewGen_EDTM_livetime_err = sqrt(term_n + term_d)`

## 8. Prescale and Trigger Bookkeeping

From `prescale_beamtime_helper.h`:

- Token parsing supports one or several `psN=setting` assignments.
- `setting = -1` means disabled.
- Example: `ps6=-1,ps3=1` selects PS3, factor 2, and `H.hTRIG3.scaler`.
- Temporary multi-trigger policy: use the first enabled assignment in token
  order, set `prescale_multiple_enabled=1`, warn, and list the run in the plot
  audit CSV for manual review.
- `trig_number` and `ps_factor` come from the same enabled assignment.
- `which_TRIG` formed as `H.hTRIG<TrigNumber>.scaler`

Prescale factor rule:

- parse integer `r` from token
- if parse fails or `r <= 0`, fallback to `1.0`
- else:
  - `ps_factor = 2^(r-1) + 1`

## 9. Beam Time Accumulation

From `accumulate_beam_time_and_scaler()`:

- Tree: `TSH`
- Required branches:
  - `H.1MHz.scalerTime`
- `H.BCM4A.scalerCurrent`
- `H.EDTM.scaler`
- `evcount`
Per scaler step:

- `dt = t_i - t_{i-1}`
- reject negative/non-finite `dt`, count as `negative_dt_intervals`
- `dEDTM = max(0, EDTM_i - EDTM_{i-1})`
- count negative raw deltas as `non_monotonic_scaler_steps`
- accumulate a time-weighted `H.S1X.scalerRate` mean and RMS over the same
  current-window and accepted-`evcount` exposure used for beam time

If in current-window:

- `scaler_edtm_total += dEDTM`

If also in accepted `evcount` range:

- `beam_time += dt`
- `scaler_edtm_total_evcount_gated += dEDTM`

## 10. Output Files and Columns

Per-kinematic CSV:

- `efficiency_<sanitized_kin>.csv`

Columns:

- `run_number`
- `kinematic_setting`
- `run_type`
- `HEL_charge_before_cut_uC`
- `HEL_neg_charge_before_cut_uC`
- `HEL_pos_charge_before_cut_uC`
- `HEL_charge_after_cut_uC`
- `HEL_neg_charge_after_cut_uC`
- `HEL_pos_charge_after_cut_uC`
- `HMS_pid_eff`
- `HMS_cal_eff_tag_cer`
- `HMS_cer_eff_tag_cal`
- `HMS_tracking_eff`
- `HMS_hodo_3of4_eff`
- `NewGen_EDTM_livetime`
- `NewGen_EDTM_num`
- `NewGen_EDTM_den`
- `NewGen_EDTM_peak`
- `NewGen_EDTM_num_good_event_gated`
- `NewGen_EDTM_den_evcount_gated_value`
- `NewGen_EDTM_den_evcount_gated`
- `prescale_token`
- `ps_factor`
- `which_TRIG`
- `beam_time`
- `file_source_used`
- `segment_count_found`
- `n_segments`
- `run_processing_status`
- `missing_branch_segments`
- `selection_failed_segments`
- `negative_dt_intervals`
- `non_monotonic_scaler_steps`
- `HMS_pid_eff_err`
- `HMS_cal_eff_tag_cer_err`
- `HMS_cer_eff_tag_cal_err`
- `HMS_tracking_eff_err`
- `HMS_hodo_3of4_eff_err`
- `NewGen_EDTM_livetime_err`

Per-kinematic selection report CSV:

- `selection_report_<sanitized_kin>.csv`

This file contains one row per run segment and carries the detailed selection diagnostics previously embedded in the run table, including:

- selection success/message
- helicity mode and quartet snap flags
- current-window and scaler stats
- hel0/hel-/hel+ charges (before and after cut)
- ready-to-use `evcount`, `evNumber`, and `g.evnum` cut strings
- accepted evcount/evNumber/g.evnum ranges

## 11. Efficiency Overlay Plots

Use `scripts/plot_efficiency_overlays.py` to make publication-quality summary
overlays from the efficiency CSVs:

```bash
python3 scripts/plot_efficiency_overlays.py \
  --input-dir output/efficiency_stuff
```

The script writes PNG and PDF plots to `output/efficiency_stuff/plots`:

- one multipanel compact-axis plot vs run for each `efficiency_*.csv`
- one multipanel plot vs current for each `efficiency_*.csv`
- one multipanel plot vs time-weighted HMS S1X rate for each `efficiency_*.csv`
- one combined multipanel plot vs current across all discovered efficiency CSVs
- one combined multipanel plot vs HMS S1X rate across all discovered efficiency CSVs
- `efficiency_multipanel_multiple_enabled_trigger_runs.csv`, listing plotted runs
  with more than one enabled `psN` assignment (header-only when none are found)
- run overlays label prescale trigger groups such as `ps4` and `ps6`; exact values like `ps6=0` and `ps6=3` are intentionally grouped together
- extreme low outliers are clipped at the panel edge and listed with run number/value so they do not collapse the useful y-axis range

Global summary CSV:

- `efficiency_runs_processed.csv`

Columns:

- `kinematic_setting`
- `total_runs_listed`
- `production_runs`
- `runs_selected_by_type`
- `runs_processed`
- `runs_not_found`
- `malformed_rows_skipped`

## 11. Run Status Definitions

- `missing_root_files`
  - no run files found in either updated or production directory
- `selection_failed`
  - no segment survived good-event selection stage
- `processed_partial`
  - run processed, but some segments had missing branches or failed selection
- `processed`
  - run processed without those partial flags

## 12. Building and Running

## 12.1 Build executable

From this directory:

```bash
g++ -O3 -std=c++17 -Wall -Wextra -pedantic `root-config --cflags` \
  compute_efficiencies_stuff.cxx -o compute_efficiencies_stuff `root-config --libs`
```

If your shell is `tcsh`, the same line works as shown above.

## 12.2 Run executable directly

Interactive kinematic selection:

```bash
./compute_efficiencies_stuff
```

Non-interactive examples:

```bash
./compute_efficiencies_stuff --kin KinC_x60_4b --types production,Production
./compute_efficiencies_stuff --all-kins --types production,Production
./compute_efficiencies_stuff --kin KinC_x60_4b --run 12345 --no-interactive
```

## 12.3 Run parallel launcher

Interactive kinematic prompt:

```bash
./run_efficiencies_stuff_parallel.sh
```

Non-interactive examples:

```bash
./run_efficiencies_stuff_parallel.sh --kin KinC_x60_4b
./run_efficiencies_stuff_parallel.sh --all-kins --jobs 64
./run_efficiencies_stuff_parallel.sh --kin KinC_x60_4b --run 12345 --jobs 8
```

## 13. Parallel Execution Behavior

The launcher:

- Builds job list as unique `(Kin_old, run_number)` pairs after type filtering.
- Caps worker count to number of queued jobs.
- Runs one job per run with output isolated in temporary directories.
- Merges partial CSVs after all workers finish.
- Provides live progress with:
  - progress bar
  - spinner (heartbeat)
  - started/running/done counters

This gives immediate visibility that workers are active even when no jobs have completed yet.

## 14. Error Handling and Diagnostics

Explicit failure cases covered by code and/or output status:

- Missing config CSV
- Missing required CSV headers
- Malformed CSV rows (counted and surfaced in summary)
- Invalid kinematic selection
- Unsupported run types relative to current filter
- Missing ROOT files
- Missing required branches in segment trees
- Selection that yields no usable segments

Per-run logs in parallel mode are written as:

- `<tmp>/jobs/<sanitized_kin>/<run>/run.log`

## 15. How to Extend the Pipeline Safely

To add a new metric:

1. Create a dedicated header with an accumulator + final function (follow existing pattern).
2. Add fields to `RunProcessingRow` in `efficiency_types.h`.
3. Extend event loop context in `hms_event_context.h` only if new branches are needed.
4. Add accumulation calls in `accumulate_segment_event_metrics()`.
5. Add output column(s) in `write_kinematic_csv()`.
6. Keep all new calculations gated by the same selected-event logic to preserve consistency.

Recommended validation strategy:

- Compare outputs for a small run subset against `compute_luminosity_scaler*.cxx` references.
- Inspect partial statuses and per-run logs before full production launch.

## 16. Notes for New Users

If you are starting from scratch, this is the fastest path:

1. Build `compute_efficiencies_stuff`.
2. Run one kinematic setting with one run to verify environment/paths.
3. Run one full kinematic setting.
4. Launch parallel mode for larger sweeps.
5. Review per-kinematic CSV + global summary.

This directory intentionally separates selection, parsing, discovery, metric modules, and execution orchestration so you can understand and modify each piece independently without destabilizing the rest of the workflow.
