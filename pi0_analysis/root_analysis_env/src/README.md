# Source Macros

The `src/` directory contains the core ROOT macros and helper headers for
per-run NPS pi0 analysis.

## Primary Entry Points

- `nps_analysis.C`: main analysis macro for skim-style input files.
- `nps_analysis_wfpi0.C`: waveform-calibrated production-data variant.
- `skim_data.C`: skim-building macro.
- `utils.C` and `utils.h`: shared utility helpers.

## Helper Headers

- `nps_helper.h`: common analysis helper functions.
- `nps_physics_var.h`: physics-variable calculations.
- `nps_comb_bg.h`: combinatorial-background handling.
- `nps_time_bg.h`: timing-background handling.
- `nps_mmiss_cor.h`: missing-mass corrections.
- `nps_charge.h`: charge/current utilities.
- `nps_livetimes.h`: livetime calculations.
- `nps_target_boiling.h`: target-boiling correction helpers.
- `nps_hms_track_eff.h`: HMS tracking-efficiency helpers.
- `nps_get_prescale.h`: prescale lookup utilities.
- `nps_yao_database_reader.h`: database reader helpers.
- `nps_pi0_curr_norm_rate.h`: current-normalized pi0 rate helpers.

## Runtime Configuration

The main analysis macros are usually driven by wrapper scripts, but they also
read these environment variables:

- `NPS_SKIM_DIR`
- `NPS_OUTPUT_DIR`
- `NPS_RUNLIST`
- `NPS_EBEAM`

Preserve these names when refactoring; the top-level `run_*.sh` scripts depend
on them.

## Output Contract

The batch wrappers expect the macro output to include a line containing:

```text
[CSV_WRITTEN] <path-to-per-run-csv>
```

If that marker changes, update the shell wrappers at the same time.

## Development Notes

- Prefer changing helper headers when logic is shared by both `nps_analysis.C`
  and `nps_analysis_wfpi0.C`.
- Be careful with branch-name changes. The original and waveform-calibrated
  data use different tree and branch names.
- Avoid writing analysis artifacts into `src/`; generated files belong under
  `output/` or another explicit output directory.
