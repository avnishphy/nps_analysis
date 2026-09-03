# NPS pi0 Analysis Workflow

## Overview
This repository contains the Hall C NPS neutral-pion analysis workflow for data processing, diagnostic output generation, optional branch-combination, optional smearing, and optional xsec extraction.

The workflow is organized so that analysis remains physics-preserving while allowing staged automation.

## Current Scope and Constraints
- Efficiency computation is frozen for this phase.
- Validated physics behavior must be preserved.
- Geant4 internals are out of scope in this phase.
- Output paths and naming contracts are explicit and must not change silently.

## Repository Layout
- `config/`
  - `nps_dvcs_all_kins_main.csv`: master run and kinematic metadata table
  - `acceptance_cuts.conf`: acceptance and timing cut definitions
  - `nps_simulation_kinematics.csv`: canonical simulation-side kinematic table (planned/currently in progress)
- `src/analysis/`
  - `nps_analysis_main.C`: main per-run analysis macro
  - `run_parallel_nps_analysis_main.sh`: main data-stage driver (parallel)
  - `combine_analysis_branches.py`: per-kin combination stage
  - `acceptance_cuts.h/.cpp`: acceptance-cuts parser and runtime API
- `src/simulation_smearing/`
  - `run_smearing_pipeline.sh`: smearing pipeline driver
  - `nps_sim_smearing_new.C`: fit/map generation
  - `simc_pi0_analysis.C`: simulation producer with smearing map usage
- `src/xsec_extract/`
  - `run_xsec_pipeline.sh`: xsec stage driver
- `scripts/`
  - `run_pipeline.sh`: top-level orchestration entrypoint
  - `generate_simc_infiles.py`: CSV generator by default; `--gen_infile` generates SIMC infiles after review
- `output/`
  - canonical production outputs by kinematic setting
- `output_phase0_baseline/`
  - frozen baseline captures for regression comparison

## Environment Setup
Use the Hall C setup script before running stages:

```bash
source /group/nps/singhav/setup.csh
cd /work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main
```

Expected tools in environment:
- ROOT
- python3
- timeout

## End-to-End Workflow
1. Select run/kinematic rows from config.
2. Run per-run analysis macro over selected runs.
3. Rebuild per-kin summary CSV from per-run rows.
4. Optionally run combine stage.
5. Optionally run smearing stage.
6. Optionally run xsec stage.

## Primary Run Commands

### One-command pipeline entrypoint (Phase 2)
```bash
./scripts/run_pipeline.sh \
  --kin KinC_x60_4b \
  --source waveform \
  --jobs 3 \
  --gevnum-cut no \
  --types production,Production \
  --no-combine
```

What this entrypoint does:
- Generates or refreshes `config/nps_simulation_kinematics.csv` (unless `--no-generate-sim-config` is used).
- Initializes new simulation offset columns to `0.0`.
- Preserves existing calibrated offsets unless `--reset-offsets` is requested directly.
- Optionally generates per-kin SIMC infiles from the simulation CSV (`--generate-simc-infiles`).
- Optionally runs analysis-side SIMC/Geant4 wrapper commands from CSV + generated infiles (`--run-sim-chain`).
- Runs the existing analysis driver with forwarded arguments.

To skip CSV generation explicitly:
```bash
./scripts/run_pipeline.sh --no-generate-sim-config -- --kin KinC_x60_4b --source waveform --gevnum-cut no --jobs 3 --no-combine
```

To generate SIMC infiles from CSV values (including explicit offset columns):
```bash
./scripts/run_pipeline.sh \
  --generate-simc-infiles \
  --simc-template /u/group/nps/singhav/simc_gfortran_updated/infiles/nps_excl_pi0_x60_4b.inp \
  --simc-outdir config/simc_infiles \
  --simc-file-prefix nps_excl_pi0_ \
  -- --kin KinC_x60_4b --source waveform --gevnum-cut no --jobs 3 --no-combine
```

Notes:
- If `--kin` is forwarded to the analysis driver, SIMC infile generation mirrors that kinematic selection by default.
- Use `--simc-kin <Kin_old>` (repeatable) to explicitly control generated SIMC infile kinematics.

To run simulation wrappers from canonical CSV values (Phase 6 integration):
```bash
./scripts/run_pipeline.sh \
  --generate-simc-infiles \
  --run-sim-chain \
  --sim-chain-run-simc \
  --sim-chain-simc-cmd "run_simc --kin {kin_old} --input {simc_infile} --output {simc_output_file}" \
  --sim-chain-run-geant4 \
  --sim-chain-geant4-cmd "run_geant4 --kin {kin_old} --input {simc_output_file} --theta {nps_theta_deg} --distance {nps_target_distance_cm} --x-offset {geant4_nps_x_offset_cm} --y-offset {geant4_nps_y_offset_cm} --z-offset {geant4_nps_z_offset_cm} --outdir {geant4_output_dir}" \
  -- --kin KinC_x60_4b --source waveform --gevnum-cut no --jobs 3 --no-combine
```

Wrapper contract notes:
- Wrapper stage entrypoint: `scripts/run_simulation_chain.py`.
- No default SIMC/Geant4 executable command is baked into the wrapper; command templates must be provided explicitly.
- The wrapper writes a per-run manifest (default: `output/simulation_chain_manifest.csv`) that records commands, offsets, and output paths for provenance.

To build a diagnostics artifact index (Phase 7 hygiene support):
```bash
./scripts/run_pipeline.sh \
  --build-diagnostics-index \
  --diagnostics-index-out output/diagnostics/diagnostics_index.csv \
  -- --kin KinC_x60_4b --source waveform --gevnum-cut no --jobs 3 --no-combine
```

Diagnostics index notes:
- Index generator: `scripts/generate_diagnostics_index.py`.
- Default index output path: `<output-base>/diagnostics/diagnostics_index.csv`.
- Indexed columns: `kin`, `stage`, `artifact_type`, `run`, `relative_path`.

To generate explicit detector-side diagnostics coverage (HMS and NPS):
```bash
./scripts/run_pipeline.sh \
  --build-diagnostics-coverage \
  --diagnostics-coverage-out output/diagnostics/diagnostics_coverage_hms_nps.csv \
  -- --kin KinC_x60_4b --source waveform --gevnum-cut no --jobs 3 --no-combine
```

Coverage report notes:
- Coverage generator: `scripts/generate_diagnostics_coverage.py`.
- Default output path: `<output-base>/diagnostics/diagnostics_coverage_hms_nps.csv`.
- For each kin, the report records whether HMS-side and NPS-side diagnostics are both present based on summary metrics and diagnostic ROOT content.

### Representative single-run analysis stage
```bash
bash src/analysis/run_parallel_nps_analysis_main.sh \
  --config output_phase0_baseline/meta/nps_dvcs_all_kins_main.trim.csv \
  --kin KinC_x60_4b \
  --run 4253 \
  --source waveform \
  --gevnum-cut no \
  --jobs 1 \
  --types production,Production \
  --no-combine \
  --output-base output_phase0_baseline
```

### Representative 3-run batch analysis stage
```bash
bash src/analysis/run_parallel_nps_analysis_main.sh \
  --config output_phase0_baseline/meta/nps_dvcs_x60_4b_batch3.trim.csv \
  --kin KinC_x60_4b \
  --source waveform \
  --gevnum-cut no \
  --jobs 3 \
  --types production,Production \
  --no-combine \
  --output-base output_phase0_baseline
```

### Enable combine stage
```bash
bash src/analysis/run_parallel_nps_analysis_main.sh \
  --config <config_csv> \
  --kin <Kin_old> \
  --source waveform \
  --gevnum-cut no \
  --jobs <N> \
  --types production,Production \
  --combine-target LH2 \
  --output-base <output_base>
```

Note:
- Combine requires per-kin efficiency CSV availability under `<output_base>/efficiency_stuff`.

### Enable smearing stage
```bash
bash src/analysis/run_parallel_nps_analysis_main.sh \
  --config <config_csv> \
  --kin <Kin_old> \
  --source waveform \
  --gevnum-cut no \
  --jobs <N> \
  --types production,Production \
  --run-smearing \
  --smearing-target LH2 \
  --output-base <output_base>
```

Smearing reproducibility:
- `src/simulation_smearing/run_smearing_pipeline.sh` supports `--smear-seed <int>`.
- The same value is forwarded as `NPS_SMEAR_RANDOM_SEED` to `simc_pi0_analysis.C` so photon smearing draws are deterministic across reruns for fixed inputs.
- The fitter caches Gaussian pulls as `double` values without changing seed or draw order, precomputes immutable event quantities (`E_safe`, `log(E_safe)`, directions, normalization membership, and weight/`Nsmear`), and uses lightweight weighted-bin accumulators during optimization.
- Zero-position-smearing evaluations reuse event directions. Shared opening-angle kinematics computes `M_gg` and `(p_target+gamma gamma)^2` together; final ROOT diagnostics retain the established path.
- Sobol points are evaluated in batches of eight while preserving each candidate's event/smear order and objective. Landau or uncached evaluations automatically retain the scalar path.
- Every fit compares scalar-fast, batched-fast, and legacy ROOT-histogram objectives at three parameter points. A mismatch automatically selects the legacy evaluator for that fit.
- Pull caching has a 1 GiB process-wide budget shared by concurrent section fits; oversized fits retain the deterministic streaming-RNG path.
- Global-prefit Sobol batches and the eight retained MIGRAD starts run in parallel with independent Minuit2 instances. Per-section fits retain outer section parallelism, so their internal Sobol/MIGRAD work remains serial to avoid nested oversubscription.
- Production compilation uses `-O3 -march=native`; `-ffast-math` is intentionally excluded to preserve floating-point/physics semantics.

### Enable xsec stage
```bash
bash src/analysis/run_parallel_nps_analysis_main.sh \
  --config <config_csv> \
  --kin <Kin_old> \
  --source waveform \
  --gevnum-cut no \
  --jobs <N> \
  --types production,Production \
  --run-xsec \
  --xsec-target LH2 \
  --output-base <output_base>
```

## Output Contracts
These paths are contract-sensitive and must not be changed silently.

Per-kin base:
- `output/<kin>/`

Core analysis outputs:
- Per-run diagnostics ROOT: `output/<kin>/root/diagnostics_run<run>.root`
- Per-kin summary CSV: `output/<kin>/summary/summary_all_runs.csv`
- Per-run logs: `output/<kin>/logs/analysis_main_run<run>.log`
- Per-run summary text: `output/<kin>/summary/summary_run<run>.txt`

Combine outputs:
- Combined ROOT: `output/<kin>/root/combined_branches_<target>.root`
- Combine log: `output/<kin>/logs/combine_<kin>.log`

Smearing outputs:
- Section map CSV and interpolated map ROOT under the selected root/output folder
- Smearing log: `output/<kin>/logs/smearing_<kin>.log`

Xsec outputs:
- ROOT/CSV/PDF products under selected output folder
- Xsec log: `output/<kin>/logs/xsec_<kin>.log`

Simulation kinematics config output:
- `config/nps_simulation_kinematics.csv`
- One row per production/LH2 `Kin_old`; malformed comment fields and physics ranges are validated.
- Existing offsets are preserved. SIMC and Geant4 offsets use separate, unit-labelled columns.
- `simc_spec_p_p_gev` defaults explicitly to `Ebeam-HMS_P` and remains user-overridable.
- SIMC event count and target offsets are explicit CSV fields, preventing template drift.

SIMC infile generation output:
- `config/simc_infiles/nps_excl_pi0_<kin_tag>.inp`
- `config/simc_infiles/nps_sidis_pi0_<kin_tag>.inp`
- `config/simc_infiles/nps_delta_pi0_<kin_tag>.inp`
- `kin_tag` is derived from `Kin_old` (for example `KinC_x60_4b -> x60_4b`).

Post-analysis data/SIMC plots:

```bash
python3 scripts/plot_publication_quality.py \
  --combined-root output/KinC_x25_3/root/combined_branches_LH2.root \
  --kin KinC_x25_3
```

- Reads combined data tree `physics` and matching exclusive, SIDIS, and delta
  Geant4 trees from the newest complete production under
  `/volatile/hallc/nps/singhav/geant4_simc/`.
- Reads each channel's `normfac` and `Ngen` from matching SIMC `.hist` file.
- Data weight is `scale*pi0_weight*(charge_uC/total_unique_run_charge_uC)`;
  simulation weight is `normfac*Weight/Ngen`. Absolute plots therefore show
  corrected yield per 1 mC while preserving channel-relative normalization.
- Missing mass is always plotted before any exclusive cut over 0--3 GeV.
- Writes 30 individual PNG/PDF comparisons, one multipage PDF, and one JSON
  normalization/selection audit under `<kin>/plots/post_analysis/`.
- By default also reads `<kin>/summary/summary_all_runs.csv` plus
  `output/efficiency_stuff/efficiency_<kin>.csv` and writes 35 stability plots
  under `<kin>/plots/post_analysis/run_stability/`: run and current trends,
  rate comparisons, normalized yield, pi0 peak/width and fit stability,
  efficiencies, beam quantities, and key quantities versus corrected rate.
- Use `--summary-csv` or `--efficiency-csv` for nonstandard layouts;
  `--ignore-runs` filters stability plots noninteractively. Use
  `--skip-run-stability` only when those inputs are intentionally unavailable.
- Use `--normalization shape` only for per-variable shape checks. Previous
  summary/efficiency run-QA plots remain available with `--legacy-run-qa`.

## How to Interpret Key Outputs
- `diagnostics_run<run>.root`:
  - event-level tree `physics`
  - per-run histograms and diagnostics canvases
- `summary_all_runs.csv`:
  - one line per successfully processed run
  - includes charge/current, selected counts, background metrics, pi0 fit metrics, and run status
- combined ROOT:
  - merged physics branches with run metadata and scale handling

## Known Issues and Debugging

### 1) Run selection mismatch from CSV spacing
Symptom:
- No runs matched selected kinematics/types.
- No Kin_old entry found for run.

Cause:
- Header/field spacing in master CSV can break strict matching in shell/awk paths.

Workaround used in baseline:
- Generate a trimmed config snapshot and run with `--config` pointing to that trimmed file.

Current status:
- SIMC CSV/infile generators normalize whitespace-padded CSV headers/fields.

### 2) AcceptanceCuts unresolved symbols in ROOT/Cling
Symptom:
- unresolved `AcceptanceCuts` symbols in run logs.

Current fix:
- driver preloads `acceptance_cuts.cpp` before executing `nps_analysis_main.C`.

### 3) Combine stage missing efficiency metadata
Symptom:
- combine log reports missing per-kin efficiency CSV.

Resolution:
- ensure `<output_base>/efficiency_stuff/efficiency_<Kin_old>.csv` exists before enabling combine.

## Phase 0 Baseline Reference
Baseline package for x60_4b:
- `output_phase0_baseline/phase0_x60_4b_20260602/README_phase0_baseline.md`

It includes:
- successful representative single-run and 3-run batch commands
- environment snapshot
- checksums and archived references

## Change Control
Any change to:
- output naming,
- output locations,
- tree/branch/schema contracts,
- or physics-sensitive flow
must be documented in `plan.md` before merge.
