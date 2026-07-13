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
  - `generate_simulation_kinematics_csv.py`: canonical simulation CSV generator
  - `generate_simc_infiles.py`: per-kin SIMC infile generator from simulation CSV
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
- Initializes all simulation offset columns to `0.0`.
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
  --sim-chain-geant4-cmd "run_geant4 --kin {kin_old} --outdir {geant4_output_dir}" \
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
- One row per `Kin_old`, with offsets explicitly present and initialized to zero in the current phase.

SIMC infile generation output:
- `config/simc_infiles/nps_excl_pi0_<kin_tag>.inp`
- `kin_tag` is derived from `Kin_old` (for example `KinC_x60_4b -> x60_4b`).

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
