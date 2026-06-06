# Repository Guidelines

## Scope

These instructions apply to the `pi0_analysis/root_analysis_env` workspace.
They are intended for Codex/agent sessions and for humans making changes in
this analysis environment.

## Project Overview

This workspace contains ROOT/C++ and Python tooling for Hall C NPS pi0
analysis. The main workflows are:

- Per-run data analysis with ROOT macros in `src/`
- Waveform-calibrated production analysis with `src/nps_analysis_wfpi0.C`
- Batch orchestration with the top-level `run_*.sh` scripts
- SIMC/data smearing and comparison tools in `scripts/`
- Cross-section extraction and publication-style plotting in `scripts/`

Most generated outputs live under `output/`, `output_pi0_xsec/`, and
`output_pi0_xsec_no_simc_model/`. Treat those as analysis artifacts unless a
task explicitly asks to update or inspect them.

## Important Files

- `README.md`: Human-facing overview and quick-start notes.
- `README_wfpi0_analysis.md`: Detailed notes for waveform-calibrated
  production input.
- `config/cuts.conf`: Simple key/value cut configuration.
- `config/branches_to_read.txt`: Branch allow-list for skimming/analysis.
- `config/runlist_*.txt`: Run lists used by orchestration scripts.
- `src/nps_analysis.C`: Main skim-data analysis macro.
- `src/nps_analysis_wfpi0.C`: Production waveform-calibrated variant.
- `scripts/simc_pi0_analysis.C`: SIMC pi0 analysis/smearing consumer.
- `scripts/nps_sim_smearing_new.C`: Section-wise smearing fit executable
  source.
- `scripts/excl_xsec_pi0_analysis*.C`: Cross-section extraction macros.

## Working Tree Rules

- Do not revert, delete, or overwrite user changes unless explicitly asked.
- This repository often contains many generated or modified files. Prefer
  targeted status/diff checks over broad `git status` when possible.
- Avoid destructive commands such as `git reset --hard`, `git checkout --`,
  or broad `rm` cleanup unless the user explicitly requests them.
- When editing manually, use `apply_patch`.

## Search And Navigation

- Prefer `rg` for text search and `rg --files` for file discovery.
- Exclude generated files when possible: ROOT files, PDFs, PNGs, compiled
  binaries, and notebook checkpoints can make searches noisy.
- Useful targeted searches:

```bash
rg "CSV_WRITTEN|NPS_SKIM_DIR|NPS_OUTPUT_DIR|NPS_RUNLIST" src scripts
rg "SetBranchAddress|NPS\\.prod|NPS\\.cal" src scripts
rg "Config::|ENABLE_POSITION_SMEARING|section_map" scripts
```

## Running Analysis

The scripts assume ROOT is available in `PATH`.

Sequential original analysis:

```bash
./run_nps_analysis.sh [runlist] [skim_dir] [output_dir] [beam_energy]
```

Sequential waveform-calibrated analysis:

```bash
./run_nps_analysis_wfpi0.sh [runlist] [production_dir] [output_dir] [beam_energy] [timeout_sec]
```

Parallel waveform-calibrated analysis:

```bash
./run_parallel_nps_analysis_wfpi0.sh [runlist] [production_dir] [output_dir] [beam_energy] [timeout_sec] [nproc]
```

Smearing pipeline:

```bash
./run_smearing_pipeline.sh
```

Do not run full production jobs during orientation or documentation-only work.
They can be long and may write many artifacts.

## Environment Variables

The ROOT analysis macros commonly read:

- `NPS_SKIM_DIR`: input ROOT file directory
- `NPS_OUTPUT_DIR`: output plot/diagnostic directory
- `NPS_RUNLIST`: run list file
- `NPS_EBEAM`: beam energy in GeV
- `NPS_SMEARING_MODE`: smearing mode used by SIMC analysis, for example
  `section`

## Output Conventions

Per-run analysis normally writes:

- `diagnostics_run<run>.root`
- diagnostic PNG/PDF files
- `summary_run<run>.txt`
- per-run CSV rows signaled in console output by `[CSV_WRITTEN]`
- `summary_all_runs.csv` reconstructed by wrapper scripts

The wrapper scripts rely on `[CSV_WRITTEN]` in ROOT output. If macro logging is
changed, preserve that marker or update the wrappers too.

## Coding Conventions

- Keep ROOT macro interfaces stable unless the calling shell scripts are also
  updated.
- Prefer explicit input/output paths over hidden global state, but preserve the
  existing environment-variable overrides.
- Keep comments focused on physics/analysis intent, not line-by-line mechanics.
- Use ASCII in new docs/code unless an existing file already uses symbols such
  as `pi0`/Greek notation and clarity benefits from them.
- For Python scripts, keep path constants near the top and document whether a
  script is production-ready or diagnostic.

## Validation

Pick validation proportional to the change:

- Documentation-only changes: no ROOT run required.
- Shell wrapper changes: run with a tiny run list if available.
- Macro changes: compile or run ROOT in batch on a single known run when
  feasible.
- Python plotting/combining changes: run the script on the smallest relevant
  input or at least import/lint syntax where dependencies permit.

## Known Cautions

- `run_parallel_nps_analysis.sh` has waveform-oriented comments/names but calls
  `src/nps_analysis.C`; verify intended behavior before changing it.
- Many scripts contain absolute Hall C filesystem paths. Check path assumptions
  before running on a different machine or scratch area.
- Notebook files may have user work embedded. Do not rewrite notebooks unless
  directly requested.
