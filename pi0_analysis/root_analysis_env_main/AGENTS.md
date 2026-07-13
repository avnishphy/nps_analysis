# Repository Guidelines

## Scope

These instructions apply to the `pi0_analysis/root_analysis_env_main`
workspace. They are intended for Codex/agent sessions and for humans making
changes in this staged Hall C NPS pi0 analysis environment.

## Token And Context Budget

- Default to concise, high-signal replies. Use fuller explanation only when the
  user asks, ambiguity would cause mistakes, or safety/irreversible actions need
  explicit wording.
- Prefer targeted reads over broad dumps. Use `rg`, `rg --files`, `sed -n`, and
  focused diffs; avoid loading generated outputs unless the task needs them.
- Summarize large command outputs instead of pasting them back. Preserve exact
  errors, commands, file paths, run numbers, fit parameters, and
  physics-relevant values.
- For broad investigation or review, prefer compressed/structured agent output
  when available, such as local `.agents/skills/cavecrew` guidance.

## Project Overview

This repository contains the newer staged Hall C NPS neutral-pion analysis
workflow. The main stages are:

- Top-level orchestration with `scripts/run_pipeline.sh`
- Per-run data analysis and branch combination in `src/analysis/`
- HMS efficiency and livetime computation in `src/efficiencies/`
- SIMC/data smearing and simulation-side production in `src/simulation_smearing/`
- Cross-section extraction in `src/xsec_extract/`
- Production outputs under `output/` and baseline captures under
  `output_phase0_baseline/`

Efficiency computation is frozen for the current phase. Preserve validated
physics behavior unless the task explicitly asks for a controlled change.
Geant4 internals are out of scope for this phase.

## Important Files

- `README.md`: Main workflow, command examples, output contracts, and cautions.
- `config/nps_dvcs_all_kins_main.csv`: Master run and kinematic metadata table.
- `config/acceptance_cuts.conf`: Acceptance and timing cut definitions.
- `config/nps_simulation_kinematics.csv`: Canonical simulation-side kinematic
  table.
- `config/dead_block_per_runs.csv`: Per-run NPS dead-block metadata.
- `src/analysis/nps_analysis_main.C`: Main per-run analysis macro.
- `src/analysis/run_parallel_nps_analysis_main.sh`: Main parallel data-stage
  driver.
- `src/analysis/combine_analysis_branches.py`: Per-kin branch combination.
- `src/efficiencies/compute_efficiencies_stuff.cxx`: HMS efficiency/livetime
  executable source.
- `src/efficiencies/run_efficiencies_stuff_parallel.sh`: Parallel efficiency
  launcher.
- `src/simulation_smearing/nps_sim_smearing_new.C`: Smearing fit/map generation.
- `src/simulation_smearing/simc_pi0_analysis.C`: Simulation producer with
  smearing-map usage.
- `src/simulation_smearing/run_smearing_pipeline.sh`: Smearing pipeline driver.
- `src/xsec_extract/run_xsec_pipeline.sh`: Cross-section stage driver.

## Working Tree Rules

- Do not revert, delete, or overwrite user changes unless explicitly asked.
- This repository often contains generated ROOT dictionaries, shared objects,
  logs, plots, and CSV outputs. Treat them as artifacts unless the task
  explicitly asks to update or inspect them.
- Prefer targeted status/diff checks over broad repository scans.
- Avoid destructive commands such as `git reset --hard`, `git checkout --`, or
  broad cleanup unless the user explicitly requests them.
- When editing manually, use `apply_patch`.

## Search And Navigation

- Prefer `rg` for text search and `rg --files` for file discovery.
- Exclude generated files when possible: ROOT files, PDFs, PNGs, compiled
  binaries, dictionary files, and notebook checkpoints can make searches noisy.
- Useful targeted searches:

```bash
rg "CSV_WRITTEN|summary_all_runs|run_status" src scripts
rg "AcceptanceCuts|acceptance_cuts|NPS_ACCEPTANCE_CUTS_CONFIG" src scripts
rg "dead_block|section_map|NPS_SMEAR_RANDOM_SEED|NPS_SMEARING_MODE" src scripts config
rg "efficiency_stuff|combine-target|run-smearing|run-xsec" src scripts README.md
```

## Running Analysis

The scripts assume ROOT is available in `PATH`. Load the Hall C/NPS ROOT
environment before running ROOT, `root-config`, compiling ROOT macros, or
running the smearing pipeline.

For interactive `csh`/`tcsh`, use:

```csh
source /group/nps/singhav/setup.csh
```

Most Codex shells are `bash`, so run ROOT commands through `csh` after sourcing:

```bash
csh -c 'source /group/nps/singhav/setup.csh; root-config --version'
csh -c 'source /group/nps/singhav/setup.csh; ./scripts/run_pipeline.sh --help'
```

In non-interactive `csh`/`tcsh` sessions, source
`/usr/share/Modules/init/csh` first if `module` is not defined.

Representative pipeline command:

```bash
./scripts/run_pipeline.sh \
  --kin KinC_x60_4b \
  --source waveform \
  --jobs 3 \
  --gevnum-cut no \
  --types production,Production \
  --no-combine
```

Representative single-run analysis command:

```bash
bash src/analysis/run_parallel_nps_analysis_main.sh \
  --kin KinC_x60_4b \
  --run 4253 \
  --source waveform \
  --gevnum-cut no \
  --jobs 1 \
  --types production,Production \
  --no-combine
```

Do not run full production jobs during orientation or documentation-only work.
They can be long and may write many artifacts.

## Output Contracts

Output paths and naming contracts are explicit and must not change silently.
Document any change to output naming, locations, tree/branch/schema contracts,
or physics-sensitive flow in `plan.md` before merge.

Key contracts include:

- Per-kin base: `output/<kin>/`
- Per-run diagnostics ROOT:
  `output/<kin>/root/diagnostics_run<run>.root`
- Per-kin summary CSV: `output/<kin>/summary/summary_all_runs.csv`
- Per-run logs: `output/<kin>/logs/analysis_main_run<run>.log`
- Per-run summary text: `output/<kin>/summary/summary_run<run>.txt`
- Combined ROOT: `output/<kin>/root/combined_branches_<target>.root`
- Efficiency CSVs under `output/efficiency_stuff/`
- Simulation kinematics config: `config/nps_simulation_kinematics.csv`
- Generated SIMC infiles under `config/simc_infiles/`

## Coding Conventions

- Keep ROOT macro and shell driver interfaces stable unless callers are updated
  together.
- Prefer explicit input/output paths over hidden global state, while preserving
  existing environment-variable overrides.
- Keep comments focused on physics/analysis intent, not line-by-line mechanics.
- Use ASCII in new docs/code unless an existing file already uses symbols and
  clarity benefits from them.
- For Python scripts, keep path constants near the top and document whether a
  script is production-ready or diagnostic.

## Validation

Pick validation proportional to the change:

- Documentation-only changes: no ROOT run required.
- Shell wrapper changes: run `--help` or a tiny/single-run command when feasible.
- Macro changes: compile or run ROOT in batch on one known run when feasible.
- Efficiency changes: preserve frozen behavior; run the smallest representative
  kinematic/run selection and compare CSV columns/values if changes are
  unavoidable.
- Python plotting/combining changes: run the script on the smallest relevant
  input or at least syntax-check/import where dependencies permit.

## Known Cautions

- Efficiency computation is frozen for this phase; avoid behavioral edits unless
  explicitly requested.
- Combine requires per-kin efficiency CSV availability under
  `<output_base>/efficiency_stuff`.
- The analysis driver preloads `acceptance_cuts.cpp` to avoid unresolved
  `AcceptanceCuts` symbols in ROOT/Cling.
- Header or field spacing in CSV inputs can affect strict shell/awk matching;
  several generators normalize whitespace, but older paths may still be strict.
- Many scripts contain absolute Hall C filesystem paths. Check path assumptions
  before running on another machine or scratch area.
- Notebook files may have user work embedded. Do not rewrite notebooks unless
  directly requested.
