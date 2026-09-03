# NPS Pi0 Analysis Refactor Plan (Living Document)

Last updated: 2026-06-02  
Status: Phase 0 and Phase 1 completed; Phase 2 complete; Phases 3, 4, 5, 6, and 7 completed for current scope with regression/interface checks

## 1. Purpose

This document defines the staged refactor and workflow hardening plan for the NPS neutral-pion analysis framework located at:

`/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main`

The goal is to preserve validated physics behavior while improving:
- code organization,
- histogram ownership clarity,
- reproducibility,
- simulation automation,
- output traceability,
- diagnostics for run-level and kinematic-level debugging,
- and one-click execution of the full analysis workflow.

This is a production analysis directory. Any change that can alter physics outputs, file layout, or workflow behavior must be explicitly documented in this file before it becomes part of the default pipeline.

---

## 2. Scope

### 2.1 In scope
- ROOT-based analysis code in this directory.
- Header/source separation and histogram lifecycle structure.
- Pipeline orchestration scripts inside this repository.
- Diagnostic plot generation and output organization.
- SIMC input-file generation from analysis-side configuration.
- Controlled automation hooks to the existing SIMC and Geant4 simulation chain.
- README and plan documentation.

### 2.2 Out of scope for this phase
- Any redesign of Geant4 internals.
- Any redesign of the SIMC physics model.
- Any change to efficiency computation.
- Any unvalidated change to physics selection thresholds or final cross-section logic.
- Any silent renaming of outputs or file paths.

---

## 3. Non-negotiable constraints

1. **Efficiency computation is frozen.**  
   Do not modify the efficiency calculation, its inputs, its algebra, or its produced values in this phase.

2. **Validated physics behavior must be preserved.**  
   If a change can affect event yields, selection, weighting, smearing, or extracted cross sections, it must be treated as physics-sensitive and compared against a baseline.

3. **Geant4 internals are not to be changed here.**  
   This repository may invoke Geant4-based workflows, but implementation changes inside Geant4 belong to a separate task.

4. **ROOT object ownership must be explicit.**  
   Histogram creation, filling, and writing must have clear ownership and must not be scattered across headers and source files.

5. **No silent contract changes.**  
   Output file names, output directories, CSV schemas, tree names, branch names, and diagnostic artifacts must not change without updating this plan and the README.

6. **Every automation entry point must be fail-fast.**  
   If required inputs, files, or environment variables are missing, scripts must exit with a clear error message.

---

## 4. Definitions

### 4.1 Main analysis file
The primary C++/ROOT analysis implementation file that owns the event loop and the analysis lifecycle. It must be the canonical place where the per-run analysis flow is readable and debuggable.

### 4.2 Header file
A header may contain:
- declarations,
- small inline helpers,
- constants,
- configuration structs,
- function prototypes.

A header must not create or own analysis histograms unless there is a very strong, documented reason and the design remains clean. The default rule is: **no histogram booking in headers**.

### 4.3 Histogram lifecycle
For each histogram, the code must make clear:
- where it is declared,
- where it is booked,
- where it is filled,
- where it is written,
- where it is deleted or allowed to go out of scope.

### 4.4 Baseline
A baseline is a frozen reference set of outputs produced before the refactor changes are merged. It includes logs, output ROOT files, summary CSVs, and representative diagnostic plots.

### 4.5 Physics regression
A physics regression is any unintended change in:
- event counts,
- selection acceptance,
- bin contents,
- weighting behavior,
- smearing behavior,
- summary yields,
- cross-section extraction,
- or any final observable that should remain stable.

### 4.6 Diagnostic plot
A diagnostic plot is any plot whose purpose is to help identify a run-quality, selection, reconstruction, weighting, smearing, or simulation mismatch problem. Every diagnostic plot must have a reason to exist.

---

## 5. Required pre-edit inventory

Before any code is changed, inspect the full tree and produce an internal map of:

1. All source files, headers, scripts, and configuration files.
2. The exact analysis entry points.
3. The exact simulation entry points.
4. The exact cross-section extraction entry points.
5. The exact locations where histograms are created, filled, and written.
6. All output directories and filename conventions.
7. All external dependencies and assumptions.
8. All places where code is duplicated or behavior is spread across headers and source files.
9. Any hidden coupling between scripts, configs, and source code.
10. Any existing diagnostic plots and how they are generated.

No refactor should start until this inventory is complete.

---

## 6. Current high-level workflow to preserve

The intended workflow is:

1. Read run/kinematic configuration.
2. Process data through the main analysis.
3. Produce run-level diagnostics.
4. Produce summary outputs.
5. Combine branches or intermediate products if required.
6. Apply simulation smearing in a controlled and physics-driven way.
7. Produce simulation inputs and outputs needed for comparison.
8. Run xsec extraction.
9. Save diagnostics, summaries, and final observables in well-defined output locations.

This workflow must remain transparent and modular.

---

## 7. Files in scope

The following files are expected to be examined and updated as needed:

- `src/analysis/nps_analysis_main.C`
- `src/analysis/run_parallel_nps_analysis_main.sh`
- `src/analysis/combine_analysis_branches.py`
- `src/analysis/acceptance_cuts.h`
- `src/analysis/acceptance_cuts.cpp`
- `src/simulation_smearing/run_smearing_pipeline.sh`
- `src/simulation_smearing/simc_pi0_analysis.C`
- `src/simulation_smearing/nps_sim_smearing_new.C`
- `src/xsec_extract/run_xsec_pipeline.sh`
- `config/acceptance_cuts.conf`
- `config/nps_dvcs_all_kins_main.csv`
- `config/nps_simulation_kinematics.csv` (new or updated canonical simulation table)
- `scripts/` (new orchestration and helper scripts)

If additional files are discovered during inventory, they must be added here before refactor work proceeds.

---

## 8. Output contract

These outputs are considered current contracts unless this plan explicitly changes them:

- Per-run diagnostics ROOT files:  
  `output/<kin>/root/diagnostics_run<run>.root`

- Combined ROOT output:  
  `output/<kin>/root/combined_branches_<target>.root`

- Summary CSV aggregate:  
  `output/<kin>/summary/summary_all_runs.csv`

- Smearing calibration artifacts:  
  `section_map.csv`, interpolated smear map ROOT file

- Smeared SIM output:  
  `simc_pi0_analysis_output_smeared.root`

Any future change to output names or output locations must be written into this plan before implementation and mirrored in the README.

---

## 9. Execution phases

## Phase 0 — Baseline freeze

### Goal
Capture reference outputs from the current workflow before refactoring.

### Required actions
1. Run the current analysis on:
   - one representative single run,
   - one representative multi-run kinematic set.
2. Save the exact command lines used.
3. Save the environment information:
   - git commit hash,
   - relevant environment variables,
   - shell version,
   - ROOT version if applicable,
   - Python version if applicable.
4. Archive the resulting outputs:
   - logs,
   - per-run diagnostic ROOT files,
   - summary CSVs,
   - combined ROOT outputs,
   - smearing outputs,
   - xsec outputs,
   - representative diagnostic PDFs/PNGs if produced.

### Exit criteria
- Baseline artifacts exist and are clearly labeled.
- The command used to generate them is documented.
- The outputs are sufficient for later comparison.

---

## Phase 1 — Documentation and contract locking

### Goal
Document the workflow and freeze the contracts before code changes.

### Required actions
1. Update `README.md` with:
   - overview of the physics goal,
   - directory layout,
   - dependencies,
   - end-to-end workflow,
   - how to run each stage,
   - how to interpret outputs,
   - how to debug common failures.
2. Update this `plan.md` as the living implementation log.
3. Lock the contracts for:
   - output file names,
   - output directories,
   - tree names,
   - branch names,
   - summary CSV schema,
   - histogram names,
   - diagnostic plot naming,
   - simulation table schema.

### Exit criteria
- README and plan both describe the same workflow.
- All contract-sensitive file names and schemas are explicitly documented.

---

## Phase 2 — Script-layer orchestration

### Goal
Create a single, transparent orchestration layer for the full workflow.

### Required actions
1. Add a top-level orchestrator, expected to live in `scripts/run_pipeline.sh`.
2. Add shared script helpers for:
   - path resolution,
   - argument validation,
   - logging,
   - stage selection,
   - environment checks,
   - fail-fast error handling.
3. Add or update scripts for:
   - SIMC input generation,
   - simulation-chain invocation,
   - simulation-kinematics table generation,
   - diagnostic PDF creation,
   - branch-combination support,
   - xsec pipeline triggering.

### Behavior requirements
- Scripts must be idempotent where possible.
- Scripts must echo the exact commands they run.
- Scripts must write logs to stage-specific log files.
- Scripts must not hide missing-file or missing-variable errors.
- Scripts must support stage toggles so users can run the full pipeline or only selected stages.
- Scripts must not hardcode per-kinematic physics constants if those values are meant to come from configuration.

### Required simulation-kinematics artifact
Create or update:

`config/nps_simulation_kinematics.csv`

This file is the canonical per-kinematic simulation input table for the analysis directory.

#### Purpose
It provides a single analysis-side source of truth for the per-kinematic values needed by SIMC and Geant4 wrappers.

#### Row structure
- One row per analysis kinematic setting.
- The row must be mapped to the corresponding `Kin_old` / `New_kin` identity used by the analysis.

#### Required physics fields
The table must contain the base kinematic values needed by the simulation chain, including:
- beam energy,
- HMS momentum,
- HMS angle,
- NPS angle,
- NPS target distance.

These values must be copied from `config/nps_dvcs_all_kins_main.csv` using the same selection conventions used by the analysis code. No alternate selection logic may be introduced silently.

#### Required offset fields
The table must contain explicit offset columns for simulation inputs. The first generated version must initialize all offset values to zero.

The default rule is:
- offsets exist explicitly in the CSV,
- offsets are not hardcoded in wrapper scripts,
- wrapper scripts read the offset values from this CSV.

#### Provenance fields
The CSV must also include provenance metadata:
- source run,
- source SIMC infile,
- selection rule,
- generation timestamp in UTC.

#### Required contract for scripts
Scripts that write SIMC or Geant4 inputs must read the CSV and must not invent or override per-kinematic offsets internally.

### Exit criteria
- One top-level pipeline script exists.
- The script layer can run the workflow in a controlled way.
- The simulation-kinematics CSV exists, is documented, and is used by the automation layer.

---

## Phase 3 — Core analysis modularization

### Goal
Make the main analysis code readable, physics-transparent, and easy to debug without changing validated behavior.

### Required actions
1. Refactor the main analysis implementation into explicit lifecycle sections:
   - initialization,
   - run loading,
   - event selection,
   - histogram booking,
   - histogram filling,
   - background subtraction,
   - pi0 weighting,
   - second-pass or weighted filling,
   - output writeout,
   - cleanup.
2. Keep the analysis logic in coherent implementation files.
3. Move histogram booking and ownership out of headers unless a specific inline helper truly requires it.
4. Keep helper code in headers limited to declarations, simple inline utilities, constants, and small pure functions.
5. Preserve existing output object names, binning, and data flow unless a change is explicitly recorded in this plan.

### Histogram rules
- Every histogram must have a single clear owner.
- Every histogram must be created in one place.
- Every histogram must be filled from one clear analysis path.
- Histograms must be written once at the appropriate output stage.
- No analysis histogram may be silently created in a header if that obscures ownership.

### Exit criteria
- The main analysis code is logically segmented.
- Histogram lifecycle is obvious.
- Output content remains compatible with baseline, unless an explicit, documented change is approved.

---

## Phase 4 — Linkage and build correctness

### Goal
Remove build or runtime issues that obscure the analysis flow.

### Required actions
1. Fix the `AcceptanceCuts` unresolved symbol issue.
2. Ensure the analysis builds and runs cleanly in the expected environment.
3. Make sure any linkage workaround is documented and local to the relevant compilation path.
4. Verify that the fix does not change physics behavior.

### Exit criteria
- No unresolved `AcceptanceCuts` symbol errors appear.
- The analysis builds and executes cleanly.
- The output physics content remains consistent with the baseline.

---

## Phase 5 — Simulation smearing cleanup

### Goal
Make the smearing stage physics-driven, non-redundant, and easy to inspect.

### Required actions
1. Identify the current smearing flow and all quantities it modifies.
2. Remove duplicated transformations or repeated corrections.
3. Make every smeared variable traceable to:
   - a source variable,
   - a smearing model,
   - a unit convention,
   - a parameter source.
4. Preserve the scientific intent of the current smearing.
5. Do not add ad hoc corrections unless they are explicitly physics-motivated and documented.
6. Keep random-number usage reproducible by controlling the seed or seed strategy.
7. Preserve the existing smearing outputs unless a change is explicitly documented.

### Optimizer performance safeguards
- Parameter-dependent photon response terms are prepared outside the `Nsmear` loop; only stochastic pulls vary inside it.
- Immutable event quantities (`E_safe`, `log(E_safe)`, directions, static normalization-bin membership, and weight/`Nsmear`) are cached once per fit context.
- When position smearing is inactive, cached directions remove repeated geometry square roots. One opening-angle calculation supplies both `M_gg` and `(p_target+gamma gamma)^2`.
- Gaussian pulls retain the deterministic seed and call ordering and may be cached as doubles under a bounded memory budget.
- Optimizer histograms use equivalent weighted bin/sumw2 accumulation; final diagnostics and persisted ROOT histograms retain the ROOT path.
- Sobol candidates are processed in fixed-size batches without changing candidate ordering, pulls, histogram fill order, or the objective; uncached/Landau modes use the scalar evaluator.
- Global-prefit Sobol batches and retained MIGRAD starts may run in parallel with independent Minuit2 instances. Section fits use their existing outer OpenMP parallelism and keep inner work serial.
- Scalar-fast, batched-fast, and legacy objectives are compared at three parameter points before each fit, with automatic legacy fallback on disagreement.
- Event selection, weights, normalization windows, chi2 definition, Sobol/MIGRAD settings, coupled sweeps, and final diagnostics remain unchanged.

### Physics requirements
- Smearing must not be applied twice to the same quantity.
- Smearing must be documented in terms of what physical effect it is approximating.
- The model should avoid hidden normalization or unit inconsistencies.
- Any interpolation or section-map generation must be reproducible and attributable to the correct input sources.

### Exit criteria
- Smearing logic is simpler, traceable, and reproducible.
- There are no redundant smear operations.
- Produced smeared outputs remain consistent with the intended physics model.

---

## Phase 6 — SIMC and Geant4 workflow integration

### Goal
Allow the analysis directory to drive the simulation chain in an automated but controlled way.

### Required actions
1. Create analysis-side wrapper scripts that:
   - read `config/nps_simulation_kinematics.csv`,
   - generate per-kin SIMC input files,
   - pass the appropriate values to the SIMC workflow,
   - hand off to the Geant4 simulation stage where required,
   - collect outputs into the analysis workflow.
2. Keep this layer separate from the internal Geant4 implementation.
3. Keep the wrappers explicit about:
   - input file location,
   - output file location,
   - kinematic setting,
   - offset values,
   - provenance of generated files.

### Required wrapper behavior
- No hidden defaults for required simulation parameters.
- No silent fallback to hardcoded constants.
- No change to Geant4 internals in this phase.
- No use of undocumented input file templates without explicit mapping.

### Exit criteria
- The analysis directory can launch the required simulation workflow.
- Inputs are generated from the canonical CSV.
- Simulation outputs are retrieved into a predictable location.

---

## Phase 7 — Diagnostics and output hygiene

### Goal
Make debugging straightforward for bad runs, poor kinematics, or simulation mismatches.

### Required actions
1. Preserve existing diagnostic plots and add missing ones where needed.
2. Organize outputs by:
   - run,
   - kinematic setting,
   - analysis stage,
   - production vs diagnostic purpose.
3. Ensure that summary outputs and diagnostic outputs are easy to distinguish.
4. Make sure every diagnostic plot has a clear failure mode it is meant to reveal.

### Minimum diagnostic categories
The analysis should be able to diagnose:
- run quality,
- event counts,
- current and livetime behavior,
- vertex distributions,
- energy and angle residuals,
- missing mass or invariant-mass stability,
- background subtraction behavior,
- kinematic agreement,
- smearing effects,
- simulation-to-data discrepancies.

### Exit criteria
- Diagnostic outputs are present and well organized.
- A bad run can be investigated without editing code.
- Diagnostic naming and placement are stable.

---

## 10. Validation matrix

Every major phase change must be checked against the following matrix.

### 10.1 Single-run smoke test
- Analysis completes successfully.
- Expected per-run diagnostic ROOT output appears.
- Expected summary row appears.
- No fatal warnings or crashes.

### 10.2 Multi-run parallel test
- Parallel execution does not cause temporary-file collisions.
- Summary regeneration is deterministic.
- Combined outputs are correct and complete.

### 10.3 Combine-stage validation
- Expected branches and trees are present.
- Output naming matches the documented contract.
- Debug-mode overrides are explicit and local.

### 10.4 Smearing-stage validation
- Section map is generated.
- Interpolated smearing map is generated.
- Smeared SIM output is generated.
- Smearing is reproducible for the same seed and input.

### 10.5 Xsec-stage validation
- Expected ROOT, CSV, and PDF outputs are generated.
- Output content is consistent with the pipeline inputs.
- No missing input dependency remains hidden.

### 10.6 Linkage/build validation
- No unresolved symbol errors.
- Code builds and runs in the documented environment.

### 10.7 Physics regression validation
Compare outputs to the baseline for:
- key spectra,
- yield trends,
- weighted distributions,
- cross-section outputs,
- diagnostic shapes.

If any difference appears, determine whether it is:
- expected and documented,
- numerical noise within tolerance,
- or an actual regression.

---

## 11. Change-control rules

1. Every meaningful code change must be traceable to a phase in this plan.
2. Every phase completion must update this document.
3. Any change in file path, output name, histogram name, CSV schema, or physics-flow ordering must be recorded here before merge.
4. If a validation fails, do not continue to later phases until the failure is explained.
5. If a change must be deferred, document the deferral explicitly in the phase status.

---

## 12. Status tracker

- [x] Phase 0 baseline freeze complete (x60_4b representative single-run 4253 and representative 3-run batch 4253/4254/4300 captured with `--gevnum-cut no`; see `output_phase0_baseline/phase0_x60_4b_20260602/README_phase0_baseline.md`).
- [x] Phase 0 x60_4b baseline snapshot package created at `output_phase0_baseline/phase0_x60_4b_20260602/`.
- [x] README authored or updated (`README.md` now documents run workflow, outputs, and debug caveats).
- [x] Contract definitions locked (output paths/names and stage contracts documented in `plan.md` and `README.md`).
- [x] Pipeline orchestration scripts added/updated (`scripts/run_pipeline.sh` and consolidated `scripts/generate_simc_infiles.py`).
- [x] `config/nps_simulation_kinematics.csv` created and validated (offset columns initialized to `0.0`, includes `KinC_x60_4b`).
- [x] Canonical CSV and SIMC infile generation consolidated in `scripts/generate_simc_infiles.py` and wired into `scripts/run_pipeline.sh`.
- [x] Consolidated generator separates review stages: default writes CSV; `--gen_infile` consumes reviewed CSV.
- [x] SIMC infile generation validated for representative kinematics (`KinC_x60_4b`, `KinC_x36_1`) with CSV-driven `Ebeam`, `spec%e%P`, `spec%e%theta`, and explicit offset fields.
- [x] Phase 2 end-to-end smoke run passed via `scripts/run_pipeline.sh` (`KinC_x60_4b`, run `4253`, `--gevnum-cut no`, `--no-combine`, output base `output_phase2_smoke`).
- [x] Pipeline SIMC generation now mirrors forwarded `--kin` selection by default (and supports explicit `--simc-kin` override), avoiding unnecessary all-kin infile emission.
- [x] Phase 3 started: extracted reusable helpers in `nps_analysis_main.C` for run kinematics resolution, cluster loading, good-cluster collection, and pair selection (no physics algorithm changes).
- [x] First and second event passes now share the same helper path for good-cluster collection and π0 pair selection to reduce duplication risk.
- [x] Post-refactor smoke validation passed (`KinC_x60_4b`, run `4253`, `--gevnum-cut no`, `--no-combine`, output base `output_phase3_smoke`).
- [x] Overlay/writeout modularization pass: repetitive 1D overlay canvas construction consolidated into a reusable helper while preserving output file names and downstream canvas write/cleanup behavior.
- [x] Post-overlay-refactor smoke validation passed (`KinC_x60_4b`, run `4253`, `--gevnum-cut no`, `--no-combine`, output base `output_phase3_smoke`).
- [x] Summary serialization modularization pass: per-run CSV row and TXT summary formatting moved to reusable helpers, preserving schema/content and output paths.
- [x] Post-summary-refactor smoke validation passed (`KinC_x60_4b`, run `4253`, `--gevnum-cut no`, `--no-combine`, output base `output_phase3_smoke`).
- [x] Event-level physics tree modularization pass: per-run `physics` tree branch declaration/fill/write extracted to a reusable helper, preserving branch names and semantics.
- [x] Post-tree-refactor smoke validation passed (`KinC_x60_4b`, run `4253`, `--gevnum-cut no`, `--no-combine`, output base `output_phase3_smoke`); ROOT output still contains the `physics` tree.
- [x] Main analysis modularization complete for current scope (run-level helper extraction completed for kinematics resolution, cluster/pair logic, overlay generation, summary serialization, and physics tree writeout; no output-contract changes).
- [x] AcceptanceCuts linkage issue fixed for driver execution path (`run_parallel_nps_analysis_main.sh` preloads `acceptance_cuts.cpp` before executing `nps_analysis_main.C`).
- [x] Phase 4 linkage/build correctness validated in runtime smoke (`KinC_x60_4b`, run `4253`): no unresolved symbol errors and successful end-to-end analysis execution.
- [x] Smearing cleanup complete for current scope: deterministic smearing-seed control added (`NPS_SMEAR_RANDOM_SEED`, `--smear-seed`) and wired through `run_smearing_pipeline.sh` to `simc_pi0_analysis.C`.
- [x] SIMC/Geant4 wrappers integrated for current scope: added `scripts/run_simulation_chain.py` and `scripts/run_pipeline.sh --run-sim-chain` stage to consume canonical simulation CSV + generated SIMC infiles and launch explicit user-provided SIMC/Geant4 command templates.
- [x] SWIF2 SIMC submission added via `scripts/submit_simc_swif2.sh` (one isolated job per infile; persists normal `worksim`, `runout`, and `outfiles` products).
- [x] Phase 6 interface validation completed (dry-run): wrapper help/options validated; `KinC_x60_4b` dry-run generated command traces and provenance manifest (`output/simulation_chain_manifest_dryrun.csv`) with predictable simulation output paths under `<output-base>/<kin_tag>/simulation/`.
- [x] Whole-production-kinematics simulation inputs hardened: malformed master-CSV comment rows repaired deterministically; LH2/good production rows and physical ranges validated; calibrated offsets preserved; SIMC NPS angle/distance/pion scale and separate SIMC/Geant4 offsets propagated explicitly.
- [x] Per-kin SIMC generation emits exclusive, SIDIS, and delta channels with explicit `doing_semi`/`which_pion` channel settings.
- [x] Diagnostics expanded and organized for current scope: added searchable diagnostics index (`scripts/generate_diagnostics_index.py`) and modular detector/experiment reports (`scripts/generate_modular_diagnostics_reports.py`) with pipeline toggles for HMS (`--build-hms-diagnostics`), NPS (`--build-nps-diagnostics`), and whole-experiment (`--build-experiment-diagnostics`) outputs.
- [x] HMS, NPS, and whole-experiment diagnostics validated on smoke outputs (`output_phase3_smoke`): modular reports confirm detector-side coverage (`has_hms_diagnostics=yes`, `has_nps_diagnostics=yes`) and experiment-level readiness for `KinC_x60_4b`.
- [x] Post-analysis publication plotting consumes combined data plus raw exclusive/SIDIS/delta Geant4 trees; absolute normalization uses charge-weighted combined-data `scale*pi0_weight` and `.hist`-derived `normfac*Weight/Ngen`, with pre-exclusive missing mass and JSON provenance. Default execution also retains all run/current/rate stability, normalized-yield, fit, efficiency, and beam-trend plots under `post_analysis/run_stability/`.
- [x] Post-analysis plot contract: 30 individual PNG/PDF comparisons, `post_analysis_<kin>_<normalization>.pdf`, and `post_analysis_<kin>_<normalization>_metadata.json` under `<kin>/plots/post_analysis/`; synthetic ROOT smoke passed for all branches and plots.
- [x] Redundant temporary update artifacts cleaned: removed non-canonical dry-run directory `output_phase3_smoke/x60_4b/` and transient script cache directory `scripts/__pycache__/`.
- [ ] Validation matrix passed
- [x] Baseline comparison references recorded (existing snapshot checksums + new phase-0 run checksums in `output_phase0_baseline/phase0_x60_4b_20260602/repro_attempt/meta/`).

---

## 13. Handoff rule

This file is a living implementation log, not a static design note.

Update it whenever:
- a new file is added,
- an output contract changes,
- a phase is completed,
- a validation fails,
- a physics-sensitive decision is made,
- or a previously unknown dependency is discovered.
