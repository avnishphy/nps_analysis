# NPS Pi0 Analysis Refactor Plan (Living Document)

## 2026-10-01: Explicit canonical raw-export authorization

- User-authorized production destination for the regenerated `KinC_x36_4`
  production-LH2 diagnostics is `output/KinC_x36_4/`.
- Preserve the raw-export safety boundary by requiring both an explicit
  `--output-base` and `--allow-canonical-raw-observation-export` when that base
  is the repository's canonical `output/` directory. Canonical raw export
  remains rejected by default.
- Regeneration is run-only: no combine, efficiency, smearing, or cross-section
  stage is authorized by this exception.

## 2026-10-01: ALG-002B approved LH2-only shadow implementation

- Approved physics increment: a shadow simultaneous timing-and-mass
  extended-Poisson model for `KinC_x36_4`.
- It splits the central true-coincidence timing component into pi0 signal
  and true-coincidence combinatorial background, retain run-specific peak
  shifts, widths, and yields through partial pooling. The complete cross-run
  replica covariance remains a required validation product.
- It does not write `pi0_weight`, alter production outputs, combine physics
  spectra, calculate efficiencies, or run cross-section extraction.
- Exact model, candidate selection, missing-run policy, validation thresholds,
  file scope, and deferred ALG-002C boundary are in
  `docs/proposals/ALG-002B_simultaneous_timing_mass_model.md`.
- User condition: only configured production-LH2 runs may enter. The launcher
  is locked to `KinC_x36_4`/`LH2`; 56/57 runs are represented and only run
  6569 is an explicit missing-input exclusion.
- Status: fitter, strict manifest, validator, and synthetic regression are
  implemented. A bounded real-data smoke fit closes exactly but is
  `NOT_PROMOTABLE`; multi-start convergence, candidate selection, profile/full
  replica covariance, factorization, leave-one-run-out, timing calibration,
  and 2,000-toy coverage remain open.
- The original production timing-box subtraction remains unchanged and is not
  an ALG-002B input. A new independent comparator reports its stored per-run
  estimate beside the joint model and exposes the known 140--160 ns legacy
  histogram versus 139--161 ns shifted-sideband support mismatch.
- Technical record:
  `docs/ALG002B_simultaneous_timing_mass_20261001.md`.

## 2026-10-01: Final-workspace migration and approval boundary

- Baseline commit: `e1e7ecc4b18d359e02140c1785df64b69b3b53aa` in the
  enclosing `/w/hallc-scshelf2102/nps/singhav/nps_analysis` repository.
- Frozen source: `pi0_analysis/root_analysis_env_main/`. No build, cache,
  diagnostic, output, or documentation write is permitted there after the
  checkpoint.
- Editable target: `pi0_analysis/root_analysis_env_final/`. The 3,315 copied
  files match their committed Git blobs; the mapping is recorded in
  `manifests/baseline_copy_manifest.tsv`.
- Runtime directories are created by `scripts/setup_final_workspace.sh` and
  kept under the target. Generated validation products live in
  `validation/runtime/`; generated publication products live in
  `publication/generated/`.
- This migration stage changes organization/provenance only. It does not
  approve or install changes to selection, background subtraction, exposure,
  corrections, smearing, response, uncertainty, fit, or publication
  eligibility. Each such update requires a concrete proposal and explicit
  user approval recorded in `docs/publication_uncertainty_worklog.md`.
- Before any migrated launcher is executed, every active absolute path and
  write destination must be classified and isolated from `MAIN`. A relocation
  that can change input membership/order, RNG consumption, normalization, or
  scientific output is approval-required.

## 2026-09-29: Optional scaled-Poisson xsec objective

The no-SIMC-model extractor now accepts `--fit-objective scaled-poisson` as an
opt-in alternative to the unchanged default Gaussian fit. It treats each
selected event's final weight, including the upstream `pi0_weight` correction,
as fixed. The SIMC forward response is fixed in this likelihood pass. Supported
zero-data reconstructed rows enter the objective with a reference weight scale
estimated from nearby selected data; no pseudocount is introduced. The new
fit uses Minuit2 and preserves continuous-angle nonnegative cross sections.

The extractor output contract adds `fit_objective`, `fit_objective_value`, and
`mc_stat_treatment` fields to CSV summaries and objective/scale/optimizer
metadata in ROOT. In scaled mode, `fit_xsec_chi2` is undefined and the deviance
is recorded as `fit_objective_value`; Gaussian outputs retain their existing
numerical values. Scaled mode also writes `scaled_poisson_rows.csv`,
`scaled_poisson_weight_groups.csv`, and
`scaled_poisson_profile_intervals.csv`, plus a signed deviance residual graph
in ROOT. The fit and its validation are documented in
`src/xsec_extract/SCALED_POISSON.md`.

## 2026-09-27: Larger LT coverage maps

- The six-panel dashboard now uses a three-row, two-column figure. W/Q2 and
  xB/Q2 occupy the entire upper row with extra height; the occupancy and three
  projections use the lower rows. Axis geometry remains measured from the
  rendered Matplotlib figure, so click coordinates follow the new positions.
  This changes only display and exported figure layout, not selected events.
  The two-setting real-data notebook validation and JupyterLab browser test pass,
  including pointer clicks on both resized maps and Apply outer limits.

## 2026-09-27: Reliable outer-limit commit in LT dashboard

- The Apply outer limits button now reads all six visible input strings in one
  browser-to-kernel message. This avoids a click overtaking a pending FloatText
  update. The kernel validates the full set, updates applied limits together,
  and acknowledges success, rejection, or an invalid selected sample.
- The existing two-setting kernel validator passes with the atomic
  string-value path. A JS component check confirms that the click sends all
  six visible strings together, and a kernel message check covers success and
  rejection. The browser test clicks directly from a focused field without
  pressing Tab; it passed with the enlarged layout using isolated Playwright.
  No event-selection formula or xsec output contract changes.

## 2026-09-27: Export the LT diamond to xsec extraction

- The notebook review block now prints the selected convex `[xB, Q2]` corners
  beside the exact bin edges for `xsec_config_*.json`. `null` means no cut.
  Existing presets use `null` until the analyst copies a finalized selection.
- The JSON generator validates and orders four physical convex corners. Both
  extractors apply the same reconstructed polygon to data and SIMC before
  accumulation; selected missing-mass slice diagnostics use it too. Generated
  truth outside the reconstructed cut remains available to response guards.
- ROOT metadata now records the configured corners. Activating a polygon
  changes selected yields and extracted cross sections, so compare both
  methods against the no-polygon baseline before production use.
- Validation: both preset JSON files render with no cut; invalid polygons are
  rejected; both C++ extractors pass syntax compilation with a real polygon;
  a C++ check covers interior, boundary, exterior and disabled selection.
  The two-setting notebook validation passes and confirms its config block
  matches the live dashboard polygon. No full xsec extraction was run.

## 2026-09-24: Combined 2D mass-cut ridge orientation

The combined ellipse now estimates the weighted pi0-mass peak in populated
missing-mass slices and rotates its semi-major axis onto a robust linear fit
of that ridge. It retains the prior covariance eigenvalues, ellipse area, and
Mahalanobis threshold; the missing-mass center moves onto the ridge. Sparse or
inconsistent ridge fits keep the original covariance geometry. The existing
ellipse branch and debug metadata contract remain intact; debug text adds the
ridge fit status, bin count, slope, and center.

This is a physics-sensitive selector change. On KinC_x36_5_407, the 60 x 80
bin combined ROOT sample gives 26 ridge slices and changes the ellipse axis
slope from about -54 to -34.84 GeV/GeV. In-memory selected data events change
from the original 9,841 to 9,521; MCD selection is independent. Regenerate the combined
ROOT and xsec metadata, then compare yields and extracted cross sections to
baseline before production use. The existing xsec PDF predates this source
change and will continue to show its old boundary until regenerated.

## 2026-09-24: Iterative SIMC pi0 ratio extraction

The SIMC-model extractor now caches `full_weight/sigcm`, reweights accepted
events at vertex kinematics with a parameterized port of
`physics_pion.f:sig_param_2021/exclfit`, and fits selected Fortran coefficients
to weighted reconstructed yields with data and finite-MC sumw2. The default
shared fit uses plus.p5, plus.p7 and plus.p9. It reports before/after
chi2, nominal ndf, covariance, parameters and residual-mismatch status.
Nonpositive trial cross sections are penalized; invalid original denominators
are rejected and counted. `--normalize_mmiss` is disallowed for an
iterative absolute fit; fixed-default mode retains the historical behavior.

This changes the SIMC-model output contract: t bins and reporting points now
use physical t; tprime remains an optional cut and diagnostic. The per-phi
CSV renames the multiplier to `model_xsec_phi_center` and adds physical-t
edges/center, before-fit ratio, bin/fit status, units and Q2 diagnostic.
Empty bins are NaN with explicit status. ROOT adds `t_edges`, model parameter
and reference trees, covariance, fit status/quality, and rejection counts.
The producer adds `epsilon_i` copied from SIMC's vertex epsilon when
available, or an explicit sentinel otherwise. Existing output without it
remains usable with a counted, explicit fixed-beam approximation warning. Plot labels and metadata distinguish
the final reweighted yields from display-only original-SIMC mass shapes.

Fortran, synthetic closure, and representative-run checks are documented in
`src/xsec_extract/ITERATIVE_SIMC_MODEL.md`. The representative fit remains
a poor statistical description (chi2 63.998/9), so its status is
`converged_residual_mismatch`; no publication-quality cross section is
claimed from that test.


## 2026-09-24: Cross-section plot clarity and flag consistency

The no-SIMC-model extractor now filters undefined ratios from the residual
panel, adds a unity guide and support count, and uses nb/GeV^2 display units
across slice, coefficient and PARTONS plots. Titles distinguish reconstructed
yield, virtual-photon angular functions, selected-MC shape diagnostics and
conditional migration fractions. Migration graphics use compact indices with
the full truth mapping kept in CSV/ROOT. These are display changes only; raw
fit and covariance products keep their existing units and names.

The wrapper forwards `--no-pdf`, `--no-png` and `--no-diagnostics` to both
extractors. In the response extractor, diagnostic plots are suppressed as a
group, while core fit and optional PARTONS plots follow the requested image
formats. Individual PDFs are merged at the end with `pdfunite`; the combined
PDF retains page order without interleaved ROOT PDF-driver warnings. Plot
tests and a real KinC_x36_5_407 flag matrix are recorded in
`src/xsec_extract/VALIDATION.md`. The full-graphics and graphics-off runs
produced byte-identical fit parameters, covariance and experimental points.

## 2026-09-24: Reference GLS and residual-corrected experimental points

The no-SIMC-model response extractor now exposes `--fit-variance data` for
Eq. 5.23 and keeps the established finite-MC iteration as the default
`--fit-variance finite-mc`. The physical `full_weight/sigcm` normalization,
target divisor, event cuts, vertex linkage, three shared truth-block
coefficients, exterior nuisance blocks, fitted curves and optional GK
comparison are retained. No eta or model-dependent bin-centering factor is
applied.

Added vertex-phi response bookkeeping and separate
`experimental_points.csv`, `truth_phi_response_cells.csv`,
`experimental_point_covariance.csv` and
`experimental_point_correlation.csv`, plus matching ROOT records. Added
`migration_correlation.csv` and ROOT parameter correlation. The slice
plot overlays the new points with distinct markers; legacy `xsec` fields
continue to mean fitted phi-center evaluations. Point covariance includes
same-data GLS dependence and, in finite-MC mode, a first-order Poissonized
response term conditional on converged variances. Target covariance remains
separate. The new products change the output schema and need downstream
readers to use their explicit names. Validation uses synthetic migration,
closure and covariance checks and a separate-output representative run.
The MC covariance also propagates the response-weighted reference epsilon
through event-level cross moments exported with the truth-phi cells.
The requested KinC_x36_5_407 pipeline with `--partons` completed under
`/tmp/nps_xsec_eq531_x36_partons`: full-rank 72-parameter fit, 235 available
experimental points, and 20 finite GK comparison predictions. Independent
fit and point checks passed; the GK convolution issued numerical warnings,
so its pointwise convergence remains to be studied. Exact commands and
results are recorded in `src/xsec_extract/VALIDATION.md`.

## 2026-09-19: SIMC model multiplier at common reference kinematics

User-authorized physics change: replace the event-yield-weighted `sigcm`
multiplier in the SIMC-model extractor with a uniform phi-bin average of the
SIMC pi0 `sig_param_2021` model at common slice Q2/xB/tprime means. Keep the
existing physical data/SIMC yields and angular-bin fit basis. Derive reference
t and epsilon consistently; use the reference t for optional PARTONS overlays.
The evaluator supports W >= 2 GeV and fails below that boundary rather than
substituting for SIMC's MAID path. Generator-model changes require revalidation.

Output contract: preserve filenames and all existing per-phi CSV fields,
retain `model_sigcm_full_weight_mean` as diagnostic-only, append the actual
`model_sigcm_phi_bin_average` and five reference-kinematics columns. Add ROOT
`model_reference` tree and PARTONS `t_ref` branch; update analysis metadata.
Validation: compare the port to original Fortran, verify analytic angular
integration, and run a synthetic 12-phi-bin extraction/fit closure at varying
reconstructed kinematics. No production outputs are overwritten.

Last updated: 2026-06-02
Status: Phase 0 and Phase 1 completed; Phase 2 complete; Phases 3, 4, 5, 6, and 7 completed for current scope with regression/interface checks

## 1. Purpose

This document defines the staged refactor and workflow hardening plan for the NPS neutral-pion analysis framework located at:

`/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final`

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

### Plot diagnostics update (2026-09-07)

- Producer-owned before/after histograms use the existing cut-debug PDF: 1D
  overlays with upper-right legends, 2D before/after panels side by side with
  shared spatial axes. NPS position maps use independent before/after color
  scales; other maps share color scales. HMS, cluster and pair stages remain explicit.
- Cluster position axes stay x=[-34,34], y=[-40,40] cm; the older `h_clustXY`
  diagnostic gains edge bins at its original 2 cm width. Outside-view labels
  are removed; flow-bin data and mass-window metadata remain available.
- Cluster x/y cut diagnostics now use 2 cm bins in 1D and 2D (34 x bins,
  40 y bins), including the after-cut clones. The existing post-HMS/pre-NPS
  multiplicity histogram is drawn as an `nclust` panel inside `cut_debug_run<RUN>.pdf`;
  the initially added standalone PNG output has been removed from the producer.
- Corrected `nclust` diagnostic filling to use input multiplicity before capping.
  Unit-width bins cover 0 through 1080, with occupancy-based display limits.
  The 20-cluster processing cap and all event selections remain unchanged.
- Added `nclust` overlay labeled "After NPS cluster cuts": recount the full
  input list using energy, position, timing and dead-block cuts; retain zero-count events and the
  same HMS population. The extra histogram is owned by plot diagnostics.
- The `nclust` legend sums actual cluster multiplicities before/after cuts;
  histogram entries and the y-axis continue to represent events.
- `nclust` uses independent y scales: after cuts on the red left axis, before
  cuts on the blue right axis. Only display copies are scaled; stored counts,
  errors, cluster totals and physics selection are unchanged.
- `nclust` also uses independent x ranges: after cuts on the red bottom axis,
  before cuts on the blue top axis. Both include zero and retain occupied bins.
- Cluster-energy and corrected missing-mass cut-debug views start at exactly
  0.4 GeV, preserving full-population legend counts. The cluster overview energy
  panel and combined missing-mass overlay use the same lower display limit.
- Fixed exit-time segmentation violation: `h_nclusters` was deleted by both
  cluster and cut-debug cleanup. It now has one cleanup owner (cut-debug).
  Full run 5237 with the good-event cut passes; its summary matches the pre-fix
  output. Event processing and plotting settings are unchanged by this fix.
- Existing physics and mass canvases compare the retained exclusivity selectors.
  A diagnostic ellipse qualification rejects undersized fit subsets and reports
  de-correlation fallback without changing existing event flags or weights.
- The -31.95 reference line is display-only. Missing-mass correction changes,
  new numerical band cuts, efficiency changes and production selection changes
  remain out of scope.
- Collector gives known multi-panel canvases more area and preserves both PNG
  and PDF representations. No new physics plot family is introduced.
- New contract: additional diagnostic histogram keys, numeric qualification
  fields and ROOT status strings. Existing filenames/tree schema remain compatible.
- User correction: `collect_kinematic_plots.py` is the only user-facing plotting
  command. Removed the standalone range-preparation helper and its CSV/env setup;
  existing analysis automatically chooses views from its pre-cut histograms.
- Details and reproducible validation: `docs/plot_diagnostic_updates.md`.

## 2026-09-16: GK convention and global migration extraction

- Replaced the monolithic extractor with one steering translation unit and
  responsibility-specific `xsec_*.h` implementation headers; detailed map and
  physics comments are in `src/xsec_extract/README.md` and the source headers.
- Corrected native GK/PARTONS conversion to divide the electron observable by
  full Hand Gamma before Fourier projection. The old extra 2 pi is removed.
  `partons_electron_flux_xbq2` consequently changes meaning to full Gamma.
- Removed both missing-mass normalization flag aliases and all fit-area
  rescaling. Data remain divided by `tgt_contam=0.584`, with 0.014 uncertainty
  exported as a separate correlated covariance. Shape-only diagnostic copies
  may still be area matched. Charge/efficiency and physical SIMC normalization
  retain their previous meaning.
- Added mandatory matched generated kinematics, generated-basis/reconstructed-row
  event response, and simultaneous unregularized full-rank SVD. Six disjoint
  outside-domain regions receive free U/LT/TT parameters when populated; this
  is an explicit coarse nuisance model requiring systematic variations.
- Reference provenance: Ali p.142, Eqs. (5.2)-(5.5), checked through indexed
  primary-PDF text; Defurne p.97 relies on prior repository reading notes because
  fresh PDF retrieval failed. Source links and this limitation are recorded
  explicitly in the extraction guide; no thesis-prescribed guard model claimed.
- Added finite-MC outer-product covariance and converged GLS iteration. Reject
  missing support, rank deficiency, invalid selected weights, inconsistent run
  charge/scale, unmatched raw events and failed convergence. Do not clip fitted
  negative cross sections, manufacture empty-row counts or apply smoothing.
- Output contract: coefficients now describe generated bins; slice chi2/ndf
  are global and labeled. New CSV/ROOT records retain the response, truth/row/
  parameter mapping, complete covariance, MC moments, singular values and
  exclusion/physics metadata. Historical filenames and leading CSV fields remain.
- Added tests under `tests/`: GK conventions, response solver/covariance and
  independent NumPy/event-split closure, plus optional native PARTONS smoke.
  Test runner and `src/xsec_extract/VALIDATION.md` preserve reproducible commands.
- KinC_x60_4b validation: 225 rows, 72 parameters, rank 72, condition 244.221,
  chi2/ndf 216.718/153, 9 MC iterations. Independent disjoint-MC input recovery:
  chi2 68.02/72; residual chi2 160.10/156. Five generated bins have negative
  angular regions with broad errors and remain flagged, not declared physical.
- Remaining approximations: coarse exterior/within-bin shape, nominal-energy
  generated epsilon, frozen adaptive boundaries, omitted zero-variance rows,
  Poissonized MC and conditional fit covariance. No detector/background/radiative
  retuning or converged GK agreement is claimed. Plot smoke completed; ROOT
  emitted object-state warnings for combined PDF export (readable 32-page file).

## 2026-09-17: Recover a global xsec fit after Q2/xB exclusions

- Try the full global migration fit first, then exclude unsupported Q2/xB
  groups (all reconstructed t'/phi rows) and globally refit the remainder.
  For unlocalized rank/convergence failures, search smaller subsets in
  decreasing retained size, with deterministic ties by excluded group index.
- Preserve every truth contribution feeding retained rows as a free fitted
  coefficient, including excluded-bin feed-in as nuisance. Never replace
  migration with independent per-bin fits, fixed backgrounds or regularization.
- Failed attempts roll back all fitted slice values; successful retained bins
  share one covariance and chi2/ndf. Binning, weights, and target correction
  remain frozen across retries. Selection is data dependent and covariance
  is conditional on the chosen subset; recovery is not physics validation.
- Append fit status/reason columns to historical CSVs; add `fit_status.csv`
  and `fit_attempts.csv` and corresponding ROOT trees. Excluded values are NaN
  with fit_ok=0; parameter mapping marks nuisance columns explicitly. All-failed
  runs save diagnostic data/response/status, skip PARTONS, then exit nonzero.
- Add deterministic recovery tests covering full success, missing MC support,
  cross-bin covariance, excluded truth feed-in, rank failure, unsupported truth,
  convergence exhaustion, and all-failed CSV/ROOT serialization.

## 2026-09-17: Whole-extraction physics review and helicity diagnostic

- Reviewed the 22-page KinC_x60_4b PDF and traced its exact smeared input to
  `/volatile/hallc/nps/singhav/nps_smearing/smear_x60_4b_test/smearing_output/KinC_x60_4b/root/simc_pi0_analysis_output_smeared.root`.
  A scratch rerun reproduces its migration matrix, coefficients and covariance;
  the workspace-default SIMC file is a different response.
- User approved correcting the diagnostic weighted-yield asymmetry variance
  using derivatives of (Yplus-Yminus)/(Yplus+Yminus). Bins without a valid
  two-helicity estimate are omitted; an entirely invalid panel says unavailable.
  The existing combined tree has no negative helicity labels. This change does
  not reconstruct missing helicity information or alter the unpolarized fit.
- Added a ROOT-independent helper and regression test for weighted yields,
  the counting limit, common scaling, signed subtraction and missing support.
  The standard test runner includes this test.
- Added explanatory comments for charge normalization, conditional errors,
  adaptive binning, SIMC response units, azimuth conventions, nominal-energy
  epsilon, exterior nuisance sensitivity, and legacy simulation-sum fields.
- Weak exterior constraints, fitted-background covariance, signed subtraction
  clipping upstream, matched final data/MC acceptance, radiative treatment and
  branching normalization require separate physics decisions/validation. No
  positivity constraint, target-factor change, or response retuning is implied.

## 2026-09-17: Thesis-based migration diagnostics

- User requested migration plots using Defurne HAL tel-01281332v1 and Ali
  OSTI 1784736. Fresh rendered-page inspection identifies Ali Fig.5.3 (p.142)
  as U/LT integrated-response matrices; TT is added using the existing basis.
  Defurne Fig.3.9 (p.51) is a DIS vertex-coverage visualization, not a matrix.
- Add `xsec_plot_migration.h`: U/LT/TT response heatmaps and separately
  normalized selected-event migration/origin fractions. Phi varies inside
  reconstructed t'/Q2/xB blocks; outside-domain truth guards are explicit.
  Signed LT/TT are never interpreted as probabilities or normalized in the fit.
- Add `xsec_plot_migration_coverage.h`: generated/reconstructed Q2-xB coverage
  for the same accepted MC. Color is count density, explicitly not unavailable
  vertex beam energy. Histograms extend to retain truth feed-in; display clones
  divide by bin area while stored ROOT histograms retain counts.
- Append three pages to the combined PDF and add `migration/` PDF/PNG files
  and ROOT objects. `--no-diagnostics` skips these additions; ROOT histograms
  remain available when only image output is disabled. Existing CSV schemas,
  event selections, response cells and fit parameters are unchanged.
- Validation compares every plotted response cell to existing CSV exports,
  both normalization directions and coverage event totals, plus baseline fit
  coefficients/covariance. Synthetic checks cover guards and empty support.

## 2026-09-17: Optional physical positivity in cross-section extraction

- User authorized adding a nonnegative-cross-section fit mode. Add
  `--positive-xsec` / `--no-positive-xsec` and `NPS_XSEC_POSITIVE_XSEC=0|1`
  to the extraction executable and wrapper; the default remains unconstrained.
- Constrain the complete U/LT/TT angular function, keeping interference terms
  signed. Apply the constraint to every active truth block, including exterior
  and excluded-group nuisance blocks. Use exact angular minima and each
  block's maximum observed generated epsilon, covering every smaller epsilon
  under the existing bin-constant coefficient model.
- Add `xsec_positive_solver.h`: full-rank SVD followed by a convex quadratic
  projection with continuous-angle constraint generation. Repeat it inside
  each finite-MC variance iteration; do not clip fitted parameters or curves.
- Add `positivity_diagnostics.csv` and ROOT tree, including epsilon maxima,
  exact angular minima, constraint flags and numerical tolerances. Add
  `migration_curvature_inverse.csv` and explicitly named ROOT curvature matrix.
  Add mode/boundary/covariance-status metadata and positive-mode fit_scope.
- When any constraint binds, global statistical/MC covariance and derived
  fitted errors are NaN and plots label central estimates without intervals.
  Unconstrained curvature is diagnostic only; boundary confidence intervals
  require dedicated profile/toy inference. Nominal ndf is descriptive only.
- Keep event selection, response construction, target normalization and
  existing unconstrained estimates unchanged. Document semantics and commands
  in `src/xsec_extract/POSITIVITY.md`. Controlled validation and the separate
  constrained PDF are under the sibling workspace's
  `audit_reports/20260917_xsec_positivity/`; no baseline PDF is overwritten.

## 2026-09-24: Constrained-fit plot errors and nanobarn display

- KinC_x36_5_407 at 0.7-1.05 GeV missing mass activates the continuous-angle
  positivity boundary; the ordinary covariance is unavailable, explaining
  missing bars in the user's 25-page PDF. Added 256 deterministic conditional
  constrained refits with the converged weighted-yield row variance and the
  full finite-MC variance iteration. Plot their sampling SD separately from
  the unavailable inverse-information covariance. Export parameter and
  Eq. 5.30 point toy covariance, including cross-bin terms, to a separate
  CSV and ROOT matrices. Toy response moments and input variances are frozen;
  this does not establish boundary-interval coverage.
- Plot structure functions and angular cross sections in nb/GeV^2 and
  nb/(GeV^2 rad). Display diagnostic response heatmaps as yield per nb/GeV^2
  coefficient without changing raw ROOT/CSV response. Show generated -t'
  bin spans alongside response-weighted vertex-mean markers; the three-term
  panel shows the shared horizontal extent once on U to reduce clutter.
- Keep full_weight/sigcm, the independent target correction, eta=1,
  vertex/reconstructed linkage, and GK as an after-fit comparison.
- Improve migration diagnostics: add a separate, readable log10(1+selected MC
  count) support page alongside the two conditional fraction directions,
  outline matching truth/reco cells,
  and label published bins with physical -t' intervals when there is one
  Q2/xB bin. Guard labels name the six exterior regions. These are display
  changes only; the fitted response and stored ROOT histogram values are raw.

This file is a living implementation log, not a static design note.

Update it whenever:
- a new file is added,
- an output contract changes,
- a phase is completed,
- a validation fails,
- a physics-sensitive decision is made,
- or a previously unknown dependency is discovered.


## 2026-09-25: Explicit cross-section bins

- Both extractors now take phi, Q2, per-Q2 xB, and tprime edges from
  `src/xsec_extract/xsec_config.h`; the iterative SIMC method additionally
  takes its physical-t edges there. Counts and selection bounds derive from
  those vectors. Data/SIMC quantiles no longer choose internal boundaries.
- The previous outer limits remain. Internal tprime and physical-t defaults
  are evenly spaced; earlier runs used sample-dependent quantiles and are
  therefore not directly comparable. Saved ROOT/CSV edge metadata and the
  extracted values can change. The legacy global xB metadata vector contains
  the first Q2 row; per-Q2 edge rows remain authoritative.
- Removed command-line bin overrides from both extractors and the wrapper.
  Invalid edge arrays or incompatible per-Q2 xB rows stop extraction before
  event processing. No output filenames or schemas changed.
- Validation: both ROOT translation units pass syntax-only compilation in
  the Hall C environment. The fixed-edge response recovery test and both
  SIMC-model synthetic closure tests pass with temporary test config headers.

## 2026-09-27: LT bin-selection notebook interaction and validation

- Preserve the original untracked notebook/module, including outputs, under
  `src/xsec_extract/backups/xsec_bin_20260927_145011/`. Keep notebook cell IDs
  and the representative-coordinate formulas; clear outdated output state.
- Replace trace-dependent Plotly editing with Matplotlib SVG plus an anywidget
  axis-coordinate click layer (`xsec_bin_canvas.js`). Stable numbered corners,
  numeric editing, queued pointer events, undo history, and explicit invalid
  states replace marker/callback patches. PDF/SVG use the same figure; SVG
  glyph outlines avoid browser symbol-font substitution (t-prime was affected).
- Apply a single convex xB/Q2 mask to all data settings and the existing SIMC
  representative sample. Project densely sampled edges into W/Q2 and check
  ROOT W (maximum allowed difference 5e-4 GeV; observed 4.44e-16 GeV).
  Retain the physical missing-mass input sample before the adjustable tprime
  window so expanding that window can recover data. No production cuts change.
- Quantiles remain data-only, pooled across settings with the original weight
  convention; positive weights only for yield quantiles. Phi stays uniform.
  Downstream cells consume the displayed edge arrays, including conditional
  manual xB arrays. Incomplete, invalid, and stale proposals cannot be exported.
  Existing single-setting/SIMC/manual-Q2-and-xB JSON safeguards remain in force.
- Full coverage stays visible, with excluded cells faded. Per-setting counts,
  all-setting overlap, full Q2/xB/tprime/phi occupancy, and sparse-setting flags
  are linked. Figure exports add user-chosen `.pdf`/`.svg` paths; the scientific
  JSON schema and ROOT/production output contracts are unchanged.
- Add reproducible kernel and real JupyterLab/Chromium test scripts and
  `README_xsec_bin_helper.md`. Test artifacts live in `output/lt_bin_validation`.
  The required Python kernel passes both supplied files and the single-setting
  configuration; independent masked counts and every full bin agree. Frontend
  tests cover creation from both views, labels 1-4, corner-4 moves, undo, clear,
  numeric edits, queued rapid clicks, and changes to modes/bins/limits.
  This does not establish behavior in VS Code, classic Notebook, or other browsers.

## 2026-09-30: Joint LT/TT with independent setting U

- `run_joint_xsec_fit.py` now fits shared sigmaLT and sigmaTT per truth block
  while refitting a separate sigmaU for each setting. It does not fit shared
  sigmaT/sigmaL or fix U to earlier single-setting estimates. The C++ solver
  uses setting-specific U columns and shared LT/TT columns in one stacked fit.
- Read U/LT/TT response moments directly from `migration_response_cells.csv`;
  if absent, select the T/LT/TT submatrix of the four-component export (T and
  U have identical constant angular basis). Keep all MC cross moments and
  per-setting normalizations. Truth-phi exports are no longer required.
- Positivity uses each setting's own event epsilon maximum and independent U
  with shared LT/TT; boundary fits continue to suppress Gaussian error claims.
  Unsupported guard U parameters are omitted for that setting; supported
  guard LT/TT are shared. Published U parameters always remain independent.
- Internal problem format changes to `joint_xsec_v3`. Output filenames stay
  unchanged. `joint_parameters.csv` now contains U/LT/TT and appends
  `setting_index` (zero-based CLI order for U, -1 for shared LT/TT).
  `joint_positivity.csv` appends `setting_index` and reports one row per
  supported setting/block. Covariance indices follow parameter_index; ROOT
  vectors use that same ordering. Manifest/ROOT/summary describe the new model.
- Add `tests/test_joint_xsec_fit.py`: temporary synthetic response inputs,
  different U values at equal epsilon, two/three settings, one-setting guard
  support, independent NumPy finite-MC fit/covariance comparison, per-setting
  epsilon constraints, boundary positivity in both variance modes, and
  unsupported-column rejection. This validates algebra, not detector bias.
- Validation setup: standalone compilation initially lacked generated
  `xsec_config.h`; the test generates it using the same renderer as the driver.
  The first regression run passed five numerical cases but used an incorrect
  expected error substring for rank failure; corrected to the actual solver
  diagnostic `unsupported truth coefficient`.
- Final validation: ROOT 6.30.04 compilation and all seven regression checks
  pass. Python syntax and changed-file whitespace checks pass. Reproduce from
  this repository with:
  `csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 tests/test_joint_xsec_fit.py'`.
  Production extractions were not rerun or overwritten.

## 2026-09-30: Fit-independent joint input preparation

- Add `--prepare-setting CONFIG VERTEX_SIMC` to the joint driver. It reads
  data/SIMC paths from each preset, prepares event yields and response moments,
  then runs only the joint LT/TT + per-setting U fit. The existing individual
  pipeline remains available and is not called by this workflow.
- The compiled extractor's `--prepare-joint-inputs` mode returns after event
  accumulation and the existing target correction, before `fit_slices()`,
  recovery, toys, PARTONS and plots. `xsec_joint_inputs.h` exports only raw fit
  inputs; no individual fit result, positivity diagnostic or fitted ROOT file
  is needed. It refuses nonempty destinations.
- New prepared-input contract: `joint_input_metadata.txt` identifies
  `prepared_joint_inputs_v1` / `fit_objective=not_run`; `joint_input_slices.csv`
  saves bin bounds, and migration row/truth/response CSVs save corrected yields,
  variances, event epsilon maxima and all U/LT/TT MC cross moments. Outputs live
  in a separate `OUT_DIR_inputs` directory, overridable with `--inputs-dir`.
  Effective config snapshots and a preparation manifest make those inputs
  reusable with `--setting`. Original configs and individual outputs are not
  overwritten. Completed prepared inputs survive a joint-fit failure.
- `--binning-config` explicitly copies all bins/diamond into effective presets;
  common selection bounds can be set with `--mmiss-select/lower/upper` during
  preparation. Without this option mismatched bins still fail before event IO.
- Remove dependence on prior extraction positivity, objective, variance and
  PARTONS flags: these do not affect the accumulated response or observations.
  The joint manifest explicitly records its own Gaussian objective and the
  upstream input stage/objective. Scaled-Poisson joint fitting is not added.
- Update the existing joint-fit README section and add an end-to-end synthetic
  ROOT test, `tests/test_joint_xsec_preparation.py`, for target normalization,
  common binning without preset mutation, no individual fit artifacts,
  coefficient recovery, refitting with different options and overwrite refusal.
- Validation passed: all seven existing joint numerical checks; complete
  synthetic ROOT preparation/joint-fit/refit test; final extractor syntax
  compilation under ROOT 6.30.04; Python syntax and tracked diff whitespace.
  The test recovered U=(6,10), LT=0.3 and TT=-0.2 with no individual-fit
  artifacts. The documented production input and ellipse-geometry paths exist;
  production processing was not run. Reproduce with
  `python3 tests/test_joint_xsec_preparation.py` and the numerical command above.

## 2026-09-30: Post-fit nominal-epsilon L/T separation

- Complete the interrupted change in `run_joint_xsec_fit.py`: after fitting
  separate U per setting and shared LT/TT, separate `U = T + epsilon_nominal L`
  independently in each truth block. Two settings use the exact intercept and
  slope; additional settings use GLS with their full within-block U covariance.
- Use preset beam energy, central HMS momentum and angle (massless electron)
  for fixed nominal epsilon, or one `--nominal-epsilon` override per setting.
  The already-added x36_5_407/x36_4 central momenta, 4.637/2.562 GeV, give
  epsilon 0.7115005066172831/0.5163181761587335. Overrides are forwarded from
  preparation into fitting. Event response/positivity epsilon is unchanged.
- Add `joint_separated_parameters.csv` with T/L/LT/TT and explicit status;
  `joint_separated_covariance.csv` carries transformed stat-plus-MC covariance,
  including cross-block and LT/TT correlations. Its indices address the new
  parameter table. `joint_lt_separation.json` and the manifest's `lt_separation`
  entry record epsilon inputs, convention and block diagnostics. Existing
  U/LT/TT CSV and ROOT contracts are preserved; no extra ROOT objects are added.
- Fixed epsilon uncertainty and target-factor uncertainty are not propagated.
  There are no additional T/L positivity constraints. Rank-deficient epsilon
  or one-setting guard blocks produce NaN T/L with status. Boundary-constrained
  two-setting fits yield central values but no Gaussian errors; more-setting
  GLS needs valid covariance. Diagnostic curvature never supplies uncertainties.
- Add `tests/test_joint_lt_separation.py` for exact recovery, correlated GLS,
  full covariance transformation, guard support, epsilon degeneracy, boundary
  errors, nominal kinematics and invalid overrides. Extend the synthetic ROOT
  preparation/refit test to check nominal separation and epsilon overrides.
- Validation passed: six post-fit numerical tests, all seven existing solver
  checks under ROOT 6.30.04, and the synthetic preparation/joint-fit/separation/
  saved-input-refit test. Python syntax, CLI help and changed-file whitespace
  checks passed. No production processing was run. Reproduce from repository
  root with `python3 tests/test_joint_lt_separation.py`,
  `python3 tests/test_joint_xsec_preparation.py`, and
  `csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 tests/test_joint_xsec_fit.py'`.

## 2026-09-30: Complete joint-fit plot report

- Treat `run_xsec_pipeline.sh` as a read-only reference. Its original 757-line
  contents were recovered from editor history after an interrupted patch,
  verified byte-for-byte, and all attempted wrapper changes were removed.
  Reference SHA256: a59ee141c7f3bff6daae1d2fe4f26b9eaf38a17a52e32b75f4dd2ecb4e67d82b.
  AGENTS.md now requires reference-only files to remain unchanged and authorized
  edits to use recoverable staged copies with atomic replacement.
- Add `joint_xsec_plots.py`: private native builds use the same plot functions
  as the reference workflow, with no edits to or invocation of its shell.
  A `--joint-plot-input` mode in the C++ extractor imports per-setting U/LT/TT
  marginal joint covariance, validates rebuilt yields/response/predictions,
  and bypasses individual fits and refit toys. It preserves standard extraction
  behavior when that mode is absent. Optional native PARTONS comparisons use
  the established library setup and GK06/GPDGK19 conventions.
- Plot by default after persisting numerical results; `--no-plots` skips it.
  `--plot-only JOINT_OUTPUT` adds a report without fitting; `--partons` enables
  model curves with configurable integration counts. Raw inputs must remain
  accessible to rebuild global/mass/epsilon/migration/helicity diagnostics.
- New artifacts: `plots/settings/<kinematic>/` contains all native PDF/PNG
  families and renderer diagnostics; `plots/joint/` contains Rosenbluth,
  separated coefficient, covariance ellipse/correlation, residual/fit-quality
  and positivity plots. `all_joint_xsec_plots.pdf` combines pages without
  duplicating per-setting combined PDFs. `plots/plot_manifest.json` records
  completeness, page inventory and uncertainty limitations.
- Joint coefficient and fitted-curve errors use the joint covariance. Native
  residual-corrected points remain explicitly central-only because the local
  influence calculation excludes other settings. Boundary-constrained errors
  are unavailable; no individual-fit toy spread is substituted. Target and
  epsilon uncertainties remain outside the joint covariance.
- Add synthetic ROOT plotting coverage test against the ordinary native
  extractor with ellipse diagnostics and optional PARTONS smoke integration.
  Check imported coefficients/errors, every native plot family, joint figures,
  PDF page count, output PNGs and unchanged reference-wrapper checksum.
- Validation passed with explicit exit code 0: `python3
  tests/test_joint_xsec_plots.py --partons` (100 warmups/1000 calls for smoke
  testing only) covers 21 native pages per setting plus each setting's combined
  PDF, eight additional joint figures and 50 combined report pages. Checks
  include finite-MC covariance import, no individual fits during joint
  rendering, ellipse diagnostics, PARTONS, boundary/degenerate-epsilon plots
  and completed-report overwrite refusal. Data-only plot-only rendering,
  six L/T numerical tests, Python syntax and final reference checksum checks
  also passed. Production processing was not run.

## 2026-09-30: Preserve diamond schema during common-binning preparation

- Fix the joint driver's `--binning-config` override: normalized diamond
  corners are tuples, but the JSON preset validator requires list pairs.
  Convert copied corners back to lists before revalidating the effective
  config; preserve all coordinates and retain null for disabled diamonds.
- The earlier synthetic fixture disabled the diamond and missed this case.
  Its preparation/refit regression now uses a nonempty diamond and verifies
  the saved common selection. The complete regression passed with exit 0.
- The user's exact x36_5_407/x36_4 command passed config/path preflight,
  stopping before ROOT execution or output creation. Both presets and the
  reference `run_xsec_pipeline.sh` passed unchanged-checksum checks.
  Python syntax passed; no production processing was run. Reproduce with
  `python3 tests/test_joint_xsec_preparation.py` from repository root.

## 2026-10-01: Independent-grid harmonic forward extraction

- Preserve existing extraction defaults and outputs; add explicit
  `--prepare-forward-inputs` cache export before rectangular cuts, retaining
  mass/diamond cuts, matched original event IDs and absolute normalization.
- Add independent reconstructed/truth/publication JSON configuration, all-row
  weighted scaled-Poisson fitting, continuous angular positivity, full nuisance
  profiles, data/MC event refits and paired common-integral binning comparisons.
- New cache and result directories refuse replacement. Persist completion
  state, effective config, input/source hashes, runtime versions, full numerical
  response and conditional uncertainty status. No publication-ready claim.
- Replace covariance Cholesky with the direct SVD factor in both constrained
  C++ solvers. Reject forced-positive-definite Minuit confidence covariance.
- Candidate x36_4 buffer is explicitly a validation configuration, not a new
  default. Background-parameter covariance was found in per-run combbg ROOT
  outputs; full event/purity-fit joint uncertainty remains outside this release.
- Technical/physics contract and commands: `src/xsec_extract/FORWARD_EXTRACTION.md`.
  Final validation results are recorded there and in the task worklog.
## 2026-10-01: ALG-002A combined timing-background shadow contract

- Approved scope: opt-in combined-statistics timing fit for `KinC_x36_4`, with
  run/mode/multiplicity identity retained and efficiencies excluded.
- `src/background_fit/` reads only ALG-001 `raw_observation` and refuses
  canonical `output/` destinations.
- Shadow outputs include manifest, run parameters, run/mass component yields,
  true-coincidence spectrum, cell residuals, curvature, provenance, run
  coverage, and validation gates.
- No `pi0_weight` is written; production analysis, combine, smearing, and
  extraction paths are unchanged.
- Validator exit 3 is required while run coverage, identifiability,
  calibration, coverage, or goodness gates remain unresolved.
- Detailed record: `docs/ALG002A_combined_timing_background_20261001.md`.
