# Smearing Pipeline Refactor Plan

Status: living document. Created before major code edits.

## Core Objective

Refactor smearing so photon energy response uses one active default mean model:

```text
E_mean(E) = a + b*E + c*ln(E_safe)
```

where `E` and `E_mean` are in GeV, `a` and `c` are GeV, `b` is dimensionless, and `E_safe=max(E,E_floor)`.
This is an additive reconstructed-energy mean model, not a multiplicative scale factor. Energy resolution is applied once around this mean.

## Current-State Audit

### Files

- `scripts/nps_sim_smearing_new_try.C`
  - Main fitter executable requested for refactor.
  - Fits per-section photon energy mean coefficients, energy resolution, optional position resolution, and optional electron momentum scale.
  - Writes `out_smear.root`, `section_map.csv`, `chi2_scans.pdf`, and `out_smear_interpolated.root`.
- `scripts/simc_pi0_analysis.C`
  - Produces unsmeared and smeared SIMC trees.
  - Consumes `section_map.csv` and interpolated ROOT maps.
- `run_smearing_pipeline.sh`
  - Orchestrates compile, fitter run, and downstream SIMC regeneration.
  - Currently compiles `scripts/nps_sim_smearing_new.C`, not requested `scripts/nps_sim_smearing_new_try.C`.

### Smearing Definitions In Fitter

- Energy mean / response:
  - Config at `scripts/nps_sim_smearing_new_try.C:300-325`.
  - Current active flag: `ENABLE_ENERGY_DEPENDENT_MU=true`.
  - Formula in `smearPhoton`: `E_sc = mu_a + mu_b*E_safe + mu_c*log(E_safe)`.
  - Legacy disabled path still exists conceptually: `E_sc=b*E`; controlled by `ENABLE_ENERGY_DEPENDENT_MU`.
  - Earlier artifact risk found: coarse/fine scans varied only scalar `mu_b` with `mu_a=0, mu_c=0`, so `a` and `c` could remain zero if the single zero-seed MIGRAD start failed, stalled, or did not improve the scalar seed.
  - Current mitigation: deterministic nonzero `a/c` MIGRAD multistarts are enabled, and each candidate is compared by chi2 so a nonzero seed can win even when the minimizer status is conservative.
- Energy resolution:
  - `computeEnergyResolution(E_scaled, sigma, res_A, res_B, res_C)`.
  - Simple model: `sigma_E=sigma*sqrt(E_mean)`.
  - 3-term model: `sigma_E=sigma*E_mean*sqrt(A^2/E + B^2/E^2 + C^2)`.
  - PDF shape can be Gaussian or Landau. Default Gaussian.
- Position smearing:
  - `sigma_pos` additive x/y smearing in cm.
  - Optional energy dependence `sigma_pos_eff=sigma_pos*sqrt(E0/E)`.
  - Current fitter default: position smearing disabled.
- Pair-energy correction:
  - Optional `ENABLE_MGG_LINEAR_ENERGY_CORRECTION`.
  - Applies extra pair-level energy correction after per-photon smearing.
  - Currently disabled. If enabled, it is a second energy-response path and must remain off unless explicitly studied and documented.
- Electron momentum scale:
  - Optional `p_e_scale` affects missing-mass calculation only.
  - Not photon smearing, but can be degenerate with photon energy response.

### Hard-Coded Parameters

- Fitter config block contains hard-coded:
  - input physics constants: beam energy, y mispointing.
  - energy model toggle, bounds, seeds, log floor.
  - sigma bounds, model choice, Gaussian/Landau choice.
  - position model toggle and bounds.
  - electron momentum scale toggle and bounds.
  - output directory and filenames.
  - optimizer controls and histogram binning.
- Downstream `simc_pi0_analysis.C` has hard-coded:
  - smearing file paths in `output/plots/x60_4b/production_wfpi0`.
  - old `ENABLE_ENERGY_DEPENDENT_MU=false`.
  - obsolete energy-dependent mu formula `mu0*((E0/(E+Eshift))^b+1)`.
  - `ENABLE_ENERGY_DEPENDENT_SIGMA_POS=true`, inconsistent with fitter default false.
  - output file paths.
- Shell pipeline has hard-coded:
  - absolute data/SIMC input paths.
  - output directory and fitter binary name.
  - grid geometry and `Nsmear`.

### Duplicate Or Conflicting Energy-Smearing Paths

- Fitter supports both legacy scalar `b*E` and required `a+bE+c ln(E)` via flag.
- Downstream consumer applies old multiplicative `mu_eff(E)` model when enabled and currently reads only one `mu` field.
- Downstream CSV parser assumes legacy header:
  - `best_mu,best_sigma,best_sigma_pos_cm,...`
  - Actual fitter writes:
  - `best_mu_a,best_mu_b,best_mu_c,best_sigma,best_sigma_pos_cm,...`
- Downstream interpolated maps load only `h_mu_interp`, but fitter also writes `h_mu_a_interp` and `h_mu_c_interp`. Consumer ignores `a,c`.
- Optional pair-level mgg energy correction is disabled but is another energy correction path if activated.
- Optional electron momentum scaling is separate physics but can mimic photon energy response in missing mass.

### Plots Depending On Smearing

- Per-section pages in `chi2_scans.pdf`:
  - data/unsmeared/smeared comparisons for `M_gg`, `M_miss`, `(p_target+gamma gamma)^2`.
  - chi2 scans for `mu_b`, `sigma`, optional `sigma_pos`.
  - pull plot for `M_gg`.
  - text summary of parameters.
- Summary pages:
  - combined all-section observable comparisons.
  - mismatch attribution.
  - parameter maps for `mu_a`, `mu_b`, `mu_c`, `sigma`, `sigma_pos`, `p_e_scale`, `chi2/ndf`.
  - `E_mean(E)` curves.
  - `E_mean(x,y)` maps at fixed energies.
  - section statistics and final text report.
- Interpolated ROOT maps:
  - `h_mu_a_interp`, `h_mu_interp`, `h_mu_c_interp`, `h_sigma_interp`, `h_sigma_pos_interp`.

### Downstream Outputs Consumed

- `section_map.csv`
  - Must contain named columns and model metadata.
  - Current consumer uses positional parsing, causing wrong columns after `best_mu_a/b/c` change.
- `out_smear_interpolated.root`
  - Must include all model parameter maps needed by downstream smearing.
  - Current consumer ignores `h_mu_a_interp` and `h_mu_c_interp`.
- `simc_pi0_analysis_output_smeared.root`
  - Produced by downstream macro.
  - Should carry metadata documenting smearing model, parameter source, and mode.
- `chi2_scans.pdf`
  - Main human diagnostic.

### Assumptions In `simc_pi0_analysis.C`

- Section smearing is only supported mode; `event` mode warning falls back to section.
- Smearing is optional; if files missing, macro silently produces only unsmeared output.
- Uses `section_map.csv` plus interpolated ROOT maps.
- CSV positional schema is legacy and wrong for current fitter output.
- Energy response is treated as multiplicative `mu*E`, with optional obsolete nonlinear factor.
- Interpolated `h_mu_interp` treated as scalar scale/linear coefficient.
- Position smearing enable/energy-dependence flags must match fitter but currently do not.
- `p_e_scale` averaged across loaded sections when enabled.

### Assumptions In `run_smearing_pipeline.sh`

- ROOT and `root-config` in `PATH`.
- Absolute production input files exist.
- Compiles `scripts/nps_sim_smearing_new.C`, not `scripts/nps_sim_smearing_new_try.C`.
- Does not remove or isolate old outputs before rerun.
- Does not verify `section_map.csv`, interpolated maps, PDF, or smeared SIMC file after each stage.
- Does not pass explicit smearing paths to downstream macro.

## Physics Requirements

- Default energy mean convention:
  - `E_mean(E) = a + b*E_safe + c*ln(E_safe)`.
  - Represents reconstructed/smeared-distribution mean energy in GeV.
  - Not a scale factor. Do not multiply original energy by it.
- Dimensional consistency:
  - `E`, `E_safe`, `E_mean`: GeV.
  - `a`: GeV.
  - `b`: dimensionless.
  - `c`: GeV because `ln(E_safe / 1 GeV)` is dimensionless by convention.
  - Code will use `ln(E_safe)` with `E_safe` numerically in GeV; docs must state reference unit is 1 GeV.
- Low-energy behavior:
  - Use `E_safe=max(E,E_floor)`.
  - Default `E_floor=0.2 GeV`.
  - Warn/count when input energy is below floor.
- Valid energy range:
  - Default model validity: `0.2 <= E <= 8.0 GeV`, pending actual SIMC/data distribution checks.
  - Outside range: clamp log input, preserve event but count warning.
- Bounds:
  - `a in [-0.5,0.5] GeV`.
  - `b in [0.95,1.12]`.
  - `c in [-0.5,0.5] GeV`.
  - `sigma >= 0`, default range from current config.
  - `sigma_pos >= 0`.
- Physical guards:
  - Clamp nonpositive `E_mean` and random energy draws to `NONPOSITIVE_CLAMP`.
  - Reject nonfinite parameters.
  - Resolution evaluator must never return negative/nonfinite sigma.
- Avoid double application:
  - Energy response once, then resolution once around that response.
  - Pair-level mgg energy correction disabled by default and documented as non-default systematic study only.
  - Electron momentum scale not part of photon energy smearing.

## Modular Architecture Plan

Introduce small ROOT-macro-friendly abstractions:

- `ModelParameter`
  - name, value, min, max, unit, description.
- `SmearingModel1D`
  - `name`, `parameters`, `evaluate(x)`, `print()`, `valid()`.
  - Supported model names conceptually:
    - `constant`
    - `linear`
    - `polynomial`
    - `logarithmic`
    - `a_plus_bE_plus_clnE`
- `SmearingConfig`
  - selected model names and toggles.
- `PhotonSmearingParameters`
  - energy mean parameters, energy sigma parameters, position parameters, `p_e_scale`, resolution constants.
- `PhotonSmearingModel`
  - `meanEnergy(E)`.
  - `sigmaEnergy(E_mean)`.
  - `sigmaPosition(E_mean)`.
  - `smearPhoton(E,x,y,rng)`.
- `SmearingStats`
  - processed/skipped events, low-energy clamps, nonpositive clamps, invalid parameter warnings.

Event loop, histogram filling, chi2, and plotting should call common helpers and never hard-code formula branches except through model construction.

Initial implementation can keep parameters inside macro structs. External config can come later after behavior is stable.

## Configuration And Parameter Storage

Phase 1:

- Keep macro-local config constants for minimal disruption.
- Write explicit metadata to ROOT output:
  - `smearing_model_energy_mean = a_plus_bE_plus_clnE`
  - `energy_mean_convention = reconstructed_energy_GeV`
  - `energy_log_floor_GeV`
  - `energy_valid_min_GeV`, `energy_valid_max_GeV`
  - `energy_sigma_model`
  - `position_sigma_model`
- Keep CSV but update schema handling:
  - parser reads by header names, not positions.
  - preserve column names `best_mu_a,best_mu_b,best_mu_c`.
  - optionally write compatibility aliases only if documented.

Phase 2 later:

- Move model selection and parameter bounds to simple text/JSON-like config if needed.

## Optimization Strategy

Current state:

- Coarse grid over `mu_b` and `sigma`, optional `sigma_pos`.
- Fine local grid around scalar best.
- MIGRAD refinement can fit `mu_a`, `mu_b`, `mu_c`, `sigma`, `sigma_pos`.
- `mu_a/mu_c` deterministic multistarts are enabled around the scalar seed to test whether zero coefficients are an optimizer artifact rather than a physical minimum.
- Deterministic RNG seed reduces stochastic objective noise.

Planned robust strategy:

- Objective:
  - default combined weighted chi2 over `M_gg`, `M_miss`, and `(p_target+gamma gamma)^2`.
  - Baker-Cousins enabled by default.
  - Normalize MC to data per observable before chi2; document ndf penalty.
- Empty bins:
  - keep current finite handling, but penalize data-positive/sim-zero bins instead of silent skip if feasible.
- Bounds:
  - enforce all parameter bounds in grid and MIGRAD.
  - reject nonfinite/negative sigma and nonphysical mean energies.
- Global search:
  - keep deterministic `a/c` multistarts active; expand to a broader coarse `a/c` grid or random restarts if zeros persist under ROOT validation.
  - random restarts within physical bounds after coarse grid.
  - keep current scalar grid as fast seed but document it is not enough by itself.
- Degeneracy diagnostics:
  - save chi2 values from restarts.
  - report best few minima.
  - report parameter covariance/correlation from Minuit2 if accessible.
  - plot `E_mean(E)` curves to expose equivalent shapes.
- Closure/systematics:
  - run unsmeared/known-smear closure when available.
  - vary fit observables/weights, binning, fit ranges, seeds, and section grid.

## Plotting Requirements

- Axis titles include units:
  - `E [GeV]`, `E_mean [GeV]`, `sigma_E [GeV]`, `sigma_pos [cm]`.
- Update all `mu_eff` labels to `E_mean(E)` or explicit `mu_b` where appropriate.
- Ratio/pull plots:
  - show denominator convention.
  - skip/flag empty histograms.
- Deterministic filenames:
  - `out_smear.root`
  - `out_smear_interpolated.root`
  - `section_map.csv`
  - `chi2_scans.pdf`
  - optional `smearing_validation.txt`.
- Avoid stale output:
  - pipeline should write into run-specific or cleaned output area.
  - verify modification time and existence after stages.
- Histogram names:
  - avoid duplicate temp names in same directory.
  - clone with unique names or keep `TH1::AddDirectory(kFALSE)`.
- Binning:
  - compare only matching histograms.
  - fail loudly on mismatch.

## Pipeline Consistency Plan

- `run_smearing_pipeline.sh`
  - compile `scripts/nps_sim_smearing_new_try.C`.
  - binary name `scripts/nps_sim_smearing_new_try`.
  - accept env overrides for data file, sim file, output dir, grid, bounds, `Nsmear`.
  - validate inputs and output directory.
  - remove or archive known stale outputs before run, or use timestamped run dir.
  - verify `out_smear.root`, `section_map.csv`, `chi2_scans.pdf`, `out_smear_interpolated.root`.
  - run downstream with explicit env paths.
  - verify `simc_pi0_analysis_output_smeared.root`.
- `simc_pi0_analysis.C`
  - consume header-based CSV.
  - load `h_mu_a_interp`, `h_mu_interp`/`h_mu_b_interp`, `h_mu_c_interp`.
  - apply `E_mean=a+bE+c ln(E_safe)` once.
  - align position-energy dependence flags with fitter metadata or explicit config.
  - fail loudly when `NPS_SMEARING_MODE=section` and required smearing files missing.

## Bias And Consistency Checks

- Compare data, unsmeared MC, smeared MC for:
  - `M_gg` peak/width.
  - `M_miss` peak/width.
  - `(p_target+gamma gamma)^2`.
- Check energy residuals:
  - `E_mean(E)-E` versus E.
  - low-energy clamp count.
- Check position residuals:
  - `sigma_pos` maps and optional energy dependence.
- Check kinematic dependence:
  - `Q2`, `W`, `xB`, `z`, `t`, `pt` if branches available downstream.
- Check acceptance/selection:
  - document cuts applied before smearing in `simc_pi0_analysis.C`.
  - verify smearing does not silently change which events are written unless intended.
- Check normalization:
  - MC normalized to data for fit plots and chi2.
  - Tree weights preserved unless documented.
- Check exclusivity weighting:
  - data and simulation exclusivity branch names are configured independently.
  - missing configured exclusivity branch is a hard error when exclusivity weighting is enabled.
  - any intentional data/sim branch-name difference must be documented.
- Document limitations:
  - parameter degeneracy between energy response, resolution, position resolution, and electron momentum scale.
  - section-wise fits with low stats may be unstable.

## Documentation Sync Rule

- `smearing.md` is now the canonical smearing documentation.
- `smearing_documentation.md` is a redirect pointer to `smearing.md`.
- Keep all relevant `.md` files updated together as code, defaults, validation, file names, parameter conventions, or known limitations change.
- At minimum, update `plan_smearing.md` and `smearing.md` whenever smearing algorithm behavior changes.

## Optimization Degeneracy Notes

- With energy mean, energy resolution, and position smearing enabled, per-section fit is at least 5D:
  - `mu_a`, `mu_b`, `mu_c`, `sigma`, `sigma_pos`.
- If electron momentum scaling is fitted too, fit becomes 6D:
  - `mu_a`, `mu_b`, `mu_c`, `sigma`, `sigma_pos`, `p_e_scale`.
- Current algorithm is now global-seeded hybrid by default:
  - deterministic ROOT/GSL Sobol seed scan across active bounded parameters.
  - keep best N seeds.
  - Minuit2 MIGRAD from each retained seed.
  - HESSE strongest-correlation diagnostics.
  - limited profile scans for strongly correlated parameters.
  - deterministic RNG/common random numbers.
- Missing robust global diagnostics:
  - true Sobol generator.
  - full HESSE covariance matrix output.
  - broad profile scans where one parameter is fixed and all others re-minimized.
  - closure tests with injected known smearing.
- Recommended future strategy:
  - bounded Sobol/Latin-hypercube seeds.
  - keep best N seeds.
  - run MIGRAD from each seed.
  - run HESSE at the selected minima.
  - profile suspicious parameters.
  - compare staged and simultaneous modes.

## Concrete Task Checklist

- [x] Inspect target files and current pipeline.
- [x] Create this living plan before code edits.
- [x] Add modular smearing model structs/helpers to fitter.
- [x] Replace direct `mu_a,mu_b,mu_c` formula use inside `smearPhoton` with model evaluation helpers.
- [x] Keep only `a_plus_bE_plus_clnE` as active default energy mean model.
- [x] Disable/document pair-level mgg energy correction as non-default only.
- [x] Improve fitter CSV/metadata output.
- [x] Update downstream CSV parsing by header name.
- [x] Update downstream interpolated map loading for `mu_a`, `mu_b`, `mu_c`.
- [x] Update downstream smearing application to use `E_mean=a+bE+c ln(E_safe)`.
- [x] Update downstream metadata strings and runtime printouts.
- [x] Update pipeline shell script to compile/run `nps_sim_smearing_new_try.C`.
- [x] Add stale-output protection and output verification.
- [x] Add smearing documentation in `smearing_documentation.md` or `docs/smearing.md`.
- [x] Add deterministic nonzero `mu_a/mu_c` optimizer starts to test zero-coefficient artifact risk.
- [x] Align `scripts/compare_mpi0_mmiss_mpgg.ipynb` defaults with smearing pipeline weighting/exclusivity conventions.
- [x] Create canonical `smearing.md` with current optimizer/degeneracy discussion.
- [x] Convert `smearing_documentation.md` to canonical-doc redirect.
- [x] Split fitter exclusivity config into `DATA_EXCLUSIVITY_BRANCH` and `SIM_EXCLUSIVITY_BRANCH`.
- [x] Add global low-discrepancy multistart default with best-N MIGRAD/HESSE diagnostics.
- [x] Add optimizer summary, seed, and profile CSV outputs.
- [x] Keep position smearing and electron momentum scale disabled by default for current work.
- [x] Remove active Pass-1/Pass-2 orchestration.
- [x] Add iterative coupled sweep workflow with nominal external response in sweep 1 and previous-sweep external response in later sweeps.
- [x] Remove final serial refinement sweep; final state is produced by barriered parallel coupled sweeps only.
- [x] Tighten default `mu_a/mu_c` bounds; HESSE/profile diagnostics are active again in the Sobol workflow.
- [x] Keep finite lower-chi2 Minuit coordinates even when Minuit reports a non-success status.
- [x] Replace Halton fallback with ROOT/GSL Sobol multistart.
- [x] Restore active HESSE/profile diagnostics in the intended Sobol workflow.
- [x] Add final-state closure consistency CSV output.
- [x] Run compile validation.
- [ ] Run smallest available pipeline validation if inputs and ROOT are available.
- [x] Update this plan with validation results and open issues.

## Validation Log

- `bash -n run_smearing_pipeline.sh`: passed.
- `git diff --check` on changed deliverables: passed after whitespace cleanup.
- `python3 -m json.tool scripts/compare_mpi0_mmiss_mpgg.ipynb`: passed after notebook edits and output clearing.
- `./run_smearing_pipeline.sh` in the bare shell stopped immediately with `[ERROR] root not found in PATH`, before archiving/touching outputs. This verifies fail-loud behavior for missing ROOT.
- ROOT is available after sourcing `/group/nps/singhav/setup.csh` through `csh`; use that setup before ROOT validation.
- Existing `output/plots/x60_4b/production_wfpi0/section_map.csv` was inspected and confirmed to use current `best_mu_a,best_mu_b,best_mu_c` schema.
- Static optimizer audit found an artificial zero-risk for `mu_a` and `mu_c`: before this update, only a single zero-seed MIGRAD start could move them after scalar `mu_b` grid scans. The new multistart printouts report number of starts and whether a nonzero `a/c` seed supplied the selected result.
- `smearing.md` created as canonical documentation and updated with current 5D/6D degeneracy risks, current optimizer status, and recommended future global-search strategy.
- `scripts/nps_sim_smearing_new_try.C` now lets data and simulation use separate exclusivity branch names while preserving hard-fail behavior when enabled branches are missing.
- `scripts/nps_sim_smearing_new_try.C` now defaults to low-discrepancy global multistart with best-N MIGRAD refinement, HESSE strongest-correlation diagnostics, and optional limited profile scans.
- New optimizer debug CSVs are written beside `section_map.csv`: `smearing_optimizer_summary.csv`, `smearing_optimizer_seeds.csv`, `smearing_optimizer_profiles.csv`.
- Static validation after global-multistart refactor:
  - `git diff --check -- scripts/nps_sim_smearing_new_try.C smearing.md plan_smearing.md smearing_documentation.md`: passed.
  - `bash -n run_smearing_pipeline.sh`: passed.
  - bare-shell `./run_smearing_pipeline.sh`: fails early with `[ERROR] root not found in PATH`; expected unless ROOT setup is sourced.
- Current optimizer correction after poor fit / long runtime report:
  - Default orchestration is now iterative coupled sweeps with no active Pass-1/Pass-2 prefit path.
  - Sweep 1 uses nominal external photon response; later sweeps use previous completed sweep results.
  - Default `mu_a` and `mu_c` bounds were tightened from `[-1,1]` GeV to `[-0.20,0.20]` GeV. The former range allowed very large low-energy nonlinearity and could drive bad minima.
  - HESSE/profile scans were temporarily disabled during debugging; they are active again in the Sobol workflow.
  - Minuit lower-chi2 coordinates are now retained even if Minuit reports a non-success status. This prevents losing improved points because of covariance/convergence-status failures.
  - Compile check passed after this change using the sourced ROOT setup:
    `csh -c 'source /group/nps/singhav/setup.csh; g++ scripts/nps_sim_smearing_new_try.C \`root-config --cflags --libs\` -O2 -std=c++17 -fopenmp -I./src -o /tmp/nps_sim_smearing_new_try_check'`.
  - `git diff --check -- scripts/nps_sim_smearing_new_try.C smearing.md plan_smearing.md`: passed.
- SIMC output metadata write fix:
  - ROOT errors like `TROOT::WriteTObject: current directory (Rint) is not associated with a file` came from writing smearing metadata/maps after 2D-mass-cut plotting changed the current ROOT directory.
  - `scripts/simc_pi0_analysis.C` now calls `outfile->cd()` and `outfile_smeared->cd()` again immediately after adding/plotting 2D mass cut branches and before writing `TNamed` metadata or smearing maps.
  - Compile check passed with sourced ROOT setup:
    `csh -c 'source /group/nps/singhav/setup.csh; g++ -c scripts/simc_pi0_analysis.C \`root-config --cflags\` -O2 -std=c++17 -I./src -o /tmp/simc_pi0_analysis_check.o'`.
- Final serial sweep removal:
  - Removed `ENABLE_FINAL_SERIAL_REFINEMENT_SWEEP` and the active `[final-serial]` fit block from `scripts/nps_sim_smearing_new_try.C`.
  - The macro now goes directly from completed coupled sweeps to final fit-status labeling, optional global `p_e_scale` Stage 4, and final chi2 recomputation.
  - Validation after removal passed: `bash -n run_smearing_pipeline.sh`, `git diff --check -- scripts/nps_sim_smearing_new_try.C plan_smearing.md smearing.md`, and ROOT compile via the repository `csh` setup command with `root-config`, `-lMathMore`, and output `/tmp/nps_sim_smearing_new_try_check`.
- Iterative coupled sweep workflow update:
  - Active code no longer runs Pass-1 or Pass-2 blocks.
  - `ENABLE_PARALLEL_SECTION_FITS=true` by default. Fits inside one sweep may run in parallel; next sweep waits for full sweep completion.
  - External photon parameters are reset to nominal before sweep 1.
  - After each sweep, cross-boundary external photon parameters refresh from completed section fits.
  - Final serial refinement sweep has been removed to avoid section-order bias.
  - The final accepted state is the last completed barriered sweep, followed by final status labeling and chi2 recomputation.
  - Compile/link check passed with direct ROOT library path and `-lMathMore`.
  - `bash -n run_smearing_pipeline.sh`: passed.
  - `git diff --check -- scripts/nps_sim_smearing_new_try.C run_smearing_pipeline.sh smearing.md plan_smearing.md`: passed.
- Sobol workflow update:
  - Halton seed generation was removed from `scripts/nps_sim_smearing_new_try.C`.
  - Multistart now uses `ROOT::Math::QuasiRandomSobol` over the active fit dimensions.
  - Active default workflow is Sobol bounded seeds -> keep best N -> MIGRAD each retained seed -> HESSE diagnostics -> profile scans for high-correlation suspicious parameters.
  - `run_smearing_pipeline.sh` now links `-lMathMore`, required for ROOT/GSL Sobol.
  - New `smearing_closure_summary.csv` records before/after final-state chi2 consistency.
  - Compile-only check passed against ROOT headers.
  - Full direct link check passed with explicit ROOT lib path and `-lMathMore`; normal `setup.csh` validation could not be repeated because the module system failed to resolve user id `11606` in this shell.

## Design Decisions

- Use additive reconstructed-energy mean convention for `a+bE+c ln(E)` because current fitter code already computes `E_sc` directly from that expression.
- Keep CSV for now, but parse by header names to avoid future schema breaks.
- Keep ROOT macro style. Avoid introducing external dependencies.
- Preserve `h_mu_interp` as compatibility alias for `mu_b`, and add explicit `h_mu_b_interp`.
- Downstream `simc_pi0_analysis.C` accepts env overrides for smearing files and output files so the pipeline can pass exact artifacts.
- Pipeline archives previous canonical outputs to `smearing_pipeline_backup/<RUN_TAG>/` instead of deleting them.
- Notebook comparison defaults now match fitter defaults: unit SIMC weights unless `USE_SIM_FULL_WEIGHT` is explicitly enabled, exclusivity applied via `is_exclusive`, and no extra cluster-energy cut beyond producer cuts.
- Canonical documentation file is `smearing.md`; `smearing_documentation.md` remains only as redirect to prevent stale duplicate content.
- Exclusivity branch config is intentionally separate for data and simulation so combined data trees and SIMC trees can use different branch naming conventions.
- Sweep 1 uses nominal cross-boundary photon response; later sweeps use the previous completed sweep. The final accepted state is therefore order-independent within each sweep and avoids a last-section advantage from serial Gauss-Seidel updates.

## Open Issues

- Need run compile and ROOT validation in an environment where `root` and `root-config` are in `PATH`.
- Need confirm production input files once ROOT environment is loaded.
- Downstream now fails when `NPS_SMEARING_MODE=section` is explicitly requested and required smearing files are missing; interactive auto mode still produces unsmeared output if smearing files are absent.
- Fitter still has legacy comments/compatibility around scalar `mu_b`, but default and downstream active energy response are `a_plus_bE_plus_clnE`.
- Need ROOT validation to determine whether `mu_a` and `mu_c` remain zero after multistart. If they do, next diagnostic is broader `a/c` chi2 surface plots and correlation/covariance output.
- Full covariance matrix output remains open before claiming complete 5D/6D degeneracy diagnostics.
- Sobol generator is now implemented using `ROOT::Math::QuasiRandomSobol`; compile-only validation passed, and full link validation passed with direct ROOT lib path plus `-lMathMore`.
- Full closure-test automation remains open.
