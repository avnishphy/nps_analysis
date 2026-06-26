# NPS Pi0 Smearing Workflow

This is the only canonical markdown document for the NPS pi0 smearing workflow
in this repository. Keep this file synchronized with:

- `scripts/nps_sim_smearing_new_try.C`
- `run_smearing_pipeline.sh`
- `scripts/simc_pi0_analysis.C`
- `scripts/compare_smearing_outputs.py`
- `scripts/summarize_smearing_run.py`

The goal of the smearing workflow is to tune detector-response parameters for
SIMC photons so that reconstructed SIMC shapes match data in the same exclusive
event selection. The tuned response is stored per calorimeter section and then
applied downstream to produce a smeared SIMC ROOT tree.

## Quick Start

Always load the Hall C/NPS ROOT environment before compiling or running:

```bash
csh -c 'source /group/nps/singhav/setup.csh; ./run_smearing_pipeline.sh'
```

The wrapper:

1. Compiles `scripts/nps_sim_smearing_new_try.C`.
2. Runs the section-wise smearing fit.
3. Writes timestamped fit artifacts and latest compatibility copies.
4. Runs `scripts/simc_pi0_analysis.C` to regenerate smeared SIMC.
5. Optionally writes reproducible data/SIMC/smeared-SIMC comparison metrics.
6. Writes compact runtime summaries and best-parameter JSON.

Useful reduced smoke test:

```bash
csh -c 'source /group/nps/singhav/setup.csh; env OMP_NUM_THREADS=8 RUN_COMPARISON=0 NX=1 NY=1 NSMEAR=1 OUT_DIR=/tmp/nps_smearing_smoke OUT_FILE=/tmp/nps_smearing_smoke/out_smear.root RUN_TAG=smoke bash ./run_smearing_pipeline.sh'
```

## Main Files

- `run_smearing_pipeline.sh`: production wrapper for compile, fit, downstream
  smeared SIMC regeneration, comparisons, and summaries.
- `scripts/nps_sim_smearing_new_try.C`: active fitter.
- `scripts/simc_pi0_analysis.C`: downstream SIMC producer and smearing consumer.
- `scripts/compare_smearing_outputs.py`: reproducible before/after comparison.
- `scripts/summarize_smearing_run.py`: compact CSV/JSON summary writer.
- `scripts/compare_smearing_strategies.py`: compares multiple smearing run
  directories.

Do not use old staged smearing notes as documentation. The active algorithm is
the coupled Sobol/MIGRAD section workflow described here.

## Inputs

The fitter takes:

```text
data.root dataTree sim.root simTree out.root nx ny x_min x_max y_min y_max overlap_frac [Nsmear]
```

The wrapper provides defaults for the current production workflow:

```bash
DATA_FILE=/.../combined_branches_LH2_wfpi0.root
SIM_FILE=/.../simc_pi0_analysis_output.root
OUT_DIR=output/plots/x60_4b/production_wfpi0
OUT_FILE=$OUT_DIR/out_smear.root
NX=7
NY=7
X_MIN=-24
X_MAX=28
Y_MIN=-34
Y_MAX=34
OVERLAP_FRAC=0.0
NSMEAR=80
NPS_SMEARING_SWEEP_ACCEPTANCE=jacobi_global_accept_rollback
```

All of these can be overridden in the shell environment.

## Required Branches

### Data

The fitter reads data branches including:

- `cluster_x_1`, `cluster_y_1`, `cluster_e_1`
- `cluster_x_2`, `cluster_y_2`, `cluster_e_2`
- `pi0_weight`
- `scale`
- `mmiss_all`
- `is_exclusive_ellipse_combined` by default

The data event weight is:

```text
w_data = pi0_weight * scale * is_exclusive_ellipse_combined
```

When exclusivity selection is enabled, the configured exclusivity branch must
exist. Missing branches are treated as errors, not warnings.

### SIMC

The fitter reads SIMC branches including:

- `clust_X`
- `clust_Y`
- `clust_E`
- `nclust`
- `full_weight`
- `siglab`
- `is_exclusive_ellipse` by default
- `sc_e_E`, `sc_e_Px`, `sc_e_Py`, `sc_e_Pz` for missing mass

SIMC smearing weights are always de-modeled from the SIMC cross-section model:

```text
w_sim_base = full_weight / siglab
w_sim = w_sim_base * is_exclusive_ellipse
```

This is intentionally consistent with
`scripts/excl_xsec_pi0_analysis_no_simc_model.C`. The point is to fit detector
response on the accepted phase-space basis, not on the generator model
cross-section shape. The branch choice is documented in `demodeling.md`.

The fitter skips SIMC events with invalid model weights:

```text
!finite(full_weight)
!finite(siglab)
abs(siglab) < 1e-20
!finite(full_weight/siglab)
```

The skipped count is printed to the log and stored in the ROOT output metadata.
If `full_weight` or `siglab` is missing, the fitter stops with a clear error.

## Exclusivity Convention

Current default:

```text
APPLY_IS_EXCLUSIVE_SELECTION = true
DATA_EXCLUSIVITY_BRANCH = "is_exclusive_ellipse_combined"
SIM_EXCLUSIVITY_BRANCH  = "is_exclusive_ellipse"
```

This is deliberate. The combined data tree uses the combined ellipse branch,
while the SIMC tree uses the SIMC ellipse branch written by
`scripts/simc_pi0_analysis.C`.

If the ellipse logic changes, regenerate both:

1. the combined data tree, and
2. the unsmeared SIMC input tree.

Otherwise the fitter can compare samples made with inconsistent exclusive cuts.

## Event Geometry And Section Assignment

The calorimeter is divided into an `nx` by `ny` base grid. Each section has:

```text
xlo, xhi, ylo, yhi
```

`OVERLAP_FRAC` can enlarge sections around their base-cell centers. Current
production default is no overlap.

Data and SIMC section buffers are filled only for events where both photons are
inside the configured calorimeter geometry. A section receives an event if
either photon lies inside that section.

For each section-specific SIMC copy, the fitter records:

```text
photon1_in_section
photon2_in_section
```

This matters because only photons owned by the fitted section receive the trial
section parameters during the section objective. The other photon uses the
external response supplied by the current global/coupled-sweep state.

The all-section summary buffer is unique-event based, not section duplicated.
It requires both photons in geometry. This prevents section overlap from
double-counting global diagnostic histograms.

## Smearing Model

### Energy Mean

The active energy response model is:

```text
E_safe = max(E, E_floor)
E_mean(E) = a + b*E_safe + c*ln(E_safe)
```

with default:

```text
E_floor = 0.2 GeV
```

Important convention:

- `E_mean` is already an energy in GeV.
- It is not a multiplicative scale factor.
- Do not compute `E * E_mean`.

Parameter units:

- `a`: GeV
- `b`: dimensionless
- `c`: GeV

The disabled scalar form is:

```text
E_mean(E) = b*E_safe
```

but the production workflow uses the full `a + bE + c ln(E)` model.

### Energy Resolution

Default simple stochastic model:

```text
sigma_E(E_mean) = sigma * sqrt(E_mean)
```

Optional three-term model:

```text
sigma_E = sigma * E_mean * sqrt((A/sqrt(E_mean))^2 + (B/E_mean)^2 + C^2)
```

The default simple model treats `sigma` as a data-driven stochastic coefficient.
The three-term model is available only when detector resolution constants are
trusted.

### Energy Smearing PDF

The configured smearing shape can be Gaussian or Landau. Gaussian mode treats
`sigma_E` as a standard deviation. Landau mode treats the fitted `sigma`
through an FWHM convention and has redraw/tail controls to avoid pathological
nonphysical tails.

Nonpositive sampled energies are protected by a small positive clamp.

### Position Smearing

Position smearing is additive:

```text
x_sm = x + Gaus(0, sigma_pos)
y_sm = y + Gaus(0, sigma_pos)
```

Optional energy dependence:

```text
sigma_pos_eff = sigma_pos * sqrt(E0 / E_mean)
```

Position smearing primarily affects opening angle and therefore `M_gg`.

### Optional Pair-Energy Correction

There is an optional linear pair-energy correction versus reconstructed
`M_gg`. It is disabled by default. Do not enable it casually because it can
absorb physics-model differences into an empirical correction.

### Optional Electron Momentum Scale

`p_e_scale` is available only when electron momentum scaling is enabled and the
objective uses `M_miss`. It is global in the final refinement mode, not a
photon smearing parameter.

Default is disabled.

## Objective Histograms

The fitter can compare:

- `M_gg`
- `M_miss`
- `(p_target + gamma gamma)^2`, called `Mpgg2` in many diagnostics

Current default objective mode uses `M_gg` and `M_miss` with weights:

```text
W_MPI0 = 1.0
W_MMISS = 1.0
W_MPGG2 = 0.0
```

`Mpgg2` is still filled and diagnosed even when its objective weight is zero.

Histogram windows are intentionally focused on useful regions:

```text
M_gg:   0.11 to 0.15 GeV/c^2, 50 bins
M_miss: 0.60 to 1.30 GeV/c^2, 80 bins
Mpgg2:  4.0 to 14.0 GeV^2, 50 bins
```

## Normalization Trick

The diagnostics normalize unsmeared simulation to the data integral first, then
apply the same scale factor to the smeared histogram.

That means:

```text
scale = Integral(data) / Integral(unsmeared_sim)
unsmeared_sim *= scale
smeared_sim   *= scale
```

This is not the same as normalizing the smeared histogram independently after
smearing. The current convention keeps smearing-induced migration into or out
of the plotted mass window visible.

Inside the objective evaluation, simulation templates are area-normalized to
the corresponding data histogram before the chi2 statistic is computed. The
fit is therefore a shape/response calibration, not a yield extraction.

## Chi2 Statistic

Baker-Cousins chi2 is enabled by default. It uses a finite floor for bins where
data are positive and simulation is empty, rather than silently skipping those
bins.

This is important because skipping data-positive/simulation-empty bins can
produce deceptively tiny chi2 values while the plot visibly disagrees.

The scalar objective is:

```text
chi2_total =
  W_MPI0  * chi2(M_gg)
+ W_MMISS * chi2(M_miss)
+ W_MPGG2 * chi2(Mpgg2)
```

Always inspect the objective breakdown, pulls, residuals, and mass plots. A
single global chi2/ndf is not enough to sign off the calibration.

## Deterministic Randomness

The smearing fit uses common random numbers:

```text
SMEAR_DETERMINISTIC_SEED = 42
```

Each parameter evaluation sees the same sequence of random draws. This makes
the chi2 surface smoother and makes repeated runs reproducible when the input
configuration is unchanged.

The per-section final output histograms use deterministic section-dependent
seeds. The point is reproducibility, not physical event-level randomness.

## Optimizer Acceleration

The optimizer does not always use every SIMC event for every trial point.
Instead, it builds deterministic, `M_gg`-stratified, weight-compensated subsets:

```text
OPT_MAX_SIM_EVENTS_PER_SECTION = 25000
OPT_MAX_SIM_EVENTS_GLOBAL_PREFIT = 60000
OPT_SUBSET_MGG_BINS = 32
```

The subset preserves weighted bin sums as well as possible while reducing
runtime. Final histograms, final chi2 recomputation, section maps, and
all-section summaries use the full SIMC buffers and command-line `Nsmear`.

The optimizer-only `Nsmear` is capped:

```text
optimizer_nsmear = min(command_line_Nsmear, OPTIMIZATION_NSMEAR)
```

Current default:

```text
OPTIMIZATION_NSMEAR = 80
```

## Optimization Workflow

The current active workflow is:

```text
global all-calorimeter prefit
-> coupled section sweeps
-> global candidate acceptance / rollback
-> optional final global p_e_scale refinement
-> final full-stat recomputation and diagnostics
```

### Global Prefit

Before section sweeps, the fitter optimizes one shared response over all
calorimeter events. This global prefit uses the same active response model and
the same Sobol/MIGRAD optimizer as the section fits.

The global prefit provides:

- the starting seed for sweep 1, and
- the out-of-section photon response for sweep 1.

If the global prefit is disabled or fails due to low statistics, sweep 1 falls
back to nominal response:

```text
a = 0
b = 1
c = 0
sigma = 0
sigma_pos = 0
p_e_scale = 1
```

### Section Fits

Each active section fit uses:

```text
Sobol bounded seeds
-> keep best GLOBAL_MULTISTART_KEEP_BEST seeds
-> Minuit2 MIGRAD from each kept seed
-> HESSE diagnostics
-> limited profile scans for strong correlations
```

The active fit parameters are:

```text
a, b, c, sigma, sigma_pos
```

If electron momentum scaling is enabled, `p_e_scale` can also participate.

### Sigma-Position Staging

To reduce over-sharpening and parameter degeneracy, section fits are staged:

1. Hold `sigma_pos` fixed at the previous/global value.
2. Optimize energy response and energy resolution.
3. Hold those energy parameters fixed.
4. Optimize `sigma_pos`.

This makes the optimizer less able to compensate energy-response mistakes by
immediately squeezing or broadening the angular response.

### Coupled Sweeps

Each sweep fits every active section. The sweep is barriered: all sections in
sweep N finish before sweep N+1 can begin.

The fitter keeps three separate states:

```text
candidate_parameters
current_accepted_parameters
best_global_parameters
```

After a sweep finishes, the candidate state is evaluated with a unique
all-calorimeter objective. If the candidate improves the accepted objective
within configured tolerances, it becomes the accepted state. If it does not,
the working state rolls back to the previous accepted state.

This rollback matters. Without it, a globally worse sweep can seed later
sweeps and make the workflow drift toward a narrow or over-optimized solution.

The fitter also detects repeated rejected candidates/cycles and stops early
instead of spending more deterministic sweeps on the same rejected solution.

There is no final serial Gauss-Seidel sweep. The final state is the best
completed barriered sweep according to the all-section objective, followed by
final status labeling and full-stat recomputation.

## Cache Replay

The fitter computes a configuration/source/input fingerprint. If the
fingerprint is unchanged and the previous optimizer CSVs exist, it can replay
the global prefit and sweep diagnostic plots without rerunning full Sobol/MIGRAD
optimization.

This is useful because Sobol sampling is deterministic and the cache can
rebuild the plots quickly.

Set:

```bash
NPS_SMEARING_FORCE_REOPT=1
```

to force a full optimization rerun.

The cache fingerprint includes:

- source hash
- input data/sim file signatures
- tree names
- grid and geometry
- overlap
- `Nsmear`

If any of those change, the optimization reruns.

## Timestamped Artifacts

`RUN_TAG` controls timestamped output names. If `RUN_TAG` is not set, the
wrapper generates one like:

```text
YYYYMMDD_HHMMSS
```

The production artifacts are timestamped:

```text
out_smear_<RUN_TAG>.root
out_smear_<RUN_TAG>_interpolated.root
chi2_scans_<RUN_TAG>.pdf
chi2_scans_progress_<RUN_TAG>/
smearing_artifacts_<RUN_TAG>.json
```

Latest compatibility copies are also maintained:

```text
out_smear.root
out_smear_interpolated.root
chi2_scans.pdf
```

Downstream tools can keep using the latest names, while provenance-sensitive
analysis can use the timestamped artifacts.

The wrapper archives previous latest outputs into:

```text
smearing_pipeline_backup/<RUN_TAG>/
```

before refreshing them.

## Main Outputs

### ROOT Outputs

`out_smear_<RUN_TAG>.root` contains:

- top-level metadata with model, optimizer, timestamp, cache, and de-modeling
  information
- per-section directories
- data histograms
- unsmeared histograms
- final smeared histograms
- `fit_params` text objects
- `diagnostic_canvases/`
- `diagnostic_maps/`

Important top-level metadata includes:

```text
smearing_model_energy_mean = a_plus_bE_plus_clnE
energy_mean_convention = reconstructed_energy_GeV
energy_sigma_model
position_sigma_model
optimizer_model
run_tag
timestamped_output_file
chi2_pdf_file
chi2_progress_dir
sim_weight_mode = full_weight_over_siglab
sim_model_xsec_branch = siglab
sim_skipped_invalid_model_weight_events
```

`out_smear_<RUN_TAG>_interpolated.root` contains:

- `h_mu_a_interp`
- `h_mu_interp`
- `h_mu_b_interp`
- `h_mu_c_interp`
- `h_sigma_interp`
- `h_sigma_pos_interp`
- response-ratio maps at fixed energies:
  - `h_response_ratio_interp_E1GeV`
  - `h_response_ratio_interp_E2GeV`
  - `h_response_ratio_interp_E3GeV`
  - `h_response_ratio_interp_E4GeV`
  - `h_response_ratio_interp_E5GeV`
- `interpolated_canvases/`

### CSV And JSON Outputs

The fit writes:

- `section_map.csv`
- `smearing_optimizer_summary.csv`
- `smearing_optimizer_seeds.csv`
- `smearing_optimizer_profiles.csv`
- `smearing_sweep_history.csv`
- `smearing_objective_breakdown.csv`
- `smearing_closure_summary.csv`
- `smearing_runtime_summary.csv`
- `best_parameters.json`
- `optimization_config.json`
- `smearing_artifacts_<RUN_TAG>.json`

`section_map.csv` is the primary tabulated parameter output. It includes:

```text
ix, iy, xlo, xhi, ylo, yhi,
n_data, n_sim,
best_mu_a, best_mu_b, best_mu_c,
best_sigma, best_sigma_pos_cm, best_p_e_scale,
res_A, res_B, res_C,
best_chi2, ndf, chi2_ndf, fit_status
```

`smearing_sweep_history.csv` is the main file for understanding how the sweep
evolved. It records candidate acceptance, previous/current/best objective,
objective components, repeated-candidate status, parameter norms, runtime, and
section parameters.

`smearing_objective_breakdown.csv` records objective contributions for each
observable. Use it to see which mass distribution is driving the fit.

`smearing_closure_summary.csv` compares optimizer-subset chi2 to final
full-event chi2.

## Diagnostic PDFs

`chi2_scans_<RUN_TAG>.pdf` is the main diagnostic PDF. `chi2_scans.pdf` is the
latest copy.

Progress pages are also written while the optimizer proceeds:

```text
chi2_scans_progress_<RUN_TAG>/
```

This lets users inspect global prefit and sweep evolution before the full job
finishes.

The PDF includes:

- global prefit mass plots
- sweep candidate mass plots
- section-wise mass plots for sweeps
- final per-section comparison pages
- pull and fractional-residual pages
- parameter-location-in-bounds diagnostics
- seed/MIGRAD diagnostics
- profile-scan diagnostics
- all-section mass comparisons
- mismatch attribution maps
- parameter maps
- response-ratio curves and maps
- final text summary

## Response-Ratio Diagnostics

Energy response plots use:

```text
R(E) = E_mean(E) / E
```

The curve x-axis is restricted to:

```text
0.5 GeV <= E <= 6.0 GeV
```

This range is deliberate. The lower tail below 0.5 GeV is not useful for the
current diagnostic purpose and can visually dominate the nonlinear response.

Fixed-energy 2D response maps show `E_mean/E` and explicitly include the
associated energy in the title. This is easier to interpret than absolute
`E_mean`, because values near 1.0 indicate no response correction at that
energy.

## Downstream Smearing Consumer

The pipeline runs:

```bash
NPS_SMEARING_MODE=section \
NPS_SMEAR_FILE="$OUT_FILE" \
NPS_SMEAR_INTERP_FILE="$INTERP_FILE" \
NPS_SECTION_MAP_FILE="$SECTION_MAP" \
NPS_SIMC_OUTPUT_FILE="$UNSMEARED_SIMC_OUT" \
NPS_SIMC_SMEARED_OUTPUT_FILE="$SMEARED_SIMC_OUT" \
root -l -b -q scripts/simc_pi0_analysis.C
```

`scripts/simc_pi0_analysis.C` consumes the section map and interpolated maps to
write:

```text
simc_pi0_analysis_output_smeared.root
```

The downstream macro must use the same energy-mean convention and position
smearing convention as the fitter. If the smearing model changes in the fitter,
verify the consumer before trusting the smeared tree.

## Comparison Metrics

When `RUN_COMPARISON=1`, the wrapper runs:

```bash
python3 scripts/compare_smearing_outputs.py
```

It writes:

```text
smearing_comparison/smearing_comparison_metrics.csv
smearing_comparison/smearing_comparison.png
```

For fast tests:

```bash
COMPARISON_MAX_ENTRIES=50000
```

The comparison script can use PyROOT or uproot. For large Hall C files, prefer
the sourced ROOT environment and PyROOT backend.

## How To Read Results

First inspect:

1. `chi2_scans_<RUN_TAG>.pdf`
2. `smearing_sweep_history.csv`
3. `smearing_objective_breakdown.csv`
4. `smearing_closure_summary.csv`
5. `smearing_comparison/smearing_comparison_metrics.csv`
6. `out_smear_<RUN_TAG>.root` canvases and maps

Questions to ask:

- Did the accepted sweep actually improve the global objective?
- Did a later candidate get rejected and rolled back?
- Are parameters stuck on bounds?
- Is one observable dominating the objective?
- Do final full-event chi2 values agree qualitatively with optimizer-subset
  chi2?
- Are response ratios physically reasonable across sections?
- Are `sigma_pos` values plausible, or is position smearing compensating for
  energy-response problems?
- Did many SIMC events get skipped due to invalid `full_weight/siglab`?
- Are data and SIMC using the intended ellipse branches?

## Common Failure Modes

### Missing ROOT Environment

Symptom:

```text
[ERROR] root not found in PATH
```

Fix:

```bash
csh -c 'source /group/nps/singhav/setup.csh; ./run_smearing_pipeline.sh'
```

### Missing `siglab`

The fitter now requires `siglab` for SIMC de-modeling. If the branch is missing,
regenerate the SIMC input with the model cross-section branch preserved.

### Missing Ellipse Branch

If exclusivity selection is enabled, the data and SIMC exclusivity branches
must exist. Regenerate combined data or unsmeared SIMC after changing ellipse
cut logic.

### Smeared Peak Too High Or Too Narrow

Check:

- whether unsmeared normalization happens before smeared comparison
- `sigma_pos` staging behavior
- response-ratio maps
- `smearing_objective_breakdown.csv`
- whether `M_gg` is over-weighted relative to `M_miss`
- whether the accepted sweep is older than a rejected later sweep

### Low Or Misleading Chi2

Check for:

- empty-simulation/data-positive bins
- excessive area normalization effects
- too few events in a section
- optimizer subset not matching full-event recompute
- pathological `full_weight/siglab` tails

### Parameters On Bounds

Bound-hugging usually means:

- the model is too restrictive,
- the section has poor statistics,
- a different observable is pulling the fit,
- the SIMC/data selection differs, or
- the optimizer found a compensation direction rather than a physical response.

## Validation Commands

Compile the fitter:

```bash
csh -c 'source /group/nps/singhav/setup.csh; g++ scripts/nps_sim_smearing_new_try.C `root-config --cflags --libs` -lMathMore -O2 -std=c++17 -fopenmp -I../src -o /tmp/nps_sim_smearing_new_try_check'
```

Check wrapper syntax:

```bash
bash -n run_smearing_pipeline.sh
```

Check documentation/code whitespace:

```bash
git diff --check -- scripts/nps_sim_smearing_new_try.C run_smearing_pipeline.sh smearing.md
```

Run reduced smoke test:

```bash
csh -c 'source /group/nps/singhav/setup.csh; env OMP_NUM_THREADS=8 RUN_COMPARISON=0 NX=1 NY=1 NSMEAR=1 OUT_DIR=/tmp/nps_smearing_smoke OUT_FILE=/tmp/nps_smearing_smoke/out_smear.root RUN_TAG=smoke bash ./run_smearing_pipeline.sh'
```

Run production-style pipeline:

```bash
csh -c 'source /group/nps/singhav/setup.csh; ./run_smearing_pipeline.sh'
```

## Maintenance Rules

- Keep this as the only smearing markdown document.
- Do not reintroduce duplicate smearing plan files.
- When changing the fitter, update this file in the same work pass.
- When changing the ellipse cut, regenerate both data and SIMC inputs.
- When changing the SIMC weighting convention, update the de-modeling section.
- When changing output names, update the timestamped artifact section.
- When changing optimizer state logic, update the coupled-sweep section.
- When changing downstream smearing application, update the consumer section.
