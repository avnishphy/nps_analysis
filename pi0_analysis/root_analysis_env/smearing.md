# NPS Pi0 Smearing Framework

## Purpose

Smearing adjusts SIMC photon energies and positions so reconstructed simulation can be compared to data with detector response effects included. Fitted parameters are stored per calorimeter section and consumed by `scripts/simc_pi0_analysis.C` to produce `simc_pi0_analysis_output_smeared.root`.

This file is canonical smearing documentation. Keep `plan_smearing.md`, `smearing.md`, and any smearing-related wrapper notes consistent as code changes. If behavior, parameter meaning, file names, defaults, validation status, or known limitations change, update relevant `.md` files in same work pass.

## Active Energy Response

Default energy mean model:

```text
E_mean(E) = a + b*E_safe + c*ln(E_safe/1 GeV)
E_safe = max(E, E_floor)
```

Convention:

- `E` is original photon energy in GeV.
- `E_mean` is reconstructed-energy mean in GeV.
- `a` has units GeV.
- `b` is dimensionless.
- `c` has units GeV.
- `E_floor` defaults to `0.2 GeV`.

This is not a scale factor. Do not compute `E * E_mean`. Resolution is applied once around `E_mean`.

## Resolution And Position Models

Energy resolution is separate from energy mean response.

Simple default:

```text
sigma_E(E_mean) = sigma * sqrt(E_mean)
```

Optional 3-term form:

```text
sigma_E = sigma * E_mean * sqrt(A^2/E_mean + B^2/E_mean^2 + C^2)
```

Position smearing is additive in x and y:

```text
x_sm = x + Gaus(0, sigma_pos)
y_sm = y + Gaus(0, sigma_pos)
```

Optional energy dependence:

```text
sigma_pos_eff = sigma_pos * sqrt(E0/E_mean)
```

Position smearing must be enabled consistently in fitter and downstream producer.

## Current Optimization

Current active tools:

- Deterministic Sobol global multistart:
  - current implementation uses `ROOT::Math::QuasiRandomSobol` backed by ROOT/GSL.
  - seed scan covers active bounded parameters.
  - best `GLOBAL_MULTISTART_KEEP_BEST` seeds are retained.
- Minuit2 MIGRAD from each retained seed through six-parameter interface:
  - `a`
  - `b`
  - `c`
  - `sigma_E`
  - `sigma_pos`
  - `p_e_scale`
- Minuit2 HESSE diagnostics run after successful improved MIGRAD minima.
- Limited profile scans run for parameters involved in strong HESSE correlations.
- Legacy staged/grid optimizer remains as fallback, but default is global multistart.
- Deterministic RNG/common random numbers so chi2 surface is reproducible and smoother.
- Baker-Cousins or Pearson chi2, depending on config.
- If Minuit returns a finite lower-chi2 point with a non-success status, the point is kept and the status flag records that Minuit did not fully converge.
- Parallel section fitting is enabled by default within each coupled sweep. Sweeps are barriered: sweep N+1 starts only after every section in sweep N finishes.
- `smearing_closure_summary.csv` compares the optimizer chi2 before final model-state recomputation with the final chi2 used for output plots/maps.

Current method is stronger than the old staged grid, but still not proof of a unique global physics solution. Use diagnostics and closure checks.

Current section orchestration:

- No Pass-1/Pass-2 prefit structure.
- Sweep 1 starts with nominal out-of-section photon response: `a=0`, `b=1`, `c=0`, `sigma=0`, `sigma_pos=0`.
- Each sweep fits every active section independently, so section fits can run in parallel.
- After all sections in sweep N finish, cross-boundary/out-of-section photon parameters are refreshed from completed sweep-N results.
- Sweep N+1 then fits in-section parameters while using sweep-N values for out-of-section photons.
- There is no final serial Gauss-Seidel sweep. The final accepted state is the last completed barriered sweep, then final status labeling and chi2 recomputation.
- This avoids section-order bias where the last serially fitted section sees fresher external coefficients than the first.

Current default bounds for the energy nonlinearity terms are conservative:

```text
a in [-0.20, 0.20] GeV
c in [-0.20, 0.20] GeV
```

Wider bounds can be used for studies, but they make it easier for the optimizer to create unphysical low-energy response and bad local minima.

## Event Weighting And Exclusivity

Fitter can optionally multiply event weights by exclusivity flags:

```text
APPLY_IS_EXCLUSIVE_SELECTION = true
```

Data and simulation can use different branch names:

```text
DATA_EXCLUSIVITY_BRANCH = "is_exclusive_ellipse_combined"
SIM_EXCLUSIVITY_BRANCH  = "is_exclusive_ellipse_combined"
```

When exclusivity weighting is enabled, both configured branches must exist. Missing data branch or missing simulation branch is a hard error. This avoids silent use of all events in one sample while gating the other.

## Dimensionality And Degeneracy

With energy and position smearing enabled, active fit can become 5D:

```text
a, b, c, sigma_E, sigma_pos
```

With electron momentum scaling enabled, it becomes 6D:

```text
a, b, c, sigma_E, sigma_pos, p_e_scale
```

Degeneracies:

- `a`, `b`, `c` change photon energy mean.
- `sigma_E` broadens energy-dependent mass shapes.
- `sigma_pos` broadens opening angle and therefore affects `M_gg`.
- `p_e_scale` shifts missing mass.
- `M_gg`, `M_miss`, and `(p_target+gamma gamma)^2` constrain overlapping combinations, not orthogonal parameters.
- Several parameter sets can produce similar chi2.

Therefore one best chi2 point is not automatically unique physics. Treat fit output as one solution inside a constrained response model, then diagnose stability.

## Algorithms Best Suited

Best practical strategy for this analysis:

```text
bounded quasi-random global seeds
-> keep best N seeds
-> Minuit2 MIGRAD local minimization from each seed
-> HESSE covariance/correlation
-> profile scans for suspicious parameters
-> closure and systematic checks
```

Recommended global seed methods:

- Sobol sequence or Latin-hypercube sampling:
  - good bounded 5D/6D coverage
  - deterministic possible
  - cheaper than full Cartesian grid
- Random restarts:
  - simple
  - less uniform than Sobol/Latin-hypercube
- Coarse Cartesian grid:
  - robust in low dimensions
  - grows too fast in 5D/6D
  - useful only for subset or diagnostic slices

Recommended local minimizers:

- Minuit2 MIGRAD:
  - already used
  - good for smooth local minimization
  - needs good seeds
- Minuit2 SIMPLEX:
  - useful if MIGRAD struggles with rough/noisy surfaces
  - slower, derivative-free
- Minuit2 HESSE:
  - covariance/correlation diagnostics after minimum
  - needed to expose degeneracy

Recommended global/rough optimizers if ROOT/tooling supports them:

- differential evolution
- CMA-ES
- particle swarm
- basin hopping with local Minuit step

These are better for nonconvex surfaces but more expensive and harder to maintain in ROOT macro style.

## What Is Missing Now

Implemented diagnostics:

- `smearing_optimizer_summary.csv`
- `smearing_optimizer_seeds.csv`
- `smearing_optimizer_profiles.csv`

Still missing:

- true Sobol generator through ROOT/GSL
- full covariance matrix output, beyond strongest correlation summary
- broad profile likelihood scans for every parameter in every section
- broad `a/c/sigma_pos/p_e_scale` degeneracy maps
- section-neighbor smoothness regularization
- MCMC/nested sampling for posterior-like degeneracy mapping
- automated closure-test mode with injected known smearing

## Stage 1 Coarse Optimization

Stage 1 coarse optimization is a seed finder, not final physics proof.

Uses:

- gives robust starting `b` and `sigma_E`
- avoids arbitrary `b=1`, `sigma=0.05` start
- helps jump over local roughness from finite smearing statistics
- separates first energy estimate from position smearing

Limits:

- freezes `a=0` and `c=0`
- freezes `sigma_pos=0`
- can bias initial seed toward scalar energy response
- partly redundant with later Stage 3/MIGRAD refinement

Current conclusion:

- keep Stage 1 as robust seed helper for now
- do not interpret Stage 1 result as final response
- test redundancy by comparing staged vs simultaneous mode on same sections
- if final chi2/parameters stable, Stage 1 is only speed/seed helper
- if final chi2/parameters differ, optimizer remains path-dependent

## Modular Model Design

Current fitter uses:

- `ModelParameter`: name, value, bounds, unit, documentation.
- `SmearingModel1D`: common `evaluate(x)` interface.
- `PhotonSmearingParameters`: photon response parameters.
- `PhotonSmearingModel`: energy mean, energy sigma, position sigma.

Event loop should depend on model interface, not math form. To add a model:

1. Add model type.
2. Implement `evaluate`.
3. Define parameter names, units, bounds.
4. Select/build model in config layer.
5. Leave event loop, histogram filling, downstream analysis unchanged.

Conceptual supported forms:

- constant
- linear
- polynomial
- `a_plus_bE_plus_clnE`

## Outputs

Fitter outputs:

- `out_smear.root`
- `section_map.csv`
- `out_smear_interpolated.root`
- `chi2_scans.pdf`
- `smearing_optimizer_summary.csv`
- `smearing_optimizer_seeds.csv`
- `smearing_optimizer_profiles.csv`

Current CSV energy columns:

```text
best_mu_a,best_mu_b,best_mu_c,best_sigma,best_sigma_pos_cm
```

Interpolated maps:

```text
h_mu_a_interp
h_mu_b_interp
h_mu_interp
h_mu_c_interp
h_sigma_interp
h_sigma_pos_interp
```

`h_mu_interp` is compatibility alias for `h_mu_b_interp`.

## Running

Default:

```bash
./run_smearing_pipeline.sh
```

Useful overrides:

```bash
DATA_FILE=/path/data.root \
SIM_FILE=/path/sim.root \
OUT_DIR=/path/output \
NX=11 NY=16 NSMEAR=80 \
./run_smearing_pipeline.sh
```

Pipeline stages:

1. Validate inputs and ROOT tools.
2. Archive previous canonical outputs into `smearing_pipeline_backup/<RUN_TAG>/`.
3. Compile `scripts/nps_sim_smearing_new_try.C`.
4. Run section smearing fit.
5. Verify ROOT, CSV, interpolated maps, PDF.
6. Run `scripts/simc_pi0_analysis.C` with explicit smearing file paths.
7. Verify smeared SIMC ROOT output.

## Plot Interpretation

Per-section PDF pages compare data, unsmeared simulation, and smeared simulation for:

- `M_gamma gamma`
- `M_miss`
- `(p_target+gamma gamma)^2`

Read plots with caution:

- 1D chi2 slice is not a full profile.
- Flat slice means weak constraint or degeneracy.
- Deep local minimum does not prove global minimum.
- Compare energy-response curves, not just coefficients.
- Check section maps for noisy/nonphysical structure.

## Validation Checklist

Minimum checks:

- compile fitter
- run pipeline or small run
- confirm `section_map.csv` current columns
- confirm `out_smear_interpolated.root` maps
- confirm downstream smeared ROOT tree has entries
- compare pi0 mass peak/width
- compare missing-mass peak/width
- compare `(p_target+gamma gamma)^2`
- inspect energy-response curves
- inspect position maps if enabled
- check kinematic distributions used in final analysis
- rerun notebook comparison after regenerating outputs

Optimizer stability checks:

- staged vs simultaneous mode
- multiple RNG seeds if stochastic choices change
- Sobol seed stability checks
- varied histogram binning/ranges
- varied observable weights
- closure with known injected smearing
- per-section neighbor smoothness check

## Common Mistakes

- treating `a+bE+c ln(E)` as scale factor
- applying pair-level `m_gg` correction together with fitted energy response without documented systematic
- fitting `sigma_E` when energy mean model already compensates bad calibration
- letting `p_e_scale` absorb photon-energy mismatch
- using different data/sim exclusivity definitions without documenting why
- comparing old notebook outputs after new pipeline run
- using stale interpolated maps after new CSV
- enabling position energy dependence downstream but not in fitter
- trusting one minimum without degeneracy diagnostics

## Current Validation Status

Local checks passed:

- `bash -n run_smearing_pipeline.sh`
- `git diff --check`
- notebook JSON validation
- ROOT compile validation with `csh -c 'source /group/nps/singhav/setup.csh; ...'`
- SIMC producer compile validation after fixing ROOT current-directory metadata writes.

Full pipeline validation still needs a small production input run after choosing the intended data/SIMC files.

## Open Work

- Add full HESSE covariance/correlation matrix output from Minuit2.
- Broaden profile scans for `a`, `b`, `c`, `sigma_E`, `sigma_pos`, `p_e_scale`.
- Add optimizer summary table to `section_map.csv` or separate diagnostics file.
- Add staged-vs-simultaneous comparison mode.
- Add closure-test mode with injected known smearing.
