# ALG-002A combined-statistics timing-background pilot

Date: 2026-10-01  
Workspace: `root_analysis_env_final`  
Validation setting: `KinC_x36_4`  
Status: **implemented shadow pilot; not promotable; no pi0 weight**

## Purpose and approval boundary

The legacy analysis derives a mass-bin `pi0_weight` independently for each run
after timing subtraction and a per-run Fermi/logistic ("PEPSI-motivated")
combinatorial-background fit. The audit and the run-6407 diagnostic showed that
low-statistics fits, clipping, incomplete timing support, and fixed geometric
transfer factors prevent those weights and conditional `sum(w^2)` errors from
representing the complete inference uncertainty.

The user approved ALG-002A on 2026-10-01 as the timing-only first step toward a
combined-statistics estimator, with these constraints:

1. Combine information within one setting without merging away run identity.
2. Preserve acquisition-mode and selected-cluster multiplicity strata.
3. Give each run timing center/width deviations and free component yields,
   constrained by the ensemble when statistics are weak.
4. Leave the true-coincidence mass spectrum unparameterized here.
5. Do not calculate or write a new `pi0_weight`.
6. Do not read, modify, or calculate efficiencies; form no cross section.
7. Keep the legacy production analysis and canonical outputs unchanged.
8. Record thesis-grade mathematical, numerical, and operational reasoning.

The implementation is under `src/background_fit/`. No production macro,
selection, pair choice, efficiency source, combiner, smearing code, or
cross-section code was edited.

## Why simultaneous fitting is not simple pooling

Simple pooling would assume one center, width, accidental rate, and
signal-to-background ratio for all runs. High-charge runs would dominate and
real calibration or rate changes could be averaged away.

ALG-002A uses a simultaneous likelihood. Timing-shape information is shared,
but yields remain run-specific. A low-statistics run borrows shape information
while retaining uncertainty in its own composition. No charge scaling is
imposed, because that would assume run-independent acceptance and efficiency.

## Observation and binning

The sole event input is ALG-001 `raw_observation`: run/source provenance,
selected-pair times, recorded pair-time cut, acquisition mode, multiplicity,
and invariant mass. Purity and efficiency-derived weights are absent.

Each observation is assigned to

```text
(run r, mode q, multiplicity m, mass bin j, timing cell a,b).
```

Fixed pilot binning:

- selected-pair time domain `[139,161)` ns;
- 1 ns by 1 ns timing cells;
- regular invariant-mass domain `[0,0.4)` GeV;
- regular mass-bin width 0.002 GeV;
- explicit mass underflow and overflow bins.

Flow bins protect provenance closure. Run 6407 contains 60 selected events
outside 0--0.4 GeV; dropping them would silently change the sample.

## Exact accepted timing support

Two-cluster events retain the production behavior: no pair-time-difference
cut. Multiplicity `>=3` is conditioned on the recorded cut (13 ns here).

Every 1 ns square is clipped against

```text
t1 - t2 <= dt_max
t2 - t1 <= dt_max
```

and the accepted polygon is triangulated. Polygon area is exact; component
densities use symmetric degree-two triangle quadrature. Exact support is

```text
two clusters:             22^2 = 484 ns^2
>=3, |dt|<=13 ns:         22^2 - (22-13)^2 = 403 ns^2.
```

A shifted 6 ns by 6 ns full box therefore has 36 ns2 for two-cluster events but
only 4.5 ns2 after the 13 ns cut. The synthetic test verifies these values to
`1e-12 ns2`.

## Extended-Poisson timing model

For one mode/multiplicity stratum,

```text
n[r,j,a,b] ~ Poisson(mu[r,j,a,b])

mu = T[r,j] G00[r,a,b]
   + H[r,j] G0[r,a] U[b]
   + V[r,j] U[a] G0[r,b]
   + R[r,j] U[a] U[b]
   + D[r,j] sum_k Gkk[r,a,b].
```

All five yields are nonnegative.

- `T`: true-coincidence timing component, still containing pi0 plus
  true-coincidence combinatorial mass background.
- `H`, `V`: one prompt-like and one random photon time.
- `R`: two random photon times.
- `D`: correlated diagonal satellite/bucket timing structure.
- `G00`: correlated bivariate Gaussian central timing core.
- `G0`: its one-dimensional marginal.
- `Gkk`: correlated Gaussian satellites at 140, 142, 144, 156, 158, 160 ns.
- `U`: positive cubic B-spline random-time density.

Clamped spline knots are 139, 143, 147, 149, 151, 153, 157, 161 ns. Log
coefficients receive a second-difference penalty, currently 10.0.
Control-region cross-validation is incomplete and remains a promotion blocker.

Every component is normalized over actual accepted support. No nominal
rectangular timing area is used as a transfer factor.

## Run-specific center and width

For run `r`,

```text
center_r = center_setting + delta_r
sigma_r  = sigma_setting * exp(eta_r).
```

Deviations have exact zero-ensemble sum and Gaussian population penalties with
fitted scales `tau_delta` and `tau_eta`. A small scale means nearly shared run
parameters; a larger scale permits heterogeneity.

### Rejected v1 and balanced v2 contrasts

V1 represented zero sum as `[free, -sum(free)]`, privileging the final run.
Run 6568 acquired a spurious `+8.074 ns` offset to compensate 55 mostly
negative coordinates. Optimizer convergence did not make that scientifically
valid.

V2 uses a balanced binary-tree orthonormal basis spanning the zero-sum
subspace. Columns have unit norm and no accumulator run. V1 is retained as
rejected evidence rather than hidden or overwritten.

## Numerical solution

Timing shapes are fitted after integrating over mass, assuming component timing
density is mass-independent. With shapes fixed, five nonnegative yields are
profiled independently in each run/mass bin.

Inner yield solution:

1. Expectation-maximization supplies a nonnegative start.
2. Slow, nearly-collinear cases use the same convex extended-Poisson objective
   with L-BFGS-B and exact analytic gradient.
3. SLSQP is an independent fallback.
4. A fallback is accepted only when its nonnegative-yield KKT residual is below
   `2e-6`; otherwise the fit fails.

The outer fit uses bounded L-BFGS-B, run-local finite differences for balanced
contrasts, global finite differences for shared shapes, and fail-closed
convergence. No tables are written until every stratum converges.

Component densities are evaluated in log space and shifted by a
component-wide maximum before exponentiation. Normalization makes this
mathematically equivalent while avoiding underflow for narrow correlated trial
Gaussians.

Inverse L-BFGS profile curvature is stored as a numerical diagnostic, not a
calibrated confidence covariance. Conditional five-yield Fisher ranks are
reported rather than silently regularized.

## Output contract

The launcher refuses paths at or below `FINAL/output`. A shadow bundle contains:

| File | Content |
|---|---|
| `FIT_STATUS.txt` | `ALG002A_SHADOW_NOT_PRODUCTION` marker. |
| `provenance.json` | Command, versions, git state, config, run list, summaries, exclusions. |
| `input_manifest.csv` | Paths, sizes, SHA-256 values, runs, entries. |
| `run_timing_parameters.csv` | Run deviations, population scales, integrated yields. |
| `component_yields_by_run_mass.csv` | Five yields, conditional covariance, rank, flow labels. |
| `true_coincidence_mass_spectrum.csv` | Projection of unparameterized `T[r,m,j]`. |
| `timing_predictions.csv` | Accepted area, observation, expectation, residuals per cell. |
| `mode*_mult*_profile_covariance.npz` | Parameters, curvature, probabilities, yields, runs. |
| `run_coverage.csv` | Expected versus represented runs. |
| `validation_gates.json` | Machine-readable pass/fail/warning/pending gates. |

No file or branch named `pi0_weight` is produced.

## Commands actually executed

Working directory:

```text
/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final
```

### Synthetic regression

```bash
python3 -m py_compile \
  src/background_fit/joint_timing_model.py \
  src/background_fit/run_joint_timing_fit.py \
  src/background_fit/validate_joint_timing_output.py \
  tests/test_joint_timing_background.py
python3 tests/test_joint_timing_background.py
```

Final result: `joint timing background synthetic test: PASS`.

This checks geometry, opposite run offsets, count/flow closure, raw ROOT
loading, noncanonical output, and absence of pi0-weight/efficiency output. It is
a mechanics regression, not a coverage ensemble.

### KinC_x36_4 raw export

```bash
csh -c 'source /usr/share/Modules/init/csh; \
source /group/nps/singhav/setup.csh; \
setenv NPS_SELECTION_REPORT_CSV \
/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/efficiency_stuff/selection_report_KinC_x36_4.csv; \
bash src/analysis/run_parallel_nps_analysis_main.sh \
  --source waveform --mode waveform --types production \
  --kin KinC_x36_4 --target LH2 --jobs 8 --run-only \
  --gevnum-cut yes --timeout 1800 --raw-observation-export \
  --output-base validation/runtime/alg002a_x36_4_raw_all'
```

This read an existing good-event-range report for unchanged selection. It did
not run or alter efficiency calculations. `--run-only` bypassed finalization.

Observed:

- 57 configured production LH2 runs;
- 56 successful files and 76,242 observations;
- 74,130 two-cluster and 2,112 `>=3`-cluster observations;
- all waveform mode with a 13 ns recorded pair cut;
- run 6569 missing because no waveform ROOT input exists.

### Corrected setting-wide fit

Thread caps remove measured OpenBLAS oversubscription without changing the
likelihood:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 src/background_fit/run_joint_timing_fit.py \
  --input 'validation/runtime/alg002a_x36_4_raw_all/KinC_x36_4/root/diagnostics_run*.root' \
  --output-dir validation/runtime/alg002a_x36_4_fit_all_v2 \
  --maxiter 300 --tolerance 2e-7
```

```bash
python3 src/background_fit/validate_joint_timing_output.py \
  validation/runtime/alg002a_x36_4_fit_all_v2 \
  --config-csv config/nps_dvcs_all_kins_main.csv \
  --kin KinC_x36_4 --target LH2 --types production
```

Validator exit status 3 is expected while blockers remain.

## KinC_x36_4 v2 measured results

These diagnose timing only; they are not pi0 yields or cross sections.

| Quantity | Multiplicity 2 | Multiplicity >=3 |
|---|---:|---:|
| Events | 74,130 | 2,112 |
| Outer iterations | 84 | 95 |
| Shared center (ns) | 149.90190 | 149.87969 |
| Shared marginal sigma (ns) | 0.92926 | 0.96858 |
| Central correlation | 0.97385 | 0.97249 |
| Curvature condition | 1.951e5 | 3.660e5 |
| Run offset range (ns) | [-0.0371,+0.0435] | [-0.00501,+0.00447] |
| Run sigma range (ns) | [0.9150,0.9406] | [0.9666,0.9706] |

Run-offset and log-width population scales sit at their lower bounds. The
likelihood supports nearly shared parameters, but variance-component inference
at a boundary is nonregular, so this is reported as a warning.

For multiplicity 2, direct prompt means and fitted centers have correlation
0.553 and count-weighted RMS difference 0.0656 ns. The sparse `>=3` stratum is
more strongly shrunk, as intended.

Count closure:

```text
raw input                76,242
mass-stratified rows     76,242
timing prediction rows  76,242
summed fitted expected   76,241.9999915
```

## Gate result and interpretation

`validation_gates.json` reports `NOT_PROMOTABLE`.

Passed:

- exact integer-count closure;
- successful outer and KKT-checked inner optimization;
- positive curvature with condition below `1e8`;
- balanced run-deviation closure;
- no pi0-weight output.

Failed:

1. Run coverage: 57 expected, 56 represented; run 6569 is missing.
2. Mass-bin identifiability: only 1,418/10,848 occupied bins have rank five
   (`13.07%`), and 4,741 occupied bins put the conditional true MLE at zero.

Pending:

- control-region selection of the spline penalty;
- synthetic-null calibration of the predeclared central timing tail;
- 2,000-toy bias, pull, coverage, and failure ensemble;
- toy-calibrated run/multiplicity goodness of fit.

The rank failure is decisive for `pi0_weight`: combined timing shapes alone do
not identify five independent yields in each individual run and 2 MeV mass
bin. A shared validated mass model is needed next. Producing weights now would
conceal non-identifiability.

Residual maps contain off-diagonal prompt-like tails. The proposal anticipated
an optional symmetric tail, but it must not be enabled merely because it
improves these data; its likelihood-ratio selection must first be calibrated
under synthetic no-tail samples.

## Uncertainty status

The output distinguishes:

1. conditional per-bin Fisher covariance with timing shapes fixed;
2. diagnostic inverse-LBFGS curvature for shared/run timing parameters;
3. missing calibrated joint uncertainty across runs and mass bins.

Only the first two exist. Neither is a publication interval. The later model
must retain shared nuisance correlations through low-rank covariance factors or
profile/toy replicas; independent event-weight errors would be incorrect.

## Retained intermediate evidence

- `validation/runtime/alg002a_x36_4_run6407`: initial single-run fit exposing
  the original correlation bound and sparse-bin rank loss.
- `validation/runtime/alg002a_x36_4_fit_all_v1`: rejected asymmetric-basis fit;
  run 6568 acquired `+8.074 ns` and this bundle must not be used.
- `validation/runtime/alg002a_x36_4_fit_all_v2`: corrected fit described here;
  still shadow-only and not promotable.

## Reference integrity

Post-work checksums:

```text
fafa5c2db7ecaceff4db369b5b2f9cd622d084039e1bea9b781dcdbc9654454d  docs/proposals/ALG-002_joint_timing_background_likelihood.md
4e929083a53781daa3a7972a6fca416a246da04e97d2077a47f5d2a9f71f2133  src/analysis/nps_raw_observation.h
5623dffa447a6bdfe16a1a17571c0f39638bb45ae19945400098cce68d310b58  src/analysis/nps_analysis_main.C
9f0df3ea13872358dc7766f7e0d75d29079e0418f846ecb841ab2ee536f64ddb  frozen pi0 uncertainty audit
```

These proposal/analysis sources remained read-only references.

## Next decision, not implemented

The next dependency is a simultaneous true-coincidence mass model with shared
signal/combinatorial shapes, run-specific yields, validated peak shifts/widths,
and full nuisance propagation. That is a separate physics update and requires
explicit approval before implementation.

