# ALG-002B proposal: simultaneous timing and mass model

Date: 2026-10-01
Workspace: `root_analysis_env_final`
Validation setting: `KinC_x36_4`
Status: **proposal awaiting explicit user approval; not implemented**

## Decision requested

Approve an isolated shadow extension of ALG-002A that splits the fitted
true-coincidence timing component into pi0 signal and true-coincidence
combinatorial background by fitting invariant mass and timing simultaneously.

This step will estimate per-run pi0 yields and their joint covariance. It will
not create an event-level `pi0_weight`, modify the production macro, combine
physics distributions, calculate efficiencies, form cross sections, or replace
any canonical output. Event-weight construction and validation remain a later
ALG-002C decision.

## Why this is the next dependency

ALG-002A deliberately left the true-coincidence mass spectrum free in every
run and 2 MeV mass bin. That protected the pi0 peak from an assumed line shape,
but the timing data alone cannot identify five component yields in most sparse
run/mass cells:

- 1,418 of 10,848 occupied cells have full conditional rank;
- 4,741 occupied cells place the true-coincidence yield at zero;
- no pi0/combinatorial separation is possible within the central timing
  component.

The legacy alternative is also inadequate. It fits a Fermi/logistic continuum
independently in each run after timing subtraction. In the 56 represented
`KinC_x36_4` runs, 34 fitted widths reach the 100 MeV upper bound. The fitted
turn mass spans 110--220 MeV and often reaches a boundary. These are symptoms
of weak per-run background-shape information, even when the minimizer reports
success.

The setting-wide data contain useful shared information:

| Quantity | Observed value |
|---|---:|
| Represented runs | 56 |
| Selected observations | 76,242 |
| Prompt-category observations | 29,477 |
| Per-run prompt observations | 41--838; median 564 |
| Legacy fitted pi0 mean | 133.521--135.247 MeV |
| Legacy fitted pi0 width | 3.249--5.581 MeV |
| Legacy mean run-to-run standard deviation | 0.443 MeV |
| Legacy width run-to-run standard deviation | 0.394 MeV |

These numbers support a simultaneous fit with run-specific peak parameters
constrained by setting-wide information. They do not support one merged
histogram with a single peak or one common event weight.

## Input and observation model

The only event input will be the existing ALG-001 `raw_observation` trees from
the `KinC_x36_4` shadow bundle. No legacy `pi0_weight`, efficiency, charge
normalization, or subtracted histogram enters the likelihood.

Each integer observation is assigned to

```text
(run r, acquisition mode q, multiplicity class k,
 mass bin j, accepted timing cell c).
```

The pilot keeps the ALG-002A support exactly:

- timing domain `[139,161)` ns in 1 ns by 1 ns cells;
- exact polygon clipping by the recorded pair-time cut;
- multiplicity classes `2` and `>=3` kept separate;
- regular mass domain `[0,0.4)` GeV in 2 MeV bins;
- explicit mass overflow category for the 3,527 observed events above 0.4
  GeV; the currently empty underflow category remains in the schema;
- every event retains run and source-segment provenance.

The signal flow probabilities are integrals of the signal density outside the
regular domain. Each background component has setting-shared, fitted
underflow/overflow probabilities; the exact mass value inside a flow category
does not enter the fit. Flow probabilities and regular-bin probabilities are
normalized together for each component, so flow events contribute to the
likelihood and closure without forcing a polynomial to describe the mass tail
up to the observed 3.59 GeV maximum.

The full binned extended-Poisson mean is

```text
mu[r,k,j,c] =
    Y_pi0[r,k]  S_pi0[r,k,j] G_true[r,k,c]
  + Y_comb[r,k] B_comb[k,j]  G_true[r,k,c]
  + Y_H[r,k]    B_H[k,j]     G_H[r,k,c]
  + Y_V[r,k]    B_V[k,j]     G_V[r,k,c]
  + Y_R[r,k]    B_R[k,j]     G_R[r,k,c]
  + Y_D[r,k]    B_D[k,j]     G_D[r,k,c].

n[r,k,j,c] ~ Poisson(mu[r,k,j,c]).
```

The timing densities `G_*` are the accepted-support-normalized ALG-002A
components. ALG-002B will initialize from the ALG-002A fit and then refit the
timing and mass parameters jointly. Freezing ALG-002A would omit timing/mass
covariance and is therefore prohibited.

All six component yields are nonnegative and free for every represented run
and multiplicity. They will not be forced to scale with charge, livetime, or
efficiency.

## Pi0 mass shape and run-specific calibration

The primary signal candidate is a normalized double-sided Crystal Ball shape.
It has a Gaussian core plus power-law low- and high-mass tails. Bin
probabilities will be integrals over bin edges, not values at bin centers.

For run `r` and multiplicity `k`,

```text
mu_pi0[r,k] = mu_setting[k] + delta_mu[r]

log sigma_pi0[r,k] = log sigma_setting[k] + delta_log_sigma[r].
```

The run deviations use the same balanced orthonormal zero-sum basis that fixed
the ALG-002A last-run defect:

```text
delta_mu[r]        ~ Normal(0, tau_mu^2)
delta_log_sigma[r] ~ Normal(0, tau_sigma^2).
```

`tau_mu` and `tau_sigma` are fitted population scales. Their boundary behavior
will be calibrated with toys. A low-statistics run therefore receives a
run-specific estimate with an appropriately broad profile uncertainty; its
mean or width is never silently replaced by the setting average.

The Crystal Ball tail parameters are shared across runs. The two multiplicity
classes may have different core means and widths, but the run shift and width
scale are common to both classes because they arise from the same run-level
calibration. A double-Gaussian signal is the predeclared comparator. The
simpler model is retained unless the tail model improves held-out predictive
performance and passes null-toy calibration.

## Mass-background shapes

The true-coincidence combinatorial component starts from three predeclared
positive candidates:

1. the existing PEPSI-motivated Fermi/logistic turn-off;
2. a positive cubic Bernstein density;
3. a positive quartic Bernstein density.

Positive coefficients prevent unphysical negative mass densities. Model
choice will use run-stratified, leave-one-run-out predictive deviance in
background-enriched timing regions and mass sidebands. The simplest candidate
within one standard error of the best predictive score is selected. The pi0
yield in the signal region is not the selection criterion.

The horizontal, vertical, two-random, and diagonal timing components each
receive a positive setting-shared mass density and run-specific normalization.
The initial model shares horizontal and vertical mass shapes, because photon
ordering should not create a physics distinction. A separate-shape alternative
is enabled only if a toy-calibrated likelihood-ratio test and held-out
prediction both reject the shared form.

Multiplicity-dependent background-shape corrections are handled in the same
gated way. This protects the `>=3` class from being forced onto the two-cluster
shape while avoiding an unconstrained model for only 2,112 observations.

No narrow background basis function or knot may be introduced in the pi0 peak
region after inspecting the fitted pi0 yield. Background flexibility is fixed
by the predeclared candidates and control-data validation.

## Factorization check

The model above factorizes each component into a mass density and timing
density. That assumption must be tested because pair selection, especially
the nominal-mass ranking for multiplicity `>=3`, can correlate mass with
timing.

The validator will compare the factorized model with predeclared coarse
mass-by-timing interaction terms. Promotion fails if the interaction model
passes synthetic calibration and gives a material held-out improvement while
changing the setting pi0 yield by more than either 0.25 fitted standard
deviations or 1%, whichever is larger. In that case, the interaction must be
part of a revised proposal rather than added silently.

## Missing run and run identity

The configured production-LH2 set has 57 runs. Run 6569 has no waveform input
in the current source locations. ALG-002B will:

- list all 57 expected runs in `run_coverage.csv`;
- fit the 56 represented runs;
- mark 6569 `excluded_missing_input` with its searched paths;
- create no yield, peak parameter, or weight for run 6569;
- fail if any other run disappears or an unexpected run enters.

An explicit, evidenced exclusion is acceptable for the shadow fit. Recovering
run 6569 later requires a new input-manifest checksum and a full rerun.

## Covariance and uncertainty contract

ALG-002B will retain correlations rather than attach independent errors to
events. It will save:

- the ordered fitted-parameter vector;
- the observed/profile curvature matrix with condition diagnostics;
- profile intervals for per-run and setting-summed pi0 yields;
- the complete cross-run pi0-yield covariance;
- parametric-replica estimates of bias, covariance, pulls, and coverage;
- correlations among timing shapes, mass shapes, run deviations, and yields.

The Hessian covariance is a numerical diagnostic until toy coverage passes.
The replica covariance is the required product for the later run-combination
decision. ALG-002B will not reduce it to independent per-run errors.

## Validation gates

Every gate is fail-closed and recorded in machine-readable JSON.

1. **Exact closure**: observed count, component expectation, mass projection,
   timing projection, and run projection agree to numerical tolerance; all
   76,242 events are accounted for.
2. **Optimizer reproducibility**: at least 20 deterministic dispersed starts
   converge to the same solution within `1e-6` relative likelihood and 0.1
   fitted standard deviations for every setting-summed yield.
3. **No hidden boundaries**: yield or population-scale boundaries are reported;
   no boundary-dependent interval is accepted without calibrated coverage.
4. **Standalone cross-checks**: high-statistics runs are also fitted alone.
   Their peak means, widths, and yields must agree with simultaneous-fit
   profiles within calibrated 95% intervals.
5. **Leave-one-run-out prediction**: every run is predicted from the other 55.
   A globally corrected predictive `p < 0.01` flags the run and blocks
   promotion until explained.
6. **Mass-model selection**: the signal and combinatorial candidates pass
   control-region predictive checks; the selected model is fixed before final
   yield extraction.
7. **Timing-model completion**: ALG-002A spline-penalty cross-validation,
   central-tail null calibration, and per-run toy goodness are completed in
   the joint model.
8. **Factorization**: the predeclared mass-by-timing interaction test passes the
   threshold above.
9. **Synthetic coverage**: at least 2,000 replicas include observed run sizes,
   fitted run shifts, sparse runs, missing run 6569, multiplicity imbalance,
   and parameter boundaries. Aggregate requirements are absolute bias below
   0.2 standard deviations, mean pull magnitude below 0.10, pull width
   0.90--1.10, 68% coverage 0.65--0.71, and 95% coverage 0.93--0.97.
10. **Model-misspecification toys**: data generated from each nonselected
    predeclared signal/background candidate must shift the setting yield by
    less than 0.5 fitted standard deviations. Larger shifts become an explicit
    model systematic or block the model.
11. **Legacy comparison**: report, but do not tune to, legacy per-run yields,
    peak means/widths, run-6407 yield, and the setting sum. Differences are
    decomposed into timing, combinatorial, clipping, and pooling effects.
12. **Output isolation**: canonical `output/` is refused, no `pi0_weight` field
    exists, and no efficiency or cross-section file is read or written.

Passing these gates establishes a validated yield-and-covariance model. It
does not by itself validate weighted phi, t, Q2, xB, or missing-mass spectra.

## Exact implementation boundary after approval

The approved implementation would add:

- `src/background_fit/joint_timing_mass_model.py`;
- `src/background_fit/run_joint_timing_mass_fit.py`;
- `src/background_fit/validate_joint_timing_mass_output.py`;
- `tests/test_joint_timing_mass_background.py`;
- `docs/ALG002B_simultaneous_timing_mass_YYYYMMDD.md` after validation.

It would update only:

- `src/background_fit/README.md`;
- `plan.md`;
- `docs/publication_uncertainty_worklog.md`.

The implementation will use the existing ALG-002A geometry and timing model
without changing ALG-002A results. Any required behavioral change to
`joint_timing_model.py` will be brought back for separate approval.

The shadow output would contain:

- `FIT_STATUS.txt`;
- `provenance.json` and `input_manifest.csv`;
- `model_selection.json`;
- `parameter_estimates.csv`;
- `run_pi0_yields.csv`;
- `pi0_yield_covariance.npz`;
- `profile_intervals.csv`;
- `mass_timing_predictions.csv`;
- `run_coverage.csv`;
- `validation_gates.json`;
- toy summaries and replica covariance under the shadow output directory.

No existing ROOT tree or production CSV schema changes in ALG-002B.

## Reproduction command planned after approval

The initial real-data fit would use the existing shadow bundle:

```bash
cd /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final

env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python3 src/background_fit/run_joint_timing_mass_fit.py \
  --input 'validation/runtime/alg002a_x36_4_raw_all/KinC_x36_4/root/diagnostics_run*.root' \
  --kin KinC_x36_4 \
  --output-dir validation/runtime/alg002b_x36_4_joint_mass_v1
```

The exact seed, model-candidate set, versions, input hashes, and command will be
written to provenance before interpreting results.

## Deferred ALG-002C decision

Only after ALG-002B passes will a separate proposal define event-level fitted
component responsibilities,

```text
w_pi0[i,r] = mu_pi0[i,r] / mu_total[i,r],
```

and test whether using those responsibilities in physics-variable spectra is
valid. Mass/time correlations with phi, t, Q2, xB, and missing mass must be
checked explicitly. Fit uncertainty will be propagated through replicas of the
complete model, not through `sum(w^2)` with fixed fitted weights.

## Method references

- ROOT `RooSimultaneous` documents category-wise simultaneous extended
  likelihoods: <https://root.cern.ch/doc/master/classRooSimultaneous.html>.
- Baker and Cousins derive likelihood-ratio fitting and goodness measures for
  Poisson histogram counts: <https://doi.org/10.1016/0167-5087(84)90016-4>.
- ROOT's generalized double-sided Crystal Ball definition is documented at
  <https://root.cern.ch/doc/master/classRooCrystalBall.html>.
- ROOT documents positive Bernstein polynomial densities at
  <https://root.cern.ch/doc/master/classRooBernstein.html>.
- Pivk and Le Diberder describe fitted component weights and the required
  relationship between discriminating and control variables:
  <https://arxiv.org/abs/physics/0402083>.

## Approval question

Approve implementation of ALG-002B under this exact boundary?
