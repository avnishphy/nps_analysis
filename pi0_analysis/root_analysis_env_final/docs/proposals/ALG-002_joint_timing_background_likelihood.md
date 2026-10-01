# ALG-002: joint two-photon timing-background likelihood

Status: **awaiting explicit user approval; not implemented**

Decision requested: **Approve ALG-002 as described?**

## Scope and efficiency exclusion

Authorize a parallel, non-default fit of integer observations in two-photon
timing and invariant mass for one setting. Its first purpose is to replace the
fixed-area accidental subtraction in a shadow analysis and quantify timing
background uncertainty without clipping or purity weights.

Efficiency calculations, definitions, files, corrections, and nuisance models
are explicitly excluded. This fit uses raw selected counts only. It neither
reads nor modifies an efficiency value.

Approval would authorize a new implementation and synthetic/one-setting shadow
validation. It would not authorize replacing the production background,
changing event/pair selection, changing the combinatorial-mass model, changing
smearing or extraction, running the full campaign, or making a publication
claim.

## Observation and fixed stratification

For one setting, retain the existing selected photon pair and form integer
counts in fixed cells:

```text
n[r,m,j,a,b]
```

where:

- `r` is run;
- `m` is multiplicity (`2` versus `>=3` accepted clusters);
- `j` is the existing 2 MeV invariant-mass bin;
- `a,b` are 1 ns bins in selected-pair `t1,t2`;
- acquisition mode is fitted separately if a setting contains both modes.

The timing domain is `[139,161] ns` so shifted waveform windows are fully
represented. Each stratum uses the exact existing pair-selection acceptance:
the two-cluster stratum has no pair-difference cut; the `>=3` stratum is
conditioned on its resolved `|t1-t2|` cut. No nominal rectangular area is used
where the accepted shape is triangular.

Phase one fits timing and mass integrated over final physics bins. Differential
phi/t/Q2/xB extension is a later proposal after the timing model passes. This
prevents a timing pilot from silently becoming a cross-section estimator.

## Proposed extended-Poisson model

For every retained cell:

```text
n[r,m,j,a,b] ~ Poisson(mu[r,m,j,a,b])

mu = T[r,m,j] * G00[r,m,a,b]
   + H[r,m,j] * G0[r,m,a] * U[r,m,b]
   + V[r,m,j] * U[r,m,a] * G0[r,m,b]
   + R[r,m,j] * U[r,m,a] * U[r,m,b]
   + D[r,m,j] * sum_k Gkk[r,m,a,b].
```

All component rates `T,H,V,R,D` are nonnegative:

- `T` is the true-coincidence mass spectrum. It intentionally includes pi0 and
  true-coincidence combinatorial mass background; ALG-002 does not yet separate
  them.
- `H,V` describe one prompt-like and one random photon time.
- `R` describes two independently random times.
- `D` describes correlated diagonal satellite/bucket structure.
- `G00` is a bivariate central prompt-time shape with a fitted correlation;
  `G0` is its marginal used in one-random components. `Gkk` are correlated
  bivariate satellite/bucket shapes, and `U` is a smooth random-time density.

Every component is normalized over the actual retained lattice cells of its
mode/multiplicity stratum after the existing pair-selection mask. This makes
the transfer a likelihood integral over accepted support rather than an area
ratio.

## Shape parameterization and constraints

- `G00` and each supported `Gkk`: bin-integrated bivariate Gaussian core with a
  fitted within-event time correlation, plus one optional symmetric tail
  parameter. A tail is retained only if its predeclared
  likelihood-ratio calibration passes synthetic null tests; otherwise the
  Gaussian-only form is used.
- `U`: positive cubic B-spline in time with knots fixed at 139, 143, 147, 149,
  151, 153, 157, and 161 ns. Coefficients use a log parameterization and a
  second-difference Gaussian smoothness constraint.
- Central timing location may vary by run. Run offsets receive a shared normal
  constraint whose width is estimated from the ensemble of run control data;
  mode/multiplicity widths are shared unless a predeclared heterogeneity test
  rejects sharing.
- Satellite centers are fixed to observed window/bucket locations during the
  pilot; amplitudes are fitted. Adding or removing a satellite after looking at
  signal-region residuals is forbidden without an amended proposal.
- Mass-bin component rates are independent nonnegative parameters in the pilot;
  no pi0 line shape or current Fermi/logistic combinatorial function is imposed.
  A weak second-difference penalty may be used only for `H,V,R,D`, with its
  strength selected by control-region cross-validation. `T` remains unpenalized
  so timing inference cannot sculpt the true-coincidence mass spectrum.

The fit uses all timing cells, including those labeled `outside` by the legacy
region scheme. Legacy region codes remain diagnostics, not the likelihood's
sampling bins.

## Alternatives considered

### Recommended: full accepted 2D timing lattice

Uses the available timing shape, zero cells, and multiplicity-dependent mask.
It directly estimates shared shape/rate nuisance covariance and does not need
negative pseudo-counts.

### Region-count on/off likelihood

Fit only prompt/diagonal/horizontal/vertical/full counts with transfer-factor
nuisances. This is simpler and will be implemented as an independent
cross-check, but it discards within-region gradients and has weaker diagnostics
for shifted-window truncation and pair-selection geometry.

### Current subtraction with full covariance

Propagate raw control counts through the existing linear formula, without
clipping, and carry its covariance. This remains a useful benchmark. It is not
the recommended primary model because one fixed transfer still cannot describe
both multiplicity masks, and the current formula mixes related control terms.

## Planned implementation if approved

- New `src/background_fit/` implementation and schema; no modification of the
  legacy subtraction helper.
- Read only ALG-001 `raw_observation` trees and segment/run identifiers.
- Separate output directory containing configuration, counts, fit parameters,
  covariance/correlation, profiles, predictions, residuals, and provenance.
- New synthetic generator and regression tests under `tests/`.
- An explicit launcher mode that remains off by default and refuses canonical
  production output paths.

No candidate code is installed before approval.

## Prospective validation and pass/fail criteria

### Exact and limiting checks

- One-component and on/off limits reproduce analytic Poisson solutions.
- Integrals over accepted 1 ns cells reproduce enumerated geometry for the
  two-cluster, HCANA `>=3`, and waveform `>=3` masks.
- Summed lattice counts equal ALG-001 selected observations exactly by run,
  mass bin, multiplicity, and mode.
- Legacy region projections reproduce the raw source histogram counts where
  histogram coverage is identical; shifted out-of-range counts are reported
  separately.

### Synthetic ensembles

- At least 2,000 successful pseudoexperiments per benchmark; all failures count.
- 68% and 95% profile-interval coverage within 3 and 2 percentage points,
  respectively, for `T` integrals and prompt-region accidental predictions.
- Median bias magnitude below `0.10` ensemble standard deviations; pull mean
  magnitude below `0.10`, pull width in `[0.90,1.10]`.
- Fit failure/nonconvergence below 1%.
- Stress samples vary timing offsets, widths, spline curvature, satellite
  strength, multiplicity mix, sparse/zero cells, and intentionally correlated
  random times.

### Model criticism

- Predeclared held-out timing stripes test interpolation of `U`.
- Pearson/deviance residual maps, marginal residuals, and run/multiplicity
  goodness-of-fit are reported with toy-calibrated p-values.
- Injection tests quantify power against an omitted correlated component. A
  failed stress test blocks promotion rather than being absorbed by uncertainty
  inflation.
- Compare full-lattice, region-count, and unclipped linear-subtraction estimates
  on the same raw observations.

### One-setting shadow

Run only after the user selects the setting/run and a separate output path.
Report component yields, prompt accidental prediction, true-coincidence mass
spectrum, covariance, fit quality, and differences from the legacy estimate.
Do not form a cross section in ALG-002.

## Expected impact and failure modes

Expected qualitative impact: remove clipping and nominal-area transfer, expose
run/multiplicity correlations, and generally widen uncertainty relative to a
fixed background weight. Correlations or strong controls could reduce some
bins. No numerical shift is predicted.

Principal risks are component non-identifiability, spline over-flexibility,
missed correlated timing structure, and insufficient per-run statistics. The
fit must report Hessian rank/condition, boundary parameters, profile failures,
and sensitivity to permitted shape variants. Persistent non-identifiability
fails ALG-002 and favors the simpler region-count/full-covariance path.

## Rollback

The implementation would be isolated, opt-in, and read-only with respect to
ALG-001 and legacy products. Rollback is disabling or reverting the new path.
No production artifact is overwritten.

## Approval request

**Approve ALG-002 as described?**
