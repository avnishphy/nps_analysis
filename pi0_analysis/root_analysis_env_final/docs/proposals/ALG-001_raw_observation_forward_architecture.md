# ALG-001: raw-observation forward-inference architecture

Status: **awaiting explicit user approval; not implemented**

Decision requested: **Approve ALG-001 as described?**

## Decision and bounded scope

Authorize a parallel, opt-in prototype for one kinematic setting whose
statistical observations are disjoint raw event-count categories, rather than
clipped/subtracted mass templates or events carrying fitted purity weights.

Approval would authorize only:

1. passive export of stable run/segment/event identifiers, photon times,
   multiplicity stratum, and one exclusive timing-category code without
   changing event selection, loop order, random-number consumption, legacy
   branches, or default outputs;
2. an explicit run/segment ledger containing charge, prescale, processing
   status, correction central values, auxiliary counts when already available,
   and zero-candidate runs;
3. a new, non-default likelihood prototype under a separate output namespace;
4. synthetic closure/coverage tests and a one-setting shadow comparison that
   cannot overwrite existing outputs.

This approval would **not** authorize a production default switch, publication
claim, efficiency-definition change, Geant4 modification, replacement timing
transfer model, replacement combinatorial-background function, final smearing
model, full campaign, or overwrite of any current result. Each of those remains
a later approval gate.

## Why this decision comes first

The current chain estimates timing background, combinatorial background,
mass-bin purity, run corrections, and detector smearing at different stages.
Downstream extraction receives only fixed event weights plus conditional
`sum(w^2)`. Once raw category membership, auxiliary counts, run exposure, and
shared fit parameters have been discarded, a downstream covariance patch
cannot reconstruct all correlations.

The architecture choice therefore precedes detailed repairs. Implementing a
new timing formula first would risk producing another summary whose uncertainty
contract is already too narrow.

## Proposed observation model

For run `r`, reconstructed analysis bin `b`, multiplicity stratum `m` (exactly
two versus at least three accepted clusters), disjoint timing category `c`, and
mass bin `j`, record the integer count

```text
n[r,b,m,c,j] ~ Poisson(mu[r,b,m,c,j]).
```

The expectation is decomposed before any subtraction:

```text
mu = mu_signal(sigma, exposure, detector_response)
   + mu_accidental(timing_controls, transfer_nuisance)
   + mu_combinatorial(mass_controls, shape_nuisance)
   + optional validated residual components.
```

The categories are mutually exclusive and exhaustive inside the prototype's
declared timing domain. Signal contributes to prompt categories through the
detector response. Sideband and full-box counts constrain accidental nuisance
parameters. Mass sidebands/control regions constrain combinatorial nuisance
parameters. No category count is negative, no subtracted bin is clipped, and
no inferred purity is treated as an observed event weight.

The first prototype keeps the current efficiency and livetime central values
unchanged. It carries the available auxiliary values and correlations as
metadata; a specific likelihood/constraint model for those quantities requires
a later proposal. Likewise, it reads the current smearing map only to reproduce
a baseline shadow prediction. Making smearing parameters joint nuisance
parameters requires a later proposal.

## Exposure and response contract

The run exposure entering the forward expectation is represented once:

```text
E_r = Q_r(mC) L_r e_r / P_r,
```

with a ledger row for every intended valid run, including zero selected
candidates. Target correction and virtual-photon/SIMC normalization enter once
in named response factors. The implementation must provide a dimensional table
showing that the predicted observation has units of events for every term.

The detector response is accumulated by original generator event. Repeated
smears from one generator event remain a correlated integration device; their
finite-MC moments are grouped at generator-event level. Exact response nuisance
treatment is deferred, but the prototype interface must not prevent it.

## Alternatives considered

### A. Recommended: simultaneous raw-observation likelihood

Advantages:

- preserves integer observations and control-region information;
- represents shared background, exposure, and response parameters explicitly;
- handles zero-count bins without invented variances or pseudocounts;
- enables coherent profiling, posterior sampling, or parametric bootstrap;
- naturally supports multiplicity- and acquisition-mode-specific timing
  acceptance.

Costs and risks:

- larger model and validation surface;
- needs passive schema additions and a run ledger;
- may be sensitive to misspecified background shapes, so goodness-of-fit and
  alternative models are mandatory;
- higher computational cost than the present weighted least-squares fit.

### B. Signed subtraction plus complete covariance

Retain the current linear subtraction but export raw component vectors and the
full covariance/Jacobian through mass, kinematics, run corrections, and
smearing. This is less invasive and can be a valuable cross-check.

It is not recommended as the primary architecture because the current pipeline
uses clipping, rejects nonpositive weights, re-estimates mass cuts, and shares
fit parameters across many event rows. A valid implementation would first have
to remove nonlinear clipping/rejection from the estimator and maintain a large,
often singular covariance across all downstream bins. Sparse signed summaries
also make the final Gaussian approximation fragile.

### C. Keep fixed purity weights and change only the final objective

Rejected as the primary solution. Gaussian, ordinary Poisson, scaled-Poisson,
or bootstrap treatment of already inferred fixed weights remains conditional on
the upstream background, corrections, and smearing. It cannot restore discarded
raw-category or auxiliary-count information.

## Planned files if approved

No file below is changed by this proposal.

- `src/analysis/nps_analysis_main.C`: passive raw-category and identifier
  branches/tree, guarded so legacy output values and ordering remain unchanged.
- `src/analysis/combine_analysis_branches.py`: retain identifiers and produce a
  complete run ledger; legacy combined tree remains available.
- New `src/xsec_extract/raw_observation/` sources: schema validator, likelihood,
  fit runner, and synthetic generators.
- `src/xsec_extract/run_xsec_pipeline.sh`: explicit opt-in mode only; current
  default unchanged.
- New tests and validation manifests under `tests/` and `validation/`.
- Documentation and worklog updates recording formulas, units, commands,
  decisions, and results.

Before editing an existing file, its checksum and recoverable pre-edit copy
will be retained as required by repository policy.

## Prospective validation contract

These criteria are part of the requested approval and will not be relaxed after
looking at results without a documented amendment.

### Engineering invariants

- Export disabled: legacy outputs match the approved baseline bit-for-bit where
  deterministic, otherwise branch-by-branch numerically with exact entry order.
- Export enabled: legacy branches, selected-entry count, event order, weights,
  and RNG trace remain identical; only new passive records appear.
- Timing categories are exclusive and exhaustive over the declared domain, and
  their sums reproduce source-loop counts by run and multiplicity.
- The run ledger matches the intended runlist and contains zero-candidate and
  partial-status rows explicitly.
- Unit tests prove the charge cancellation and prevent a duplicate `1000`,
  target, flux, or bin-width factor.

### Statistical tests

- Analytic one-bin Poisson and on/off limiting cases agree with independent
  calculations.
- Synthetic closure includes zero counts, boundary solutions, timing-transfer
  variation, background-shape variation, migrations, finite MC, zero-candidate
  runs, partial segments, and correlated normalization nuisances.
- At least 2,000 successful pseudoexperiments per declared benchmark point;
  every failure is counted, not silently dropped.
- Nominal 68% interval coverage must be within 3 percentage points and nominal
  95% coverage within 2 percentage points for every primary fitted component.
  With 2,000 trials these bands are about three or more binomial standard
  errors at the nominal coverages.
- Median fitted bias magnitude must be below 0.10 of the ensemble statistical
  standard deviation; pull mean magnitude below 0.10 and pull width within
  `[0.90,1.10]`.
- Fit failure/nonconvergence below 1%; otherwise the benchmark fails and no
  coverage claim is made.
- A deliberately misspecified timing or mass model must be detected by a
  predeclared goodness-of-fit/stress diagnostic often enough to expose, rather
  than hide, the tested discrepancy. Exact power targets require the later
  background-model proposal.

### One-setting shadow comparison

For one user-selected setting, report—not tune away—differences in selected
counts, exposure, background components, cross-section central values,
uncertainties, correlations, fit quality, and support. A central-value change
is neither automatic failure nor evidence of improvement; it must be traced to
a named estimator difference. No production output is overwritten.

## Expected impact and unknowns

Expected direction: uncertainties should become more complete because shared
background/exposure/response information can contribute. Their numerical size
can increase or decrease after correlations and control constraints are
included. Central values may move because clipping, nonpositive rejection, and
mass-only purity assignment are absent from the new estimator.

No magnitude is claimed. The production impact is unresolved until the
one-setting shadow run and pseudoexperiment suite pass. Runtime and memory are
also unresolved; profiling will be reported before any campaign proposal.

## Reproducibility, rollback, and failure conditions

- New mode is opt-in and writes to a distinct directory containing config,
  commit, input/output hashes, environment, run ledger, fit diagnostics, and
  applied RNG seeds.
- The legacy default remains unchanged. Rollback is disabling the new mode or
  reverting its isolated commits; historical outputs are never overwritten.
- Failure of schema invariants, closure, coverage, unit checks, convergence, or
  source-to-ledger reconciliation blocks promotion.
- A later change to timing/background parameterization, efficiency behavior,
  smearing likelihood, response statistics, or publication default requires a
  new explicit proposal and approval.

## Approval request

**Approve ALG-001 as described?**
