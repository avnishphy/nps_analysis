# Optional nonnegative angular cross section

The default extraction remains unconstrained. Add `--positive-xsec` to require
the fitted unpolarized angular cross section to be nonnegative. Use
`--no-positive-xsec` to select the original mode explicitly. Both the compiled
executable and `run_xsec_pipeline.sh` accept these switches. The environment
variable `NPS_XSEC_POSITIVE_XSEC` accepts `0` or `1`; the last explicit switch
overrides it. Other effective environment values are rejected.

## What is constrained

In the extractor's convention, the fitted angular cross section is

```text
2*pi * d2sigma/(dt dphi)
  = U + sqrt(2*epsilon*(1+epsilon))*LT*cos(phi)
      + epsilon*TT*cos(2*phi).
```

The entire right-hand side must be >= 0 for every phi. Zero is allowed: this
is nonnegativity, not strict positivity. LT and TT remain signed interference
coefficients. Requiring each coefficient to be positive would impose a
different, unjustified physical restriction.

The constraint applies to all active truth blocks: published bins, exterior
guard regions and any excluded-group truth bins retained as feed-in nuisance
parameters. It uses the maximum generated epsilon among selected MC entries
contributing to that block. The event response continues to use each event's
epsilon and phi. This choice guarantees angular nonnegativity for every
smaller epsilon in the same bin-constant U/LT/TT model, including the displayed
mean-epsilon curve. It does not establish validity beyond the observed epsilon
range or resolve U's physical within-bin epsilon dependence.

For the analytic check, put `z=cos(phi)` and minimize

```text
U - epsilon*TT + sqrt(2*epsilon*(1+epsilon))*LT*z
  + 2*epsilon*TT*z*z,   -1 <= z <= 1.
```

Both endpoints and an interior minimum, when present, must be tested. Checking
only phi-bin centers can miss a negative region between centers. The exact
minimum is nonincreasing with epsilon for fixed U/LT/TT. In the interior branch
(`TT>0`) it is `U-epsilon*TT-(1+epsilon)*LT*LT/(4*TT)`; in the endpoint branch it
is `U+epsilon*TT-abs(LT)*sqrt(2*epsilon*(1+epsilon))`. The endpoint-minimum
condition makes the latter nonincreasing as well. Therefore the largest
observed epsilon is sufficient for the full lower-epsilon interval.

This mode changes the fit problem. It is not clipping a plotted curve or
replacing a negative coefficient after fitting. Finite-MC variance is evaluated
at the fitted coefficients during the existing iterative variance treatment.
The response normalization, data selection, migration matrix and signed
interference convention retain their existing meanings.

## Uncertainties and fit diagnostics

When a positivity constraint binds, the fitted estimator lies on a boundary.
The ordinary inverse-information matrix is not a confidence covariance for
that estimator. The Gaussian/asymptotic regularity assumptions and their
possible failure at parameter boundaries are discussed in the
[PDG 2025 Statistics review, Sections 40.2.2.1 and 40.4.2.3](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-statistics.pdf),
pages 5 and 29. The current mode therefore returns central values while
marking the ordinary statistical/MC covariance and derived statistical errors
unavailable: those numerical products contain NaN. Plots use a separately
named conditional refit-toy sampling SD when it can be computed. NaN means
unavailable, not zero uncertainty.

An unconstrained curvature matrix is saved separately for numerical diagnosis
in `migration_curvature_inverse.csv` and the ROOT matrix
`migration_unconstrained_curvature_inverse`.
It must not be relabeled as constrained covariance, used as a confidence
interval, or substituted for NaN errors in downstream fits. Constrained
profile intervals or calibrated repeated-sample refits would be needed for
boundary uncertainty inference. Target normalization remains a separate
correlated scale uncertainty; its existence does not supply the missing
statistical/MC intervals.

For plot diagnostics, 256 deterministic Gaussian pseudo-observations are
drawn around the constrained forward prediction with the converged row
variances. Every pseudo-observation repeats the constrained fit and finite-MC
variance iteration. The displayed bars are sample standard deviations for
the coefficients, forward yields, angular functions and residual-corrected
Eq. 5.30 points. They are not calibrated confidence intervals. The full
parameter and point toy covariances, including cross-bin terms, are in
`positivity_refit_toy_covariance.csv` and separately named ROOT matrices;
`fit_positive_refit_toys_successful` records the usable count. The original
inverse-information and Eq. 5.30 covariance products remain NaN at an active
boundary. These toys freeze the response, binning and input variance estimate,
using Gaussian row noise to represent the finite-MC term. They do not sample
event-level response changes, the target scale or model systematics.

The constraint and uncertainty records are:

- `positivity_diagnostics.csv` and the ROOT `positivity_diagnostics` tree:
  each active truth block's observed epsilon maximum, exact angular minimum
  of the response bracket, its minimizing cos(phi), and constraint status.
  The bracket is `2*pi*d2sigma/(dt dphi)`; it is not the angular cross section
  itself. Exported solver tolerances distinguish roundoff from a material
  violation.
- ROOT `fit_positive_xsec`, `fit_positivity_boundary_active` and
  `fit_positivity_iterations`: selected mode and final solver status.
- ROOT `migration_covariance_status` and `analysis_metadata`: interpretation
  of covariance, uncertainty availability and nominal degrees of freedom.
- Existing `migration_covariance.csv`, `migration_parameters.csv`, slice and
  phi summaries retain their usual names. Their statistical/MC covariance and
  derived error fields become NaN when the global fit has an active boundary.
  The separate target covariance remains a scale diagnostic.

The positive-mode `fit_scope` suffix is `_positive`: a full fit uses
`global_migration_positive`, and a recovered subset uses
`global_migration_subset_positive`. Existing unconstrained scope names remain.
Because coefficients are globally correlated, one active boundary suppresses
statistical/MC intervals for the whole fit, not just that block.

If positivity mode is enabled and the unconstrained solution already meets
the constraints, the ordinary conditional covariance can remain available;
consult the exported constraint and uncertainty status. In either mode the
existing finite-MC covariance approximations and unmodeled background,
acceptance and radiative systematics still apply.

The reported design rank and condition describe the measured response before
the constraints. They do not count independent physical directions after
applying a boundary. `rows - parameters` is only the nominal residual degree
count in a constrained fit, not a calibrated chi-square reference
distribution. Constraint activation is not evidence that background or
migration modeling has been corrected, and nonnegative curves alone do not
validate the extraction. Compare unconstrained and constrained results with
the same inputs, binning and cuts.

## Run and verify

From `/work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main`,
run the audited data/SIMC inputs with a separate positive-mode output directory:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash src/xsec_extract/run_xsec_pipeline.sh --kin KinC_x60_4b --target LH2 --root-dir output/KinC_x60_4b/root --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x60_4b_test/smearing_output/KinC_x60_4b/root/simc_pi0_analysis_output_smeared.root --vertex_simc_file output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x60_4b.root --mmiss-lower 0.6 --mmiss-upper 1.1 --positive-xsec --out-dir output/KinC_x60_4b/xsec_positive'
```

For an unconstrained comparison in its own directory:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash src/xsec_extract/run_xsec_pipeline.sh --kin KinC_x60_4b --target LH2 --root-dir output/KinC_x60_4b/root --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x60_4b_test/smearing_output/KinC_x60_4b/root/simc_pi0_analysis_output_smeared.root --vertex_simc_file output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x60_4b.root --mmiss-lower 0.6 --mmiss-upper 1.1 --no-positive-xsec --out-dir output/KinC_x60_4b/xsec_unconstrained'
```

These wrapper examples use the current configured binning; compare the saved
bin edges before comparing coefficients. Add `--partons` only to request the
optional native GK point projection. Positivity does not extend the GK model's
forward domain or convert that point projection into a bin-folded comparison.

Read `analysis_metadata` in the output ROOT file and the fit-status/parameter
tables before using a constrained result. Check mode, active constraints,
epsilon domain, solver and MC convergence, and uncertainty status. The combined
PDF retains the historical filename inside each selected output directory.
Do not interpret an absent error bar as an exact measurement or a toy SD as a
calibrated coverage interval.
