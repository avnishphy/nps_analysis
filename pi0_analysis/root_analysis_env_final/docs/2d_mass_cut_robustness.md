# Robust 2D Mass-Cut Design, Operation, and Validation

## 1. Purpose and scope

This document describes the two-dimensional cut built from
`mpi0_all` and `mmiss_all`, explains the observed failure in
`KinC_x25_1`, defines the revised algorithm, and gives reproducible validation
and tuning procedures.

Two implementations must remain behaviorally equivalent:

- `src/analysis/nps_2d_mass_cut.h`: C++/ROOT implementation used while each
  run is analyzed and by SIMC processing.
- `src/analysis/combine_analysis_branches.py`: NumPy implementation used after
  per-run trees are combined and normalized.

The combined diagnostic PDF is produced by the Python implementation, not by
the header. Therefore a header-only change cannot fix
`combined_branches_LH2_combined_2d_mass_cut_canvas.pdf`.

The cut adds selector branches. It does not delete events:

- per-run: `is_exclusive_ellipse`, `is_exclusive_mcd`;
- combined: `is_exclusive_ellipse_combined`,
  `is_exclusive_mcd_combined`.

## 2. Inputs and physics meaning

Each accepted input point is

```text
(x, y, w) = (mpi0_all, mmiss_all, positive event weight).
```

Per-run C++ weight is `pi0_weight`. Combined Python weight is
`pi0_weight * scale`; `scale` contains run normalization. Non-finite values,
non-positive weights, and points outside the configured histogram range are
ignored.

Default physical anchor:

```text
mpi0  = 0.13498 GeV
mmiss = 0.93800 GeV
```

These are an anchor and validation reference, not a forced fitted mean. Data
may move within configured tolerances.

`mmiss_all` is the uncorrected missing mass. Existing `is_exclusive` uses
`mmiss_all_corr`, so its blue contour in the combined plot is a useful
cross-check but is not mathematically identical to either new cut.

## 3. Measured original failure

The old `KinC_x25_1` combined debug output reported:

```text
seed peak       = (0.133875, 0.935250) GeV
core fraction   = 56.3494%
ellipse center  = (0.128883, 1.257484) GeV
MCD center      = (0.128463, 1.226266) GeV
ellipse fraction= 78.6536%
MCD fraction    = 32.1952%
```

The seed was physically sensible. Failure happened afterward: as density
threshold decreased, eight-neighbor flood-fill connected the seed component
to the long anti-correlated background ridge. The auto-selector chose a small
relative jump of only `1.086`, despite the resulting component containing more
than half the weighted histogram. Ordinary covariance and MCD then inherited
the contaminated location.

Additional per-run weakness existed in C++: all bins tied with the global
maximum were inserted as seeds. Sparse histograms often contain many tied
maxima, so disconnected background islands could enter the same logical core.

The former MCD target was `0.30 * total_weight`. That forced MCD to collect a
large amount of global background even when the physically relevant local
candidate population was much smaller.

## 4. Required invariants

The revised implementation enforces these invariants:

1. Peak discovery occurs only near the physical mass pair.
2. Exactly one peak bin seeds connectivity.
3. Smoothing may stabilize topology but may not alter fitted weights.
4. A core may broaden only while its weight and centroid remain safe.
5. MCD may use only a bounded candidate region around the physical anchor.
6. Covariance must remain invertible for sparse or collinear occupied bins.
7. Ellipse and MCD validity are independent; one failure must not erase the
   other valid selector.
8. Every safety decision must be visible in debug parameters.

## 5. Revised algorithm

### 5.1 Input filtering and histogram

Build a weighted `160 x 200` histogram over:

```text
0.11 <= mpi0_all  < 0.15 GeV
0.60 <= mmiss_all < 1.50 GeV
```

The raw histogram remains the source of all weights, means, covariances, and
reported fractions.

### 5.2 Density smoothing

Create a second, temporary density grid with a separable triangular kernel.
Default radius is one bin, equivalent to the normalized kernel

```text
1 2 1
2 4 2
1 2 1
```

At histogram boundaries normalization includes only in-range neighbors.
Smoothing is used only for peak choice, threshold comparison, and connected
component topology. It prevents one high-weight bin or a one-bin gap from
dominating the decision. It does not fabricate weight for covariance fitting.

### 5.3 Anchored, unique seed

Search smoothed density only inside:

```text
|mpi0  - 0.13498| <= 0.010 GeV
|mmiss - 0.93800| <= 0.030 GeV
```

Choose the maximum smoothed bin. If maxima tie numerically, choose the bin with
smallest normalized squared distance to the anchor. Start flood-fill from this
single bin only.

If the seed window contains no positive density, the result is invalid and all
selector flags remain zero. The code does not silently switch to a remote
global maximum.

### 5.4 Threshold scan and connected core

For threshold fraction `f`, accept eight-connected bins satisfying

```text
smoothed_density(bin) >= f * smoothed_density(seed).
```

Scan `f` downward from `0.60` to `0.05` in steps of `0.005`. For each connected
component, calculate from raw occupied-bin weights:

- occupied-bin count;
- weighted core weight and total-weight fraction;
- weighted `mpi0` and `mmiss` means;
- weight and bin-count growth from previous threshold.

A component is safe only when all conditions hold:

```text
occupied bins >= 3
core weight > 0
core fraction <= 0.30
|mean mpi0  - seed mpi0 | <= 0.008 GeV
|mean mmiss - seed mmiss| <= 0.100 GeV
```

A safe component is qualified when it also has at least eight occupied bins
and at least 0.5% of total weight. The broadest qualified component is chosen.
If statistics never reach qualification, the broadest safe component is used.
If a formerly safe component violates the fraction or centroid bounds, the
scan stops and `core_leak_rejected=1` is recorded.

For an explicitly configured positive `peak_fraction`, the same safety checks
still apply. Manual configuration cannot bypass anchoring and leak protection.

### 5.5 Regularized covariance

For weighted points `z_i=(x_i,y_i)` with weights `w_i`, location and covariance
are

```text
mu = sum(w_i z_i) / sum(w_i)
C  = sum(w_i (z_i-mu)(z_i-mu)^T) / sum(w_i).
```

Sparse histograms can give zero variance in one direction. The revised code
adds a diagonal regularizer equal to half a histogram bin in each coordinate,
squared. Correlation is capped at `|rho|=0.995`. This is small compared with a
resolved peak but prevents singular inversion caused only by binning.

The squared Mahalanobis distance is

```text
d2 = (z-mu)^T C^-1 (z-mu).
```

The density-connected component is only an anchored seed. It is intentionally
small because the physical exclusive population forms a continuous
anti-correlated band rather than an isolated density island.

### 5.6 Iterative ellipse-fit-subset growth

Starting from the seed covariance, the ellipse fitting subset grows by iterative
Mahalanobis inclusion. At each iteration, local candidate bins satisfying

```text
d2 <= chi2_quantile_df2(0.990)
```

are retained from the full analysis histogram and their raw-weight covariance
is recomputed. No rectangular candidate gate is applied, avoiding an abrupt fit
edges. Growth stops at convergence or after 20 iterations. Every update must
remain inside the same centroid tolerances and the grown subset may not exceed
30% of total weight. Thus the fit can follow the diagonal event band but
cannot flood into the remote global background.

This subset is an estimator-training detail. Its boundary is determined partly
by the 30% contamination guard and therefore is not displayed as a physics
core.

### 5.7 Reported core contour

After the ordinary covariance is fixed, the reported red core is recomputed as
the bins satisfying

```text
d^2 <= chi2_quantile_df2(core_quantile)
```

with `core_quantile=0.6827`. It is therefore a smooth inner Mahalanobis contour,
strictly nested in and aligned with the final 99% ellipse. The core does not
participate in fitting or event-branch selection; changing `core_quantile`
changes only this diagnostic definition. The former grown population remains
available as `h_ellipse_fit_subset` and through `fit_subset_*` parameters.
Final ordinary ellipse uses the 99% two-dimensional chi-square contour with
`1.05^2` padding.

### 5.8 Local MCD candidate population

MCD may examine occupied bins only inside:

```text
|mpi0  - 0.13498| <= 0.015 GeV
|mmiss - 0.93800| <= 0.180 GeV
```

Each C-step sorts those candidates by current `d2` and retains 50% of local
candidate weight. It recomputes covariance and repeats up to eight times. A
candidate covariance is accepted only if its mean satisfies the same anchor
tolerances as the ordinary model.

This differs critically from retaining 30% of the entire analysis histogram.
Remote background cannot satisfy the local target merely by carrying large
global weight.

### 5.9 MCD consistency correction and contour

Covariance calculated from the central fraction of a Gaussian is biased low.
For two dimensions and retained fraction `h`, the implementation uses

```text
c       = -2 ln(1-h)
v_trunc = 1 - c(1-h)/(2h)
C_full  = C_trimmed / v_trunc.
```

The final MCD contour uses the exact two-degree-of-freedom chi-square quantile

```text
chi2_quantile(q) = -2 ln(1-q)
```

with `q=0.975`, followed by `1.05^2` padding. This gives MCD a defined
probabilistic target instead of taking a 99.5th percentile from an already
truncated subset.

### 5.10 Independent validity

`Params` now contains:

```text
valid
ellipse_valid
mcd_valid
core_leak_rejected
```

`valid` means at least the anchored ordinary model succeeded. MCD may fail
independently. In that case ellipse flags are still applied and MCD flags stay
zero. `apply_mass_cuts` follows these per-model validity bits.

## 6. Default configuration reference

| Parameter | Default | Meaning |
|---|---:|---|
| `seed_mpi0` | 0.13498 GeV | physical pi0 anchor |
| `seed_mpi0_half_width` | 0.010 GeV | peak-search x half-width |
| `seed_mmiss` | 0.938 GeV | proton missing-mass anchor |
| `seed_mmiss_half_width` | 0.030 GeV | peak-search y half-width |
| `smoothing_radius_bins` | 1 | triangular smoothing radius |
| `core_quantile` | 0.6827 | reported inner-core chi-square probability |
| `ellipse_quantile` | 0.990 | fit growth/final ellipse chi-square probability |
| `ellipse_grow_iterations` | 20 | maximum ordinary-core growth steps |
| `auto_min_core_total_fraction` | 0.005 | qualification floor |
| `auto_max_core_total_fraction` | 0.30 | leakage ceiling |
| `auto_min_core_bins` | 8 | qualification occupancy |
| `max_model_mpi0_offset` | 0.008 GeV | fitted-center x tolerance |
| `max_model_mmiss_offset` | 0.100 GeV | fitted-center y tolerance |
| `mcd_candidate_mpi0_half_width` | 0.015 GeV | MCD local x gate |
| `mcd_candidate_mmiss_half_width` | 0.180 GeV | MCD local y gate |
| `mcd_keep_candidate_fraction` | 0.50 | local trimmed-weight fraction |
| `mcd_iterations` | 8 | maximum C-steps |
| `mcd_ellipse_quantile` | 0.975 | final chi-square probability |
| `covariance_regularization_bins` | 0.50 | diagonal floor in bin units |

Defaults are identical in C++ and Python. When changing one, change the other
and rerun both regression tests.

## 7. Reading diagnostics

Important fields in `*_params.csv` and combined `*_debug.txt`:

- `peak_mpi0`, `peak_mmiss`: anchored density maximum;
- `peak_weight`: raw content of peak bin;
- `smoothed_peak_weight`: density used by threshold scan;
- `peak_fraction`: chosen density threshold fraction;
- `core_total_fraction`: raw weight inside the final 68.27% inner contour divided by total weight;
- `fit_subset_*`: iteratively grown population used to estimate the ordinary covariance;
- `seed_core_*`: density-connected seed statistics before Mahalanobis growth;
- `ellipse_growth_steps`: accepted ordinary-core growth iterations;
- `core_leak_rejected`: later threshold crossed a safety bound;
- `mean_*`: ordinary ellipse center;
- `mcd_mean_*`: robust center after consistency scaling;
- `ellipse_valid`, `mcd_valid`: whether each selector is trustworthy;
- `ellipse_total_fraction`, `mcd_total_fraction`: final accepted weighted
  fractions over the configured analysis range;
- `mcd_candidate_weight`, `mcd_candidate_bins`: size of local MCD population.

Red core bins in the canvas must form a smooth contour nested inside the
magenta ellipse. They are not the covariance-training subset.
Neither fitted center should violate configured anchor tolerances. A
`core_leak_rejected=1` value is expected when the scan successfully detects the
background connection; it is not itself a fit failure.

## 8. Validation performed for this change

Run from repository root:

```bash
cd /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main
```

### 8.1 Python regression tests

Actually run:

```bash
MPLCONFIGDIR=/tmp/nps-mpl-cache \
PYTHONPYCACHEPREFIX=/tmp/nps-pycache \
python3 -m unittest -v tests/test_2d_mass_cut.py
```

Expected result:

```text
Ran 2 tests
OK
```

The main synthetic test contains a physical anti-correlated signal plus a
larger, stronger remote ridge. Measured Python result:

```text
core fraction                = 18.81%
ellipse-fit-subset fraction  = 28.92%
signal ellipse efficiency    = 99.65%
remote-background acceptance = 14.80%
```

The second test verifies covariance remains invertible for collinear occupied
bins once bin-scale regularization is applied.

### 8.2 C++ regression test

Actually run:

```bash
csh -c 'source /group/nps/singhav/setup.csh; \
  g++ -std=c++14 `root-config --cflags` tests/test_nps_2d_mass_cut.cpp \
      `root-config --libs` -o /tmp/test_nps_2d_mass_cut; \
  /tmp/test_nps_2d_mass_cut'
```

Expected result starts with `PASS`. Measured result:

```text
peak=(0.134625,0.95325)
core=0.187498
fit_subset=0.287529
mean=(0.134975,0.937992)
mcd_valid=1
```

The command sources the required Hall C environment, compiles only the small
test executable, writes that executable under `/tmp`, and does not modify
analysis outputs.

### 8.3 Existing `KinC_x25_1` histogram

The stored combined histogram was read without modifying its ROOT file. Each
occupied bin center and content was passed through the revised Python fitter.
Measured revised result:

```text
peak                  = (0.135875, 0.957750) GeV
core fraction         = 18.9156%
ellipse-fit subset    = 28.7072%
ordinary model center = (0.131086, 1.024421) GeV
core leak rejected    = 1
ellipse fraction      = 43.0461%
MCD center            = (0.130626, 1.035787) GeV
MCD fraction          = 42.3050%
existing de-corr frac = 40.6552%
```

The 68.27% core is now a smooth inner contour of the ordinary model. The
28.7072% grown subset still estimates that model but is no longer mislabeled or
drawn as the core. Ordinary ellipse and MCD acceptance are unchanged by this
semantic/display correction.

This histogram replay reproduces fit geometry and weighted fractions at bin
resolution. It does not replace a complete event-level regeneration of the
combined tree.

The same stored histogram was also replayed through C++. C++ and Python agreed
at printed precision on peak, threshold, model center, fit-subset/core
fractions, ellipse fraction, MCD center, and MCD fraction. This directly checks
behavioral parity for the reported failure dataset.

## 9. Regenerating outputs

The following command reruns the combined stage and overwrites its configured
combined ROOT/debug outputs. Run it only when regeneration is intended:

```bash
cd /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main
MPLCONFIGDIR=/tmp/nps-mpl-cache \
PYTHONPYCACHEPREFIX=/tmp/nps-pycache \
python3 src/analysis/combine_analysis_branches.py \
  --kin KinC_x25_1 \
  --target LH2 \
  --no-analysis-plots
```

Key outputs:

```text
output/KinC_x25_1/root/combined_branches_LH2.root
output/KinC_x25_1/root/combined_branches_LH2_combined_2d_mass_cut_debug.txt
output/KinC_x25_1/root/combined_branches_LH2_combined_2d_mass_cut_canvas.pdf
output/KinC_x25_1/root/combined_branches_LH2_combined_2d_mass_cut_canvas.png
```

Per-run products require rerunning the main analysis for those runs because
the flags are generated inside `nps_analysis_main.C`. Follow repository run
wrappers and use a one-run list before launching a full production run.

## 10. Cross-kinematics tuning protocol

Do not tune from only the final accepted fraction. For every kinematic setting:

1. Confirm seed peak lies inside the intended physical window.
2. Confirm red core bins form a smooth contour nested in the magenta ellipse.
3. Confirm ordinary and MCD centers satisfy anchor limits.
4. Compare accepted and rejected 1D projections.
5. Compare with `mmiss_all_corr`/legacy selection as an independent diagnostic,
   not as exact ground truth.
6. Inspect per-run stability; combined success can hide a weak individual run.
7. Compare data and SIMC using the same configuration.
8. Record efficiency/background tradeoffs before changing a default.
9. Rerun Python and C++ regression tests after every parameter or algorithm
   change.

If valid signal calibration genuinely moves outside default anchor tolerances,
configure anchors per kinematic setting rather than widening all bounds. Broad
global windows recreate the original failure mode.

## 11. Failure behavior and troubleshooting

### No positive density in seed window

Meaning: no usable weighted events near physical anchor, incorrect units, bad
calibration, or wrong kinematic input. Result is invalid; selectors stay zero.

### No safe anchored core

Meaning: seed component has too few occupied bins or immediately violates
fraction/centroid bounds. Inspect statistics and binning. Do not automatically
fall back to global maximum.

### `ellipse_valid=1`, `mcd_valid=0`

Meaning: ordinary core succeeded but local MCD lacked enough candidates or an
iteration escaped anchor limits. Ellipse remains usable; MCD flags stay zero.

### `core_leak_rejected=1`

Meaning: descending threshold found the connection to background and retained
the preceding safe component. Usually desirable. Investigate only if chosen
core is visibly too small.

### Fit changes strongly with one histogram bin

Check effective statistics, unusually large normalization weights, and whether
per-run results should be combined before fitting. Increase smoothing only
after validating that it does not merge physical and background components.

## 12. Maintenance rules

- Keep C++ and Python defaults synchronized.
- Preserve existing public function names and output branch names.
- Keep smoothing separate from raw fit weights.
- Never remove anchoring/leak checks to make a difficult dataset return valid.
- Add a regression dataset for every newly observed failure mode.
- Validate on smallest representative input before production regeneration.
