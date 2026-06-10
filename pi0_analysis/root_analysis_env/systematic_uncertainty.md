# Systematic Uncertainty Notes

This file collects uncertainty sources and implementation details for the pi0
analysis. Add future sources as separate sections with the same structure:

- physics or reconstruction question being tested
- nominal implementation
- stored parameters needed to reproduce the cut/correction
- variation knobs
- outputs to compare
- recommended uncertainty prescription

## 2D Mass Cut: `mmiss_all` vs `mpi0_all`

### Motivation

The current missing-mass correction in `src/nps_mmiss_cor.h` intentionally
removes the correlation between proton missing mass and pi0 invariant mass.
As an alternate exclusivity selection, `scripts/test_2D_mass_cut.C` tests a
direct 2D cut on the weighted `mmiss_all:mpi0_all` distribution.

The input test file is:

```text
output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root
```

The tree is:

```text
physics
```

The variables are:

```text
x = mpi0_all
y = mmiss_all
w = pi0_weight
```

### Current Test Macro

Main implementation:

```text
scripts/test_2D_mass_cut.C
```

ROOT environment:

```csh
source /group/nps/singhav/setup.csh
```

In non-interactive `csh`/`tcsh` sessions, source module init first if needed:

```csh
source /usr/share/Modules/init/csh
source /group/nps/singhav/setup.csh
```

Run:

```bash
root -l -b -q 'scripts/test_2D_mass_cut.C()'
```

### Auto Peak-Core Selection

The macro first builds a connected high-density core around the highest-weight
2D bin. A bin belongs to the core if:

```text
bin_weight >= peak_fraction * peak_bin_weight
```

and if it is connected to the peak bin through 8-neighbor connectivity.

By default:

```text
peak_fraction = -1.0
```

means auto mode. Auto mode scans `peak_fraction` from high to low:

```text
auto_peak_max = 0.60
auto_peak_min = 0.05
auto_peak_step = 0.005
```

For each scanned threshold, the macro records:

```text
selected_bins
selected_weight
selected_total_fraction
weight_ratio_to_previous
bin_ratio_to_previous
```

If the connected component grows suddenly, that is interpreted as a sideband
joining the signal component. The selected `peak_fraction` is the previous
threshold before the largest growth jump.

For the current test file, the sideband join is visible around:

```text
peak_fraction  selected_bins  selected_weight
0.130          1465           52281.1
0.125          4878           128131
```

The auto-selected nominal value is:

```text
peak_fraction = 0.13
auto_jump_ratio = 3.32969
```

### Standard Mahalanobis Ellipse

After the connected core is selected, the macro computes the weighted mean and
covariance from the core bins:

```text
mu = (mean_mpi0, mean_mmiss)
C  = [[cov_mpi0_mpi0, cov_mpi0_mmiss],
      [cov_mpi0_mmiss, cov_mmiss_mmiss]]
```

The squared Mahalanobis distance is:

```text
d2 = [x - mu]^T C^-1 [x - mu]
```

The ellipse cut is:

```text
d2 <= ellipse_d2_cut
```

The current macro chooses `ellipse_d2_cut` from the weighted quantile of core
bin distances, then pads the radius:

```text
ellipse_d2_cut = quantile(core_d2, ellipse_core_quantile) * ellipse_padding^2
```

Nominal test settings:

```text
ellipse_core_quantile = 0.995
ellipse_padding = 1.05
```

### Minimum Covariance Distance Method

The macro also tests a weighted minimum-covariance-distance (MCD) style robust
ellipse. The current implementation is a histogram-bin MCD approximation:

1. Start from the auto-selected core covariance.
2. Compute Mahalanobis distance for all positive-weight bins.
3. Keep the closest fixed fraction of total weighted content.
4. Recompute weighted mean and covariance.
5. Repeat for `mcd_iterations`.
6. Keep the covariance model with the smallest determinant.
7. Build a final MCD ellipse from a weighted distance quantile plus padding.

Nominal test settings:

```text
mcd_keep_total_fraction = 0.30
mcd_iterations = 8
mcd_ellipse_quantile = 0.995
mcd_padding = 1.05
```

The MCD cut equation is the same as the standard Mahalanobis ellipse, but with
robust MCD parameters:

```text
d2_mcd = [x - mu_mcd]^T C_mcd^-1 [x - mu_mcd]
accept if d2_mcd <= mcd_d2_cut
```

### Parameters To Save For Systematics

For the standard Mahalanobis ellipse, save:

```text
mean_mpi0
mean_mmiss
cov_mpi0_mpi0
cov_mpi0_mmiss
cov_mmiss_mmiss
ellipse_d2_cut
ellipse_axis1
ellipse_axis2
ellipse_angle_deg
peak_fraction
auto_peak_fraction
auto_peak_min
auto_peak_max
auto_peak_step
auto_jump_ratio
ellipse_core_quantile
ellipse_padding
```

For the MCD ellipse, save:

```text
mcd_mean_mpi0
mcd_mean_mmiss
mcd_cov_mpi0_mpi0
mcd_cov_mpi0_mmiss
mcd_cov_mmiss_mmiss
mcd_det
mcd_d2_cut
mcd_keep_total_fraction
mcd_target_weight
mcd_subset_weight
mcd_subset_bins
mcd_iterations
mcd_ellipse_quantile
mcd_padding
```

For all variations, also record:

```text
accepted weighted yield
accepted unweighted event count
mpi0_all projection
mmiss_all projection
Q2, xB, t, phi distributions after cut
final cross section or asymmetry observable
```

### Current Output Files

The test macro writes:

```text
output/plots/x60_4b/production_wfpi0/test_2D_mass_cut.root
output/plots/x60_4b/production_wfpi0/test_2D_mass_cut.pdf
output/plots/x60_4b/production_wfpi0/test_2D_mass_cut.png
output/plots/x60_4b/production_wfpi0/test_2D_mass_cut_mask_bins.csv
output/plots/x60_4b/production_wfpi0/test_2D_mass_cut_ellipse_bins.csv
output/plots/x60_4b/production_wfpi0/test_2D_mass_cut_mcd_bins.csv
output/plots/x60_4b/production_wfpi0/test_2D_mass_cut_ellipse_params.csv
output/plots/x60_4b/production_wfpi0/test_2D_mass_cut_peak_scan.csv
```

The parameter CSV is the main reproducibility file:

```text
test_2D_mass_cut_ellipse_params.csv
```

The scan CSV diagnoses the auto `peak_fraction` choice:

```text
test_2D_mass_cut_peak_scan.csv
```

### Suggested Variation Knobs

Standard Mahalanobis ellipse:

```text
ellipse_padding = 0.95, 1.00, 1.05, 1.10
ellipse_core_quantile = 0.990, 0.995, 0.998
peak_fraction = auto, 0.12, 0.13, 0.14
```

MCD ellipse:

```text
mcd_keep_total_fraction = 0.25, 0.30, 0.35
mcd_padding = 0.95, 1.00, 1.05, 1.10
mcd_ellipse_quantile = 0.990, 0.995, 0.998
mcd_iterations = 5, 8, 12
```

Histogram/binning sensitivity:

```text
n_mpi0_bins = 120, 160, 200
n_mmiss_bins = 150, 200, 250
mpi0 range variations around 0.11-0.15 GeV
mmiss range variations around 0.6-1.5 GeV
```

### Recommended Uncertainty Evaluation

Use one nominal cut definition for the production analysis. Then rerun the
final analysis with a controlled set of variations. For each final observable,
compare variation results to nominal:

```text
delta_i = observable_i - observable_nominal
```

Possible prescriptions:

```text
envelope = max(abs(delta_i))
RMS      = sqrt(mean(delta_i^2))
```

Use the prescription that matches the collaboration/analysis convention. Keep
standard Mahalanobis and MCD variations separate at first; combine only after
checking whether they represent the same cut-shape uncertainty or distinct
model choices.

### Open Items

- Decide which cut becomes nominal in production: standard covariance ellipse
  or MCD ellipse.
- Add event-level branch output once this cut enters `src/nps_analysis*.C`.
- Add per-run or per-kinematics stability checks for saved ellipse parameters.
- Study sensitivity to `pi0_weight` definition and possible `scale` weighting.
- Compare final cross-section impact against existing `mmiss_all_corr` method.

## Future Source: Placeholder

Use this template for later uncertainty sources.

### Motivation

Describe physics/reconstruction reason.

### Nominal Implementation

Describe nominal code path and parameters.

### Parameters To Save

List reproducibility quantities.

### Variation Knobs

List controlled changes.

### Outputs To Compare

List plots, yields, and final observables.

### Recommended Prescription

State envelope/RMS/model-comparison approach.
