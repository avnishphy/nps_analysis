# Neutral-pion cross-section extraction

Two extraction methods are available. `excl_xsec_pi0_analysis_no_simc_model.C`
steers the response-fit method described below. The iterative data/SIMC ratio method is implemented by
`excl_xsec_pi0_analysis_simc_model.C` with the parameterized Fortran pi0
model in `simc_pi0_reweight.h`. See
[ITERATIVE_SIMC_MODEL.md](ITERATIVE_SIMC_MODEL.md) for the equations,
physical-t and phi-center conventions, fit diagnostics, limitations, and
reproducible validation. Physics defaults come from a JSON preset in
`xsec_config/`. The wrapper generates `xsec_config.h` in its temporary build
directory and uses that header for either extractor.

For independent reconstructed/truth/publication binning, finite-MC event
resampling and paired binning comparisons, see
[Independent-grid forward extraction](FORWARD_EXTRACTION.md). This explicit
mode includes x36_4 validation results and the remaining publication checks.

Use `--xsec_config xsec_config_x36_4.json` (or a path to a JSON file in
`src/xsec_extract/xsec_config/`) to select a preset. If the flag is omitted,
`--kin KinC_x36_4` selects `xsec_config_x36_4.json` when that file exists.
For a new setting, copy a JSON preset to `xsec_config_<setting>.json`, edit
its values, then pass it with `--xsec_config`. The JSON file contains
paths, beam energy, bin edges, cuts, and fit defaults; the generator rejects
unknown keys, nonfinite numbers, and inconsistent Q2/xB bins. An unmatched
`--kin` needs an explicit preset. `--kin` still selects paths only; it never changes the
physics values in a selected JSON file. Explicit CLI values override
non-binning defaults. The extractor warns if the selected preset's
`configured_kinematic` differs from `--kin`.

## Bin edges

Set `phi_bins` (uniform 0 to 2pi), or replace it with a numeric
`phi_bin_edges` array. Set `tprime_bin_edges`, `t_bin_edges`,
`q2_bin_edges`, and `xb_bin_edges_by_q2` in the JSON preset. Supply one
xB edge row per Q2 interval; each row must have the same count and outer
limits. Bin counts and
selection bounds follow the vectors automatically. The response fit uses
tprime; the SIMC-model ratio fit uses physical t. Both use the same Q2, xB,
and phi edges. Invalid or unordered vectors stop extraction. CLI bin-count
and range options have been removed so these vectors remain authoritative.
Changing these defaults changes the extracted physics bins; rerun both methods
and compare their saved edge metadata after editing the JSON preset.

Set `diamond_xb_q2_vertices` to `null` for no polygon cut, or four
`[xB, Q2]` corners copied from the bin-selection notebook. The config
generator checks physical, distinct, convex corners and orders them
counterclockwise. Both extractors intersect this polygon with the bin-edge
limits for reconstructed data and SIMC, including their per-slice missing-mass
diagnostics. The W/Q2 boundary is derived from these xB/Q2 corners; no separate
W vertices are supplied. The selected vertices are recorded in ROOT metadata.
The two existing presets default to `null`, preserving their previous selection.

The SIMC-model method extracts U, LT and TT. It reports LT' as unavailable in
ROOT metadata and the slice CSV because helicity-specific luminosities and beam
polarization are absent from its input contract. It does not fit or publish a
helicity-difference amplitude as LT'.

The response-fit method orders input checks, frozen binning, response
accumulation, the global fit, optional GK comparison, plotting, and
serialization. Its headers are implementation modules, not separately compiled
libraries.

## Organization

| Files | Responsibility |
| --- | --- |
| `excl_xsec_pi0_analysis_simc_model.C`, `simc_pi0_reweight.h` | Iterative data/SIMC ratio extraction and parameterized Fortran pi0 model |
| `xsec_root.h`, `xsec_config_template.h.in`, `generate_xsec_config.py`, `xsec_types.h`, `xsec_analysis.h` | Dependencies, JSON-to-header configuration, records, analysis state/interfaces |
| `xsec_cli.h`, `xsec_input.h` | Options, mandatory branches, matched generated-event lookup |
| `xsec_physics.h`, `xsec_binning.h` | Kinematic conventions and frozen bin boundaries |
| `xsec_accumulation.h` | Selected data/MC event loops and diagnostic histograms |
| `xsec_response.h` | Generated-to-reconstructed response, exterior regions, MC moments |
| `xsec_linear_solver.h`, `xsec_fit.h` | Full-rank SVD, global fit, finite-MC uncertainty |
| `pi0_response_conventions.h`, `partons_pi0_projection.h`, `xsec_model.h` | Hand-flux/Fourier conventions and optional native GK comparison |
| `xsec_plot_style.h`, `xsec_plot_global.h`, `xsec_plot_slices.h` | Global and fitted-slice displays |
| `xsec_plot_migration.h`, `xsec_plot_migration_coverage.h` | Component response, selected-MC fractions and vertex/reco coverage |
| `xsec_output.h`, `xsec_migration_output.h` | ROOT/CSV records, complete covariance, reproducibility metadata |

## Reconstructed mass selection

The response-fit extractor defaults to the shared 0.80-1.10 GeV missing-mass
window. Use `--mmiss_select mcd` or `--mmiss_select ellipse` to
select data with `is_exclusive_mcd_combined` or
`is_exclusive_ellipse_combined`, respectively. Both modes apply the exact
combined-data covariance shape to reconstructed SIMC `mpi0` and `mmiss`;
SIMC remains restricted to the generated-exclusive channel.

Combined-data production now writes covariance and fit-domain fields to
`combined_branches_<target>_combined_2d_mass_cut_debug.txt`. For older data
files, the wrapper exports those fields to the xsec output directory after
recomputing the cut and checking both stored flags event by event. Direct
extractor runs can use `--mmiss-cut-file <path>`; to export metadata manually:

```bash
python3 src/analysis/export_combined_mass_cut_metadata.py \
  output/<kin>/root/combined_branches_LH2.root \
  --out /tmp/combined_mass_cut_debug.txt
```

The `mmiss_selection_2d.pdf/png` figure shows all and selected data events
beside all and selected generated-exclusive SIMC events, with the chosen
geometric boundary when available. The same four histograms are saved in
the output ROOT file. The next PDF page, `mmiss_selection_1d.pdf/png`,
keeps corrected data Mx (`mmiss_all_corr`) alone, with full and peak-zoom
overlays of all and selected events. The following page,
`mmiss_reconstructed_comparison.pdf/png`, shows the same two views for
reconstructed data `mmiss_all` and SIMC `mmiss` in aligned rows. SIMC's
full sample is generated-exclusive. These six unweighted histograms are also
saved in the output ROOT file. `--no-diagnostics` omits these figures.

## Plot output and switches

The wrapper and compiled extractor accept `--no-pdf`, `--no-png` and
`--no-diagnostics`. The first two select image formats independently; with
both disabled no plot files or plot subdirectories are created. The diagnostic
switch omits global QA, missing-mass, epsilon and migration figures while
retaining slice fits, coefficient figures and optional `--partons` comparisons.
Numeric ROOT/CSV fit products remain available in every mode. `--partons`
still computes its comparison points when graphics are disabled.

Individual PDF and PNG filenames are retained. The combined PDF is assembled
from the completed individual PDFs with `pdfunite` when PDF output is enabled;
the Hall C analysis host provides this utility. This avoids ROOT warnings from
interleaved PDF streams and gives the combined file the same page order as the
individual plots. No combined PDF is written for PNG-only or graphics-off runs.

Slice angular and coefficient figures display nb/GeV^2 (and nb/(GeV^2 rad)
for the angular ordinate), matching the PARTONS comparison figures. Diagnostic
response heatmaps display weighted yield per nb/GeV^2 coefficient; their ROOT
matrices and CSV retain raw SIMC response units. These are display conversions
only. Coefficient and PARTONS markers sit at accepted vertex means of -t',
with horizontal whiskers showing the generated-bin extent (the U term alone
shows the shared bin span on the three-term figure).
Undefined reconstructed-yield ratios are omitted as points and counted in the
panel annotation. They remain NaN in CSV/ROOT. Ratio error bars use data sumw2
only and are same-data residual diagnostics. Experimental points, fitted
curves, target-systematic separation and positivity-boundary interval status
retain their explicit legends and output meanings. Global data/SIMC shape
plots use display-only area scaling, never a fit normalization.

## Physics convention and normalization

The fitted unpolarized virtual-photon differential cross section is

```text
d2sigma/(dt dphi) = [U + sqrt(2 epsilon (1+epsilon)) LT cos(phi)
                        + epsilon TT cos(2 phi)] / (2 pi).
```

`U` denotes `sigma_T + epsilon sigma_L`; a single setting does not separate
T and L. The fit treats U, LT and TT as constants within each generated bin.
The legacy name `TL` in CSV fields means LT. Coefficients retain SIMC `sigcm`
units, microbarn/MeV^2; plots multiply by 1e9 for nb/GeV^2. Phi is in radians.

Data retain their existing efficiency/charge weights and are divided by
`tgt_contam`, default `0.584 +/- 0.014`. Its common uncertainty is a separate
rank-one covariance, not independent noise in each bin. SIMC uses its physical
normalization through `full_weight / sigcm`. This quotient already retains
generation volume, luminosity convention, flux, Jacobian and simulated radiation.
Do not multiply another flux or divide by reconstructed bin widths.

Both `--normalize_mmiss` and `--normalize-mmiss` now fail explicitly. The fit
never matches the SIMC integral to measured data. Area matching remains only
in copies used for missing-mass shape diagnostics and never changes the response.

GK's native electron observable is divided by the **full Hand flux Gamma**.
The GK06 process prefactor is Gamma/(2 pi); confusing these introduced the old
extra factor of 2 pi. Fourier projection then recovers U/LT/TT in the convention
above. Corrected GK coefficients are the old adapter's coefficients divided by
2 pi at identical kinematics/integration settings. GK points use response-weighted
generated-bin means and are comparisons, not bin-centering corrections or inputs
to the fit. Native integration precision and the experiment's LT azimuth-sign
convention need their own checks.

## Migration response and fit

### Defurne Eqs. 5.12-5.31 implementation (2026-09-24)

Indices are explicit: reconstructed row `r=(it,iq,ix,ip_reco)`;
truth block `b=(it,iq,ix)` or one of six exterior regions; truth angular
cell `v=(b,ip_vertex)`; component `n=(U,LT,TT)`. The fitted design is
`D[r,b,n]=sum_selected_in_(r,b) (full_weight/sigcm)*basis_n(phi_vertex,epsilon_vertex)`.
`truth_phi_response_cells.csv` retains the finer `D[r,v,n]` whose sum over
vertex phi gives the original design. All `v` in one `b` share three fitted
coefficients. The corresponding pair for experimental points uses identical
grid indices at reconstruction and vertex: `r=(b,ip)`, `v=(b,ip)`.
Migration from other truth cells, including exterior regions, remains in
`mu_r` and is subtracted. A missing diagonal or denominator small relative
to the row's absolute contributions is explicitly unavailable.

The basis `[1,sqrt(2 epsilon (1+epsilon)) cos(phi),epsilon cos(2phi)]/(2pi)`
implements the helicity-even virtual-photon structure functions of Eq. 5.19.
`U=T+epsilon L` is fitted as one bin-constant coefficient while epsilon
varies event by event in LT/TT. This approximates the true within-bin
`T+epsilon L`; it is not a Rosenbluth separation. The inputs do not justify
an LT' fit. The displayed `PhiBin::xsec` is the historical fitted function
at the phi-bin center using the truth-block mean epsilon. It is not an
experimental point or a phi-bin average.

The new `experimental_points.csv` and ROOT tree use the vertex-phi cell's
accepted-response-weighted Q2/xB/t' and epsilon, and its phi-bin center, as
the reference point. `sigma_fit_reference` evaluates the fitted angular
function there. For `M[r,v]=sum_n D[r,v,n] X[b(v),n]` and
`mu_r=sum_v M[r,v]`, the exported `subtracted_yield=y_r-(mu_r-M[r,v])`,
`correction=subtracted_yield/M[r,v]`, and
`sigma_exp=correction*sigma_fit_reference`. Exact forward closure gives
`correction=1`; a no-migration row reduces to `y_r/mu_r`. These points
remain model-dependent through the bin-constant fitted angular response,
but receive no GK/PARTONS bin-centering factor and no eta correction
(`eta=1`). Separate columns compare observed reconstructed Q2/xB/t'
means with the forward-predicted reconstructed means, using the same selected
events and fitted basis weights. Truth reference means are separately named
because data do not directly measure vertex coordinates.

`--fit-variance data` selects the fixed-variance Eq. 5.23 reference:
`V_data,rr=sum_selected_data w_event^2/tgt_contam^2`, with zero-variance
rows excluded and recorded. The default `--fit-variance finite-mc` adds
iterated Poissonized MC outer-product variance as an extension. The wrapper
accepts the same option; `NPS_XSEC_FIT_VARIANCE` supplies its default.
Both export exact last-solve variances, fit rows, convergence and support.
These data variances are conditional weighted-yield variances: they do not
include uncertainty of the upstream `pi0_weight` model. With fixed
variances, `A=D'V^-1D`, `B=D'V^-1y`, `AX=B` and `C=A^-1`;
the implementation uses a full-rank scaled SVD for the solve and covariance.
`migration_correlation.csv` and ROOT
`migration_parameter_correlation` contain
`rho_ij=C_ij/sqrt(C_ii*C_jj)`, never a normalized information matrix.
The independent target-divisor uncertainty stays in separate rank-one
covariance products. A positivity boundary leaves confidence covariance and
correlation unavailable; inverse curvature remains a distinct diagnostic.
For display, a boundary fit draws 256 deterministic Gaussian weighted-yield
toys around its fitted forward prediction, using the converged data-plus-MC
row variance. Each toy repeats the constrained fit and finite-MC variance
iteration. Vertical bars show the sample standard deviations of fitted
coefficients, forward yields, angular functions and Eq. 5.30 points. These
are conditional sampling spreads, not 68% coverage intervals.
`positivity_refit_toy_covariance.csv` and separately named ROOT covariance
and correlation matrices retain cross-bin correlations. Ordinary covariance
and point-error fields remain NaN on a boundary. The toy calculation freezes
response moments, bin edges, target correction and the final input variance
estimate; MC response fluctuations are approximated as Gaussian row noise.
Event-level MC response changes and model/target systematics need separate
study. No toy result changes the fitted central values.

For experimental-point uncertainty, put `z=y_r-D_r X`, `f=g_v X`,
`m=d_rv X`, and `E=f(1+z/m)`. At fixed response and row variances,
`K=C D' V^-1` and
`h=(1+z/m)g_v-(f/m)D_r-(fz/m^2)d_rv`, giving
`dE/dy_s=(f/m)delta_rs+h'K_s`. The exported full point covariance
contracts these gradients with the same data's `V_data`; it therefore keeps
numerator, fit and denominator correlation and cross-bin terms. In finite-MC
mode it additionally differentiates the same fitted estimator with respect
to every vertex-phi response cell and contracts those derivatives with its
stored 3x3 event outer products. This is a first-order Poissonized-MC
extension conditional on final GLS weights and fixed configured bin edges.
The reference epsilon's MC ratio uncertainty is included using its
cross moments with the same vertex-phi response events.
It does not propagate changes in the iterated variance estimator, fixed-Ngen
multinomial correlations, upstream signal-weight uncertainty, or detector
and radiation systematics. `experimental_point_covariance.csv`,
`experimental_point_correlation.csv` and matching ROOT matrices expose
cross-bin dependence; target covariance is a separate matrix.

The output is the virtual-photon `d2sigma/(dt dphi)` in the raw SIMC
`sigcm` convention (microbarn/MeV^2/radian). Eq. 5.30's printed four-fold
electron quantity would require a separately specified Hand flux and unit
conversion; no such display factor is inserted into the normalized response.

For normalization, the producer writes
`full_weight=Weight*normfac/Ngen`. The generator factorization recorded in
`docs/physics_audit_demodelled_smearing_xsec_20260911.md` is
`Weight=(generation/radiation/angular-Jacobian weight)*siglab` and
`siglab=sigcm*(hadron Jacobian)*(virtual-photon flux)*(Fermi factor)`.
Thus `full_weight/sigcm` retains the luminosity and generation phase-space
normalization, radiation, flux and Jacobians needed to map the fitted
virtual-photon cross section into weighted reconstructed yield per mC.
`normfac=L*genvol*Naccepted/Ntried` in the audited SIMC convention; dividing
by the histogram's `Ngen=Naccepted` leaves `L*genvol/Ntried`. On the data
side the combiner's `scale=PS/(Qrun_mC*LT*eff)` is multiplied by
`Qrun/Qtotal` and `pi0_weight`, giving corrected signal yield per total mC.
The target divisor is then applied once. No second luminosity, generation
volume, flux, Jacobian or bin-width factor enters the response.

For reconstructed row r, generated bin b and component a:

```text
A[r,b,a] = sum(events reconstructed in r, generated in b)
             (full_weight / sigcm) * [1, sqrt(2 eps (1+eps)) cos(phi),
                                      eps cos(2 phi)]_a / (2 pi)
expected_yield[r] = sum(b,a) A[r,b,a] sigma[b,a].
```

Cuts and row indices use reconstructed coordinates. Basis functions and column
indices use generated coordinates. `ti` stores positive -t; it is negated and
combined with the exact two-body forward limit to obtain signed t' = t-t_min.
Generated Q2, W, t and phi are mandatory. The no-simc-model extractor also
requires the matching original exclusive h10 via `--vertex_simc_file`:
`event_id` selects its entry, `sigcm` is checked event by event, and vertex
epsilon is calculated from `Q2i`, `Wi`, `hsxptari`, `hsyptari` and the HMS
central angle in the matching SIMC `.hist` (set as `hms_theta_deg` in the JSON
preset). The original h10 `epsilon` and a smeared `epsilon_i` are reconstructed
values and are not used in the vertex response.

The configured Q2/xB/t' boundaries are shared between
generated and reconstructed coordinates. Generated events outside these bounds
can migrate into the measured region. They enter six disjoint exterior regions:
t' below/above, then Q2 below/above, then xB below/above. This priority assigns
corners once. Each populated exterior region has three free nuisance coefficients
and contributes to the full covariance. This coarse exterior parameterization is
an implementation choice, not a prescription taken from either thesis. It must
be varied for a physics systematic study; it does not guarantee an arbitrary
outside-domain cross-section shape is represented.

The solver whitens rows by uncertainty and scales columns for numerical units,
then uses full SVD. It rejects unsupported or rank-deficient published bins;
the default has no smoothing, singular-value truncation, prior or positivity
constraint. Negative angular regions are counted and reported, never clipped. Defaults are
`--svd-rank-tolerance 1e-10`, `--mc-max-iterations 30`,
`--mc-fit-tolerance 1e-6` (see `--help`).

`--positive-xsec` optionally requires the complete angular cross section to be
nonnegative over all phi in every active truth block, including exterior
nuisances. LT and TT remain signed and zero cross section is allowed. The
constraint uses each block's largest observed generated epsilon, which covers
all smaller epsilon values in the same bin-constant model. The default remains
off; `--no-positive-xsec` explicitly selects the unconstrained fit. When
constraints bind, fitted central values are returned but statistical/MC
covariance and error bars are unavailable (NaN in numerical output, central
values only in plots). Separate unconstrained curvature is a diagnostic, not
a confidence interval. See [POSITIVITY.md](POSITIVITY.md) for scope, uncertainty
limits, and a command that writes to `output/KinC_x60_4b/xsec_positive`.

Data variance is the sum of squared event weights. MC variance uses event-level
outer products of the three basis terms, preserving their correlations. GLS
iterates this parameter-dependent variance to convergence. The unconstrained
covariance is conditional on final weights and is not multiplied by chi2/ndf. This uses
Poissonized MC counts; fixed-generated-count multinomial correlations and
detector/radiative systematics are not included. Same-fit yield uncertainty
accounts, to first order, for correlation between the fitted coefficients and
their response MC. The plotted data/prediction ratio is a residual diagnostic;
its error bar shows only data error divided by the fitted prediction.

Rows with zero observed data variance are explicitly excluded and recorded;
no pseudocount is invented. Sparse-bin and bin-choice effects require
closure with varied samples/binnings. Epsilon at generated Q2/xB currently uses
the nominal beam energy: event-dependent radiative beam energy is not available
in the reduced input contract. U's within-bin epsilon dependence is likewise
not separately resolved. These are stated approximations, not corrected physics.

## Thesis basis and provenance

The requested references are [Defurne, HAL tel-01281332v1](https://theses.hal.science/tel-01281332v1)
and [Ali, OSTI 1784736](https://www.osti.gov/biblio/1784736).

Ali's printed p.142, Eqs. (5.2)-(5.5) and Fig. 5.3, explicitly distinguish
generated/vertex bins from reconstructed bins and integrate component-dependent
angular factors over simulated events linking them. In particular the LT term
uses the vertex azimuth. The implementation above follows that forward-response
construction, extending it across Q2 and xB as well as t'. The basis-weighted
matrix is not a unit-sum stochastic matrix even though the thesis calls it a
migration probability matrix. Printed p.140 discusses an additional t' bin for
migration checks. Ali uses positive t_min-t; this code retains signed t-t_min.
These pages and Fig. 5.3 were checked again in the complete local copy of the
[exact OSTI-linked PDF](https://misportal.jlab.org/sti/publications/16620/attachments/7016/SFALI_THESIS_OCT2020.pdf#page=180);
the page was rendered to verify axis orientation and the signed LT color scale.

For Defurne, the complete local copy of the
[JLab thesis mirror](https://hallaweb.jlab.org/experiment/DVCS/documents/results/m_defurne.pdf)
was freshly inspected at printed pp.51, 76-77 and 97 (PDF pp.57, 82-83 and 103).
Fig. 3.9 illustrates DIS migration using vertex Q2/xB, with vertex incoming
electron energy as color and the reconstructed HRS acceptance contour. It is
not a normalized response-matrix plot. Eqs. (5.12)-(5.22) describe integrated
generated-to-reconstructed response; the pi0 extraction adopts it on p.97.
Fresh web retrieval of HAL/JLab failed; byte identity of the existing mirror
copy with the requested HAL v1 remains unverified.

The free exterior-region model, SVD rank policy, finite-MC GLS iteration and
empty-row policy are our implementation choices, not asserted prescriptions of
these theses. Their limitations and closure checks are documented separately.

## Migration plots

With diagnostics enabled (default), the extractor appends four pages to
`all_generated_plots_no_simc_model.pdf`, writes individual PDF/PNG files under
`migration/`, and saves the numerical histograms in the ROOT `migration/`
directory. `--no-diagnostics` skips these new products; with only `--no-pdf
--no-png`, their numerical ROOT histograms are still written.

- `response_components`: Ali Fig. 5.3-style U, LT and TT integrated responses.
  Horizontal coordinate is generated truth block; vertical coordinate is
  reconstructed row `((it*n_q2+iq)*n_xb+ix)*n_phi+ip`. Phi varies fastest.
  In a single Q2/xB grid, the `b` and `s` labels also show the positive -t'
  bin intervals in GeV^2; `g` labels identify exterior feed-in regions.
  Larger grids retain compact indices, whose bounds are in the output tables.
  The heatmap display converts the raw response to weighted yield per
  nb/GeV^2 coefficient; ROOT and CSV keep native response units.
  LT/TT use signed linear color scales; negative entries are interference
  coefficients, not negative migration probabilities. No matrix is normalized
  or multiplied by fitted coefficients for this display.
- `selected_migration_support`: aggregate unweighted selected-event counts
  over reconstructed phi first and show log10(1+selected MC events), so sparse
  support remains visible. Zero-event cells are blank.
- `selected_migration_fractions`: the left panel divides each truth column by
  its total; the right divides each reconstructed row by its total, including
  all guards. Outlined cells mark matching truth/reconstructed bin indices;
  other populated cells show migration or exterior feed-in.
  These answer where selected truth events reconstruct, and where a selected
  reconstructed sample originated. They depend on the sampled MC population.
  Detector-lost/rejected events are absent, so neither is an efficiency.
  Unsupported denominators are stored as zero and described in the caption.
- `vertex_reco_q2_xb_coverage`: Defurne-inspired generated/reconstructed Q2-xB
  views of exactly the MC entries entering the response. Vertex incoming energy
  is unavailable, so color is selected-event density per xB-Q2 bin area, not
  energy. Dashed lines are our reconstructed analysis boundaries, not an HRS
  contour. Display-only density scaling handles different extended bin widths;
  the ROOT histograms retain raw counts, including outside-range truth.

The extraction uses signed `t'=t-tmin`, opposite to the positive `tmin-t`
abscissa used in Ali. No event or fit convention is changed for these plots.
Empty guard columns may be hidden in figures; ROOT keeps all six guards in
their original index order. All reconstructed rows are included, even those
excluded by fit recovery. The figures diagnose the available response and do
not establish acceptance closure or constrain exterior coefficients.

New ROOT objects: `response_U`, `response_LT`, `response_TT`, `selected_counts`,
`p_reco_given_truth_selected`, `p_truth_given_reco_selected`,
`vertex_q2_xb_selected`, `reco_q2_xb_selected`, plus `definitions` and
`coverage_semantics`. All live under `migration/`. Existing CSV schemas remain.

## Output contract

### Global fit recovery

The full global fit is always attempted first. If a fit fails, the extractor
excludes Q2/xB groups with data outside MC support or unsupported generated
bins and retries a **single global fit of the remaining groups**. Exclusion
removes all reconstructed t'/phi rows in that Q2/xB group. Binning and target
correction are computed only once; retries do not rescale the data again.

Every truth block contributing to retained rows stays in the response, including
truth from excluded Q2/xB groups and exterior regions. These have free nuisance
coefficients and appear in the full covariance with `is_nuisance=1`. Recovery
can still fail if those contributions are not identifiable. They are never
fixed to zero to force a fit.

For remaining rank or convergence failures, try progressively smaller global
subsets, stopping at the first successful largest subset. Ties use increasing
excluded group index (`iq*n_xb+ix`), not fit chi2. There are at most `2^N-1`
nonempty attempts for N Q2/xB groups (15 at the default N=4); very fine Q2/xB
grids can make recovery expensive. This selection is data dependent; reported
covariance is conditional on the retained subset, and closure/selection-bias
checks are still required before using recovered results for physics.

No new command-line flag is needed. New `fit_status.csv` (one row per Q2/xB
group) and `fit_attempts.csv` (retained groups and exact error for every attempt)
also have ROOT trees named `fit_status` and `fit_attempts`. Excluded bins have
`fit_xsec_ok=0`, NaN fitted cross sections/predictions, and `fit_failure_reason`.
Their slice plots say "Excluded from global fit"; PARTONS skips these bins.
`fit_scope=global_migration_subset` identifies recovery. All successful slice
chi2/ndf values still refer to the same final global fit and must not be summed.
Reconstructed rows excluded with their group have `exclusion=excluded_q2_xb`;
raw response cells and data remain exported for diagnosis. Predictions for
excluded rows are NaN because their full truth model may no longer be fitted.

Partial success exits zero with a warning and saves the retained results.
If no subset succeeds, diagnostics are saved and the executable exits nonzero;
there is no fitted-parameter ROOT matrix in that case. Input/normalization/file
errors remain fatal. Run the standalone recovery checks from the repository:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash tests/run_xsec_bin_recovery_tests.sh /tmp/nps_xsec_recovery_check'
```

### Joint fit across kinematic settings

The joint workflow prepares data yields and SIMC response matrices directly
from each setting's raw ROOT inputs, then fits shared sigmaLT/sigmaTT and
independent sigmaU per setting in each truth block. No individual extraction
fit is required. From the repository root:

```bash
python3 src/xsec_extract/run_joint_xsec_fit.py \
  --prepare-setting xsec_config_x36_5_407.json output/simc/simc_x36_5_407/worksim/ \
  --prepare-setting xsec_config_x36_4.json output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/ \
  --binning-config xsec_config_x36_4.json \
  --mmiss-select ellipse --mmiss-lower 0.6 --mmiss-upper 1.1 \
  --positive-xsec --fit-variance finite-mc --partons \
  --out-dir output/joint_x36_5_407_x36_4_LH2
```

Each `--prepare-setting CONFIG VERTEX_SIMC` reads the data and smeared SIMC
paths from CONFIG. VERTEX_SIMC can be the original exclusive h10 ROOT file or
its worksim directory. The driver loads the Hall C environment automatically.
It uses production event selection, vertex matching, target correction and
response accumulation, then returns before individual fitting, fit recovery,
toys, PARTONS evaluation or plots. The individual extraction workflow remains
available through `run_xsec_pipeline.sh`.

`--binning-config` explicitly copies all bins and the diamond selection from
one preset into new effective config files. In this example it chooses the
x36_4 first tprime edge (-0.75) for both settings. Original presets remain
unchanged. Without this option, all setting bin definitions must already agree.
Missing-mass overrides apply only when preparing raw settings. Matching
ellipse/MCD geometry must exist beside each combined data ROOT file.

Preparation is saved separately in `OUT_DIR_inputs` (override with
`--inputs-dir`). Both the preparation and fit directories must be new.
Effective configs live in `configs/`; each kinematic subdirectory contains
`joint_input_metadata.txt`, `joint_input_slices.csv`, `migration_reco_rows.csv`,
`migration_truth_blocks.csv`, and `migration_response_cells.csv`. No individual
fit files are needed or produced. Completed inputs remain available even if
the subsequent joint fit fails its rank or convergence checks.

To refit those inputs without repeating preparation, use:

```bash
python3 src/xsec_extract/run_joint_xsec_fit.py \
  --setting output/joint_x36_5_407_x36_4_LH2_inputs/configs/KinC_x36_5_407.json output/joint_x36_5_407_x36_4_LH2_inputs/KinC_x36_5_407 \
  --setting output/joint_x36_5_407_x36_4_LH2_inputs/configs/KinC_x36_4.json output/joint_x36_5_407_x36_4_LH2_inputs/KinC_x36_4 \
  --positive-xsec --fit-variance finite-mc \
  --out-dir output/joint_x36_5_407_x36_4_LH2_refit
```

`--setting` also accepts previous no-SIMC-model extraction exports, provided
bins, selection, normalization and input provenance match. Earlier fitted
coefficients, fit objective, positivity and variance choices are not inputs
to the new joint fit. Full U/LT/TT event covariance is retained; a four-component
response export can supply the equivalent T/LT/TT submatrix when necessary.
Legacy epsilon maxima can come from positivity diagnostics; prepared inputs
store the maxima directly, without solving positivity first.

The joint objective currently remains **Gaussian chi-square**, with optional
iterated finite-MC variance; `--fit-objective scaled-poisson` is not supported.
`--fit-variance data` uses only observed data variance. Positivity applies to
each setting's U and shared LT/TT using that setting's event epsilon maximum.
After the joint fit, the driver separates `sigmaU = sigmaT + epsilon*sigmaL`
in each truth block using one fixed nominal epsilon per setting. For two
settings this is the exact intercept/slope solution; for three or more it is
generalized least squares using the joint U covariance within that block.
The numerical joint fit does not propagate target-factor uncertainty. Optional
PARTONS comparisons are evaluated afterward by the plotting stage.

Nominal epsilon is calculated from preset `ebeam`, `hms_p_gev`, and
`hms_theta_deg`, with the ultra-relativistic electron convention `E' = p`:
`Q2 = 4 E p sin(theta/2)^2`, `nu = E-p`, and
`epsilon = 1/[1 + 2(1 + nu^2/Q2) tan(theta/2)^2]`.
The x36_5_407 and x36_4 presets give **0.7115005066172831** and
**0.5163181761587335**, respectively. This nominal epsilon is used only in
the post-fit separation; event epsilon still defines the response and
positivity constraints. No epsilon uncertainty is propagated.

To override the nominal values, append `--nominal-epsilon VALUE` once per
setting, in command-line setting order. For example, with x36_5_407 first:
`--nominal-epsilon 0.7115005066172831 --nominal-epsilon 0.5163181761587335`.
Both preparation and saved-input modes accept this option. Presets without
`hms_p_gev` require these overrides or a central momentum added to the preset.
When reusing prepared inputs, repeat any desired overrides; otherwise the
saved preset's central kinematics are used.

Results include `joint_parameters.csv`, `joint_covariance.csv`, `joint_rows.csv`,
`joint_setting_chi2.csv`, `joint_positivity.csv`, `joint_summary.txt`,
`joint_xsec_output.root`, and `joint_manifest.json`. Parameter `setting_index`
is the zero-based setting order for U, and -1 for shared LT/TT. Covariance
indices follow `parameter_index`. A positivity boundary makes ordinary
Gaussian covariance/errors unavailable (NaN); curvature remains diagnostic.

Additional outputs are `joint_separated_parameters.csv` (T/L/LT/TT per block),
`joint_separated_covariance.csv` (indices refer to that parameter table), and
`joint_lt_separation.json` (nominal epsilon provenance, per-block status,
chi-square and degrees of freedom). The JSON is also embedded under
`lt_separation` in `joint_manifest.json`. The covariance propagates correlations
between U values, across truth blocks, and with shared LT/TT. The original
joint parameter/covariance tables and ROOT objects retain U/LT/TT; the added
separation products are CSV/JSON only.

Separation does not impose positivity on T or L. Bins supported by fewer than
two settings, or with indistinguishable nominal epsilons at `--rank-tolerance`,
receive NaN T/L and an explicit status, while shared LT/TT remain available.
At a positivity boundary, two distinct epsilons still give T/L central values
but errors remain NaN. With more settings and unavailable joint covariance,
GLS separation is unavailable; diagnostic curvature is never used as covariance.

Reproduce the synthetic raw-input/preparation/refit check with
`python3 tests/test_joint_xsec_preparation.py`. The numerical regression is
`csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 tests/test_joint_xsec_fit.py'`.
These use temporary artifacts and do not run production data.
Post-fit numerical checks (no ROOT required) are
`python3 tests/test_joint_lt_separation.py`.

#### Joint plot report

Joint runs now generate PDF and PNG plots by default. Add `--partons` for the
same native GK06/GPDGK19 comparison plots requested from the individual
pipeline; `--partons-warmups` and `--partons-calls` retain its defaults of
10000 and 100000. `--no-plots` runs only input preparation, fitting and L/T
separation. All ROOT/PARTONS compilation occurs in temporary build directories.
`run_xsec_pipeline.sh` is a read-only reference and is not changed or invoked
by the joint renderer.

The report includes all native plot families for each setting, when their
corresponding inputs/options apply:

| Plot family | Contents |
| --- | --- |
| Global kinematics | 1D distributions, Q2/xB, tprime and phi binning checks, data/SIMC occupancy |
| Mass selection | Missing-mass comparisons per slice, 2D ellipse/MCD geometry, selected/unselected 1D spectra and zooms |
| Epsilon | Generated reference epsilon versus tprime and slice index |
| Phi slices | Forward yields, residual ratios, angular cross sections and components, helicity diagnostics when available |
| Cross-section terms | U/LT/TT versus minus tprime for each Q2/xB group |
| Migration | Component responses, selected-event support/fractions, vertex/reconstruction coverage |
| PARTONS (`--partons`) | Native GK06/GPDGK19 U/LT/TT comparisons at each setting's generated reference kinematics |

Additional joint plots show U versus nominal epsilon (Rosenbluth lines and
covariance bands), separated T/L and shared LT/TT versus minus tprime,
independent U comparisons, T/L and LT/TT Gaussian confidence ellipses, full
joint and separated correlation matrices, per-setting residuals/chi-square
contributions/MC variance, and positivity envelopes with nominal/event epsilon.
Guards are included in correlation and positivity diagnostics; published bins
receive the separation and coefficient plots. Plot units match the native
pipeline: a factor of 1e9 converts raw SIMC coefficients to nb/GeV^2.

Outputs reside under `plots/settings/<kinematic>/` and `plots/joint/`.
`all_joint_xsec_plots.pdf` collects every individual PDF page once.
`plots/plot_manifest.json` records page paths, count, options and completion
status. Per-setting `render.log`, `joint_plot_slices.csv`,
`joint_plot_diagnostics.root` and `joint_plot_input.txt` support auditing.
A completed report is not overwritten. An interrupted report can be retried
with `--plot-only` and the same PARTONS option.

Native plots reaccumulate the original data/SIMC/vertex ROOT inputs, import
the setting's marginal joint coefficients/covariance, and verify yields,
variances, responses and predictions against the saved fit. No individual
fit or individual constrained-refit toys run. Thus old saved inputs still
need accessible raw files and recorded `vertex_source` for the complete
diagnostic suite. Numerical fit outputs are saved before plotting; a plot
failure does not discard the joint fit.

Curve and coefficient errors use joint covariance. Residual-corrected angular
points are explicitly central-only: the individual extractor's point-error
formula omits observations from the other joint settings and is not reused.
At a positivity boundary, Gaussian errors/ellipses remain unavailable rather
than substituting an individual fit's toy errors. PARTONS curves are point
predictions at setting-specific reference kinematics, not joint bin averages.

To add plots to an existing joint fit without refitting, run from repository
root (include `--partons` for theory comparisons):

```bash
python3 src/xsec_extract/run_joint_xsec_fit.py \
  --plot-only output/joint_x36_5_407_x36_4_LH2 --partons
```

Reproduce the full synthetic comparison against the native individual
extractor, including ellipse diagnostics and a small PARTONS smoke test:

```bash
python3 tests/test_joint_xsec_plots.py --partons
```

The test checks native plot coverage, imported coefficient values/errors,
absence of individual fits during joint rendering, PNG/PDF completeness,
combined PDF page count and the reference shell's unchanged checksum.

### Result files

Historical output filenames remain. Slice coefficients now belong to generated
bins. `fit_xsec_chi2` and `fit_xsec_ndf` refer to the one global fit, repeated
across slice rows; do not sum them. `fit_scope` is `global_migration` for a full
fit or `global_migration_subset` after exclusions.
Positive mode appends `_positive` to these scope names; see
[POSITIVITY.md](POSITIVITY.md) for its additional diagnostic outputs and
boundary-uncertainty status.
Reconstructed means remain diagnostics; generated reference means are explicit.
`partons_electron_flux_xbq2` now stores full Gamma, twice pi its former meaning.

Additional CSV files make the fit reproducible; the standard response exports
also have ROOT counterparts:

- `migration_design.csv`: integrated response matrix, all reconstructed rows.
- `migration_reco_rows.csv`: measured yields, observed and final solve variances,
  inclusion map, prediction and conditional residuals.
- `migration_parameters.csv`: complete U/LT/TT column mapping, including nuisances.
- `migration_covariance.csv`: full statistical+MC and separate target covariance.
- `migration_response_cells.csv`: per-cell U/LT/TT basis sums and all nine MC moments.
- `migration_joint_response_cells.csv`: per-cell T/L/LT/TT basis sums and
  all 16 event covariance moments for the joint fit.
- `migration_truth_blocks.csv`: generated reference moments and exterior labels.
- `migration_singular_values.csv`: singular values of the scaled fit matrix.

ROOT `analysis_metadata` records conventions and limitations. Use the complete
available covariance for subsequent fits, not only marginal slice errors.
Do not substitute the separate unconstrained curvature for unavailable
boundary-constrained covariance. No background
subtraction or detector-acceptance retuning is performed by this update.

## Reproduce the checks

From the repository root, the following runs standalone convention/solver tests,
compiles the extractor, checks failure paths, fits the existing KinC_x60_4b data,
and independently reconstructs the response and split-sample closure:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash tests/run_xsec_tests.sh /tmp/nps_xsec_check'
```

The optional second argument selects the repository containing fixture outputs.
An explicit SKIP means event fixtures were unavailable, not that closure passed.
See [VALIDATION.md](VALIDATION.md) for checks actually run and their limits.

To run the wrapper using the audited setting, with outputs in a new directory:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash src/xsec_extract/run_xsec_pipeline.sh --kin KinC_x60_4b --target LH2 --sim-file output/KinC_x60_4b/root/simc_pi0_analysis_output_smeared.root --vertex_simc_file output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x60_4b.root --out-dir /tmp/KinC_x60_4b_migration'
```

Add `--partons` only when a native model comparison is needed; it adds numerical
integration work and requires the installed PARTONS environment.
