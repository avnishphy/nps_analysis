# Validation performed 2026-09-16

These checks test implementation and stated statistical approximations. They
do not establish that background, detector modeling, radiation or point GK
predictions describe the measured cross section.

## Inputs and commands actually run

Code was staged under `/tmp/xsec_update_20260916`. With the Hall C/NPS environment
loaded, the complete suite was run as:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash /tmp/xsec_update_20260916/tests/run_xsec_tests.sh /tmp/xsec_update_20260916/validation/final'
```

It completed with exit status zero. The event fixtures, opened read-only, were
under `root_analysis_env_main/output/`:

- `KinC_x60_4b/root/combined_branches_LH2.root`
- `KinC_x60_4b/root/simc_pi0_analysis_output_smeared.root`
- `simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x60_4b.root`

The defaults used 0.6 <= missing mass <= 1.1 GeV, 5.1 <= Q2 <= 6.4 GeV^2,
0.45 <= xB <= 0.65, -1 <= t' <= 0 GeV^2, and target factor 0.584 +/- 0.014.
Selected counts were 37,666 data events and 210,270 MC events. Every selected
MC event had a unique matched raw identifier and a consistent model cross section.

## Results

| Check | Observed result |
| --- | --- |
| Analytic GK conversion | Full-Hand/native-prefactor, units, signed LT/TT and multiple phi/epsilon checks passed |
| SVD solver | Identity/mixed-bin closure, full covariance, extreme unit scaling and rejection paths passed |
| MC covariance | NaN, material negative variance and overflow rejected; roundoff behavior checked |
| Input contract | Missing generated input and both obsolete normalization aliases rejected |
| Real-data fit | 225 fitted reconstructed rows, 72 parameters, rank 72 |
| Real-data fit diagnostics | Scaled condition 244.221; chi2/ndf = 216.718/153; 9 MC iterations |
| Independent NumPy solve | Parameter relative discrepancy 3.35e-13; covariance 4.01e-14 |
| Manufactured sloped U, signed LT/TT and exterior terms | Relative algebra closure below 5e-13 |
| Independent event reconstruction of response | Relative discrepancy 7.99e-17 |
| Independent event reconstruction of MC variance | Relative discrepancy 4.27e-15 |
| Target factor | Parameter/covariance scaling and separate correlated covariance passed |

The 20 published generated bins give 60 parameters; four populated exterior
regions give 12 nuisance parameters. Three MC-supported rows had zero observed
data variance and were explicitly excluded. Other omitted rows had neither data
nor response. Full cross-bin covariance is substantial and must be retained.

For a reproducible, disjoint MC split (seed 20260916), 105,367 events built the
response and 104,903 independently generated pseudodata under prescribed binwise
U/LT/TT coefficients. The recovered coefficient vector gave chi2 = 68.02 for
72 dimensions; reconstructed residual chi2 = 160.10 for 156 degrees of freedom.
The largest marginal pull was 2.61. This statistical check is stronger than
same-matrix algebra closure but still uses the same detector simulation and a
representable binwise-constant injected cross section.

An intentionally incorrect synthetic fit omitting exterior response produced
up to 61.49 nb/GeV^2 coefficient bias. This demonstrates why feed-in cannot be
silently dropped; it does not bound the systematic from coarse exterior shapes.

Five real-data generated bins have a negative fitted angular region. These are
flagged without clipping or positivity constraints. At the most negative angle,
the magnitudes are roughly 0.11, 1.66, 0.11, 0.70 and 1.54 marginal standard
errors for flattened bins 0, 4, 6, 9 and 12. They remain diagnostic issues;
the extraction is not declared physically validated.

## Native model and plots

An optional native PARTONS test at Q2=5.8 GeV^2, xB=0.58, t=-1 GeV^2,
E=10.538 GeV linked successfully and reconstructed its native electron observable
with zero measured round-trip error. Full Hand Gamma was 0.00015297672532529311
and epsilon was 0.757607321260502. Low-statistics convolution warnings occurred
at 1,000 warmup / 10,000 integration calls. This is an adapter/link smoke test,
not a converged GK numerical prediction.

A separate actual-data PDF plotting run completed and produced a readable
32-page combined PDF. ROOT emitted PDF object-state warnings while interleaving
individual and combined exports; these are a remaining plotting issue. Numerical
outputs and the standalone fit tests were unaffected.

## Remaining physics work

Thesis verification: Ali's migration equations (5.2)-(5.5), printed p.142, were
checked through indexed primary-PDF text. Fresh full-PDF retrieval failed;
Defurne's p.97 guidance is supported only by the repository's previous reading
notes. See the extraction guide for links and explicit provenance.

Vary exterior-region parameterization, configured bin edges and the
missing-mass selection; study low-occupancy rows, within-bin shapes and nominal
beam-energy epsilon in radiative events. Validate phi/LT sign conventions,
converge native GK integration, and evaluate detector/background uncertainties.
The update adds no smoothing, positivity constraint, external model prior or
data-driven absolute normalization to conceal these effects.

Re-run after installation using the repository-relative command in
[README.md](README.md#reproduce-the-checks). The test runner writes the exported
fit, detailed log and machine-readable `closure.json` to its chosen scratch path.

## 2026-09-24: KinC_x36_5_407 Eq. 5.30-5.31 and PARTONS check

The updated extractor was compiled under ROOT 6.30.04. The standalone xsec
suite and bin-recovery suite passed. The x60 event-closure fixture named by
the existing runner was unavailable, so real-event output checks used the
KinC_x36_5_407 inputs below. All output was isolated under `/tmp`.

The requested pipeline command was run with only `--out-dir` added:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash src/xsec_extract/run_xsec_pipeline.sh --kin KinC_x36_5_407 --target LH2 --root-dir output/KinC_x36_5_407/KinC_x36_5/root/ --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x36_5_407/smearing_output/KinC_x36_5_407/root/simc_pi0_analysis_output_smeared.root --vertex_simc_file output/simc/simc_x36_5_407/worksim/ --mmiss-lower 0.6 --mmiss-upper 1.1 --partons --out-dir /tmp/nps_xsec_eq531_x36_partons'
```

It exited zero: 236 fitted rows, 72 parameters, full rank, scaled condition
72.4067, chi2/ndf 162.63/164, and nine finite-MC variance iterations.
The native GK06/GPDGK19 comparison returned finite, nonzero coefficients in
all 20 published generated bins. Its 240 convolution warnings included NaN
terms replaced by zero, so numerical convergence of the theory points has
not been established. PARTONS did not alter extraction: parameter and
covariance CSV files were byte-identical to a separate finite-MC run without
PARTONS on the same input. ROOT also reported PDF object-state warnings; the
combined 67-page PDF was readable by `pdfinfo`.

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 tests/validate_migration.py /tmp/nps_xsec_eq531_x36_partons --json /tmp/nps_xsec_eq531_x36_partons/migration_validation.json'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 tests/validate_experimental_points.py /tmp/nps_xsec_eq531_x36_partons'
```

Both passed. Independent NumPy parameter/covariance/correlation relative
discrepancies were 3.49e-14, 3.04e-15, and 6.98e-15. Of 240 corresponding
experimental-point cells, 235 were available and five had an explicit
`no_diagonal_response` status. The separate finite-MC point-covariance check
reproduced nine selected matrix entries with 1.02e-13 relative discrepancy.
The full data-only reference-mode point covariance was independently checked
on the same input (235 x 235 available cells; 1.18e-13 relative discrepancy).

A final point-uncertainty refinement after the `--partons` run added the MC
fluctuation of the reference epsilon, with event cross moments in
`truth_phi_response_cells.csv`. The final source was compiled and the xsec
and bin-recovery suites rerun successfully. A separate final-source
finite-MC run, with the same input and `--kin KinC_x36_5_407` but without
plots or PARTONS, wrote to `/tmp/nps_xsec_eq531_x36_finite_latest`.
`tests/validate_experimental_points.py` passed there with 1.02e-13
relative discrepancy for the nine independently calculated covariance
entries. Relative to the completed `--partons` run, fit coefficients and
their covariance were byte-identical, as were all experimental-point
central values; the largest relative point-error change was 0.00265.
The native PARTONS run was not repeated because this change only affects
the separately exported point covariance and reference-epsilon moments.
The exact final source was subsequently compiled and run in both fit modes
under `/tmp/nps_xsec_eq531_x36_final` (`finite-mc`) and
`/tmp/nps_xsec_eq531_x36_data_final` (`data`), with `--kin
KinC_x36_5_407 --no-pdf --no-png --no-diagnostics --quiet` and the same
data/SIMC/vertex paths and missing-mass range above. Both independent
validators passed. The data-only fit had 236 rows, rank 72, condition
85.3242, chi2/ndf 203.602/164, and zero MC iterations; its entire
235 x 235 available point covariance and correlation agreed with NumPy
to 1.18e-13 and 1.09e-13 relative discrepancy. The final finite-MC fit
had the same diagnostics as the PARTONS run and the nine checked point
covariances agreed to 1.02e-13.

## 2026-09-24: Plot cleanup and flag matrix

With the required Hall C/NPS setup and ROOT 6.30.04, the extractor compiled:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ -O2 -std=c++17 src/xsec_extract/excl_xsec_pi0_analysis_no_simc_model.C `root-config --cflags --libs` -o /tmp/nps_xsec_plot_cleanup_extractor'
```

The existing xsec, bin-recovery, migration-plot and migration-coverage tests
all passed. The x60 event fixture in the xsec runner was unavailable, so the
same KinC_x36_5_407 input paths and missing-mass window specified above were
used for representative flag runs. The current user configuration has one Q2
bin and one xB bin: the finite-MC fit had 60 rows, 27 parameters, full rank,
chi2/ndf 33.4556/33 and 11 MC iterations.

| Flags / output directory | PDFs | PNGs | Result |
| --- | ---: | ---: | --- |
| Full graphics, `/tmp/nps_xsec_plot_cleanup_final_full` | 23 (22 pages plus combined) | 22 | Exact final source; combined PDF readable; no ROOT TPDF object-state warnings |
| `--no-pdf --no-png --no-diagnostics`, `/tmp/nps_xsec_plot_cleanup_nographics` | 0 | 0 | Fit parameters, covariance, points and point covariance byte-identical to full graphics |
| `--fit-variance data --positive-xsec --no-png --no-diagnostics`, `/tmp/nps_xsec_plot_cleanup_pdf_only` | 7 (6 pages plus combined) | 0 | Active boundary labeled as having unavailable intervals |
| `--no-pdf --no-diagnostics`, `/tmp/nps_xsec_plot_cleanup_png_only` | 0 | 6 | Only slice and coefficient PNGs |
| Wrapper `--partons --partons-warmups 1000 --partons-calls 10000 --no-pdf --no-diagnostics`, `/tmp/nps_xsec_plot_cleanup_partons_png` | 0 | 9 | Five GK predictions; slice/coefficient/model PNGs only |

The wrapper also forwarded all three disabling switches together and selected
`--fit-variance data` successfully in `/tmp/nps_xsec_plot_cleanup_wrapper`.
The reduced-call PARTONS check tests output routing and layout; its native
convolution warnings mean it does not establish model convergence.

## 2026-09-24: Positive-fit error bars and generated t-prime axis

The user's KinC_x36_5_407 PDF had an active continuous-angle positivity
boundary. Its five coefficient errors were NaN by design, so the prior plots
showed no vertical bars. With `--mmiss-lower 0.7 --mmiss-upper 1.05`, the
full-rank fit has 60 reconstructed rows, 27 parameters, condition 87.981,
chi2/nominal ndf 66.5826/33 and four finite-MC iterations. The central U,
LT and TT values are exactly unchanged from the original PDF's slice CSV.

The exact user flag combination, with only a separate output directory added,
completed:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; ./src/xsec_extract/run_xsec_pipeline.sh --kin KinC_x36_5_407 --target LH2 --root-dir output/KinC_x36_5_407/KinC_x36_5/root/ --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x36_5_407/smearing_output/KinC_x36_5_407/root/simc_pi0_analysis_output_smeared.root --vertex_simc_file output/simc/simc_x36_5_407/worksim/ --mmiss-lower 0.7 --mmiss-upper 1.05 --partons --positive-xsec --out-dir /tmp/nps_xsec_positive_plot_full'
```

Its 25-page combined PDF and 25 individual PNGs were produced. All five GK
comparison points were present. The 256 refits yielded 253 converged toys;
all 27 parameter and 60 available experimental-point toy marginal variances
were positive. Both covariance matrices were symmetric with positive
eigenvalues (parameter min/max 4.06e-18/1.97e-10; point
2.19e-19/1.45e-15 in raw squared units). Visual inspection confirmed
vertical toy bars, -t' bin-span whiskers, centered axis labels and nb-based
cross-section units on coefficient, angular and GK comparison pages.

`bash tests/run_xsec_tests.sh` passed after sourcing the ROOT environment;
its x60 event closure fixture was unavailable as before. A data-only positive
run with `--fit-variance data --no-pdf --no-png --no-diagnostics` completed
under `/tmp/nps_xsec_positive_plot_data_nographics`: 60 rows, rank 27,
zero MC iterations, 256/256 converged refits, 60 positive point toy
variances, and no PDF/PNG files. Its exported `variance_used` equals
`data_variance` and `mc_variance_used=0` in every fit row. The plot-only toy
covariance remains separate from the NaN constrained confidence covariance.

The subsequent bin-migration refinement adds a separate selected-MC support
page. The final exact-flag run used the same command above with output changed
to `/tmp/nps_xsec_plots_all_flags_final`. It produced a readable 26-page PDF
and 26 PNGs; five GK comparison points and 253/256 toy refits were present.
All 27 parameter and 60 point toy diagonal variances were positive, and the
five bins' U/LT/TT central coefficients were exactly equal to the original
production CSV. Visual inspection confirmed the support heading, -t' ranges,
guard labels, and diagonal outlines on the support and fraction pages.

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ -O2 -std=c++17 tests/test_xsec_migration_plots.C `root-config --cflags --libs` -o /tmp/nps_xsec_migration_plots_test'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; /tmp/nps_xsec_migration_plots_test /tmp/nps_xsec_migration_plot_fixture_final'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 tests/validate_xsec_migration_plots.py /tmp/nps_xsec_plots_all_flags_final --baseline /tmp/nps_xsec_positive_plot_full'
```

The synthetic fixture passed support/fraction normalization, guard and
zero-support handling, unchanged raw response/fit, image switches and palette
restoration. The independent final ROOT/CSV validator passed with 131572
selected MC entries and matching numerical records for design, parameters,
covariance, reconstructed rows, response cells and truth blocks. Native GK
convolution emitted its existing null-error/range warnings; this run confirms
plot wiring and five model points, not model-integration convergence.
