# Cross-section production modes

Run from the repository root with the Hall C ROOT environment loaded. The
wrapper compiles a standalone C++ executable using the selected JSON preset;
it does not invoke a ROOT macro through `root -q` and does not combine data.
Python numpy/uproot are required for artifact validation; model reports also
use matplotlib. PDF assembly uses `pdfunite`, as do the common ROOT plots.

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash'
```

Inside that initialized bash shell, these are production examples (not claims
that the current combined input passes validation):

```bash
common=(--kin KinC_x36_4 --xsec_config xsec_config_x36_4.json
  --data-file output/KinC_x36_4/root/combined_branches_LH2.root
  --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x36_4/smearing_output/KinC_x36_4/root/simc_pi0_analysis_output_smeared.root
  --vertex_simc_file output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x36_4.root)
bash src/xsec_extract/run_xsec_pipeline.sh --mode no_simc_model "${common[@]}"
bash src/xsec_extract/run_xsec_pipeline.sh --mode simc_model "${common[@]}"
```

No-model remains the default. `--xsec-method no-simc-model|simc-model`,
`--no-simc-model`, `--simc-model`, and `NPS_XSEC_METHOD` remain supported.
Both modes share the same input, selection, raw-event matching, target,
objective, positivity, finite-MC, and plotting controls. Model-specific
minimizer flags are forwarded to the model entry point. Model parameter names
and formulas are not interpreted by the wrapper or production plotter.

The JSON must be an `xsec_config*.json` preset in `src/xsec_extract/xsec_config/`.
It supplies physics/binning defaults. `--kin` selects canonical paths and
infers the preset name. Explicit CLI paths override environment defaults;
without explicit data/SIMC paths, the wrapper resolves them beneath `root-dir`.
The raw SIMC file or worksim directory must be supplied. No reconstructed
epsilon fallback exists. `--fit-objective scaled-poisson` enforces positivity.
The geometric selectors `--mmiss_select ellipse|mcd` use verified combined cut
metadata, exporting it with the existing helper when necessary.

## Output contract

```text
output/<kin>/xsec/                # no-model; historical names preserved
output/<kin>/xsec_simc_model/     # model
  excl_xsec_pi0_analysis_<mode>_output.root
  excl_xsec_pi0_analysis_<mode>_summary.csv
  excl_xsec_pi0_analysis_<mode>_slice_summary.csv
  all_generated_plots_<mode>.pdf
  global/, slices/, sigma_vs_tprime/, migration/, ...  # shared ROOT pages
  model/                        # new report pages, PDF/PNG
  model_*.csv, model_fit.root, model_vertex_epsilon.root
  logs/pipeline_<mode>_<UTC timestamp>_<pid>.log
  pipeline_config.json, pipeline_config.h
  pipeline_status.csv, pipeline_artifacts.json, model_plot_manifest.json
```

Explicit output overrides remain available. Use a separate directory and
separate file paths for each kinematic/mode. Directory ownership markers reject
cross-mode/cross-kin reuse; a directory lock rejects simultaneous writers.
Same-mode reruns replace their products, retain unique logs, and perform no
automatic resume. Old optional files may remain when flags disable their
generation; only artifacts listed in the current manifest belong to that run.
`pipeline_status.csv` records mode, stage, exit code, and start time. Consumers
must require `stage=complete`, `exit_code=0`, and a matching `started_ns` in
`pipeline_artifacts.json`; an old manifest is not evidence of a successful rerun.

`SUCCESS` requires extractor exit zero plus nonempty, fresh ROOT/CSV/covariance
and enabled PDF artifacts. ROOT directories must be readable and nonempty.
Model report PNG/PDF files listed in its manifest are also checked. Artifacts
are hashed in `pipeline_artifacts.json`. `--no-pdf` disables individual and
combined PDF requirements; `--no-png` disables PNG; `--no-diagnostics` retains
core model curves and run summary. `--prepare-forward-inputs` remains a distinct
no-model contract, checked for both event CSVs and the forward-cache manifest.

Failures identify kin, mode, stage, exit code and log. Missing input files or
the `analysis_runs` manifest fail in `upstream-input-preparation`, before
compilation/extraction. The extractor retains authoritative run exposure,
current run-status, matching, sigcm, and covariance checks. A failed run must
be repaired by the analysis producer and then recombined. No manifest is
invented, and no failed run is silently dropped.

The upstream `scripts/run_pipeline.sh` delegates to
`src/analysis/run_parallel_nps_analysis_main.sh`: run analysis is parallel;
per-kin combine/smearing/extraction is serial, with existing failure handling.
Its extraction subprocess inherits `NPS_XSEC_METHOD=simc_model`; no second batch
implementation or changed batch failure policy is introduced.

## Plot and CSV interpretation

All common ROOT plot families remain available in model mode. Global SIMC
shape overlays are diagnostic area matches, not absolute model predictions.
Migration fractions describe selected MC and must not be read as acceptance
efficiencies. The added model pages explicitly separate Born structure
functions, folded detector yields, and raw response coverage.

The combined model PDF opens with the run summary, retains the common ROOT
pages, then adds input statistics,
coverage, row and phi residual/pull/component panels, structure curves and
optional overlays, correlations, parameter tables and identifiability.
Standalone model report pages live in `model/`; raw covariance and original
canvases remain in `model_fit.root`. The four legacy standalone model preview
plots are retained but omitted from the combined PDF to avoid repeating the
new structure-function and detector pages.

The production model `sigparam2021_pi0` uses a provisional charged-inspired
baseline with a smooth non-pole neutral longitudinal term, folded at generated
event level. Global U/LT/TT normalizations and one pivoted U slope are fitted;
LT/TT shapes remain fixed. The response-weighted tau0 is fixed and exported in
model_context.csv. Use --model-fix-u-slope to reproduce the three-normalization
model, or --model-before-dir for a diagnostic-only before/after comparison.
There is no physics-model selector: choose `--mode simc_model` or
`--mode no_simc_model`. See [SIGPARAM2021_PI0.md](SIGPARAM2021_PI0.md) for the
model/fit contract, current shape diagnostics and reproducible checks, and
[the validation report](../../validation/xsec_pipeline_20261004/REPORT.md)
for this implementation's measured results and remaining production blocker.

The optional comparison defaults to sibling
`xsec/excl_xsec_pi0_analysis_no_simc_model_slice_summary.csv` and can be selected
with `--reference-slice-csv` or `NPS_XSEC_REFERENCE_SLICE_CSV`. Absence is normal
and skips comparison. These points are read only after fitting; they never
enter the likelihood. A user-supplied comparison must have corresponding
kinematics/selection; the plot is diagnostic, not a compatibility certification.

Existing CSV columns are preserved; additions are appended:

| File | Stable identifiers and meaning |
|---|---|
| `model_parameters.csv` | index/name/role; value, confidence error, raw Hesse width, lower/upper bounds, at-bound and relative-error flags |
| `model_fit_status.csv` | convergence/status/EDM/objective, rows/parameters/nominal DOF, MC convergence, confidence validity; model identifier, statistic, physical/tprime-nuisance/fixed-feed-in counts |
| `model_covariance.csv` | i,j,role; confidence covariance, correlation, raw Hessian covariance |
| `model_structure_functions.csv` | truth_block,tprime,tau; U/LT/TT and confidence errors at physical reporting points |
| `model_structure_curves.csv` | tprime,tau,component,value,error; values/Jacobian propagation evaluated in C++ |
| `model_structure_covariance.csv` | pairs of truth-block/component ids; confidence and raw propagated covariance |
| `model_reconstructed_yields.csv` | row/included; data,sumw2,prediction,residual,objective pull/variance; it/iq/ix/ip, phi bounds in radians, signed folded U/LT/TT plus fixed Q2/xB feed-in, prediction error |
| `model_identifiability.csv` | parameter name/role/value/error/relative error, bound and weak flags |
| `model_strong_correlations.csv` | named pairs with absolute raw Hessian correlation above 0.9 |
| `migration_truth_blocks.csv` | original block id, active index, region, selected events, response strength and generated means; inactive blocks remain present |
| `model_plot_manifest.json` | schema version, kin/mode, generated page artifacts, confidence validity and reference-overlay path |

Internal structure-function units remain microbarn/MeV^2; tprime/tau are
GeV^2; yields are per mC after the configured target divisor. Gaussian pulls
use the actual objective variance; scaled-Poisson pulls are signed square roots
of row deviances. Excluded rows have no objective pull. Prediction bands are
conditional parameter-covariance bands, not independent data errors.
Signed component sums reproduce the full detector prediction; components are
drawn separately, never positive-stacked. Zero-response blocks have no measured
generated mean and are listed explicitly instead of plotted at fabricated tprime.

Boundary/invalid-covariance states retain central curves, signed predictions,
status and raw Hessian diagnostics. Confidence errors remain NaN, error bars
and bands are omitted, and every model page carries a warning. Raw correlation
plots are labeled diagnostic. Nominal objective/DOF is not a calibrated p-value.

## Reproduce validation

The current U-slope model, closure results, three-versus-four-parameter comparison,
unit checks, and complete PDF inspection are documented in
[the current model report](../../validation/xsec_uslope_20261004/REPORT.md).

In the initialized shell, choose a new output directory:

```bash
bash src/xsec_extract/tests/run_sigparam_validation.sh \
  src/xsec_extract/xsec_config/xsec_config_x36_4.json /tmp/my_sigparam_validation
```

This compiles both modes and runs charged Fortran reproduction, event/tprime-nuisance
Jacobian and MC checks, Gaussian/scaled-Poisson four-parameter closure, 17 Python
tests including display units and ratio guards, matching-rejection and empty-row fixtures.
The accompanying phase report records actual two-mode wrapper runs, reference
CSV comparisons, artifact-failure probes and rendered PDF inspection.
