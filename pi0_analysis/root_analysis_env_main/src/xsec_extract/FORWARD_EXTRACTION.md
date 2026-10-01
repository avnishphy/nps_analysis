# Independent-grid forward extraction

The new mode keeps the vertex U/LT/TT response and separates reconstructed
observations, truth coefficients and publication bins. It is an explicit
alternative to the existing shared-grid plotting pipeline. Existing presets
and default extraction remain unchanged.

The event cache is selected with the same reconstructed mass criterion and
Q2/xB diamond as the original extractor, but **before** rectangular t', Q2 and
xB cuts. It preserves original event IDs, vertex kinematics and normalization.
Changing a common mass/diamond cut requires a new cache. Widening a subsequent
rectangular range cannot restore events removed by the common diamond or
upstream reconstruction.

## Prepare events

From repository root, load ROOT and create a new output directory. This command
matches the supplied x36_4 extraction's ellipse selection; edit input paths for
other production samples. The existing xsec directory is not overwritten.

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash src/xsec_extract/run_xsec_pipeline.sh --kin KinC_x36_4 --data-file output/KinC_x36_4/root/combined_branches_LH2.root --sim-file /lustre24/expphy/volatile/hallc/nps/singhav/nps_smearing/smear_x36_4/smearing_output/KinC_x36_4/root/simc_pi0_analysis_output_smeared.root --vertex_simc_file output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x36_4.root --mmiss_select ellipse --mmiss-cut-file output/KinC_x36_4/root/combined_branches_LH2_combined_2d_mass_cut_debug.txt --mmiss-lower 0.6 --mmiss-upper 1.1 --prepare-forward-inputs --out-dir output/KinC_x36_4/forward_cache'
```

New files are `data_events.csv`, `mc_events.csv`, the effective generated
`forward_export_config.h`, `forward_export_provenance.json`, and
`forward_cache_manifest.json`. The completion manifest is committed last;
existing/partial caches are rejected. ROOT UUIDs, source/configuration hashes,
selection, event counts and normalization conventions are saved.

Data weights are yield/mC after exactly one division by the target factor.
MC retains `full_weight/sigcm`; do not normalize its columns, divide by another
bin width, or multiply an additional luminosity/flux. Every fitted coefficient
is in nb/GeV^2: the Python harmonic kernel includes `1e-9/(2*pi)`. All angles
are radians and t'=t-t_min is negative in the physical region. Phi edges may
start at any angle but must span one complete period.

## Fit independently configured domains

`forward_x36_4_reference.json` reproduces the rectangular grid in the supplied
fit. `forward_x36_4_buffer_candidate.json` is a **validation candidate**:
reconstructed t' extends to -1 GeV^2 and retains six bins; the truth buffer
[-1,-0.55] has one independent harmonic triple; only [-0.55,0] is reported.
It imposes a binwise-constant shape over that broader buffer, which needs
independent shape stress tests. No candidate is automatically certified.

```bash
python3 -s src/xsec_extract/run_forward_xsec.py \
  --cache output/KinC_x36_4/forward_cache \
  --config src/xsec_extract/xsec_config/forward_x36_4_buffer_candidate.json \
  --output output/KinC_x36_4/forward_buffer \
  --bootstrap 200 --seed 20261001
```

The `-s` disables user-site Python packages. On the inspected host this selects
compatible system NumPy 1.23.5 and SciPy 1.9.3; the user-site NumPy 1.26.4
conflicts with that SciPy's declared range. It does not alter the installation.
Use a supported NumPy/SciPy environment on other hosts. Runtime versions and
exact Python source snapshots are included with each result.

Configuration separates `reco`, `truth` and `publication_truth_blocks`.
Interior truth IDs iterate t', Q2, xB; xB counts may vary between Q2 bins.
Publication selection does not remove the other supported truth coefficients.
All supported exterior regions remain explicit nuisance triples, with corner
priority t', Q2, xB. No automatic rank truncation, region deletion, edge
clamping, exterior prior or regularization is introduced. A required unsupported
publication bin or a rank-deficient fit fails explicitly.

The fit uses the full harmonic prediction and exact continuous-phi positivity;
LT and TT remain signed. Its scaled-Poisson objective includes empty rows,
borrowing an explicitly recorded pooled weight scale. This is a conditional
two-moment approximation for positive weighted events, not an integer-Poisson
likelihood for efficiency-corrected counts. Signed weights are rejected.

The event bootstrap regenerates the data sums and MC response, then refits
every coefficient, including exterior nuisance terms. Copies of the same
original event share their multiplier. Fixed full-cache event catalogs make
the same seed a genuinely paired comparison across different reconstructed
partitions. Accepted-event MC normalization is kept fixed (Poissonized MC).
No failed fit is silently dropped: failures, successful indices, and conditional
success-only covariance/quantiles are all recorded. Failed replicas set a
requires-reevaluation status; their surviving covariance is not a usable
publication error matrix by itself.

## Read outputs and compare binnings

- `fit.json`: coefficients, predictions, full objective, optimizer information,
  positivity boundaries and empty-row scales.
- `coefficients.csv`: nb/GeV^2 coefficients and explicit publication flags.
- `response_diagnostics.json`: geometric SVD modes, support/effective MC counts,
  exterior contributions and information after nuisance profiling.
- `variance_diagnostics.json`: data/MC variance decomposition and, when all rows
  have positive variance, local Gaussian information. This is diagnostic only.
- `conditional_bootstrap.json`: all successful samples, failed indices/reasons,
  covariance, quantiles and the precise statistical limitations.
- `target_covariance.json`: separate common target-normalization covariance.
- `response_problem.npz`, `bin_metadata.json`, `provenance.json`, source snapshots:
  numerical problem and reproduction metadata.
- `status.json`: completion, fit/refit status, and remaining publication checks.

Optional `--profiles requests.json` takes a list of
`{"name":"U_first","functional":[1,0,0,...],"values":[0,10,20,...]}`.
Supply complete numeric arrays; ellipses above are notation only. Each point
refits all other harmonics and nuisance coefficients. These are fixed-response
profile curves; standard deviance thresholds are not coverage-calibrated here.

Use `compare_forward_xsec.py` to compare **common bin integrals**, retaining
the paired covariance. Its observable JSON is a list such as

```json
[{"name":"U_mid","component":"U","q2":[3.3,4.7],
  "xb":[0.29,0.44],"tprime":[-0.55,-0.40]}]
```

```bash
python3 -s src/xsec_extract/compare_forward_xsec.py \
  output/KinC_x36_4/forward_reference output/KinC_x36_4/forward_buffer \
  --observables common_integrals.json --output paired_comparison.json
```

Both fits need the same full event cache, statistical source and bootstrap seed.
The comparison rejects unaligned partial truth bins and unavailable publication
regions; it never silently adds a bin-centering model. t'-integrated harmonic
units are nb at the specified Q2/xB bin. Paired covariance accounts for the
same measured events entering both fits. Failed-refit selection remains a
limitation even when common successful replicas can be compared.

## Validation and remaining physics

The original direct-SVD solve is retained. Both C++ positivity solvers now use
the SVD inverse-information factor directly, avoiding Cholesky of a numerically
degraded covariance. The legacy scaled-Poisson path rejects Minuit covariance
status 2 (forced positive definite), requiring status 3 for that calculation.

Tests can be run from repository root:

```bash
PYTHONNOUSERSITE=1 python3 -m unittest discover -s tests -p 'test_forward_xsec_*.py'
python3 -s -m unittest discover -s src/xsec_extract/tests -p 'test_forward_xsec_problem.py'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 tests/test_forward_export.py'
```

The synthetic ROOT test checks mass cuts, deferred rectangular cuts,
normalization/units, matched truth and overwrite refusal. Python tests include
independent grids, shifted phi origin, exact angular minima, signed interference,
empty cells, repeated event IDs, paired fluctuations, all-nuisance profiles,
rank failures, yield-unit invariance and end-to-end CLI output.

### Measured x36_4 validation, 2026-10-01

The 34 Python tests, compiled ROOT exporter test, direct-SVD regressions and
existing five-case joint-fit synthetic test passed. The real selected cache
contains 162414 matched exclusive MC events. With the reference binning, the
cached response, measured yields and sumw2 reproduce the archived CSVs to
relative differences of 3.44e-15, 1.17e-15 and 1.42e-15, respectively. This checks
normalization and response construction; it is not independent physics closure.

Real fits and 200 paired data/MC refits are saved under
`/scratch/singhav/pi0_forward_validation_20261001`:

| Result directory | Scaled condition number | Successful refits |
| --- | ---: | ---: |
| `reference_200` | 38.1423 | 86/200 |
| `buffer_200` | 20.3861 | 200/200 |
| `buffer_phi24_200` | 18.2674 | 195/200 |

These are the new objective's numerical condition numbers in `fit.json`:
columns are normalized after weighting rows by
`1/sqrt(s_r * max(y_r,s_r))`, with `s_r` the recorded effective weight scale.
They are not the geometric SVD or the legacy Gaussian condition numbers.
The reference loses rank in many refits because of its sparsely sampled
exterior block. Improved buffer conditioning comes with a broader binwise
constant exterior description; validate that shape assumption independently.

The 24-phi-bin check uses the same buffer truth grid and a phi origin of pi/24.
Its five failed replicas (20, 37, 128, 157, 158) all lose MC support in the row
with t' in [-1,-0.75] and phi in [1.1780972451,1.4398966329] radians. That row
has three independent MC events and one data event. All three MC multipliers
are zero while the data multiplier is positive. The bootstrap probability is
`exp(-3)*(1-exp(-1)) = 0.03147`, consistent with five observed failures.
This is an observed MC-support limit, not a solver failure. More MC or a
justified coarser reconstructed partition is needed before using that setup.

For the common t' integral [-0.55,0], the 12-to-24-phi shifts in U/LT/TT are
[0.02584,0.29836,0.26340] nb. Paired standard deviations of the differences are
[1.06552,1.06902,2.90481] nb, conditional on the 195 shared successful refits.
See `paired_phi_binning.json`; `paired_common_integrals.json` compares the
reference and buffer on 86 shared successful refits. These conditional
comparisons do not establish coverage or justify omitting failed replicas.

The existing cache can be reused in the fit command above by setting
`--cache /scratch/singhav/pi0_forward_validation_20261001/cache`; always give
a fresh output directory. The shifted-phi configuration is saved alongside
these results as `buffer_phi24_shifted_config.json`. The repository's
`xsec_config/forward_x36_4_common_integrals.json` supplies the common U/LT/TT
observables to the comparison command.

The workflow always records `publication_ready=false`: conditional bootstrap
spreads and nominal fixed-response profile scans are not calibrated confidence
intervals. Originally empty empirical data rows remain empty in this bootstrap;
unsampled MC support cannot be regenerated. Fixed-N generator correlations,
normalization uncertainty, detector/radiative systematics and the joint
upstream signal-weight uncertainty require separate treatment.

Actual background parameter covariance exists at
`output/KinC_x36_4/plots/run_<run>/combbg_run<run>_order4_results.root`, under
`run_<run>/background_fit_covariance`. It alone does not include the dependence
of fitted purity weights on the same event data or all accidental-subtraction
uncertainty. Refit the upstream yield estimation in independent ensembles, or
provide its complete joint covariance; do not insert the 3x3 matrix as an
independent event-noise term.

Publication additionally needs independent within-bin and exterior-shape closure,
detector/radiation validation, and calibrated intervals with failures included.
The physical harmonic basis and migration buffer are supported by
[Dlamini et al.](https://arxiv.org/pdf/2011.11125); weighted-count approximations
by [Bohm and Zech](https://arxiv.org/pdf/1309.1287); the need to test coverage
separately from apparent smoothness by
[Brenner et al.](https://arxiv.org/pdf/1910.14654).
