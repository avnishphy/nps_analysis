# ALG-002B simultaneous timing/mass shadow implementation

Date: 2026-10-01
Setting: `KinC_x36_4`
Run scope: configured `production` runs with target exactly `LH2`
Status: implemented shadow fitter; validation gates remain open

## Approval and boundary

The user approved the model in
`docs/proposals/ALG-002B_simultaneous_timing_mass_model.md` and added the
requirement that only LH2 runs be used. The implementation therefore locks the
launcher and validator to `KinC_x36_4`, selects `Type=production` and
`target=LH2` from `config/nps_dvcs_all_kins_main.csv`, and rejects every input
run outside that exact set.

The configured set contains 57 runs. The current ALG-001 bundle represents 56;
run 6569 remains the only allowed missing-input exclusion. A different missing
run or any unexpected run is a hard error before model construction.

ALG-002B reads only ALG-001 `raw_observation` trees. It does not read an
efficiency, charge normalization, legacy purity weight, or cross-section
product. It does not write `pi0_weight` or modify a production tree. Output
under the canonical `output/` tree is refused.

## Implemented likelihood

Every selected event keeps its run, acquisition mode, multiplicity class,
invariant-mass bin, and accepted timing cell. The extended likelihood has six
nonnegative components per represented run and multiplicity stratum:

```text
pi0, true-coincidence combinatorial, horizontal, vertical, random, diagonal.
```

The pi0 and combinatorial components share the central-coincidence timing
density but have distinct mass densities. Horizontal and vertical components
share a mass density while retaining distinct timing densities. All timing
densities use the exact clipped ALG-002A timing support.

The signal candidates are a bin-integrated double-sided Crystal Ball and a
double Gaussian. The combinatorial candidates are the legacy logistic form and
positive cubic or quartic Bernstein densities. Other background components use
positive cubic Bernstein densities. Underflow, regular-range, and overflow
probabilities are normalized together.

Run peak-mean and log-width deviations use the balanced orthonormal zero-sum
basis already validated in ALG-002A. Population scales are fitted. Every
softmax mass shape fixes one reference logit, removing the otherwise exact
common-shift degeneracy.

At each shape point, the six run yields are profiled with a nonnegative
expectation-maximization update. Slow profiles use analytic-gradient L-BFGS-B
and SLSQP fallbacks and must satisfy a KKT residual below `2e-6`. The outer
SciPy fit alternates a mass block
with a complete ALG-002A timing-parameter block for each stratum. Timing is
therefore refitted rather than frozen. The existing ALG-002A covariance NPZ
files supply deterministic initial states and exact parameter/run-order checks.

## Interpreter and runtime

The system Python combines NumPy 1.26.4 with SciPy 1.9, which emits SciPy's
unsupported-version warning. The user-provided interpreter is the supported
execution environment:

```text
/group/nps/singhav/software/python/bin/python
Python 3.12.13
NumPy 2.5.0
SciPy 1.18.0
uproot 5.7.4
```

No additional package was required. On the full 76,242-observation bundle, a
cached objective evaluation took 0.22--0.32 seconds with numerical-library
threads capped at one, compared with 0.53--0.69 seconds under the system
interpreter. The mass block has 160 parameters and the two timing blocks have
252 parameters in total, so a converged multi-start campaign remains a batch
job even with the compatible interpreter.

## Output contract

One successful write creates a new, initially empty shadow directory with:

- `FIT_STATUS.txt`;
- `provenance.json` and `input_manifest.csv`;
- `run_coverage.csv` and `run_pi0_yields.csv`;
- `parameter_estimates.csv` with the complete ordered optimizer vectors and
  derived run quantities;
- `mass_timing_predictions.csv`;
- `conditional_intervals.csv`;
- `pi0_yield_covariance.npz`;
- `model_selection.json`.

Files are first written to a hidden sibling staging directory. The staging
directory is renamed to the requested path only after every artifact and the
status marker succeeds; an existing destination is refused.

The interval and covariance products are deliberately labeled
`SHAPES_FIXED_CONDITIONAL_NOT_COVERAGE_CALIBRATED`. They are diagnostics, not
the proposal's final profile intervals or complete cross-run replica
covariance. The model-selection file names the current candidates as
provisional until candidate comparison is complete.

`validate_joint_timing_mass_output.py` adds `validation_gates.json` and exits
with code 3 while any promotion gate fails or remains pending. It checks the
exact LH2 manifest, closure, optimizer finiteness and convergence, multi-start
agreement, boundaries, model selection, full covariance, toy coverage,
leave-one-run-out prediction, factorization, timing calibration, absence of a
pi0-weight product, and exclusion of efficiencies and cross sections.

## Validation completed

The focused regression
`tests/test_joint_timing_mass_background.py` verifies:

- production-LH2 filtering is case-insensitive only for the configured run
  type and exact for the target;
- LD2, fan-test, wrong-kinematic, arbitrary, and unapproved-missing runs fail;
- all six mass probabilities, including flow bins, are nonnegative and
  normalized;
- a two-run end-to-end fit closes exactly, writes only shadow products, and
  creates no `pi0_weight`.

The real-input read-only check found:

| Quantity | Result |
|---|---:|
| Configured production-LH2 runs | 57 |
| Represented runs | 56 |
| Missing runs | 6569 only |
| Observations | 76,242 |
| Strata | mode 2, multiplicity 2 and >=3 |
| Initial objective | 153817.237547 |
| Initial expected count | 76,242 |

A deliberately bounded real-data smoke optimization used one start, one
coordinate cycle, three mass iterations, and two timing iterations per
stratum. The final KKT-audited pass took about ten minutes under the compatible
managed Python. It produced:

| Diagnostic | Result |
|---|---:|
| Final smoke objective | 150879.268285 |
| Expected count | 76241.99999991708 |
| Absolute closure difference | 8.29e-8 |
| Run/stratum rows | 112 |
| Parameters at bounds | 4 timing population scales |
| Component yields at zero | 37 |
| Validator status | `NOT_PROMOTABLE` |

Both mass and timing blocks correctly report that they reached their small
iteration limits. All 112 inner profiles pass with maximum KKT residual
`2.58e-8`. The smoke pi0-yield sum is therefore not a scientific result
and must not be quoted or compared with the legacy analysis.

The smoke bundle is
`validation/runtime/alg002b_x36_4_joint_mass_smoke_v2/`. Runtime products are
ignored and are not part of the source commit.

## Reproduction

From the FINAL workspace, the focused regression is:

```bash
/group/nps/singhav/software/python/bin/python \
  tests/test_joint_timing_mass_background.py
```

The bounded real-data smoke command is:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/group/nps/singhav/software/python/bin/python \
  src/background_fit/run_joint_timing_mass_fit.py \
  --input 'validation/runtime/alg002a_x36_4_raw_all/KinC_x36_4/root/diagnostics_run*.root' \
  --config-csv config/nps_dvcs_all_kins_main.csv \
  --kin KinC_x36_4 \
  --output-dir validation/runtime/alg002b_x36_4_joint_mass_smoke_v2 \
  --timing-initial-dir validation/runtime/alg002a_x36_4_fit_all_v2 \
  --starts 1 --coordinate-cycles 1 \
  --mass-maxiter 3 --timing-refit-maxiter 2 --seed 20261001
```

The validator is:

```bash
/group/nps/singhav/software/python/bin/python \
  src/background_fit/validate_joint_timing_mass_output.py \
  validation/runtime/alg002b_x36_4_joint_mass_smoke_v2 \
  --config-csv config/nps_dvcs_all_kins_main.csv \
  --kin KinC_x36_4 --target LH2
```

Exit 3 is expected for the bounded smoke result.

## Open gates

The implementation is not promotable until the approved campaign completes:

1. converged dispersed starts and 20-start reproducibility;
2. calibrated treatment of yield and population-scale boundaries;
3. predeclared signal/background candidate comparison;
4. standalone high-statistics and leave-one-run-out checks;
5. mass/timing factorization checks;
6. completion of the ALG-002A timing calibration;
7. at least 2,000 coverage and misspecification replicas;
8. complete cross-run replica covariance and profile intervals.

ALG-002C event-weight construction remains unapproved and unimplemented.
