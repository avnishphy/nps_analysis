# ALG-002B simultaneous timing/mass shadow implementation

Date: 2026-10-01
Setting: `KinC_x36_4`
Run scope: configured `production` runs with target exactly `LH2`
Status: implemented shadow fitter; validation gates remain open
Implementation commit: `ec833cf614c0d3acbc48b0a57a9f3887acf58f28`
Legacy-comparison commit: `7575aa7d275f765102987a95484e3ee4673b0702`

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

## Relationship to the original timing-background method

The production macro still uses the original timing estimator and was not
changed. For each run it calculates

```text
N_acc = D + 0.5*(H + V) - 0.5*(F1 + F2),
```

where each control-box count is scaled to the nominal prompt-box area. It then
subtracts the resulting mass template before fitting the true-coincidence
combinatorial background.

ALG-002B does not use `N_acc`, its uncertainty, or the subtracted mass
histogram. It uses the same selected events and measured photon times, and its
horizontal, vertical, two-random, and diagonal components represent the same
physical timing-background mechanisms. Their shapes and normalizations are
instead inferred simultaneously over the accepted timing support.

`compare_legacy_timing_background.py` preserves the original estimator as an
independent validation comparator. It reads the authoritative
`accidental_est` and `accidental_err` values stored by production, reconstructs
the nominal box formula from raw region masks as a separate support diagnostic,
and compares both with the ALG-002B prompt-accidental prediction. Nothing from
that comparison is passed back into the fit.

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
pi0-weight product, exclusion of efficiencies and cross sections, and the
independent legacy-method comparison.

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

The independent comparison for this unconverged smoke result reports:

| Timing-background diagnostic | Setting sum |
|---|---:|
| Stored production legacy estimate | 1528.67 |
| Same nominal box formula on exported 139--161 ns support | 1771.72 |
| Raw-support formula minus stored legacy | 243.06 |
| ALG-002B prompt-accidental prediction | 2088.28 |
| Run-level ALG-002B/legacy correlation | 0.98945 |

The shifted production boxes extend to 139 and 161 ns, while the legacy timing
histogram spans 140--160 ns. The raw-support difference therefore confirms the
known clipping mismatch. The ALG-002B difference cannot be interpreted until
the outer fit converges and the shared-data covariance is obtained through
replicas.

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
  --nproc "$nproc" \
  --starts 1 --coordinate-cycles 1 \
  --mass-maxiter 3 --timing-refit-maxiter 2 --seed 20261001
```

After running the independent comparison command below, the validator is:

```bash
/group/nps/singhav/software/python/bin/python \
  src/background_fit/validate_joint_timing_mass_output.py \
  validation/runtime/alg002b_x36_4_joint_mass_smoke_v2 \
  --config-csv config/nps_dvcs_all_kins_main.csv \
  --kin KinC_x36_4 --target LH2 \
  --legacy-comparison \
    validation/runtime/alg002b_x36_4_joint_mass_smoke_v2_legacy_timing_v2
```

Exit 3 is expected for the bounded smoke result.

## 2026-10-01 legacy support and CPU update

The production legacy timing histogram now spans `[139,161)` ns with 220 bins,
retaining its 0.1 ns resolution. This includes the complete shifted waveform
sidebands. Diagnostics produced before this source change still contain the old
`[140,160)` histogram and must be regenerated before the updated legacy
comparison gate can pass.

The same validation exposed an older argument-order defect: the scalar legacy
estimator counted the first complete-accidental rectangle twice instead of
using the reflected `F2` rectangle. The helper contract and production call now
use the natural `(full2_t1, full2_t2)` order. On regenerated production-LH2 run
6407, the stored estimate is 30.111111 events and the independent raw-region
formula is 30.111111 events (difference `3.55e-15`); the stored histogram axes
are exactly 139 and 161 ns. This single-run check validates implementation, but
the complete 56-run diagnostics still need regeneration before comparison.

ALG-002B accepts `--nproc N`. With SciPy 1.16 or newer, L-BFGS-B distributes its
finite-difference objective calls over `N` forked CPU workers. Fork provides
copy-on-write access to the read-only fit dataset and avoids Python's thread
interpreter lock. Keep BLAS libraries at one thread as shown above to avoid
nested oversubscription. The configured worker count is stored in
`provenance.json` through the fit configuration.

On eight representative real-data objective calls, eight forked workers took
0.279 s versus 1.850 s serial, a 6.62x throughput gain with identical objective
values. A bounded 56-run production-LH2 CLI smoke with eight workers completed
and recorded `nproc=8`, `target=LH2`, and `lh2_only_enforced=true` in provenance.

The user subsequently authorized regeneration under the canonical setting
directory `output/KinC_x36_4/`. The analysis driver retains its fail-closed
default and requires the explicit
`--allow-canonical-raw-observation-export` option in addition to
`--raw-observation-export` and `--output-base <repo>/output`. The exception is
restricted to run-only diagnostics; it does not authorize combination,
efficiency, smearing, cross-section, or event-weight stages.

The authorized canonical regeneration completed on 2026-10-01 with eight
parallel ROOT jobs. It produced 56 diagnostic ROOT files, 56 per-run summaries,
and 57 logs under `output/KinC_x36_4/`; only input-missing run 6569 failed.
The new files contain 76,242 raw observations, and every `h_t1_t2` axis spans
exactly 139--161 ns with 0.1 ns bins.

`compare_raw_observation_inputs.py` hashes the exact Awkward buffers for every
branch and run. The regenerated bundle matches the earlier ALG-002A input for
all 56 runs and all 76,242 entries, so the existing timing initialization is
reusable. The machine-readable report is
`output/KinC_x36_4/validation/raw_observation_equivalence.json`.

The regenerated independent legacy comparison has exact setting closure:

| Diagnostic | Value |
|---|---:|
| Stored legacy accidental sum | 1771.722222 |
| Independent raw-mask formula sum | 1771.722222 |
| Raw-mask minus stored sum | 0.0 |
| Timing histogram ranges | 139--161 ns only |

An eight-worker, one-iteration canonical integration smoke at
`output/KinC_x36_4/alg002b/smoke_139_161_v1/` has objective 151259.147278 and
count-closure difference `5.59e-8`. It remains intentionally unconverged and
`NOT_PROMOTABLE`.

`run_joint_timing_mass_campaign.py` runs deterministic absolute start IDs as
independent atomic jobs, limits the total CPU budget, resumes matching complete
starts, and aggregates the approved 20-start thresholds. A two-start real-data
smoke used two concurrent starts with four gradient workers each. Both outputs
completed; rerunning the identical command reused them in two seconds. Its
one-iteration objective and yield spreads are deliberately nonpromotable.
The completed aggregate is supplied to the scientific validator with
`--start-campaign`; absent or nonpassing campaign evidence keeps the gate
pending or failed, respectively.

## 2026-10-02 converged autograd mass-block baseline

The 160-parameter finite-difference mass block remained iteration-bound after
300 steps. The managed analysis Python now installs the pinned, fork-safe
`autograd==1.9.1` dependency. The fitter differentiates the existing bounded
mass probability, penalty, and profiled-yield likelihood exactly; it does not
change the statistical model or output schema. Timing gradients retain the
validated forked SciPy workers.

Validation of the gradient backend found:

- NumPy/autograd mass probabilities agree within `2.36e-16`;
- mass penalties agree exactly;
- representative real gradients agree with central finite differences within
  `8.4e-7`;
- five real mass iterations took 2.78 s with autograd versus 28.98 s with
  16-worker finite differences, with objective difference `1.16e-4` after the
  same iteration count;
- warning-as-error mechanics and synthetic fit regressions pass.

With L-BFGS history 40, the continuation converged at
`output/KinC_x36_4/alg002b/autograd_probe_dscb_bernstein3_v3/`:

| Diagnostic | Value |
|---|---:|
| Objective | 147713.081597 |
| Final relative coordinate change | `9.42e-7` |
| Count-closure difference | `9.55e-6` |
| Maximum yield-profile KKT residual | `2.67e-7` |
| Sum of diagnostic pi0 yields | 34985.468217 |
| Optimizer convergence gate | pass |

The yield sum remains diagnostic until model choice, boundaries, replicas, and
coverage pass. The converged independent comparison reports legacy accidentals
1771.722222, joint prompt accidentals 2045.220786, difference 273.498564, and
run-level correlation 0.991843. The legacy value is not used by the fit.

The independent legacy comparison is:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/group/nps/singhav/software/python/bin/python \
  src/background_fit/compare_legacy_timing_background.py \
  validation/runtime/alg002b_x36_4_joint_mass_smoke_v2 \
  validation/runtime/alg002b_x36_4_joint_mass_smoke_v2_legacy_timing_v2 \
  --config-csv config/nps_dvcs_all_kins_main.csv \
  --kin KinC_x36_4
```

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
