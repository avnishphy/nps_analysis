# Scaled-Poisson extraction

Select this approximation explicitly with `--fit-objective scaled-poisson`.
The default `gaussian` objective retains the historical weighted least-squares
fit. The wrapper accepts the same option. Scaled mode requires continuous-angle
positivity; it enables positivity unless explicitly disabled, in which case it
fails.

The upstream `pi0_weight` correction and every other event-weight factor are
fixed inputs. No raw-count or background likelihood is constructed. For each
reconstructed row, `Y=sum(w)`, `V=sum(w²)`, and the fixed reference scale is
`s=V/Y` when `Y>0`. The fitted objective is

```text
D = sum_rows 2/s [mu - Y + Y log(Y/mu)]  (Y > 0)
D = sum_rows 2 mu/s                   (Y = 0)
mu = A theta, theta = (U, LT, TT) for every active truth block
```

The real-valued effective count `Y/s=Y²/V` is a diagnostic, not an integer
Poisson observation. The response `A` is held fixed. Its saved event second
moments are retained for later finite-SIMC inference. The old finite-MC
Gaussian variance never enters this objective.

For an empty supported row, `s` is borrowed in this order: the same
reconstructed t'/Q²/xB block pooled over phi if it has at least 20 effective
events; otherwise the same Q²/xB slice pooled over t' and phi if it has at
least 20; otherwise all selected data. Use `--scaled-empty-scale slice` or
`global` only for controlled sensitivity checks. Unsupported rows are
classified separately. A supported zero row remains in the objective and
contributes `2 mu/s`.

The optimizer is Minuit2 Migrad with `ErrorDef=1`. Each active truth block
is parameterized by signed LT and TT plus a nonnegative squared gap above the
exact continuous-angle minimum required U. Thus every trial cross section is
physical. Nonpositive reconstructed predictions are rejected. Interior errors
come from the transformed Minuit Hessian. At an active boundary, symmetric
Hessian covariance is unavailable as an interval; the transformed Hessian is
saved separately as a diagnostic. MINOS profiles are saved for published LT
and TT. A gap profile is not a U profile, so U intervals are not reported at
the boundary.

The output directory contains `scaled_poisson_rows.csv`,
`scaled_poisson_weight_groups.csv`, and
`scaled_poisson_profile_intervals.csv`. The ROOT file stores the signed
deviance-residual graph and fit metadata. A PNG/PDF residual plot is emitted
when diagnostics and the corresponding image format are enabled. Existing
weighted-yield and ratio diagnostics remain available. Residual-corrected
experimental points are post-fit central diagnostics; their historical
Gaussian analytic errors are not applied in scaled mode.

For `mu=Y+delta` with small `delta`,
`D=delta²/(sY)+O(delta³)`. Since `sY=V`, this approaches the historical
weighted Gaussian row term `delta²/V`.

Validation commands:

```bash
g++ -std=c++17 -O2 src/xsec_extract/tests/test_scaled_poisson_stat.cpp -o /tmp/test_scaled_stat
/tmp/test_scaled_stat
python src/xsec_extract/tests/validate_scaled_poisson_toys.py FIT_OUTPUT_DIR COMBINED_DATA_ROOT 500
```

The toy script tests one published U/LT/TT block using the saved forward
response and empirical selected-event weights, while holding the other truth
blocks fixed. It is a diagnostic for estimator shifts, not a full coverage
calibration. The scaled-Poisson objective approximates the compound-Poisson
distribution; upstream `pi0_weight` uncertainty, target/efficiency
systematics, and finite-SIMC response uncertainty remain separate.
