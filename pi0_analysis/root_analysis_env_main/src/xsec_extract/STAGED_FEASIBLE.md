# Opt-in constrained Gaussian central fit

`joint_minuit` remains the default. `staged_feasible` is an explicitly selected
diagnostic central-fit strategy for `simc_model`, `--positive-xsec`, and the
Gaussian objective. It rejects scaled-Poisson and positivity-off requests;
those continue to use the unchanged joint-Minuit path. No boundary confidence
errors or covariance are supplied by the staged strategy.

From this repository, after loading the Hall C/NPS environment:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash src/xsec_extract/run_xsec_pipeline.sh --mode simc_model --kin KinC_x36_4 --positive-xsec --fit-strategy staged_feasible --help'
```

Supply the same explicit input, configuration, and selection options as the
joint-Minuit extraction. The full recovered-input command, measured validation
results, scripts, and figures are in
`validation/xsec_staged_fit_20261004/REPORT.md` and `reproduce.sh`.

## Mathematical scope

For an exterior block the production constraint is the exact minimum of
`U + B LT z + epsilon TT (2 z^2 - 1)` on `[-1,1]`, where
`B = sqrt(2 epsilon (1+epsilon))` and epsilon is the block's production epsilon.
The independent slack representation `U = U_required(LT,TT) + S`, `S >= 0`,
is exact. The implemented Gaussian solver instead solves each three-variable
convex quadratic directly over that same closed cone.

Its boundary consists of two planar wedges and the curved surface

```text
(U,LT,TT) = t * (epsilon*(1+2*z*z), -4*epsilon*z/B, 1)
t >= 0, -1 <= z <= 1.
```

The planar wedges are spanned by `(3*epsilon,+/-4*epsilon/B,1)` and
`(epsilon,0,-1)` with nonnegative coefficients. Their small quadratic problems
are solved exactly. On the curved surface, eliminate `t` analytically and
isolate all stationary roots of a polynomial of degree at most five. There is
no angular grid or softened penalty. The unconstrained interior solution and
the cone apex are also candidates. Positive-definite block curvature makes
this enumeration a global block solve, up to floating-point precision.

Boundary construction uses `minimum_response()` itself and at most eight
`nextafter` steps on U to represent a nonnegative margin. This is an ulp-scale
rounding correction, not a positive physical floor. Zero slack is allowed.

Physics updates use an exact-Hessian Newton direction with numerical damping
of the search direction when needed, and backtracking that tests the original
event-level positivity before evaluating chi-square. Damping does not enter
the objective. Full row predictions are algebraically regrouped for speed;
the final objective is evaluated by `ProxyProblem::objective()`.

## Acceptance and outputs

Each cycle updates every nuisance block and refits physics. Two consecutive
complete cycles must satisfy all block-success flags, strict feasibility,
objective change, curvature-scaled parameter change, and constrained KKT
stationarity. With `r = min(1,sqrt(chi2))`, the thresholds are respectively
`1e-13 + 1e-9*r`, `1e-10 + 1e-6*r`, and `3e-11 + 3e-7*r`. Physics uses scaled
gradient tolerance `1e-12 + 1e-8*r`. The tighter absolute limits near zero
resolve injected synthetic boundary solutions. Minuit settings are unchanged.

Nonnegative KKT multipliers are fit to all active angular normals, including
both endpoints at nonsmooth intersections. The cone apex is checked against
the exact dual cone. Physics boundary solutions that do not pass the physics
gradient test are rejected; this is not a general optimizer for every possible
active physics constraint or every kinematic setting.

The existing finite-MC event variance and outer threshold are unchanged. With
no meaningful HESSE errors, parameter changes are divided by `abs(parameter)`
(with the existing `1e-30` numerical scale protection). This is a stricter
sufficient test than dividing by `max(abs(parameter), HESSE_error, 1e-30)`.

`staged_solver_history.csv` records all parameter updates and convergence
checks. `model_fit_strategy.txt` distinguishes constrained-solver status 0
from ROOT/Minuit status. Status 10 means invalid seed, 11 means cycle limit,
and 12 means the final exact production evaluator rejected the point. Call
limits and singular subproblems throw and fail extraction. Failed points are
never accepted. EDM and all covariance/error entries are unavailable (NaN).

The real-data tests support an opt-in central minimum, not a global-uniqueness
proof or a default-backend change. Profile chi-square, constrained toys, or
another boundary-aware uncertainty method needs separate validation.
