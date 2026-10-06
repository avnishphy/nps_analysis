# Provisional SigParam2021-inspired neutral-pion model

`--mode simc_model` uses one production model, `sigparam2021_pi0`:

- Component-wise average of the pi+ and pi- parameterizations.
- Charged-pion pole removed from L.
- Smooth non-pole neutral-pion L placeholder retained.
- Global U, LT and TT normalizations, plus exactly one pivoted U slope fitted.

The physics parameters are `(N_U,DeltaB_U,N_LT,N_TT)`, followed by `(eta_U,eta_LT,eta_TT)`
for each populated exterior/excluded-group truth block in ascending active-block
order. The model is `U=N_U*exp(-DeltaB_U*(tau-tau0))*(Tbar+epsilon_vertex*Lbar)`,
`LT=N_LT*LTbar`, `TT=N_TT*TTbar`. Tau is -tprime; tau0 is the positive
model-independent-response-weighted mean over accepted physical events. It is
fixed throughout each fit, recorded in model_context.csv and model_fit.root.
Positive DeltaB_U suppresses large tau relative to the pivot. Bounds are
[-20,20] GeV^-2; the sign is not restricted. No LT/TT slopes or independent T/L
normalizations are fitted. Initial values are `(1,0,1,1)`; NU has lower bound
1e-12, while NLT/NTT are signed and unbounded. Six starts cover both interference
signs, modest NU variants, and slope seeds 0,+0.5,-0.5. Use
`--model-initial NU,DeltaB_U,NLT,NTT`, `--model-starts`
and existing minimizer controls. The removed `--model` selector is rejected.
JSON `model_identifier` is identifying metadata, not a physics-model selector.
`--model-fix-u-slope` fixes DeltaB_U=0 for three-normalization validation. The
four-slot vector retains a zero covariance row/column for the fixed slope;
nominal DOF and conditioning use only free parameters. `--model-before-dir`
adds a diagnostic-only comparison with a prior model output using identical data
and response. Neither comparison points nor ratios enter the objective.

`xsec_sigparam2021_pi0_model.h` contains the unchanged pp/pm tables and baseline.
The original Fortran is `/u/group/nps/singhav/simc_gfortran_updated/physics_pion.f`.
Its `s_gev` is W2, not W. Neutral L retains `(Q2/1 GeV2)*G(Q2)^2`, where
`G=1/(1+p1*Q2+p2*Q2^2)` is empirical damping, not a pi0 electromagnetic or
transition form factor. The charged pole is absent; L is not forced zero.
Both charged sets are evaluated before averaging components. Every component
receives `8.539/(W2-0.938^2)^2 / 1e6` once. Components are ub/MeV2; only the
detector response contains `1/(2*pi)`.

Accepted generated events are cached and individually folded.
`vertex_epsilon_from_exclusive_simc` uses matched raw events, never smeared
`epsilon_i`. Generated theta_cm comes from exclusive two-body invariants.
The U Jacobian is `(exp(-DeltaB_U*(tau-tau0))*U0, -(tau-tau0)*U)`;
LT/TT derivatives are their baseline values in columns 2/3. Full
detector-induced parameter/nuisance covariance is retained. Weighted Gaussian,
finite-MC sumw2 iteration, scaled-Poisson fixed-response treatment, row masks,
event matching and normalization are unchanged. Negative values are not clamped.

| Product | Interpretation |
| --- | --- |
| model_parameters.csv, model_fit.root | Four physics parameters plus nuisances; full covariance |
| model_context.csv | Fixed pivot, weighting, fixed-slope flag, synthetic flag, unit contract |
| model_U_correction.csv | Dimensionless U shape factor over accepted tau |
| model_shape_comparison.csv, model_before_after.json | Optional before/after residuals, ratios and descriptive slopes |
| model_reconstructed_yields.csv | Exact event-folded means, components, residuals and pulls |
| model_event_cache.csv, model_row_jacobian.csv | Reproducible event fold and covariance |
| model_structure_functions.csv | Response-base-weighted event averages by truth block |
| model_structure_curves.csv | Fixed generated weighted mean Q2/W2/epsilon context, declared in baseline diagnostics |
| model_tprime_shape.csv | Residual sums, mean pulls and objective contributions versus reconstructed tprime |
| model_shape_ratios.csv, model_shape_summary.csv | No-model/model ratios, ranges and sign discrepancies |
| model_charged_spread.csv, model_charged_spread_summary.csv | D=(plus-minus)/[0.5*(abs(plus)+abs(minus))] at accepted events; not an uncertainty |
| model_longitudinal_summary.csv | MODEL-DEPENDENT T/L DECOMPOSITION and L fraction |
| model_event_positivity.csv | Analytic full-phi minimum and 720-angle scan before/after normalization |

Ratios never enter the objective. Denominators below 0.001 of the component's
largest bin magnitude are flagged and omitted. Ratio bars show numerator
uncertainty with the model denominator fixed; shared-data cross-covariance is
not evaluated, so they are not independent pulls or significance tests. The
legacy constant-cell residual-corrected experimental points remain explicitly
unavailable for the varying-event model.

U gains one slope; LT and TT differ only by normalization. Detector phi plots, row
residuals, tprime summaries and no-model ratios expose shape deficiencies.
Future changes should respond to observed discrepancies with minimal
physics-motivated modifications, not automatically add flexibility.

All structure-function plots use nb/GeV2 (internal microbarn/MeV2 multiplied by
1e9, with no extra 2pi). Existing CSV physics units remain microbarn/MeV2.
Event-averaged model points are primary; reference-kinematic curves are thin
dashed diagnostics. Shading marks the accepted generated range; the exact forward
limit outside it is extrapolation. Synthetic headers derive from the explicit
ROOT `fixture_provenance` metadata, never the input filename. Real data without
that synthetic metadata are not labeled synthetic.

This remains provisional: the baseline is charged-pion-inspired, averaging is
heuristic, L is phenomenological, and T/L are not experimentally separated.
W<2 GeV is unsupported. Mathematical guards do not certify a validity domain.

Reproduce checks in a fresh directory:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash src/xsec_extract/tests/run_sigparam_validation.sh src/xsec_extract/xsec_config/xsec_config_x36_4.json /tmp/sigparam4_check'
```

This compiles both entries and tests, checks charged Fortran agreement,
normalization/nuisance closure, matching guards, empty rows, both objectives,
and exported event replay. Full pipeline commands and the historical comparison
are in `validation/xsec_uslope_20261004/REPORT.md`. Historical artifacts are
preserved and are not selectable production models.
