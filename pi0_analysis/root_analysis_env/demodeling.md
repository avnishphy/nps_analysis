# SIMC Weighting Conventions

The current workflow uses two related SIMC weighting conventions:

- Smearing calibration removes the full SIMC event-level cross-section factor:

```text
w_sim_base = full_weight / siglab
```

- The no-SIMC-model cross-section extraction fits a CM-response basis:

```text
w_xsec_base = full_weight / sigcm
```

Do not interchange these two denominators without also changing the fitted
quantity and the response factors.

## Common SIMC Factorization

The processed SIMC tree branch `full_weight` is built from the input SIMC
`Weight` branch multiplied by the run/sample normalization in
`scripts/simc_pi0_analysis.C`:

```cpp
// pi0_analysis/root_analysis_env/scripts/simc_pi0_analysis.C
out_full_weight = Weight * weight_factor;
```

The processed tree writes `siglab` as the cross-section branch used by the
downstream de-modeling code:

```cpp
// pi0_analysis/root_analysis_env/scripts/simc_pi0_analysis.C
outtree->Branch("siglab", &out_siglab, "siglab/F");
out_siglab = in_siglab;
```

In SIMC, `Weight` is written from `main%weight`. The cross-section factor inside
`main%weight` is `main%sigcc`:

```fortran
! /group/nps/singhav/simc_gfortran_updated/event.f:1558
main%weight = main%SF_weight*main%jacobian*main%gen_weight*main%sigcc
```

The ntuple branch `siglab` is also written from `main%sigcc`:

```fortran
! /group/nps/singhav/simc_gfortran_updated/results_write.f:144
ntu(44) = main%sigcc              ! d5sig

! /group/nps/singhav/simc_gfortran_updated/results_write.f:203
ntu(39) = main%sigcc              ! d5sig
```

Therefore `full_weight / siglab` removes the same event-level cross-section
model factor that SIMC multiplied into `Weight`, while preserving generator,
phase-space, acceptance, and normalization factors.

Use `siglab` for smearing de-modeling across exclusive, SIDIS pi0, and delta
samples. The SIMC source assigns `main%sigcc` by channel:

```fortran
! /u/group/nps/singhav/simc_gfortran_updated/event.f
doing_pion  -> main%sigcc = peepi(vertex,main)
doing_delta -> main%sigcc = peedelta(vertex,main)
doing_semi  -> main%sigcc = peepiX(vertex,vertex0,main,survivalprob,.FALSE.)
```

The ntuple writes that same `main%sigcc` as `siglab`:

```fortran
! pion/kaon/delta ntuple
ntu(44) = main%sigcc              ! d5sig

! semi/rho ntuple
ntu(39) = main%sigcc              ! d5sig
```

## Smearing

Smearing should fit detector response, not the SIMC model cross-section shape.
The active smearing fitter is configured as:

```cpp
// scripts/nps_sim_smearing_new_try.C
const bool USE_SIM_MODEL_XSEC_DEMODELING = true;
const char* SIM_MODEL_XSEC_BRANCH = "siglab";
```

The smearing event weight is therefore:

```text
w_sim = (full_weight / siglab) * is_exclusive_ellipse
```

This is the correct common choice for exclusive, SIDIS, and delta rows because
`siglab` is the exact `main%sigcc` factor in `Weight`.

## Cross-Section Extraction

The no-SIMC-model exclusive cross-section macro intentionally uses:

```text
w_xsec_base = full_weight / sigcm
```

For exclusive pion production, SIMC forms:

```fortran
! /u/group/nps/singhav/simc_gfortran_updated/physics_pion.f
ntup%sigcm = sigma_eepi
peepi = sigma_eepi*jacobian*(gtpr*fac)

! /u/group/nps/singhav/simc_gfortran_updated/event.f
main%sigcc = peepi(vertex,main)
```

Thus:

```text
full_weight / sigcm
  = normalization * SF_weight * jacobian * gen_weight
    * (siglab / sigcm)
  = normalization * SF_weight * jacobian * gen_weight
    * davejac * gtpr * fac
```

This removes the SIMC CM model cross section while keeping SIMC's lab/CM
Jacobian, flux convention, and Fermi-motion factor in the response matrix. That
is why the cross-section macro must not multiply an additional gamma/flux factor
on top of this basis.

## Current Consumers

- `scripts/nps_sim_smearing_new_try.C`: `full_weight/siglab`
- `scripts/excl_xsec_pi0_analysis_no_simc_model.C`: `full_weight/sigcm`
- `scripts/excl_xsec_pi0_analysis.C`: model-ratio legacy path; check
  `model_xsec_mode` before using outputs

The practical rule is: use `siglab` when removing the SIMC model from detector
smearing weights; use `sigcm` only for the exclusive CM-response cross-section
fit that is designed around that basis.
