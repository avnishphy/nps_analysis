# SIMC De-Modeling Convention

This workflow removes the SIMC event-level cross-section model from SIMC event
weights with:

```text
w_sim_base = full_weight / siglab
```

Use `siglab` for both exclusive and SIDIS pi0 samples.

## Reason

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

## Exclusive Pi0

For exclusive pion production, SIMC returns the lab five-fold cross section to
`event.f`:

```fortran
! /group/nps/singhav/simc_gfortran_updated/physics_pion.f:191
peepi = sigma_eepi*jacobian*(gtpr*fac)

! /group/nps/singhav/simc_gfortran_updated/event.f:1459
main%sigcc = peepi(vertex,main)
```

Since `siglab` is written from `main%sigcc`, it is the correct denominator for
exclusive de-modeling.

## SIDIS Pi0

For SIDIS, the returned cross section is `sigma_eepiX`:

```fortran
! /group/nps/singhav/simc_gfortran_updated/semi_physics.f:562
sigma_eepiX = sigsemi*jacobian/1.e6

! /group/nps/singhav/simc_gfortran_updated/semi_physics.f:584
sigma_eepiX = sigma_eepiX*fac

! /group/nps/singhav/simc_gfortran_updated/semi_physics.f:595
peepiX = sigma_eepiX

! /group/nps/singhav/simc_gfortran_updated/event.f:1521
main%sigcc = peepiX(vertex,vertex0,main,survivalprob,.FALSE.)
```

Again, `siglab` is written from `main%sigcc`, so `full_weight / siglab` removes
the event-level cross-section factor in `Weight`.

## Pipeline Usage

The active de-modeling consumers are configured to use `siglab`:

- `scripts/nps_sim_smearing_new_try.C`
- `scripts/excl_xsec_pi0_analysis_no_simc_model.C`
- `scripts/excl_xsec_pi0_analysis.C`

In the smearing fitter this is enforced by:

```cpp
// pi0_analysis/root_analysis_env/scripts/nps_sim_smearing_new_try.C
const char* SIM_MODEL_XSEC_BRANCH = "siglab";
double sim_base_weight = static_cast<double>(s_full_weight) / static_cast<double>(s_model_xsec);
TNamed("sim_weight_mode", "full_weight_over_siglab").Write();
```
