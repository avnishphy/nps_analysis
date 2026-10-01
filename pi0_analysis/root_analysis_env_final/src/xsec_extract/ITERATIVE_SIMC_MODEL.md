# Iterative SIMC pi0 model extraction

This path implements the yield ratio method of [Blok et al., Sec. V A,
Eq. (13)](https://arxiv.org/abs/0809.3161), with the local SIMC neutral-pion
model. The cited paper's charged-pion coefficients are not used.

## Weight and model

For each selected reconstructed exclusive event, the original input branches
give `base_weight = full_weight/sigcm` once, after finite-input and
`sigcm > 0` checks. At every trial parameter vector `p`,

```text
event_weight(p) = simc_yield_scale * base_weight * sigma_SIMC(vertex,p)
Ysim_b(p) = sum(reconstructed events in bin b) event_weight(p)
Vsim_b(p) = sum(reconstructed events in bin b) event_weight(p)^2
```

`simc_yield_scale` defaults to 1 and is applied once. The base retains the
original generation, luminosity, acceptance, radiative, flux and Jacobian
factors. No extra flux or bin-width factor is applied.

The parameterized implementation is in `simc_pi0_reweight.h`, with a small
model registry for future alternatives. Its default 34 coefficients exactly
match `/u/group/nps/singhav/simc_gfortran_updated/physics_pion.f`,
`sig_param_2021/exclfit` for `doing_pizero`: the arithmetic mean of the
pi+ and pi- fits, with `fpifact=0`. For each 17-coefficient fit,

```text
T  = p5/Q2 * exp(p6*Q2^2) / ((W^2)^p12 + W^p16) * exp(p14*|t|)
LT = p7/(1+p10*Q2) * exp(p8*|t|) * sin(theta*) / (W^2)^p13
TT = p9/(1+Q2) * exp(-7*|t|) * sin(theta*)^2
sigma = [T + sqrt(2*epsilon*(1+epsilon))*LT*cos(phi)
           + epsilon*TT*cos(2*phi)]
        * 8.539/(W^2-0.938^2)^2 / (2*3.1415928*10^6)
```

The source model's longitudinal term vanishes for this pi0 prescription;
that is a property of this provisional generator model, not an L/T
measurement. The model returns microbarn/MeV^2/radian. The separate SIMC
MAID mixture below W=2 GeV is unsupported; such selected events are counted
as rejected.

`sigcm` is evaluated at the SIMC interaction vertex. The producer copies
`Q2i` (GeV^2), `Wi` (GeV), `ti` (positive `-t`, GeV^2) and `phipqi`
(radians); the extractor evaluates at `t=-ti`. The updated producer also
copies the generator `epsilon` as `epsilon_i` when its input carries it.
Older files lack this branch, and some producer inputs can carry no epsilon,
leaving explicit `-999` sentinels. For those events epsilon is reconstructed
from `Q2i`, `Wi`, and configured beam energy. The fallback count, warning
and metadata state that reweighting is approximate. The default model's discrepancy from stored `sigcm` is
reported. Exact vertex reweighting requires a regenerated producer output
with `epsilon_i`, or another verified event-level vertex epsilon source.

## Fit and reporting

Data bin edges and yields are fixed before optimization. Physical `t`
bins use `t_bin_edges` from `xsec_config.h`. The `tprime_bin_edges`
span remains a phase-space cut and diagnostic axis; it does not set the
reported physical-t bins. Q2, xB, and phi edges are also read directly from
that header. Every event is assigned by reconstructed `Q2,xB,t,phi`.

By default, the shared selected-bin fit floats the actual Fortran
coefficients `plus.p5,plus.p7,plus.p9` (T, LT, TT amplitude controls).
`--model-free` selects other named coefficients; all others retain their
Fortran defaults. The Minuit2 objective is

```text
chi2(p) = sum_supported_bins
    [Ydata_b - Ysim_b(p)]^2 / [Vdata_b + Vsim_b(p)] .
```

`Vdata` is the data weight sumw2, scaled by the target correction squared.
This is a yield-space Gaussian approximation for weighted or subtracted
yields; it does not substitute Poisson errors. Nonfinite and nonpositive
trial event predictions receive a fit penalty. Nonconvergence is fatal.
The fit stores parameter errors and covariance where Minuit2 supplies
them, plus before/after chi2, nominal ndf, p-value and per-bin ratios.
A converged fit with p-value below 0.01 is marked
`converged_residual_mismatch`, so a poor model description is visible.
One parameter vector covers all selected bins. No separate L/T extraction is
claimed from a single epsilon setting.

For each physical-t/Q2/xB slice, accepted SIMC events are averaged with
immutable base weights to obtain one common `W_ref` and `xB_ref` for all
phi bins. `Q2_ref=xB_ref*(W_ref^2-Mp^2)/(1-xB_ref)`; directly averaged Q2
is a diagnostic. Epsilon follows from the derived reporting point. The
reported `t_center` and `phi_center` are geometric bin centers.
The final per-bin result is

```text
R_b = Ydata_b / Ysim_b(p_best)
sigma_data_b = R_b * sigma_SIMC(W_ref,xB_ref,t_center,phi_center,p_best).
```

Simulation yields still integrate **event-level** model weights.
`xsec_err` combines data sumw2 and final finite-MC sumw2 conditional on
the fitted model. Fitted-parameter errors are not added independently to
the same-data extraction error. Empty or unsupported bins carry NaN and
specific `bin_status` values (`empty_data`, `empty_sim`,
`unsupported_reference`, or `nonphysical_reporting_model`) in CSV. ROOT contains the model parameter
tree, covariance, reference tree, status, rejected-event counters, and
metadata. Reconstructed SIMC means in CSV/ROOT are recomputed with final
model weights; reporting W and xB retain immutable base-weight means. The CSV states cross-section units and before/after ratios.

`--normalize_mmiss` is rejected during iterative extraction: deriving an
extraction scale from the same data removes absolute model sensitivity. The
historical option remains available with `--fixed-default-model`.
Area-normalized global and missing-mass plot clones are display only; the
missing-mass shape plot explicitly uses original SIMC weights.

## Reproduce focused validation

From the repository root, load the required environment first:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; root-config --version'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 tests/compare_simc_fortran.py'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 tests/iterative_simc_closure.py /tmp/nps_iterative_closure'
```

Reproduced on 2026-09-25 with explicit temporary test bins: four physical points match the original Fortran
evaluator to 2.22e-16 relative difference. Synthetic truth
`(plus.p5,p7,p9)=(20.47049,1.044876,-222.29116)` fitted to
`(20.4705,1.04488,-222.291)` with chi2 6.6e-12/9. All 12 post-fit
ratios were within 9.5e-7 conditional standard deviations of unity; the
largest cross-section truth difference at reporting centers was 1.6e-8 relative.
Independent event sums agreed with output MC sumw2 to 3.1e-7 relative,
model center values to 1.3e-7, and Q2 references to 2.3e-10.
Two invalid `sigcm` values and one nonfinite `full_weight` were rejected;
one missing `epsilon_i` sentinel used the counted fixed-beam fallback.
The test also checks normalization conflict, optimizer failure, and empty
bin reporting.

The representative command actually run after the build above was:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; /tmp/nps_iterative_xsec --data-file output/KinC_x36_5_407/KinC_x36_5/root/combined_branches_LH2.root --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x36_5_407/smearing_output/KinC_x36_5_407/root/simc_pi0_analysis_output_smeared.root --out-dir /tmp/nps_iterative_representative_final2 --mmiss-lower 0.6 --mmiss-upper 1.1 --no-diagnostics --no-png --no-pdf'
```

A representative KinC_x36_5_407 extraction with two t bins and six phi
bins selected 138784 SIMC events. The default model gave chi2 2011.64;
the three-parameter fit gave 63.998/9 and p=2.26e-10. It completed and
was correctly marked `converged_residual_mismatch`. The older SIMC file
lacks `epsilon_i`, so all 138784 selected events used the fixed-beam
epsilon approximation; 703 differed from stored `sigcm` by over 1%,
with maximum 8.23%. These outputs are
diagnostic, not a validated publication cross section. Rerun the
simulation producer with `epsilon_i` and investigate remaining residual
structure before interpreting those fitted cross sections physically.
