# Physics audit: demodelled smearing and pi0 cross sections

Date: 2026-09-11. Scope: the requested `root_analysis_env_main/src`, `simc_gfortran_updated`, and the active `geant4_simc` application. No analysis code, efficiencies, or production outputs were changed. This is a source-level physics audit, supported by selected SIMC normalization records and analytic mechanism checks; it is not an experimental closure result.

**Assessment:** the inverse-model response approach is viable, and its central weight algebra is consistent. The present implementation cannot yet establish an artifact-free detector response or a validated Born cross section. Several selections and response operations differ between fitting, production, and extraction. Those differences should be resolved before interpreting disagreements as a need to change the physics model.

## 1. Scope and current configuration

Source keys used below; numbers following a key are source line numbers at audit time:

| Key | File |
|---|---|
| S | [Smearing fitter](../src/simulation_smearing/nps_sim_smearing_new.C) |
| P | [Simulation producer](../src/simulation_smearing/simc_pi0_analysis.C) |
| X | [Cross-section extraction](../src/xsec_extract/excl_xsec_pi0_analysis_no_simc_model.C) |
| D | [Data analysis](../src/analysis/nps_analysis_main.C) |
| C | [Run combination](../src/analysis/combine_analysis_branches.py) |
| A | [Acceptance configuration](../config/acceptance_cuts.conf) |
| F | `/group/nps/singhav/simc_gfortran_updated/` |
| G | `/work/hallc/nps/singhav/geant4_simc/HallC_NPS/DVCS_evt_gen/DVCS/src/` |
| N | `/work/hallc/nps/singhav/geant4_simc/NPS_SOFT/` |

The active chain is SIMC weighted events -> Geant4 block deposits -> reconstructed clusters -> analysis acceptance/pair selection -> additional response -> cross-section fit. The producer uses `clust_E/X/Y`, not individual truth-photon deposit branches. Unused Geant4 variants and the efficiency implementation were outside this investigation.

| Item | Current source |
|---|---|
| Smearing demodelling | **Off**, `USE_SIM_MODEL_XSEC_DEMODELING=false`; optional divisor is `siglab` (S:484-490) |
| Extraction demodelling | `full_weight/sigcm` (X:890) |
| Smearing observables | `Mgg` and `(p_target+gamma gamma)^2`, weights 1 and 1; missing-mass weight 0 (S:248-250) |
| Photon response | `Emean=a+bE+c ln(E/1 GeV)`, Gaussian extra width `sigma*sqrt(Emean)` |
| Position smearing | Fitter **on** (S:372), producer **off** (P:111) |
| Calibration selection | Independently determined data and unsmeared-MC ellipses |
| Extraction selection | Data reconstructed exclusivity; MC generated-channel label |

The user confirmed this is the intended source and that the purpose includes comparing `sigcm` with `siglab`. The user also identified `tgt_contam=0.584 +/- 0.014` as the yield-reduction factor from helium contamination of the hydrogen target; division is intended to recover the lost hydrogen yield. The user requested documenting the inherited hadronic radiation treatment as unresolved; no validated prescription was supplied. Existing output files may predate current source; no claim is made that every archived production used these settings.

## 2. What is already physically consistent

**Weight factorization.** For stationary-proton exclusive production, the inspected code has

```text
Weight = generation/radiation weight * angular Jacobian * siglab
siglab = sigcm * hadron Jacobian * virtual-photon flux * Fermi-flux factor
full_weight = normfac/Ngen * Weight
```

The Fermi factors reduce appropriately for hydrogen. `sigcm` is the vertex hadronic `d^2sigma/(dt dphi)` in microbarn/MeV^2/radian, not `d sigma/dOmega_cm` (F/physics_pion.f:164-191; F/event.f:1558). Thus `full_weight/sigcm` retains the flux, radiative weight, and transformation of integration variables needed for a unit hadronic cross section. The extractor correctly avoids multiplying another flux (X:1106-1111).

The `1/(2pi)` Fourier convention is consistent with phi-integrated U, LT and TT structure functions. Numerical units `microbarn/MeV^2 = barn/GeV^2` explain the existing plot label; conversion to `nb/GeV^2` requires multiplication by `10^9`. No additional bin-width divisor belongs in an already integrated Monte Carlo response. At fixed Q2,xB, `t'=t-tmin` has unit Jacobian.

**Charge normalization.** C:595-596 uses `PS/[Qrun_mC * LT * efficiency]`; X:828-832 multiplies by `Qrun/Qtotal`. Together these produce corrected counts/mC. This algebra is sound; it does not validate every efficiency input. SIMC's `normfac=L*genvol*Naccepted/Ntried` gives the appropriate accepted-event factor `normfac/Naccepted`, provided the histogram and event population match (F/simc.f:368-396).

The actual staged `nps_simc_20260824_135058/.../nps_excl_pi0_x60_4b.hist` contains `Ntried=33694137`, requested/contributing events `1000000`, charge `1 mC`, `L=2.67882e9 /microbarn`, `genvol=2.67137`, and `normfac=2.12385e8`: event factor **212.385**. This confirms the normalization convention for that sample, not its complete Geant4 processing history.

**Detector foundation.** Geant4 transports the two photons from the SIMC vertex, preserves scalar event weights, and normally writes zero/one-cluster events too (G/PrimaryGeneratorAction.cc:188-234; G/EventAction.cc:598,650). The isotropic two-photon decay is appropriate for spin-zero pi0 (F/pizero_decay.f:46-74). Additional smearing can represent missing detector fluctuations; it must represent the residual beyond Geant4's existing shower fluctuations and leakage.

## 3. Findings and remedies

### 3.1 Demodelling does not isolate detector resolution

**Established distinction:** optional `full_weight/siglab` removes the hadronic model **and** virtual-photon flux and hadron Jacobian. It describes a unit lab cross section under the generator/radiative measure; it does not flatten Q2,xB,t,phi. `full_weight/sigcm` defines a different measure. Neither produces a detector-only event population.

The current fit pools mass distributions across each section. In particular,

```text
M(p+gamma gamma)^2 = Mp^2 + Mgg^2 + 2 Mp (E1+E2).
```

Its shape depends strongly on the physical pion-energy population. Changing model weights, channel mixture, or accepted kinematics can change this shape without changing the detector. Fitting that difference with `a,b,c,sigma` can distort energy scale and resolution. Merely switching the divisor from `siglab` to `sigcm` does not solve this identifiability problem. Current weighted fitting is also sensitive to an inaccurate model population.

**Remedy:** calibrate conditionally using independently constrained electron kinematics/energy proxies and geometry; or include smooth physics-population parameters as separate nuisances. If binning in reconstructed photon energy, model the response-dependent selection/truncation: narrow bins alone do not establish identifiability. Hold the detector response fixed while varying plausible truth populations and require recovery of the same response. Keep missing mass as a held-out diagnostic under current zero-weight configuration. Do not treat matching the energy-sensitive distribution as an independent calibration proof.

### 3.2 The response fitted is not the response subsequently applied

**Established mismatches:** position smearing is enabled in S and suppressed in P. Any fitted nonzero angular broadening is therefore absent from the final sample. Fitter summaries use section coefficients first and interpolate only as fallback (S:7270-7345); P prefers interpolated maps (P:1484-1536).

The interpolation additionally treats coefficients stored at section centers as though anchored at lower edges: S:2650-2675 uses `(x-xmin)/cell_width` without the center offset. At a first-section center, adjacent coefficients `b=(1.00,1.10)` give `b=1.05`, not the fitted `1.00`. This is a concrete mechanism for shifting spatial calibration patterns. It vanishes for a constant map. Unfitted cells start with default response, including `sigma=0.05` (S:2722-2726), so interpolation can also blend an assumed width into fitted regions.

**Remedy:** use one shared response evaluator and configuration in fitter and producer. First apply exactly the fitted section response; introduce interpolation only with center-consistent coordinates, explicit treatment of unmeasured cells, and a fresh final-response validation. Fit or validate the interpolated representation that will actually be used. Allow an identity/null response: the fitter's `sigma>=0.01` currently excludes zero extra energy broadening (S:338). Do not force extra width where Geant4 already matches or exceeds the data width.

### 3.3 Selection must act on the final reconstructed event

**Established:** P:2181-2334 applies cluster energy/fiducial cuts, retains up to four good clusters, chooses a pair, and only then smears that pair. P:2336-2438 writes the result without reapplying those cuts or reconsidering the pair. Consequently it retains events smeared below threshold, cannot admit initially rejected clusters migrating upward, and cannot change pair identity. Position smearing would similarly affect fiducial migration once enabled.

**Remedy:** start from all available reconstructed Geant4 clusters, apply the residual response, then run the data-equivalent cluster cuts, masking, multiplicity handling and pairing. Retain a sufficiently loose precursor sample. Simply adding a final cut to the existing preselected two-cluster tree fixes outward losses but cannot restore inward migration. Check the effective threshold and photon-separation dependence, not only mass peaks.

### 3.4 Data and MC do not yet implement the same acceptance

**Established:** D:1541-1547 removes dead-region clusters before pairing. P:2293-2302 only records `passes_dead_block_mask`; neither S nor X consumes it. This can create spatial deficits in data relative to MC, which correlate with t and phi. A single representative-run mask would also be insufficient if the combined runs have different masks.

**Established:** X requires `is_exclusive` in both samples, but in data it means the corrected missing-mass cut (D:3307-3317), whereas in MC it means the generated exclusive channel (P:2142 and 2428). An exclusive generator still has reconstructed radiative and resolution tails; channel identity does not implement the data cut. The configured `KinC_x60_4b` data vertex range is also +/-7 cm, while simulation inherits +/-8 cm (A).

**Remedy:** apply the same reconstructed exclusivity transformation and cut to signal MC, with identical fiducial/vertex cuts unless a separate efficiency correction is explicitly demonstrated. Apply run-dependent masks before MC pairing and combine predictions with measured run luminosity/efficiency conventions. Keep generated channel labels separate from reconstructed selection. This has priority over tuning a model to match data/MC yields.

### 3.5 Calibration ellipses condition the distributions being fitted

**Established:** S:477-479,4493-4500,4625-4645 preselects data using its combined ellipse and MC using its own unsmeared ellipse. P:915-965 trains the MC ellipse on all three model-weighted channels. Those event gates are then fixed during the smearing optimization; final MC receives a newly trained ellipse. Thus calibration fits truncated, differently selected populations. Resolution tails removed before the fit cannot be recovered by broadening the retained events. Changes to the channel model can also change the training gate.

**Remedy:** use a broad common calibration domain and explicitly model backgrounds; alternatively define a calibrated common gate and include its response-dependent efficiency by re-evaluating it after every trial response. Train any data-derived selection independently and freeze its definition for the measurement. Validate tails outside the calibration gate and the actual final exclusivity efficiency. Matching two separately fitted ellipses is not acceptance closure.

### 3.6 Mass-purity weights do not guarantee differential signal yields

**Established estimator:** D:3191-3216 subtracts accidentals and a combinatorial fit in Mgg, then sets `pi0_weight(Mgg)=Nsignal/Nall`. D:3311-3317 retains positive weights and applies the corrected-missing-mass gate. This can reproduce a marginal mass yield without reproducing the signal distribution in t,phi,Q2 or in the selected missing-mass subset: background and signal need not have identical conditional kinematics at fixed Mgg. Clipping negative signal estimates can additionally bias sparse bins upward.

**Remedy:** subtract/fit backgrounds within the final physics slices or use a joint discriminant/kinematic model, with correlations and subtraction uncertainties. Test the current estimator using injected backgrounds whose phi,t and missing-mass correlations differ from signal. Preserve signed estimates if using subtraction. This is a validation requirement, not a measured numerical bias in the present data.

### 3.7 The smearing objective is a useful discrepancy, not the stated count likelihood

**Established:** S:470,1475-1496 applies a Poisson/Baker-Cousins deviance to efficiency/charge/purity-weighted yields, ignoring their weight variances. These are not independent Poisson counts. The two mass observables share events; section samples share photons/events too. Hence the reported objective/ndf and HESSE errors do not have their nominal independent-count interpretation. Rescaling yield units by 1000 rescales this deviance by 1000, affecting absolute quality thresholds and inferred errors without changing the physics.

S:1990-2048 separately normalizes each observable using its unsmeared in-range integral. This preserves a penalty for response moving events across histogram boundaries, but it is neither a freely normalized shape fit nor a full common-rate likelihood.

**Remedy:** choose the estimand explicitly: use a joint raw-count likelihood with background/efficiency nuisances, or a covariance-aware fit to weighted distributions. For a shape-only fit, profile the appropriate normalization including truncation; for a rate fit, use one physical normalization. Bootstrap original events/runs and background fits. The 80 smearing replicas are integration draws, not 80 independent generated events; retain their common-event correlation. Include effective sample size when judging sparsely populated sections.

### 3.8 Extraction needs a vertex-to-reconstruction response

**Established approximation:** X:871-900 assigns events to reconstructed Q2,xB,t',phi bins and evaluates the Fourier harmonics at reconstructed phi. Independent fits in each reconstructed slice do not associate migrated events with their originating physics bin. Dividing out the original vertex model does not undo that migration.

For example, true bin cross sections `(1,3)` with symmetric 20% migration produce reconstructed yields `(1.4,2.6)`. Dividing those by a flat-model response returns `(1.4,2.6)`, not the vertex cross sections. This is an analytic example, not an estimate of NPS migration.

X:1037-1052,1103-1140 also uses epsilon evaluated at one slice mean outside the event sum. Generally `sum(w*epsilon*cos(2phi))` differs from `epsilon(mean)*sum(w*cos(2phi))`.

**Remedy:** retain vertex and reconstructed quantities together; bin the predicted yield by reconstructed coordinates but evaluate physics coefficients/harmonics and epsilon at the interaction vertex. Include off-diagonal migration and feed-in from beyond the reported range. Fit smooth vertex functions or a truth-bin response matrix. Use a diagonal approximation only after independent sloped/modulated closure shows its error is acceptable. At one epsilon setting, report U=T+epsilon L; do not infer independent T and L without additional information. Define how any epsilon variation transports U to its quoted reference point.

### 3.9 Helicity handling currently affects even the unpolarized fit

**Established:** P:1438 assigns alternating MC helicity signs. Since both trees have helicity branches, X:1112 activates four parameters. X:1144-1174 then fits total, plus and minus yields as independent observations although the total contains the other two. That double-counts information. Alternating MC signs also do not encode measured helicity luminosities; the basis lacks beam polarization.

**Remedy:** begin with the three unpolarized terms after checking that helicity-odd contributions cancel through balanced/corrected exposures or are negligible. For LT', fit only independent plus/minus observations with measured exposures, beam polarization and the verified phi/helicity sign convention. Use totals for display; retaining this redundant observation requires explicit handling of the singular covariance. Propagate finite-MC, background, response and normalization uncertainties; X currently uses conditional data sumw2, and its fitted-prediction variance is not MC counting variance (X:1238-1243).

### 3.10 Neutral-pion radiation needs a dedicated physics decision

**Established implementation; physical validity unresolved:** configured `using_rad=1, one_tail=0` enables three tails. The staged x60_4b and x36_5 histograms report nonzero third-tail parameters (`lambda3=0.016`, central proposal fractions 0.126 and 0.130). These fractions are not cross-section corrections.

F/init.f:693-711 supplies `vertex%p` (the generated pion leg) to the hadronic radiator and virtual-correction call. F/radc.f:768-801 implements the inherited charged-hadron/proton-tail expression; :495 removes third-tail energy from that pion leg. Pi0 photons are generated from the vertex pion earlier (F/event.f:899-900). A neutral pion is not a charged radiating external leg; the recoil proton is a different particle. This treatment therefore needs explicit justification for this channel. Demodelling retains it, so it cannot fix a wrong radiator assignment.

**Remedy:** validate an exclusive-pi0 radiative prescription, including which charged legs/interference terms it retains. Compare the present configuration with a consistent electron-tail-only diagnostic (`one_tail=-3`) and an appropriate exclusive-pi0 calculation; do not promote that flag change directly to the final correction. Check missing-mass-cut dependence and regenerated detector yields if the event kinematics change. [EXCLURAD's documented formalism](https://www.jlab.org/RC/exclurad/) provides an independent non-peaking internal-radiation reference, but its physics model, kinematic coverage and external-radiation treatment must match the intended comparison.

### 3.11 Geant4 response and absolute normalization need explicit boundaries

G/EventAction.cc:278 scales block deposits by **1.04**, then calls `TriggerSim(0.5)` before clustering. N/TCaloEvent.cxx:244-343 removes blocks not belonging to passing overlapping 2x2 sums. This changes shower energy and clustering as well as event efficiency; it is not merely an event trigger. Blocks have no simulated timing spread in this path. A post-cluster response cannot restore deleted energy or missing clusters.

**Remedy:** establish equivalence with the data readout, suppression and clustering thresholds, or demonstrate the selected sample lies on the relevant efficiency plateau. Compare zero/one/two/multiple-cluster fractions and photon-separation dependence. Match vertex/shower-depth conventions: N/TCaloCluster.cxx:232-305 uses uncorrected logarithmic centroids, while both analysis helpers reconstruct photon directions from the nominal origin. This common approximation needs residual checks versus vertex and energy, not an automatic declaration of error.

SIMC already requires both photons inside its calorimeter envelope before Geant4 (F/simc.f:1499-1537). That does not double-count acceptance when Ntried is retained, but can remove potential inward migration. Widen the generator envelope with final cuts fixed and require convergence. Check Geant4 processing completeness: a requested first-N subset is permitted; normalization cannot silently assume the entire SIMC sample was transported.

Forced gamma-gamma decay has unit branching probability; no downstream branching factor was found in the audited chain. Apply the physical branching ratio once for a pi0-production cross section, or explicitly report the gamma-gamma-channel cross section.

**Helium contamination clarified by the user:** `0.584` represents reduction of hydrogen yield, so dividing by it is physically consistent when SIMC uses nominal uncontaminated hydrogen luminosity. This is not an inverted signal-purity correction. The inspected x60_4b input uses `targ%abundancy=100`, `rho=0.0723 g/cm^3` and `thick=717.80 mg/cm^2` (input:38-41). Verify the correction's measured run/target-period coverage and that the adopted thickness has not already incorporated the same reduction. Its 2.40% relative uncertainty is a correlated normalization contribution. Correcting hydrogen luminosity does not itself subtract helium-induced events or establish the mixture's energy-loss/radiation thickness; check those separately where relevant. The yield-loss factor alone should not be interpreted as a helium atom or mass fraction.

## 4. Conventional Hall C extraction: concrete implementation path

The conventional prescription is

```text
sigma_exp(v0) = [Ydata_i / YMC_i(model)] * sigma_model(v0),
```

where both yields use the same reconstructed acceptance and charge convention, and the model is evaluated at a specified Born reference point. The model supplies the finite-acceptance/bin-centering dependence; iteration tests its adequacy. This approach and model-variation studies are described in [Basnet et al., Phys. Rev. C 100, 065204, Eq. (7) and Sec. IV](https://link.aps.org/accepted/10.1103/PhysRevC.100.065204). That paper supplies the method, not a pi0 model for these kinematics.

### A. Establish the common response first

Resolve Sections 3.2-3.5 and the radiation prescription. Freeze detector response, selection, run exposures and background treatment while fitting the cross-section model. A physics-model update must not compensate for an acceptance hole or a different smearing operator.

### B. Add a dedicated pi0 hadronic model

The current high-W pi0 prior is an average of charged-pion parameterizations with the longitudinal pion-pole contribution disabled (`fpifact=0`; F/physics_pion.f:738-825). This is explicitly a provisional prescription; it does **not** demonstrate physical `sigma_L(pi0)=0`.

Add a pi0-specific evaluator at `peepi`'s model dispatch (F/physics_pion.f:126) or within the pi0 branch of `sig_param_2021`. Preserve charged-pion behavior and explicitly handle the low-W transition. Fit a compact, smooth parameterization of U, LT and TT versus Q2,W,t, with nonnegative total cross section over supported phi and appropriate forward-angle behavior. Introduce T/L separation only with sufficient epsilon coverage or stated external assumptions.

Keep flux, Jacobians and radiation outside this hadronic evaluator. A model supplied in `nb/GeV^2/rad` must return `sigcm` in `microbarn/MeV^2/rad` by multiplying by `10^-9`. Preserve the `1/(2pi)` and phi-sign conventions exactly.

### C. Preserve sufficient event information

Carry through `(sample_id,event_id)`, original weights/model version, vertex `Q2,W,t,phi,epsilon,Ein`, and reconstructed variables/clusters. SIMC already exports `Q2i,Wi,ti,phipqi` (F/results_write.f:149-153), but P does not retain these in its reduced output. Explicit vertex epsilon/beam energy are not currently exported. Reconstructed electron vectors are not substitutes in radiative events. Recovering vertex epsilon for old samples may be possible under validated hydrogen/collinear assumptions; otherwise regenerate the needed sample.

Keep Ntried, contributing/generated counts, genvol, luminosity, processed fraction, radiation settings and decay convention as sample metadata. Never substitute the number surviving final cuts for the generation denominator.

### D. Iterate through forward predictions

For unchanged generator support and radiative/detector prescription:

```text
full_weight_new(e) = full_weight_old(e)
                   * sigcm_new(vertex_e)/sigcm_old(vertex_e)

YMC_i(theta) = BR * sum_e [full_weight_old(e)/sigcm_old(e)]
                         * sigma_theta(vertex_e) * I_i(reconstructed_e).
```

BR is included only if not already accounted for. This also shows how conventional iteration and a corrected demodelled basis fit can share the same response. Binning and all cuts act on reconstruction; the cross section acts at the vertex.

The inspected radiative sampling factors do not themselves call the hadronic model, so changing only the model permits event reweighting within that approximation. Recalculate model-dependent radiative/acceptance integrals from the new weights. No new Geant4 transport is necessary for weight-only changes; changes to support, radiator kinematics, geometry, thresholds or reconstruction require the corresponding regeneration/reconstruction.

Fit all relevant reconstructed slices jointly; inspect Q2,W,t,phi and mass distributions. Do not replace the truth model with jagged reconstructed data/MC bin ratios. Re-extract at fixed reference points and iterate until both predictions and reported cross sections stabilize. Repeat with plausible alternative starting shapes; treat the spread as model dependence rather than tuning it away.

## 5. Minimal validation sequence

Choose one well-populated kinematic setting and a representative run/mask mixture first. No full-production campaign is needed to establish the method.

| Order | Controlled check | Required evidence |
|---|---|---|
| 1 | Identity and known-response injection | Same evaluator in fit/producer; recover zero and nonzero energy/position response; reproduce final selections |
| 2 | Threshold, pairing, mask, exclusivity closure | Cut flows and selected yields agree for identical reconstruction; loose precursor sample includes inward migration |
| 3 | Model-population variation at fixed detector response | Recovered response stable when Q2,t,phi/channel populations change |
| 4 | Independent cross-section pseudo-data | Recover constant, sloped and azimuthally modulated vertex cross sections; show migration matrix/feed-in and covariance pulls |
| 5 | Correlated background and radiation studies | Signal estimator closes after final cuts; mass-tail/cut dependence understood under validated radiation |
| 6 | Real-data model iteration | Differential yield ratios and extracted values stable across iterations and plausible starting models |

Set numerical tolerances against the intended uncertainty budget before these checks. No tolerances or bias magnitudes are inferred from visual mass agreement. Existing `smearing_closure_summary.csv` compares optimizer/full-sample objective evaluations (S:6887-6915); it is useful numerical consistency checking, not the independent physics closure in this table.

## 6. Reproduction and work performed

Performed: targeted source tracing, two staged exclusive SIMC histogram checks, primary-reference verification, and four analytic mechanism checks. No ROOT or Geant4 jobs were run. The analytic checks are saved in `/tmp/pi0_physics_audit_checks_20260911.py`; working notes are `/tmp/pi0_physics_audit_20260911.md` and `/tmp/pi0_{simc,geant,xsec}_audit_20260911.md`.

```bash
cd /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main
rg -n 'USE_SIM_MODEL_XSEC_DEMODELING|SIM_MODEL_XSEC_BRANCH|ENABLE_POSITION_SMEARING' src/simulation_smearing/*.C
rg -n 'passes_dead_block_mask|pass_nps_cluster|out_is_exclusive_sm' src/simulation_smearing/simc_pi0_analysis.C
rg -n 'base_w|sim_base_cos|row_total|row_plus|row_minus|tgt_contam' src/xsec_extract/excl_xsec_pi0_analysis_no_simc_model.C
python3 /tmp/pi0_physics_audit_checks_20260911.py
```

Analytic outputs: center-map example `1.00 -> 1.05`; a 0.61 GeV cluster with threshold 0.60 GeV and 0.05 GeV Gaussian width has pass probability `0.579260`, whereas frozen selection retains it; the migration example returns `(1.4,2.6)` from truth `(1,3)`; rescaling weighted yields by 1000 multiplies the Poisson deviance by 1000. These demonstrate mechanisms, not measured NPS effects.
