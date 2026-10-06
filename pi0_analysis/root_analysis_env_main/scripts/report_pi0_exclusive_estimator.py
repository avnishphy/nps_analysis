#!/usr/bin/env python3
"""Generate the blocked exclusive-estimator report from measured diagnostics."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot
from audit_pi0_exclusive_estimator import DATA, NOMINAL, REPO, read_csv, save_json, write_csv

FLOW = r'''# Actual event estimator: KinC_x36_4

```mermaid
flowchart TD
 A[HMS selection, accepted NPS clusters, best photon pair] --> B[Per-run all-candidate mass histogram]
 B --> C[Event timing-region subtraction]
 C --> D[Per-run Fermi combinatorial sideband fit]
 D --> E[Signed residual / all-candidate mass-bin count]
 E --> F[pi0 signal estimate: exclusive + SIDIS + delta]
 F --> G[Combined data-derived mpi0 versus missing-mass ellipse]
 G --> H[Exclusive-enriched sample; detector tprime x phi rows]
 H --> I[CURRENT: normalize and fit as exclusive]
 H --> J[REQUIRED: conditional background control and non-exclusive subtraction]
 J --> K[Exclusive detector yield - presently unavailable]
```

## Exact current definitions

For run r and 2 MeV mass bin m in [0,0.4) GeV:

```
c_e = I_C - I_D/6 - (I_H + I_V)/12 + (I_F1 + I_F2)/18
H_rm = sum_(e in r,m) c_e
Bhat_rm = A_r / (1 + exp((mass_bin_center - mt_r)/width_r))
S_rm = H_rm - Bhat_rm
p_e = pi0_weight_e = S_rm / N_all,rm   (0 when N_all,rm = 0)
```

The fit is per run, inclusive over tprime/phi and missing mass, AFTER upstream
HMS/NPS acceptance, photon-pair selection and timing subtraction, BEFORE any
exclusivity selection. It fits the sidebands [0.01,0.11] and [0.15,0.40] GeV.
The subsequent Gaussian is a mass-summary diagnostic, not the definition of
the weights. Weights are signed yield redistribution, not probabilities,
purities, sWeights, or exclusive/SIDIS classification. Empty denominator bins
receive no event weight even if a fit extrapolates residual content there.

Evidence: `src/analysis/nps_analysis_main.C:2988,3181,3216,3244,3344` and
`src/analysis/nps_comb_bg_pepsi.h:270,535,586`. Event cuts are at
`nps_analysis_main.C:2840-2988`: HMS electron/acceptance cuts, cluster packing,
dead-block mask, cluster energy/position/time cuts and eligible best pair.
No final cross-section diamond, tprime or phi subdivision enters this mass fit.
Timing coefficients come from `nps_time_bg.h` and are independently implemented
in `scripts/bootstrap_pi0_data.py:59`.

With exclusive, SIDIS, delta/non-exclusive real pi0, combinatorial and other
components, H estimates E + S + D + C + O after accidental subtraction. Hence
the intended mass residual is E + S + D + residual(C + O). Real SIDIS pi0 has
the same invariant-mass identity as exclusive pi0. Wrong photon combinations
in a SIDIS event are a separate reconstruction background; generator channel
labels alone do not certify correct photon pairing.

## Geometry and detector rows

`src/analysis/combine_analysis_branches.py:952` fits the common ellipse from
the combined target/kinematic sample using POSITIVE `pi0_weight*scale` within
0.11 <= mpi0 < 0.15 and 0.6 <= Mmiss < 1.5 GeV. It uses a smoothed 60 x 80
histogram, a peak-connected core, an automatic core-fraction scan, weighted
covariance, growth, and ridge alignment. Signed weights are still used in
the final yield. The geometry is data-derived, not a fixed proton-mass window.
Per-run ellipse/decorrelation flags also exist; the model uses the COMBINED
ellipse flag. The combiner applies that ellipse to every finite candidate in
the geometry domain, including negative-weight events (`:1269`).

For x=(mpi0,Mmiss), selection is `(x-mu)^T Cov^-1 (x-mu) <= d2cut`, within the
domain above. Actual stored parameters are:

```
mu = (0.13456952206348707, 0.94384662572514622) GeV
Cov = [[8.9892855956150908e-06, -0.00016704095020567288],
       [-0.00016704095020567288, 0.0051211477719482808]] GeV^2
d2cut = 10.154400260103738
```

`src/xsec_extract/xsec_analysis.h:241` uses the stored flag for data and
reapplies the stored data geometry to SIMC through `xsec_mass_cut.h:13`.
The Q2/xB diamond and rectangular ranges then select rows: Q2 [3.3,4.7],
xB [0.29,0.44], signed tprime=t-tmin edges [-0.75,-0.55,-0.40,-0.25,-0.13,0],
12 equal wrapped phi bins, row=12*tprime_bin+phi_bin. Kinematic float32
conversions are reproduced by the audit. The positive variable tau=-tprime
belongs to the physics model, not a different detector binning.

Let A_ei denote ellipse membership AND detector-row cuts. The current row is

```
scale_r = PS_r / ((Q_r,uC/1000) * LT_r * efficiency_r)
a_r = float32(scale_r) * float32(Q_r,uC) / Q_total,uC / 0.584
y_current,i = sum_e a_r * A_ei * p_e
v_conditional,i = sum_e (a_r * A_ei * p_e)^2
```

Exposure is the sum of 56 accepted serialized run charges:
880666.5126953125 uC = 880.6665126953125 mC. `scale` incorporates prescale,
livetime and efficiencies. 0.584 is the existing external hydrogen-yield
factor; it is not a SIDIS estimate. The result is yield/mC, not a cross section
or a bin-width-normalized density. The ordinary sum of squares freezes fitted
weights and geometry; it is not the full data variance.
See `combine_analysis_branches.py:635`, `xsec_accumulation.h:159`, and
`xsec_fit.h:13`. No such normalization constants were changed here.

## Required exclusive definition; not yet an executable corrected estimator

On the same normalization and detector selection, a valid subtraction form is

```
y_excl,i = sum_r a_r * [sum_(e in r) A_ei*c_e
                      - Bhat_comb,r,i(A) - Bhat_other,r,i(A)]
           - bhat_SIDIS,i(A) - bhat_delta,i(A)
```

The selected combinatorial estimate must condition on selection/kinematics;
it cannot simply transport the inclusive mass-bin mean. SIDIS and delta terms
here are already normalized in yield/mC. Equivalent additive forward fitting
would retain those terms in `f_i = f_excl,i(theta) + b_SIDIS,i + b_delta,i`.
There is no division by ellipse efficiency at the data stage: selected
exclusive acceptance is already handled in the folded exclusive response.

Current expectation need not equal an exclusive yield:
`E[y_current] = selected exclusive + selected SIDIS + selected delta
               + residual combinatorial/other + weight-transport bias`.
Mass subtraction cannot determine the individual real-pi0 components.

The existing fixed-window selected-sample fitter in `bootstrap_pi0_data.py`
does not validate this ellipse estimator. In particular, the ellipse's
0.11--0.15 GeV mass domain removes both existing combinatorial sidebands.
Refitting the unchanged sideband model only inside the ellipse cannot identify
its continuum. An independently constrained conditional continuum or validated
joint treatment is needed. No such ellipse treatment was found.

## Available non-exclusive treatment and uncertainty ownership

`scripts/generate_simc_infiles.py:82` creates exclusive/SIDIS/delta channels.
`src/simulation_smearing/simc_pi0_analysis.C:1050` reads all three and stores
their labels and `full_weight = Weight * normfac / Ngen` (`:2030`). The existing
post-analysis plotter compares these channels (`scripts/plot_publication_quality.py:1239`).
These are simulation/plotting capabilities, not a validated data subtraction.
The active model loader/selector rejects generated non-exclusive events
(`xsec_analysis.h:247`). There is no SIDIS fit, subtraction, sideband-derived
normalization, additive template or normalization nuisance in its objective.
Absence of treatment is not evidence that the contamination is negligible.

Use the existing channel templates as candidates for a future small additive
background treatment only after validating their conditional shapes and
normalization against data and resolving combinatorial transport. Do not
promote their generator normalization to a measured background. Data-derived
background amplitudes and geometry belong in the physical-event DATA ensemble;
finite channel-template statistics belong in independently resampled exclusive,
SIDIS and delta MC ensembles. Shared generator provenance must be checked before
asserting independent samples. Geometry from a data replica must also reselect
the response and all background templates. The 56 `signal_events_run*.root`
caches retain event IDs; the combined tree drops event_id.
'''


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);a=p.parse_args()
    out=a.input
    if (out/'REPORT.md').exists():p.error('REPORT.md already exists; use a fresh reproduction directory')
    summ=read_csv(out/'audit/sidis_summary_by_tprime.csv')
    closure=read_csv(out/'toys/estimator_closure.csv')
    transport=json.loads((out/'audit/weight_transport.json').read_text())
    # Independently verify the CURRENT per-run histogram ratio for every event.
    checks=[]
    with uproot.open(DATA) as f:events=f['physics'].arrays(['run_number','mpi0_all','pi0_weight'],library='np')
    for run in np.unique(events['run_number']):
        with uproot.open(DATA.parent/f'diagnostics_run{run}.root') as f:
            den=f[f'h_mpi0_all_run{run}'].values();sig=f['h_pi0_final'].values()
            ratio=np.divide(sig,den,out=np.zeros(200),where=den>0)
            np.testing.assert_allclose(ratio,f[f'h_pi0_weight_run{run}'].values(),rtol=0,atol=1e-14)
            take=events['run_number']==run;m=events['mpi0_all'][take]
            mb=np.searchsorted(np.linspace(0,.4,201),m,side='right')-1
            valid=(mb>=0)&(mb<200);expected=np.zeros(len(m));expected[valid]=ratio[mb[valid]]
            err=float(np.max(abs(expected-events['pi0_weight'][take])))
            assert err<1e-12
            checks.append(dict(run=int(run),max_weight_error=err,events=int(take.sum())))
    write_csv(out/'per_run_weight_verification.csv',checks)
    report='''# Exclusive-pi0 central-estimator decision - KinC_x36_4

**EXCLUSIVE-YIELD ESTIMATOR NOT VALIDATED.** The current mass-bin residual is a
true-pi0 estimator candidate, not an exclusive yield. It redistributes timing
and combinatorial subtraction before a correlated, data-fitted ellipse. The
model then fits the selected mixture with an exclusive-only response. Existing
SIDIS/delta MC is available, but no validated selected combinatorial treatment
or constrained residual real-pi0 subtraction was found. Production covariance,
M0 correction, nested models and final U/LT/TT errors were not run.

This audit covers model extraction only. No no-model comparison, charged-pion
diagnostic, hadronic-model redesign, git commit, or production-input change was
performed. Read [event_estimator_flow.md](event_estimator_flow.md) for source
locations, all event cuts, normalization and component equations.

## What the current weight means

Per run and 2 MeV mass bin: pi0_weight = (timing-subtracted candidates minus
fitted combinatorial continuum) / all-candidate count. It is signed yield
redistribution, not a posterior purity or sWeight. It is computed before the
ellipse; true SIDIS and delta pi0 remain part of the mass signal. The fit is
inclusive across missing mass/tprime/phi within each run. A later Gaussian
summarizes the mass peak; it does not calculate the weights.

Independent checks reproduced every combined-event weight against the original
per-run residual/all histogram for all 56 runs. Data row yields reproduce the
old nominal export; simulation exclusive row membership also matches the
model event cache. This establishes the actual estimator, not its validity.

## Composition: predictions available; data decomposition unavailable

The smeared file has 238727 exclusive, 104668 SIDIS, and 201128 delta events.
Templates use their stored physical full_weight normalization, independently
of the exclusive model fit. The table below is explicitly NOT a measurement
of the SIDIS fraction in data. Data-derived exclusive/SIDIS/delta components
and combinatorial residuals are NaN in the machine-readable tables.

| Signed reco tprime [GeV2] | Current selected /mC | MC SIDIS /mC | SIDIS / all MC | (SIDIS+delta) / all MC | MC SIDIS / current data |
|---|---:|---:|---:|---:|---:|
'''
    for r in summ:
        report+=f"| [{r['tprime_lo']}, {r['tprime_hi']}] | {float(r['ellipse_selected_pi0']):.6f} | {float(r['template_SIDIS']):.6f} | {100*float(r['template_SIDIS_fraction']):.3f}% | {100*float(r['template_nonexclusive_fraction']):.3f}% | {100*float(r['template_SIDIS_over_observed']):.3f}% |\n"
    report+='''
Raw exclusive simulation exceeds the current selected data by factors of about
1.5--2.1, so template-mixture fractions cannot be transported to the fitted
data composition without validating relative normalization. No arbitrary
percentage cut declares SIDIS negligible. Its U/LT/TT impact is unknown until
the selected true-pi0 estimator and background normalization are established.
This is an unresolved release blocker, not a finding of dominant SIDIS.

The search inventory is `sidis_search_evidence.txt` in the original audit
directory. Relevant capabilities are SIMC infile generation, three-channel
smearing/production, and post-analysis plots; none implements an active
exclusive-yield SIDIS subtraction. Delta is kept separate from SIDIS throughout.

## Correct interpretation of the 10.43% diagnostic

'''
    report+=f"For {transport['zero_background_runs']} nominal-zero-continuum runs, the SAME ellipse-selected events give {transport['mass_weight']:.12f}/mC with mass weights and {transport['timing']:.12f}/mC with direct timing coefficients: difference {transport['difference']:.12f}/mC, {transport['percent_relative_mass_weight']:.5f}% relative to mass weights.\n"
    report+='''
This is a conditional **timing-weight redistribution** mismatch. The direct
timing sum has not subtracted selected combinatorial background and is not
known exclusive truth. A fitted zero continuum is not proof of absent physical
combinatorial background. The comparison does not change ellipse membership,
does not label either real-pi0 source, and does not measure a composition shift.
Correlations of timing coefficients with selection can explain a difference
without any SIDIS. Source composition could affect such correlations but
cannot be inferred or apportioned from this diagnostic. Calling the entire
10.43% a measured combinatorial bias or a SIDIS fraction would be unjustified.

The exact zero-combinatorial-branch toy with no SIDIS and independent selection
closes for both estimators (biases about +0.13 events with MC errors 0.17--0.19).
With correlated accidental acceptance, mass redistribution recovers about
333.23 of 799.55 selected injected pions; direct timing gives 799.45. These are
10000 trials per condition, seed 20261005. This demonstrates the statistical
failure mechanism, not its numerical magnitude in data.

## Current-chain component toys

Seed 20261004; 12 requested pseudoexperiments per case; 56 synthetic runs per
pseudoexperiment. Each reruns the production C++ mass fitter per run, forms
the same residual/all weights, refits the actual combined ellipse geometry,
and accumulates the same detector rows and normalization. Its SIDIS treatment
is the current one: NONE. Truth is the injected exclusive yield inside each
replica's own ellipse, not the pre-selection yield.

True-pi0 shapes and relative rates come from existing exclusive/SIDIS simulation.
They are channel-origin proxies, not validated photon-ancestry templates.
Combinatorial mass follows the existing Fermi form, at a declared 10% of the
exclusive count with uniform missing mass; that rate is a stress assumption,
not a data measurement. Case 5 uses SIDIS x1.2 and missing mass +30 MeV. These
component toys have prompt timing; the separate transport test covers timing
accidentals. They are diagnostic tests of the current chain, NOT full validation
of a corrected estimator or calibrated realistic data ensembles.

| Toy | Components | Completed/requested | Injected exclusive /mC | Current output /mC | Bias /mC | Selected SIDIS /mC |
|---:|---|---:|---:|---:|---:|---:|
'''
    labels=['exclusive','exclusive + combinatorial','exclusive + SIDIS','exclusive + SIDIS + combinatorial','varied SIDIS + exclusive + combinatorial']
    for case in range(1,6):
        r=[v for v in closure if int(v['case'])==case and v['scope']=='tprime']
        nums=[sum(float(v[k]) for v in r) for k in ['injected_exclusive','recovered_current_chain','bias','injected_selected_SIDIS']]
        report+=f"| {case} | {labels[case-1]} | {r[0]['successful']}/{r[0]['requested']} | "+' | '.join(f'{v:.6f}' for v in nums)+' |\n'
    report+='''
Toys 1/2 each had ellipse-fit failures; toy 5 had mass-fit failures. Conditional
means over completed replicas are shown only diagnostically; failures were
saved, not discarded to claim closure. `estimator_closure.csv` supplies biases,
empirical residual spreads and MC errors for every row and tprime slice. The
full exclusive estimator does not close: no corrected estimator exists, and
the successful exclusive+SIDIS case retains approximately its surviving SIDIS
yield. The geometry itself also changes with mixture. These tests are sufficient
to reject release; they are not a coverage or uncertainty study.

## Rows and M0 decision

Rows 0,5,6,11 have respectively 289,318,317,299 events before the ellipse,
zero observed selected events, and 22,3,9,32 selected exclusive MC events.
Their old conditional variances are zero. Corrected exclusive yield and full
bootstrap variance remain unavailable, not zero. No data-dependent exclusion
is endorsed. Neither reinstatement with an invented variance nor retention of
the old exclusion is a validated solution. A supported zero-count statistical
treatment must accompany the repaired estimator.

Old diagnostic M0 Q = 54.182524936522192. Old detector/nuisance/physics ranks
were 16/12/4; old predictions obeyed the checked positivity constraints and
six starts converged through 15 variance iterations. These are inherited
diagnostic results, not checks of a corrected fit. The corrected fit, ranks,
positivity, convergence and U/LT/TT shifts are unavailable. The old values and
NaN corrected values are in `central_M0_before_after.csv`; no no-model result
is used. No LT/TT slopes were added.

## Concrete decision and next prerequisite

The available estimator cannot yet distinguish/control selected SIDIS/delta
and conditional combinatorial background sufficiently to report an unbiased
exclusive yield. Stop here. The next work must establish a selected true-pi0
estimator retaining usable continuum constraints and validate the existing
non-exclusive channel shapes/normalization, then repeat complete closure.
This is the same central-estimator task, not another model diagnostic campaign.

The data bootstrap cannot now be called authoritative. Once repaired it must
resample physical run/event IDs and rerun mass/background fits, weights,
geometry, any data-derived channel normalization, yields and geometry-dependent
response/template selection. Exclusive/SIDIS/delta finite MC must be separately
resampled as appropriate. No production command is valid at present. The new
`central` command intentionally exits 2 with this reason and writes no result.

## Artifacts and reproduction

- [Event flow and equations](event_estimator_flow.md)
- [Per-row composition](audit/exclusive_composition_by_row.csv)
- [tprime composition summary](audit/sidis_summary_by_tprime.csv)
- [Current-chain closure](toys/estimator_closure.csv)
- [Former zero-variance rows](audit/corrected_row_variances.csv)
- [M0 before/after availability](audit/central_M0_before_after.csv)
- [Flow figure](audit/event_estimator_flow.pdf)
- [Mass before ellipse](audit/mpi0_before_exclusivity.pdf)
- [Missing mass and separate templates](audit/missing_mass_templates.pdf)
- [Template SIDIS fractions](audit/sidis_fraction_vs_tprime.pdf)
- [Closure bias](toys/exclusive_yield_closure_bias.pdf)
- [Supported empty rows; bootstrap unavailable](audit/formerly_zero_variance_rows.pdf)

Figures include the requested available diagnostics. No figure claims a
corrected bootstrap or measured SIDIS fraction. NaN in CSV means unavailable.

From the repository root, the complete tested reproduction driver is:

```bash
cd /work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main
bash validation/exclusive_estimator_20261004/reproduce.sh validation/exclusive_estimator_rerun_20261004
```

The driver refuses an existing destination, builds the unchanged C++ fitter
under the Hall C environment, runs audit and toys, and generates this report.
It contains exact independent audit/toy/report commands, not placeholders.
For the corrected central-estimator availability check (expected exit 2):

```bash
python3 scripts/audit_pi0_exclusive_estimator.py central --output validation/exclusive_estimator_20261004/central
```

No source/config file that existed before this task was modified. New source:
`scripts/audit_pi0_exclusive_estimator.py` and
`scripts/report_pi0_exclusive_estimator.py`. The reproduction driver and all
other new files are under `validation/exclusive_estimator_20261004/` (plus the
independent reproduction directory). The initial git HEAD, status and diff-stat
and reference checksums were recorded before additions. Final checksum/status
evidence is retained alongside this report in the original audit directory.
No commit was made.

**EXCLUSIVE-YIELD ESTIMATOR NOT VALIDATED**
'''
    (out/'event_estimator_flow.md').write_text(FLOW)
    (out/'REPORT.md.tmp').write_text(report);(out/'REPORT.md.tmp').replace(out/'REPORT.md')
    save_json(out/'release_decision.json',dict(validated=False,production_allowed=False,
        reason='Conditional true-pi0 estimator and residual SIDIS/delta control not validated',
        production_data_replicas=0,production_MC_replicas=0,corrected_M0_refit=False,
        new_source_hashes={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [REPO/'scripts/audit_pi0_exclusive_estimator.py',Path(__file__).resolve()]}))
    print(out/'REPORT.md')


if __name__=='__main__':main()
