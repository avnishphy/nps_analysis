# Pi0 workflow evidence and dependency map

Date: 2026-10-01

Workspace: `root_analysis_env_final`

Scientific status: source map only; no publication extraction is selected or
certified.

## Evidence labels and boundaries

- **Source-confirmed** means the cited committed code implements the stated
  behavior. It does not establish the size of an effect in production data.
- **Synthetic-check-only** means a deterministic calculation or toy study was
  run without experimental data. Exact commands are in
  `validation/formula_checks_20261001.md`.
- **Production-confirmed** would require a versioned real-data/SIMC run with a
  complete manifest. No item in this map has that label.
- **Unresolved** marks missing data, provenance, or validation that prevents a
  scientific conclusion.

The read-only audit used as an input is
`/work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/docs/pi0_uncertainty_audit_2026-10-01.md`
(SHA-256
`9f0df3ea13872358dc7766f7e0d75d29079e0418f846ecb841ab2ee536f64ddb`).
Its findings were rechecked against the copied `FINAL` source; the audit itself
was not edited.

## Executable call graph and active defaults

```text
scripts/run_pipeline.sh
  |-- optional SIMC config/infile generation and simulation wrapper
  `-- src/analysis/run_parallel_nps_analysis_main.sh
        |-- one src/analysis/nps_analysis_main.C invocation per run
        |-- src/analysis/combine_analysis_branches.py       [default: yes]
        |-- src/simulation_smearing/run_smearing_pipeline.sh [default: no]
        `-- src/xsec_extract/run_xsec_pipeline.sh            [default: no]
```

The top launcher defaults runtime products to `FINAL/output`. The analysis
driver defaults to combining runs, but smearing and cross-section extraction
remain opt-in. The extraction wrapper defaults to the no-SIMC-model response
method, Gaussian objective, finite-MC variance iteration, and no positivity
constraint. The alternative C++ SIMC-model path is Gaussian only. The separate
Python forward-event bootstrap declares `publication_ready=False`; the joint
multi-setting solver is Gaussian and has a different parameterization. These
are distinct estimators, not interchangeable output backends.

The active acceptance configuration is `config/acceptance_cuts.conf`, SHA-256
`1f2ebf451c8ffd22f3d3f06140a839cf83cb729c9d88d8da0f9ce317136bca05`.
Its timing-sideband switch is `auto`, so the acquisition-mode defaults described
below are active unless a launcher/environment override is recorded.

## Stage map

| Stage | Observations and sampling unit | Selection | Estimator / weights / units | Estimated versus fixed inputs | Outputs and downstream use | Uncertainty treatment and missing information |
|---|---|---|---|---|---|---|
| Raw production input | ROOT tree entries, split into run segments | Branch availability and macro-level event cuts | Unweighted event records in detector conventions | Run configuration fixed; detector response embodied in input | Per-run analysis input | Segment completeness is not carried into the final combined event tree. |
| Efficiency/livetime | Trigger/event counts by run or segment | Tracking: electron/hodoscope sample; hodoscope: reconstructed-track acceptance; EDTM counters | Tracking and hodoscope are binomial plug-in ratios; EDTM livetime is `P*N/D` | Central corrections estimated; behavior explicitly frozen | Efficiency CSV consumed by combiner | Tracking/hodoscope numerator and denominator are not exported by the main CSV. EDTM propagation treats numerator and denominator independently. Partial segments may contribute charge while a missing branch suppresses a metric. |
| Per-run event analysis | One selected event; candidate photon pair | Detector/kinematic cuts; exactly two clusters bypass pair-time-difference selection; multiplicity >=3 chooses the allowed pair closest to nominal pi0 mass | Raw histograms, timing-region counts, then a mass-bin signal fraction | Timing mode and cut values fixed by configuration; timing and combinatorial backgrounds estimated from the same run | `diagnostics_run<run>.root`, summary row, selected-event tree | Output event tree retains `event_id` but not both photon times or a disjoint timing-category label. The learned background parameters and their covariance are not propagated with each event. |
| Timing accidental estimate | Counts in prompt, diagonal, horizontal, vertical, and two full-box timing regions | Region boundaries plus the pair selector | `A = D + (H+V)/2 - (F1+F2)/2`, after each aggregate is area-scaled to the 4 ns2 prompt box | Region counts estimated; geometric transfer factors fixed | Accidentals-subtracted mass template | Raw-count Poisson terms are propagated once in the scalar estimate, then related errors enter the template again through `Sumw2` plus a normalization term. Full-box geometric areas do not equal the accepted areas after multiplicity-dependent pair selection. |
| Combinatorial mass background | Prompt/accidental-subtracted mass bins | Configured mass fit range | Fermi/logistic fit and subtracted mass spectrum | Fit parameters estimated; invalid fit can continue with current parameters and zero covariance | Final mass histogram used to form purity weights | Fit covariance is diagnostic only. Negative subtracted bin content is clipped to zero. |
| Event purity assignment | Selected event in mass bin `j` | Same selected event pool | `p_j = N_final,j / N_all,j`; event retained downstream only for `p_j > 0` | `p_j` estimated from shared data but then treated as fixed | `pi0_weight` in per-run and combined event trees | Conditional `sum(w^2)` omits shared uncertainty and bin correlations. A mass marginal does not demonstrate differential closure in phi, t, Q2, xB, helicity, or exclusivity. This is not an sWeight construction. |
| Run combination | One selected event row | Skip missing ROOT, missing efficiency, or empty trees | `scale_r=P_r/[Q_r(mC)L_r e_r]`; copied to every event with an approximate error | Run charge, livetime, tracking and hodoscope corrections estimated; prescale fixed | `combined_branches_<target>.root`, combined mass-cut metadata | `event_id` is omitted. No run ledger is emitted, so valid zero-candidate runs cannot reach extraction. Processing status and auxiliary efficiency counts are not required. A combined 2D mass cut is re-estimated with `pi0_weight*scale`. |
| Smearing calibration | Data and SIMC events; each generator event is repeated `Nsmear` times | Exclusive and geometry cuts; event contributes to every detector section touched by either photon | Data weight `pi0_weight*scale*charge_fraction`; simulation normalized to data integral per observable; default objective is ordinary Baker-Cousins deviance | Section smearing parameters estimated; data weights and normalization treated as fixed | Section map CSV/ROOT and smeared SIMC products | Ordinary count deviance is applied to weighted, arbitrarily normalized bins. Overlapping section objectives share events. Repeated-smear `sumw2` is copy-wise, not generator-event-wise. Full joint parameter covariance is not exported. RNG seeding is time-dependent; requested producer seed is recorded but not applied by the legacy producer. |
| No-SIMC-model response extraction | Combined weighted data rows and accepted SIMC events | Reconstructed selection for data; truth/reconstructed support for SIMC | Data `sumw` and conditional `sumw2`; response uses `full_weight/sigcm` times U/LT/TT basis | Cross sections fitted; purity, run corrections, target divisor, smearing map and response generally fixed or only partly modeled | Cross-section ROOT/CSV/plots | Gaussian mode can exclude rows with nonpositive data variance and adaptively drop groups; finite-MC is an iterative row-variance approximation. Scaled-Poisson mode treats positive event weights as fixed, borrows a scale for empty rows, fixes MC response, and uses boundary-dependent Minuit errors. |
| SIMC-model ratio extraction | Weighted data and simulated reconstructed yields | Similar reconstructed selection | Data/SIMC ratio times generator model; optional historical missing-mass area match | Model/reweight parameters may be fitted; normalization mode differs | Separate ROOT/CSV/plots | Gaussian objective only; model dependence and conditional upstream weights remain. |
| Python forward bootstrap | Exported weighted data events and SIMC response events | Rejects negative weights | Conditional scaled-Poisson fit; Poisson event bootstrap for accepted data and MC | Final exported event weights fixed inside resampling | Diagnostic forward-fit products | Explicitly not publication-ready. Cannot recover raw timing/background or auxiliary-efficiency fluctuations absent from exports. |
| Joint multi-setting fit | Per-setting prepared response/yield products | Configured common blocks and setting support | Gaussian joint fit with epsilon-dependent separation | Shared cross sections fitted; response and upstream weights conditionally fixed | Joint separation products | Optional finite-MC/positivity modes exist, target uncertainty is omitted, and a nominal setting epsilon enters the separation. |

## Selection and background dependencies

The current pair and timing logic is multiplicity dependent:

- With exactly two accepted clusters, the sole pair is accepted without the
  pair-time-difference gate.
- With three or more clusters, candidate pairs must satisfy the active
  `|t1-t2|` threshold and are ranked by closeness to nominal pi0 mass, with
  total energy as a tie-breaker.
- HCANA-like defaults use a 10 ns pair threshold and unshifted windows;
  waveform defaults use 13 ns and shifted windows.
- The timing histogram spans 140--160 ns, whereas shifted waveform sidebands
  extend to 139 and 161 ns. The histogram-based raw estimator therefore does
  not represent the entire configured shifted region.

The nominal aggregate areas are 4 ns2 for prompt, 24 ns2 independently for
each of diagonal/horizontal/vertical, and 36 ns2 for each full box. Within one
full box, the accepted area is 8 ns2 for the unshifted 10 ns gate and 4.5 ns2
for the shifted 13 ns gate; it is 36 ns2 for the two-cluster bypass. Thus one
fixed full-box transfer cannot simultaneously describe these strata. This is a
geometric statement, not a measured bias.

## Exposure and normalization derivation

Define, for run `r`:

- `Q_r` = accepted run charge in microcoulombs;
- `q_r = Q_r/1000` = the same charge in millicoulombs;
- `P_r` = prescale factor;
- `L_r` = livetime;
- `e_r = e_track,r * e_hodo,r` = correction efficiency;
- `p_i` = inferred pi0 mass-bin weight for event `i`;
- `Q_tot = sum_r Q_r` over runs represented in the combined event rows.

The combiner stores

```text
scale_r = P_r / (q_r L_r e_r)                         [1/mC]
```

and extraction multiplies by the dimensionless charge fraction `Q_r/Q_tot`:

```text
w_ri = p_i scale_r (Q_r/Q_tot)
     = p_i P_r / [(Q_tot/1000) L_r e_r]               [1/mC].
```

The individual run charge cancels algebraically; inserting another factor of
1000 would be a unit error. This cancellation also shows why the run set used
to form `Q_tot` is part of the estimator. At present that set is discovered by
scanning event rows, so a valid selected run with zero candidate rows is absent
from the exposure.

The dimensionless hydrogen target divisor is applied once after accumulation.
The response fit uses

```text
A_rbn = sum accepted SIMC (full_weight/sigcm) * basis_n(phi_vertex, epsilon)
y_r   = sum_b,n A_rbn sigma_bn.
```

`full_weight/sigcm` is an integrated yield-per-cross-section factor, not a
probability-normalized migration matrix. Consequently no additional Q2, xB,
t-prime, or phi bin-width divisor belongs in this equation. Fitted coefficients
retain SIMC `sigcm` units (microbarn/MeV2); display plots multiply by `1e9` to
show nb/GeV2. This unit chain is source-confirmed, while its external SIMC
normalization and data/MC closure are not production-confirmed here.

## Supported extraction paths and statistical meaning

| Path | Current default? | Data likelihood / variance | MC statistics | Positivity and intervals | Publication status |
|---|---:|---|---|---|---|
| No-SIMC-model Gaussian response | Yes when xsec is enabled | Weighted least squares using conditional data `sumw2` | Default iterative diagonal row-variance approximation | Optional constrained solve; covariance is local/conditional | Not established |
| No-SIMC-model scaled Poisson | No | Scaled-Poisson objective for positive fixed weights | Fixed response in objective | Positivity forced; Hesse and selected MINOS profiles | Not established |
| SIMC-model ratio/iterative | No | Gaussian | Method-specific fixed/iterative response | Method-specific | Not established |
| Python forward bootstrap | No | Conditional scaled Poisson | Poisson accepted-event bootstrap | Diagnostic profiles/bootstrap | Source explicitly says false |
| Joint multi-setting | No | Gaussian | Optional approximation | Optional positivity | Not established |

ROOT's histogram documentation defines `Sumw2` as accumulated squared weights;
that is the conditional variance convention used here, not propagation of an
estimated shared weight. ROOT's fitting guide distinguishes Poisson likelihood
for count histograms from weighted-histogram treatment. Baker--Cousins concerns
Poisson counts. Barlow--Beeston introduces nuisance treatment for finite
simulation counts. Those references motivate tests but do not, by themselves,
make any method applicable to these signed or inferred weights.

Primary/official references:

- ROOT `TH1` reference: <https://root.cern.ch/doc/master/classTH1.html>
- ROOT fitting guide: <https://root.cern.ch/manual/fitting/>
- ROOT Minuit2 guide: <https://root.cern.ch/root/htmldoc/guides/minuit2/Minuit2.html>
- Baker and Cousins (1984): <https://doi.org/10.1016/0167-5087(84)90016-4>
- Barlow and Beeston (1993): <https://doi.org/10.1016/0010-4655(93)90005-W>
- Pivk and Le Diberder, sPlot: <https://arxiv.org/abs/physics/0402083>
- Bohm and Zech, weighted Poisson histograms: <https://arxiv.org/abs/1309.1287>
- Particle Data Group statistics review: <https://pdg.lbl.gov/2025/reviews/rpp2025-rev-statistics.pdf>

## Revalidated finding-to-source map

| Finding | Evidence | Principal source locations | Consequence still to establish |
|---|---|---|---|
| BKG-001 clipping/nonpositive rejection | source-confirmed | `src/analysis/nps_time_bg.h:489-502`; `nps_analysis_main.C:3208-3218,3312` | Production shift by final bin |
| BKG-002 multiplicity-dependent timing acceptance | source-confirmed + synthetic geometry | `nps_analysis_main.C:1570-1573,2960-2968`; `nps_helper.h:917-953`; `nps_time_bg.h:109-235` | Rate and kinematic dependence in data |
| BKG-003 duplicated control-count uncertainty | source-confirmed | `nps_time_bg.h:224-235,270-290,489-502` | Exact covariance inflation/omission by bin |
| BKG-004 shifted-window histogram truncation | source-confirmed | `nps_analysis_main.C:2603-2607,2698-2707` | Entries affected in waveform data |
| BKG-005 invalid mass fit can continue | source-confirmed | `nps_comb_bg_pepsi.h:267-278,576` | Frequency and affected runs |
| EXP-001 zero-candidate exposure omission | source-confirmed mechanism; affected runs unresolved | `combine_analysis_branches.py:585-647`; `xsec_accumulation.h:408-430` | Whether any intended valid run is absent |
| EXP-002 partial processing status ignored | source-confirmed | `compute_efficiencies_stuff.cxx:818-822`; combiner required-column logic | Which production rows are partial |
| EFF-001 auxiliary counts unavailable downstream | source-confirmed | efficiency headers and CSV export schema | Valid nuisance model and correlations |
| SMEAR-001 ordinary deviance on weighted normalized bins | source-confirmed + synthetic unit check | `nps_sim_smearing_new.C:436,1314-1335,1434-1473,1939-1948` | Parameter and cross-section change |
| SMEAR-002 copy-wise finite-MC moments | source-confirmed + synthetic variance check | `nps_sim_smearing_new.C:1412-1417,1965,2175-2186` | Correct event-level covariance magnitude |
| SMEAR-003 overlapping section objectives/no joint covariance | source-confirmed | `nps_sim_smearing_new.C:4528-4566,7192` | Parameter cross-correlation and coverage |
| SMEAR-004 requested seed not applied | source-confirmed | `run_smearing_pipeline.sh:129`; `simc_pi0_analysis.C:2655-2659`; time seeds in smearing source | Reproducibility of a real rerun |
| XSEC-001 Gaussian row exclusion/recovery | source-confirmed | `xsec_fit.h:81-145,227-244` | Dataset-dependent retained support |
| XSEC-002 scaled-Poisson is conditional | source-confirmed | scaled-Poisson solver/output metadata | Coverage under inferred purity and finite MC |
| XSEC-003 bootstrap remains conditional | source-confirmed | `run_forward_xsec.py`; `forward_xsec_statistics.py` | Coverage once upstream auxiliaries vary |
| XSEC-004 joint-fit omissions | source-confirmed | `run_joint_xsec_fit.py`; `xsec_joint_solver.C` | Multi-setting separation coverage |

## Output provenance status

`FINAL/output` contains no production result established by this work. Existing
publication material in the snapshot predates the new frozen baseline unless a
separate manifest proves otherwise. No real run, smearing fit, extraction, or
end-to-end pipeline was executed during this evidence pass. A raw production
file for run 4398 was located in the external waveform directory, but it was
not opened or processed. Therefore all numerical production impacts remain
unresolved.

Any future publication candidate needs, at minimum, the Git commit, complete
configuration hashes, run/segment ledger including zero-candidate runs,
software environment, RNG seeds and actual seed application, input/output
checksums, fit status, covariance/nuisance products, validation results, and an
explicit approval record for every physics-sensitive change.

## Dependency order implied by the evidence

1. Decide whether inference starts from disjoint raw observation categories or
   from signed subtracted summaries with a complete covariance.
2. Preserve identifiers, run ledgers, timing categories, and auxiliary counts
   needed by the selected statistical model.
3. Specify timing and combinatorial-background nuisance models by acquisition
   mode and multiplicity.
4. Specify efficiency/livetime nuisance models without silently changing the
   frozen central-value behavior.
5. Rebuild smearing inference and finite-MC treatment with reproducible seeds
   and joint covariance.
6. Validate one bounded setting, then multi-setting separation and systematic
   combinations.
7. Only after coverage, closure, stress tests, and explicit approval may a path
   become a publication default.

The first approval package, `docs/proposals/ALG-001_raw_observation_forward_architecture.md`,
addresses step 1 only. No part of it is implemented.
