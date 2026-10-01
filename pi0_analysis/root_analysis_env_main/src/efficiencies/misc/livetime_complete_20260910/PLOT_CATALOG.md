# Plot catalog: scope and provenance

The nineteen current figures below are in the report. Eighteen have both
vector PDF and PNG; the applied NewGen reference is the unchanged supplied PNG.
Eleven figure pairs were made for this documentation edition from existing
results. Seven pairs were copied unchanged from the full-cache analysis.
No ROOT/event reread was used. The presentation adds an editable TikZ routing
schematic from the dated user-supplied logbook text; it is not a measured timing
or gate diagram.

Current numerical population: all 43 production runs / 183 segments, except
figures explicitly restricted to a run or cohort. No proposed ratio is truth.
Raw/tight labels and red counter-defect flags must survive any later relabeling.
For scatter plots without errors, points are count ratios, not uncertainty-free
measurements. Rate associations are descriptive, not correction fits.

## Current figure files

| Figure | What it establishes or illustrates | Source |
|---|---|---|
| [applied_NewGen.png](figures/applied_NewGen.png) | Original applied multipanel, preserved unchanged; saved coverage can differ from refreshed files | User path / full-cache snapshot |
| [refresh_livetime_comparison.pdf](figures/refresh_livetime_comparison.pdf) | Broad EDTM, tight timing, all-trigger and proposed raw-subtracted ratios; retains values above one | Full-cache `make_refresh_figures.py` |
| [refresh_normalization_sensitivity.pdf](figures/refresh_normalization_sensitivity.pdf) | Separates changed input files from same-file method sensitivity; fixed yield/charge illustration | Same |
| [refresh_correlated_difference.pdf](figures/refresh_correlated_difference.pdf) | Raw/tight differences with paired statistical scales; no truth baseline | Same / acceptance tests |
| [refresh_phase_timing.pdf](figures/refresh_phase_timing.pdf) | Relative/corrected timing versus independent TI-clock phase | Same / phase bridge |
| [refresh_rejected_tags.pdf](figures/refresh_rejected_tags.pdf) | Tight-rejected broad tags retain independent pulse association | Same / phase checks |
| [refresh_counter_exposure.pdf](figures/refresh_counter_exposure.pdf) | Recorded timestamp exposure versus scaler clock exposes 4303/4305 deficits | Same / event-clock checks |
| [raw4305_exposure_step.pdf](figures/raw4305_exposure_step.pdf) | Persistent raw L1/event deficit and duplicate unassigned pulse-like channels | Same / single raw4305 file |
| [cohort_ti6.pdf](figures/cohort_ti6.pdf) | Five TI6 coincidence runs; four p=1 and 4259 p=5 | New `make_comparisons.py` |
| [cohort_ti4_early.pdf](figures/cohort_ti4_early.pdf) | First 17 unprescaled HMS EL-REAL runs; counter flags retained | Same |
| [cohort_ti4_late.pdf](figures/cohort_ti4_late.pdf) | Remaining 16 unprescaled HMS EL-REAL runs | Same |
| [cohort_ti4_prescaled.pdf](figures/cohort_ti4_prescaled.pdf) | Five HMS EL-REAL runs with p=2 | Same |
| [method_decomposition_current.pdf](figures/method_decomposition_current.pdf) | Sequential coverage/current/joins changes and net change, in percentage points; order-dependent | Same / per-run JSON |
| [run4398_refinement_steps.pdf](figures/run4398_refinement_steps.pdf) | Same six files: E112586 ->112263 ->112658 ->112658, D112746 | Same / verified refinement JSON |
| [rate_comparison_by_trigger.pdf](figures/rate_comparison_by_trigger.pdf) | p=1 cohorts versus S/exposure; paired-difference error bars; defective exposure can affect both axes | Same / acceptance tests |
| [current_intervals_4303.pdf](figures/current_intervals_4303.pdf) | S-pN, N-A and D-pEraw on selected intervals; anomalous counter exposure | Same / current interval CSV |
| [current_intervals_4305.pdf](figures/current_intervals_4305.pdf) | Corresponding raw-verified counter defect; no repair adopted | Same |
| [current_intervals_4398.pdf](figures/current_intervals_4398.pdf) | Real loss bursts with N=A retained; distinct from N>A inconsistency | Same |
| [current_intervals_4551.pdf](figures/current_intervals_4551.pdf) | Prescaled interval history; S-pN can fluctuate from prescale sampling | Same |

The new cohort plots show the original estimator evaluated on the same files,
matched broad raw EDTM, matched tight timing and the proposed raw-subtracted
ratio. Their nominal count formulas are unchanged. The full-cache comparison
also supplies the all-trigger ratio. A common selected interval does not by
itself complete physics-yield/good-helicity/charge matching.

## Historical plots, explicitly not the current sample

`figures/historical_168_segments/` retains all seventeen figures from the old
deck: applied_NewGen; comparison_all; comparison_early; comparison_late;
comparison_unprescaled; coverage; decomposition; intervals_4303;
intervals_4305; intervals_4398; intervals_4551; normalization_effect;
prescales; sensitivity_all; timing_peaks; timing_sensitivity; uncertainties.
All except the unchanged applied PNG are PDFs. Their numerical sample was the
earlier 168-segment audit. Do not silently substitute these for current figures.

The complete old deck and generators are under `evidence/historical_beamer/`;
the original numerical archive is under `evidence/historical_168_run_audit/`.
Run4398-only figures (endpoint, timing, current and retained-loss diagnostics)
remain under `evidence/run4398/`. Other source-paper diagrams remain in the
original PDFs with their original conditions and attribution.

`artifact_inventory.json` enumerates every packaged figure, table, source file
and script with a checksum. It complements the narrative catalog rather than
assigning current applicability to every historical file.
