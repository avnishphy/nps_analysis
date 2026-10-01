# Run 4398 livetime investigation: current entry point

Updated 2026-09-10. The six-segment check is complete and is now part of the
full 43-run investigation. Read the [complete report](livetime_complete_20260910/REPORT.md),
[presentation](livetime_complete_20260910/KinC_x60_4b_LH2_livetime_presentation.pdf)
and [research log](LIVETIME_RESEARCH_LOG.md). Earlier detailed material is
preserved under `livetime_4398_check_20260909/` and in the new evidence package.

## Exact comparison, with input populations distinguished

| Quantity | E | D | Ratio |
|---|---:|---:|---:|
| Saved original one-segment NewGen | 19341 | 19310 | 1.0016053858 |
| Same segment, unsupported tail removed | 19292 | 19310 | 0.9990678405 |
| Original NewGen on all six updated files | 112586 | 112746 | 0.9985808809 |
| Same six files, per-file coverage only / old current | 112263 | 112746 | 0.9957160343 |
| Ending-current aligned, final stitched broad raw diagnostic | 112658 | 112746 | 0.9992194845 |
| Matched tight +/-2 ns diagnostic | 112650 | 112746 | 0.9991485286 |

The core pE/D estimator was already the user's method. The net same-file change
is +72 numerator tags, with unchanged denominator. It is not simply dropping
tails. The new [step comparison](livetime_complete_20260910/figures/run4398_refinement_steps.pdf)
separates the contributions. Raw-subtracted proposed ratio: 0.9991656330;
the legacy tight-subtracted result differs. CLT_physics is never the baseline.

Full4398 has 3,222,394 contiguous events and 2233 real scaler snapshots.
The earlier fixed diagnostic band 33.3625--45.1375 uA gives N=A=2893284,
S=2895694, D=112746, Eraw=112658, Etight=112650, 2821.898170 s and
109930.299545 uC. Its alternative boundary inclusion adds two tight tags;
this local boundary sensitivity is not an all-run systematic confidence limit.

Segment0's terminal row repeats the counters of boundary541321 while changing
the event boundary to542585. The unsupported1264 events contain49 tags.
The initial snapshot already includes878 sent pulses before common event
coverage. Adjacent segments can recover exposure; do not discard each file
tail mechanically. The retained loss burst near820.125614 s has S2082,
N=A1799, D80, E69 and is distinct from the 4303/4305 counter defects.

The TSHelH.evcount/TSH.evcount mismatch is a separate diagnostic issue: record
indices of the two streams are not interchangeable. It did not invalidate
the saved main NewGen D19310. Do not promote the mismatched gated diagnostic
18970/6550 to a livetime.

## Physical interpretation and provenance

4398 uses TI4 HMS EL-REAL singles under the user-confirmed August2023 routing.
An NPS hardware cluster condition is not required by that selected input.
NPS readout/reconstruction remain separate considerations. The detailed HMS
EDTM path and sent-counter tap still need evidence; no total LT is certified.
PRE widths/rates do not identify their upstream trigger.

The database CPU_LT0.9992 and master metadata Computer_live_time0.992 have
different, incompletely resolved provenance. Neither is a truth reference.
Replay job/tar evidence is scoped to the recorded4398 job; it is not proof of
every pass2 executable. Exact source snapshots are preserved.

The historical top-level memory is saved in
`livetime_complete_20260910/previous_docs/LIVETIME_4398_INVESTIGATION.md`.
The earlier staging request completed; no staging or production change is
needed for this documentation. The one later raw-file allowance was used for
4305, not4398, and no additional raw staging is authorized.
