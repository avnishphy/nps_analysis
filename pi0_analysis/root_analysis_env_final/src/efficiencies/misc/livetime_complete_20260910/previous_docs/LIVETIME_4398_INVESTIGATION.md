# Run 4398 livetime investigation memory

## Latest: editable technical presentation and source extension

See [Beamer PDF](livetime_beamer_20260909/KinC_x60_4b_LH2_livetime_Beamer.pdf),
[LaTeX](livetime_beamer_20260909/livetime.tex) and
[research/definitions](livetime_beamer_20260909/RESEARCH.md).
The 56-slide deck retains the frozen 168-segment numerical results, adds exact
stage/weighting definitions and historical formulas, and explicitly cautions
that tight corrected-time EDTM cuts need an acceptance test. New NPS sources
point to beginning-of-run EVIO configuration records; PRE routing remains open.
User staging is in progress; later availability is archived separately, not
folded into the prior ratios. No jcache or production correction change.
Temporary exported columns may no longer exist; regenerate from the frozen
manifest if needed. Full new-cache reanalysis remains a separate next step.

## Current continuation: all KinC_x60_4b LH2 runs

The 2026-09-09 cache-only extension is complete: see
[`LIVETIME_KinC_x60_4b_LH2.md`](LIVETIME_KinC_x60_4b_LH2.md) and
[`livetime_KinC_x60_4b_lh2_20260909/REPORT.md`](livetime_KinC_x60_4b_lh2_20260909/REPORT.md).
46 cached runs (all 43 production), 168 selected segments, 70,155,126 events;
all 51 LH2 metadata entries accounted for. No jcache was invoked in this extension.
**The current user instruction is no jcache; historical staging instructions
below are superseded.** New counter defects in 4303/4305 and counter/clock
discontinuities in 4350 are documented. Prescaled residuals remain unresolved.
Production efficiency code/CSVs remain unchanged; PRE attribution is still unknown.

Updated 2026-09-09. Phase: SOURCE SET APPROVED; ALL SIX UPDATED SEGMENTS CHECKED;
TOTAL-LIVETIME HARDWARE ATTRIBUTION UNRESOLVED. No final total-livetime prescription selected.

## Latest full-run checkpoint

- User clarified PRE40/100/150/200 issue: UNKNOWN TRIGGER ATTRIBUTION, not missing
  branches. Nonzero values or rate correlation are not proof of physical routing.
  Do not use a PRE formula until input trigger, cable/module routing, widths,
  and trigger-population applicability are established.
- Tape request **100230264 completed**. All six updated segments inspected:
  3222394 contiguous events (1--3222394), 2233 real scaler snapshots. Primary
  EDTM/TRIG4/L1/clock/charge counters are monotonic across stitched segments.
- Fixed diagnostic current window 33.3625--45.1375 uA: N=A=2893284,
  S=2895694, D=112746, E_raw=112658, E_tight(+/-2 ns)=112650.
  EDTM=0.9991485286; physics CLT=0.9991685076. These are conditional component
  measurements, not a certified total-livetime correction.
- Primary convention [b_i,b_next); using (b_i,b_next] changes E_tight to112652,
  EDTM to0.9991662675 and CLTP to0.9991677890. Preserve this latch-boundary
  ambiguity. One terminal non-EDTM event is excluded by the primary convention.
- With current >2uA, CLTP=0.9992068726, consistent with database CPU_LT=0.9992
  at its stored precision; does not validate total electronic coverage.
- Full timing counts +/-2/3/5/10ns:112650/112651/112652/112652. Conditional
  binomial scales EDTM8.69e-5,CLTP1.73e-5; 10/20/60s time-block resampling
  spreads ~1.3e-4/~1.1e-4. These are diagnostics, not final uncertainties.
- Largest observed production-current loss burst near820.125614s: S2082,N1799,
  D80,E69 in~2s. Retained in averages. Shows time variation.
- Deliverables: `livetime_4398_check_20260909/REPORT.md`,
  `Run4398_livetime_investigation.pdf`, full/segment count tables and figures,
  reproducible scripts/macros, code/DAQ/replay provenance and checksums.
- Remaining: run-specific EDTM/PRE routing, TI buffer/readout/prescale evidence,
  original executable version, matched production charge/yield/good-helicity
  integration and final uncertainty/total-loss scope. User says logbook requires
  sign-in; no login attempted. No production code or output CSV changes made.
- `/tmp/nps4398_resume/full_run` has diagnostic binary column caches; scripts
  recreate them. Do not substitute `production` ROOT files or rerun tape retrieval.

## Earlier segment-0 checkpoint (full-run checkpoint above supersedes staging status)

- User explicitly approved resumption: "yeah, resume". Source-approval gate satisfied.
- User reconfirmed updated ROOT glob and said Hall C DAQ logbook requires sign-in.
- Results/scripts/figures: `livetime_4398_check_20260909/REPORT.md` and
  `segment0_diagnostic_figures.pdf`. These are interim diagnostics, not the final presentation.
- Only updated segment 0 was cached. Tape catalog lists 0--5. Requested missing
  segments 1--5 with `env -u DEBUG jcache get ...`; request **100230264**, active
  with all five pending at checkpoint. Check status; do not resubmit blindly.
- T=542584 events (g.evnum 1--542584), TSH=355 rows, TSHelH=19435 rows.
  PRE40/100/150/200 branches exist and are populated in T and TSH.
- Reproduced production E/D=19341/19310=1.0016053858. Final TSH row has updated
  evNumber but identical scaler values/time. Uncovered [541321,542585) contains
  1264 events and 49 raw-window EDTM events. Covered production-current ratio is
  19292/19310=0.9990678405. +/-2 ns corrected EDTM timing removes 2 outliers:
  19290/19310=0.9989642672. These are segment-0 pulse survival measurements only.
- Matched production-current S=495679, N=A=495221. Physics CLT (tight EDTM)
  (N-E)/(S-D)=0.9990805447. Current window 33.3625--45.1375 uA; covered time
  483.315153 s. Stored event current is from PRECEDING scaler interval.
- Current code's +/-500 "ns" is actually raw channels (~+/-48.83 ns).
- Found TSHelH evcount ranges applied to TSH evcount: different counters
  (19435 versus 354 maximum). Saved beam_time=163.950216 and gated denominator
  6550 therefore do not describe the intended common good intervals.
- Existing file discovery accepts any cached updated subset; no completeness check.
- Source job JSON found: `/group/nps/nps-ana/hcswif/jsons/NPS_PASS2_x60_4b_x60_2.json`.
  Input tar `/group/nps/nps-ana/nps_replay_pass2_github.tar.gz`, archived git HEAD
  2d462a85b860c54785d8287e644619a1cb1acc19. Job's stdout/stderr absent.
- Run_Data stores trigger-4 factor 1, but prescale set/read flags are 0: not
  independently decoded TI prescale proof. Ten embedded FADC/VTP configs saved;
  no TI buffer/ROC-lock setting located. Exact original executable still unproven.
- Independent ROOT/uproot equality verified: 542584 x 92 T values and 355 x 696
  TSH values; no mismatches. Current installed hcana source commit
  2f93d89f1ef482f36666a4e91821493ce4c22fc5 supports terminal-row behavior.
- Production efficiency source/CSV files untouched. Diagnostic caches (.npz/.bin)
  remain under `/tmp/nps4398_resume`; compact archive omits them.
- Next: when staging completes, inspect/stitch ALL updated segments using
  absolute event/scaler boundaries. Check duplicated initialization/end rows,
  resets, charge/yield matching. Establish actual TI mode, EDTM/PRE signal paths,
  NPS coverage, timing/background and covariance/systematics before choosing
  total livetime or making the requested final presentation. Do not simply
  multiply EDTM survival by CLT or transplant a historical unbuffered correction.
- Environment: `DEBUG=release` breaks jcache Java launch; `env -u DEBUG jcache ...`
  works. ROOT module and Python data reads required execution outside sandbox.

## Historical source-collection checkpoint

## User contract

- Reference glob: `/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_4398_*_1_-1.root`.
- Determine authoritative analysis livetime despite reportedly missing pretrigger branches and EDTM problems.
- FIRST collect internet documents and obtain user approval. ONLY THEN conduct complete ROOT/code/physics check.
- Archive relevant docs and compact memory in `/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc`.
- Very detailed presentation at investigation end, once conclusions supported. Ask about ambiguities; use tokens efficiently.

## Prior work supplied by user

- Presentation: `https://indico.jlab.org/event/946/contributions/16510/attachments/12602/20074/Luminosity_nps_collaboration_meeting_2025.pdf` (archived as `singh_nps_luminosity_202505.pdf`).
- Current implementation: `/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies`.
- Earlier implementation: `/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/src/nps_livetimes.h`.
- At the original source-collection checkpoint code inspection was deferred.
  Approval and subsequent inspection are recorded above.

## Completed

- Read applicable ancestor/destination AGENTS.md: efficiency behavior frozen; no physics edits made.
- Applied teach-back and token-saving skills; no agents spawned.
- Source catalog: `SOURCES_LIVETIME_4398.md`. Originals, URL/SHA-256 manifests, collection scripts and PDF text derivatives in `livetime_4398_sources_20260909/`.
- No ROOT file opened, branch inventory, event counting, scaler calculation, replay, or production job run.
- Missing pretriggers/EDTM issues remain user-reported, not independently verified.
- NPS wiki Computer and Electronic Deadtime page is a redlink (`page does not exist`), not evidence.
- Several web-tool DocDB opens failed; public metadata/PDF downloads succeeded without authentication. Do not describe successfully archived PDFs as inaccessible.

## Evidence rules

- Hierarchy: actual run hardware/DAQ configuration + pass2 provenance; measured file contents reconciled with code/maps; DAQ-expert and NPS documents; other-experiment precedents/models.
- Distinguish observed fact, inference, assumption, unresolved. Every result needs source/page/revision and reproducible command.
- Historical unbuffered EDTM correction is conditional; no transplantation into NPS without checking mode and signal path.
- Distinguish detector pretrigger, legacy width scalers, trigger-input counts and accepted events; branch names alone do not establish semantics.

## Full-check plan (approval now granted; progress recorded above)

1. Inspect prior work. Establish actual run trigger types/overlaps, prescale register versus factor, TI event types/masks, DAQ buffer/readout mode, EDTM injection and NPS VTP coverage, run-start flags/logs, exact pass2 code/config revision.
2. Inventory all matching segments, keys, trees, branches, metadata and event/scaler linkage. Distinguish missing, renamed, unexported, zero/unfilled channels. Do not assume T/TSH or any branch exists.
3. Define desired population and trigger-stage boundary. Separate intentional prescaling, detector/trigger efficiency, electronics loss, DAQ busy and offline cuts; live fraction versus reciprocal correction.
4. Map candidate numerator/denominator to hardware taps before/after busy, prescale and gating. Determine EDTM decode/rejection, timing windows, hit multiplicity/accidentals and trigger attribution.
5. Align good-beam intervals for event counts, scalers and charge. Check cumulative versus incremental counters, reset/rollover, split-file offsets/duplicates and incomplete first/last intervals.
6. Evaluate only identifiable estimators with explicit EDTM inclusion/subtraction and trigger/prescale conditions. Compare independent estimators. Never silently call CLT total livetime or invent missing counters.
7. Quantify finite EDTM statistics, periodic-pulser/beam-ramp bias, covariance, rate dependence, timing/current sensitivity and uncovered stages. No clipping unphysical results to [0,1].
8. If evidence insufficient, identify measurable quantities and exact missing evidence/additional data needed. Deliver signal-flow and branch-provenance tables, formulas, counts/plots, uncertainties and limitations.
9. Detailed presentation only at the end, with references and reproducibility appendix.

## Resume / reproduce

Read the resumed checkpoint and `livetime_4398_check_20260909/REPORT.md` first.
The original source catalog and source archive remain valid. Collection scripts
write to `/tmp/nps_livetime_4398_sources_20260909`; their manifests separate
download success from scientific applicability.

Before later ROOT/root-config/compile commands (NOT run yet):

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; root-config --version'
```

## Environment

Requested archive is a sibling workspace outside original writable roots; installing it requires filesystem escalation, separate from scientific approval. Initial sandbox curl failed DNS; subsequent sandbox shell launches failed `Failed to create unified exec process: No such file or directory (os error 2)`. Elevated commands worked. Patch tool can add staged files but could not read existing /tmp files for Update; complete replacement of our own staged notes used. No user files overwritten.
