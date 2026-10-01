# KinC_x60_4b LH2 livetime: active calculation

## FINAL UPDATE: one authorized raw file fully analyzed

Raw4305 arrived and completed:13,341,141,036bytes, CRC32=53e8a0c0,
Adler32=756ab5b0, both match tape stub. No other file staged/read.
`RAW4305_REPORT.md` is the latest detailed result; `REPORT.md` covers all ROOT runs.

- Raw extraction561307records,560647physics events,660scaler banks; all event
  numbers/timestamps/trigger masks exactly match ROOT. All660 raw clock/EDTM/L1/
  HMS-trigger scaler snapshots match ROOT. Final ROOT extra row repeats counts
  at boundary560648 (rawfinal560647), not additional exposure.
- Anomaly is in rawbanks125/126 (records105247/106294), events105080/106127.
  Slot12 clock increment5481, EDTM0; slot10 L1=2; slot11 TRIG4=2. Hodoscope
  slots6--9 continue counting. No ROOT-only counter/export failure.
- L1−recorded number afterdefect is−1046 in524snapshots,−1045 in11; before0/1.
  Persistent deficit, not recoverednextread. Still cannot distinguishhardware
  inhibit vsROC/FIFO/accumulation/readout failure withoutactualROC5 code.
- Slot7/ch31 ==slot9/ch31 all660snapshots; equalmappedEDTM first125, then+80
  inall535remaining. MaplabelsEmpty_12/Empty_24: NO attribution or replacement
  denominator adopted. Same cautionasPRE. `raw4305/validation.json` exactdata.
- Rawtextconfigs:5FADC+5VTP, noevent120/133prescalerecords. VTPinternalprescales
  areNOTHallCtriggerprescales. Types135/136areuuencodedanglephotos, notlogic.
- Userlogs4303/4305 confirmTI slaveblocklevel1/broadcastbuffer5busyenabled;
  notfullDAQmodeproof. NoEDTMrouting/masterreadback. Raw4303notread.
- AllROOTdata/figurescomplete. Newrawfigure`figures/raw4305_exposure_step.png`.
  NeedpersistcompactarchiveandfinalMarkdownlinks. Oldarchives/TeXunchanged.
- CLT_physics remains a candidate only. No baseline/truth/validationstandard.
  ExactEDTM/totalLTstillunestablished. Next evidence:actualROC5 scaleraccumulation
  code plus run-periodEDTM/masterwiring. No furtherstagingauthorized.

The active/pending raw status below is historical and superseded here.

## Latest state: ROOT calculation complete; one-file raw investigation active

- All 43 production runs / 183 updated segments / 77,838,384 events complete.
  Including controls/junk: 46 runs / 187 segments / 79,564,296 events.
  `finalize.done`, `validation.json`, `refresh_validation.json` verify completion.
- Final results/report: `/tmp/nps_livetime_refresh_20260909/REPORT.md`,
  `production_summary.csv`, six focused figures, interval/phase/statistics JSONs.
  Not yet copied into a new persistent archive; do that after raw investigation.
- **User says CLT_physics is NOT a baseline or truth: just a candidate definition.**
  All mutual comparisons must retain that interpretation; no estimator is
  established as a reference by agreement. Production correction stays unchanged.
- **New user authorization supersedes no-jcache for ONE FILE ONLY:**
  `/mss/hallc/c-nps/raw/nps_coin_4305.dat.0` (13,341,141,036 bytes).
  Submitted once via `env -u DEBUG jcache get ...`, request **100232289**.
  Last state active / file pending, cache placeholder zero bytes. Do not submit
  another get, and do not stage any other file. `DEBUG=release` breaks jcache's
  Java wrapper; unset DEBUG for status as well. Read-only status queries permitted.
- User supplied logbook logs for 4305 and 4303; exact bytes + provenance in
  `logbook/run4305_user_log.txt`, `run4303_user_log.txt`, `manifest.json`.
  Both explicitly identify their run; end counts 560,647 and 3,310,432 agree
  with corresponding ROOT event totals. Ten TI slave status blocks per log
  (prestart/end), block level 1 and broadcast block-buffer level 5, busy enabled.
  No EDTM/PRE routing named. Slave buffer depth does not alone establish whole
  DAQ buffering/deadtime: master/HMS busy can still constrain effective behavior.
- Local decoder source now found:
  `/u/group/halla/apps/analyzer/1.7.12/src/hana_decode/CodaDecoder.cxx` and
  `THaCodaFile.cxx`; corresponding installed headers in el9/RelWithDebInfo/include.
  ROOT/hcana environment still requires Hall C setup. Need inspect raw run-start
  config/prescale banks and scaler banks surrounding events105080--106127.
- Final pulse result: 3,650,941 broad-window tags; 51 lack edge interpolation.
  All 3,650,890 modeled tags pass 25 ticks. Tight cut rejects 264 of them;
  all264 pass25 ticks,226 pass10. Genuine pulse association is NOT causal-trigger
  proof. 9,395 positive tags outside raw window; 5,515 within500 ticks follow
  relative timing slope -3.9987366 ns/tick, r=-0.9997145: nearby pulse capture.
- All37 unprescaled run differences <1.7 paired-block sigma. Do not claim
  equivalence: common-mean diagnostic over35 counter-clean runs gives raw
  difference (4.3734 +/-1.4349)e-5 and tight (-3.0850 +/-1.5160)e-5. Shared
  systematics excluded; sign changes with timing. No baseline truth assumed.
- 52,025 exact-boundary selected intervals: only4303/4305 differ by >0.1s
  between event/scaler clocks. Their actual spans ~2s vs tiny scaler increments.
- 21 production sources unchanged; old402-file numerical and148-file Beamer
  archive hashes all pass. Presentation untouched. Independent uproot checks
  passed outside sandbox; stalled inside-sandbox duplicate was stopped.

The remainder below records earlier progress and is superseded by this section.

Updated 2026-09-09 evening (UTC 2026-09-10). User: resume full cached production
sample; pin down EDTM or alternative method; leave LaTeX untouched; maintain
this Markdown checkpoint. No jcache. No production correction edits.

## Frozen inputs and execution

- Work directory: `/tmp/nps_livetime_refresh_20260909`.
- `inventory.py` refreshed current source/CSV snapshots and ROOT size/mtime.
- 51 LH2 entries; 46 cached; 187 selected updated segments, including
  **183 segments across all 43 production runs**. Catalog comparison finds no
  production gaps. Remaining catalog gaps concern junk 4483/4552 only.
- `inventory_initial.json` and refreshed inventory have identical run/file
  records. Updated-first per run, never mix replay variants.
- `run_exports.py` runs three ROOT readers after sourcing Hall C setup;
  4397/4399/4402 last, with fresh size/mtime checks before dispatch.
- `.done` written only after T and TSH exports; `dispatch.jsonl` records phases.
- `analyze_ready.py` computes completed runs while export continues; logs and
  `analysis_progress.json` are the live progress record.
- Audit compares saved selections only for the same replay path identity;
  changed replay variants cannot be treated as exact saved-file controls.

Commands actually run:

```bash
cd /tmp/nps_livetime_refresh_20260909
python3 inventory.py > inventory_refresh.log 2>&1
python3 run_exports.py > dispatch.log 2>&1
python3 catalog_coverage.py
python3 analyze_ready.py > analysis.log 2>&1
```

The long commands run in background tool sessions. Do not launch duplicate
export workers on resume. Inspect logs, `.done` markers and active processes.

## Questions to resolve with the new export

1. **Coverage/current:** repeat old/current and common-interval raw/tight ratios;
   compare against prior 168-segment baseline with replay changes labeled.
2. **EDTM acceptance:** corrected-time tails can contain genuine pulses in other
   configurations. Inspect raw/corrected times, multiplicity, trigger timing,
   and independent event-timestamp periodicity. A narrow peak is not efficiency.
3. **Prescaling:** quantify pE/D fluctuations under stated models; compare to
   block variation and all-trigger/physics ratios. Do not force unity or claim
   periodic bias without evidence. Test time/rate dependence.
4. **Counters:** inspect 4303/4305 flagged intervals and 4350 restart/overflow
   evidence in the updated replay. Separate N>A inconsistency from real S>N
   losses with N=A. No inferred counter repair or unmatched exclusions.
5. **Alternative:** matched recorded/L1 versus configured trigger scalers gives
   a DAQ-stage estimator subject to counter integrity, pulser routing and
   prescale conditions. It does not measure upstream electronics/NPS losses.

## Fixed scientific conditions

### Timestamp pilot (new measured evidence)

`python3 timestamp_probe.py 4253 4255 4259 > timestamp_pilot.log 2>&1`
fits the pulser period from g.evtime and alternates core events between local
interpolation anchors and held-out validation. Period ~6257191.7 event-clock
ticks. Held-out absolute residual 99.9th percentile ~7--11 ticks. An initial
20-second straight-line model failed to follow oscillator wander; replaced by
neighboring-anchor interpolation before drawing conclusions.

- Run 4253: all 13 raw-window candidates rejected by +/-2 ns lie within 10
  timestamp ticks of predicted pulser phase. Run 4255: all 3 do too. Examples
  have corrected times 265--328 ns despite the ~245-ns core. This supports
  genuine accepted pulser tails; tight timing is not an unbiased classifier.
- Run 4259: 6632 tight tags remain; none of 7 additional positive-raw tags is
  within 125 ticks. Its pE/D=1.03042168 survives this discrimination check.
- Prescaled EDTM fluctuations have large conditional errors. Previous-baseline
  run 4259 excess is 2.73 binomial sigma relative to L=1, not an automatic
  hardware-bias diagnosis. Final refreshed sample/statistical tests pending.
- Phase windows remain diagnostic; account for interpolation support, held-out
  validation and accidental events before using any phase-derived E.

N=recorded; E=accepted EDTM; S=selected trigger input; D=sent EDTM; A=L1.
On identical intervals: pE/D, pN/S, p(N-E)/(S-D), pA/S. Accepted E is subtracted
from N; sent D from S. Shared-count identity prevents treating these as
independent measurements. PRE attribution remains unknown, so no PRE correction.
Common current cuts do not establish physics/charge weighting or full
good-helicity/yield matching. No certified total-LT prescription yet.

Earlier immutable archives and research: sibling
`livetime_KinC_x60_4b_lh2_20260909` and `livetime_beamer_20260909` in `_main/misc`.
The existing LaTeX/PDF will not be updated during this investigation.

## Continuation findings (full-sample jobs still running)

- Exporter now reads selected TBranch objects directly. Exact binary agreement
  with the original TTree reader: 58,600 x 19 values in 4259 and 450,326 x 19
  values in 4253 segment 0. Completed exports retained; own readers restarted
  with append-only logs. No physics definition changed.
- At 29 completed production runs / 123 segments, all 2,558,067 broad raw-window
  candidates with an interpolation model are within 25 event-clock ticks of
  pulse phase; 31 have no model at file edges. All 195 raw candidates rejected
  by corrected +/-2 ns pass this independent phase check (164 within 10 ticks).
  This is direct NPS evidence against adopting the tight corrected cut.
- Positive raw hits outside the broad window include near-phase candidates:
  they require a separate acceptance/background study. A uniform-phase
  sideband subtraction is unjustified: among 46,419,270 raw-zero events, none
  lies within 25 ticks; 506 lie within 125 ticks. Trigger/DAQ correlations matter.
- Subtracting the selected trigger's corrected time does not fix the tails:
  relative +/-2 ns recovers only 1/14, 1/3, 0/0, 1/6, 0/4, 3/12 whole-file
  raw timing tails for 4253,4255,4259,4301,4305,4350 respectively.
- `event_clock_checks.py`: exact event-number matches at scaler boundaries;
  median healthy-interval clock ratio, not a fit across counter discontinuities.
  Of 34,938 nominal intervals checked so far, only 4303 and 4305 disagree by
  >0.1 s. Their event spans are 2.000258 and 2.001508 s while scaler spans are
  0.000869 and 0.005481 s. This establishes missing scaler exposure relative to
  event timestamps, not whether hardware inhibition or decoding caused it.
  A merely absent cumulative snapshot would normally recover counts on the
  next read: do not assert a missing bank as the specific cause.
- Comparing EDTM to the scaler clock alone misses these defects because both
  lose exposure together. No counter repair or unmatched yield exclusion.
- `acceptance_tests.py`: compare pE/D and p(N-E)/(S-D) using an explicitly
  exchangeable-label hypergeometric null plus paired 20-second block bootstrap.
  Run 4259 p=5: difference 0.0351962, hypergeometric z=2.97, block z=2.61.
  A single prescaled >1 estimate is not proof of a hardware bias. Periodic
  pulser labels need not satisfy the exchangeability model; report that caveat.
- Raw cache directory checked for 4259,4303,4305,4350,4398,4493: no EVIO files.
  No staging requested. Run-start EVIO/configuration and authenticated logbook
  remain the external evidence needed for routing/mode and counter root cause.
