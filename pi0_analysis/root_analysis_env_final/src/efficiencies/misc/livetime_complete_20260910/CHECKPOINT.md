# KinC_x60_4b LH2 livetime: resume checkpoint

Updated UTC 2026-09-10. Full ROOT calculation and one authorized raw-file
investigation complete. **Exact EDTM/total-livetime prescription still unresolved.**

## User constraints

- **CLT_physics is only a proposed definition. Never use it as baseline, truth,
  calibration target or validation standard.** Mutual agreement proves neither
  estimator. Production correction/code/CSVs unchanged.
- The user now explicitly authorized all documentation and LaTeX updates.
  The earlier hold is superseded; it is not evidence that the method is final.
  Earlier immutable archives are preserved. New edition and checks below.
- User authorized **one raw file only**, already used:
  `/mss/hallc/c-nps/raw/nps_coin_4305.dat.0`, request100232289.
  No more staging authorized. File is now cached and verified. Do not resubmit.
- Work independently, keep notes, use tokens efficiently. No subagents.

## Current documentation edition

- [Index/edit/build instructions](README.md)
- [Report](REPORT.md), [references](REFERENCES.md), [plot catalog](PLOT_CATALOG.md)
- [Report PDF](KinC_x60_4b_LH2_livetime_report.pdf), [slides PDF](KinC_x60_4b_LH2_livetime_presentation.pdf)
- Editable report.tex, sections/ and presentation.tex; run `bash build.sh` in
  the edition directory. No ROOT is required. Do not regenerate Markdown-to-TeX
  sections over direct user LaTeX edits without preserving them.
- All43 current count/ratio rows, nineteen current figures, complete frozen
  evidence/source packages, old editions and original top-doc backups included.
- Eleven new plots read existing summaries/intervals; seven current diagnostic
  pairs and the supplied NewGen PNG copied unchanged. No event reread or staging.
- Core pE/D method was already the user's. Same-file refinements concern
  exposure/current/joins. For4398 E112586 ->112263 ->112658 ->112658, D112746.
- No exact total-LT prescription, production replacement or new raw allowance.
- Build, numerical and preservation checks: `validation.json`.

## Artifacts / execution

- User requested detailed ongoing audit: `LIVETIME_RESEARCH_LOG.md` is the
  living research record; append dated evidence/commands/decisions as work
  proceeds. Short checkpoint does not replace the detailed record.

- Persistent compact archive: `evidence/full_cache/` beside this note.
  Start with `REPORT.md`, `RAW4305_REPORT.md`, `production_summary.csv`.
  Seven standalone diagnostic figures, all interval/phase/raw tables, source
  snapshots, user logs, scripts and SHA256SUMS included. No full raw copy.
- Scratch: `/tmp/nps_livetime_refresh_20260909`; no jobs need restarting.
  Large ROOT column caches remain there and are omitted from the compact archive.
  Recreate from frozen manifest if absent; commands in REPORT.md.
- ROOT environment: `csh -c 'source /usr/share/Modules/init/csh; source
  /group/nps/singhav/setup.csh; hcana ...'`. jcache wrapper requires `env -u DEBUG`.

## Established numerical findings

- All43 production runs /183 updated segments /77,838,384 events.
  With controls/junk:46 cached runs /187 segments /79,564,296 events.
  No production catalog gaps; all187 input size/mtimes stable. All11 exact
  saved-file controls reproduce old NewGen counts. Changed replay/coverage
  comparisons are explicitly separated from method changes.
- Shared scaler intervals, ending-snapshot current, no gap/restart bridging.
  N=recorded, E=EDTM tags, S=configured input, D=sent EDTM, A=L1, p=prescale.
  Compare pE/D, pN/S, p(N-E)/(S-D), pA/S conditionally, without a truth baseline.
  Exact mixture identity verified. Yield/good-helicity/charge matching remains.
- 3,650,941 broad raw-window tags:51 lack edge phase interpolation; all other
  3,650,890 pass25 event-clock ticks. Tight corrected +/-2ns rejects264 tags;
  all264 pass25ticks,226 pass10. Held-out core validates the local phase model.
  Pulse association does NOT prove the pulse caused the recorded trigger.
- 9,395 positive hits outside raw window;5,515 within500ticks follow relative
  timing slope-3.9987366ns/tick (r=-0.9997145), consistent with nearby pulse
  capture. Don't count all positive hits or apply a uniform phase sideband.
- 4259,p5: EDTM1.03042168, proposed subtracted ratio0.99522546; difference
  z=2.97 under exchangeable-label hypergeometric model,2.61 paired20s bootstrap.
  Not a hardware-bias proof. Other prescaled candidate differences <1sigma.
- Per-run unprescaled differences <1.7 paired sigma, but common-mean diagnostic
  over35 counter-clean runs gives raw (4.3734+/-1.4349)e-5, tight
  (-3.0850+/-1.5160)e-5. Timing changes sign; do not claim precise equivalence.
- 52,025 exact-boundary intervals: only4303/4305 have >0.1s event/scaler clock
  disagreements. Both about2s versus scaler0.000869/0.005481s. No repair adopted.
- Fixed-yield/charge old-to-matched-raw sensitivity spans-0.76662% to+0.69107%.
  Illustration of potential artifacts, not measured cross-section bias.

## Raw 4305 / logbook evidence

- One raw file,13,341,141,036bytes; CRC32=53e8a0c0, Adler32=756ab5b0 match tape.
  561307records,560647physics events,660scaler banks. All physics event numbers,
  timestamps/masks and all relevant raw scaler snapshots match ROOT exactly.
- Rawbanks125/126 (records105247/106294), events105080/106127: clock increment
 5481ticks,EDTM0,L1=2,TRIG4=2 while1047events/80tags recorded. Hodoscope
  slots6--9 advance. Deficit already in raw, not created by ROOT/export.
- L1−event number afterdefect is-1046 in524snapshots,-1045 in11, before0/1.
  Persistent count deficit; cause could involve hardware/ROC accumulation/FIFO.
- Slot7/ch31 ==slot9/ch31 all660snapshots: equalmappedEDTM first125, then+80
  inall535remaining. Maps label Empty_12/Empty_24. **Not replacement denominators.**
  Attribution remains unknown, as for PRE scalers.
- User logs4303/4305 under `logbook/` confirm their event totals and TI slave
  blocklevel1/broadcastbuffer5/busyenabled. Not effective whole-DAQ-mode proof.
  No explicit EDTM/PRE wiring or master prescale readback. Raw4303 not read.
- Raw contains5FADC+5VTPtextconfigs; VTPinternalprescales are not HallCprescales.
  No types120/133prescale records. Types135/136 hold angle photos, not logic.

## Next evidence / work

NEW dated wiring evidence (supersedes blanket "EDTM route unknown"):
`evidence/routing/REPORT.md`, excerpts and43-run mapping CSV.
User supplied log4170042(2023-08-23) and4171556(2023-08-29), confirmed applicable.
Use as working map; web401 token_expired prevented independent retrieval.
ROC1 TI1=NPS1cluster OR2cluster ORdelayedEDTM;TI3=HMS3/4;TI4=HMS ELREAL;
TI5=1AND3;TI6=1AND4. Local NPS/cable T6 (2cluster) is NOT ROC1 TI6(coincidence).
EDTM copy joins NPS at CH after~3.5us delay, bypassing NPS cluster formation;
VLD traverses FADC/VTP/V1495 instead. HMS discriminator fanout/ELREAL pulser
path and sent-scaler tap still not fully specified; PRE remains unmapped.
Metadata check:5 TI6 runs4253/54/55/56/59 (2706910events),38 TI4 runs
(75131474events); all43 sole enabled PS_N factors match frozen summary.
Interpret pE/D by cohort: injected coincidence vs HMS singles. No automatic
NPS hardware-trigger correction for TI4; no EDTM validation of bypassed NPS
formation for TI6. Exact total LT unresolved;4303/4305counterdefect unchanged.

Focused NPSlib follow-up: `evidence/npslib/REPORT.md`, also appended
to living research log. /group and /u/group NPSlib paths same directory.
THcNPSCoinTime.h:78 says T1=NPS VTP||EDTM,T3=HMS3/4; comment from2023-09-12,
not verified wiring. Timing implementation uses configured entries0,1,3.
Event137 parser uses unique-key map emplace: first occurrence retained per
instance;4305 raw text has32 VTP prescale entries per VTP and differing
latencies. Keep original per-record/indexed configs; no full table from one
GetInfo value. Source VTP bits0--5 labeled singles/cosmic-scin/cosmic-column/
different-crate-pair/same-crate-pair/VLD, not global TI bits or EDTM causality.
Decision-time mask0x7F conflicts with11-bit comment: verify firmware before
using VTP decision timing. No correction/code change or extra raw read.

Apps lead checked: `/u/group/nps/apps` has hcana/NPSlib/panguin, no online ROC
code found in bounded2562-file search. Evidence/source hashes/commands:
`evidence/apps_search/`. Three hcana1.0.0/1.0.1/1.0.2 scaler handlers
identical and clean; these decode raw banks offline, not hardware FIFO code.
Neighboring nps_replay module references `/cdaqfs1/apps`, absent on this node.
Both logs distinguish `ROC5` from `npsvme5`; five NPS TI dumps must NOT be
assumed to describe the HMS scaler ROC. ROC5's own end events match totals:
4303:3310432 (log8025),4305:560647 (log5173). End event counts do not certify
per-channel scaler values. Exact online source/build/config remains missing.

Official CODA follow-up: `evidence/coda/REPORT.md`, hashed
public sources and reproducible `check_logs.py`. TI3v11.3 header matches logged
firmware label: all20 register pairs select broadcastbuffer5, busy enabled;
local blockBuffer1 is NOT the selected threshold. All five end TI slave
Readout/Ack/L1A triplets match recorded totals in each run (4303:3310432;
4305:560647). Accepted/readout agreement is not a trigger-loss denominator.
Initial/end timers are not certified cumulative exposure: nonzero initial
timers with L1A0, busy decreases for matching board0x7101150e in both runs.
Two interleaved4305 timer chunks deliberately unassigned. Need latch/reset
history; no TI timer ratio adopted. Current TI hardwarePDF updated2026-08-27,
so later features cannot be assumed for2024. Exact tiStatus implementation
not obtained. Published docs readable at /group/da/distribution/coda.
CODA configuration example identifies saved COOL_HOME configuration ROC rol1/
rol2 paths: retrieve actual coin_sparse ROC5 entries and matching scaler
source/build. No exact NPS EDTM/PRE mapping found on CODA site.

Actual ROC5 scaler readout/accumulation code for these runs; run-period EDTM
injection/trigger wiring and master/TS settings. Investigate the raw deficit
without asserting its hardware cause. Excluding invalid intervals requires
the identical physics-yield and charge selection; not implemented here.
Keep broad/tight EDTM and other definitions diagnostic until their physical
endpoints and acceptance are established. Do not force unity, fit a correction
to CLT_physics, use unassigned/PRE channels, or multiply overlapping factors.
