# KinC_x60_4b LH2 livetime: resume checkpoint

Updated UTC 2026-09-10. Full ROOT calculation and one authorized raw-file
investigation complete. **Exact EDTM/total-livetime prescription still unresolved.**

## User constraints

- **CLT_physics is only a proposed definition. Never use it as baseline, truth,
  calibration target or validation standard.** Mutual agreement proves neither
  estimator. Production correction/code/CSVs unchanged.
- No LaTeX/presentation update until the method is established. Old Beamer
  archive's148 hashes and numerical archive's402 hashes verified unchanged.
- User authorized **one raw file only**, already used:
  `/mss/hallc/c-nps/raw/nps_coin_4305.dat.0`, request100232289.
  No more staging authorized. File is now cached and verified. Do not resubmit.
- Work independently, keep notes, use tokens efficiently. No subagents.

## Artifacts / execution

- Persistent compact archive: `livetime_refresh_20260909/` beside this note.
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

Actual ROC5 scaler readout/accumulation code for these runs; run-period EDTM
injection/trigger wiring and master/TS settings. Investigate the raw deficit
without asserting its hardware cause. Excluding invalid intervals requires
the identical physics-yield and charge selection; not implemented here.
Keep broad/tight EDTM and other definitions diagnostic until their physical
endpoints and acceptance are established. Do not force unity, fit a correction
to CLT_physics, use unassigned/PRE channels, or multiply overlapping factors.
