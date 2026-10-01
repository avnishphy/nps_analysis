# KinC_x60_4b LH2 livetime investigation

## Detailed research record and apps search

[Living research log and audit trail](LIVETIME_RESEARCH_LOG.md) records prior
tests, definitions, evidence, uncertainty and dated continuation decisions.
[Apps-search evidence](nps_roc_apps_research_20260910/README.md) preserves exact
commands, source snapshots and Git revisions. `/u/group/nps/apps` contains
replay/monitoring software; online ROC5 readout code remains unlocated.
ROC5 and npsvme5 are distinct components. ROC5's own end counts agree with
recorded totals but do not establish integrity of its scaler payload.

## CODA documentation follow-up

UTC2026-09-10: [DAQ/ROC evidence and definitions](coda_livetime_research_20260910/REPORT.md).
TI3v11.3 documentation confirms broadcast buffer5 is selected; local setting1
is not the active threshold. All five end TI slave Readout/Ack/L1A counts
match each run's recorded total. This does not measure pre-acceptance losses.
Logged timers lack established common exposure/reset/latch history; no new
ratio adopted. Saved coin_sparse ROC5 rol1/rol2 entries and matching readout
source are concrete next artifacts. Exact EDTM/PRE routing still unknown.
Public sources, register checks and provenance archived; LaTeX unchanged.

## Latest continuation: full cache, one raw file and logbook checks

UTC 2026-09-10: [current checkpoint](LIVETIME_REFRESH_CHECKPOINT.md),
[full ROOT report](livetime_refresh_20260909/REPORT.md),
[raw/logbook findings](livetime_refresh_20260909/RAW4305_REPORT.md),
[43-run table](livetime_refresh_20260909/production_summary.csv).
All43 production runs /183 updated segments /77,838,384 events calculated.
User authorized exactly one raw file,4305.dat.0; it was staged and validated.
The 4305 scaler deficit is present in raw banks, while every recorded physics
event matches ROOT. Both user logs are preserved with provenance.

**CLT_physics is a proposed definition, never a baseline or truth.** Pulse
phase distinguishes associated EDTM hits but does not identify causal triggers.
PRE and newly identified pulse-like unassigned channels remain unused. Exact
EDTM routing/total LT remains unresolved; production correction unchanged.
LaTeX/PDF presentation untouched, as requested. Earlier archives remain frozen.

## Previous continuation: LaTeX presentation and source research

2026-09-09: [technical Beamer PDF](livetime_beamer_20260909/KinC_x60_4b_LH2_livetime_Beamer.pdf),
[editable LaTeX](livetime_beamer_20260909/livetime.tex),
[research/definitions](livetime_beamer_20260909/RESEARCH.md) and
[build instructions](livetime_beamer_20260909/README.md).
56 slides: 27 main and 29 backup. The initial 168-segment results are preserved;
the new user-staged coverage has not yet been recalculated in these figures.

The tight corrected-time EDTM window is explicitly diagnostic: a PionLT source
provides a warning about genuine timing tails, not proof of NPS acceptance.
NPS-period trigger documentation identifies beginning-of-run EVIO configuration
records as a useful next lead. PRE routing, TI mode, prescale sampling and total
loss coverage remain unresolved. New public source bytes/hashes and a later
cache filename/size snapshot are archived. No jcache or production edit performed.

Several formerly missing updated segments have appeared during user staging.
Future work must freeze a fresh per-run replay manifest, retain raw/tight timing
variants and old/same-cache controls, and repeat interval/counter audits. Do not
mix updated and production replays. This Beamer package supersedes the earlier
presentation for discussion; the original numerical archive stays immutable.
Temporary column caches were absent in this resumed session; regenerate them
from the frozen manifest when reproducing the ROOT calculation.

Updated 2026-09-09. Cache-only calculation, report and presentation complete.
**No certified total-livetime prescription selected. Production code/CSVs unchanged.**

## Deliverables

- [Report](livetime_KinC_x60_4b_lh2_20260909/REPORT.md)
- [Presentation](livetime_KinC_x60_4b_lh2_20260909/KinC_x60_4b_LH2_livetime_investigation.pdf)
- [All-run CSV](livetime_KinC_x60_4b_lh2_20260909/run_summary.csv)
- Archive includes interval tables, source/baseline snapshots, per-run JSONs,
  standalone figures, reproducible scripts, read/metadata checks and SHA256SUMS.

## Scope and constraints

- All 51 LH2 metadata entries: 43 production, six junk, two efficiency.
- 46 runs calculated from 168 selected cached segments: 70,155,126 events.
- No cache for junk runs 4307, 4486, 4552, 4553, 4554; no LT invented for them.
- **Do not invoke jcache or retrieve tape data.** This overrides historical
  staging instructions in the run-4398 memory.
- Existing updated-first / production-fallback policy retained per run. 253
  cached files include 85 alternate replay versions; do not combine versions.
- Preserve explicit partial coverage and separate components at event gaps,
  counter restarts and clock restarts. Catalog availability is not exposure.
- PRE40/100/150/200 trigger attribution is unknown. No PRE correction applied.

## Verified findings

- All 42 production runs with unchanged saved coverage reproduce original
  NewGen E, D and raw peak exactly. All 168 selection/file/metadata records
  checked; 40 stored prescale records populated and agreeing, 128 all-zero
  placeholders, all set/read flags zero. Every event includes its configured
  trigger bit. None is independent run-start TI/routing proof.
- Run 4398: saved one-segment NewGen=1.00160539; same-cache six-segment
  NewGen=0.99858088; matched +/-2-ns EDTM=0.99914853; physics CLT=0.99916851.
- Run 4259, p=5: matched EDTM=1.03042168; all-trigger CLT=0.99925323;
  subtracted physics CLT=0.99522546. Sampling/routing residual unresolved.
- Run 4493, p=2: old=0.98521430, matched=0.99282550. Do not force unity.
- Selected counter defects: 4303 has N=1915,A=S=3,D=1,E=79 in 0.000869 s;
  4305 has N=1047,A=S=2,D=0,E=80 in 0.005481 s. An independent uproot
  reader verifies representative ROOT values and the 4305 event mismatch.
- Diagnostic exclusion of those intervals is tabulated but NOT adopted.
  Real S-N loss bursts with N=A are retained. Yield/charge consistency remains
  mandatory if an interval is eventually removed.
- 4350: zero-clock and near-2^32 EDTM counter anomalies occur below the current
  window; EDTM cumulative restart at segments3->4, clock restart at4->5.
  These boundaries are not stitched; raw evidence is retained.
- Same-cache old/new normalization illustration spans -0.76662% (4493) to
  +0.69710% (4308), for fixed yield/charge and changing ONLY L. This is not an
  observed cross-section bias or a certification of the new total correction.

## Remaining scientific conditions

Run-period EDTM/PRE routing, TI mode and decoded prescales, exact replay
executable, anomalous scaler-block cause, physical latch convention,
production good-helicity/yield/charge matching and final uncertainty/total-loss
coverage. A tighter pulse window, plausible ratio or flatter yield does not
establish those conditions. Do not multiply overlapping EDTM and computer LT.

## Reproduction

Exact commands and prerequisites are in the report. Working column caches
were written to `/tmp/nps_x60_4b_lh2_livetime_20260909/columns`; they were absent
in the resumed session. The compact archive omits
the binaries. The archived macro recreates them from the frozen cache manifest
without retrieval. ROOT commands source `/group/nps/singhav/setup.csh`.
The independent uproot reader required execution outside the sandbox after its
sandboxed read stalled; normal diagnostic ROOT exports completed in sandbox.
