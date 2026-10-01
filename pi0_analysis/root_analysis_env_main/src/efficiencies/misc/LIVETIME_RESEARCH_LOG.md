# KinC_x60_4b LH2 livetime: research log and audit trail

Started UTC 2026-09-10 at the user's request for a detailed, auditable research
document. This is a living record; the linked numerical/source archives remain
frozen. Entries below distinguish results recovered from previous reports from
checks executed in this continuation. **The exact total-livetime correction is
still unresolved. No production correction has been changed.**

## Current documentation edition

[Complete report, presentation, references and plots](livetime_complete_20260910/README.md)
are now updated at the user's explicit request. The final entry below records
this documentation build; earlier dated entries retain their historical scope.

## Reading and maintaining this record

The [resume checkpoint](LIVETIME_REFRESH_CHECKPOINT.md) is a compact handoff.
This document records questions, methods, evidence, alternatives and decisions.
For every continuation, add the date, exact inputs/source revisions, commands
or executable script, outputs, interpretation and remaining uncertainty.
When a conclusion changes, add an explicit correction rather than silently
rewriting the old evidence. Distinguish source-code behavior from behavior
verified in the run. Filenames, directory dates and current module defaults
do not prove the software/configuration loaded during a particular run.

Initial evidence packages (later follow-ups and current edition appear below):

| Package | Role and present status |
|---|---|
| [Initial numerical report](livetime_KinC_x60_4b_lh2_20260909/REPORT.md) | Earlier cache coverage; retained for comparison, not current all-run inventory |
| [Presentation research](livetime_beamer_20260909/RESEARCH.md) | Earlier source research and presentation; frozen historical edition |
| [Full-cache report](livetime_refresh_20260909/REPORT.md) and [raw report](livetime_refresh_20260909/RAW4305_REPORT.md) | Current all-production calculation, timing/exposure tests, raw4305 and user logs |
| [CODA follow-up](coda_livetime_research_20260910/REPORT.md) | Version-qualified documentation, register checks, timer limitations |

## Non-negotiable analysis constraints

- **CLT_physics is a proposed definition, never truth, a baseline, calibration
  target or validation standard.** Agreement between candidates is only a
  consistency test. Historical output labels do not override this restriction.
- The dated logbook-routing entry below now supplies a working EDTM route;
  detailed injection fanout/taps and PRE40/100/150/200 attribution remain
  unresolved. Do not infer wiring from names, rates or near-identical counts.
- User authorized one raw file; the allowance was used for
  `/mss/hallc/c-nps/raw/nps_coin_4305.dat.0`, request100232289. No further
  staging/resubmission. Cached ROOT production coverage has been processed.
- The user has now authorized all documentation and LaTeX updates. The earlier
  hold is superseded. No numerical repair, event exclusion or production
  code/CSV change has been adopted.
- Any future exclusion must apply identically to physics yield, charge and
  livetime exposure. A plausible livetime curve alone is insufficient.

## Previous work: reconstructed from the frozen reports

This section summarizes established prior work; it is not a claim that the
full ROOT or raw analysis was rerun during today's source search.

### Inputs and counting method

43 production runs,183 updated segments,77,838,384 events. Including available
controls/junk:46 runs,187 segments,79,564,296 events. The frozen file manifest
and source snapshots identify the calculation inputs. All11 old-result
controls with exactly matching file coverage reproduced original NewGen counts.
Other old/new differences can include replay or coverage changes; they cannot
be attributed solely to changing the livetime method.

The comparison uses shared scaler intervals `[event_i,event_(i+1))`, ending
snapshot BCM4A current, the existing per-file current-peak band (+/-15%), and
joins only across compatible counter/event continuity. Heads, tails, gaps and
restarts are handled explicitly. This defines diagnostic exposure; complete
physics-yield, good-helicity and charge selection still needs matching.

Let N=recorded events, E=EDTM-tagged subset, S=selected trigger input counts,
D=sent EDTM counts, A=L1 counts and p=prescale factor. Compared quantities are
`pE/D`, `pN/S`, `p(N-E)/(S-D)` and `pA/S`. Their algebraic relation

```text
pN/S = (1-D/S) * p(N-E)/(S-D) + (D/S) * pE/D
```

was checked numerically. They share counts and are not independent estimates.
To interpret the subtraction physically requires verified pulse contribution
to S and correct identification of the accepted pulse subset E. Prescaling,
overlap and trigger causality remain conditions, not established assumptions.

### Attempt, evidence and resulting decision

| Attempt/question | Evidence in full-cache report | Decision/limitation |
|---|---|---|
| Use original NewGen across refreshed coverage | Exact-coverage controls reproduce; other inputs changed | Preserve same-cache comparisons and distinguish coverage from method |
| Match event/scaler exposure | Shared boundaries/current and compatible joins are implemented | Necessary counting correction; does not certify total physics survival |
| Tight corrected EDTM peak +/-2ns | Rejects264 broad-window tags; all264 have pulse-phase association within25 event-clock ticks | Narrow timing does not establish acceptance; keep timing sensitivity |
| Recover tails with EDTM-minus-trigger timing | Does not universally recover the raw-window tails | Not a demonstrated replacement acceptance |
| Count all positive EDTM hits | Thousands show timing consistent with a nearby pulse captured in another trigger's TDC window | Do not equate positive TDC hits with independently accepted pulses |
| Uniform phase sideband | Background is depleted around pulse phase | Flat-background subtraction is not justified |
| Compare EDTM scaler with scaler clock | Agreement can survive simultaneous exposure loss | Cross-check against recorded event timestamps and accepts |
| Interpret prescaled EDTM greater than one as hardware failure | 4259,p5 gives1.03042168; discrepancy depends on explicit statistical model | Sampling can exceed one; no forced-unity correction |
| Use PRE/unassigned pulse-like scalers | Physical connections not verified | Leave unused as denominators |
| Repair missing scaler counts or exclude intervals | Raw deficit confirmed; cause and matched-yield treatment unresolved | Neither action adopted |

There are3,650,941 broad raw-window tags:51 lack edge phase interpolation; all
remaining3,650,890 pass25ticks. This establishes association, not which signal
caused the readout. Raw-window width is +/-500 raw channels (about +/-48.83ns),
not +/-500ns. A periodic pulser samples availability at pulse times; physics
acceptance weights the actual eligible-event population. Their equivalence
requires physical justification.

### Raw4305 localization

The single cached raw file was read completely and verified against the tape
stub (13,341,141,036bytes; CRC32=53e8a0c0, Adler32=756ab5b0).
All560,647 physics event numbers/timestamps/trigger masks and all660 relevant
scaler snapshots match ROOT. One extra final ROOT scaler row repeats counter
values; it does not supply additional exposure.

Across events105080 to106127,1047 events and80 EDTM tags were recorded over
about2.001508s, while mapped scaler increments are clock5481ticks,EDTM0,L1=2,
TRIG4=2. The raw module-header sequence is intact; slots6--9 advance while the
mapped clock/trigger group scarcely advances. The L1 count deficit persists
through the remaining raw snapshots. It is already in the raw record.

This rules out a ROOT-only read/export conversion as the cause of that defect.
It does not distinguish hardware inhibit, lost FIFO contents, online
accumulation/readout loss or another cause. A simply omitted cumulative
snapshot would normally recover its counts at the next snapshot; that did not
happen here. Raw4303 was not read; its corresponding ROOT clock discrepancy
is2.000258s versus0.000869s. Avoid extending raw4305 verification to4303.

## 2026-09-10: CODA documentation follow-up

The [CODA evidence report](coda_livetime_research_20260910/REPORT.md) contains
public URLs, exact downloaded source hashes and reproducible log checks.
TI3v11.3 header definitions agree with the firmware label in the logs.
All20 initial/end register pairs select broadcast buffer5 with busy enabled;
local blockBuffer1 is not the selected threshold. Other crates/master/busy
paths can limit effective whole-DAQ acceptance.

All five end TI status dumps per run have Readout=Ack=L1A matching the recorded
totals. They do not measure requests lost before acceptance. The initial/end
timer snapshots do not have verified common reset/latch exposure: initial
timers are nonzero with L1A0, and busy decreases for the same unambiguously
paired board in both logs. Do not subtract or use endpoint timer ratios.
The current TI hardware manual is dated2026-08-27; later features are not proof
of2024 behavior. Exact online tiStatus implementation remains unlocated.

CODA examples identify ROC `rol1`/`rol2` paths in the saved configuration.
These are the concrete artifacts to retrieve for `coin_sparse`, together with
their source/build and the actual scaler driver. Generic examples are not NPS
configuration evidence.

## 2026-09-10: search of /u/group/nps/apps

**Question:** does the user-suggested installation contain the online ROC5
readout list or a path to it?

**Execution/evidence:** [source-search archive](nps_roc_apps_research_20260910/README.md).
`evidence/commands.json` records exact argv, shell equivalent, UTC, exit codes
and output sizes. Full stdout/stderr, file inventory, Git revisions and six
source snapshots are preserved. No remote login, DAQ command or ROOT job ran.

Observed contents: `hcana`, `NPSlib`, `panguin` and setup scripts. The bounded
inventory contains2562 files, including hidden/ignored paths but excluding
Git internals and not following directory symlinks. The precise text search
found no coin_sparse, COOL_HOME, rol1/rol2, SIS3801 online-driver, tiLib or
/site/coda references in the searched non-build text. An earlier broad tiLib
pattern matched compiler `multilib` text; those false positives were inspected
and the boundary-aware query is recorded separately. Absence is limited to
this search scope; it does not prove no other installation/archive exists.

`THcScalerEvtHandler.cxx` is available under hcana1.0.0/1.0.1/1.0.2. The three
files are byte-identical and clean relative to their respective Git HEADs.
HEADs are207fd147dc4b4a460d27c458da002a28d89285bb,
6ff4b5e21d8dc4062efee6e9b22bbcced4fa50db and
96c3ac979592e6397b759c39abf5364ce6d90f10, respectively. Exact hashes and source
copies are in `evidence/source_manifest.json`.

Focused reading: lines278--354 parse already-recorded banks and dispatch
`Scaler3801` decoding; lines371--374 extend a wrapped clock; lines495--501
interpret a decreasing count as a32-bit wrap. This is **offline replay**
behavior. It does not operate the physical FIFO or establish how online ROC
counts were accumulated. The raw4305 agreement remains stronger evidence
against a ROOT-only origin than merely inspecting this source. A32-bit-wrap
branch in code is not evidence that the observed persistent deficit is a wrap.
No candidate source version is asserted to be the exact pass2 build.

The neighboring modulefile `nps_replay/11.11.23` explicitly sets an offline
replay environment and references hcana/1.0.1 and NPSlib/bf8ec57. Its comment
mentions `/cdaqfs1/apps`. The2.20.24 module references hcana/1.0.2. These are
location/version clues, not evidence which module the run or replay loaded.
`/cdaqfs1`, `/u/cdaqfs1`, `/home/coda`, `/u/home/coda`, `/adaqfs` are absent on
this node. The saved exact errors say No such file or directory. Top-level
inspection of neighboring NPS/HallC directories found offline software,
simulation or older documents, not a retrieved coin_sparse configuration.

### Explicit ROC identity check

Both logs independently list event-builder inputs named `ROC5` and `npsvme5`
(lines25/27 and15/17). These are distinct component names. Do not identify the
five detailed NPS TI dumps with the HMS scaler ROC simply because a name ends
in5. Attributing interleaved register text by nearest log prefix is unsafe.

ROC5's own end messages give:

| Run | ROC5 blocks / events | Source line in frozen user log |
|---|---:|---:|
| 4303 | 3310432 /3310432 | 8025 |
| 4305 | 560647 /560647 | 5173 |

These agree with recorded totals. **A ROC event count is not a per-channel
scaler-integrity check.** A correctly counted event can contain stale or
incompletely accumulated scaler values. This reinforces why neither the ROC
end count nor TI acceptance agreement fixes the scaler deficit. The wording
"five TI slaves" in the prior report must not be read as verification of all
DAQ components or ROC5's hardware configuration.

**Decision:** the apps lead establishes useful replay-source provenance and
an online-filesystem clue, but has not supplied the needed ROC code. Do not
change any correction based on it.

## Outstanding hypotheses and discriminating evidence

| Hypothesis/question | Present evidence | Needed discriminator |
|---|---|---|
| ROOT introduced4305 deficit | All relevant raw counters and physics events agree with ROOT | Ruled out for this defect within the checked channels/file |
| A single cumulative bank was omitted | Missing L1 count persists through later banks | Online storage/accumulation model; simple omission alone insufficient |
| Scaler group inhibited or cleared | Only selected module groups lose exposure | Run-period input/inhibit wiring, control registers and transition/readout code |
| FIFO/online accumulation lost counts | Compatible with persistent deficit; not demonstrated | Actual driver, buffer/overflow handling, accumulator update and sync paths |
| EDTM is representative total-physics probe | Timing establishes pulse association only | Injection point, trigger overlap/causality, prescale and loss-stage coverage |
| PRE or unassigned channels can replace D | Counts correlate but mapping unknown | Dated wiring/configuration and independent channel identification |

No hypothesis is accepted solely because it produces a smoother curve or
agrees with CLT_physics. The normalized yield is proportional to1/L when other
factors are fixed, so a run-dependent error in L creates the opposite change
in normalization. The previous same-cache sensitivity of-0.76662% to+0.69107%
is an illustration with fixed yield/charge, not a measured cross-section bias.
Beam-current, timing and denominator exposure choices must remain visible.

Next useful input: an accessible copy/path for the run-period Hall C
`coin_sparse` configuration, ROC5 rol1/rol2 libraries and sources, scaler
driver, and master/TS setup. Preserve revision, date and provenance when found;
inspect its initialization, readout, sync and end behavior before designing a
repair. EDTM wiring remains a separate requirement even after ROC diagnosis.

## 2026-09-10: focused NPSlib review

User requested `/group/nps/apps/NPSlib` specifically. It is the same filesystem
directory as `/u/group/nps/apps/NPSlib`; this review inspected source/history
beyond the earlier online-driver keyword search. [Detailed report and evidence](npslib_roc_research_20260910/REPORT.md)
contain20 hashed source snapshots, commands, clean source status and a model
of configuration-key handling using the already archived raw text.

**New clue:** `THcNPSCoinTime.h:78` in all three installed trees describes T1
as `NPS VTP||EDTM`, T3 as `HMS 3/4`. Blame traces this comment to2023-09-12,
commit34c3ae03. Timing implementation uses configured trigger-detector entries
0,1,3, described as pTRIG1_ROC1,pTRIG3_ROC1,pTRIG4_ROC1. These comments support
a routing hypothesis; they are not proof of run-period wiring or a map from
timing labels to TS/scaler inputs. If injection is only after VTP formation,
EDTM would bypass earlier VTP formation losses, but that condition is unverified.

**Config preservation limitation:** the event137 handler, added2024-01-11,
parses text into a unique-key std::map using emplace. Within one instance it
keeps the first value per key, without crate/index qualification. Archived4305
text has32 VTP_NPS_TRIG_PRESCALE rows per VTP, and five VTP latency values
2560,2556,2556,2536,2536. A first-key model retains only `0 1` and `2560`.
This is a source-semantics model, not an execution of pass2 or proof which
handler build/lifecycle was used. Original text and all repeated rows remain
preserved. Do not infer a full configuration from one GetInfo value, or equate
VTP internal prescales with Hall C trigger prescales.

**VTP decode clue:** source labels bits0--5 as singles, cosmic scintillator,
cosmic column, different-crate pair, same-crate pair, VLD. These are not the
global TI bits and do not identify EDTM causality. A decision-time mask0x7F
conflicts with an eleven-bit comment; matching firmware must be checked before
any new VTP timing analysis. Previously used TI timestamps/EDTM TDC values are
different fields. No decoder change or physics conclusion follows from the
comment mismatch alone.

Decision: preserve these as evidence leads; exact ROC5 accumulation/FIFO code,
PRE mapping, master settings and EDTM loss coverage remain unresolved. No new
ROOT/raw job, staging, correction or LaTeX edit. Next targeted check is the
run-period trigger-detector setup/map behind these timing labels plus the actual
online coin_sparse readout sources. The latest installed source HEAD predates
the runs, but its use in the run/pass2 has not been established.

## 2026-09-10: dated logbook routing supplied by user

**Update to the earlier "routing unknown" status:** the user provided the
firmware/output mapping entry [4170042,2023-08-23](https://logbooks.jlab.org/entry/4170042)
and counting-house/ROC1 routing entry [4171556,2023-08-29](https://logbooks.jlab.org/entry/4171556),
and confirmed these should apply to these runs. Adopt that mapping as the
working configuration. Original pages were not retrieved because web returned
401 token_expired; supplied text/dates/URLs and applicability are preserved.
This improves on the source-comment hypothesis without pretending to be an
independent hardware survey. No need to request applicability confirmation again.

[Detailed routing report](nps_trigger_routing_20260910/REPORT.md) includes the
two excerpts, provenance, timing schematic, all43 run labels, exact input
hashes and reproducible metadata checks. The new join skipped aggregate
run-range rows and asserts exactly one matching row per production run.

Reported ROC1 inputs:1=NPS one-cluster OR two-cluster OR delayed EDTM;
2=NPS cosmic OR LED;3=HMS h3/4;4=HMS EL-REAL;5=1 AND3;6=1 AND4.
Local NPS/cable T6 is the two-cluster output; ROC1 TI6 is the resulting
coincidence with HMS EL-REAL. VTP internal decision bits are a third namespace.
PRE mappings are still not established by the stated pulse widths.

The direct EDTM copy is delayed about3.5us and OR'ed with NPS at the counting
house; HMS pulses also traverse about3.5us loops. NPS gets100ns before the OR
and50ns afterward; at coincidence HMS40ns-wide pulses follow NPS100ns-wide
pulses by about40ns. Singles are further delayed relative to coincidence.
VLD, unlike that direct EDTM copy, is explicitly routed FADC->VTP->V1495->TS.
Firmware capabilities described in the first post are not active settings.

**Physical consequence:** direct EDTM bypasses NPS cluster-trigger formation.
For TI6 it probes an injected coincidence only if the pulser reaches the
selected HMS branch and overlaps correctly. It cannot certify the bypassed
NPS trigger-formation efficiency. For TI4, NPS hardware cluster formation
is not a requirement of the configured input; its reconstruction/readout
efficiencies remain separate questions. Exact HMS discriminator/EL-REAL
fanout and sent-scaler tap still need validation, as do causality and sampling.

Cross-check of the frozen summary with effective PS_N metadata:5 TI6 runs
(4253,4254,4255,4256,4259),2,706,910 whole-file events;38 TI4 runs,75,131,474
events. Every run has a single positive PS_N agreeing with the summary.
TI6 has four p1 and one p5;TI4 has33 p1 and five p2. Do not infer trigger logic
from the metadata coin_status string or equate local and ROC1 numbering.

Decision: future comparisons must retain this physical cohort distinction.
No total-LT number is certified by these excerpts; no normalization changes,
LaTeX changes or new ROOT/raw job. The4303/4305 scaler defect still requires
ROC investigation. CLT_physics is not a baseline, and overlapping or
population-inappropriate correction factors must not be multiplied.

## 2026-09-10: precise meaning of "refinements"

The user correctly pointed out that pE/D with the broad raw EDTM timing window
was already their method. Calling it a new estimator would be misleading.
Re-read the frozen original sources and diagnostic analyze.py to distinguish
implementation differences from validation and interpretation:

1. Original numerator selects events by the current stored in T; denominator
   selects scaler increments by the ending TSH current. The matched diagnostic
   assigns events to `[TSH.evNumber_i,TSH.evNumber_(i+1))` and applies the same
   ending-snapshot current decision to both. The +/-15% band itself is unchanged.
2. Original numerator can include event heads/tails outside the paired scaler
   increments. Matched counts use only intervals fully covered by available
   events, counting the accepted tags and scaler increments over common exposure.
   This is not a new physics/PID selection or completed yield/charge matching.
3. Original denominator accumulates each file separately. The diagnostic joins
   adjacent segments only with exact event continuity and nondecreasing relevant
   counters/clock, retaining eligible cross-file intervals; gaps/restarts split
   components. Identical counter/clock snapshots are handled explicitly.

Same six updated4398files, p1: original E112586,D112746,L0.9985808809.
Restricting numerator to per-file coverage with its old current gives E112263;
aligning to ending-interval current gives E112658. Final matched raw result
E112658,D112746,L0.9992194845, identical to that intermediate value in this
particular run. Thus removal alone is not the final change; net numerator
change is+72. This example does not establish that every run has the same
decomposition or certify the matched result as physical truth.

Formula pE/D, prescale conversion, run-peak definition and broad +/-500raw-channel
window were retained. Timing/phase tests, counter-integrity flags, prescale-aware
statistical comparisons and the TI4/TI6 routing interpretation are validation
and interpretation, not changes to that point-estimator formula. No repair of
4303/4305 was made; the production implementation was not replaced. Its original
independent-Poisson-style error function was inspected but not rewritten here.

Exact4398 numbers and input hash: [comparison record](nps_livetime_refinements_20260910/run4398_comparison.json).
Source references in the frozen full-cache archive: snapshot/
compute_efficiencies_stuff.cxx:530--540, snapshot/prescale_beamtime_helper.h:182--226,
snapshot/newgen_edtm_livetime.h, analyze.py:55--119. Reproduction: inspect these
sources and results/run4398.json; no ROOT/event reread is needed.

## 2026-09-10: complete documentation and LaTeX edition

The user explicitly requested all findings, references, relevant plots and
comparisons, EDTM capabilities/limitations, Hall C treatment, trigger-setting
dependence and updated TeX. This supersedes the earlier hold on presentation
updates. It does not establish a total-livetime prescription or authorize a
production correction change.

New edition: [livetime_complete_20260910/README.md](livetime_complete_20260910/README.md).
Its editable report and Beamer sources compile offline. The report includes
all43 production counts and ratios, nineteen current figures, a complete
reference catalog, failed alternatives, source/version limits and the detailed
reason for separating TI4 singles from TI6 coincidence interpretation. The
slides have a main discussion and technical backup with all-run tables.

No ROOT/raw data were reread and no staging was performed. The new eleven
comparison figure pairs read frozen per-run JSON and interval CSV inputs.
The seven full-cache figure pairs and the applied NewGen PNG are unchanged
copies. Historical168-segment figures and full earlier archives remain separately
labeled. The five pre-update top-level Markdown documents are backed up.

The core method was already the user's broad raw-window pE/D. The nominal
changes are common event/scaler exposure, the same ending-interval current
decision, and compatible segment joins. For4398, same six files:
E112586 ->112263 ->112658 ->112658, D112746 throughout; final raw0.9992194845.
The net+72 tags demonstrate why this is not merely endpoint removal. Narrow
timing remains diagnostic, with pulse-phase-supported rejected tags. CLT_physics
remains only a proposed definition and never supplies a truth baseline.

Working configuration: dated entries4170042 and4171556, user-confirmed applicable.
Five TI6 runs /38 TI4 runs, with matching sole-enabled metadata factors. EDTM
joins the NPS branch at CH downstream of cluster formation. ForTI6 it cannot
measure bypassed NPS formation losses; forTI4 an NPS hardware cluster is not
the selected-input requirement. HMS fanout, sent tap, causal acceptance and
population weighting still need evidence. PRE attribution is still missing.

The documentation explicitly distinguishes conventional electronic/computer/
total LT endpoints, block level from buffer threshold, periodic-pulser from
eligible-event weighting, and single-trigger prescaling from possible
multi-trigger OR/overlap behavior. Actual metadata has one enabled input per
run; mask overlap is not proof of several independent enabled accept paths.

4303/4305 counter defects remain visible. Raw4305 proves the checked deficit
already exists in raw; it does not localize the hardware/ROC mechanism or prove
raw4303 behavior. No repair, missing-pulse addition, yield/charge exclusion,
PRE correction, rate fit, forced unity or total-loss multiplication was adopted.
Same-file old-to-matched-raw fixed-yield/charge normalization sensitivity spans
-0.7666206% (4493) to+0.6910719% (4308), not a measured cross-section bias.

Reproduce the documentation from the edition directory:

```bash
bash build.sh
```

This compiles editable TeX without replacing prose. Optional
`python3 make_comparisons.py` regenerates the eleven new plots from frozen
outputs. `python3 generate_report.py` deliberately replaces generated article
sections from Markdown/CSV and therefore must not be run over unpreserved
direct LaTeX edits. No ROOT environment is needed for document builds.

Build/numeric/link/preservation results are in
[validation.json](livetime_complete_20260910/validation.json). Issued checksums
and file inventory cover the editable source, PDFs and evidence. PDF layouts
were rendered for visual inspection; current table rows and refinements were
checked against the frozen CSV/JSON. Old archive copies are byte-checked against
their originals. Later edits will intentionally change issued checksums.

Remaining scientific work: run-period ROC5 readout/build/configuration, detailed
injection/sent taps, justified pulse/physics weighting, complete yield/charge/
good-helicity matching, and a final uncertainty/total-loss prescription. No
job needs restarting. No additional raw staging is authorized.
