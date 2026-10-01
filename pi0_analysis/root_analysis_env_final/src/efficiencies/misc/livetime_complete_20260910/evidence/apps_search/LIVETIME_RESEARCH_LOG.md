# KinC_x60_4b LH2 livetime: research log and audit trail

Started UTC 2026-09-10 at the user's request for a detailed, auditable research
document. This is a living record; the linked numerical/source archives remain
frozen. Entries below distinguish results recovered from previous reports from
checks executed in this continuation. **The exact total-livetime correction is
still unresolved. No production correction has been changed.**

## Reading and maintaining this record

The [resume checkpoint](LIVETIME_REFRESH_CHECKPOINT.md) is a compact handoff.
This document records questions, methods, evidence, alternatives and decisions.
For every continuation, add the date, exact inputs/source revisions, commands
or executable script, outputs, interpretation and remaining uncertainty.
When a conclusion changes, add an explicit correction rather than silently
rewriting the old evidence. Distinguish source-code behavior from behavior
verified in the run. Filenames, directory dates and current module defaults
do not prove the software/configuration loaded during a particular run.

The investigation has four existing evidence packages:

| Package | Role and present status |
|---|---|
| [Initial numerical report](livetime_KinC_x60_4b_lh2_20260909/REPORT.md) | Earlier cache coverage; retained for comparison, not current all-run inventory |
| [Presentation research](livetime_beamer_20260909/RESEARCH.md) | Earlier source research and presentation; frozen pending a justified method |
| [Full-cache report](livetime_refresh_20260909/REPORT.md) and [raw report](livetime_refresh_20260909/RAW4305_REPORT.md) | Current all-production calculation, timing/exposure tests, raw4305 and user logs |
| [CODA follow-up](coda_livetime_research_20260910/REPORT.md) | Version-qualified documentation, register checks, timer limitations |

## Non-negotiable analysis constraints

- **CLT_physics is a proposed definition, never truth, a baseline, calibration
  target or validation standard.** Agreement between candidates is only a
  consistency test. Historical output labels do not override this restriction.
- EDTM routing and PRE40/100/150/200 trigger attribution are unknown. Do not
  infer wiring from channel names, rate correlations or near-identical counts.
- User authorized one raw file; the allowance was used for
  `/mss/hallc/c-nps/raw/nps_coin_4305.dat.0`, request100232289. No further
  staging/resubmission. Cached ROOT production coverage has been processed.
- No LaTeX revision until the physical method is established. No numerical
  repair, event exclusion or production code/CSV change has been adopted.
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
