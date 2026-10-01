# KinC_x60_4b LH2: full-cache livetime continuation

Latest extension: [one authorized raw file and both logbook logs](RAW4305_REPORT.md).
The 4305 defect is now confirmed directly in raw scaler banks; all recorded
physics events and relevant scaler values agree with ROOT. Exact EDTM routing
and a final correction remain unestablished. CLT_physics is not a baseline.

Completed 2026-09-09 evening (UTC 2026-09-10). **All 43 production runs,
183 updated-replay segments, 77,838,384 events.** No production gaps relative
to the selected-source tape catalog. All 51 LH2 metadata entries were checked:
46 cached runs / 187 segments / 79,564,296 events including controls and junk.
The ROOT calculation used cache only. Subsequently, the user authorized
staging exactly one raw file: nps_coin_4305.dat.0, request 100232289. That file
and both supplied logs are analyzed in RAW4305_REPORT.md. No production
correction or presentation edit.

## Decision supported by this calculation

**User clarification: CLT_physics is one candidate definition, never the
baseline, truth, calibration target or validation standard.** Neither
agreement with it nor agreement with pN/S establishes the true livetime.
Statistical differences below compare candidate definitions only.

Use identical scaler-covered intervals for event and scaler counts. Retain
the existing broad raw EDTM window as the reference diagnostic; do not promote
the corrected +/-2 ns cut or all-positive-hit counting to a certified EDTM
acceptance. Report pN/S and pulser-subtracted physics ratios as conditional
DAQ-stage checks. The exact total-livetime correction is **not established**:
run-period injection/routing, trigger causality and eligible physics/charge
exposure still matter. No PRE-based correction is justified by these files.

The remaining timing ambiguity is small compared with the original exposure
and prescale effects, but must not be silently erased. Preserve uncertainty,
counter failures and meaningful rate dependence instead of forcing unity.

## What was tried, and what each test established

| Method/test | Finding and limitation |
|---|---|
| Original NewGen on the refreshed files | Reproduces all 11 runs with identical saved file coverage exactly. Other saved comparisons mix changed coverage/replay with method effects. |
| Restrict event heads/tails, use ending-interval current, join compatible segments | Removes event/scaler exposure mismatch without changing the raw timing definition. Decomposition remains in each run JSON. A current cut alone is not full good-helicity/yield/charge matching. |
| Corrected EDTM +/-2 ns | Rejects 264 raw-window tags that all match independent pulse phase within 25 ticks. A narrow peak does not prove acceptance efficiency; genuine pulse association also does not prove which signal caused the trigger. |
| EDTM time minus selected-trigger time | Does not universally restore raw timing tails. Whole-file checks recover only 1/14, 1/3, 0/0, 1/6, 0/4, 3/12 with relative +/-2 ns for 4253,4255,4259,4301,4305,4350. |
| All positive EDTM hits / wider phase windows | Includes nearby pulses consistent with capture in an earlier physics event. Not an independently accepted EDTM count. |
| Uniform pulse-phase sideband subtraction | Unjustified: raw-zero background is strongly depleted immediately around the pulse phase. Trigger/DAQ correlations violate that simple model. |
| EDTM scaler versus scaler clock | No selected inconsistency greater than two pulses, yet this misses the 4303/4305 defects because both counters lose exposure together. |
| Event timestamp versus scaler clock | Identifies about two seconds of missing scaler exposure in each of 4303 and 4305. Does not determine whether hardware inhibition, decoding or another failure caused it. |
| Infer missing scaler counts / exclude bad intervals | Not adopted. A repair needs raw evidence; an exclusion must remove corresponding yield and charge exposure as well. |
| Prescale-corrected pE/D versus unity | Sampling fluctuations can give values above one. 4259 is a modest statistical discrepancy under explicit models, not automatic evidence of bad hardware. |
| PRE40/100/150/200 | Addresses/counts do not identify the upstream trigger. Unused. |

## Exact quantities and exposure

On the same selected intervals, define N = recorded events, E = EDTM-tagged
subset, S = configured trigger-input scaler increments, D = sent-EDTM scaler
increments, A = L1-accept scaler increments, p = configured prescale factor.

```text
L_EDTM = p E / D
L_all  = p N / S
L_phys = p (N - E) / (S - D)
L_L1   = p A / S
L_all  = (1 - D/S) L_phys + (D/S) L_EDTM
```

The last line is an exact shared-count identity, verified numerically for all
43 production runs; the ratios are not independent measurements. Subtract
accepted E from N and sent D from S. Calling L_phys a physics acceptance
requires one sent pulse to contribute to the chosen input, E to identify the
corresponding accepted subset, and representative prescale sampling. Trigger
overlap, accidental tagging and upstream loss require separate validation.

Physical endpoints are distinct. For eligible ideal physics population G,
formed trigger input B, prescale pass P, and recorded event R:

```text
L_elec = Pr(B | G)
L_DAQ  = Pr(R | P,B,G)
Pr(R | G) = L_elec * Pr(P | B,G) * L_DAQ
```

For representative Pr(P|B,G)=1/p, prescale-removed total LT is
L_elec*L_DAQ. An EDTM measures only its injection path. A periodic pulse samples
availability at its pulse times; physics weights availability by eligible
input rate, and yield normalization needs the appropriate charge/exposure
weighting. These measures can differ. Do not multiply overlapping EDTM and
DAQ corrections. Detector threshold/PID/tracking efficiencies need their own
reference populations.

The implementation uses [TSH.evNumber_i, TSH.evNumber_(i+1)) and the ending
snapshot's BCM4A current, with the existing per-file helper's current band
(0.5-uA histogram peak, +/-15%). It joins adjacent files only across compatible
event/counter/clock continuity. Initial cumulative baselines, uncovered heads
and tails, gaps and restarts are not silently included. Duplicate snapshots
are removed only when clock and all primary counters are unchanged. Current,
boundary and timing alternatives remain diagnostic. Prescale setting zero
maps to p=1; positive setting s maps to 2^(s-1)+1 in this sample.

Raw E uses raw>1 and +/-500 raw channels about the run peak (about +/-48.83 ns,
not +/-500 ns). Tight E uses raw>1 and the corrected-time peak +/-2 ns;
on these files every tight candidate also falls within the broad raw window.
Both have identical nominal event/scaler exposure. The tables retain ratios
above one and counter flags without clipping or inferred counter repair.

## Independent pulse-phase result

Pulse spacing is about 6,257,191.7 g.evtime ticks. Fit alternate core events as
local interpolation anchors; validate with the held-out half. No extrapolation
at file edges. Healthy event/scaler clock comparisons are consistent with
about 250 MHz; phase results remain in ticks rather than claiming an absolute
clock calibration. Held-out per-segment 99.9% absolute residuals span
6.452--10.516 ticks. The earlier 20-second straight-line model could not track
oscillator wander and was replaced before classification conclusions.

Production totals on nominal intervals:

| Selection | Count | Within 25 ticks | No interpolation support |
|---|---:|---:|---:|
| Broad raw window | 3,650,941 | 3,650,890 | 51 |
| Corrected +/-2 ns core | 3,650,677 | 3,650,626 | 51 |
| Raw window rejected by tight cut | 264 | 264 | 0 |
| Positive raw hits outside broad window | 9,395 | 153 | 1 |

All 1,825,258 held-out core events with a model pass 25 ticks; 14 are outside
10 ticks. All 264 rejected raw tags pass 25 ticks, 226 pass 10 ticks. This is
strong pulse-association evidence, **not a causal trigger label**. Half of the
full core sample is used as anchors, so use held-out validation to assess the
model, not the forced-zero anchor residuals.

For 5,515 positive hits outside the raw window and within 500 phase ticks,
relative EDTM-trigger time versus event phase has slope -3.9987366 ns/tick and
correlation -0.9997145. This supports nearby pulses captured in physics-event
TDC windows; it does not authorize adding them to accepted E. Among 69,906,458
raw-zero events, none passes 25 ticks, two pass 50 ticks, and 738 pass 125 ticks.
A flat phase-background assumption is therefore demonstrably inadequate.

Changing matched raw E to tight E moves any production ratio by at most
0.0001595202 (0.015952 percentage points, run 4254). Retain this sensitivity;
it is not a validated confidence interval or the complete timing uncertainty.

## Counter exposure defects

52,025 selected intervals across 46 calculated runs have exact recorded event
numbers at both scaler boundaries. Only these two disagree by more than 0.1 s:

| Run | Event interval | Event-clock span (s) | Scaler span (s) | N | A=S | D | Eraw |
|---|---|---:|---:|---:|---:|---:|---:|
| 4303 | [455000,456915) | 2.000258 | 0.000869 | 1915 | 3 | 1 | 79 |
| 4305 | [105080,106127) | 2.001508 | 0.005481 | 1047 | 2 | 0 | 80 |

The clock scale uses a median of healthy positive intervals; a fit across a
counter step could absorb the defect. A merely omitted cumulative snapshot
would normally recover exposure at the next read, so do not assert that a
missing scaler bank is the specific cause. Charge increments are also tiny.
Ratios for these runs retain visible flags and must not be treated as sound
normalization values. Raw-bank/configuration evidence is needed for repair;
matched yield/charge exclusion is an alternative to investigate, not applied.

4350 retains its updated-replay restart and near-2^32/zero-clock evidence.
Affected local anomalies are below the current band; component boundaries
remain split. Real input-to-accept loss intervals with N=A are retained, not
mistaken for the N>A exposure defects.

## Prescaling and residual agreement

4259 has p=5, N=56,200, E=6,632, D=32,181, S=281,210, A=56,200:

```text
EDTM                    1.03042168
all-trigger ratio       0.99925323
subtracted physics      0.99522546
EDTM - physics          0.03519623
```

Under independent representative thinning with true LT=1, the pE/D standard
deviation is sqrt((p-1)/D)=0.011149, making the excess about 2.73 sigma. For a
comparison with physics acceptance, condition on N,S,D and assume pulser
labels exchangeable among input trials:

```text
E ~ Hypergeometric(S,D,N)
mean(E) = N D / S
var(E) = N (D/S) (1-D/S) (S-N)/(S-1)
```

The conditional z is 2.97; paired 20-second block resampling gives difference
z=2.61. The exact two-sided hypergeometric probability is retained in the CSV.
Periodic pulses need not satisfy exchangeability; block errors are not a
systematic uncertainty, and multiple run comparisons weaken an isolated
outlier claim. Other prescaled runs' paired residuals are below one sigma.

Every unprescaled run's raw-minus-physics residual is below 1.7 paired block
sigma, but **that does not establish absence of a collective offset**.
For the 35 counter-clean unprescaled runs, an inverse-variance common-mean
diagnostic gives raw-minus-physics (4.3734 +/- 1.4349)e-5 (3.05 sigma), while
tight-minus-physics gives (-3.0850 +/- 1.5160)e-5 (-2.03 sigma). These assume
independent run errors and exclude shared systematics. The sign changes with
timing selection; do not turn per-run agreement into a precision calibration.

## Why this can create artificial normalization structure

4398's saved one-segment NewGen is 1.00160539. Original NewGen recomputed on
all six updated segments is 0.99858088. On the same six files, matched raw
EDTM is 0.99921948, tight EDTM 0.99914853, all-trigger ratio 0.99916773, and
raw-subtracted physics ratio 0.99916563. File and method changes are distinct.

Holding yield and charge fixed, replacing same-cache old LT with matched raw
LT changes normalization by Lold/Lmatched-1: -0.7666206% in 4493 to +0.6910719%
in 4308. This is a sensitivity illustration, not an observed cross-section
bias or a certified replacement correction. Prescale changes, coverage,
current/endpoints and counter failures can otherwise mimic run dependence.

## Files, checks and reproduction

- [Production table](production_summary.csv): Eraw and raw-subtracted physics
  explicitly paired. `run_summary.csv` is the broader legacy audit table;
  its CLT_physics column uses Etight. Do not mix the two definitions.
- [Raw/tight statistical tests](acceptance_tests.csv),
  [combined residual diagnostic](combined_residual_diagnostic.json),
  [phase totals](phase_summary.json), [event-clock checks](event_clock_checks.json).
- [All-run comparison](figures/refresh_livetime_comparison.png),
  [timing/phase relationship](figures/refresh_phase_timing.png),
  [rejected tags](figures/refresh_rejected_tags.png),
  [correlated differences](figures/refresh_correlated_difference.png),
  [counter exposure](figures/refresh_counter_exposure.png),
  [normalization sensitivity](figures/refresh_normalization_sensitivity.png).
  Standalone PDFs accompany these diagnostic figures; the presentation is untouched.
- Original applied NewGen plot and source/CSV snapshots are in `snapshot/`.
  Run intervals, JSONs, phase diagnostics, selection records and logs retained.
- All 187 selected file size/mtimes stable; 187 stored prescale records checked:
  45 populated and agreeing, 142 placeholders. All set/read flags zero; all
  events have the configured trigger bit. These are not independent run-start
  TI/routing proof.
- Independent uproot checks reproduce representative T and TSH scalar values
  exactly. It completed outside the sandbox after the sandboxed read stalled;
  the stalled validation process was then stopped. ROOT exports ran in sandbox.
- 21 production source/snapshot files remain unchanged. Previous numerical
  archive's 402 hashes and Beamer archive's 148 hashes all verify unchanged.
- Reader optimization reproduces 4259 and 4253-seg0 binary columns exactly.
  Faster interval-mask selection reproduces run4253 phase JSON byte-for-byte.

For a **fresh writable copy** of this compact archive, recreate the omitted
column cache from the frozen manifest (requires these same cached files):

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; hcana -l -b -q "export_cached.C(0,0,1)"' > export_batch_reproduction.log 2>&1
python3 analyze.py
python3 timestamp_ready.py
python3 audit.py
python3 acceptance_tests.py
python3 clock_scaler_checks.py
python3 event_clock_checks.py
python3 phase_timing_bridge.py
python3 final_consistency.py
python3 make_refresh_figures.py
```

Do not run inventory.py/run_exports.py for an exact reproduction: those are
live-cache refresh/dispatch tools. Use the archived manifest and catalog.
Audit's size/mtime checks detect changed inputs. `final_consistency.py` also
checks original local source/archive paths; distinguish changed provenance
from a numerical reproduction failure if rerunning later elsewhere.

## Specific remaining evidence

At the initial scan no production raw EVIO files were cached
(`raw_cache_availability.json`). The user subsequently authorized just
`/mss/hallc/c-nps/raw/nps_coin_4305.dat.0`; request 100232289 supplied the file,
which was fully read and verified against both tape checksums.
Do not stage any additional file. User-provided run4303/4305 logbook logs are
preserved in `logbook/`, with attachment provenance and SHA256 hashes.
Need run-period EDTM injection/cabling and trigger logic,
TI buffering/prescale readback, and raw scaler banks/configuration around the
4303/4305 intervals. The ROOT scalar branches do not recover unrecorded events
or uniquely identify a pulse as the cause of a trigger. The supplied logbook
logs do not contain those wiring/master details. Existing source context remains in the earlier Beamer
`RESEARCH.md`; its historical or other-experiment evidence is not run proof.

Next: use that evidence to define the EDTM's physical endpoints and timing
acceptance, resolve/exclude defective exposure consistently with physics
yield and charge, then validate the final correction. Do not update the
presentation or production correction merely because one ratio looks flatter.
