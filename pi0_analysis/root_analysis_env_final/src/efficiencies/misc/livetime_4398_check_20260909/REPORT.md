# Run 4398 livetime investigation: six updated segments

2026-09-09. Source set approved by the user in the resumed chat. The six-segment
file/code investigation is complete. **Total-livetime attribution remains
conditional on hardware evidence**, not a finalized normalization prescription.
Production source code and efficiency CSVs have not been modified.

## Full-run result

All six updated segments were retrieved and inspected. They contain exactly
3222394 contiguous event numbers, 1--3222394. Joining real scaler snapshots
across segments gives 2233 snapshots, spanning clock 21.963170--4160.053707 s.
No reset or rollover occurs in EDTM, TRIG4, L1-accept, clock, or BCM4A charge
over this stitched sequence. Six repeated terminal rows are excluded.

For the **fixed diagnostic current window 33.3625--45.1375 uA**:

| Quantity | Count / value |
|---|---:|
| Covered time | 2821.898170 s |
| TSH BCM4A charge over those intervals | 109930.299545 uC |
| Recorded events N / L1 accepts A | 2893284 / 2893284 |
| Trigger-4 input count S | 2895694 |
| Delivered EDTM count D | 112746 |
| Accepted EDTM, raw window | 112658 |
| Accepted EDTM, corrected time +/-2 ns | 112650 |
| EDTM survival, raw window | 0.9992194845 |
| EDTM survival, +/-2 ns | **0.9991485286** |
| Computer LT, physics, +/-2 ns | **0.9991685076** |

These figures use the stated half-open event-boundary convention. It excludes
one final non-EDTM event from the full data set; its effect is negligible at
the displayed diagnostic precision but boundary/latch timing is not asserted
as exact. Switching to `(b_i,b_(i+1)]` changes the +/-2-ns pulser count from
112650 to 112652: EDTM becomes 0.9991662675 and physics CLT 0.9991677890. This
measured boundary sensitivity is 0.000017739 for EDTM and 0.000000719 for CLT.
The segment-0 unsupported tail is recovered by subsequent real
snapshots in the complete run. Do not discard every file's tail and first
snapshot independently when a later segment supplies the missing coverage.

With a 2-uA current threshold the measured physics computer LT is 0.9992068726,
consistent with the stored database value 0.9992 at its reported precision.
That is consistency of a computer-livetime quantity, not proof of total
electronic/DAQ coverage. The fixed window above is an explicit diagnostic
selection, not an emulation of every per-segment production helicity selection.

The user clarified that PRE channels have **uncertain trigger attribution**.
No PRE-based correction is adopted. Establish their physical inputs and actual
widths before assigning them an electronic-livetime meaning.

See the [detailed presentation](Run4398_livetime_investigation.pdf),
[full-run figures](full_run_diagnostic_figures.pdf),
[full-run counts](full_run_estimators.csv), and
[full-run interval records](full_run_intervals.csv).

## Segment-0 result that identified the endpoint problem

The recorded EDTM ratio above one is reproducible and has an identifiable
endpoint mismatch in the available updated segment 0. The numerator includes
49 EDTM events in an event interval with **no new scaler counts**. Removing
that uncovered interval from the numerator changes the existing ratio from
`19341/19310 = 1.0016053858` to `19292/19310 = 0.9990678405`.

A reference-subtracted EDTM timing cut of +/-2 ns about 245.11508 ns removes two
additional outliers, yielding `19290/19310 = 0.9989642672`. This is a measured
pulser survival ratio over the covered intervals. Whether it is the correct
total livetime for physics still requires signal-path and DAQ-mode evidence.

The PRE40/100/150/200 branches are present **and populated** in both scaler and
TDC data. Missing pretrigger branches are not established for this file.

See [diagnostic figures](segment0_diagnostic_figures.pdf),
[exact counts](matched_estimators.csv), and [machine-readable results](results.json).

## Coverage and provenance

Input inspected:
`/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_4398_0_1_-1.root`.
Only segment 0 was cached at the start. Tape catalog and archived run list both
contain segments 0 through 5. The existing efficiency CSV reports one segment,
so it cannot be described as a complete-run result.

The five missing updated segments were requested with `jcache get`, request
**100230264**, and all completed before the full-run analysis. The older
`production` directory also has all six segments, but those were not substituted
for the requested updated replay. The full-run calculation examines continuity
and cumulative scaler offsets across all updated segments; it does not average
per-file ratios.

Segment-0 contents: `T` has 542584 entries, `TSH` 355, `TSHelH` 19435, and `E` 313.
The event numbers in this segment are exactly 1 through 542584, with no gaps or
duplicates. Every recorded physics event has `g.evtyp=1`; trigger masks vary,
and every mask includes bit 3 (value 8). Event type alone therefore cannot
separate EDTM from physics. Multiple bits do not establish multiple enabled
prescaled triggers.

The archived job `NPS_PASS2_x60_4b_x60_2_nps_coin_4398.dat.0` explicitly names the
updated output and `/group/nps/nps-ana/nps_replay_pass2_github.tar.gz` as its replay
input. That tarball contains git HEAD
`2d462a85b860c54785d8287e644619a1cb1acc19`, with relevant source snapshots in
`replay_snapshot/`. Its working-tree contents, rather than the git identifier
alone, are the archived configuration evidence. The job's original stdout and
stderr paths are no longer present at the recorded locations.

`Run_Data` records February 13, 2024, 17:03:18 and stored prescale factors
`-1/-1/-1/1/-1/-1` for triggers 1--6. **Its prescale set/read flags are both zero**:
these stored parameters are not independent proof of a decoded DAQ prescale
event. The master CSV also says `ps4=0`, interpreted by the analysis as factor 1.
The nearly equal accepted/trigger counts support this interpretation, but the
run-start TI configuration is still desirable.

Ten embedded DAQ configuration strings were extracted to `daq_config.txt`.
They contain FADC/VTP settings, including VTP firmware 65536 and trigger width 20.
No TI/TS buffer-level, block-level, ROC-lock or busy configuration was found in
those strings. VTP prescales must not be confused with TI prescales.

The installed hcana source inspected for decoder behavior is commit
`2f93d89f1ef482f36666a4e91821493ce4c22fc5`; it is supporting implementation
evidence, not a recovered hash of the executable that produced these files.
`code_provenance.json` records the source paths and SHA-256 hashes actually read.

## What each counter measures

The archived `MAPS_db_HScalevt.dat` maps the following channels. Hardware labels
are map evidence; historical wiring still needs run-period confirmation.

| Quantity | Branch | Archived scaler mapping / interpretation |
|---|---|---|
| Sent pulser count, D | `H.EDTM.scaler` | ROC 5, slot 12, channel 14; labeled HMS EDTM |
| Copy of pulser count | `H.EDTM_CP.scaler` | ROC 5, slot 10, channel 6; not an accepted-pulser monitor |
| Trigger-4 input count, S | `H.hTRIG4.scaler` | ROC 5, slot 11, channel 13 |
| L1 accepted count, A | `H.hL1ACCP.scaler` | ROC 5, slot 10, channel 16 |
| PRE100 / PRE150 | `H.hPRE100.scaler`, `H.hPRE150.scaler` | ROC 5, slot 10, channels 9 / 10 |
| Recorded event count, N | `g.evnum` / rows of `T` | Count only events belonging to covered scaler intervals |
| Accepted pulser count, E | `T.hms.hEDTM_tdcTimeRaw`, `...tdcTime` | Event-level timing classification |
| Common event boundary | `TSH.evNumber`, `TSHelH.evNumber`, `T.g.evnum` | Must distinguish event number from each tree's record counter |

In this segment EDTM and EDTM_CP counters agree, while the two 1 MHz clocks
agree to a few ticks. These are duplicate monitors, not independent measures
of busy-gated versus ungated livetime.

Conceptually, the relevant path is:

```text
detector/analog response -> trigger-forming electronics -> TRIG4 scaler
                                                    -> TI prescale/busy -> L1 scaler -> T
EDTM injection -> only the stages downstream of its actual injection points
NPS VTP / waveform acceptance -> coverage must be established separately
```

The archived September 2023 NPS run plan, printed page 5, identifies trigger 4
as HMS EL-REAL and trigger 6 as NPS-trigger AND HMS EL-REAL. This is historical
context, not run-4398 wiring proof. A TRIG4-based computer-livetime ratio does
not, by itself, certify NPS clustering/readout efficiency or analog losses
upstream of the trigger scaler.

## Endpoint and current alignment

The first real TSH record is at event 1, clock 21.963170 s. It already contains
878 EDTM pulses, 23020 trigger-4 inputs, and one L1 accept. These counts include
time before the first recorded physics event. Blindly including the first
cumulative count as a denominator for events starting at event 1 would create
a different mismatch. Consecutive differences remove this baseline.

The last real snapshot has event boundary 541321 and clock 674.895349 s. A final
row has event boundary 542585 but exactly the same scaler values and time.
Thus the half-open interval `[541321,542585)` contains 1264 recorded events,
including 49 EDTM candidates, with zero new trigger, EDTM, or clock counts.

The inspected `THcScalerEvtHandler::End()` increments the event number and fills
a final scaler row even after processing zero delayed scaler updates. This is
consistent with the observed duplicated endpoint. Reject intervals with no
positive clock increment for this diagnostic; later segments may supply real
snapshots spanning a segment boundary, so stitch before discarding coverage.

For each consecutive real scaler boundary pair `[b_i,b_(i+1))`, count events
in that same interval, use the difference of cumulative counters, and apply
the current reported at its ending scaler snapshot. The event-stored clock
matches the **preceding** TSH snapshot for every event. Its current therefore
describes an earlier interval, not necessarily the interval containing that
event. This matters at beam trips and ramps.

For the production current window 33.3625--45.1375 uA, the matched covered
duration is 483.315153 s. Recorded events and L1 scaler counts both total
495221. Individual intervals differ by up to two events near latch boundaries;
their sum agrees exactly. Changing inclusive/exclusive endpoints provides a
one-event boundary sensitivity, not a substitute for hardware latch timing.

The cumulative-count diagnostic shows the endpoint jump directly:

![Endpoint mismatch](endpoint_mismatch.png)

## Measured ratios and assumptions

These formulas assume the expected single enabled trigger and factor 1, and
that the pulser contributes one count to the trigger-input denominator.
Let N be all recorded trigger events and E the accepted EDTM subset. Then:

```text
pulser survival       = E / D
computer LT, all      = N / S
computer LT, physics  = (N - E) / (S - D)
```

These are live fractions; an inverse correction would be `1/L`. If EDTM covers
both electronics and computer losses, multiplying its ratio by computer
livetime again double counts the covered computer loss.

| Segment-0 covered selection | N | E, raw window | D | S | EDTM raw | Computer LT, physics, raw |
|---|---:|---:|---:|---:|---:|---:|
| All currents | 541320 | 26070 | 26087 | 541790 | 0.99934833 | 0.99912159 |
| Current >2 uA | 538363 | 23270 | 23286 | 538833 | 0.99931289 | 0.99911938 |
| Production current | 495221 | 19292 | 19310 | 495679 | 0.99906784 | 0.99907635 |

With the +/-2 ns corrected-time cut, production-current E is 19290, pulser
survival is 0.99896427 and physics computer LT is 0.99908054. The latter changes
very little because EDTM is a small part of the trigger population.

Subtracting **sent** pulses from the accepted L1 numerator gives
`(A-D)/(S-D)=0.99903856`, which differs from subtracting the measured accepted
pulses. The database label saying EDTM is subtracted does not establish which
quantity was subtracted in its original producer.

## Timing, uncertainty, and PRE monitors

`tdcTimeRaw` stores raw channels. The archived conversion is 0.09766 ns/channel;
the current code's +/-500 raw-channel window is +/-48.83 ns, although its
comment calls it 500 ns. The corrected-time pulse peak is near 245.11508 ns.
The measured numerator is unchanged for corrected half-widths 2, 3, 5 and 10 ns.
Two raw-window candidates lie farther away. They are timing outliers; their
origin is not established as an accidental background by this test alone.
All raw-window EDTM candidates in the covered sample have multiplicity one.

The measured pulser rate in the production intervals is about 39.95 Hz.
The conditional independent-Bernoulli standard-error scale is 0.000220 for the
raw-window pulser fraction and 0.0000432 for all-trigger computer livetime.
These are **not final uncertainties**: periodic sampling, time correlations,
boundary assignment, and signal-path coverage have not been certified. The
current implementation's independent-Poisson propagation of numerator and
denominator gives 0.01018935 for the original ratio; it omits their covariance
and cannot be justified merely by calling this a ratio of counts.

Tightening timing changes pulser survival by 0.00010357. Restricting to intervals
whose previous current also passes the production window changes the raw
ratio to 0.99906005. The lower-current comparison in the table is diagnostic,
not proof of a rate law. No unbuffered 250-us correction from a different
experiment has been applied; actual buffering and readout conditions remain
unresolved.

For the matched production-current intervals, HMS PRE100 and PRE150 contain
612581944 and 606870963 counts, respectively: ratio 0.99067720. PRE200/PRE100 is
0.94158061; PRE40 is essentially equal to PRE100. The corresponding pPRE150 /
pPRE100 ratio is 0.95547866. None is interchangeable with the measured EDTM
ratio. A branch label or historical width formula alone does not establish
which trigger population these channels sample in run 4398, their true widths,
or whether a linear deadtime approximation applies. Do not choose an arm or
formula just because it produces a plausible correction.

## Existing implementation findings

1. `newgen_edtm_livetime.h` computes `E*ps_factor/D`. The mismatch originates in
   event/scaler accumulation and coverage, not that division alone.
2. `compute_efficiencies_stuff.cxx` accepts event-current EDTM candidates without
   a check for a subsequent real scaler snapshot. This reproduces the 49-pulse
   uncovered tail. Its timing-unit comment is also incorrect.
3. `good_event_selection_helper.h` builds `accepted_evcount_ranges` from
   **TSHelH** records. `prescale_beamtime_helper.h` applies them to **TSH.evcount**.
   The two counters span 1--19435 versus 1--354 and do not identify the same
   intervals. This explains why the saved diagnostic denominator is only 6550
   versus 19310 and its reported beam time is 163.950216 s versus 483.315153 s
   for covered current intervals. Use common event/time boundaries when
   connecting these trees. Do not promote `18970/6550` to a livetime.
4. `root_file_discovery.h` chooses the updated directory as soon as it finds
   any matching file, without an expected-segment completeness check. This
   permits the one-segment result to be marked processed.
5. The active older `nps_livetimes.h` reads database `CPU_LT`; its old calculator
   is commented out. The database gives run 4398 `CPU_LT=0.9992`, six segments,
   and a 2-uA cut. It is not an independent segment-0 total-livetime measurement.
   The master metadata CSV's `Computer_live_time=0.992` has different unresolved
   provenance and should not be substituted silently.
6. `compute_luminosity_scaler.cxx` counts TDC events over the whole available
   event stream, then applies beam-on fractions. Its EDTM expression simplifies
   to all-event accepted EDTM divided by all-interval scaler EDTM, rather than
   a directly matched beam-on estimator. It shares endpoint exposure.

These are documented findings, not edits to the frozen efficiency calculation.

## Full-run timing and time-variation checks

At fixed production current, the corrected-time numerator is 112650 for +/-2 ns,
112651 for +/-3 ns, and 112652 for both +/-5 and +/-10 ns. The latter ratio is
0.9991662675. Timing outliers are observable; they have not been conclusively
classified as accidentals. The raw-window result is higher, 0.9992194845.

The conditional independent-binomial scales are 0.0000869 for EDTM and
0.0000173 for physics CLT. Time-block resampling with 10/20/60-s groups and
2000 replicates gives spreads of about 0.00013 for EDTM and 0.00011 for CLT.
This is a diagnostic of time variation and correlations, **not an adopted
systematic error or a missing-hardware-coverage uncertainty**. See
`full_run_block_bootstrap.json`; the seed is 4398.

Losses are not uniform in time. In the ~2-s interval ending at clock 820.125614 s,
current is 39.220992 uA, trigger inputs are 2082, accepted events are 1799, and
accepted/delivered EDTM is 69/80. This is 283 lost trigger inputs and 11 lost
EDTM pulses in the same interval. Such real loss bursts are retained in the
average. Their correlation supports the pulser's sensitivity to measured DAQ
losses, but cannot certify upstream electronics or NPS coverage.

## Remaining requirements before a final total-livetime prescription

- The six-segment staging, continuity, primary-counter reset and endpoint checks
  are complete. Confirm the exact physical latch convention before production
  integration, and preserve the run-level joining of snapshots.
- Reconcile the chosen event sample and charge integration with the same good
  beam/helicity intervals. Excluding a tail from the livetime study does not
  authorize keeping unmatched physics yield or charge in the production result.
- Recover run-4398 TI prescale/readout/buffer configuration and actual EDTM / PRE
  signal routing, including whether the desired NPS path is sampled. The user
  confirms the Hall C logbook requires sign-in; no login was attempted.
- Establish the pass2 executable/module version and any deviations from the
  archived replay source; confirm scaler latch timing and trigger-mask meaning.
- Only then select a total physics correction and finalize its systematic/
  periodic-pulser treatment. The presentation supplied here documents measured
  results and unresolved attribution. If wiring
  evidence remains absent, report identifiable component ratios and the
  unresolved coverage rather than inventing total livetime.

## Reproduce

Run from this directory. Python requires uproot, NumPy and Matplotlib.
The scripts read the explicitly named updated segment 0 and write diagnostic
files locally; they do not run the production pipeline. Generated NPZ and binary
column caches were left under `/tmp/nps4398_resume`, not installed in this archive.

```bash
python3 inventory.py
python3 analyze.py
MPLCONFIGDIR=/tmp/nps4398-mpl python3 results.py
```

`inventory.py` records all keys/branch types and caches only scalar diagnostic
columns. `results.py` generates the matched count tables and four-page figure
PDF. The separate `analyze.py` deliberately includes naive ratios to demonstrate
their failure; use `matched_estimators.csv` for covered-interval measurements.

Commands actually used for independent ROOT verification (the macros retain
their `/tmp/nps4398_resume` output paths):

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; hcana -l -b -q /tmp/nps4398_resume/inspect_run.C'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; root -l -b -q /tmp/nps4398_resume/export_columns.C'
```

ROOT and uproot agreed exactly on all 542584 rows of the 92 exported event
columns and all 355 rows of the 696 scaler columns. See
`reader_verification.json`. Repeating that optional cross-check elsewhere
requires updating the macros' output directory or placing them in the shown
`/tmp` directory. No full production job was run.

The full-run diagnostic commands actually used were:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; hcana -l -b -q /tmp/nps4398_resume/export_full_run.C' > /tmp/nps4398_resume/full_export.log
MPLCONFIGDIR=/tmp/nps4398-mpl python3 /tmp/nps4398_resume/full_results.py
MPLCONFIGDIR=/tmp/nps4398-mpl python3 /tmp/nps4398_resume/make_presentation.py
```

To repeat outside `/tmp/nps4398_resume`, first set the output directory in the
ROOT macros to this archive directory. Run the Python scripts from here;
`full_results.py` uses the sibling `full_run/` exports and verifies that all
18 tree exports completed. The binary column caches were omitted from the
installed archive, so recreate them before rerunning full-run calculations.

```bash
env -u DEBUG jcache status 100230264
sha256sum -c SHA256SUMS
```

`DEBUG=release` inherited in this environment interferes with the Java cache
client; unsetting it for this command fixes the launcher. In this session the
ROOT environment and Python file-data reads needed execution outside the
sandbox; the sandbox module could not resolve the user ID and data reads stalled.

## References used

- Source catalog: `../SOURCES_LIVETIME_4398.md`; original PDFs/text are in the
  adjacent `livetime_4398_sources_20260909/` directory.
- D00, Avnish's May 2025 luminosity presentation: prior findings and motivation.
- D01, Pooser, Live Time Calculations, 1022-v1: slides 7--10 and 14--15 for
  decoding, timing, interval selection and component definitions.
- D02, Mack, 1001-v2, September 20, 2019: slides 2 and 4--7; conditional
  unbuffered periodic-pulser correction, not adopted here.
- D03, Hall C Trigger Electronics, 1028-v5, `trigger_v3.pdf`, printed page 14
  and following: historical trigger/scaler/EDTM topology. Page 14 was checked
  visually against the PDF; NPS run-period applicability remains conditional.
- Archived NPS Run Plan, September 11, 2023, printed page 5: planned trigger
  numbering. D07, Crafts, May 2025, slide 13: a legacy PRE-width convention,
  not validation for this run.
- Local replay/code snapshots, `code_provenance.json`, `run_metadata.log`,
  `daq_config.txt`, `inventory.json`, and original run-specific CSV extracts.
- Supplemental public decoder documentation:
  <https://hallcweb.jlab.org/hcana/docs/THcScalerEvtHandler_8h_source.html>.
  Local source/measurements support the endpoint finding; web retrieval of
  additional raw-source/PDF pages failed and is not treated as evidence.
