# KinC_x60_4b LH2: livetime, coverage and artificial structure

2026-09-09. **Diagnostic measurements, not a certified total-livetime correction.**
Production efficiency code, applied efficiency CSVs and physics yields were not
modified. PRE monitors were not used as a correction. No jcache command or tape
retrieval was performed in this extension.

## Results and coverage

Completed **46 cached runs**, including all 43 production runs, from **168 segments and 70,155,126 recorded events**. Five junk runs have no selected cache files: 4307, 4486, 4552, 4553, 4554. The same-coverage control reproduces all 42 saved production count pairs and peaks exactly. 2 selected intervals in runs 4303, 4305 have recorded-event increments exceeding L1 increments by more than four. These defects remain visibly flagged; no ratio was forced to one.

The master CSV selects **51 LH2 entries** by stripped `Kin_old=KinC_x60_4b`
and case-insensitive `target=lh2`: 43 production, six junk, two efficiency.
All 51 have a status in [run_summary.csv](run_summary.csv). Main comparison plots
contain the 43 production runs; efficiency/junk entries remain labeled in the
tables. A missing run is not assigned zero or unit livetime.

The frozen cache inventory has 253 files across both replay variants. The
existing updated-first, production-fallback policy selects **168 files** for
46 runs; the other 85 are alternate replay versions, not additional independent
events. Never combine these versions within a run. Several partial updated
replays have a fuller production alternative: that is recorded, not silently
substituted. `inventory.json` contains every candidate and the selected paths,
size and modification time. `catalog_coverage.json` uses read-only `/mss`
directory listings to identify absent segments; it reads no payload or staging
service. Catalog completeness means file availability, not usable exposure.

| Run | Source | Cached | Missing catalogued segments |
| --- | --- | --- | --- |
| 4253 | updated | 0,2 | 1 |
| 4300 | updated | 0,1,2,3,5 | 4 |
| 4351 | updated | 0,1,2,3,4,5,7 | 6 |
| 4352 | updated | 7 | 0,1,2,3,4,5,6 |
| 4354 | updated | 1 | 0 |
| 4356 | updated | 1,3,4,5,6 | 0,2 |
| 4483 | updated | 1 | 0 |
| 4488 | updated | 4 | 0,1,2,3,5 |
| 4490 | updated | 1 | 0 |
| 4552 | updated | none | 0 |

The supplied [original multipanel PNG](snapshot/efficiency_multipanel_vs_run_KinC_x60_4b_lh2.png)
is copied unchanged. Numeric saved values come from the accompanying snapshotted
`efficiency_KinC_x60_4b.csv`, not pixel extraction. Recomputing NewGen on current
cache files is the control for changing the method. All 42 production runs with
the same saved file coverage reproduce **E, D and the raw peak exactly**. The
43rd, run 4398, now has six segments rather than one. Canonical cache aliases
(`/cache` versus `/lustre24/expphy/cache`) are compared by source and filename.
File identity is additionally checked by frozen size/mtime; full multi-terabyte
ROOT payload hashes were not computed.

## What was tried previously, and what each comparison establishes

The original May 2025 presentation is archived in the sibling
`livetime_4398_sources_20260909/singh_nps_luminosity_202505.pdf`.
Its slide 23 states four definitions:

| Previous method | Actual population / implementation | Limitation established by this investigation |
|---|---|---|
| Total EDTM | Whole-file accepted raw-TDC EDTM / accumulated EDTM scalers; no current cut in that slide | Whole-file events need not have matching initial/final scaler coverage; periodic pulses need representative sampling |
| Computer LT all, TSH | Current-cut L1 accepts / trigger-input scalers | Measures the interval between those taps; not upstream electronic/NPS coverage, and scaler integrity must be checked |
| Computer LT all, TDC | Whole-file trigger-TDC count multiplied by scaler beam-on accept fraction, divided by beam-on trigger inputs | A scaler beam-on fraction is a proxy for event selection; cannot establish direct event/scaler matching or correct missing intervals |
| Computer LT physics, TDC | Trigger-TDC events without EDTM, beam-on-fraction weighted, divided by trigger inputs minus sent EDTM | Same proxy/coverage problem; accepted EDTM must be removed from accepted events and sent EDTM from input counts; prescale sampling assumptions matter |

The historical plots concern luminosity scans, not a rerun of this kinematic
setting. They reported EDTM excesses and current trends, as well as DAQ/logbook
issues in runs 1526--1534. Those historical hardware explanations are **not
established for KinC_x60_4b**. A flat carbon yield is not a criterion for choosing
a correction. The inspected `compute_luminosity_scaler.cxx` applies whole-file
TDC counts and beam-on fractions; its EDTM expression algebraically reduces to
the all-event/all-interval ratio. Source evidence is preserved in the earlier
run-4398 archive, not reinterpreted as a direct interval calculation.

The active older `nps_livetimes.h` reads a database CPU_LT; its old calculator is
commented out. Run 4398's database CPU_LT=0.9992 agrees at stored precision with
the earlier matched >2-uA physics CLT=0.9992068726. This checks a component and
selection, not total loss coverage. The master CSV's separate
`Computer_live_time=0.992` has different unresolved provenance.

**Current NewGen:** the current helper computes `p*E_old/D_old`. For each valid
segment, the selection helper infers a TSHelH current peak using 0.5-uA bins
and sets a +/-15% window. The event numerator uses the current stored on `T`,
requires finite `hEDTM_tdcTimeRaw>1`, and counts within +/-500 **raw channels**
of the run-level 10-channel histogram maximum. The denominator sums consecutive
TSH EDTM differences passing the ending-row current window, separately per
file. The divisor itself is not the endpoint bug.

Current `accepted_evcount_ranges` come from **TSHelH** but the beam-time/gated
diagnostic code applies them to **TSH.evcount**. These are different record
counters. In 4398 segment 0 the saved gated denominator was 6550 versus current-
only 19310. This is a separate diagnostic-gate error: the main NewGen denominator
does not use that gate. The new diagnostic remains current-only too; it does
not claim full good-helicity/yield/charge matching.

## New calculation: an explicit common coverage

1. Freeze the run list, cache paths, source code, selection CSV and applied PNG.
   Source the Hall C/NPS environment. Export only scalar branches needed here.
   Run the snapshotted current selection helper for each segment; failed
   selections and incomplete exports cannot silently enter a result.
2. Check event numbers within each file. Join files only when segment numbers
   and event numbers are contiguous **and** cumulative clock, trigger, EDTM,
   L1 and charge values do not decrease across their boundary. A clock or
   counter restart starts a new component; never bridge missing coverage.
3. Remove a duplicate snapshot only if its clock **and all primary counters**
   are unchanged. A zero-clock record with changing counters is an anomaly,
   not an ordinary duplicate; retain it in the audit. Internal decreasing
   counters would stop that run for investigation rather than be clamped.
4. For boundaries `b_i=TSH.evNumber`, count recorded events in
   `[b_i,b_(i+1))` and difference the corresponding cumulative counters. The
   entire interval must lie inside the component's available event range.
   Initial cumulative counts are a baseline, not a denominator for unrecorded
   pre-Go events. Uncovered heads/tails are reported. A later compatible file
   can recover an earlier segment tail.
5. Apply the ending snapshot's current window to that same interval. At a file
   boundary use the ending file's selection window. `T`-stored scaler current
   normally belongs to the preceding snapshot, hence to an earlier interval.
   Its lag matters at beam ramps/trips. Current is an interval average, not
   an event-by-event instantaneous measurement.
6. Find the corrected EDTM peak for each run from a 0.1-ns maximum and median
   within 1 ns. Require raw TDC>1 and corrected time within +/-2 ns of that
   peak. Retain raw-window and +/-0.5/1/3/5/10/20/50-ns results as diagnostics.
   The historical conversion 0.09766 ns/channel makes +/-500 raw channels
   approximately +/-48.83 ns; the production comment's "500 ns" is incorrect.
7. Sum counts over selected intervals, then divide. Do not average individual
   interval ratios, clip ratios to one, or remove real loss bursts to flatten
   a run. Audit recorded-versus-L1 closure and abnormal scaler jumps separately.

The physical latch convention remains unresolved at single-event precision.
The alternative `(b_i,b_(i+1)]` is evaluated with the same scalers and cuts.
Neither convention is certified by the field name alone. Timing, boundary and
current variants appear in the presentation and all-run CSV.

## Formulas and what they do not establish

Let N be recorded events, E their accepted EDTM subset, S the configured trigger
input scaler, D sent EDTM pulses, A L1 accepted scaler increments, and p the
intentional prescale factor. The diagnostic ratios are:

```text
L_EDTM        = p E / D
L_all         = p N / S
L_physics     = p (N - E) / (S - D)
L_L1-input    = p A / S
```

For each run, select TRIG3, TRIG4 or TRIG6 from the sole enabled metadata token.
The current convention is p=1 for setting 0; otherwise
`p=2^(setting-1)+1`. Thus `ps4=1` means 2, and `ps6=3` means 5.
These ratios remove intentional prescaling *in expectation*. Prescale weighting
of actual physics yields remains a separate production operation.

The physics subtraction assumes one sent EDTM pulse is represented in the
chosen trigger-input population, E classifies the corresponding accepted
subset, and the same prescale applies representatively. Subtract accepted E
from accepted N; subtract sent D from input S. At p=1, `(A-D)/(S-D)` is only
equivalent to accepted-pulser subtraction when E=D (and N=A). Branch names,
trigger masks and stored prescales support configuration checks but cannot
independently prove all of this routing.

The three ratios share counts and obey the algebraic identity
`L_all = (1-D/S)*L_physics + (D/S)*L_EDTM`. Agreement is not independent
validation. A prescaled EDTM estimate above one can be a sampling fluctuation
or bias, rather than literally more accepted than sent pulses: E/D can remain
below one while pE/D exceeds one. Periodic pulser/prescale correlation is a
possibility, not an established explanation for the residuals.

| Run | Trigger / p | pE/D | pN/S | p(N-E)/(S-D) |
| --- | --- | --- | --- | --- |
| 4259 | 6 / 5 | 1.03042168 | 0.99925323 | 0.99522546 |
| 4308 | 4 / 2 | 0.99538264 | 0.99689449 | 0.99695611 |
| 4356 | 4 / 2 | 0.99453103 | 0.99610108 | 0.99616445 |
| 4404 | 4 / 2 | 0.99304482 | 0.99624596 | 0.99637490 |
| 4493 | 4 / 2 | 0.99282550 | 0.99797624 | 0.99818275 |
| 4558 | 4 / 2 | 0.99887773 | 0.99722025 | 0.99715272 |

## Counter defects can survive endpoint matching

| Run / component / interval | dt (s) | N | A | S | D | E |
| --- | --- | --- | --- | --- | --- | --- |
| 4303 / 0 / 314 | 0.000869 | 1915 | 3 | 3 | 1 | 79 |
| 4305 / 0 / 124 | 0.005481 | 1047 | 2 | 2 | 0 | 80 |

| Run | Matched EDTM | If suspect intervals excluded | Matched physics CLT | If excluded |
| --- | --- | --- | --- | --- |
| 4303 | 0.99978196 | 0.99915204 | 0.99975006 | 0.99914534 |
| 4305 | 1.00154745 | 0.99955074 | 1.00159674 | 0.99965856 |

The supplementary `closure_exclusion_*` columns remove selected intervals with
**N-A>4** as a sensitivity study. Four is a diagnostic tolerance above the
ordinary one/two-event latch differences, not a calibrated hardware cut.
This variant is **not adopted** as the corrected answer. Negative N-A can
represent loss after L1 and is retained. In particular, the criterion does
not remove S-N losses when N agrees with A. The removed interval counts and
charge are explicitly tabulated. Removing these events only from livetime
would still leave a biased production result unless physics yield and charge
were treated consistently. A repair requires raw/replay evidence, not the
numerical attraction of a result near one.

Run 4350 includes a zero-clock record with D=35540 at current 0, followed by an
EDTM increase of 4294931836 at approximately 1 uA. The latter produces a
cumulative value near 2^32; this is an observed replay-counter anomaly, not a
measured physical pulser burst. Both fall outside the selected beam window.
The 3->4 segment boundary resets the EDTM cumulative offset; 4->5 restarts
the scaler clock. These boundaries are not bridged. All raw values and
selection membership remain in the audit.

Real busy periods remain in the results. The prior full 4398 check found
S=2082, N=A=1799, D=80, E=69 in the interval ending at 820.125614 s. This
coincident loss of trigger and pulser accepts is different from N>A with
almost no elapsed scaler time. Time variation is not automatically an artifact.

## Comparison plots and why normalization can acquire false structure

![All production runs](figures/comparison_all.png)
![Factor-one production runs](figures/comparison_unprescaled.png)
![Sequential method changes](figures/decomposition.png)
![Normalization sensitivity](figures/normalization_effect.png)

For fixed physics count n, charge Q and other corrections,
`Y=n/(Q*L)`, so replacing only L gives
`Y_new/Y_old-1=L_old/L_new-1`. The plot uses **old and new calculations on
the same cache**. It illustrates normalization sensitivity; it is not an
observed cross-section bias or a claim that the new ratio is certified total
livetime. Endpoint errors scale with file boundaries and run duration, current
lag with trips/ramps, and sampling with prescale/trigger configuration. These
quantities correlate with run groups and rates and can therefore create false
current trends, steps or kinematic differences. Missing coverage is not a
simple number-of-files correction. A plausible global ratio can conceal
counter defects that partially cancel real losses.

## Uncertainty and coverage limits

The saved NewGen error propagates independent Poisson numerator and denominator
errors. Its E and D are not automatically independent. Conditional on D and
independent identical pulse outcomes, with `q=E/D`, a Bernoulli scale is
`sigma(pE/D)=p*sqrt(q*(1-q)/D)`. The analogous physics scale uses
`q=(N-E)/(S-D)` and denominator S-D. No binomial error is assigned when those
raw fractions are outside [0,1]. They are not clipped to obtain an error bar.

The archive also has 10/20/60-s block-resampling spreads: 1000 replicates,
run-number seed, blocks grouped within each continuous component. These measure
observed time variation under a resampling model; they do not establish
hardware, boundary or periodic-sampling systematics. One block length is not
chosen because it gives a preferred error. The current-stability comparison
requires both adjacent current records to pass their windows; that changes
the population and is not solely a measurement-error estimate.

Across all calculated runs, the maximum absolute EDTM boundary-convention change is 0.00020601566; the maximum absolute +/-2 to +/-10-ns change is 6.0663047e-05. These are observed variations, not adopted systematic errors. Every per-run variation is in `run_summary.csv` and the JSON/timing tables.

PRE40/100/150/200 are populated, but the trigger monitored by each h/p channel,
its actual pulse width and its relevance to the selected trigger remain
unknown. Nonzero values, correlations and a historical width formula do not
establish attribution. Earlier segment-0 hPRE150/hPRE100=0.99067720 and
pPRE150/pPRE100=0.95547866 are not interchangeable. Neither was applied here.
EDTM_CP is a copy monitor, not an accepted-pulser denominator. An unbuffered
250-us bias model from another DAQ configuration was not transferred to these
runs without TI buffering/readout evidence.

TRIG-input -> L1 -> recorded-event comparisons address the stages between
those taps. EDTM samples only the paths downstream of its actual injection.
NPS VTP/cluster/waveform response, analog/detector behavior and upstream trigger
losses require separate coverage evidence. If EDTM already includes computer
losses, multiplying EDTM survival by CLT counts those losses twice.

All recorded events contain the configured trigger bit. Of 168 embedded
Run_Data records, 40 have populated prescales agreeing with the metadata;
128 have all-zero placeholders. Every prescale set/read flag is zero. These
are not independently decoded TI run-start settings, and an all-zero later-
segment record is not evidence that its trigger was disabled. Exact original
replay executables and run-period
EDTM/PRE routing remain unresolved. Earlier source/job evidence for **4398**
must not be asserted as byte-level provenance for every run or both sources.
The current-only ratio must also be reconciled with the actual physics
good-helicity intervals and charge integration before adoption. A single
unweighted run-average LT is not a replacement for exposure-consistent
normalization.

## Reproduce and verify

The following commands were executed from the working directory
`/tmp/nps_x60_4b_lh2_livetime_20260909`; the compact archive omits the large
binary column caches. Copy the archive to a writable scratch directory before
recreating them. Python needs NumPy and Matplotlib; ROOT/hcana comes from the
Hall C environment. The current manifest is frozen, so a later missing cache
file must be reported, never fetched automatically.

```bash
cd /tmp/nps_x60_4b_lh2_livetime_20260909
# The original inventory.py invocation snapshots inputs; do not overwrite an
# archived snapshot merely to rerun on the existing frozen manifest.
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; hcana -l -b -q "export_cached.C(0,0,3)"' > export_batch0.log 2>&1 &
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; hcana -l -b -q "export_cached.C(0,1,3)"' > export_batch1.log 2>&1 &
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; hcana -l -b -q "export_cached.C(0,2,3)"' > export_batch2.log 2>&1 &
wait
python3 analyze.py > final_analysis.log 2>&1
python3 audit.py
python3 independent_checks.py
python3 make_figures.py
python3 make_deliverables.py
sha256sum -c SHA256SUMS
```

`*.done` appears only after both T and TSH column exports finish. The audit
requires all manifest runs to have final statuses, verifies current-window
reproduction, checks the 42 identical-coverage saved counts, and rechecks ROOT
size/mtime. Independent uproot comparisons check representative T and TSH
values, including the 4305 mismatch and the 4350 anomalous counters;
`independent_reader_checks.json` records these comparisons. `validation.json`
gives the exact check totals. No full production
efficiency or cross-section job was run. `run_summary.csv`, per-run JSONs,
interval CSVs, timing histograms, selection records, logs and the source hashes
make each reported number traceable. `SHA256SUMS` verifies installed artifacts;
it is not a hash of the original ROOT payloads.

## Source trail

- Current implementation: snapshotted `compute_efficiencies_stuff.cxx`,
  `newgen_edtm_livetime.h`, `prescale_beamtime_helper.h`,
  `good_event_selection_helper.h`, `NPS_selection_helper.h`,
  `root_file_discovery.h`; hashes in `source_provenance.json`.
- Applied baseline: user's multipanel PNG and saved numerical/selection CSVs
  under `snapshot/`. Run metadata: snapshotted master CSV.
- Prior 4398 findings and code/replay evidence:
  `../livetime_4398_check_20260909/REPORT.md`, especially endpoint/current,
  counter mappings, embedded metadata and full-run sections.
- Original May 2025 talk: sibling source archive, slides 3--11, 20 and 23;
  [original PDF](https://indico.jlab.org/event/946/contributions/16510/attachments/12602/20074/Luminosity_nps_collaboration_meeting_2025.pdf).
- Source archive D01 Pooser 1022-v1, slides 7--10 and 14--15 (definitions);
  D02 Mack 1001-v2, slides 2, 4--7 (conditional non-Poisson model);
  D03 trigger electronics 1028-v5, printed p14 (topology); NPS Sept-2023 run
  plan, printed p5 (historical trigger labels). These provide context, not
  run-specific routing proof. See `../SOURCES_LIVETIME_4398.md` for URLs/hashes.

## All-run ratio table

| Run | Type | Source | Segments | p | Saved NewGen | Same-cache NewGen | Matched EDTM | Physics CLT | Flags |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 4253 | production | updated | 0;2 | 1 | 1.00039076 | 1.0003907585600547 | 0.9997801983099692 | 0.9997665813952755 | multiple_coverage_components;partial_catalog_coverage |
| 4254 | production | updated | 0;1;2 | 1 | 1.00039931 | 1.0003993099923332 | 0.9994895354772844 | 0.999667176995274 | CALCULATED |
| 4255 | production | updated | 0 | 1 | 1.0 | 1.0 | 0.9999186969836581 | 0.9999732075876112 | CALCULATED |
| 4256 | production | updated | 0 | 1 | 0.99998163 | 0.9999816279475662 | 0.9999540698689154 | 0.9999409785752228 | CALCULATED |
| 4257 | efficiency | updated | 0 | 1 |  | 0.9997334517992004 | 0.9993558418480675 | 0.9994104249062324 | coverage_differs_from_saved |
| 4258 | efficiency | updated | 0;1 | 1 |  | 1.002607675257508 | 0.9995444487830275 | 0.9996790856190095 | coverage_differs_from_saved |
| 4259 | production | updated | 0 | 5 | 1.03197539 | 1.0319753892048102 | 1.0304216773872783 | 0.9952254556698216 | EDTM_above_one;prescale_sampling_conditional |
| 4300 | production | updated | 0;1;2;3;5 | 1 | 0.99884534 | 0.9988453376607165 | 0.9987111794786513 | 0.9987491457503327 | multiple_coverage_components;partial_catalog_coverage |
| 4301 | production | updated | 0;1;2;3;4;5 | 1 | 0.9976526 | 0.9976526033531659 | 0.9990790982385497 | 0.9990894709565477 | CALCULATED |
| 4302 | production | updated | 0;1;2;3;4;5;6;7 | 1 | 1.00165841 | 1.0016584099024715 | 0.9991635093244108 | 0.9992163875523277 | CALCULATED |
| 4303 | production | updated | 0;1;2;3;4;5;6 | 1 | 0.99766461 | 0.997664608720949 | 0.9997819556158543 | 0.9997500641326204 | recorded_exceeds_L1_interval |
| 4304 | production | updated | 0;1 | 1 | 0.99917579 | 0.999175794382388 | 0.9991339179387247 | 0.9991813632296042 | CALCULATED |
| 4305 | production | updated | 0 | 1 | 1.0016972 | 1.0016971996206259 | 1.0015474467129237 | 1.001596744649399 | EDTM_above_one;CLT_above_one;recorded_exceeds_L1_interval |
| 4306 | production | updated | 0 | 1 | 0.99988414 | 0.9998841373499052 | 0.999852538445334 | 0.9998212002523438 | CALCULATED |
| 4307 | junk | none | none |  |  |  |  |  | NO_CACHED_FILES |
| 4308 | production | updated | 0;1;2;3;4;5 | 2 | 1.00232144 | 1.0023214447224817 | 0.9953826447019571 | 0.9969561094537874 | prescale_sampling_conditional |
| 4349 | production | production | 0;1;2;3;4;5 | 1 | 0.99975407 | 0.999754073163234 | 0.9993149180975802 | 0.9993811252436242 | CALCULATED |
| 4350 | production | production | 0;1;2;3;4;5 | 1 | 1.00132499 | 1.0013249921920102 | 0.9989888394333721 | 0.9991332036830061 | multiple_coverage_components;cross_segment_counter_or_clock_restart;scaler_jump_or_clock_anomaly |
| 4351 | production | updated | 0;1;2;3;4;5;7 | 1 | 0.99856876 | 0.9985687575491834 | 0.9989327729702003 | 0.9990937209908162 | multiple_coverage_components;partial_catalog_coverage |
| 4352 | production | updated | 7 | 1 | 0.9979351 | 0.9979351032448378 | 0.9979351032448378 | 0.9986903202762234 | missing_start_segment;partial_catalog_coverage |
| 4353 | production | production | 0;1 | 1 | 1.00091745 | 1.0009174494904167 | 0.9995420060136602 | 0.9995354659989878 | CALCULATED |
| 4354 | production | updated | 1 | 1 | 0.99986072 | 0.9998607242339833 | 0.9998607242339833 | 0.9996135969706003 | missing_start_segment;partial_catalog_coverage |
| 4355 | production | production | 0 | 1 | 0.9998164 | 0.9998163995960792 | 0.999795999551199 | 0.9998432624675012 | CALCULATED |
| 4356 | production | updated | 1;3;4;5;6 | 2 | 0.99829727 | 0.9982972714813738 | 0.9945310251243883 | 0.9961644490399505 | multiple_coverage_components;missing_start_segment;prescale_sampling_conditional;partial_catalog_coverage |
| 4397 | production | production | 0;1;2;3;4;5;6 | 1 | 1.00087975 | 1.0008797542755299 | 0.9993196250377986 | 0.9993179319207556 | CALCULATED |
| 4398 | production | updated | 0;1;2;3;4;5 | 1 | 1.00160539 | 0.9985808809181701 | 0.9991485285509021 | 0.9991685076400997 | coverage_differs_from_saved |
| 4399 | production | production | 0;1;2;3;4;5;6 | 1 | 1.00130973 | 1.001309732615406 | 0.9991529368670874 | 0.9991735712709365 | CALCULATED |
| 4400 | production | production | 0;1;2;3;4;5;6 | 1 | 1.00212966 | 1.0021296595482185 | 0.9991877177062904 | 0.999205414062634 | CALCULATED |
| 4401 | production | production | 0;1 | 1 | 0.99767994 | 0.9976799440950385 | 0.9994409503843467 | 0.9995293884898017 | CALCULATED |
| 4402 | production | production | 0 | 1 | 0.99821995 | 0.9982199506520972 | 0.9995417694747973 | 0.9996535284465821 | CALCULATED |
| 4403 | production | production | 0 | 1 | 0.99985178 | 0.9998517825992772 | 0.9997947759066914 | 0.9997224565121795 | CALCULATED |
| 4404 | production | production | 0;1;2;3;4;5 | 2 | 0.99891162 | 0.998911617782491 | 0.9930448222565688 | 0.996374898422675 | prescale_sampling_conditional |
| 4483 | junk | updated | 1 | 1 |  | 1.0103875349580504 | 0.999001198561726 | 0.9992002326595899 | missing_start_segment;coverage_differs_from_saved;partial_catalog_coverage |
| 4484 | production | production | 0;1;2;3;4;5;6 | 1 | 1.00052066 | 1.0005206648732654 | 0.9988518582595429 | 0.9990029074219635 | CALCULATED |
| 4485 | production | production | 0;1;2;3;4;5 | 1 | 1.00052449 | 1.0005244875498693 | 0.9987729956068094 | 0.9988388998654093 | CALCULATED |
| 4486 | junk | none | none |  |  |  |  |  | NO_CACHED_FILES |
| 4487 | production | production | 0;1;2 | 1 | 0.9983163 | 0.9983163005599279 | 0.9990406828771683 | 0.9990231813251176 | CALCULATED |
| 4488 | production | updated | 4 | 1 | 1.00100923 | 1.0010092344956352 | 0.9988898420548015 | 0.9992362817362818 | missing_start_segment;partial_catalog_coverage |
| 4489 | production | production | 0;1;2;3;4;5;6 | 1 | 0.99947196 | 0.9994719577566206 | 0.9990337778487582 | 0.9991470422596819 | CALCULATED |
| 4490 | production | updated | 1 | 1 | 0.99979398 | 0.9997939843428101 | 0.9997939843428101 | 0.9994804046210823 | missing_start_segment;partial_catalog_coverage |
| 4491 | production | production | 0 | 1 | 0.99839398 | 0.9983939840532063 | 0.9997732683369233 | 0.9997254853985087 | CALCULATED |
| 4492 | production | production | 0 | 1 | 0.99995268 | 0.9999526772827296 | 0.999940846603412 | 0.9998297034248123 | CALCULATED |
| 4493 | production | production | 0;1;2;3;4;5 | 2 | 0.9852143 | 0.9852142990829122 | 0.9928255037744089 | 0.9981827527855911 | prescale_sampling_conditional |
| 4551 | production | production | 0;1;2;3;4;5 | 1 | 0.9995637 | 0.9995636963003062 | 0.9992011297155537 | 0.999295001887861 | CALCULATED |
| 4552 | junk | none | none |  |  |  |  |  | NO_CACHED_FILES |
| 4553 | junk | none | none |  |  |  |  |  | NO_CACHED_FILES |
| 4554 | junk | none | none |  |  |  |  |  | NO_CACHED_FILES |
| 4555 | production | production | 0;1;2;3;4;5 | 1 | 0.99952725 | 0.9995272542391023 | 0.9990202191146343 | 0.9990444147891844 | CALCULATED |
| 4556 | production | production | 0;1 | 1 | 1.00270595 | 1.0027059465357337 | 0.9992157717024289 | 0.9992493108908569 | CALCULATED |
| 4557 | production | production | 0 | 1 | 0.99993474 | 0.9999347364986132 | 0.999923859248382 | 0.9997971947000215 | CALCULATED |
| 4558 | production | production | 0;1;2;3;4;5 | 2 | 0.99844555 | 0.9984455484775518 | 0.9988777336285601 | 0.9971527223385245 | prescale_sampling_conditional |
