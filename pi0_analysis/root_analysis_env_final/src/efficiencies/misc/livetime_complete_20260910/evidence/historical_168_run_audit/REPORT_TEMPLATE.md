# KinC_x60_4b LH2: livetime, coverage and artificial structure

2026-09-09. **Diagnostic measurements, not a certified total-livetime correction.**
Production efficiency code, applied efficiency CSVs and physics yields were not
modified. PRE monitors were not used as a correction. No jcache command or tape
retrieval was performed in this extension.

## Results and coverage

@@SUMMARY@@

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

@@COVERAGE_TABLE@@

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

@@PRESCALES@@

## Counter defects can survive endpoint matching

@@COUNTER_AUDIT@@

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

@@SENSITIVITIES@@

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

@@ALL_RUN_TABLE@@
