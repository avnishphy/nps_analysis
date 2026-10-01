# KinC_x60_4b LH2 livetime investigation: complete technical record

Edition: 10 September 2026. Editable synthesis of the source, ROOT, raw-EVIO,
logbook and code investigations. This edition supersedes the earlier
168-segment presentation. The earlier evidence remains frozen and accessible.
The user explicitly authorized this documentation and LaTeX update; the previous
instruction to hold presentation updates is superseded. Production corrections,
physics yields and efficiency CSVs have not been changed.

## 1. Present conclusion and the question still open

The core EDTM estimator is the user's existing method: $pE/D$, with a broad
raw-TDC window. The investigation did not discover a fundamentally different
EDTM estimator. Its numerical refinements align event/scaler exposure, current
assignment and segment joins. Its other contributions are tests, failure
localization, source provenance and a more precise physical interpretation.

The exposure-matched broad-window ratio is our preferred current EDTM
diagnostic. It is not yet a certified universal total-physics correction.
Runs 4303/4305 have scaler exposure defects. The direct counting-house EDTM copy
bypasses NPS cluster-trigger formation. Injection fanout, the sent-counter tap,
accepted-pulse causality, periodic sampling and matched yield/charge exposure
still set limits on interpretation.

**CLT_physics is only a proposed definition. It is never the baseline, truth,
calibration target or validation standard.** All comparisons with it test
mutual consistency; neither agreement nor a flatter run dependence validates
either quantity. Historical variable names are preserved in data files for
reproducibility, without endorsing their physical labels.

The full-cache calculation covers 43 LH2 production runs,183 updated segments,
77,838,384 events. The five TI6 runs are coincidence-triggered;38 TI4 runs use
HMS EL-REAL singles under the supplied working routing. This distinction is
essential to deciding which losses belong in normalization. No extra NPS
hardware-trigger factor should be imposed on the TI4 population merely because
it is needed for a coincidence-triggered population.

## 2. Evidence hierarchy, chronology and scope

We distinguish measured file contents, inspected implementation, reported
hardware configuration, interpretation under assumptions, and unresolved
questions. A source revision is not proof of the executable loaded online.
A source comment is not a wiring measurement. A run plan is not a run-start
readback. An original log entry plus the user's applicability confirmation is
stronger routing evidence than a guessed meaning of a branch name.

| Stage | Evidence added | Present interpretation |
|---|---|---|
| User's May 2025 work | Four livetime definitions and unresolved luminosity behavior | Starting point; historical samples differ |
| 4398 segment0 | Reproduced E19341/D19310; uncovered49-pulse tail | Endpoint mismatch demonstrated |
| 4398 all six segments | Contiguous events and matched counters | Full-run raw ratio 0.99921948; no total-LT certification |
| First all-run audit | 168 segments,70,155,126 events including controls | Historical coverage; superseded for current numbers |
| Full-cache continuation | 183 production segments; timing/phase and statistics | Current numerical evidence |
| One raw file4305 | All relevant counters and events agree with ROOT | Deficit already in raw, cause unresolved |
| User logs4303/4305 | TI slave status and ROC endpoint totals | Acceptance/readout checks, not sent-trigger denominators |
| CODA and local source review | Buffer fields, timer semantics, replay code | Version-qualified support, not exact ROC5 readout source |
| Dated August 2023 entries | Trigger namespaces and EDTM route | Adopted working map after user confirmation |
| This edition | Consolidated documentation and new comparison layouts | No new event calculation or correction change |

All 51 LH2 metadata entries were inventoried. With cached controls/junk,46 runs
and 187 segments contain 79,564,296 events. Missing junk-cache runs 4307, 4486, 
4552, 4553, 4554 were not assigned livetimes. No production catalog gaps remain
relative to the selected updated-source catalog. ROOT data from different
replay variants were not mixed within a nominal run calculation.

The user later authorized one raw file only,4305.dat.0, and that allowance was
used. No additional raw staging is authorized or performed in this edition.
Earlier4398 ROOT staging and its request 100230264 are historical actions;
they are not permission to submit new requests.

## 3. Definitions: the population and endpoints come first

Let G be the eligible physics population with fixed detector/trigger criteria,
B the formation of the selected trigger input, P passage through intentional
prescaling, A an issued accept, and R a recorded event. Define the upstream
survival and downstream DAQ survival by their endpoints:

$$
L_{\rm upstream}=\Pr(B\mid G),\qquad
L_{\rm DAQ}=\Pr(R\mid P,B,G).
$$

$$
\Pr(R\mid G)=\Pr(B\mid G)\Pr(P\mid B,G)\Pr(R\mid P,B,G).
$$

This product is a conditional-probability chain, not an independence
assumption. If prescale sampling is representative, the middle probability
is1/p. Calling the prescale-removed product total livetime requires a precise
definition of G and which losses have been excluded from it. Threshold
acceptance, detector response, tracking/PID, NPS readout and reconstruction
efficiencies require their own populations. They must not be counted twice.

| Quantity | Operational definition | Why it differs from another LT |
|---|---|---|
| Time availability | Live duration divided by live plus busy duration at a specified decision point | Requires valid timer exposure and busy-source definition |
| Trigger acceptance | Accepted requests divided by eligible input requests | Depends on prescale, overlap and which input tap is counted |
| Readout survival | Recorded accepted events divided by issued accepts | Cannot see losses before accept generation |
| EDTM survival | Recorded injected probes divided by sent probes, with justified prescale treatment | Samples the injected route and pulse times |
| Total physics survival | Survival of the defined physics population through all corrected stages | Not supplied automatically by any branch named livetime |

The labels electronic LT, computer LT, DAQ LT and total LT are not interchangeable
without endpoints. In the conventional factorization, electronic LT describes
survival through the specified upstream trigger electronics; computer/DAQ LT
describes acceptance/readout downstream of the specified trigger-request tap.
The total is their conditional product only when the factors partition the
same eligible population without overlap. Deadtime fraction is 1 minus the
corresponding livetime, not necessarily a duration per event. The local variable
`CLT_physics` does not become a measurement of conventional computer LT by name.
Likewise, a scaler L1/input ratio is a computer-LT candidate only after those
two taps, their prescale placement and their exposure are established.

With ROC lock, a readout busy can inhibit the next accept while readout proceeds.
With buffering, some events can be accepted during earlier events' readout;
finite buffer occupancy and enabled external busy still limit acceptance.
Block level (events per block), buffer threshold (queued blocks) and individual
event readout duration are different quantities. A slave setting of block level
1 therefore does not prove an unbuffered whole DAQ. A formula that assigns one
fixed readout dead interval after every accepted event requires validation for
the actual mode; the archived unbuffered model is not applied to these runs.

For availability a(t), eligible input rate r(t), and beam current I(t):

$$
L_t=\frac{\int a(t)dt}{T},\qquad
L_r=\frac{\int r(t)a(t)dt}{\int r(t)dt},\qquad
L_Q=\frac{\int I(t)L(t)dt}{\int I(t)dt}.
$$

The last expression assumes the relevant ideal yield scales with charge.
A periodic pulser samples availability at its own times. A shared current cut
helps match exposure but does not prove equality of these weightings. Good
helicity, physics-yield and charge selection must be matched before applying
an effective production factor. Mean per-run or per-file ratios should not
replace a count/exposure combination without justified weights.

## 4. How Hall C uses EDTM, and what transfers to NPS

Pooser's Hall C formalism describes artificial pulses injected into trigger
legs near the detector electronics, counted in scalers and recognized in
recorded TDC data. The intended advantage is to sample the acquisition chain
with a known probe. That design principle does not imply every experimental
implementation injects upstream of every loss. Pooser also explicitly
distinguishes ROC-lock from buffered behavior. The title date says 2018 while
body/metadata indicate 2019; the source catalog records this ambiguity. [D01]

Mack's periodic-pulser discussion concerns an explicitly unbuffered model and
beam-ramp effects. It supplies a warning about non-Poisson sampling, not a
ready-made correction or universal readout time for NPS. His electronics
examples distinguish lost leading edges, extended logic gates, recovery and
reshaping. A scaler's edge count and a coincidence's acceptance can respond
differently. [D02,D05]

Murphy's2022 PionLT study found meaningful corrected-time tails and
multi-trigger/prescale competition. That motivated our independent tests of
tight timing and causality. Its trigger delays, windows and conditions were
not imported numerically into NPS. In particular, a narrow peak is not by
itself proof of complete EDTM acceptance. [N01]

The NPS February/March/July2024 studies report earlier timing, scaler-repair
and rate-dependence investigations. The March workflow used an automatic raw
peak, a2-uA threshold, inferred missing-pulser counts and removal of terminal
readouts. The inspected later DT_Analyzer_new.C also adds an inferred scaler
contribution. Those choices are historical evidence to assess, not adopted
defaults for the current replay. [N02,N04,N05,N06]

The 2024 SHMS instrument draft provides another description of EDTM and DAQ
buffering, scoped to its instrument and draft version. Crafts' May 2025 NPS
presentation gives a PRE100/PRE150 convention and TRIG6 component expression.
Neither establishes our actual PRE wiring. Raydo's NPS presentations separate
trigger thresholds, readout thresholds and embedded configuration data.
Waveform/trigger agreement in accepted data cannot reveal all events that
failed to enter that sample. [D04,D06,D07,N03]

## 5. Dated NPS trigger logic and timing

The user supplied entry 4170042 dated 2023-08-23 and entry 4171556 dated 2023-08-29,
then confirmed they should apply to these runs. The URLs and supplied text
are preserved in evidence/routing. Web retrieval returned401 token_expired;
the original pages and any attachments were not independently retrieved.
We adopt the reported mapping as the working configuration without claiming
an independent hardware survey or ruling out every possible later change.

Three namespaces must remain distinct: VTP internal bits, local NPS TS/cable
outputs, and counting-house ROC1 TI inputs after combinational logic.

| Number | Local NPS TS / outgoing cable, Aug 23 | ROC1 TI input, Aug 29 |
|---|---|---|
| 1 | NPS single cluster | Single cluster OR two clusters OR delayed EDTM |
| 2 | Cosmic scintillator | NPS cosmic OR LED |
| 3 | Cosmic-column OR of crates | HMS h3/4 |
| 4 | Cosmic-column AND of crates | HMS hEL-REAL |
| 5 | VLD | TI1 AND TI3 |
| 6 | At least two clusters | TI1 AND TI4 |

Local NPS T6 is not ROC1 TI6. The latter is a coincidence with HMS EL-REAL.
NPSlib's internal VTP bits 0--5 are labeled cluster singles, cosmic scintillator,
cosmic column, different-crate cluster pair, same-crate cluster pair and VLD.
Those labels are a third mapping, not a contradiction resolved by renumbering
the global TI bits. The selected bit in the raw TI mask is separate evidence.

The second entry reports the following physical arrangement:

- HMS h3/4 and hEL-REAL are delayed about 3.5us using three counting-house to
  SHMS-hut loops and an additional hall-floor loop. Their pretrigger timing
  difference was40 ns; an additional27 ns was applied to h3/4.
- HMS EDTM enters discriminators, about 52 ns ahead of h3/4 at the counting
  house, and follows the HMS delay. The precise fanout and hEL-REAL formation
  from injected pulses are not fully specified in the excerpt.
- A separate EDTM copy is delayed about 3.5us and OR'ed at the counting house
  with NPS single-cluster OR two-cluster outputs.
- NPS outputs receive100 ns delay before that OR and 50 ns afterward. At
  coincidence formation, HMS pulses are40 ns wide and NPS100 ns; HMS arrives
  about 40 ns later and sets the coincidence timing.
- Singles are delayed further relative to coincidence before reaching TI:
  NPS about 40 ns, HMS about 24 ns. The actual enabled-input arbitration cannot
  be reconstructed from these approximate delays alone.

The first entry describes firmware capabilities: cluster hit coincidence
extended from +/-16ns to as much as +/-48ns;5x5 or7x7 waveform regions; full
readout every N events up to 65535; a separate two-cluster threshold; and VLD
traversing FADC->VTP->V1495->TS. It also says driver/config support was pending.
Capabilities are not active run settings. Waveform-pattern N is not
automatically the trigger prescale p. These approximate timing statements
must not be combined into a calibrated propagation model without common
reference points. A40 ns pulse width does not identify the PRE40 scaler.

## 6. What EDTM can and cannot establish for these triggers

The direct delayed EDTM copy joins the NPS branch at the counting house.
It therefore bypasses NPS FADC/VTP/V1495 cluster-trigger formation. Its
survival cannot certify the efficiency of those bypassed stages. VLD's
reported full electronics route is different from that direct EDTM route.
Do not equate a VLD test or VTP bit 5 with an EDTM acceptance measurement.

| Cohort or question | What the EDTM diagnostic can probe, conditionally | What it cannot establish alone |
|---|---|---|
| TI4 HMS EL-REAL singles | Survival of the sampled HMS/injection/DAQ path | NPS reconstruction/readout efficiency or a nonexistent NPS trigger requirement |
| TI6 NPS AND HMS EL-REAL | Survival of an injected coincidence when both pulse branches reach and overlap | Upstream NPS physics-cluster formation that EDTM bypasses |
| Recorded TDC pulse association | Whether timing agrees with an injected pulse | Which signal caused the global accept |
| Ratio near unity | Numerical agreement of counted populations | Correct wiring, complete loss coverage or sound scaler exposure |
| Agreement with proposed physics ratio | Conditional consistency on shared counts | Independent validation or a truth baseline |
| Cross-section normalization | A possible factor after population/exposure validation | Automatic correction of trigger thresholds, PID, readout or missing data |

Sent-counter tap, injected branch fanout, pulse formation, timing acceptance,
prescale sampling and overlap remain explicit conditions. Electronic rate
losses and fixed threshold/trigger acceptance should not all be called DAQ
deadtime. For TI6 an upstream NPS trigger-formation factor may need independent
study for the target population. For TI4 that hardware cluster condition is
not a requirement of the configured input; adding such a factor simply by
analogy with TI6 would be inappropriate. NPS readout/selection still matters.

Frozen metadata agrees with all 43 configured trigger/prescale choices. TI6
runs 4253, 4254, 4255, 4256, 4259 (2024-02-08) contain 2,706,910 whole-file events;
four have p1 and 4259 p5. The 38 TI4 runs start 2024-02-10 and contain 75,131,474
events;33 have p1 and five p2. These are inventory counts, not selected E/N.
The descriptive coin_status field is not the enabled-input definition:
4259 says off there but explicitly selects ps6=3, effective factor 5.

For a sole enabled trigger with p=1, EDTM survival needs no intentional-prescale
factor, although tag acceptance, exposure and path assumptions still apply.
For a sole enabled prescaled trigger, multiplying E/D by p removes intentional
thinning only under representative sampling of pulse opportunities. Values
above one are possible and must not be clipped. If several triggers are enabled,
one pulse can populate several inputs while one event is recorded after their
combined acceptance logic. A single input's p is then not automatically the
correct factor for the pooled E. Even the familiar OR probability
1 minus the product of individual rejection probabilities assumes independent
decisions; synchronous pulses, shared counters and correlated prescale states
can violate that assumption. Trigger priority, timing, overlap and busy must be
established from configuration and data. Multiple bits in a recorded mask do
not by themselves establish multiple enabled, independently accepting triggers.
Our 43-run configuration inventory supports one enabled input per run; it does
not recover the complete hardware state of every prescaler.

## 7. Earlier methods and why they did not settle the issue

The user's May 2025 luminosity work and inspected older code compared four
quantities. Let H0 be the whole-file nonzero selected-trigger TDC count, F0
the corresponding raw-EDTM-zero count, E0 positive EDTM, and subscript c
current-selected scaler values. Define fA=Ac/A0 and fD=Dc/D0 where valid.

$$
L_{\rm EDTM,old}=\frac{pE_0}{D_c/f_D}=\frac{pE_0}{D_0},\quad
L_{\rm TSH,old}=\frac{pA_c}{S_c},\quad
L_{\rm TDC,old}=\frac{pH_0f_A}{S_c},\quad
L_{\rm phys,old}=\frac{pF_0f_A}{S_c-D_c}.
$$

The simplification assumes a valid, nonzero, unclipped fraction. The inspected
source clamps beam fractions; pathological fractions cannot be simplified
blindly. Multiplying a whole-file TDC count by a beam-on fraction is a proxy,
not direct event-to-scaler exposure assignment. The scaler L1/input ratio
measures only its tap-to-tap population. EDTM subtraction also requires the
right accepted and sent populations. These are separate shortcomings, not
evidence that one alternate definition is the true LT.

The older active nps_livetimes.h reads database CPU_LT; a former calculator
is commented out. For 4398 that database gives 0.9992 with six segments and a
2-uA condition, whereas the master metadata Computer_live_time is 0.992 with
unresolved differing provenance. Agreement of a diagnostic with a rounded
database value supports at most component consistency. Neither value should
silently replace the current NewGen reference.

Historical NPS event-spacing fits, PRE formulas and scaler repairs remain
conditional alternatives. No shifted-exponential fit, width-ratio correction,
unbuffered 250-us correction, empirical rate slope or missing-pulse repair was
adopted in this calculation. No ratio is clipped to unity.

## 8. Original NewGen versus the actual refinements

The original NewGen implementation already computes $pE/D$, uses a run raw
peak and a +/-500 raw-channel window, and uses the existing per-file current
helper. The width is about +/-48.83ns at the historical0.09766ns/channel
conversion, not +/-500ns. The helper finds a current peak in 0.5-uA bins and
uses +/-15%; its threshold definition was retained.

Original numerator: current stored with each T event, finite raw EDTM greater
than 1, and the broad window. Original denominator: per-file consecutive TSH
EDTM differences passing the ending-row current window. That asymmetry matters
because event-stored current normally reflects the preceding scaler snapshot.

The matched diagnostic changes only the following nominal count bookkeeping:

1. Assign events to `[TSH.evNumber_i,TSH.evNumber_(i+1))` and apply the same
   ending-snapshot current decision to events and scaler increments.
2. Use fully event-covered scaler intervals for both numerator and denominator;
   do not count unsupported event heads/tails or initial cumulative exposure.
3. Join adjacent segments only with exact event continuity and nondecreasing
   relevant counters/clock; keep gaps/restarts separate and handle identical
   counter/clock snapshots explicitly. Storage boundaries are not automatically
   exposure boundaries.

This is not a new physics/PID cut or a completed good-helicity/yield/charge
integration. Counter integrity checks and timing studies are validation, not
a replacement estimator. The original independent-Poisson-style error
function was inspected, not rewritten in production.

For 4398 on the same six updated files, p1:

| Step | E | D | Ratio |
|---|---:|---:|---:|
| Original NewGen | 112586 | 112746 | 0.9985808809 |
| Restrict numerator to per-file coverage; old current | 112263 | 112746 | 0.9957160343 |
| Align to ending-interval current | 112658 | 112746 | 0.9992194845 |
| Final stitched diagnostic | 112658 | 112746 | 0.9992194845 |

Net E change is+72, not merely removal of edge tags. These sequential changes
are order-dependent and this run's decomposition is not universal. The saved
one-segment NewGen ratio 1.00160539 is a different input comparison. It must not
be presented as the same-file method change to 0.99921948.

![4398 refinement steps](figures/run4398_refinement_steps.png)

![Sequential changes on the same files](figures/method_decomposition_current.png)

## 9. Exact common-interval quantities and boundary checks

On identical intervals define N=recorded events, E=broad-window EDTM-tagged
subset, D=sent EDTM increments, S=configured input increments, A=L1 increments,
and p=effective prescale. The four displayed ratios are:

$$
\widehat L_E=\frac{pE}{D},\quad \widehat L_{\rm all}=\frac{pN}{S},\quad
\widehat L_{\rm phys}=\frac{p(N-E)}{S-D},\quad \widehat L_{A}=\frac{pA}{S}.
$$

$$
\widehat L_{\rm all}=(1-D/S)\widehat L_{\rm phys}+(D/S)\widehat L_E.
$$

The exact shared-count identity was verified for all 43 production runs.
Accepted E belongs in the numerator subtraction; sent D in the denominator.
Physical interpretation of the subtracted quantity assumes each relevant sent
pulse contributes to S and E identifies its accepted subset, with suitable
prescale sampling. The quantities are statistically correlated and not
independent proofs of each other. Zero/invalid denominators remain invalid.

The nominal half-open convention has a checked alternative using the opposite
boundary inclusion. A recorded event at a latch boundary does not prove an
exact hardware latch order. Initial baselines and duplicated terminal rows
are treated separately. The 4398 first snapshot at event 1 already has 878 EDTM,
23020 trigger inputs and 1 L1 at clock 21.963170s; these are not counts from
an observed event interval beginning at the first event.

In4398 segment0, the last real boundary 541321 and terminal boundary 542585
have identical scaler values/time. The unsupported1264 events include49
EDTM tags. The original19341/19310 ratio becomes19292/19310=0.9990678405
when that tail is removed; tight timing yields19290/19310=0.9989642672.
Later segments supply real snapshots and recover much of the tail exposure;
discarding every file tail independently is not the final stitched method.

Full 4398 has 3,222,394 contiguous events and 2233 real scaler snapshots. The
fixed diagnostic window33.3625--45.1375 uA gives 2821.898170s and 109930.299545 uC,
N=A=2893284, S=2895694, D=112746, Eraw=112658, Etight=112650. Switching boundary
inclusion makes tight E112652, a measured LT change0.000017739. That local
boundary study has its own fixed window and must not be mislabeled as a new
all-run systematic. The one final non-EDTM event has no nominal extra exposure.

TSHelH.evcount and TSH.evcount are different record indices. In segment0 their
ranges reach 19435 and 354, respectively. Applying ranges of one to the other
produced the saved diagnostic denominator6550 and beam time 163.950216s rather
than common covered-current exposure483.315153s. This did not define the main
NewGen denominator, which remains19310. Do not promote18970/6550 into an LT.
Use common event/time coordinates for future helicity-to-scaler matching.

## 10. Timing, pulse association and rejected alternatives

Broad E requires finite raw TDC greater than 1 and +/-500 raw channels around
the run peak. Tight E requires raw greater than 1 and corrected time within
2 ns of its run peak. Every tight candidate on this sample also lies in the
broad window. Both nominal variants use the same exposure; tight is a
sensitivity diagnostic, not a demonstrated acceptance correction.

The independent pulse-phase analysis uses event-clock timestamps. A first
20-second straight-line model failed to track oscillator wander. It was
replaced by local interpolation between alternating core pulses, with the
other half held out for validation. No unsupported edge extrapolation is
used. Pulse spacing is about 6,257,191.7 ticks. Healthy scaler-clock comparisons
are consistent with roughly250MHz, but phase results stay in ticks rather
than claiming an independently calibrated absolute clock.

| Nominal production selection | Count | Within25 ticks | No interpolation support |
|---|---:|---:|---:|
| Broad raw tags | 3650941 | 3650890 | 51 |
| Tight corrected core | 3650677 | 3650626 | 51 |
| Broad tags rejected by tight | 264 | 264 | 0 |
| Positive raw outside broad window | 9395 | 153 | 1 |

All 1,825,258 modeled held-out core events pass 25 ticks;14 are outside10 ticks.
Per-segment99.9% held-out absolute residuals span6.452--10.516 ticks. Of the 264
rejected tags,226 pass 10 ticks. Using interpolation anchors themselves to claim
model precision would be circular; their fitted residuals are forced small.

For 5515 positive hits outside the broad window but within 500phase ticks,
EDTM-minus-selected-trigger timing versus phase has slope-3.9987366ns/tick
and correlation-0.9997145. This supports nearby pulses captured in a physics
event's TDC window, not adding all positive hits to accepted E. Among69,906,458
raw-zero events, none passes 25 ticks, two pass 50 ticks and 738 pass 125 ticks.
A uniform-phase background subtraction is therefore not justified.

Whole-file relative-timing +/-2ns recovers only1/14,1/3,0/0,1/6,0/4,3/12 of
the raw-versus-tight tails in 4253,4255,4259,4301,4305,4350. It does not
universally repair tight acceptance. Broad-to-tight changes any production
ratio by at most0.0001595202 (0.015952 percentage points,4254). This is a
measured sensitivity, not a complete systematic confidence interval.

![Pulse phase and corrected/relative timing](figures/refresh_phase_timing.png)

![Broad tags rejected by tight timing](figures/refresh_rejected_tags.png)

The central limitation is unchanged: pulse association is not causal trigger
identification. An injected pulse and a physics trigger may coexist in the
same recorded event. The documented downstream OR makes the distinction
physically relevant, rather than removing it.

## 11. Counter exposure and the raw 4305 investigation

An EDTM/clock comparison alone missed4303/4305 because both counters lose
exposure together. Among52,025 selected intervals with exact recorded events
at both scaler boundaries, only these two disagree by more than 0.1s:

| Run | Event boundaries | Event-clock seconds | Scaler seconds | N | A=S | Eraw | D |
|---|---|---:|---:|---:|---:|---:|---:|
| 4303 | 455000 to 456915 | 2.000258 | 0.000869 | 1915 | 3 | 79 | 1 |
| 4305 | 105080 to 106127 | 2.001508 | 0.005481 | 1047 | 2 | 80 | 0 |

The clock scale comes from healthy positive intervals; fitting across a
counter step could absorb the defect. Charge increments are also tiny.
Values remain flagged in diagnostic tables. Neither inferred repair nor a
production interval exclusion has been adopted; exclusion requires the same
yield and charge exposure treatment.

![Independent event-clock exposure test](figures/refresh_counter_exposure.png)

The one authorized file nps_coin_4305.dat.0 is 13,341,141,036 bytes, CRC32
53e8a0c0 and Adler32 756ab5b0, matching the tape stub. Request100232289 completed;
the raw file was read to EOF. There are561307 raw records,560647 physics
events and 660 relevant scaler banks. All physics event numbers/timestamps/
six-bit masks and all corresponding clock, EDTM, L1 and six trigger counters
match ROOT. ROOT's extra final scaler row repeats the last counter values
with boundary 560648 instead of 560647; it is bookkeeping, not exposure.

Raw records105247/106294, scaler banks125/126, have the same anomalous
increments. Module headers remain intact. Slots6--9 continue counting while
mapped clock/trigger counters in slots 10--12 scarcely advance. Before the
defect L1 minus event number is 0 or1; afterward it is-1046 in 524 snapshots
and-1045 in 11 through run end. The deficit does not recover at the next read.
A simply omitted cumulative snapshot would ordinarily recover its counts;
that alone is not the observed pattern.

Slot7/channel31 and slot9/channel31 are equal in all 660 snapshots. They equal
mapped EDTM for the first 125 and exceed it by80 in every remaining 535.
The maps call them Empty_12 and Empty_24. They localize a discrepancy, but
their identity is not established and they are not replacement denominators.
Hardware inhibit, FIFO loss, online accumulation/readout or another upstream
cause remain possibilities. ROOT-only conversion is ruled out for this checked
defect. Raw 4303 was not examined, so raw verification must not be extended to it.

![Persistent raw scaler deficit](figures/raw4305_exposure_step.png)

Run 4350 has a separate updated-replay restart and near-2^32/zero-clock evidence.
Its affected local anomalies are below the current band; components remain
split. This should not be conflated with4303/4305 selected exposure defects.
Ordinary losses with N=A below S remain in the LT calculation. For 4398 near
820.125614s, S2082,N=A1799,D80,E69 over about 2s is a retained loss burst, not
grounds for flattening the curve by removal.

![Run 4398 retained losses](figures/current_intervals_4398.png)

![Run 4303 interval diagnostics](figures/current_intervals_4303.png)

![Run 4305 interval diagnostics](figures/current_intervals_4305.png)

![Run 4551 interval diagnostics](figures/current_intervals_4551.png)

## 12. Prescales, covariance and uncertainty

The analysis convention maps setting0 to p1 and positive setting s to
2^(s-1)+1; setting-1 disables an input. The effective PS_N fields independently
agree with all 43 summary choices. Of 187 stored Run_Data prescale records,
45 are populated and agreeing and 142 are placeholders; set/read flags are
zero throughout. Recorded events carry the configured trigger bit, but these
facts are not a recovered hardware prescale readback or proof of multiple
enabled triggers when several mask bits are present.

Prescaling reduces the accepted pulser sample. Under independent thinning
with q=E/D, a conditional binomial scale is:

$$
\sigma(pE/D)=p\sqrt{q(1-q)/D}.
$$

For true LT1 under representative p-fold thinning, the scale becomes
sqrt((p-1)/D). Thus a corrected estimate greater than 1 need not mean accepted
pulses outnumber sent pulses. Independent Poisson propagation for E and D
ignores their shared population; neither it nor the conditional binomial
model automatically covers phase correlation, beam structure or hardware
systematics. The production error function has not been replaced.

Run 4259,p5 has N56200,E6632,D32181,S281210,A56200. Its ratios are EDTM
1.03042168, all-trigger0.99925323, proposed raw-subtracted0.99522546; difference
0.03519622. E/D itself is below 1. Under LT1/thinning, sigma0.011149 gives
an excess about 2.73 sigma. Conditional on N,S,D with exchangeable pulse labels:

$$
E\sim{\rm Hypergeometric}(S,D,N),\qquad
\mathrm{E}[E]=\frac{ND}{S},\qquad
\mathrm{Var}(E)=N\frac{D}{S}(1-D/S)\frac{S-N}{S-1}.
$$

The comparison z is 2.97; paired20-second block resampling gives 2.61. Exact
two-sided probabilities are retained in acceptance_tests.csv. Periodic pulses
need not satisfy exchangeability. Block spreads are diagnostics of dependence,
not a final systematic error; isolated outliers also require multiple-testing
context. Other prescaled comparisons have paired residuals below 1 sigma.

Every unprescaled per-run raw-minus-physics residual is below 1.7 paired sigma.
That does not exclude a collective difference. For 35 counter-clean unprescaled
runs, an inverse-variance common-mean diagnostic gives raw-minus-physics
(4.3734 +/-1.4349)e-5,3.05sigma; tight-minus-physics gives
(-3.0850 +/-1.5160)e-5,-2.03sigma. It assumes independent run errors and omits
shared systematics. Timing changes the sign; no precision equivalence or
calibration against CLT_physics follows.

The earlier fixed-window4398 study used 2000 bootstrap replicas with seed 4398
and10/20/60-second blocks, obtaining spreads near 1.3e-4 for EDTM and 1.1e-4 for
the proposed ratio. Other archive stages have different replica counts and
variants; their scripts and JSON are retained, not silently harmonized into
one claimed uncertainty. The current paired raw/tight comparison is the
full-cache acceptance_tests output. No final hardware-path uncertainty is assigned.

![Paired statistical comparisons](figures/refresh_correlated_difference.png)

## 13. DAQ, ROC and source-code evidence

The user logs for4303/4305 identify coin_sparse and ten initial/end TI status
blocks each, configured as slaves with firmware 11.3 (tip113.svf), blocklevel 1,
broadcast buffer 5 and busy enabled. The CODA TI3v11.3 header gives:

| Register | Value | What is established |
|---|---|---|
| dataFormat | 0x05050006 | Broadcast buffer level5 in bits 31:24 |
| vmeControl | 0x00800010 | Local selector bit 22 clear; busy-on-buffer bit 23 set |
| blockBuffer | 0x00000001 | Local value, not the selected threshold |
| busy | 0x00000003 | Switch-slot A/B sources enabled, not their busy duration |
| Block level | 1 | One event per configured block |

These settings support buffered operation in those slaves; another crate,
master setting or busy path can constrain effective whole-DAQ acceptance.
The current TI hardware PDF is updated 2026-08-27, after the runs. Its later
features are not assumed for 2024; the matching firmware-labelled driver
documentation is the relevant cross-check. The archived TS PDF is a2017
Version 4 document; it does not identify the actual NPS master model.

All five end TI Readout/Ack/L1A triplets equal3310432 in 4303 and 560647 in 4305.
ROC5 itself reports those same totals. ROC5 and npsvme5 are distinct event-builder
components; the five NPS TI dumps cannot be assumed to configure the HMS
scaler ROC. Agreement of event totals does not prove integrity of scaler
payloads inside those events or measure lost pre-acceptance requests.

Timer registers are not validated run-boundary exposure. Initial L1A is zero
while timers are nonzero; same-board busy can decrease. For board 0x7101150e:
4303 busy0x260a to 0x1afd;4305 busy0x670 to 0x29a. TI getters require latching;
readout APIs offer no latch, latch and latch/reset options. Exact online reset
history and tiStatus implementation are unknown. Two interleaved4305 timer
chunks were deliberately left unassigned. No timer ratio was adopted.

Raw 4305 contains five FADC and five VTP text configurations. No event120/133
prescale records were found. Types135/136 hold uuencoded angle photographs,
not trigger logic despite generic decoder labels. Literal VTP trigger widths,
thresholds, latency and internal prescales require firmware units/meaning;
they are not a global TI configuration.

The NPSlib event 137 handler, added 2024-01-11, parses a unique-key std::map with
emplace. Within one instance it retains the first value per key, without
crate/index qualification. Each archived4305 VTP config has 32 prescale rows;
latencies in record order are2560,2556,2556,2536,2536. A first-key model retains
only prescale0 1 and latency 2560. This is an inspected-source/model limitation,
not a reproduced pass 2 bug claim. Original indexed raw text remains preserved.
Do not infer a full configuration from one GetInfo value.

NPSlib also contains the September 2023 T1=NPS VTP OR EDTM timing comment, now
supported by stronger dated logbook routing. Its VTP decision-time mask 0x7F
conflicts with an eleven-bit comment; firmware must be checked before using
that field. Previously analyzed TI timestamps and EDTM TDC values are different
fields and were not modified. No decoder patch was made.

The searches of /u/group/nps/apps and /group/nps/apps/NPSlib found replay and
monitoring software, not the online ROC5 readout list. The two NPSlib paths
are the same filesystem directory. The 2562-file apps inventory, precise
negative searches, source snapshots and Git revisions are retained. The
three hcana scaler handlers1.0.0/1.0.1/1.0.2 are identical and clean; they
decode recorded banks and extend wraps offline. They do not reveal the
physical FIFO/online accumulator. A /cdaqfs1/apps modulefile clue was found,
but that filesystem and /home/coda are absent here.

CODA configuration examples identify rol1/rol2 shared-library paths. The
concrete missing source artifact is the run-period coin_sparse ROC5
configuration and matching source/build/driver. Source snapshots identify
what was inspected, not an executable proven to produce every pass 2 file.
For 4398 the recorded replay job names the updated output and an input tar
with HEAD2d462a85b860c54785d8287e644619a1cb1acc19; its original job logs are
absent. The tar working contents are archived. This provenance is not
silently extended to every run or alternate replay.

## 14. PRE and other fallback models: what remains unlicensed by evidence

PRE40/100/150/200 exist and are populated. The unresolved issue is their
upstream trigger attribution, actual pulse widths, module response and gating.
Mapped ROC5 PRE100/150 use slot 10 channels9/10. In the original4398 segment0
selected diagnostic, PRE100=612581944 and PRE150=606870963, ratio 0.99067720;
PRE200/PRE100=0.94158061; PRE40 is almost PRE100. These are not full-run totals.
The pPRE150/pPRE100 ratio 0.95547866 is another channel population. None can
be chosen by closeness to EDTM or by matching a desired correction.

For an explicitly hypothetical nonparalyzable Poisson monitor of rate R and
dead width w over time T, Mw=RT/(1+Rw). For w2 greater than w1:

$$
r=\frac{M_{w_2}}{M_{w_1}}=\frac{1+Rw_1}{1+Rw_2},\qquad
1-L_{\rm path}\simeq \frac{\tau}{w_2-w_1}(1-r).
$$

This low-loss approximation requires same input, known widths, a validated
nonparalyzable response and a relevant path dead interval tau. A historical
60/50 coefficient encodes assumptions, not a universal calibration. Updating
gates and analog recovery can invalidate it. No PRE model is fitted here.

Even identical plane marginal availabilities need not imply identical trigger
survival. A toy independent3-of-4 trigger has product(Li) plus the sum of each
single-failed-plane combination. Four independent Li=0.9 give0.9477; perfectly
shared availability0.9 gives 0.9. The toy illustrates correlation sensitivity,
not a measured HMS correction.

For stationary Poisson input with fixed nonparalyzable tau, event-spacing
density is R exp[-R(delta-t minus tau)] above tau, and L=1/(1+R tau).
A minimum recorded spacing alone does not establish total busy/readout loss.
Buffering, bursts, trigger mixing, periodic input and changing rates break
the simple model. These teaching models remain separate from observations.

## 15. Comparisons and artificial normalization structure

The original user-supplied NewGen multipanel is retained unchanged; the saved
CSV supplies exact numerical comparisons. Do not read precision values from
the plotted pixels. All 11 exact saved-file controls reproduce original
NewGen counts. Other saved comparisons mix changed file/replay coverage with
method changes. Earlier claims of 42 matching controls belong to the old
168-segment snapshot and are not the full-cache control count.

![Original applied NewGen multipanel](figures/applied_NewGen.png)

![All 43 current diagnostic ratios](figures/refresh_livetime_comparison.png)

![TI6 cohort](figures/cohort_ti6.png)

![TI4 early unprescaled cohort](figures/cohort_ti4_early.png)

![TI4 late unprescaled cohort](figures/cohort_ti4_late.png)

![TI4 prescaled cohort](figures/cohort_ti4_prescaled.png)

These ratio plots retain4303/4305 flags and values above one. No curve is a
truth reference. The proposed raw-subtracted ratio uses Eraw; the legacy
run_summary.csv CLT_physics column uses Etight. Those cannot be mixed by
name alone. The authoritative current publication table is `production_summary.csv`.

With fixed yield and charge, replacing Lold by Lnew changes normalization by
Lold/Lnew minus1. On the same refreshed files, old to matched broad raw spans
-0.7666206% in 4493 to+0.6910719% in 4308. This is a sensitivity illustration,
not measured cross-section bias or a recommended replacement. Earlier tight
comparisons have slightly different extrema and are retained as historical.

![File changes versus method changes](figures/refresh_normalization_sensitivity.png)

![Rate associations by trigger](figures/rate_comparison_by_trigger.png)

The rate comparison uses S divided by selected interval exposure for p1 runs.
It is descriptive only; invalid scaler exposure can distort both axes. No
slope is fitted or adopted. Zhang's July 2024 source shows a KinC_60_4b row
with 19 runs, maximum TRIG1 rate1200 kHz, and slope entry-2.13 under the literal
header10^-5 %/kHz. Those original units and physical cause require verification;
no conversion or correction is inferred from that historical table.

Artificial structures can arise from unmatched beam trips/endpoints, different
replay coverage, timing acceptance, prescale sampling, counter loss, unmatched
yield/charge cuts, double-counted stages or applying a TI6-only trigger factor
to TI4. Such errors can track run, current or trigger changes rather than
appearing as a common scale. A smooth curve is not proof of validity.

## 16. Validation, unresolved decisions and reproduction

The frozen analysis verified all 187 input size/mtimes, representative independent
ROOT/uproot values, optimized-reader agreement, selected trigger bits and
the shared-count identity. The raw 4305 comparison is exhaustive for the
specified physics fields and relevant scaler channels. Prior source snapshots
and archival checksums protect provenance. No full ROOT rerun was needed
to produce this documentation edition.

Still required before a production prescription:

1. Identify the sent EDTM tap and detailed HMS discriminator/EL-REAL pulser
   fanout, including whether the selected coincidence reliably receives both
   pulses. Preserve timing and causality conditions rather than forcing unity.
2. Locate actual ROC5 readout/accumulator/FIFO code and run-period configuration
   to diagnose4303/4305, or design a justified matched yield/charge exclusion.
3. Establish effective whole-DAQ master/busy/prescale behavior as needed for
   any model, rather than substituting slave buffer settings for the whole DAQ.
4. Match the intended physics, good-helicity and charge exposure, and separately
   establish any NPS trigger/readout/reconstruction factors for the appropriate
   trigger cohort. Do not use CLT_physics as the target.
5. Quantify statistical dependence, timing/boundary/current sensitivity and
   remaining path/systematic uncertainty with explicit assumptions.

For document editing, change presentation.tex or report.tex/sections/report_body.tex
and run build.sh. These are standalone editable LaTeX sources; the metadata
and plot generator need not run for prose edits. REPORT.md is the parallel
human-readable source. generate_report.py regenerates the article body from
Markdown and tables, so run it only when deliberately replacing LaTeX prose
edits; it is not called by build.sh.

For plot regeneration from the frozen evidence, run make_comparisons.py.
The older diagnostic plots are copied with their generating scripts in
evidence/full_cache. The historical_168_segments figure directory is explicitly
historical and is not used as current numerical evidence. The accompanying
plot catalog identifies every copied/new figure and its scope.

For exact ROOT reproduction, use a fresh writable copy of evidence/full_cache,
the archived manifest and identical cached inputs. Load the HallC environment
before ROOT/hcana. The original export/analyze/timestamp/audit/acceptance/clock/
phase/final-consistency scripts are retained. Do not use live inventory or
dispatch scripts for an exact reproduction, and do not stage more raw data.
The raw extractor requires the already cached4305 file and the specified
column exports for full validation. Original commands and caveats remain in
the evidence reports; local absolute-path checks can fail after relocation
without implying numerical disagreement.

## 17. References and complete supporting record

The reference catalog following the main report lists every source family
used, its URL, scientific role and limitation. Original source PDFs/HTML,
searchable text and retrieval/checksum manifests are packaged in references.
The CODA source snapshots, source-code searches, supplied logbooks and dated
routing excerpts are in evidence. Missing or unverified source access is
identified explicitly; no login-protected page is claimed retrieved when
only the user text was available.

Appendices include all 43 production counts and ratio comparisons and the complete
reference catalog. The detailed research log is supplied as a companion Markdown file. Historical reports,
all current per-interval/per-run outputs, raw scalar arrays, analysis scripts
and original references are supplied as editable companion files rather than
thousands of printed rows. This preserves detail without making the presentation
unreadable. The presentation has its own main discussion and technical backup.
