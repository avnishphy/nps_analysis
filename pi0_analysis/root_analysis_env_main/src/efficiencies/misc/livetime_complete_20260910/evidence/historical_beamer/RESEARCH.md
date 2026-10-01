# Livetime definitions and new source evidence

2026-09-09. This note extends the frozen KinC_x60_4b LH2 calculation. It does
not certify a replacement production correction. PRE attribution remains
unresolved. The editable technical presentation is `livetime.tex`.

## Findings that change how we present the result

1. The matched +/-2 ns ratio must remain a **timing-selection diagnostic**.
   A narrow peak and a small width sensitivity do not measure the efficiency
   for accepting genuine EDTM pulses. The next check should classify rejected
   candidates using raw versus reference-subtracted time, available multiple
   hits, trigger decisions and interval/rate dependence. The matched raw-window
   estimate is therefore retained alongside the tight estimate in every
   production numeric table. No production timing selection is changed.
2. FADC/VTP beginning-of-run EVIO configurations are a concrete evidence lead.
   They may establish settings even when logbook access is unavailable, but
   they do not necessarily encode PRE cabling, analog module modes or the TI
   prescale/readback configuration. Those require their own records.
3. Historical NPS scaler repairs need explicit validation before transfer.
   The local later analyzer `sources/local/DT_Analyzer_new.C` adds
   `EDTM_Scaler_correction` to its denominator near lines 575--588. Its current
   bytes are archived; they are not identified as the executable that produced
   every historical plot. Our existing diagnostic has no such inferred repair.

## New primary sources and their limits

| ID | Source and location | Relevant evidence / scope |
|---|---|---|
| N01 | [Murphy, EDTM Study Report, 2022-06-28](https://redmine.jlab.org/attachments/download/1499/EDTM_Study_Report.pdf), sections 1--4 | PionLT observed valid corrected-time tails and trigger-priority/prescale competition. This cautions against assuming universal timing cuts or a single-factor multi-trigger correction; it does not prove either effect in NPS. |
| N02 | [Zhang, NPS Deadtime and Efficiency, 2024-07-18](https://indico.jlab.org/event/866/contributions/14918/subcontributions/283/attachments/11529/17843/Deadtime%20and%20Efficiency.pdf), slides 3, 5--13 | Period trigger definitions, earlier LT discrepancies and an LH2 rate correlation. Slide 13's KinC_60_4b row has 19 runs, maximum TRIG1 rate 1200 kHz, slope entry -2.13; the literal header is `10^-5 %/kHz`. Original fit units and physical cause require verification. No slope correction adopted. |
| N03 | [Raydo, NPS Trigger/DAQ, 2024-07-17](https://indico.jlab.org/event/866/contributions/14914/subcontributions/234/attachments/11505/17805/NPS_July17_2024_Raydo.pdf), slides 2, 7--8, 11--12 | Separate trigger/readout thresholds; EVIO configuration records and VTP decisions support future run-specific checks. Waveform/trigger agreement tests the observed sample, not automatically the events that failed to trigger. |
| N04 | [Zhang, 2024-02-19 report](https://redmine.jlab.org/attachments/download/2343/Dead%20time%20analysis_2024-02-19.pdf) | Historical raw-TDC/scaler investigation; examples are not measurements of all runs in our snapshot. |
| N05 | [Zhang, 2024-03-03 update](https://redmine.jlab.org/attachments/download/2368/Updates_2024-03-03.pdf), slide 1 | Earlier workflow: automatic raw-time cut, 2-uA threshold, pulser-scaler repair and terminal-readout removal. These choices are not silently imported. |
| N06 | [NPS issue 837](https://redmine.jlab.org/issues/837) | Public source of the chronological attachments above. A historical resolved claim does not certify our different replay sample. |
| N07 | [NPS issue 866](https://redmine.jlab.org/issues/866), [DAQ issue 836](https://redmine.jlab.org/issues/836) | Electronics-deadtime and event-blocking task records; neither supplies run 4398's hardware state. |
| N08 | [Hall C Trigger History](https://hallcweb.jlab.org/wiki/index.php/Trigger_History) | Links [entry 3532998](https://logbooks.jlab.org/entry/3532998), the 2018-02-16 electronic-deadtime cabling lead. Logbook contents were not read; historical cabling cannot establish 2024 PRE attribution. |
| N09 | [Hall C EDTM Pulser](https://hallcweb.jlab.org/wiki/index.php?title=Hall_C_EDTM_Pulser) | Points to [entry 3893007](https://logbooks.jlab.org/entry/3893007) for the cdaqpi1 pulser change. It does not establish any run's actual rate/phase. |

Downloaded bytes, retrieval timestamps, redirects and SHA256 values are in
`sources/web_manifest.json`. Public HTML was retrieved without credentials.
Web-tool fetch failures were followed by ordinary public downloads, not login
bypass. The user-protected logbook was not accessed. Prior primary documents
used in the deck were copied from the approved archive and are enumerated in
`sources/local_manifest.json`; the full earlier catalog remains in the sibling
`SOURCES_LIVETIME_4398.md`.

## Exact operational definitions

Fix a reference population G that would satisfy the chosen detector/trigger
criteria absent rate losses. Let B mean the selected trigger input forms,
P that it passes intentional prescaling, and R that it is recorded.

```text
L_elec = Pr(B | G)
L_DAQ  = Pr(R | P,B,G)
Pr(R | G) = L_elec * Pr(P | B,G) * L_DAQ
```

If representative prescaling gives Pr(P|B,G)=1/p, the prescale-removed total
livetime is p*Pr(R|G)=L_elec*L_DAQ. This is a conditional chain, not an
independence assumption. Detector inefficiency, thresholds and tracking/PID
efficiencies require separately stated reference populations. A branch called
Computer_live_time does not define the endpoints of a probability.

On identical selected intervals define N=recorded, E=accepted EDTM subset,
S=selected trigger input, D=sent EDTM, A=L1 accepts:

```text
L_EDTM_hat = p E / D
L_all_hat  = p N / S
L_phys_hat = p (N-E) / (S-D)
L_L1_hat   = p A / S
```

The physical labels are conditional: one sent pulse must contribute to the
chosen input, E must identify that accepted subset, and prescale sampling must
be representative. Accepted E belongs in the numerator subtraction; sent D
belongs in the denominator subtraction. The estimators share counts and obey
the exact mixture identity; they are not independent tests.

For availability a(t), time livetime is integral(a)/T, while physics livetime
weights by the eligible input rate r(t): integral(r*a)/integral(r). A periodic
pulser samples availability at its pulse times; a common current cut alone
does not make that measure equal to the physics-weighted one. If yield per
charge is stable, an effective factor is integral(I*L)/integral(I). Matching
good-helicity selection, yield and charge is still required before production
use. Multiplying EDTM and a DAQ factor can double count common losses.

## PRE and model assumptions

PRE40/100/150/200 have counts but unknown run-period trigger attribution.
The archived scaler map establishes addresses, not the upstream signal source.
Correlation is insufficient to assert HMS hEL-REAL, NPS, or a particular
coincidence population. No PRE correction has been applied.

The deck derives the width-ratio formula for an explicitly hypothetical
nonparalyzable Poisson monitor, rather than declaring that model describes the
hardware. It also contrasts leading-edge counts with updating logic acceptance
and shows why a 3-of-4 trigger can have different losses from its planes.
Those are teaching examples; no plane-rate model is fitted to these data.
Likewise, a shifted-exponential event-spacing fit describes a fixed
nonparalyzable model, not arbitrary buffered TI/readout operation. The
unbuffered periodic-pulser literature does not establish the actual TI mode.

## Cache and continuation

All ratios retain the initial 168-segment audit. The separate
`cache_availability_only.json` contains a later filename/size snapshot during
user-initiated staging; it is not a ROOT payload read or a completeness test.
Several formerly absent segments have appeared. Some selected updated replays
still have gaps even when fuller production alternatives exist. A new run must
freeze and label its per-run replay choice; never mix updated and production
copies or silently substitute them in a same-cache method comparison.

The next calculation must preserve original and refreshed controls, verify
selected file stability, repeat the same interval/counter audits, retain the
matched raw and timing variants, and report coverage gained. No jcache command
was issued. The current production efficiency code, CSVs and yields were not
changed by this presentation/research extension.
