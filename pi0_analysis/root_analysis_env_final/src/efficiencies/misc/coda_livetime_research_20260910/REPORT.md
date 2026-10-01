# CODA DAQ / ROC evidence for the livetime investigation

Updated UTC 2026-09-10. Official CODA sources requested by the user were read,
downloaded and hashed. **No total-livetime prescription established.**
CLT_physics remains a proposed definition only. No production correction,
PRE attribution, scaler repair, extra staging or LaTeX change is made.

## What CODA documentation establishes

CODA connects hardware readout controllers (ROCs), event building and output.
Experiment-specific hardware readout is configured through readout lists.
The official CODA 3 example stores the configuration under `COOL_HOME` and gives
each ROC `rol1` and `rol2` shared-library paths. Thus the concrete next artifact
is the saved **coin_sparse ROC5 configuration**, including those paths and the
matching source/build, not another ROOT replay. Example paths on the website
are illustrative; they are not NPS paths.
[CODA configuration example](https://coda.jlab.org/drupal/content/example-single-crate-configuration).
The [3.10 walkthrough](https://coda.jlab.org/drupal/system/files/coda/3.10/walkthrough/index.html)
also identifies ROL1 as the compiled readout-list setting. Neither page supplies
the February 2024 NPS readout list or wiring.

The [CODA driver index](https://coda.jlab.org/drupal/content/vme-module-drivers)
links TI documentation labeled **3v11.3**, matching firmware `11.3 (tip113.svf)`
printed in both user logs. This is version-relevant documentation, not proof
of the actual compiled library revision used online.
The [TI header](https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/tiLib_8h_source.html)
defines the broadcast buffer mask, local/broadcast selector and busy enable.
The [configuration API](https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/group__Config.html)
distinguishes buffer-threshold busy from the external busy-source mask.

The [TI hardware manual](https://coda.jlab.org/drupal/system/files/pdfs/HardwareManual/TI/TI.pdf)
downloaded here was updated **August 27, 2026**, later than the runs. Its register
discussion describes live/busy timers and a required scaler latch. Do not
silently apply its later features to the 2024 firmware. The saved
[TS manual](https://coda.jlab.org/drupal/system/files/pdfs/HardwareManual/TS/TS.pdf)
is dated March 16, 2017 and describes a Version 4 TS. It distinguishes trigger
inputs, prescaling/trigger decisions, distributed accepts and returned busy.
The NPS master model, firmware and settings have not been established by these
slave logs; the TS manual is architectural context only.

## Verified log observations

`check_logs.py` checks all 20 initial/end TI register pairs in runs 4303/4305.
Exact source hashes and values are in `log_checks.json` and
`ti_log_registers.csv`. The previously archived user logs remain the source.

| Field | Logged value | Interpretation using TI header |
|---|---|---|
| dataFormat | 0x05050006 | Bits 31:24 give broadcast buffer level 5 |
| vmeControl | 0x00800010 | Bit 22 clear: broadcast selected; bit 23 set: buffer-threshold busy enabled |
| blockBuffer | 0x00000001 | Local setting; not the selected threshold |
| busy | 0x00000003 | Switch-slot A/B busy sources enabled |
| Block Level | 1 | One event per configured block |

Inference: these settings support buffered operation in the reported slaves.
They do not prove five events can always be accepted by the complete DAQ:
other participating crates, master settings and busy paths can constrain it.
The busy enable mask does not measure how long any source was busy.

All five end-status blocks in each run have equal Readout/Ack/L1A counts:
**3,310,432 for 4303; 560,647 for 4305**, matching the recorded ROOT totals
(and the independently read raw 4305 total). This establishes accepted/readout
count consistency at those endpoints, not the number of triggers lost before
acceptance. Do not turn this agreement into a livetime of one.

Initial TI L1A counts are zero but initial timers are nonzero. For the
unambiguously paired boardID `0x7101150e`:

| Run | Initial live / busy hex | End live / busy hex |
|---|---|---|
| 4303 | 1d39dad5 / 260a | 202dbeaf / 1afd |
| 4305 | 0a1a58d9 / 0670 | 0a3e6c71 / 029a |

The busy value decreases. These are not demonstrated common, cumulative
run-boundary snapshots. Reset, latch and rollover history must be established
before subtraction or endpoint ratios. Two interleaved timer chunks in4305
are deliberately left unassigned by the script; nearby log prefixes do not
reliably identify the owner of mixed output.

The [TI status API](https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/group__Status.html)
documents timer units of 7.68 microseconds, `tiLive(sflag>0)` as integrated and
its return scale as tenths of a percent. It also documents `tiReadScalers`
options without latch, with latch, and with latch/reset. The
[timer-latch API](https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/group__MasterConfig.html)
requires latching before timer getters. The available API descriptions do not
establish whether the actual `tiStatus` call latched or which reset calls the
NPS readout list executed. Consequently no TI timer correction is calculated.

## Definitions and consequences for this analysis

These definitions are mathematical distinctions, not calibrated estimates:

- **Time availability at a specified busy decision:** live-duration divided
  by live-plus-busy duration over a verified common exposure. The decision
  point and included busy sources must be stated. Time availability need not
  equal event survival when event rate and busy are correlated.
- **Trigger acceptance:** accepted requests divided by eligible requests at
  explicitly defined points before/after prescaling and busy. Overlapping
  triggers, multiplicity and prescale selection affect count correspondence.
- **Readout survival after acceptance:** recorded accepted events divided by
  issued accepts over the same exposure, accounting for blocks/event types.
- **EDTM survival:** recorded injected pulses divided by sent injected pulses,
  adjusted for a verified sampling/prescale rule. It probes the actual pulse
  route. A timing-associated TDC hit need not have caused the recorded trigger.
- **Total physics survival:** probability a physics opportunity survives all
  stages intended to be corrected. A product requires nonoverlapping,
  conditional stages for the same population; arbitrary livetime products
  can count one loss twice.

For 4305, the earlier raw investigation already found 1047 recorded events
and 80 EDTM tags while the suspect interval's mapped clock advanced only
5481 ticks, with L1/TRIG4 increments2 and EDTM0. Other scaler channels advanced.
The endpoint TI agreement adds independent accepted-count evidence. Inference:
ordinary upstream deadtime alone cannot explain an interval containing those
recorded accepts while its mapped counters scarcely advance. Hardware inhibit
of scaler groups, lost FIFO samples or ROC accumulation/readout behavior remain
possibilities. These observations do not distinguish those causes.

No CODA page inspected identifies NPS EDTM injection points or PRE40/100/150/200
connections. No alternative denominator is certified. Changing denominator
exposure, tag acceptance or correction factors inconsistently can generate
run-dependent normalization structure; agreement with CLT_physics is not a
validity test. Invalid-interval exclusion must match physics yield and charge.

## Evidence needed next

1. Saved run-period `coin_sparse` configuration: ROC5 rol1/rol2, shared-library
   build/source and scaler driver. Inspect FIFO draining, accumulation, latch,
   clear/reset/inhibit, overflow/error handling and sync/end bank handling.
2. Actual master/TS configuration: model/firmware, trigger lookup/rules,
   prescales, active busy inputs, buffer settings and timer latch/reset calls.
3. EDTM fanout/injection and discriminator/trigger routing; exact PRE mappings.

These are existing configuration/source artifacts to retrieve read-only. No
DAQ control operations or messages to operators were performed. Broad local
source search was stopped when unproductive; the exact readout list remains
unlocated. The public ROC landing page alone contains no implementation detail.

## Reproduce the completed checks

From this report's directory, Python3 only; no ROOT setup or staging needed:

```bash
python3 check_logs.py
sha256sum -c SHA256SUMS --quiet
rg -n 'TI_DATAFORMAT_BCAST_BUFFERLEVEL_MASK|TI_VMECONTROL_USE_LOCAL_BUFFERLEVEL|TI_VMECONTROL_BUSY_ON_BUFFERLEVEL' tiLib_8h_source.html.txt
```

Expected: 20 verified register pairs and ten end count triplets matching the
recorded totals; hash verification has no output on success. `check_logs.py`
reads the adjacent frozen analysis archive and rewrites only its two small
diagnostic outputs. `fetch_sources.py URL ...` reproduced public downloads;
the manifest records exact URLs, UTC timestamps, bytes, SHA256 and failures.
The guessed `tiLib_8c.html` returned404; function documentation and header were
available, not the matching implementation. No dependency was installed.

CODA's [downloads page](https://coda.jlab.org/drupal/content/downloads) identifies
`/group/da/distribution/coda`; its local LibraryManual/tiLib directory is readable.
`local_sources.json` records the local MasterConfig copy and equality checks
between downloaded documentation and local published copies. Source hashes
freeze what was inspected; neither a current website nor filesystem mtime
establishes the library loaded in February2024.
