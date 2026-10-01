# Raw 4305 and user-provided logbook evidence

Continuation of [full-cache ROOT investigation](REPORT.md), UTC 2026-09-10.
**CLT_physics remains a proposed definition, not a baseline or true livetime.**
No final EDTM/total-LT prescription has been adopted.

## Authorized input and integrity

User authorized staging exactly one file. Submitted once:

```bash
env -u DEBUG jcache get /mss/hallc/c-nps/raw/nps_coin_4305.dat.0
```

Request 100232289. `DEBUG=release` otherwise breaks the Java wrapper.
13,341,141,036 bytes; CRC32 `53e8a0c0`, Adler32 `756ab5b0`, both match the
archived tape stub. File fully read to EOF. No other raw file staged or read.
The earlier no-cache inventory is an initial snapshot, not the final state.

User attachments for runs4303/4305 are copied byte-for-byte under `logbook/`;
`manifest.json` retains source paths, line references and SHA256 hashes.
The raw file was chosen because it contains both run-start data and the
4305 scaler anomaly. Run4303 raw data were NOT examined.

## What the raw file establishes

- 561,307 raw records, 560,647 physics events, 660 relevant scaler banks.
  Event number, event-clock timestamp and six-bit trigger mask agree exactly
  with ROOT for **every physics event**. There are no event-number gaps.
- All 660 raw scaler snapshots match ROOT's clock, EDTM, L1 and six HMS
  trigger branches. The additional final ROOT row duplicates final counter
  values but increments the event boundary from560647 to560648. This is
  bookkeeping, not extra scaler exposure; do not double count it.
- Raw records105247 and106294, scaler banks125 and126, delimit physics
  events105080 and106127. The same small increments seen in ROOT are present
  directly in the raw banks: clock5481 ticks=0.005481s, EDTM0, L1=2,
  configured HMS TRIG4=2, versus1,047 recorded events and80 EDTM tags.
  Event timestamps span about2.0015s.
- The fixed scaler-module header sequence is intact. Slots6--9 hodoscope
  channels continue advancing; clock/EDTM/trigger scaler counts in slots10--12
  have the anomalously small exposure. This rules out this analysis's ROOT
  read/export or a ROOT-only counter conversion as the source of the defect.
- Before the anomaly, raw L1 minus recorded event number is0 or1. Afterwards
  it is-1046 in524 snapshots and-1045 in11 snapshots, through run end.
  The deficit does not recover at the next cumulative snapshot.

This establishes a persistent deficit in the raw scaler record. It does not
distinguish scaler hardware inhibition, FIFO/accumulation/readout software
loss, or another upstream cause. A common delayed latch is also not an
automatic explanation for a nearly fixed missing event count through changing
rates. Need the actual ROC5 scaler readout/accumulation logic and configuration.
Do not infer or add missing counts without that evidence.

## Unassigned pulse-like channels: evidence, not a replacement

Slot7/channel31 and slot9/channel31 are identical in all660 snapshots.
Both equal slot12/channel14 (mapped EDTM) in the first125 snapshots, and
exceed it by exactly80 in every subsequent snapshot. The sole step occurs
across the anomalous interval. Multiple other unused hodoscope channels show
the same pulse-like increment pattern there.

The archived scaler map calls these two channels `Empty_12` and `Empty_24`.
Their behavior supports the localization of the scaler discrepancy, but
**count agreement does not establish wiring or authorize an EDTM denominator**.
This is the same attribution problem as PRE40/100/150/200. No replacement
denominator, inferred repair, or production interval exclusion was applied.

[Raw exposure/counter step plot](figures/raw4305_exposure_step.png).
Exact per-slot/channel increments, raw/ROOT comparisons and the unassigned
channel step are in `raw4305/validation.json`; full module arrays are in
`raw4305/scaler_modules.npz`, with a compact named-counter CSV alongside.

## What the logbook and embedded configurations establish

Both user logs explicitly identify their runs and end with event counts
matching ROOT (4305:560647; 4303:3310432). Both show `coin_sparse`, and ten
TI slave status blocks (prestart/end) with block level1 and broadcast block
buffer level5, busy enabled. For example4305 lines455--515 contain a status
block, and line4065 records ROC1's final event count.

This proves those slave settings, **not the effective buffering/deadtime of
the entire acquisition**. Master/TS and HMS busy behavior can still constrain
it. Neither log names EDTM/PRE routing or gives a master prescale readback.
Generic strings such as `ignore data errors = false`, `Bus Errors Disabled`
and error-status table headings are not observed error messages. No explicit
ERROR/FATAL message explaining the anomalous interval was found. Interleaved
ROC log output should not be assigned to a crate merely by nearest line.

The raw file contains ten text configuration banks: five FADC and five VTP,
preserved under `raw4305/config_text/`. Literal VTP fields include
TRIG_WIDTH20, latency2560/2556/2536, internal TRIG_PRESCALE entries and cluster
threshold1600/pair threshold800. These are module settings; **do not identify
the VTP internal prescale with the Hall C trigger prescale**, or infer physical
threshold/time units without the matching firmware/source definition.

There are no event-type120/133 prescale records in this raw file. Types135/136
contain uuencoded SHMS/HMS angle photographs, despite generic analyzer labels
suggesting detector/trigger files; they are not an EDTM trigger-logic diagram.
Thus absence of a decoded prescale readback in Run_Data is unsurprising here.
The embedded text and EPICS search did not identify EDTM or PRE wiring.

## Reproduce without another staging request

From a fresh writable copy of this archive with the already-cached raw file:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; hcana -l -b -q "raw_extract.C(0)"' > raw_extract.log 2>&1
python3 raw_validate.py
python3 raw_config_probe.py
python3 raw_figure.py
```

`raw_validate.py` additionally requires the ROOT scalar column cache described
in REPORT.md. The extraction uses EVIO's byte handling and the installed CODA3
trigger-bank decoder, without replaying detector reconstruction. Exported bank
words are host-endian uint32, not original-file byte order. Decoder sources
and map bytes are snapshotted with hashes in `raw4305/source_evidence/`.

Next evidence needed: ROC5's actual scaler readout/accumulation code and
run-period EDTM/trigger wiring/master settings. These results explain the
4305 invalid count ratios, but do not establish a total livetime or turn any
candidate definition into a reference.
