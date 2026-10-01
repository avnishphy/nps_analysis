# NPSlib source clues for EDTM / ROC investigation

UTC2026-09-10. Read-only source/history and archived configuration-text audit.
No ROOT replay, raw staging, correction change or LaTeX edit.

## Scope and provenance

`/group/nps/apps/NPSlib` and `/u/group/nps/apps/NPSlib` are the same filesystem
directory (`os.path.samefile` checked). This focused review goes beyond the
earlier apps search by inspecting timing, VTP and config-parser implementations.
Versions0269b41,bf8ec57,fddb8e9 were searched. `evidence/commands.json` records
exact commands, times and exit codes; `evidence/sources.json` hashes20 source
snapshots. Relevant source/README status checks are clean in all three trees.
The newest tree's HEAD isfddb8e9d5355212510b8e2f5b5e2fffc8ce0f3bd,
dated2024-02-09. This does not prove it was loaded in the run or pass2 replay.
NPSlib's README describes a reconstruction add-on for hcana. Its modules here
decode recorded data, not a retrieved ROC5 hardware readout list.

## 1. EDTM signal interpretation in coincidence source

`src/THcNPSCoinTime.h:78` in all three versions describes T1 as
`NPS VTP||EDTM`, and T3 as `HMS 3/4`. Git blame dates the line to2023-09-12,
commit34c3ae03. This is a pre-run-period software comment, not a dated wiring
measurement. It supports investigating an OR of the VTP and pulser at that
named timing signal. It cannot establish the physical injection point relative
to all losses or certify the selected trigger scaler's meaning.

`THcNPSCoinTime.cxx:167--169` gets timing entries0,1,3 from the configured
`THcTrigDet`; comments identify pTRIG1_ROC1 (NPS VTP), pTRIG3_ROC1 (HMS3/4),
pTRIG4_ROC1 (HMS ELREAL). These are timing-channel descriptions. Mapping them
to actual TDC channels and TS input/scaler channels requires the replay map
and run-period hardware configuration. Do not equate the suffix or array
index with a verified trigger-input number.

Conditional physics implication: if the EDTM is injected only at the indicated
OR after VTP trigger formation, it would bypass earlier VTP formation losses.
That is a question to resolve from routing, not an established total-LT model.
This comment does not define PRE40/100/150/200 connections.

## 2. Event137 parser is not a complete configuration archive

`THcNPSConfigEvtHandler` was added2024-01-11 (history saved). It defaults to
event137 and reads `DAQInfoExtra::strings`. Its comment says separate parsing
is needed because VTP keys are missing from `THaRunBase::GetDAQInfo`.
The default exported strings are FADC sparsification and three VTP
cluster/readout/pair thresholds; additional keys can be requested.

The header declares `std::map<string,string> fNPSConfigData`. Parsing tokenizes
each line into first-token key and remaining-token value, then uses `emplace`
without replacing an existing entry. Thus **one handler instance retains the
first occurrence of a key**, including indexed rows sharing the same key.
No per-crate or index-qualified key is constructed in this code. We have not
established the exact handler lifecycle in pass2; this is inspected-source
behavior, not a reproduced pass2 metadata result or a production bug claim.

The already archived4305 raw text contains five VTP configs, each with32
`VTP_NPS_TRIG_PRESCALE` lines. It also contains latency values2560,2556,2556,
2536,2536 in raw-record order43,44,45,46,49. A first-key model over the archived
text retains only prescale value `0 1` and latency `2560`. The other entries
remain in our original raw-text archive. `evidence/config_key_audit.json`
records all line numbers/values and hashes for the ten input config records.
The model demonstrates information loss from unique-key insertion; it does
not execute the ROOT handler. No new raw read was needed.

Decision: retain indexed rows and raw-record/crate provenance; do not use a
single `GetInfo` value as a full multi-crate/per-channel configuration. The
VTP internal prescale table still is not the Hall C trigger prescale table.

## 3. VTP decision bits are another conditional diagnostic

`VTPModule.cxx:187--217` decodes bits0--5 of the VTP decision pattern.
`THcNPSCalorimeter.cxx:399--406` gives these output descriptions:

| Bit | Source label |
|---|---|
| 0 | Cluster singles |
| 1 | Cosmic scintillator |
| 2 | Cosmic column |
| 3 | Cluster pair, different crate |
| 4 | Cluster pair, same crate |
| 5 | VLD |

VLD is not labeled EDTM. These VTP bits are distinct from the six-bit global
TI trigger mask already compared between ROOT and raw. Source descriptions
do not identify which signal caused the final recorded global trigger.
The calorimeter code keeps VTP crate tags and checks event-number agreement;
any future branch-based study must verify that mapping and its error flags.

Additional caution before using VTP decision timing: line201 masks with
`0x7F` (seven bits), while its comment says eleven bits. No matching firmware
format was validated here; do not silently change the mask or assume a timing
unit from this comment. This does not affect the previously used TI event
timestamps or EDTM TDC analysis, which did not use this VTP decision-time field.

## Outcome and remaining work

Useful new evidence: pre-run-period EDTM OR comment, VTP output labels and a
specific configuration-parser limitation. Still no ROC5 readout source, scaler
FIFO/accumulation implementation, verified PRE mapping or exact EDTM path.
CLT_physics remains a proposed definition only. No new denominator is adopted.

The next targeted source check is the run-period trigger detector map and
setup associated with these timing labels, alongside the saved coin_sparse
ROC5 rol1/rol2 source and master/TS configuration. VTP decision-bit comparisons
could test conditional correlations once fields/mapping are validated; they
cannot recover unrecorded triggers or replace missing scaler exposure alone.

## Reproduce the completed checks

From this archive, copy the script to a fresh writable directory:

```bash
mkdir -p /tmp/npslib_roc_audit_repeat
cp collect.py /tmp/npslib_roc_audit_repeat/
python3 /tmp/npslib_roc_audit_repeat/collect.py
sha256sum -c SHA256SUMS --quiet
```

Requires Python3, rg, git and access to the named NPSlib and frozen raw-text
paths. The script reads those sources and writes only beside itself. Expected:
same-directoryTrue,20 source snapshots,10 config records,32 prescale rows in
each of five VTP records. Hash verification checks this frozen archive and
is silent on success. Exact Git commands and negative-search scope are saved;
no implication that every historical object or external filesystem was searched.
