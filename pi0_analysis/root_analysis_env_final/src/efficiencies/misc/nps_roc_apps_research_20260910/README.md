# NPS apps / ROC source search evidence

UTC2026-09-10. See the living [research log](../LIVETIME_RESEARCH_LOG.md),
especially the dated apps-search entry. This directory freezes that entry's
inputs and exact command outputs; its copy of the research log is a dated
snapshot. Later changes belong in the living parent document.

The search found replay/monitoring software, an offline scaler decoder and a
modulefile reference to `/cdaqfs1/apps`. It did not locate the online ROC5
readout list. ROC5 and npsvme5 are distinct log components. ROC5 end counts
agree with recorded totals but do not certify scaler-channel integrity.

## Reproduce

Completed commands and their outputs are in `evidence/commands.json` and
`evidence/*.stdout.txt` / `*.stderr.txt`. To repeat discovery into a fresh
writable directory without editing this archive:

```bash
mkdir -p /tmp/nps_roc_apps_audit_repeat
cp collect_evidence.py /tmp/nps_roc_apps_audit_repeat/
python3 /tmp/nps_roc_apps_audit_repeat/collect_evidence.py
```

Requires Python3, rg, git and read access to the documented JLab filesystem
paths. No ROOT, network, DAQ login or staging needed. The script reads the
previous frozen user logs at their explicit path. It writes only beside itself.
Expected:2562 inventoried paths for this snapshot; precise ROC-text query exit1
with empty stderr (no match); three decoder status checks exit0 with empty
output (clean for that file); online location ls exit2 with ENOENT messages.
Future inventory changes are possible. Broad-query compiler multilib matches
are false positives; do not interpret exit0 as an online-driver discovery.

`source_manifest.json` includes original and resolved paths, byte sizes and
SHA256. All three THcScalerEvtHandler source snapshots have the same hash;
their repository revisions differ. These are inspected source files, not an
assertion of exact run/replay build provenance.

Verify the completed archive from its directory:

```bash
sha256sum -c SHA256SUMS --quiet
```

No output with exit0 means all archived file hashes match. No production
analysis, raw staging, correction or LaTeX changes were made in this search.
