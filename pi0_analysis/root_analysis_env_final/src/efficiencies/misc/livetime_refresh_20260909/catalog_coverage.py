"""Recreate the read-only catalog listing; no payload reads or staging."""
from pathlib import Path
import json
P=Path(__file__).resolve().parent
inv=json.loads((P/'inventory.json').read_text())['runs'];rows=[]
for r in inv:
    source=r['source'] if r['source']!='none' else 'updated'
    d=Path('/mss/hallc/c-nps/analysis/pass2/replays')/source
    segs=sorted(int(f.name.split('_')[4]) for f in d.glob(f"nps_hms_coin_{r['run']}_*_1_-1.root"))
    rows.append(dict(run=r['run'],source=source,catalog_segments=segs,cached_segments=[f['segment'] for f in r['files']],missing_catalog_segments=sorted(set(segs)-{f['segment'] for f in r['files']})))
(P/'catalog_coverage.json').write_text(json.dumps(dict(method='Read-only directory listing of /mss catalog; no file payload reads and no staging commands',runs=rows),indent=2))
