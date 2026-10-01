"""Independent uproot read validates representative ROOT-exported values."""
from pathlib import Path
import json
import numpy as np
import uproot
from analyze import P,read
inv={r['run']:r for r in json.loads((P/'inventory.json').read_text())['runs']}
checks=[]
for run,seg,tree in [(4259,0,'T'),(4305,0,'T'),(4350,3,'TSH'),(4350,4,'TSH'),(4350,5,'TSH')]:
    path=next(f['path'] for f in inv[run]['files'] if f['segment']==seg)
    exported=read(run,seg,tree)
    keys=['g.evnum','g.evtime','g.trigbits','T.hms.hEDTM_tdcTimeRaw','T.hms.hEDTM_tdcTime'] if tree=='T' else list(exported)
    with uproot.open(path) as f:a=f[tree].arrays(keys,library='np')
    for k in keys:
        if not np.array_equal(a[k],exported[k],equal_nan=True):raise RuntimeError('Reader mismatch '+str((run,seg,tree,k)))
    checks.append(dict(run=run,segment=seg,tree=tree,rows=len(a[keys[0]]),columns=len(keys),all_values_identical=True))
    if run==4305:
        m=(a['g.evnum']>=105080)&(a['g.evnum']<106127)
        E=np.sum(m&(a['T.hms.hEDTM_tdcTimeRaw']>1)&(abs(a['T.hms.hEDTM_tdcTime']-245.11508)<=2))
        assert m.sum()==1047 and E==80
        assert np.all(np.diff(a['g.evtime'][m])>0)
        checks[-1].update(suspect_interval_N=int(m.sum()),suspect_interval_E=int(E),event_timestamps_strictly_increasing=True)
(P/'independent_reader_checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
