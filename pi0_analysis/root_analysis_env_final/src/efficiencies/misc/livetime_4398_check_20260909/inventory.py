"""Read-only run-4398 diagnostic; outputs only into this script's directory."""
from pathlib import Path
import csv, json, re
import numpy as np
import uproot

OUT = Path(__file__).resolve().parent
INPUT = '/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_4398_0_1_-1.root'
f = uproot.open(INPUT)
inventory = {'input': INPUT, 'size': Path(INPUT).stat().st_size,
             'uuid': str(f.file.uuid), 'keys': f.classnames(), 'trees': {}}
for key, cls in f.classnames().items():
    if cls == 'TTree':
        t = f[key]
        inventory['trees'][key] = {'entries': t.num_entries, 'branches': t.typenames()}
(OUT/'inventory.json').write_text(json.dumps(inventory, indent=2))
s = f['TSH'].arrays(library='np')
np.savez_compressed(OUT/'scalers.npz', **s)
h = f['TSHelH'].arrays(library='np')
np.savez_compressed(OUT/'helicity_scalers.npz', **h)
branches = [b for b in f['T'].keys() if re.fullmatch(
    r'g\..*|T\.hms\.(hEDTM|[hp]PRE\d+|hTRIG\d|npsTRIG\d)_(tdcTimeRaw|tdcTime|tdcMultiplicity)|H\.(BCM4A|1MHz|EDTM|hTRIG4|hL1ACCP)\.[^.]+', b)]
e = f['T'].arrays(branches, library='np')
np.savez_compressed(OUT/'events.npz', **e)
rows=[]
for name,a in s.items():
    if not re.search(r'PRE|TRIG|L1ACCP|EDTM|1MHz|BCM4A|EL_CLEAN|evcount|evNumber',name): continue
    rows.append(dict(branch=name,first=float(a[0]),last=float(a[-1]),min=float(np.min(a)),
                     max=float(np.max(a)),nonzero=int(np.count_nonzero(a)),
                     negative_steps=int(np.sum(np.diff(a)<0))))
with (OUT/'scaler_summary.csv').open('w') as out:
    w=csv.DictWriter(out,fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)
print('trees', {k:v['entries'] for k,v in inventory['trees'].items()})
for r in rows:
    if r['branch'].endswith('.scaler') or r['branch'] in ['evcount','evNumber','H.1MHz.scalerTime']:
        print(r)
for b in ['g.evnum','g.evtyp','g.trigbits','T.hms.hEDTM_tdcTimeRaw','T.hms.hTRIG4_tdcTimeRaw']:
    a=e[b]; u,n=np.unique(a,return_counts=True)
    print(b,'range',a.min(),a.max(),'nonzero',np.count_nonzero(a),'top',sorted(zip(n,u),reverse=True)[:8])
