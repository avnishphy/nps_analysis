"""Apply user-confirmed working routing to frozen run metadata; no event reread."""
import csv,hashlib,json
from pathlib import Path
P=Path(__file__).resolve().parent
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
A=M/'livetime_refresh_20260909'
inputs=[A/'production_summary.csv',A/'snapshot/nps_dvcs_all_kins_main.csv']
summary=list(csv.DictReader(inputs[0].open()))
metadata={}
for r in csv.DictReader(inputs[1].open()):
 d={k.strip():v.strip() for k,v in r.items() if k is not None}
 # Aggregate metadata rows such as 2499-2512 are not individual runs.
 # Every requested production run must still have one exact match below.
 if d['run_number'].isdigit():metadata.setdefault(int(d['run_number']),[]).append(d)
rows=[]
for r in summary:
 n=int(r['run']);t=int(r['trigger']);p=int(r['p']);matches=metadata[n]
 assert len(matches)==1,(n,len(matches))
 d=matches[0]
 assert float(d[f'PS_{t}'])==p,(n,t,p,d[f'PS_{t}'])
 enabled=[i for i in range(1,7) if float(d[f'PS_{i}'])>0]
 assert enabled==[t],(n,enabled,t)
 rows.append(dict(run=n,date=d['day'],trigger=t,p=p,prescale_token=d['prescale'],metadata_coin_status=d['coin_status'],events=int(r['events']),
  working_trigger_description={4:'HMS EL-REAL singles',6:'NPS cluster-or-delayed-EDTM AND HMS EL-REAL coincidence'}[t]))
with (P/'run_trigger_interpretation.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
groups={str(t):dict(runs=[r['run'] for r in rows if r['trigger']==t],events=sum(r['events'] for r in rows if r['trigger']==t)) for t in [4,6]}
assert len(rows)==43 and sum(r['events'] for r in rows)==77838384
result=dict(inputs=[dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in inputs],groups=groups,all_single_enabled=True,
 applicability_basis='User confirmation; logbook4171556 dated2023-08-29 and4170042 dated2023-08-23. Physical routing not independently measured here.',new_events_read=0)
(P/'mapping_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS:43 runs; single enabled trigger and p agree with frozen metadata; TI6=5 runs/2706910 events; TI4=38 runs/75131474 events.')
