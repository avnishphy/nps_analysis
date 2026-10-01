"""Compare scaler exposure with recorded-event timestamps; never repair counts.

Use only exact event-number matches at both scaler boundaries. Infer relative
clock frequency from the median of healthy positive intervals; short event
latency and boundary jitter remain in the residual. Large deficits identify
exposure disagreement, not its hardware or decoding cause.
"""
from pathlib import Path
import csv,json,numpy as np
from timestamp_probe import read
P=Path(__file__).resolve().parent
checks=[];bad=[]
for f in sorted((P/'results').glob('run*.json')):
    r=json.loads(f.read_text())
    if r['status']!='CALCULATED':continue
    run=r['run'];ii=list(csv.DictReader((P/f'results/run{run}_intervals.csv').open()))
    for c in r['components']:
        parts=[read(run,s,'T') for s in c['segments']]
        ev=np.concatenate([x['g.evnum'] for x in parts]);ts=np.concatenate([x['g.evtime'] for x in parts])
        if np.any(np.diff(ev)<=0):raise ValueError((run,c['component'],'unordered events'))
        rows=[x for x in ii if int(x['component'])==c['component'] and x['nominal']=='1']
        if not rows:continue
        b=np.array([[int(x['ev_start']),int(x['ev_end'])] for x in rows]);j=np.searchsorted(ev,b)
        inside=np.all(j<len(ev),axis=1);safe=np.minimum(j,len(ev)-1)
        exact=inside&np.all(ev[safe]==b,axis=1)
        dt=np.array([float(x['dt']) for x in rows]);ticks=np.diff(ts[safe],axis=1)[:,0]
        good=exact&(dt>1)&(ticks>0)
        freq=float(np.median(ticks[good]/dt[good]));event_dt=ticks/freq;delta=event_dt-dt
        flagged=exact&(abs(delta)>.1)
        checks.append(dict(run=run,component=c['component'],segments=c['segments'],selected_intervals=len(rows),exact_boundary_intervals=int(exact.sum()),relative_ticks_per_scaler_second=freq,event_minus_scaler_seconds_quantiles=np.quantile(delta[exact],[0,.01,.5,.99,1]).tolist(),disagreement_gt_0p1s=int(flagged.sum())))
        for i in np.flatnonzero(flagged):
            x=rows[i]
            bad.append(dict(run=run,run_type=r['run_type'],component=c['component'],segment_start=int(x['segment_start']),segment_end=int(x['segment_end']),ev_start=int(b[i,0]),ev_end=int(b[i,1]),scaler_seconds=float(dt[i]),event_seconds=float(event_dt[i]),difference_seconds=float(delta[i]),N=int(x['N']),A=float(x['A']),S=float(x['S']),D=float(x['D']),E_raw=int(x['E_raw']),charge_uC=float(x['charge_uC'])))
(P/'event_clock_checks.json').write_text(json.dumps(dict(components=checks,disagreements=bad),indent=2)+'\n')
print('Components',len(checks),'exact boundaries',sum(x['exact_boundary_intervals'] for x in checks),'disagreements',len(bad))
print(json.dumps(bad,indent=2))
