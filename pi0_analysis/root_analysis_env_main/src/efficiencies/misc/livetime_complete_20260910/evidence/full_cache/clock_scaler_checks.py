"""Compare periodic-pulser counts, event clock and scaler clock without repairing counters."""
from pathlib import Path
import json,csv,numpy as np
from timestamp_probe import read
P=Path(__file__).resolve().parent;rows=[];anomalies=[]
for f in sorted((P/'results').glob('run*.json')):
    r=json.loads(f.read_text())
    if r['status']!='CALCULATED':continue
    run=r['run'];phase=P/f'timestamp/run{run}.json'
    if not phase.exists():continue
    pp={s['segment']:s for s in json.loads(phase.read_text())['segments']}
    for d in r['segment_details']:
        seg=d['segment'];e=read(run,seg,'T');s=read(run,seg,'TSH');b=s['evNumber'];sc=s['H.1MHz.scalerTime']
        j=np.searchsorted(e['g.evnum'],b);m=(j>0)&(j<len(e['g.evnum']))
        x=sc[m];y=e['g.evtime'][j[m]];co=np.polyfit(x-x[0],y-y[0],1)
        err=(y-y[0])-np.polyval(co,x-x[0]);good=abs(err-np.median(err))<max(1e5,10*np.median(abs(err-np.median(err))))
        co=np.polyfit((x-x[0])[good],(y-y[0])[good],1)
        freq=float(co[0]);period=pp[seg]['period_ticks']/freq;dt=np.diff(sc);D=np.diff(s['H.EDTM.scaler']);delta=D-dt/period
        current=s['H.BCM4A.scalerCurrent'][1:];beam=(current>=d['low_uA'])&(current<=d['high_uA']);valid=beam&(dt>0)
        row=dict(run=run,segment=seg,event_ticks_per_scaler_second=freq,pulser_period_seconds=period,pulser_rate_hz=1/period,fit_points=int(good.sum()),beam_intervals=int(valid.sum()),D_minus_expected_quantiles=np.quantile(delta[valid],[0,.01,.5,.99,1]).tolist(),selected_abs_excess_gt2=int((valid&(abs(delta)>2)).sum()))
        rows.append(row)
        for i in np.flatnonzero(beam&(abs(delta)>2)):
            anomalies.append(dict(run=run,segment=seg,ev_start=int(b[i]),ev_end=int(b[i+1]),dt=float(dt[i]),D=float(D[i]),expected=float(dt[i]/period),excess=float(delta[i]),current=float(current[i])))
(P/'clock_scaler_checks.json').write_text(json.dumps(dict(segments=rows,selected_anomalies=anomalies),indent=2)+'\n')
print('Clock-scalers checked',len(rows),'selected pulse-count inconsistencies',len(anomalies))
print(json.dumps(anomalies[:10],indent=2))
