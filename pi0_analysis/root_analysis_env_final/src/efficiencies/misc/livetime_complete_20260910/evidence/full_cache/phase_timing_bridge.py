"""Link pulser phase to EDTM/trigger timing without assigning causal triggers.

Near-phase positive hits can be a pulser recorded in an earlier physics event.
Pulse association alone does not establish that the pulser caused the event.
Only production runs; identical nominal scaler intervals. Core plot uses held-
out events only and deterministic thinning; all noncore near-phase rows saved.
"""
from pathlib import Path
import json,csv,numpy as np
from timestamp_probe import read
P=Path(__file__).resolve().parent;rows=[];summary=[]
for f in sorted((P/'results').glob('run*.json')):
    r=json.loads(f.read_text())
    if r['status']!='CALCULATED' or r['run_type']!='production':continue
    run=r['run'];ii=[x for x in csv.DictReader((P/f'results/run{run}_intervals.csv').open()) if x['nominal']=='1']
    for sd in r['segment_details']:
        seg=sd['segment'];e=read(run,seg,'T');t=e['g.evtime'];ev=e['g.evnum'];raw=e['T.hms.hEDTM_tdcTimeRaw'];tc=e['T.hms.hEDTM_tdcTime'];tr=e[f'T.hms.hTRIG{r["trigger"]}_tdcTime']
        core=(raw>1)&(abs(tc-r['corrected_peak_ns'])<=2);rawmask=(raw>1)&(abs(raw-r['raw_peak_channels'])<=500)
        ct=t[core];dt=np.diff(ct);dt=dt[dt>0];u,n=np.unique(np.rint(dt/1000).astype(np.int64),return_counts=True);p0=float(u[n.argmax()]*1000)
        k=np.rint(dt/p0);near=(k>=1)&(abs(dt-k*p0)<.02*p0);period=float(np.median(dt[near]/k[near]));coef=np.polyfit(np.rint((ct-ct[0])/period),ct-ct[0],1);period=float(coef[0]);origin=float(ct[0]+coef[1])
        ci=np.flatnonzero(core);train=ci[::2];test=ci[1::2];nn=np.rint((t[train]-origin)/period);tt=t[train]-origin;good=np.r_[True,np.diff(nn)>0];nn=nn[good];tt=tt[good]
        residual=t-origin-np.interp(np.rint((t-origin)/period),nn,tt,left=np.nan,right=np.nan)
        selected=np.zeros(len(t),bool)
        for x in ii:
            if int(x['segment_start'])<=seg<=int(x['segment_end']):
                lo,hi=np.searchsorted(ev,[int(x['ev_start']),int(x['ev_end'])]);selected[lo:hi]=True
        relcenter=float(np.median((tc-tr)[core]));held=np.zeros(len(t),bool);held[test]=True
        for name,mask in [('core_heldout',held),('raw_tail',rawmask&~core),('positive_outside_raw',(raw>1)&~rawmask)]:
            m=selected&mask;finite=m&np.isfinite(residual)
            summary.append(dict(run=run,segment=seg,category=name,total=int(m.sum()),model_missing=int((m&~np.isfinite(residual)).sum()),phase10=int((finite&(abs(residual)<=10)).sum()),phase25=int((finite&(abs(residual)<=25)).sum()),phase125=int((finite&(abs(residual)<=125)).sum()),phase500=int((finite&(abs(residual)<=500)).sum())))
            ix=np.flatnonzero(finite&(abs(residual)<=500))
            if name=='core_heldout':ix=ix[::max(1,len(ix)//50)]
            for i in ix:rows.append(dict(run=run,segment=seg,event=int(ev[i]),category=name,phase_ticks=float(residual[i]),raw_delta_channels=float(raw[i]-r['raw_peak_channels']),corrected_delta_ns=float(tc[i]-r['corrected_peak_ns']),relative_delta_ns=float(tc[i]-tr[i]-relcenter),multiplicity=int(e['T.hms.hEDTM_tdcMultiplicity'][i])))
(P/'phase_timing_bridge.json').write_text(json.dumps(dict(summary=summary,points=rows),indent=2)+'\n')
for name in ['core_heldout','raw_tail','positive_outside_raw']:
    a=[x for x in rows if x['category']==name];s=[x for x in summary if x['category']==name]
    print(name,{k:sum(x[k] for x in s) for k in ['total','model_missing','phase10','phase25','phase125','phase500']})
    if name!='core_heldout' and len(a)>2:
        x=np.array([q['phase_ticks'] for q in a]);y=np.array([q['relative_delta_ns'] for q in a]);fit=np.polyfit(x,y,1)
        print('relative timing vs phase: slope ns/tick, intercept ns',fit.tolist(),'correlation',float(np.corrcoef(x,y)[0,1]))
