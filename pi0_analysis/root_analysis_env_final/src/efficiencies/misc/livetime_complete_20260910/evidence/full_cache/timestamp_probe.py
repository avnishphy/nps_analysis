"""Independent event-clock pulse-phase diagnostic; no correction adopted.

Fit pulse period from tight-tagged events. Cross-check phase with held-out
core events and classify raw/tight timing tails on identical selected intervals.
All phase units below are g.evtime ticks; scale to ns only after clock validation.
"""
from pathlib import Path
import json,sys,numpy as np
P=Path(__file__).resolve().parent;O=P/'timestamp';O.mkdir(exist_ok=True)
def read(run,seg,tree):
    stem=P/f'columns/run{run}_seg{seg}_{tree}'
    names=Path(str(stem)+'_columns.txt').read_text().splitlines()
    a=np.memmap(str(stem)+'.bin',mode='r',dtype='float64').reshape(-1,len(names))
    return {n:a[:,i] for i,n in enumerate(names)}
def probe(run):
    r=json.loads((P/f'results/run{run}.json').read_text());peak=r['corrected_peak_ns'];rawpeak=r['raw_peak_channels']
    import csv
    ii=[x for x in csv.DictReader((P/f'results/run{run}_intervals.csv').open()) if x['nominal']=='1']
    rows=[];examples=[];hist=np.zeros((4,400),dtype=np.int64)
    for sd in r['segment_details']:
        seg=sd['segment'];e=read(run,seg,'T');s=read(run,seg,'TSH')
        t=e['g.evtime'];ev=e['g.evnum'];raw=e['T.hms.hEDTM_tdcTimeRaw'];tc=e['T.hms.hEDTM_tdcTime']
        core=(raw>1)&(abs(tc-peak)<=2);rawmask=(raw>1)&(abs(raw-rawpeak)<=500)
        ct=t[core];dt=np.diff(ct);dt=dt[dt>0]
        # Most common coarse spacing, not median (prescaling skips pulses).
        u,n=np.unique(np.rint(dt/1000).astype(np.int64),return_counts=True);p0=float(u[n.argmax()]*1000)
        k=np.rint(dt/p0);near=(k>=1)&(abs(dt-k*p0)<.02*p0)
        period=float(np.median(dt[near]/k[near]))
        # One absolute pulse index; interpolate alternate core anchors to allow drift.
        ncore=np.rint((ct-ct[0])/period);coef=np.polyfit(ncore,ct-ct[0],1)
        period=float(coef[0]);origin=float(ct[0]+coef[1])
        selected=np.zeros(len(t),bool)
        for x in ii:
            if int(x['segment_start'])<=seg<=int(x['segment_end']):
                lo,hi=np.searchsorted(ev,[int(x['ev_start']),int(x['ev_end'])]);selected[lo:hi]=True
        ci=np.flatnonzero(core);train=ci[::2];test=ci[1::2]
        nn=np.rint((t[train]-origin)/period);tt=t[train]-origin
        good=np.r_[True,np.diff(nn)>0];nn=nn[good];tt=tt[good]
        ni=np.rint((t-origin)/period)
        expected=np.interp(ni,nn,tt,left=np.nan,right=np.nan)
        residual=t-origin-expected
        validations=residual[test];validations=validations[np.isfinite(validations)]
        fitted=np.diff(tt)/np.diff(nn)
        # A 125-tick window corresponds to 0.5 us ONLY if 250 MHz is verified.
        nearphase=abs(residual)<=125;side=(abs(residual)>=500)&(abs(residual)<1500)
        row=dict(run=run,segment=seg,period_ticks=period,period_block_median=float(np.median(fitted)),core_clock_residual_quantiles_ticks=np.quantile(np.abs(validations),[.5,.9,.99,.999]).tolist(),events=int(selected.sum()),phase_model_events=int((selected&np.isfinite(residual)).sum()),scaler_clock_span=float(s['H.1MHz.scalerTime'][-1]-s['H.1MHz.scalerTime'][0]),event_clock_span_ticks=float(t[-1]-t[0]))
        for name,mask in [('core',core),('raw',rawmask),('raw_not_core',rawmask&~core),('positive_not_core',(raw>1)&~core),('zero_raw',raw<=1),('all',np.ones(len(t),bool))]:
            m=selected&mask;row[name]=dict(total=int(m.sum()),phase125=int((m&nearphase).sum()),side500_1500=int((m&side).sum()),model_missing=int((m&~np.isfinite(residual)).sum()))
            row[name]['phase_windows']={str(w):int((m&(abs(residual)<=w)).sum()) for w in [10,25,50,125,250,500]}
        for j,m in enumerate([core,rawmask&~core,(raw>1)&~rawmask,raw<=1]):
            hist[j]+=np.histogram(residual[selected&m],bins=np.linspace(-2000,2000,401))[0]
        interesting=selected&((rawmask&~core)|((raw<=1)&nearphase)|(core&~nearphase))
        for i in np.flatnonzero(interesting)[:150]:
            examples.append(dict(run=run,segment=seg,event=int(ev[i]),ticks=float(t[i]),raw=float(raw[i]),corrected=float(tc[i]),multiplicity=float(e['T.hms.hEDTM_tdcMultiplicity'][i]),residual_ticks=float(residual[i]),raw_window=bool(rawmask[i]),core=bool(core[i])))
        rows.append(row)
    (O/f'run{run}.json').write_text(json.dumps(dict(run=run,segments=rows,examples=examples),indent=2)+'\n')
    np.savez_compressed(O/f'run{run}_phase_hist.npz',edges=np.linspace(-2000,2000,401),counts=hist)
    names=['core','raw_not_core','positive_not_core','zero_raw','all']
    summary={n:{k:sum(r[n][k] for r in rows) for k in ['total','phase125','side500_1500','model_missing']} for n in names}
    print(run,json.dumps(summary),flush=True)
for run in map(int,sys.argv[1:]):probe(run)
