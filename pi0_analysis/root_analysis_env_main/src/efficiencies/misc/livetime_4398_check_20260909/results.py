"""Segment-0 measurements and figures; no production correction is applied.

Prerequisite: run inventory.py to read the explicitly named ROOT file.
All paths written below are relative to this script's directory.
"""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

P=Path(__file__).resolve().parent
s=np.load(P/'scalers.npz'); e=np.load(P/'events.npz')
ev=e['g.evnum']; raw=e['T.hms.hEDTM_tdcTimeRaw']; corrected=e['T.hms.hEDTM_tdcTime']
edges=s['evNumber']; time=s['H.1MHz.scalerTime']; dt=np.diff(time)
cur=s['H.BCM4A.scalerCurrent'][1:]; current_e=e['H.BCM4A.scalerCurrent']
lo,hi=33.3625,45.1375
good=(cur>=lo)&(cur<=hi)&(dt>0)
raw_ed=(raw>1)&(abs(raw-1765)<=500)
tight_ed=(raw>1)&(abs(corrected-245.11508)<=2)
counts=np.diff(np.searchsorted(ev,edges,side='left'))
edcounts=np.diff(np.searchsorted(ev[raw_ed],edges,side='left'))
tightcounts=np.diff(np.searchsorted(ev[tight_ed],edges,side='left'))
de=np.diff(s['H.EDTM.scaler']); ds=np.diff(s['H.hTRIG4.scaler']); da=np.diff(s['H.hL1ACCP.scaler'])
index=np.searchsorted(edges,ev,side='right')-1
valid=(index>=0)&(index<len(dt))
event_interval_good=valid&good[np.clip(index,0,len(dt)-1)]
rows=[]
for label,mask in [('all covered',dt>0),('I > 2 uA',(dt>0)&(cur>2)),
                   ('production current',good),
                   ('stable production current',good&(s['H.BCM4A.scalerCurrent'][:-1]>=lo)&(s['H.BCM4A.scalerCurrent'][:-1]<=hi))]:
    N=int(counts[mask].sum()); E=int(edcounts[mask].sum()); Et=int(tightcounts[mask].sum())
    S=int(ds[mask].sum()); D=int(de[mask].sum()); A=int(da[mask].sum())
    lt=E/D; p=N/S
    rows.append(dict(selection=label,time_s=float(dt[mask].sum()),N=N,E_raw=E,E_tight=Et,S=S,D=D,A=A,
                     EDTM_raw=lt,EDTM_tight=Et/D,CLTA=p,CLTP_raw=(N-E)/(S-D),CLTP_tight=(N-Et)/(S-D),
                     hardware_subtract_sent=(A-D)/(S-D),
                     conditional_binomial_sigma_EDTM=float(np.sqrt(lt*(1-lt)/D)),
                     conditional_binomial_sigma_CLTA=float(np.sqrt(p*(1-p)/S))))
with (P/'matched_estimators.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)

sens=[]
for halfwidth in [0.5,1,2,3,5,10,20,40,50]:
    n=int((event_interval_good&(raw>1)&(abs(corrected-245.11508)<=halfwidth)).sum())
    sens.append(dict(halfwidth_ns=halfwidth,accepted=n,ratio=n/de[good].sum()))
with (P/'timing_sensitivity.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=sens[0]);w.writeheader();w.writerows(sens)

windows=[]
for low,high in [(2,100),(25,45.1375),(30,45.1375),(33.3625,45.1375),(35,45.1375),(37,41),(38,40)]:
    m=(dt>0)&(cur>=low)&(cur<=high)
    N=int(counts[m].sum());E=int(edcounts[m].sum());D=int(de[m].sum());S=int(ds[m].sum())
    windows.append(dict(current_min=low,current_max=high,N=N,E=E,D=D,S=S,EDTM=E/D,CLTA=N/S))
with (P/'current_sensitivity.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=windows[0]);w.writeheader();w.writerows(windows)

summary={
 'scope':'Only cached updated segment 0; no final run-level or total-livetime recommendation.',
 'input':str(P/'inventory.json'), 'window_uA':[lo,hi],
 'legacy_EDTM_numerator':int((raw_ed&(current_e>=lo)&(current_e<=hi)).sum()),
 'legacy_EDTM_denominator':int(de[(cur>=lo)&(cur<=hi)].sum()),
 'uncovered_tail':dict(start_event=int(edges[-2]),end_event_exclusive=int(edges[-1]),
                       events=int(counts[-1]),edtm_raw=int(edcounts[-1]),dt=float(dt[-1]),
                       scaler_EDTM=int(de[-1]),scaler_TRIG4=int(ds[-1])),
 'measured_EDTM_rate_Hz':float(de[good].sum()/dt[good].sum()),
 'event_clock_matches_previous_scaler':bool(np.all(e['H.1MHz.scalerTime']==time[np.clip(index,0,len(time)-1)])),
 'matched_estimators':rows,
 'uncertainty_limit':'Binomial values are conditional scales only; periodic-pulser sampling, covariance and hardware coverage are not certified.',
}
(P/'results.json').write_text(json.dumps(summary,indent=2))

# Independent reader verification against ROOT's binary export, when available.
verification={}
for name,data in [('T',e),('TSH',s)]:
    fn=P/(name+'.bin'); cn=P/(name+'_columns.txt')
    if fn.exists() and cn.exists():
        columns=cn.read_text().splitlines(); a=np.fromfile(fn,dtype=np.float64).reshape(-1,len(columns))
        shared=[c for c in columns if c in data.files]
        bad=[c for c in shared if not np.array_equal(a[:,columns.index(c)],data[c],equal_nan=True)]
        verification[name]=dict(rows=len(a),columns_checked=len(shared),mismatches=bad)
        if bad: raise RuntimeError('ROOT/uproot mismatch: '+str(bad))
(P/'reader_verification.json').write_text(json.dumps(verification,indent=2))

plt.rcParams.update({'font.size':11,'axes.grid':True,'grid.alpha':0.25,'figure.dpi':150})
with PdfPages(P/'segment0_diagnostic_figures.pdf') as pdf:
    fig,axs=plt.subplots(3,1,figsize=(10,8),sharex=True,layout='constrained')
    axs[0].plot(time[1:],cur,lw=1);axs[0].axhspan(lo,hi,color='green',alpha=.1);axs[0].set_ylabel('BCM4A current (uA)')
    axs[1].plot(time[1:],np.cumsum(edcounts-de),lw=1.3);axs[1].scatter(time[-1],np.cumsum(edcounts-de)[-1],c='red',zorder=3)
    axs[1].set_ylabel('Cumulative EDTM\naccepted - scaler')
    axs[2].plot(time[1:],np.cumsum(counts-da),lw=1.3);axs[2].set_ylabel('Cumulative events\nrecorded - L1 scaler');axs[2].set_xlabel('Scaler clock (s)')
    fig.suptitle('Run 4398, updated segment 0: endpoint mismatch')
    pdf.savefig(fig);fig.savefig(P/'endpoint_mismatch.png');plt.close(fig)

    fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    q=event_interval_good&(raw>1)
    axs[0].hist(raw[q],bins=np.arange(1200,2301,10),histtype='step');axs[0].set_xlabel('EDTM raw TDC channels');axs[0].set_ylabel('Events / 10 channels')
    axs[0].axvspan(1265,2265,color='green',alpha=.1)
    axs[1].hist(corrected[q],bins=np.arange(240,290,.1),histtype='step');axs[1].set_yscale('log');axs[1].set_xlabel('Reference-subtracted EDTM time (ns)');axs[1].set_ylabel('Events / 0.1 ns')
    axs[1].axvspan(243.11508,247.11508,color='green',alpha=.1)
    fig.suptitle('Covered production-current intervals: raw and corrected timing')
    pdf.savefig(fig);fig.savefig(P/'edtm_timing.png');plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,4.5),layout='constrained')
    labels=['Existing event-current\nEDTM result','Matched raw-window\nEDTM','Matched +/-2 ns\nEDTM','Matched physics\ncomputer livetime']
    prod=rows[2];v=[summary['legacy_EDTM_numerator']/summary['legacy_EDTM_denominator'],prod['EDTM_raw'],prod['EDTM_tight'],prod['CLTP_tight']]
    ax.scatter(np.arange(4),v,s=65);ax.axhline(1,color='black',ls='--');ax.set_xticks(np.arange(4),labels);ax.set_ylabel('Measured count ratio')
    for i,y in enumerate(v):ax.annotate(f'{y:.8f}',(i,y),xytext=(0,10),textcoords='offset points',ha='center')
    ax.set_xlim(-.5,3.5);ax.set_ylim(min(v)-.0005,max(v)+.0005);fig.suptitle('Segment-0 ratios; no total-livetime prescription implied')
    pdf.savefig(fig);fig.savefig(P/'estimator_comparison.png');plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,4.5),layout='constrained')
    for field,label in [('EDTM','EDTM, raw window'),('CLTA','Computer livetime, all triggers')]:
        ax.plot(range(len(windows)),[r[field] for r in windows],'o-',label=label)
    ax.set_xticks(range(len(windows)),[f"{r['current_min']:g}--{r['current_max']:g}" for r in windows],rotation=25);ax.set_xlabel('Current interval (uA)');ax.set_ylabel('Matched count ratio');ax.legend();fig.suptitle('Current-cut sensitivity, segment 0')
    pdf.savefig(fig);plt.close(fig)
print(json.dumps(summary,indent=2));print('Reader verification',verification)
