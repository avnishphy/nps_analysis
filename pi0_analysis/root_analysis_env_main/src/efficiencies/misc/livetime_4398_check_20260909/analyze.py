"""Diagnostic counts, without choosing or applying a production correction."""
from pathlib import Path
import json,csv
import numpy as np
P=Path(__file__).resolve().parent
s=np.load(P/'scalers.npz'); e=np.load(P/'events.npz'); h=np.load(P/'helicity_scalers.npz')
lo,hi=33.3625,45.1375
en=e['g.evnum']; n=len(en)
ed=e['T.hms.hEDTM_tdcTimeRaw']; ec=e['H.BCM4A.scalerCurrent']
st=s['H.1MHz.scalerTime']; se=s['evNumber']; sc=s['H.BCM4A.scalerCurrent']
d={k:np.diff(s[k]) for k in s.files}
is_ed=(ed>1)&(np.abs(ed-1765)<=500)
peakbin=int(np.bincount((ed[(ed>1)&(ec>=lo)&(ec<=hi)]/10).astype(int)).argmax())
print('peakbin center',10*peakbin+5,'EDTM width is in raw channels')
print('event current distribution',np.quantile(ec,[0,.1,.5,.9,1]))
print('edtm raw >1',int((ed>1).sum()),'peak',int(is_ed.sum()),'out_of_window',int(((ed>1)&~is_ed).sum()))
print('all event types',np.unique(e['g.evtyp'],return_counts=True))
print('g.evnum steps',np.unique(np.diff(en),return_counts=True))
print('TSH evNumber',se[:6],se[-6:],'time',st[:6],st[-6:],'evcount',s['evcount'][-6:])
print('TSHelH evNumber',h['evNumber'][:6],h['evNumber'][-6:],'evcount',h['evcount'][-6:])
print('FIRST scalers',{k:float(s[k][0]) for k in ['H.1MHz.scalerTime','H.BCM4A.scalerCurrent','H.hTRIG4.scaler','H.hL1ACCP.scaler','H.EDTM.scaler']})
# A record at evNumber b is decoded before the physics event b in current hcana.
# Evaluate both half-open conventions to expose one-event boundary sensitivity.
results=[]
for bound in ['left','right']:
    counts=np.diff(np.searchsorted(en,se,side=bound))
    edcounts=np.diff(np.searchsorted(en[is_ed],se,side=bound))
    for label,sm,em in [
        ('all',np.ones(len(st)-1,dtype=bool),np.ones(n,dtype=bool)),
        ('current_gt2',sc[1:]>2,ec>2),
        ('production_current',(sc[1:]>=lo)&(sc[1:]<=hi),(ec>=lo)&(ec<=hi)),
        ('stable_current', (sc[1:]>=lo)&(sc[1:]<=hi)&(sc[:-1]>=lo)&(sc[:-1]<=hi), (ec>=lo)&(ec<=hi)),
    ]:
        a=float(d['H.hL1ACCP.scaler'][sm].sum()); tr=float(d['H.hTRIG4.scaler'][sm].sum()); de=float(d['H.EDTM.scaler'][sm].sum())
        row=dict(selection=label,boundary=bound,intervals=int(sm.sum()),time=float(d['H.1MHz.scalerTime'][sm].sum()),
                 charge=float(d['H.BCM4A.scalerCharge'][sm].sum()),scaler_edtm=de,scaler_trig4=tr,scaler_l1=a,
                 event_current_N=int(em.sum()),event_current_EDTM=int((em&is_ed).sum()),
                 aligned_N=int(counts[sm].sum()),aligned_EDTM=int(edcounts[sm].sum()),
                 edtm_event_current=float((em&is_ed).sum()/de),edtm_aligned=float(edcounts[sm].sum()/de),
                 hw_all=a/tr,hw_subtract_sent=(a-de)/(tr-de),
                 event_clt_all=float(counts[sm].sum()/tr),
                 event_clt_physics=float((counts[sm].sum()-edcounts[sm].sum())/(tr-de)))
        results.append(row);print(row)
    if bound=='left':
        rows=[]
        for i in range(len(st)-1):
            rows.append(dict(index=i+1,ev_start=se[i],ev_end=se[i+1],time_end=st[i+1],dt=st[i+1]-st[i],
                             current=sc[i+1],event_N=int(counts[i]),event_EDTM=int(edcounts[i]),
                             scaler_EDTM=d['H.EDTM.scaler'][i],scaler_TRIG4=d['H.hTRIG4.scaler'][i],
                             scaler_L1=d['H.hL1ACCP.scaler'][i],
                             production_current=bool(lo<=sc[i+1]<=hi)))
        with (P/'intervals.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
(P/'estimators.json').write_text(json.dumps(results,indent=2))
print('event scaler clock match fraction',np.mean(np.isin(e['H.1MHz.scalerTime'],st)))
idx=np.searchsorted(se,en,side='right')-1
idx=np.clip(idx,0,len(se)-1)
print('event clock matches preceding TSH by evNumber',np.mean(e['H.1MHz.scalerTime']==st[idx]))
for tag,mask in [('all',np.ones(n,bool)),('EDTM',is_ed),('non_EDTM',~is_ed)]:
    u,c=np.unique(e['g.trigbits'][mask],return_counts=True)
    print('trigger masks',tag,dict(zip(u.astype(int).tolist(),c.tolist())))
print('EDTM raw-window sensitivity production event-current')
for width in [50,100,200,300,400,500,750,1000]:
    print(width,int(((ed>1)&(np.abs(ed-1765)<=width)&(ec>=lo)&(ec<=hi)).sum()))
print('PRE width monitor count deltas')
for k in ['H.hPRE40.scaler','H.hPRE100.scaler','H.hPRE150.scaler','H.hPRE200.scaler']:
    print(k,float(d[k].sum()))
print('event pretrigger nonzero')
for k in e.files:
    if 'PRE' in k and k.endswith('tdcTimeRaw'): print(k,int(np.count_nonzero(e[k])))
