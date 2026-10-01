"""Diagnostic ratios on real, contiguous scaler coverage; no production edits.

Run: python3 analyze.py [run ...]. Reads export_cached.C output. A missing
export is pending, never silently accepted as complete. See REPORT.md for
conditional trigger/prescale and physical-latch assumptions.
"""
from pathlib import Path
import csv,json,sys,re,math
import numpy as np
P=Path(__file__).resolve().parent;B=P/'columns';O=P/'results';O.mkdir(exist_ok=True)
inventory=json.loads((P/'inventory.json').read_text())['runs']
def read(run,seg,name):
    stem=B/f'run{run}_seg{seg}_{name}'
    cols=Path(str(stem)+'_columns.txt').read_text().splitlines()
    a=np.fromfile(str(stem)+'.bin',dtype=np.float64).reshape(-1,len(cols))
    return {k:a[:,j] for j,k in enumerate(cols)}
def writecsv(p,rows):
    if not rows:return
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def div(a,b):return float(a/b) if b>0 else float('nan')
def peak(a,width):
    bins,counts=np.unique(np.floor(a[np.isfinite(a)]/width).astype(np.int64),return_counts=True)
    return float((bins[np.argmax(counts)]+.5)*width) if len(bins) else float('nan')
def joined(parts,key):return np.concatenate([q[key] for q in parts])
def run_one(q):
    run=q['run'];files=q['files'];flags=[]
    if not files:return dict(run=run,status='NO_CACHED_FILES',run_type=q['run_type'])
    missing=[f['segment'] for f in files if not (B/f"run{run}_seg{f['segment']}.done").exists()]
    if missing:return dict(run=run,status='PENDING_EXPORT',missing=missing)
    trig,setting=map(int,re.fullmatch(r'ps([1-6])=(-?\d+)',q['prescale_token']).groups());p=1 if setting<=0 else 2**(setting-1)+1
    sk=f'H.hTRIG{trig}.scaler';dk='H.EDTM.scaler';ak='H.hL1ACCP.scaler';clk='H.1MHz.scalerTime';curk='H.BCM4A.scalerCurrent';qk='H.BCM4A.scalerCharge'
    rawk='T.hms.hEDTM_tdcTimeRaw';timek='T.hms.hEDTM_tdcTime';multk='T.hms.hEDTM_tdcMultiplicity'
    ee=[];ss=[];segments=[];excluded=[];scaler_anomalies=[]
    for f in files:
        seg=f['segment'];sel=dict(line.split('\t',1) for line in (B/f'run{run}_seg{seg}_selection.tsv').read_text().splitlines())
        if sel['ok']!='1':excluded.append(dict(segment=seg,reason=sel['message']));continue
        e=read(run,seg,'T');s=read(run,seg,'TSH');lo,hi=float(sel['low']),float(sel['high'])
        for k in [sk,dk,ak,clk,curk,qk,'evNumber']:
            if k not in s:raise ValueError('Missing scaler '+k)
        for k in ['g.evnum','g.trigbits','g.evtyp',rawk,timek,multk,curk,clk]:
            if k not in e:raise ValueError('Missing event '+k)
        if len(e['g.evnum'])==0:excluded.append(dict(segment=seg,reason='empty T'));continue
        if np.any(np.diff(e['g.evnum'])!=1):raise ValueError(f'Internal event discontinuity {seg}')
        if np.any(np.diff(s[clk])<0):raise ValueError(f'Clock reset {seg}')
        zero=np.diff(s[clk])==0
        unchanged=np.ones(len(zero),bool)
        for k in [sk,dk,ak,qk]:
            unchanged &= np.diff(s[k])==0
            if np.any(np.diff(s[k])<0):raise ValueError('Counter reset '+k)
        duplicate=np.r_[False,zero&unchanged]
        anomalous=(zero&~unchanged)|(np.diff(s[dk])>1000+100*np.diff(s[clk]))
        for j in np.flatnonzero(anomalous):
            scaler_anomalies.append(dict(segment=seg,row_end=int(j+1),ev_start=float(s['evNumber'][j]),ev_end=float(s['evNumber'][j+1]),dt=float(s[clk][j+1]-s[clk][j]),D=float(s[dk][j+1]-s[dk][j]),S=float(s[sk][j+1]-s[sk][j]),A=float(s[ak][j+1]-s[ak][j]),current=float(s[curk][j+1]),in_current_window=bool(lo<=s[curk][j+1]<=hi)))
        oldgood=(s[curk][1:]>=lo)&(s[curk][1:]<=hi)
        oldD=float(np.diff(s[dk])[oldgood].sum())
        e['oldcurrent']=(e[curk]>=lo)&(e[curk]<=hi)
        e['segment']=np.full(len(e['g.evnum']),seg)
        s={k:v[~duplicate] for k,v in s.items()};s['segment']=np.full(len(s[clk]),seg);s['low']=np.full(len(s[clk]),lo);s['high']=np.full(len(s[clk]),hi)
        if len(s[clk])<2:raise ValueError(f'Insufficient real snapshots {seg}')
        idx=np.searchsorted(s['evNumber'],e['g.evnum'],side='right')-1
        inside=(idx>=0)&(idx<len(s[clk])-1)
        newcur=(s[curk][1:]>=lo)&(s[curk][1:]<=hi)
        e['perfilecovered']=inside;e['perfilealigned']=inside&newcur[np.clip(idx,0,len(newcur)-1)]
        staleidx=np.clip(idx,0,len(s[clk])-1)
        segments.append(dict(run=run,segment=seg,source=q['source'],events=len(e['g.evnum']),first_event=int(e['g.evnum'][0]),last_event=int(e['g.evnum'][-1]),real_snapshots=len(s[clk]),duplicate_rows=int(duplicate.sum()),first_boundary=int(s['evNumber'][0]),last_boundary=int(s['evNumber'][-1]),low_uA=lo,high_uA=hi,old_D=oldD,hel_charge_before=float(sel['hel_charge_before']),hel_charge_after=float(sel['hel_charge_after']),event_clock_matches_previous=int(np.sum(e[clk]==s[clk][staleidx])),event_current_matches_previous=int(np.sum(e[curk]==s[curk][staleidx]))))
        ee.append(e);ss.append(s)
    if not ee:return dict(run=run,status='NO_VALID_SELECTION',excluded=excluded,run_type=q['run_type'])
    e={k:joined(ee,k) for k in ee[0]};raw=e[rawk];tc=e[timek];en=e['g.evnum'];oldmask=e['oldcurrent']&(raw>1)&np.isfinite(raw)
    rawpeak=peak(raw[oldmask],10);rawed=(raw>1)&np.isfinite(raw)&(abs(raw-rawpeak)<=500)
    # Run-specific corrected-time peak: dense 0.1-ns bin, then median within 1 ns.
    tcrough=peak(tc[oldmask&rawed],.1)
    tcnear=tc[oldmask&rawed&(abs(tc-tcrough)<=1)]
    tcpeak=float(np.median(tcnear)) if len(tcnear) else float('nan')
    tight=(raw>1)&np.isfinite(tc)&(abs(tc-tcpeak)<=2)
    oldE=int((rawed&e['oldcurrent']).sum());oldD=sum(x['old_D'] for x in segments)
    perfileE=int((rawed&e['oldcurrent']&e['perfilecovered']).sum());perfileAlignedE=int((rawed&e['perfilealigned']).sum())
    # Components must be neighboring file numbers AND exact event continuity.
    groups=[];component_breaks=[]
    for j in range(len(ss)):
        reasons=[]
        if j:
            if segments[j]['segment']!=segments[j-1]['segment']+1:reasons.append('missing_segment')
            if segments[j]['first_event']!=segments[j-1]['last_event']+1:reasons.append('event_gap')
            for k in [clk,sk,dk,ak,qk]:
                if ss[j][k][0]<ss[j-1][k][-1]:reasons.append('decreasing_'+k)
        if j==0 or reasons:
            groups.append([])
            if reasons:component_breaks.append(dict(before_segment=segments[j]['segment'],reasons=reasons))
        groups[-1].append(j)
    intervals=[];timing=[];components=[];hist=[];allmasks=[]
    for component,g in enumerate(groups):
        s={k:joined([ss[j] for j in g],k) for k in ss[0]}
        eg={k:joined([ee[j] for j in g],k) for k in ee[0]};ev=eg['g.evnum'];rt=eg[rawk];tt=eg[timek]
        unchanged=np.ones(len(s[clk])-1,bool)
        for k in [sk,dk,ak,qk]:
            unchanged &= np.diff(s[k])==0
        dup=np.r_[False,(np.diff(s[clk])==0)&unchanged]
        s={k:v[~dup] for k,v in s.items()}
        b=s['evNumber'];dt=np.diff(s[clk])
        if np.any(dt<0) or np.any(np.diff(b)<0):raise ValueError('Decreasing component snapshots')
        for k in [sk,dk,ak,qk]:
            if np.any(np.diff(s[k])<0):raise ValueError('Cross-file counter reset '+k)
        # Require the entire snapshot interval to be covered by recorded events.
        covered=(b[:-1]>=ev[0])&(b[1:]<=ev[-1]+1)
        current=(s[curk][1:]>=s['low'][1:])&(s[curk][1:]<=s['high'][1:])
        previouscurrent=(s[curk][:-1]>=s['low'][:-1])&(s[curk][:-1]<=s['high'][:-1])
        rb=(rt>1)&(abs(rt-rawpeak)<=500);tb=(rt>1)&(abs(tt-tcpeak)<=2)
        def counts(mask=None,side='left'):
            return np.diff(np.searchsorted(ev if mask is None else ev[mask],b,side=side))
        n=counts();er=counts(rb);et=counts(tb);na=counts(side='right');ea=counts(tb,side='right')
        d=np.diff(s[dk]);tr=np.diff(s[sk]);a=np.diff(s[ak]);charge=np.diff(s[qk])
        where=np.searchsorted(b,ev,side='right')-1;inside=(where>=0)&(where<len(current))
        nominal=inside&(current&covered)[np.clip(where,0,len(current)-1)]
        for width in [.5,1,2,3,5,10,20,50]:
            timing.append(dict(component=component,width_ns=width,E=int((nominal&(rt>1)&(abs(tt-tcpeak)<=width)).sum()),D=float(d[current&covered].sum())))
        bins=np.arange(-60,60.00001,.1);h,_=np.histogram((tt-tcpeak)[nominal&(rt>1)],bins=bins);hist.append(h)
        for i in range(len(dt)):
            intervals.append(dict(run=run,component=component,index=i,segment_start=int(s['segment'][i]),segment_end=int(s['segment'][i+1]),ev_start=int(b[i]),ev_end=int(b[i+1]),t_start=float(s[clk][i]),t_end=float(s[clk][i+1]),dt=float(dt[i]),current=float(s[curk][i+1]),low=float(s['low'][i+1]),high=float(s['high'][i+1]),covered=int(covered[i]),nominal=int(current[i]&covered[i]),stable=int(current[i]&previouscurrent[i]&covered[i]),beam2=int(s[curk][i+1]>2 and covered[i]),N=int(n[i]),E_raw=int(er[i]),E_tight=int(et[i]),S=float(tr[i]),D=float(d[i]),A=float(a[i]),N_alt=int(na[i]),E_alt=int(ea[i]),charge_uC=float(charge[i]),s1x_rate=float(s['H.S1X.scalerRate'][i+1])))
        components.append(dict(component=component,segments=[segments[j]['segment'] for j in g],first_event=int(ev[0]),last_event=int(ev[-1]),first_boundary=int(b[0]),last_boundary=int(b[-1]),head_events=int((ev<b[0]).sum()),tail_events=int((ev>=b[-1]).sum()),head_raw_EDTM=int((rb&(ev<b[0])).sum()),tail_raw_EDTM=int((rb&(ev>=b[-1])).sum()),cross_duplicate_rows=int(dup.sum())))
    if len(groups)>1:flags.append('multiple_coverage_components')
    if any(any(x.startswith('decreasing_') for x in b['reasons']) for b in component_breaks):flags.append('cross_segment_counter_or_clock_restart')
    if files[0]['segment']!=0:flags.append('missing_start_segment')
    if excluded:flags.append('selection_exclusions')
    if scaler_anomalies:flags.append('scaler_jump_or_clock_anomaly')
    def sums(selection):
        ii=[x for x in intervals if x[selection]]
        z={k:sum(x[k] for x in ii) for k in ['N','E_raw','E_tight','S','D','A','N_alt','E_alt','dt','charge_uC']};z['intervals']=len(ii)
        N,E,S,D=z['N'],z['E_tight'],z['S'],z['D']
        z.update(EDTM_raw=div(p*z['E_raw'],D),EDTM_tight=div(p*E,D),CLT_all=div(p*N,S),CLT_physics=div(p*(N-E),S-D),L1_over_input=div(p*z['A'],S),EDTM_alternative=div(p*z['E_alt'],D),CLT_alternative=div(p*(z['N_alt']-z['E_alt']),S-D))
        qr=div(E,D);qc=div(N-E,S-D)
        z['EDTM_conditional_sigma']=p*math.sqrt(qr*(1-qr)/D) if D>0 and 0<=qr<=1 else float('nan')
        z['CLT_conditional_sigma']=p*math.sqrt(qc*(1-qc)/(S-D)) if S>D and 0<=qc<=1 else float('nan')
        return z
    nominal=sums('nominal');stable=sums('stable');beam2=sums('beam2')
    if nominal['EDTM_tight']>1:flags.append('EDTM_above_one')
    if nominal['CLT_physics']>1:flags.append('CLT_above_one')
    if p>1:flags.append('prescale_sampling_conditional')
    bits=e['g.trigbits'].astype(np.int64)
    triggerbit_fraction=float(np.mean((bits&(1<<(trig-1)))!=0))
    if triggerbit_fraction!=1:flags.append('configured_trigger_not_in_every_mask')
    if np.any(e['g.evtyp']!=1):flags.append('non_type1_events')
    if not q['same_saved_files']:flags.append('coverage_differs_from_saved')
    timingtot=[dict(width_ns=w,E=sum(x['E'] for x in timing if x['width_ns']==w),D=nominal['D'],ratio=div(p*sum(x['E'] for x in timing if x['width_ns']==w),nominal['D'])) for w in [.5,1,2,3,5,10,20,50]]
    z=dict(run=run,status='CALCULATED',run_type=q['run_type'],source=q['source'],prescale_token=q['prescale_token'],prescale_factor=p,trigger=trig,triggerbit_fraction=triggerbit_fraction,events=len(en),segments=len(segments),components=components,excluded=excluded,flags=flags,same_saved_files=q['same_saved_files'],saved=q['saved'],raw_peak_channels=rawpeak,corrected_peak_ns=tcpeak,corrected_peak_seed_fraction=div(len(tcnear),int(oldmask.sum())),old_E=oldE,old_D=oldD,old_ratio=div(p*oldE,oldD),perfile_covered_oldcurrent_E=perfileE,perfile_covered_oldcurrent_ratio=div(p*perfileE,oldD),perfile_covered_aligned_E=perfileAlignedE,perfile_covered_aligned_ratio=div(p*perfileAlignedE,oldD),nominal=nominal,stable=stable,beam2=beam2,timing=timingtot,multiplicity_counts={str(int(k)):int(v) for k,v in zip(*np.unique(e[multk][oldmask&rawed],return_counts=True))},segment_details=segments)
    # Full interval records permit checking each count and locating real bursts.
    writecsv(O/f'run{run}_intervals.csv',intervals)
    np.savez_compressed(O/f'run{run}_timing_hist.npz',centers=(bins[1:]+bins[:-1])/2,counts=np.sum(hist,axis=0))
    # Conditional time-block resampling; not a final systematic uncertainty.
    rng=np.random.default_rng(run);boot=[]
    for width in [10,20,60]:
        blocks={}
        for x in intervals:
            if x['nominal']:
                key=(x['component'],int(x['t_end']//width));blocks.setdefault(key,np.zeros(4));blocks[key]+=np.array([x['N'],x['E_tight'],x['S'],x['D']])
        values=np.array(list(blocks.values()))
        if len(values)>1:
            sample=values[rng.integers(0,len(values),size=(1000,len(values)))].sum(axis=1)
            ed=p*sample[:,1]/sample[:,3];cl=p*(sample[:,0]-sample[:,1])/(sample[:,2]-sample[:,3])
            boot.append(dict(block_seconds=width,blocks=len(values),replicates=1000,seed=run,EDTM_sigma=float(ed.std(ddof=1)),CLT_sigma=float(cl.std(ddof=1))))
    z['block_bootstrap']=boot
    z['scaler_anomalies']=scaler_anomalies
    z['component_breaks']=component_breaks
    return z
if __name__=='__main__':
    chosen=set(map(int,sys.argv[1:]))
    for q in inventory:
        if chosen and q['run'] not in chosen:continue
        try:z=run_one(q)
        except Exception as e:z=dict(run=q['run'],status='ERROR',error=str(e))
        (O/f"run{q['run']}.json").write_text(json.dumps(z,indent=2))
        if z['status']=='CALCULATED':print(z['run'],z['events'],'old',z['old_ratio'],'matched',z['nominal']['EDTM_tight'],'CLT',z['nominal']['CLT_physics'],'flags',','.join(z['flags']),flush=True)
        else:print(z,flush=True)
