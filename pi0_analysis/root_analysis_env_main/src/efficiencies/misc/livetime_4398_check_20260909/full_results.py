"""Join real cumulative snapshots across all six updated run-4398 segments.

Prerequisite: export_full_run.C. Output remains diagnostic: hardware attribution
and total-livetime coverage require independent run-period evidence.
"""
from pathlib import Path
import csv,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
P=Path(__file__).resolve().parent
B=P/'full_run'
def read(seg,name):
    cols=(B/f'seg{seg}_{name}_columns.txt').read_text().splitlines()
    data=np.fromfile(B/f'seg{seg}_{name}.bin',dtype=np.float64).reshape(-1,len(cols))
    return {k:data[:,j] for j,k in enumerate(cols)}
log=(P/'full_export.log').read_text()
if len(re.findall('EXPORTED segment=',log))!=18:
    raise RuntimeError('Need all 18 completed tree exports before calculating')
ss=[];ee=[];segrows=[];removed=[];naive=[]
lo,hi=33.3625,45.1375
for seg in range(6):
    s=read(seg,'TSH');e=read(seg,'T');ee.append(e)
    t=s['H.1MHz.scalerTime'];idx=np.r_[True,np.diff(t)>0]
    if np.any(np.diff(t)<0):raise RuntimeError('Clock reset inside segment')
    for row in np.flatnonzero(~idx):
        for k in ['H.EDTM.scaler','H.hTRIG4.scaler','H.hL1ACCP.scaler']:
            if s[k][row]!=s[k][row-1]:raise RuntimeError('Zero clock interval with changed '+k)
        removed.append(dict(segment=seg,row=int(row),event_boundary=float(s['evNumber'][row]),clock=float(t[row])))
    ss.append({k:v[idx] for k,v in s.items()})
    en=e['g.evnum']; raw=e['T.hms.hEDTM_tdcTimeRaw'];cur=e['H.BCM4A.scalerCurrent']
    current_ok=(cur>=lo)&(cur<=hi);raw_ed=(raw>1)&(abs(raw-1765)<=500)
    mask=(s['H.BCM4A.scalerCurrent'][1:]>=lo)&(s['H.BCM4A.scalerCurrent'][1:]<=hi)
    naive.append(dict(segment=seg,E=int((current_ok&raw_ed).sum()),D=float(np.diff(s['H.EDTM.scaler'])[mask].sum())))
    segrows.append(dict(segment=seg,events=len(en),first_event=int(en[0]),last_event=int(en[-1]),
                        first_real_scaler_event=int(s['evNumber'][idx][0]),last_real_scaler_event=int(s['evNumber'][idx][-1]),
                        first_clock=float(t[idx][0]),last_clock=float(t[idx][-1]),
                        scaler_rows=len(t),real_scaler_rows=int(idx.sum()),
                        event_gaps=int(np.sum(np.diff(en)!=1))))
s={k:np.concatenate([q[k] for q in ss]) for k in ss[0]}
e={k:np.concatenate([q[k] for q in ee]) for k in ee[0]}
en=e['g.evnum'];bound=s['evNumber'];time=s['H.1MHz.scalerTime'];dt=np.diff(time)
if np.any(np.diff(en)!=1):raise RuntimeError('Event continuity is not exact')
if np.any(np.diff(bound)<=0) or np.any(dt<=0):raise RuntimeError('Non-increasing snapshot boundaries')
reset_counts={k:int(np.sum(np.diff(s[k])<0)) for k in ['H.EDTM.scaler','H.hTRIG4.scaler','H.hL1ACCP.scaler','H.BCM4A.scalerCharge']}
if any(reset_counts.values()):raise RuntimeError('Cumulative reset needs explicit treatment: '+str(reset_counts))
cur=s['H.BCM4A.scalerCurrent'][1:];raw=e['T.hms.hEDTM_tdcTimeRaw'];tc=e['T.hms.hEDTM_tdcTime']
raw_ed=(raw>1)&(abs(raw-1765)<=500)
tight_ed=(raw>1)&(abs(tc-245.11508)<=2)
counts=np.diff(np.searchsorted(en,bound,side='left'))
edcounts=np.diff(np.searchsorted(en[raw_ed],bound,side='left'))
tightcounts=np.diff(np.searchsorted(en[tight_ed],bound,side='left'))
ds=np.diff(s['H.hTRIG4.scaler']);de=np.diff(s['H.EDTM.scaler']);da=np.diff(s['H.hL1ACCP.scaler'])
rows=[]
for label,m in [('all covered',np.ones(len(dt),bool)),('I > 2 uA',cur>2),
                ('fixed production current',(cur>=lo)&(cur<=hi)),
                ('stable production current',(cur>=lo)&(cur<=hi)&(s['H.BCM4A.scalerCurrent'][:-1]>=lo)&(s['H.BCM4A.scalerCurrent'][:-1]<=hi))]:
    N=int(counts[m].sum());Er=int(edcounts[m].sum());Et=int(tightcounts[m].sum());D=int(de[m].sum());S=int(ds[m].sum());A=int(da[m].sum())
    rows.append(dict(selection=label,intervals=int(m.sum()),time_s=float(dt[m].sum()),
                     charge_uC=float(np.diff(s['H.BCM4A.scalerCharge'])[m].sum()),
                     N=N,A=A,S=S,D=D,E_raw=Er,E_tight=Et,EDTM_raw=Er/D,EDTM_tight=Et/D,
                     CLTA=N/S,CLTP_raw=(N-Er)/(S-D),CLTP_tight=(N-Et)/(S-D),
                     HW_all=A/S,HW_subtract_sent=(A-D)/(S-D),
                     conditional_sigma_EDTM=float(np.sqrt((Et/D)*(1-Et/D)/D)),
                     conditional_sigma_CLTP=float(np.sqrt(((N-Et)/(S-D))*(1-(N-Et)/(S-D))/(S-D)))))
def writecsv(name,rows):
    with (P/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
writecsv('full_run_estimators.csv',rows);writecsv('full_run_segments.csv',segrows)
intervals=[]
for i in range(len(dt)):
    intervals.append(dict(index=i,ev_start=bound[i],ev_end=bound[i+1],time_start=time[i],time_end=time[i+1],
                         dt=dt[i],current=cur[i],N=int(counts[i]),A=int(da[i]),S=int(ds[i]),D=int(de[i]),
                         E_raw=int(edcounts[i]),E_tight=int(tightcounts[i])))
writecsv('full_run_intervals.csv',intervals)
tail=(en>=bound[-1]);head=(en<bound[0])
summary=dict(scope='All six updated segments; matched real-snapshot coverage, not a certified total-livetime correction.',
             total_events=len(en),first_event=int(en[0]),last_event=int(en[-1]),
             snapshots=len(bound),removed_duplicate_rows=removed,counter_resets=reset_counts,
             boundary_range=[float(bound[0]),float(bound[-1])],clock_range=[float(time[0]),float(time[-1])],
             uncovered_head_events=int(head.sum()),uncovered_tail_events=int(tail.sum()),
             uncovered_tail_EDTM_raw=int((tail&raw_ed).sum()),
             fixed_current_window_uA=[lo,hi],naive_per_segment=naive,
             naive_sum_ratio=sum(q['E'] for q in naive)/sum(q['D'] for q in naive),
             estimators=rows,segments=segrows)
(P/'full_run_results.json').write_text(json.dumps(summary,indent=2))

# Correlation suggests candidate routing; it does not establish trigger identity.
pre=[]
m=(cur>=lo)&(cur<=hi)
for arm in ['h','p']:
    base=np.diff(s[f'H.{arm}PRE100.scaler'])[m].sum()
    for width in [40,100,150,200]:
        v=np.diff(s[f'H.{arm}PRE{width}.scaler'])[m].sum()
        pre.append(dict(branch=f'H.{arm}PRE{width}.scaler',counts=float(v),ratio_to_PRE100=float(v/base)))
writecsv('full_run_PRE_diagnostics.csv',pre)

windows=[]
for low,high in [(2,100),(25,45.1375),(30,45.1375),(33.3625,45.1375),(35,45.1375),(37,41),(38,40)]:
    m=(cur>=low)&(cur<=high);N=counts[m].sum();E=tightcounts[m].sum();D=de[m].sum();S=ds[m].sum()
    windows.append(dict(low=low,high=high,N=int(N),E=int(E),D=int(D),S=int(S),EDTM=float(E/D),CLTP=float((N-E)/(S-D))))
writecsv('full_run_current_sensitivity.csv',windows)
where=np.searchsorted(bound,en,side='right')-1;inside=(where>=0)&(where<len(cur))
prod=inside&((cur>=lo)&(cur<=hi))[np.clip(where,0,len(cur)-1)]
timing=[]
for halfwidth in [.5,1,2,3,5,10,20,40,50]:
    n=int((prod&(raw>1)&(abs(tc-245.11508)<=halfwidth)).sum())
    timing.append(dict(halfwidth_ns=halfwidth,E=n,D=rows[2]['D'],EDTM=n/rows[2]['D']))
writecsv('full_run_timing_sensitivity.csv',timing)
prod_intervals=(cur>=lo)&(cur<=hi)
alt_N=int(np.diff(np.searchsorted(en,bound,side='right'))[prod_intervals].sum())
alt_E=int(np.diff(np.searchsorted(en[tight_ed],bound,side='right'))[prod_intervals].sum())
boundary_check=dict(primary='[b_i,b_(i+1))',alternative='(b_i,b_(i+1)]',
                    alternative_N=alt_N,alternative_E=alt_E,D=rows[2]['D'],S=rows[2]['S'],
                    alternative_EDTM=alt_E/rows[2]['D'],
                    alternative_CLTP=(alt_N-alt_E)/(rows[2]['S']-rows[2]['D']))
(P/'full_run_boundary_sensitivity.json').write_text(json.dumps(boundary_check,indent=2))

# Time-block resampling diagnoses nonstationarity; it is not a hardware
# coverage uncertainty and is not adopted as a final correction error.
rng=np.random.default_rng(4398);bootstrap=[]
good=(cur>=lo)&(cur<=hi)
for width in [10,20,60]:
    groups=np.floor((time[1:]-time[0])/width).astype(int)
    values=np.array([np.bincount(groups[good],weights=v[good]) for v in [counts,tightcounts,ds,de]]).T
    values=values[values[:,3]>0]
    sample=values[rng.integers(0,len(values),size=(2000,len(values)))].sum(axis=1)
    ed=sample[:,1]/sample[:,3];cl=(sample[:,0]-sample[:,1])/(sample[:,2]-sample[:,3])
    bootstrap.append(dict(block_seconds=width,blocks=len(values),replicates=2000,seed=4398,
                          EDTM_bootstrap_sigma=float(ed.std(ddof=1)),CLTP_bootstrap_sigma=float(cl.std(ddof=1)),
                          EDTM_percentile95=np.percentile(ed,[2.5,97.5]).tolist(),CLTP_percentile95=np.percentile(cl,[2.5,97.5]).tolist()))
(P/'full_run_block_bootstrap.json').write_text(json.dumps(bootstrap,indent=2))

plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.25})
with PdfPages(P/'full_run_diagnostic_figures.pdf') as pdf:
    fig,axs=plt.subplots(3,1,figsize=(11,8),sharex=True,layout='constrained')
    axs[0].plot(time[1:],cur,lw=.8);axs[0].axhspan(lo,hi,alpha=.1,color='green');axs[0].set_ylabel('Current (uA)')
    axs[1].plot(time[1:],np.cumsum(tightcounts-de),lw=1);axs[1].set_ylabel('Cumulative EDTM\naccepted - sent')
    axs[2].plot(time[1:],np.cumsum(counts-da),lw=1);axs[2].set_ylabel('Cumulative events\nrecorded - L1');axs[2].set_xlabel('Scaler clock (s)')
    for r in segrows[1:]:
        for ax in axs:ax.axvline(r['first_clock'],color='gray',ls=':',lw=.7)
    fig.suptitle('Run 4398: real scaler snapshots joined across six updated segments')
    pdf.savefig(fig);fig.savefig(P/'full_run_alignment.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    q=prod&(raw>1)
    axs[0].hist(tc[q],bins=np.arange(240,290,.1),histtype='step');axs[0].set_yscale('log');axs[0].set_xlabel('EDTM reference-subtracted time (ns)');axs[0].set_ylabel('Events / 0.1 ns');axs[0].axvspan(243.11508,247.11508,color='green',alpha=.1)
    axs[1].plot([q['halfwidth_ns'] for q in timing],[q['EDTM'] for q in timing],'o-');axs[1].set_xscale('log');axs[1].set_xlabel('Timing half-width (ns)');axs[1].set_ylabel('EDTM count ratio')
    fig.suptitle('Fixed production-current window: timing sensitivity')
    pdf.savefig(fig);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.5),layout='constrained')
    ax.plot(range(len(windows)),[r['EDTM'] for r in windows],'o-',label='EDTM, +/-2 ns');ax.plot(range(len(windows)),[r['CLTP'] for r in windows],'o-',label='Physics computer LT')
    ax.set_xticks(range(len(windows)),[f"{r['low']:g}--{r['high']:g}" for r in windows],rotation=25);ax.set_xlabel('Current interval (uA)');ax.set_ylabel('Matched ratio');ax.legend();fig.suptitle('Full covered run: current-cut sensitivity')
    pdf.savefig(fig);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(10,6),sharex=True,layout='constrained')
    axs[0].plot(time[1:][good],(ds-counts)[good],'.',ms=3);axs[0].set_ylabel('Trigger inputs - recorded\nper interval')
    axs[1].plot(time[1:][good],(de-tightcounts)[good],'.',ms=3);axs[1].set_ylabel('EDTM sent - accepted\nper interval');axs[1].set_xlabel('Scaler clock (s)')
    fig.suptitle('Losses vary with time: retain real busy intervals in the average')
    pdf.savefig(fig);fig.savefig(P/'full_run_loss_bursts.png',dpi=150);plt.close(fig)
print(json.dumps(summary,indent=2))
