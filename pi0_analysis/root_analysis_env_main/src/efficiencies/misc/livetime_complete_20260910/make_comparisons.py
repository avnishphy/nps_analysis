"""Regenerate comparisons from frozen summaries/intervals, without ROOT reads."""
from pathlib import Path
import csv,json,os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parent/'build/mpl'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;F=P/'figures';B=P/'evidence/full_cache'
R=[json.loads(f.read_text()) for f in sorted((B/'results').glob('run*.json'))]
R=[r for r in R if r.get('run_type')=='production' and r.get('status')=='CALCULATED']
assert len(R)==43
A={(r['run'],r['variant']):r for r in json.loads((B/'acceptance_tests.json').read_text())}
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.17,'savefig.dpi':170,'pdf.fonttype':42})
colors={'old':'#b56517','raw':'#007d86','tight':'#686e76','phys':'#6650a3','flag':'#b62f46'}
def save(fig,name):
 for ext in ['pdf','png']:fig.savefig(F/f'{name}.{ext}',bbox_inches='tight')
 plt.close(fig)
def ticks(ax,rs):ax.set_xticks(range(len(rs)),[str(r['run']) for r in rs],rotation=60,ha='right',fontsize=10)
def compare(rs,name,title):
 fig,ax=plt.subplots(figsize=(12.4,4.8),layout='constrained');x=np.arange(len(rs))
 for offset,vals,label,color,marker in [(-.23,[r['old_ratio'] for r in rs],'Original method, same files','old','x'),(-.08,[r['nominal']['EDTM_raw'] for r in rs],'Matched broad raw EDTM','raw','o'),(.08,[r['nominal']['EDTM_tight'] for r in rs],'Matched tight (+/-2 ns), diagnostic','tight','+'),(.23,[A[r['run'],'raw']['physics_CLT'] for r in rs],'p(N-Eraw)/(S-D), proposed ratio','phys','^')]:
  ax.scatter(x+offset,vals,label=label,color=colors[color],marker=marker,s=35)
 for i,r in enumerate(rs):
  if A[r['run'],'raw']['counter_suspect_intervals']:ax.axvspan(i-.43,i+.43,color=colors['flag'],alpha=.12)
 ax.axhline(1,color='black',lw=.8);ticks(ax,rs);ax.set_ylabel('Count ratio');ax.set_xlabel('Run (ordered positions)');ax.set_title(title,loc='left',fontsize=15)
 ax.legend(fontsize=9,ncol=2,loc='best');save(fig,name)
compare([r for r in R if r['trigger']==6],'cohort_ti6','TI6 coincidence cohort: five runs; 4259 has p=5')
ti4=[r for r in R if r['trigger']==4 and r['prescale_factor']==1]
compare(ti4[:17],'cohort_ti4_early','TI4 HMS EL-REAL singles, p=1: early runs; shaded = counter defect')
compare(ti4[17:],'cohort_ti4_late','TI4 HMS EL-REAL singles, p=1: later runs')
compare([r for r in R if r['trigger']==4 and r['prescale_factor']>1],'cohort_ti4_prescaled','TI4 HMS EL-REAL singles: five runs with p=2')
fig,ax=plt.subplots(figsize=(13,5.3),layout='constrained');x=np.arange(len(R));bottomp=np.zeros(len(R));bottomn=bottomp.copy()
arr=np.array([[r['old_ratio'],r['perfile_covered_oldcurrent_ratio'],r['perfile_covered_aligned_ratio'],r['nominal']['EDTM_raw']] for r in R])
diff=100*np.diff(arr,axis=1)
for j,(label,color) in enumerate([('Restrict to paired per-file coverage','#b56517'),('Align event current to ending snapshot','#007d86'),('Join compatible segments / full coverage','#6650a3')]):
 v=diff[:,j];base=np.where(v>=0,bottomp,bottomn);ax.bar(x,v,bottom=base,color=color,label=label,width=.78);bottomp+=np.maximum(v,0);bottomn+=np.minimum(v,0)
ax.scatter(x,100*(arr[:,-1]-arr[:,0]),color='black',s=10,label='Net change');ax.axhline(0,color='black',lw=.8);ticks(ax,R);ax.set_ylabel('Change in livetime ratio (percentage points)');ax.set_title('Actual refinements to the existing pE/D estimator',loc='left');ax.legend(fontsize=9,ncol=2);save(fig,'method_decomposition_current')
fig,ax=plt.subplots(figsize=(11.5,4.7),layout='constrained');r=next(r for r in R if r['run']==4398)
vals=[r['old_ratio'],r['perfile_covered_oldcurrent_ratio'],r['perfile_covered_aligned_ratio'],r['nominal']['EDTM_raw']]
labels=['Original, same six files','Per-file coverage only','Ending-current aligned','Final stitched diagnostic']
ax.plot(range(4),vals,'o-',color=colors['raw'],lw=2);ax.set_xticks(range(4),labels,fontsize=10);ax.set_ylabel('pE/D');ax.set_title('4398: E changes; D stays at 112,746',loc='left');ax.set_ylim(min(vals)-.00065,max(vals)+.00065)
for i,(v,n) in enumerate(zip(vals,[112586,112263,112658,112658])):ax.annotate(f'E={n:,}\nL={v:.8f}',(i,v),xytext=(0,13),textcoords='offset points',ha='center',fontsize=10)
save(fig,'run4398_refinement_steps')
fig,axs=plt.subplots(1,2,figsize=(12.2,5),layout='constrained')
for t,marker,c in [(4,'o',colors['raw']),(6,'s',colors['phys'])]:
 rs=[r for r in R if r['trigger']==t and r['prescale_factor']==1];xx=[r['nominal']['S']/r['nominal']['dt']/1000 for r in rs]
 axs[0].scatter(xx,[r['nominal']['EDTM_raw'] for r in rs],marker=marker,color=c,label=f'TI{t}, p=1',s=35)
 axs[1].errorbar(xx,[100*A[r['run'],'raw']['difference'] for r in rs],yerr=[100*A[r['run'],'raw']['block20_difference_sigma'] for r in rs],fmt=marker,color=c,ms=4,capsize=2,label=f'TI{t}, p=1')
 for r,rate in zip(rs,xx):
  if r['run'] in [4303,4305]:axs[0].annotate(str(r['run']),(rate,r['nominal']['EDTM_raw']),fontsize=9,color=colors['flag'])
for ax in axs:ax.set_xlabel('Selected-input scaler rate S / exposure (kHz)');ax.legend(fontsize=10)
axs[0].set_ylabel('Matched broad EDTM ratio');axs[1].set_ylabel('EDTM minus proposed physics ratio (pp)');axs[0].axhline(1,color='black',lw=.8);axs[1].axhline(0,color='black',lw=.8)
fig.suptitle('Observed rate associations only: no fit and no correction inferred\nDefective scaler exposure can affect both rate and livetime.',fontsize=13);save(fig,'rate_comparison_by_trigger')
for run in [4303,4305,4398,4551]:
 rs=list(csv.DictReader((B/'results'/f'run{run}_intervals.csv').open()));r=next(r for r in R if r['run']==run);p=r['prescale_factor']
 rs=[q for q in rs if int(q['nominal'])];tt=np.array([float(q['t_end']) for q in rs]);fig,axs=plt.subplots(3,1,figsize=(11.7,5.8),sharex=True,layout='constrained')
 ys=[np.array([float(q['S'])-p*float(q['N']) for q in rs]),np.array([float(q['N'])-float(q['A']) for q in rs]),np.array([float(q['D'])-p*float(q['E_raw']) for q in rs])]
 for ax,y,label,c in zip(axs,ys,['S - pN','N - A','D - pEraw'],[colors['phys'],colors['flag'],colors['raw']]):ax.plot(tt,y,'.',ms=2,color=c);ax.set_ylabel(label);ax.axhline(0,color='black',lw=.6)
 axs[-1].set_xlabel('Ending scaler clock (s; selected intervals only)');fig.suptitle(f'Run {run}, p={p}: distinguish lost accepts from inconsistent counters',fontsize=14);save(fig,f'current_intervals_{run}')
(P/'data/method_decomposition.json').write_text(json.dumps([dict(run=r['run'],coverage_pp=float(v[0]),current_pp=float(v[1]),joins_pp=float(v[2]),net_pp=float(v.sum())) for r,v in zip(R,diff)],indent=2)+'\n')
print('Created 11 new comparison figures; no ROOT/event reread.')
