"""Standalone comparison figures; all values are unclipped diagnostics."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
P=Path(__file__).resolve().parent;F=P/'figures';F.mkdir(exist_ok=True)
R=json.loads((P/'audited_results.json').read_text());valid=[r for r in R if r['status']=='CALCULATED'];prod=[r for r in valid if r['run_type']=='production'];byrun={r['run']:r for r in valid}
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2,'savefig.dpi':180,'axes.titlesize':16,'axes.labelsize':12})
blue='#176b9b';orange='#d66a16';green='#23764a';red='#b3294d';gray='#797e85'
def save(fig,name):
    fig.savefig(F/(name+'.png'));fig.savefig(F/(name+'.pdf'));plt.close(fig)
def ticks(ax,rs,every=1):
    indices=sorted(set(range(0,len(rs),every))|{i for i,r in enumerate(rs) if r['run']==4493})
    ax.set_xticks(indices,[rs[i]['run'] for i in indices],rotation=90)
    ax.tick_params(axis='x',labelsize=10);ax.set_xlabel('Run (ordered positions; gaps in run numbering are compressed)')
def comparison(rs,name,title):
    fig,ax=plt.subplots(figsize=(12.4,5.1),layout='constrained');x=np.arange(len(rs))
    ax.scatter(x-.12,[float(r['saved']['NewGen_EDTM_livetime']) for r in rs],facecolors='none',edgecolors=orange,s=45,label='Saved NewGen')
    ax.scatter(x,[r['old_ratio'] for r in rs],color=orange,marker='x',s=25,label='NewGen recomputed, same cache')
    ax.scatter(x+.12,[r['nominal']['EDTM_tight'] for r in rs],color=blue,s=27,label='Matched EDTM, +/-2 ns')
    ax.scatter(x+.20,[r['nominal']['CLT_physics'] for r in rs],color=green,marker='^',s=28,label='Physics CLT (conditional)')
    for i,r in enumerate(rs):
        if r['prescale_factor']>1:ax.axvspan(i-.45,i+.45,color=gray,alpha=.1)
        if r['closure_exclusion_diagnostic']['removed_intervals']:ax.scatter(i+.12,r['nominal']['EDTM_tight'],s=90,facecolors='none',edgecolors=red,lw=1.5)
    ax.axhline(1,color='black',lw=.8);ax.set_ylabel('Count ratio (no clipping)');ticks(ax,rs,2 if len(rs)>25 else 1)
    ax.set_title(title);ax.legend(fontsize=9,ncol=2,loc='best');save(fig,name)
comparison(prod,'comparison_all','All 43 production runs: saved and cache-matched comparisons')
comparison(prod[:22],'comparison_early','Production runs 4253--4397: changes are not all in one direction')
comparison(prod[22:],'comparison_late','Production runs 4398--4558: preserve real structure and anomalies')
rs=[r for r in prod if r['prescale_factor']==1]
comparison(rs,'comparison_unprescaled','Factor-one runs: endpoint/current effects and residual defects')

fig,axs=plt.subplots(2,1,figsize=(12.4,6.2),sharex=True,layout='constrained');x=np.arange(len(prod))
for ax in axs:ax.axhline(0,color='black',lw=.8)
for key,previous,label,color in [('perfile_covered_oldcurrent_ratio','old_ratio','Restrict to per-file covered events',orange),('perfile_covered_aligned_ratio','perfile_covered_oldcurrent_ratio','Assign ending-interval current',blue)]:
    axs[0].plot(x,[100*(r[key]-r[previous]) for r in prod],'o-',ms=3,lw=.8,label=label,color=color)
axs[1].plot(x,[100*(r['nominal']['EDTM_raw']-r['perfile_covered_aligned_ratio']) for r in prod],'o-',ms=3,label='Join compatible neighboring snapshots',color=green)
axs[1].plot(x,[100*(r['nominal']['EDTM_tight']-r['nominal']['EDTM_raw']) for r in prod],'o-',ms=3,label='Replace raw window with +/-2 ns',color=red)
for ax in axs:ax.set_ylabel('Change (percentage points)');ax.legend(fontsize=10)
ticks(axs[1],prod,2);fig.suptitle('Sequential changes from the original calculation on identical cache files');save(fig,'decomposition')

fig,ax=plt.subplots(figsize=(12.4,4.9),layout='constrained')
delta=np.array([100*(r['old_ratio']/r['nominal']['EDTM_tight']-1) for r in prod])
ax.bar(x,delta,color=[red if r['prescale_factor']>1 or r['closure_exclusion_diagnostic']['removed_intervals'] else blue for r in prod]);ax.axhline(0,color='black',lw=.8)
ticks(ax,prod,2);ax.set_ylabel('Hypothetical yield change (%)');ax.set_title('Fixed yield and charge: changing only L changes the normalization')
ax.text(.02,.97,r'$Y_{new}/Y_{old}-1=L_{old}/L_{new}-1$',transform=ax.transAxes,va='top',fontsize=15)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=red,label='Prescaled or counter-flagged'),Patch(color=blue,label='Other production runs')],fontsize=9,loc='lower left')
save(fig,'normalization_effect')

cat=json.loads((P/'catalog_coverage.json').read_text())['runs'];fig,axs=plt.subplots(1,2,figsize=(12.4,8.2),layout='constrained')
for ax,rr in zip(axs,[cat[:26],cat[26:]]):
    a=np.zeros((len(rr),8))
    for i,r in enumerate(rr):
        for j in r['catalog_segments']:a[i,j]=1
        for j in r['cached_segments']:a[i,j]=2
    ax.imshow(a,aspect='auto',interpolation='none',cmap=ListedColormap(['white','#edc5ce',blue]),vmin=0,vmax=2)
    ax.set_xticks(range(8));ax.set_xlabel('Segment number');ax.set_yticks(range(len(rr)),[str(r['run']) for r in rr],fontsize=10);ax.grid(False)
fig.suptitle('Selected replay coverage: blue = cached; rose = catalogued but absent\nWhite = no entry in the selected-source catalog; no tape staging',fontsize=15);save(fig,'coverage')

fig,axs=plt.subplots(2,1,figsize=(12.4,6),sharex=True,layout='constrained')
axs[0].plot(x,[r['raw_peak_channels'] for r in prod],'o',color=orange);axs[0].set_ylabel('Raw peak (channels)')
axs[1].plot(x,[r['corrected_peak_ns'] for r in prod],'o',color=blue);axs[1].set_ylabel('Corrected peak (ns)');ticks(axs[1],prod,2)
fig.suptitle('Find each run\'s pulse peak; raw channels and corrected ns are different');save(fig,'timing_peaks')

fig,axs=plt.subplots(1,2,figsize=(12.4,5.1),layout='constrained')
for run in [4259,4305,4398,4551]:
    r=byrun[run];h=np.load(P/f'results/run{run}_timing_hist.npz');axs[0].step(h['centers'],h['counts'],where='mid',label=str(run))
    axs[1].plot([t['width_ns'] for t in r['timing']],[100*(t['ratio']-r['nominal']['EDTM_tight']) for t in r['timing']],'o-',label=str(run),ms=4)
axs[0].set_xlim(-12,12);axs[0].set_yscale('log');axs[0].set_xlabel('Corrected EDTM time - run peak (ns)');axs[0].set_ylabel('Candidates / 0.1 ns');axs[0].axvspan(-2,2,color=blue,alpha=.08);axs[0].legend(fontsize=10)
axs[1].set_xscale('log');axs[1].set_xlabel('Timing half-width (ns)');axs[1].set_ylabel('Change from +/-2 ns (percentage points)');axs[1].legend(fontsize=10)
fig.suptitle('Timing is a sensitivity study, not a background-subtraction proof');save(fig,'timing_sensitivity')

fig,axs=plt.subplots(3,1,figsize=(12.4,7),sharex=True,layout='constrained')
axs[0].bar(x,[100*(r['nominal']['EDTM_alternative']-r['nominal']['EDTM_tight']) for r in prod],color=blue);axs[0].set_ylabel('Boundary\nchange (pp)')
axs[1].bar(x,[100*(next(t['ratio'] for t in r['timing'] if t['width_ns']==10)-r['nominal']['EDTM_tight']) for r in prod],color=orange);axs[1].set_ylabel('2 to 10 ns\nchange (pp)')
axs[2].bar(x,[100*(r['stable']['EDTM_tight']-r['nominal']['EDTM_tight']) for r in prod],color=green);axs[2].set_ylabel('Stable current\nchange (pp)')
for ax in axs:ax.axhline(0,color='black',lw=.8)
ticks(axs[2],prod,2);fig.suptitle('Boundary, timing and current variations have different meanings');save(fig,'sensitivity_all')

ps=[r for r in prod if r['prescale_factor']>1];fig,ax=plt.subplots(figsize=(12,5.1),layout='constrained');xx=np.arange(len(ps))
for key,shift,label,color,m in [('EDTM_tight',-.18,'pE/D',blue,'o'),('CLT_all',0,'pN/S',gray,'s'),('CLT_physics',.18,'p(N-E)/(S-D)',green,'^')]:ax.scatter(xx+shift,[r['nominal'][key] for r in ps],label=label,color=color,marker=m,s=65)
ax.axhline(1,color='black',lw=.8);ax.set_xticks(xx,[f"{r['run']}\np={r['prescale_factor']} / TRIG{r['trigger']}" for r in ps]);ax.set_ylabel('Prescale-corrected count ratio');ax.legend(ncol=3,fontsize=11);ax.set_title('Prescaled pulses need representative sampling; these ratios share counts');save(fig,'prescales')

for run in [4303,4305,4398,4551]:
    r=byrun[run];ii=[{k:float(v) for k,v in a.items()} for a in csv.DictReader((P/f'results/run{run}_intervals.csv').open())];g=[a for a in ii if a['nominal']]
    t=np.array([a['t_end'] for a in g]);fig,axs=plt.subplots(2,1,figsize=(12.4,5.8),sharex=True,layout='constrained')
    axs[0].plot(t,[a['S']-r['prescale_factor']*a['N'] for a in g],'.',label='S - pN',color=green,ms=3)
    axs[0].plot(t,[a['N']-a['A'] for a in g],'.',label='N - A',color=red,ms=3);axs[0].set_ylabel('Counts / interval');axs[0].legend(fontsize=10)
    axs[1].plot(t,[a['D']-r['prescale_factor']*a['E_tight'] for a in g],'.',color=blue,ms=3);axs[1].set_ylabel('D - pE / interval');axs[1].set_xlabel('Scaler clock (s; separate components must not be bridged)')
    for ax in axs:ax.axhline(0,color='black',lw=.6)
    fig.suptitle(f'Run {run}: separate real input-to-accept losses from inconsistent counters');save(fig,f'intervals_{run}')

fig,ax=plt.subplots(figsize=(12.4,5.1),layout='constrained')
ax.semilogy(x,[float(r['saved']['NewGen_EDTM_livetime_err']) for r in prod],'o',color=orange,label='Saved independent-Poisson propagation')
ax.semilogy(x,[r['nominal']['EDTM_conditional_sigma'] for r in prod],'s',ms=4,color=green,label='Conditional Bernoulli scale')
ax.semilogy(x,[next(t['EDTM_sigma'] for t in r['block_bootstrap'] if t['block_seconds']==20) for r in prod],'^',ms=4,color=blue,label='20-s block-resampling spread')
ticks(ax,prod,2);ax.set_ylabel('EDTM ratio uncertainty / diagnostic spread');ax.legend(fontsize=9);ax.set_title('Different statistical assumptions; none covers unknown hardware routing');save(fig,'uncertainties')
print('Wrote',len(list(F.glob('*.pdf'))),'standalone figures')
