"""Focused diagnostic figures for the refreshed sample; no presentation edit."""
from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;F=P/'figures';F.mkdir(exist_ok=True)
R=[r for r in json.loads((P/'audited_results.json').read_text()) if r['status']=='CALCULATED' and r['run_type']=='production']
B=json.loads((P/'phase_timing_bridge.json').read_text());A={r['run']:r for r in json.loads((P/'acceptance_tests.json').read_text()) if r['variant']=='raw'}
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18,'savefig.dpi':180})
blue='#176b9b';orange='#cf6a18';green='#23764a';red='#b3294d';gray='#71777e'
def save(fig,name):
    fig.savefig(F/(name+'.png'));fig.savefig(F/(name+'.pdf'));plt.close(fig)
def ticks(ax,rs):
    ax.set_xticks(range(len(rs)),[str(r['run']) for r in rs],rotation=90,fontsize=9)
    ax.set_xlabel('Run (ordered positions)')
fig,axs=plt.subplots(2,1,figsize=(13,8),layout='constrained',gridspec_kw={'height_ratios':[1.35,1]})
for ax,rs,title in zip(axs,[[r for r in R if r['prescale_factor']==1],[r for r in R if r['prescale_factor']>1]],['Unprescaled production runs','Prescaled production runs: EDTM fluctuations need their own uncertainty']):
    x=np.arange(len(rs))
    for dx,key,label,col,mark in [(-.18,'EDTM_raw','Matched raw EDTM',orange,'o'),(-.06,'EDTM_tight','Corrected +/-2 ns diagnostic',blue,'x'),(.06,'CLT_all','pN/S',gray,'s')]:
        ax.scatter(x+dx,[r['nominal'][key] for r in rs],s=30,label=label,color=col,marker=mark)
    ax.scatter(x+.18,[A[r['run']]['physics_CLT'] for r in rs],s=30,label='p(N-Eraw)/(S-D), conditional',color=green,marker='^')
    for i,r in enumerate(rs):
        if r['closure_exclusion_diagnostic']['removed_intervals']:ax.axvspan(i-.42,i+.42,color=red,alpha=.13)
    ax.axhline(1,color='black',lw=.8);ticks(ax,rs);ax.set_ylabel('Count ratio; no clipping');ax.set_title(title,loc='left')
axs[0].legend(ncol=2,fontsize=9);fig.suptitle('All 43 production runs / 183 cached updated segments\nShared exposure; red bands flag counter defects. These are not certified total livetimes.',fontsize=14)
save(fig,'refresh_livetime_comparison')

fig,axs=plt.subplots(1,2,figsize=(12.5,5.4),layout='constrained')
for ax,key,title in zip(axs,['corrected_delta_ns','relative_delta_ns'],['Corrected EDTM time','EDTM time minus selected-trigger time']):
    for name,label,col,alpha,size in [('core_heldout','Core: held-out, thinned',gray,.15,5),('raw_tail','Raw-window tails: all',orange,.8,18),('positive_outside_raw','Positive hits outside raw window',blue,.35,9)]:
        a=[r for r in B['points'] if r['category']==name]
        ax.scatter([r['phase_ticks'] for r in a],[r[key] for r in a],s=size,alpha=alpha,color=col,label=label,rasterized=True)
    ax.set_xlim(-500,500);ax.set_xlabel('Event-clock residual from predicted pulse (ticks)');ax.set_ylabel('Time minus core center (ns)');ax.set_title(title,loc='left')
    ax.plot([-500,500],[2000,-2000],'k--',lw=.8,label='-4 ns/tick reference')
axs[0].legend(fontsize=8);fig.suptitle('Pulse association does not identify which signal caused the trigger\nTiming/phase correlation is consistent with nearby pulses entering a physics-event TDC window.',fontsize=14)
save(fig,'refresh_phase_timing')

fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
tail=[r for r in B['points'] if r['category']=='raw_tail']
ax.hist([r['phase_ticks'] for r in tail],bins=np.arange(-30,31,2),color=orange,edgecolor='white')
ax.axvspan(-25,25,color=blue,alpha=.07);ax.set_xlabel('Event-clock residual from predicted pulse (ticks)');ax.set_ylabel('Raw-window candidates rejected by corrected +/-2 ns')
ax.set_title('Independent timestamp test of rejected EDTM tags',loc='left')
ax.text(.02,.97,f'{len(tail)} pulse-associated tails\nAnchors fitted from alternate core events; no extrapolation',va='top',transform=ax.transAxes)
save(fig,'refresh_rejected_tags')

fig,axs=plt.subplots(2,1,figsize=(12.7,7),layout='constrained')
for ax,rs in zip(axs,[[r for r in R if r['prescale_factor']==1],[r for r in R if r['prescale_factor']>1]]):
    x=np.arange(len(rs));a=[A[r['run']] for r in rs]
    ax.errorbar(x,[100*q['difference'] for q in a],yerr=[100*q['block20_difference_sigma'] for q in a],fmt='o',color=blue,capsize=3,ms=4)
    for i,q in enumerate(a):
        if q['counter_suspect_intervals']:ax.axvspan(i-.42,i+.42,color=red,alpha=.13)
    ax.axhline(0,color='black',lw=.8);ticks(ax,rs);ax.set_ylabel('EDTM - physics ratio\n(percentage points)')
fig.suptitle('Compare two proposed count ratios; neither is a livetime baseline\nError bars: paired 20-second block spread. Agreement does not validate either definition.',fontsize=14)
save(fig,'refresh_correlated_difference')

C=json.loads((P/'event_clock_checks.json').read_text())['disagreements']
fig,axs=plt.subplots(2,2,figsize=(11.5,7),layout='constrained')
for row,c in zip(axs,C):
    ax=row[0];v=[c[k] for k in ['N','A','E_raw','D']];bars=ax.bar(range(4),v,color=[blue,gray,orange,green]);ax.set_xticks(range(4),['Recorded N','L1 A','EDTM tags E','Sent D']);ax.bar_label(bars,padding=3,fontsize=10);ax.set_ylim(0,max(v)*1.2);ax.set_ylabel('Counts');ax.set_title(f"Run {c['run']}: events {c['ev_start']} to {c['ev_end']}",loc='left')
    ax=row[1];v=[c['event_seconds'],c['scaler_seconds']];bars=ax.bar(range(2),v,color=[blue,gray]);ax.set_xticks(range(2),['Event clock','Scaler clock']);ax.bar_label(bars,labels=[f'{x:.6f} s' for x in v],padding=3);ax.set_ylim(0,2.5);ax.set_ylabel('Exposure (s)')
fig.suptitle('Two intervals lose almost 2 seconds of scaler exposure\nExact event boundaries. Counter/clock inconsistency is established; its cause is not.',fontsize=14)
save(fig,'refresh_counter_exposure')

fig,axs=plt.subplots(2,1,figsize=(12.7,7),sharex=True,layout='constrained');x=np.arange(len(R))
axs[0].bar(x,[100*(r['old_ratio']/float(r['saved']['NewGen_EDTM_livetime'])-1) for r in R],color=gray);axs[0].set_ylabel('Ratio change (%)');axs[0].set_title('Recompute the original method: refreshed files versus saved NewGen',loc='left')
axs[1].bar(x,[100*(r['old_ratio']/r['nominal']['EDTM_raw']-1) for r in R],color=blue);axs[1].set_ylabel('Hypothetical normalization\nchange (%)');axs[1].set_title('Same refreshed files: old method to matched raw EDTM (Lold/Lmatched - 1)',loc='left');ticks(axs[1],R)
for ax in axs:ax.axhline(0,color='black',lw=.8)
fig.suptitle('Separate input/replay changes from calculation changes\nLower panel holds yield and charge fixed; it is not a measured cross-section bias.',fontsize=14)
save(fig,'refresh_normalization_sensitivity')
print('Created',len(list(F.glob('refresh_*.png'))),'focused diagnostic figures')
