#!/usr/bin/env python3
"""KinC_x36_4 estimator diagnostics, NOT a corrected estimator or release tool.

Inputs are read-only. Each command creates a fresh output directory. Template
predictions are never substituted for an identified data component. ROOT is
needed only for `toys`, which reuses the production mass fitter and geometry.
"""
import argparse
import contextlib
import csv
import hashlib
import io
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'src/analysis'))
from combine_analysis_branches import _fit_combined_2d_mass_cut
from bootstrap_pi0_data import Fitter, timing_components, MASS_EDGES, DIAMOND, TP_EDGES

DATA = REPO / 'validation/xsec_input_recovery_20261004/output/KinC_x36_4/root/combined_branches_LH2.root'
SIM = Path('/volatile/hallc/nps/singhav/nps_smearing/smear_x36_4/smearing_output/KinC_x36_4/root/simc_pi0_analysis_output_smeared.root')
NOMINAL = REPO / 'validation/model_release_20261004/nominal'
CHANNELS = ('exclusive', 'sidis', 'delta')


def write_csv(path, rows):
    with path.with_suffix('.csv.tmp').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    path.with_suffix('.csv.tmp').replace(path)


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save_json(path, value):
    tmp = path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n'); tmp.replace(path)


def load():
    with uproot.open(DATA) as f:
        e = f['physics'].arrays(library='np')
        runs = f['analysis_runs'].arrays(library='np')
    with uproot.open(SIM) as f:
        s = f['simulation'].arrays(['mpi0','mmiss','Q2','xB','t','tmin','phi','full_weight',
                                    'is_exclusive','is_sidis','is_delta'], library='np')
        provenance = {k.split(';')[0]: f[k].member('fTitle') for k in f.keys()
                      if any(t in k for t in ('normfac','ngen','normalization_source','producer_source_identity'))}
    return e, runs, s, provenance


def rows(e):
    q, x = [e[k].astype('float32').astype(float) for k in ('Q2','xB')]
    tp = e['t'].astype('float32').astype(float) - e['tmin'].astype('float32').astype(float)
    it = np.searchsorted(TP_EDGES, tp, side='right') - 1
    phi = np.mod(e['phi'], 2*np.pi)
    ip = np.floor(np.where(np.isfinite(phi), phi, 0)*12/(2*np.pi)).astype(int)
    cross = np.array([(b[0]-a[0])*(q-a[1])-(b[1]-a[1])*(x-a[0])
                      for a,b in zip(DIAMOND,np.roll(DIAMOND,-1,axis=0))])
    diamond = np.all(cross>=-1e-12,axis=0) | np.all(cross<=1e-12,axis=0)
    ok = diamond & (q>=3.3) & (q<=4.7) & (x>=.29) & (x<=.44) & (it>=0) & (it<5) & np.isfinite(phi)
    return it*12+ip, ok


def geometry():
    path = DATA.with_name(DATA.stem+'_combined_2d_mass_cut_debug.txt')
    d = {}
    for line in path.read_text().splitlines():
        if '=' in line:
            k,v = line.split('=',1)
            try: d[k] = float(v)
            except ValueError: pass
    return d


def ellipse(m, x, d):
    dx, dy = m-d['mean_mpi0'], x-d['mean_mmiss']
    a,b,c = d['cov_mpi0_mpi0'],d['cov_mpi0_mmiss'],d['cov_mmiss_mmiss']
    return ((c*dx*dx-2*b*dx*dy+a*dy*dy)/(a*c-b*b)<=d['ellipse_d2_cut']) & \
        (m>=d['mpi0_min']) & (m<d['mpi0_max']) & (x>=d['mmiss_min']) & (x<d['mmiss_max'])


def counts(row, ok, weight):
    return np.bincount(row[ok], weights=weight[ok], minlength=60)


def figsave(out, name):
    plt.tight_layout(); plt.savefig(out/(name+'.pdf')); plt.savefig(out/(name+'.png'), dpi=140); plt.close()


def audit(out):
    e,runs,s,provenance = load()
    dr, dk = rows(e); sr, sk = rows(s); d = geometry()
    de = e['is_exclusive_ellipse_combined'] != 0
    assert np.array_equal(de, ellipse(e['mpi0_all'],e['mmiss_all'],d))
    se = ellipse(s['mpi0'],s['mmiss'],d)
    charge = float(runs['charge_uC'].astype(float).sum())
    factor = e['scale'].astype('float32').astype(float)*e['charge_uC'].astype('float32').astype(float)/charge/.584
    w = e['pi0_weight']*factor
    c = np.array([1.,-1/6,-1/12,-1/12,1/18,1/18]) @ timing_components(e)
    c *= (e['mpi0_all']>=0) & (e['mpi0_all']<.4)
    pre,sel,timing = counts(dr,dk,w),counts(dr,dk&de,w),counts(dr,dk&de,c*factor)
    old = read_csv(NOMINAL/'model_reconstructed_yields.csv')
    np.testing.assert_allclose(sel,[float(v['data']) for v in old],rtol=2e-12,atol=1e-14)
    templates = {ch:counts(sr,sk&se&(s['is_'+ch]!=0),s['full_weight']) for ch in CHANNELS}
    before = {ch:counts(sr,sk&(s['is_'+ch]!=0),s['full_weight']) for ch in CHANNELS}
    # Verify exclusive MC membership against the actual model event cache.
    cache = np.genfromtxt(NOMINAL/'model_event_cache.csv',delimiter=',',names=True)
    support = np.bincount(cache['row'].astype(int),minlength=60)
    np.testing.assert_array_equal(counts(sr,sk&se&(s['is_exclusive']!=0),np.ones(len(sr))),support)
    composition = []
    for i in range(60):
        total = sum(t[i] for t in templates.values())
        composition.append(dict(row=i,tprime_bin=i//12,phi_bin=i%12,
            total_pi0_signal=pre[i],ellipse_selected_pi0=sel[i],exclusive_component=np.nan,
            SIDIS_component=np.nan,SIDIS_fraction=np.nan,combinatorial_residual=np.nan,
            delta_component=np.nan,raw_selected_events=int(np.sum(dk&de&(dr==i))),
            timing_selected_before_combinatorial=timing[i],
            template_exclusive=templates['exclusive'][i],template_SIDIS=templates['sidis'][i],
            template_delta=templates['delta'][i],
            template_SIDIS_fraction=templates['sidis'][i]/total if total else np.nan,
            statistical_source='observed mass-weight sum; frozen input geometry; SIMC full_weight',
            estimation_method='UNVALIDATED mass transport; template fractions are model-only, not fitted data composition',
            units='per_mC; data includes target divisor 0.584'))
    write_csv(out/'exclusive_composition_by_row.csv',composition)
    summary = []
    for it in range(5):
        sl=slice(12*it,12*(it+1)); t={ch:float(v[sl].sum()) for ch,v in templates.items()}; total=sum(t.values())
        summary.append(dict(tprime_bin=it,tprime_lo=TP_EDGES[it],tprime_hi=TP_EDGES[it+1],
            total_pi0_signal=float(pre[sl].sum()),ellipse_selected_pi0=float(sel[sl].sum()),
            SIDIS_fraction=np.nan,exclusive_component=np.nan,SIDIS_component=np.nan,
            template_exclusive=t['exclusive'],template_SIDIS=t['sidis'],template_delta=t['delta'],
            template_SIDIS_fraction=t['sidis']/total,
            template_nonexclusive_fraction=(t['sidis']+t['delta'])/total,
            template_SIDIS_over_observed=t['sidis']/sel[sl].sum(),
            status='DATA_COMPOSITION_UNIDENTIFIED; absolute generator prediction only'))
    write_csv(out/'sidis_summary_by_tprime.csv',summary)
    zero = []
    for run in runs['run_number']:
        status=read_csv(DATA.parent/f'analysis_status_run{run}.csv')[0]
        if status['classification']=='zero_background':zero.append(run)
    z=np.isin(e['run_number'],zero)&dk&de
    transport = dict(zero_background_runs=len(zero),mass_weight=float(w[z].sum()),timing=float((c*factor)[z].sum()))
    transport['difference']=transport['timing']-transport['mass_weight']
    transport['percent_relative_mass_weight']=100*transport['difference']/transport['mass_weight']
    transport['interpretation']='Same selected events, different timing redistribution; no exclusive/SIDIS labels. Neither is known exclusive truth.'
    save_json(out/'weight_transport.json',transport)
    write_csv(out/'weight_transport_by_row.csv',[dict(row=i,mass_weight=sel[i],timing_before_comb=timing[i],difference=timing[i]-sel[i]) for i in range(60)])
    variance_rows=[]
    for i in (0,5,6,11):
        variance_rows.append(dict(row=i,observed_events_before_ellipse=int(np.sum(dk&(dr==i))),
            observed_events=int(np.sum(dk&de&(dr==i))),pi0_signal_estimate=pre[i],
            ellipse_selected_pi0=sel[i],SIDIS_contribution=np.nan,exclusive_yield=np.nan,
            bootstrap_variance=np.nan,old_conditional_variance=float(old[i]['data_sumw2']),
            MC_response_support=int(support[i]),final_decision='PENDING; no data-dependent exclusion authorized',
            reason='Full exclusive estimator unavailable; empirical fixed-geometry resampling cannot populate an empty observed row'))
    write_csv(out/'corrected_row_variances.csv',variance_rows)
    central=read_csv(NOMINAL/'model_structure_functions.csv')
    q=float(read_csv(NOMINAL/'model_fit_status.csv')[0]['objective'])
    write_csv(out/'central_M0_before_after.csv',[dict(tprime_bin=i,old_Q=q,corrected_Q=np.nan,
        **{'old_'+k:float(v['sigma_'+k])*1e9 for k in ('U','LT','TT')},
        corrected_U=np.nan,corrected_LT=np.nan,corrected_TT=np.nan,
        shift_U=np.nan,shift_LT=np.nan,shift_TT=np.nan,units='nb/GeV2',
        status='BLOCKED: no validated exclusive detector yield') for i,v in enumerate(central)])
    save_json(out/'input_provenance.json',dict(data=str(DATA),simulation=str(SIM),charge_mC=charge/1000,
        simulation_channels={ch:int(np.sum(s['is_'+ch])) for ch in CHANNELS},simulation_metadata=provenance,
        target_divisor=.584,geometry=d,nominal_parity='passed data yields and exclusive simulation row counts',
        mass_distribution='200 bins [0,0.4) GeV per run; residual/all-event count',
        templates='full_weight=Weight*channel_normfac/channel_Ngen; no data normalization fitted'))
    fig,axs=plt.subplots(1,3,figsize=(13,3.8))
    for ax,it in zip(axs,(0,2,4)):
        take=dk&(dr//12==it)
        ax.hist(e['mpi0_all'][take],bins=np.linspace(0,.3,101),weights=(c*factor)[take],histtype='step',label='Timing-subtracted candidates')
        ax.hist(e['mpi0_all'][take],bins=np.linspace(0,.3,101),weights=w[take],histtype='step',label='Mass pi0 estimate')
        ax.set(title=f"t' [{TP_EDGES[it]}, {TP_EDGES[it+1]}]",xlabel='mpi0 [GeV]',ylabel='Yield / mC / bin');ax.legend(fontsize=7)
    figsave(out,'mpi0_before_exclusivity')
    fig,axs=plt.subplots(1,3,figsize=(13,3.8))
    for ax,it in zip(axs,(0,2,4)):
        take=dk&(dr//12==it);bins=np.linspace(.4,2.,65)
        ax.hist(e['mmiss_all'][take],bins=bins,weights=w[take],histtype='step',color='k',label='Data mass-weight estimate')
        for ch in CHANNELS:
            take=sk&(sr//12==it)&(s['is_'+ch]!=0)
            ax.hist(s['mmiss'][take],bins=bins,weights=s['full_weight'][take],histtype='step',label=ch+' MC (unfitted)')
        ax.set(title=f"Reco t' bin {it}",xlabel='Missing mass [GeV]',ylabel='Yield / mC / bin');ax.legend(fontsize=7)
    figsave(out,'missing_mass_templates')
    plt.figure(figsize=(7,4))
    plt.plot(range(5),[r['template_SIDIS_fraction'] for r in summary],'o-',label='SIDIS / total MC')
    plt.plot(range(5),[r['template_nonexclusive_fraction'] for r in summary],'s-',label='(SIDIS + delta) / total MC')
    plt.title('Ellipse-selected template fractions; data fractions UNKNOWN');plt.xlabel("Reco t' bin");plt.ylabel('Fraction');plt.legend()
    figsave(out,'sidis_fraction_vs_tprime')
    plt.figure(figsize=(8,3))
    plt.bar([str(v['row']) for v in variance_rows],[v['MC_response_support'] for v in variance_rows])
    plt.title('Supported rows: corrected bootstrap NOT RUN');plt.xlabel('Row');plt.ylabel('Exclusive MC events')
    figsave(out,'formerly_zero_variance_rows')
    plt.figure(figsize=(12,3));plt.axis('off')
    labels=['Selected candidates','Mass pi0 estimate\nexclusive + SIDIS + delta','Data-fitted ellipse\nexclusive-enriched','Exclusive estimate\nBLOCKED: residual backgrounds']
    for j,label in enumerate(labels):
        plt.text(.12+.255*j,.5,label,ha='center',va='center',bbox=dict(boxstyle='round',facecolor='lightgray'),fontsize=9,transform=plt.gca().transAxes)
        if j<3:plt.annotate('',xy=(.255+.255*j,.5),xytext=(.22+.255*j,.5),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'))
    figsave(out,'event_estimator_flow')
    print(json.dumps(dict(transport=transport,tprime=summary),indent=2))


def transport_toys(out, seed):
    """Exact zero-combinatorial production branch, one mass bin, two regions.

    A prompt pion has c=1. Six diagonal accidental windows have c=-1/6.
    Different ellipse acceptances model timing/kinematic correlations. This
    deliberately has NO SIDIS; any difference cannot be source composition.
    """
    rng=np.random.default_rng(seed);table=[]
    for name,sa,aa in [('independent_selection',.5,.5),('correlated_selection',.8,.2)]:
        means=np.array([1000*sa,1000*(1-sa),500*aa,500*(1-aa),3000*aa,3000*(1-aa)])
        n=rng.poisson(means,size=(10000,6));c=np.array([1,1,1,1,-1/6,-1/6])
        p=(n@c)/n.sum(axis=1);selected=n[:,[0,2,4]].sum(axis=1)
        mass=p*selected;direct=n[:,0]+n[:,2]-n[:,4]/6;truth=n[:,0]
        for method,recovered in [('mass_redistribution',mass),('event_timing',direct)]:
            residual=recovered-truth
            table.append(dict(case=name,method=method,injected_selected_mean=float(truth.mean()),
                recovered_mean=float(recovered.mean()),bias=float(residual.mean()),
                residual_sd=float(residual.std(ddof=1)),bias_mc_se=float(residual.std(ddof=1)/100),
                pseudoexperiments=10000,SIDIS_injected=0,seed=seed,
                scope='exact zero-combinatorial branch, fixed selection; not full data estimator validation'))
    write_csv(out/'transport_controlled_toys.csv',table)


def toys(out, library, seed, experiments):
    """Full current mass/geometry chain on labelled component pseudo-events.

    56 synthetic runs; known prompt components; no timing accidentals here
    (separate exact transport toys test that failure). Shapes/relative true-pi0
    rates use existing MC. Combinatorial 10% is a declared stress assumption,
    not a measured rate. No SIDIS correction exists in the current chain.
    """
    fit=Fitter(library);e,runs,s,_=load();rng=np.random.default_rng(seed)
    coeff=np.array([1.,-1/6,-1/12,-1/12,1/18,1/18])@timing_components(e)
    signal_count=max(1.,float(coeff.sum()))
    run_ids=runs['run_number'];prob=np.array([np.sum(e['run_number']==r) for r in run_ids],float);prob/=prob.sum()
    charge=float(runs['charge_uC'].astype(float).sum())
    scales=runs['scale'].astype(float);fac=scales*runs['charge_uC'].astype(float)/charge/.584
    pool={};rate={}
    for ch in ('exclusive','sidis'):
        take=np.flatnonzero((s['is_'+ch]!=0)&np.isfinite(s['full_weight'])&(s['full_weight']>0)&(s['mpi0']>=0)&(s['mpi0']<.4))
        weights=s['full_weight'][take].astype(float);pool[ch]=(take,weights/weights.sum());rate[ch]=weights.sum()
    ratio=rate['sidis']/rate['exclusive']; n_excl=signal_count/(1+ratio)
    raw=[]; failures=[]
    for case,use_sidis,use_comb,variation in [(1,0,0,1.),(2,0,1,1.),(3,1,0,1.),(4,1,1,1.),(5,1,1,1.2)]:
        for rep in range(experiments):
            parts=[]
            for label,ch,nmean in [(0,'exclusive',n_excl),(1,'sidis',n_excl*ratio*use_sidis*variation)]:
                idx=rng.choice(pool[ch][0],size=rng.poisson(nmean),p=pool[ch][1])
                part={k:s[k][idx].copy() for k in ('Q2','xB','t','tmin','phi')}
                part.update(mpi0_all=s['mpi0'][idx].copy(),mmiss_all=s['mmiss'][idx].copy(),label=np.full(len(idx),label))
                if case==5 and label==1:part['mmiss_all']+=.03
                parts.append(part)
            nb=rng.poisson(.1*n_excl*use_comb)
            idx=rng.integers(len(e['mpi0_all']),size=nb)
            part={k:e[k][idx].copy() for k in ('Q2','xB','t','tmin','phi')}
            centers=(MASS_EDGES[:-1]+MASS_EDGES[1:])/2
            cb=1/(1+np.exp((centers-.145)/.025));cb/=cb.sum()
            part.update(mpi0_all=rng.choice(centers,size=nb,p=cb)+rng.uniform(-.001,.001,nb),
                mmiss_all=rng.uniform(.6,1.5,nb),label=np.full(nb,2));parts.append(part)
            a={k:np.concatenate([v[k] for v in parts]) for k in parts[0]}
            n=len(a['label']);run=rng.choice(len(run_ids),size=n,p=prob)
            mb=np.searchsorted(MASS_EDGES,a['mpi0_all'],side='right')-1
            w=np.zeros(n);failed=False
            for ir in range(len(run_ids)):
                take=run==ir;hist=np.bincount(mb[take],minlength=200).astype(float)
                code,final,info=fit.fit(hist,hist)
                if code:
                    failures.append(dict(case=case,replica=rep,run=int(run_ids[ir]),stage='mass_fit',status=code));failed=True
                else:w[take]=np.divide(final,hist,out=np.zeros(200),where=hist>0)[mb[take]]
            if failed:continue
            df=pd.DataFrame(dict(mpi0_all=a['mpi0_all'],mmiss_all=a['mmiss_all'],pi0_weight=w,scale=scales[run]))
            with contextlib.redirect_stdout(io.StringIO()):debug=_fit_combined_2d_mass_cut(df)
            if debug is None or not debug['params'].get('ellipse_valid',0):
                failures.append(dict(case=case,replica=rep,run=-1,stage='ellipse_fit',status=1));continue
            rr,ok=rows(a);ok &= df['is_exclusive_ellipse_combined'].to_numpy()!=0
            truth=counts(rr,ok&(a['label']==0),fac[run]);est=counts(rr,ok,w*fac[run])
            sid=counts(rr,ok&(a['label']==1),fac[run]);comb=counts(rr,ok&(a['label']==2),fac[run])
            for i in range(60):raw.append(dict(case=case,replica=rep,row=i,tprime_bin=i//12,
                injected_exclusive=truth[i],current_chain_output=est[i],bias=est[i]-truth[i],
                injected_SIDIS=sid[i],injected_combinatorial=comb[i]))
        print(f'Toy {case}: {experiments} requested; completed {len(set(v["replica"] for v in raw if v["case"]==case))}',flush=True)
    write_csv(out/'toy_failures.csv',failures or [dict(case=-1,replica=-1,run=-1,stage='none',status=0)])
    summary=[]
    for case in range(1,6):
        for scope,groups in [('row',range(60)),('tprime',range(5))]:
            for group in groups:
                selected=[v for v in raw if v['case']==case and v['row' if scope=='row' else 'tprime_bin']==group]
                replicas=sorted(set(v['replica'] for v in selected))
                arr=np.array([[sum(v[k] for v in selected if v['replica']==rep) for k in
                    ('injected_exclusive','current_chain_output','bias','injected_SIDIS')] for rep in replicas])
                vals=arr.mean(axis=0) if len(arr) else [np.nan]*4
                sd=float(arr[:,2].std(ddof=1)) if len(arr)>1 else np.nan
                summary.append(dict(case=case,scope=scope,index=group,injected_exclusive=vals[0],
                    recovered_current_chain=vals[1],bias=vals[2],residual_spread=sd,
                    bias_mc_se=sd/np.sqrt(len(arr)) if len(arr) else np.nan,injected_selected_SIDIS=vals[3],
                    requested=experiments,successful=len(arr),corrected_exclusive_yield=np.nan,
                    status='CURRENT_CHAIN_DIAGNOSTIC_ONLY; no validated SIDIS correction',units='per_mC'))
    write_csv(out/'estimator_closure.csv',summary)
    if raw:write_csv(out/'toy_replica_rows.csv',raw)
    save_json(out/'toy_configuration.json',dict(seed=seed,experiments_per_case=experiments,
        runs=len(run_ids),mean_exclusive_before_ellipse=n_excl,SIDIS_exclusive_ratio=ratio,
        combinatorial_fraction_of_exclusive=.1,combinatorial_mass='Fermi turn-off mt=.145,width=.025 GeV',
        SIDIS_variation='case 5 normalization x1.2 and mmiss +.03 GeV',
        timing='prompt-only; independent transport_controlled_toys test correlated accidentals',
        geometry='production combined geometry refitted each pseudoexperiment',
        mass='production C++ fitter separately for every synthetic run',
        limitations='MC channel labels are generator origin, not photon ancestry; rates not fitted to data; failed replicas are reported, never validation passes'))
    transport_toys(out,seed+1)
    plt.figure(figsize=(10,4))
    for case in range(1,6):
        v=[r for r in summary if r['case']==case and r['scope']=='tprime']
        plt.errorbar(range(5),[r['bias'] for r in v],yerr=[r['bias_mc_se'] for r in v],marker='o',label=f'Toy {case}')
    plt.axhline(0,color='k',lw=.7);plt.xlabel("Reco t' bin");plt.ylabel('Current-chain minus injected exclusive [/mC]');plt.legend();plt.title('Diagnostic closure; bars = Monte Carlo error on bias')
    figsave(out,'exclusive_yield_closure_bias')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['audit','toys','central'])
    p.add_argument('--output',type=Path,required=True);p.add_argument('--library',type=Path)
    p.add_argument('--seed',type=int,default=20261004);p.add_argument('--experiments',type=int,default=12)
    a=p.parse_args()
    if a.mode=='central':
        p.exit(2,'EXCLUSIVE-YIELD ESTIMATOR NOT VALIDATED: no conditional combinatorial estimator or constrained residual SIDIS/delta contribution. No production command is enabled.\n')
    if a.mode=='toys' and not a.library:p.error('--library is required for production-fitter toys')
    a.output.mkdir(parents=True,exist_ok=False)
    if a.mode=='audit':audit(a.output)
    else:toys(a.output,a.library,a.seed,a.experiments)


if __name__=='__main__':main()
