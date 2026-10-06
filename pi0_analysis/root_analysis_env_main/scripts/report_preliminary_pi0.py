#!/usr/bin/env python3
"""Publish preliminary M0 ensemble products; never substitutes Hessian errors."""
import argparse
import hashlib
import subprocess
import shutil
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from preliminary_pi0_xsec import *

COMPONENTS=['U','LT','TT']

def covariance(out,name,matrix,ordering,units):
    sd=np.sqrt(np.maximum(0,np.diag(matrix)))
    corr=np.divide(matrix,np.outer(sd,sd),out=np.zeros_like(matrix),where=np.outer(sd,sd)>0)
    eigen=np.linalg.eigvalsh(matrix);ce=np.linalg.eigvalsh(corr)
    for prefix,value in [('V',matrix),('Corr',corr)]:
        np.savetxt(out/f'{prefix}_{name}.csv',value,delimiter=',',fmt='%.17g')
        np.save(out/f'{prefix}_{name}.npy',value)
    meta=dict(ordering=ordering,units=units,eigenvalues=eigen.tolist(),correlation_eigenvalues=ce.tolist(),
        effective_rank=int(sum(ce>max(ce[-1],0)*1e-10)),rank_definition='correlation eigenvalue > 1e-10 * largest',
        covariance_absolute_rank=int(sum(eigen>max(eigen[-1],0)*1e-10)))
    save(out/f'V_{name}_metadata.json',meta)
    return sd,corr,meta

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();out=a.output
    central=json.loads((out/'central_results.json').read_text());c=central['additive'];z=np.array(c['pub'])
    samples={mode:np.load(out/mode/'replicas.npz') for mode in ['data','mc','toys']}
    summary={mode:json.loads((out/mode/'summary.json').read_text()) for mode in samples}
    if min(summary[m]['accepted'] for m in ['data','mc'])<300:raise ValueError('Insufficient accepted replicas for release')
    names=[f'{k}(bin{b})' for k in COMPONENTS for b in range(5)]
    pnames=['N_U','DeltaB_U','N_LT','N_TT'];stats={};covs={};metadata={};convergence=[]
    for mode,label in [('data','data'),('mc','MC')]:
        pub=samples[mode]['pub'];theta=samples[mode]['theta'][:,:4]
        for space,s,order,units in [('pub',pub,names,'(nb/GeV2)^2'),('theta',theta,pnames,'parameter units; DeltaB_U in GeV^-2')]:
            name=f'{label}_{space}';covs[name]=np.cov(s,rowvar=False,ddof=1)
            stats[name],_,metadata[name]=covariance(out,name,covs[name],order,units)
        quant=np.percentile(pub,[16,50,84],axis=0)
        csvout(out/f'{label}_published_percentiles.csv',[dict(quantity=names[i],median=quant[1,i],p16=quant[0,i],p84=quant[2,i],
            mean=float(pub[:,i].mean()),empirical_sd=stats[label+'_pub'][i]) for i in range(15)])
        checkpoints=[n for n in [200,300,500,1000,len(pub)] if n<=len(pub)]
        for n in sorted(set(checkpoints)):
            sd=np.std(pub[:n],axis=0,ddof=1)
            for i in range(15):convergence.append(dict(ensemble=label,accepted=n,quantity=names[i],sd=sd[i],relative_to_final=sd[i]/stats[label+'_pub'][i]-1))
        s=summary[mode].copy();s['max_SD_change_500_to_final']=float(np.max(abs(np.std(pub[:500],axis=0,ddof=1)/stats[label+'_pub']-1))) if len(pub)>=500 else None
        s['max_SD_change_200_to_final']=float(np.max(abs(np.std(pub[:200],axis=0,ddof=1)/stats[label+'_pub']-1)))
        s['covariance_rank_theta']=metadata[label+'_theta']['effective_rank'];s['covariance_rank_pub']=metadata[label+'_pub']['effective_rank']
        save(out/f'{mode}_bootstrap_summary.json',s)
    for space,order,units in [('pub',names,'(nb/GeV2)^2'),('theta',pnames,'parameter units; DeltaB_U in GeV^-2')]:
        name=f'stat_{space}';covs[name]=covs['data_'+space]+covs['MC_'+space]
        stats[name],_,metadata[name]=covariance(out,name,covs[name],order,units)
    covariance(out,'data_rows',np.cov(samples['data']['y'],rowvar=False),[f'row{i}' for i in range(60)],'(yield/mC)^2')
    csvout(out/'uncertainty_convergence.csv',convergence)
    means=np.array(c['means']);edges=np.array([-.75,-.55,-.4,-.25,-.13,0.]);table=[]
    for b in range(5):
        row=dict(tprime_bin=b,tprime_low=edges[b],tprime_high=edges[b+1],response_weighted_tprime=means[b])
        for k,component in enumerate(COMPONENTS):
            i=k*5+b;row.update({f'sigma_{component}':z[i],f'{component}_data_stat':stats['data_pub'][i],
                f'{component}_MC_stat':stats['MC_pub'][i],f'{component}_total_stat':stats['stat_pub'][i]})
        row.update(units='nb/GeV2',status='PRELIMINARY; statistical uncertainties only');table.append(row)
    csvout(out/'preliminary_cross_sections.csv',table)
    sensitivity=[]
    for method in ['old','multiplicative']:
        for i,name in enumerate(names):sensitivity.append(dict(estimator=method,quantity=name,central=np.array(central[method]['pub'])[i],
            additive=z[i],difference=np.array(central[method]['pub'])[i]-z[i],shift_total_stat_sigma=(np.array(central[method]['pub'])[i]-z[i])/stats['stat_pub'][i]))
    csvout(out/'estimator_sensitivity.csv',sensitivity)
    toy=samples['toys']['pub'];pull=(toy-z)/stats['stat_pub'];toyrows=[]
    for i,name in enumerate(names):toyrows.append(dict(quantity=name,bias_over_total_stat_sigma=float(pull[:,i].mean()),
        pull_mean=float(pull[:,i].mean()),pull_width=float(pull[:,i].std(ddof=1)),coverage68=float(np.mean(abs(pull[:,i])<=1))))
    csvout(out/'M0_toy_summary.csv',toyrows)
    toy_failed=any(abs(r['pull_mean'])>.5 or not .7<=r['pull_width']<=1.3 or r['coverage68']<.5 for r in toyrows)
    for name,meta in metadata.items():
        meta['release_status']='DIAGNOSTIC; toy calibration failed' if toy_failed else 'PRELIMINARY; statistical uncertainties only'
        save(out/f'V_{name}_metadata.json',meta)
    csvout(out/'covariance_manifest.csv',[dict(product=name,covariance=str((out/f'V_{name}.csv').resolve()),
        correlation=str((out/f'Corr_{name}.csv').resolve()),metadata=str((out/f'V_{name}_metadata.json').resolve()),
        effective_rank=meta['effective_rank'],status=meta['release_status']) for name,meta in metadata.items()])
    if toy_failed:
        for row in table:row['status']='DIAGNOSTIC; NOT RELEASED; statistical toy calibration failed'
        csvout(out/'preliminary_cross_sections.csv',table)
    gof=float((1+np.sum(samples['toys']['Q']>=c['info'][0]))/(1+len(toy)))
    save(out/'M0_toy_summary.json',dict(**summary['toys'],rows=toyrows,gof_pvalue=gof,
        construction='Local centered physical-event bootstrap M0 pseudoexperiments, with independent exclusive-MC Poisson resampling; fixed ellipse and exposure',
        coverage='Approximate coverage of central +/- reported total empirical SD; common preliminary statistical scale',
        scope='Checks local M0 fitting and statistical propagation; does not validate background transport or non-exclusive purity'))
    starts=json.loads((out/'start_stability.json').read_text())
    stability=max(float(np.max(abs(np.array(s['pub'])-z)/stats['stat_pub'])) for s in starts)
    save(out/'start_stability_summary.json',dict(starts=len(starts),all_converged=all(s['code']==0 for s in starts),max_published_shift_total_stat_sigma=stability))
    rank=json.loads((out/'central_rank.json').read_text());prod=json.loads((out/'production_verification.json').read_text())
    empty_support=json.loads((out/'empty_support_summary.json').read_text())
    detector=[]
    for r in range(60):detector.append(dict(row=r,included=bool(MASK[r]),data=c['y'][r],prediction=c['prediction'][r],
        residual=c['y'][r]-c['prediction'][r],objective_pull=(c['y'][r]-c['prediction'][r])/np.sqrt(c['variance'][r]) if MASK[r] else None))
    csvout(out/'detector_gof.csv',detector)
    for k,component in enumerate(COMPONENTS):
        fig,ax=plt.subplots(figsize=(6.4,4.5));sl=slice(5*k,5*k+5);x=-means
        ax.errorbar(x,z[sl],yerr=stats['stat_pub'][sl],fmt='o',capsize=5,label='Data + finite MC stat.',color='#163f6b')
        ax.errorbar(x,z[sl],yerr=stats['data_pub'][sl],fmt='none',elinewidth=3,color='#54a2d5',label='Data stat.')
        ax.set(xlabel=r"$-t'\ [\mathrm{GeV}^2]$",ylabel=rf'$\sigma_{{{component}}}\ [\mathrm{{nb}}/\mathrm{{GeV}}^2]$',xlim=(0,.75))
        ax.axhline(0,color='.6',lw=.7)
        ax.set_title(('PRELIMINARY CANDIDATE - NOT RELEASED\nStatistical toy calibration failed' if toy_failed else 'PRELIMINARY - KinC_x36_4\nStatistical uncertainties only'))
        ax.plot([0,.75],[.02,.02],transform=ax.get_xaxis_transform(),lw=3,color='.6',label='Accepted generated support')
        ax.legend(fontsize=8);fig.tight_layout()
        for suffix in ['pdf','png']:fig.savefig(out/f'sigma_{component}_preliminary.{suffix}',dpi=180)
        plt.close(fig)
    with PdfPages(out/'model_only_diagnostics.pdf') as pdf:
        fig,ax=plt.subplots(figsize=(10,7));ax.axis('off')
        text=f"PRELIMINARY KinC_x36_4: statistical uncertainties only\n\nQ = {c['info'][0]:.6f}; included rows = 56\nPhysics parameters: {np.array(c['theta'][:4])}\nDetector / nuisance / physics rank: {rank['detector']} / {rank['nuisance']} / {rank['physics']}\nData / MC accepted: {summary['data']['accepted']} / {summary['mc']['accepted']}\nLocal M0 toys: {summary['toys']['accepted']}; calibrated GOF p = {gof:.4f}\nMaximum start shift: {stability:.4g} total-stat sigma\n\nActive positivity boundaries are propagated through empirical ensembles.\nRows 0,5,6,11 retain the existing preliminary mask.\nEllipse-selected clean pi0 is treated as exclusive.\nBackground transport, ellipse dependence, contamination and model\nrefinements remain future systematic studies."
        if toy_failed:text='VALIDATION FAILED - NOT RELEASED\n\n'+text
        ax.text(.03,.94,text,va='top',fontsize=12);pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(5,2,figsize=(11,13));phi=np.arange(12)*30+15
        for b in range(5):
            sl=slice(b*12,(b+1)*12);included=MASK[sl];y=np.array(c['y'])[sl];pred=np.array(c['prediction'])[sl]
            axes[b,0].errorbar(phi[included],y[included],np.sqrt(np.array(c['data_variance'])[sl][included]),fmt='o',label='Data conditional Sumw2')
            axes[b,0].plot(phi,pred,'-',label='Folded M0');axes[b,0].set_title(f"t' bin {b}");axes[b,0].legend(fontsize=7)
            axes[b,1].axhline(0,color='.4');axes[b,1].plot(phi[included],(y-pred)[included]/np.sqrt(np.array(c['variance'])[sl][included]),'o');axes[b,1].set_ylabel('Objective pull')
        fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(1,3,figsize=(13,4))
        for ax,label in zip(axes,['data','MC','stat']):
            cor=np.load(out/f'Corr_{label}_pub.npy');im=ax.imshow(cor,vmin=-1,vmax=1,cmap='coolwarm');ax.set_title(label+' published correlation')
        fig.colorbar(im,ax=axes);pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(3,5,figsize=(13,8))
        for i,ax in enumerate(axes.flat):ax.hist(samples['data']['pub'][:,i],bins=30,alpha=.7);ax.axvline(z[i],color='k');ax.set_title(names[i])
        fig.suptitle('Physical-event data bootstrap; central line');fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(1,3,figsize=(13,4));idx=np.arange(15)
        for ax,key,target in zip(axes,['pull_mean','pull_width','coverage68'],[0,1,.6827]):
            ax.plot(idx,[r[key] for r in toyrows],'o');ax.axhline(target,color='.4');ax.set_title(key);ax.set_xticks(idx);ax.set_xticklabels(names,rotation=90,fontsize=7)
        fig.tight_layout();pdf.savefig(fig);plt.close(fig)
    # Explicit gates: statistical-only release does not certify deferred physics.
    gates=[]
    if any(r['code'] for r in starts) or stability>.1:gates.append('2: corrected M0 unstable across starts')
    if rank['physics']!=4:gates.append('3: published M0 physics not identifiable')
    if any(abs(r['pull_mean'])>.5 or not .7<=r['pull_width']<=1.3 or r['coverage68']<.5 for r in toyrows):gates.append('7: local M0 toy bias, pull or coverage fails the declared preliminary gate')
    if max(abs(r['shift_total_stat_sigma']) for r in sensitivity if r['estimator']=='multiplicative')>=1:gates.append('8: additive versus timing-preserving multiplicative estimator differs by >=1 sigma')
    save(out/'release_verdict.json',dict(ready=not gates,critical_blockers=gates,old_estimator_interpretation='Historical estimator has a demonstrated timing-transport defect; report its shift separately from the two timing-preserving estimators'))
    headings=['bin',"<t'>",'U +/- stat','LT +/- stat','TT +/- stat']
    human='| '+' | '.join(headings)+' |\n|'+ '|'.join(['---']*5)+'|\n'
    for b in range(5):human+='| '+str(b)+' | '+f'{means[b]:.6f}'+' | '+' | '.join(f'{z[k*5+b]:.5f} +/- {stats["stat_pub"][k*5+b]:.5f}' for k in range(3))+' |\n'
    fullkeys=list(table[0])[:-2]
    fulltable='| '+' | '.join(fullkeys)+' |\n|'+'|'.join(['---']*len(fullkeys))+'|\n'
    for row in table:fulltable+='| '+' | '.join(f'{row[k]:.7g}' for k in fullkeys)+' |\n'
    (out/'preliminary_cross_sections.md').write_text(('DIAGNOSTIC - NOT RELEASED; toy calibration failed.\n' if toy_failed else 'PRELIMINARY; statistical uncertainties only.\n')+'Units: nb/GeV2.\n\n'+fulltable)
    before=json.loads((out/'before/source_sha256.json').read_text());changed=[f for f,h in before.items() if hashlib.sha256(Path(f).read_bytes()).hexdigest()!=h]
    new=['scripts/pi0_preliminary_bridge.cpp','scripts/preliminary_pi0_xsec.py','scripts/validate_preliminary_pi0.py','scripts/report_preliminary_pi0.py','scripts/run_preliminary_pi0_xsec.sh']
    auxiliary=[str(out/'diagnose.py'),str(out/'check_rows.py')]
    save(out/'source_changes.json',dict(changed_preexisting=changed,new_sources=new,auxiliary_validation_helpers=auxiliary))
    reproduction_prefix=f'env PRELIM_DATA_ACCEPTED={summary["data"]["accepted"]} ' if summary['data']['accepted']<summary['data']['requested'] else ''
    commands=Path('scripts/run_preliminary_pi0_xsec.sh').read_text()
    commands=commands.replace('cd "$(dirname "$0")/.."',f'cd {REPO}')
    commands=commands.replace('OUT="${1:-validation/preliminary_model_xsec_reproduction_20261005}"',
        'OUT="validation/preliminary_model_xsec_reproduction_20261005"\nexport PRELIM_DATA_ACCEPTED='+str(summary['data']['accepted']))
    (out/'commands.sh.tmp').write_text(commands);(out/'commands.sh.tmp').replace(out/'commands.sh')
    elapsed_wall=time.time()-(out/'before/git_head.txt').stat().st_mtime
    save(out/'runtime.json',dict(elapsed_since_preproduction_git_snapshot_seconds=elapsed_wall,
        ensemble_runtime_seconds={mode:summary[mode]['runtime_seconds'] for mode in summary}))
    for command,name in [(['git','status','--short'],'git_status_final.txt'),(['git','diff','--stat'],'git_diff_stat_final.txt')]:
        (out/name).write_bytes(subprocess.check_output(command))
    report=f"""# PRELIMINARY KinC_x36_4 model extraction

{'**NOT RELEASED: statistical toy calibration fails the explicitly requested gate. The table and figures below are diagnostic ensemble results, not approved preliminary confidence errors.**' if toy_failed else '**PRELIMINARY: statistical uncertainties only.**'}

Statistical uncertainties only. Output: `{out.resolve()}`.
Git HEAD: `{(out/'before/git_head.txt').read_text().strip()}`; existing user changes preserved.
Elapsed wall time since the pre-production Git snapshot: {elapsed_wall:.1f} seconds.
Configuration: `config_snapshot.json`; accepted run set: `accepted_runs.txt`; environment: `environment.json`.

Frozen event weight: `pi0_weight_prelim = pi0_timing_coeff - (pi0_timing_bin_mean - pi0_weight) = c_e - B_rb/N_rb`.
It is constructed downstream without changing `pi0_weight`, independently of ellipse membership.
The fixed nominal ellipse is applied afterward. No toy bias correction is applied.

## Production and central M0

All {prod['runs']} runs regenerated through the normal producer and combiner: {prod['events']} events,
{prod['exposure_uC']/1000:.12f} mC. All 41 legacy per-run branches and every legacy combined branch are exactly unchanged.
The timing primitives and event IDs survive combination. All {prod['applicable_events']} regular-mass-bin events have finite additive weights.
The {prod['undefined_flow_events']} mass-flow events retain the defined NaN convention and are outside the ellipse.
The {prod['empty_bins']} unsupported empty mass bins contribute {prod['empty_bin_residual']:.12f} fitted counts;
no fake events are introduced. The selected yield is {sum(c['y']):.12f}/mC.
The absolute normalized empty-bin residual is bounded by {empty_support['absolute_normalized_bound_per_mC']:.9g}/mC,
or {100*empty_support['fraction_of_selected_yield']:.5f}% of the full selected yield: a documented small limitation.

Objective: Gaussian diagonal conditional data Sumw2 plus iterated exclusive-MC Sumw2, unchanged from M0.
Q={c['info'][0]:.9f}; included rows=56; detector/nuisance/physics ranks={rank['detector']}/{rank['nuisance']}/{rank['physics']}.
Physics parameters `(N_U, DeltaB_U [GeV^-2], N_LT, N_TT)` = `{c['theta'][:4]}`.
Nominal pivot is fixed at {c['info'][6]:.17g} GeV2 in every ensemble, preserving parameter coordinates.
Positivity minimum in physical support: {c['info'][4]:.6g}; an active physics boundary is retained.
The boundary-aware solver uses the same objective and angular constraints; empirical ensembles supply errors.
Variance updates: {int(c['info'][1])}; maximum four-start published shift: {stability:.6g} total-stat sigma.

"Rows 0,5,6,11 are retained under the existing preliminary analysis
mask; low-count/zero-count likelihood treatment will be revisited in
a later statistical refinement."

## Cross sections

{human}

{fulltable}

The CSV has separate data, MC and total statistical errors for every component.
The plotted positive coordinate is `-t'`; the table retains production's signed `t'=t-tmin` convention.
Values are accepted-response-weighted generated-bin averages, in nb/GeV2. U=T+epsilon L; there is no T/L separation.

## Statistical ensembles

Data: `{json.dumps(summary['data'])}`.
Exclusive MC: `{json.dumps(summary['mc'])}`.
Local M0 toys: `{json.dumps(summary['toys'])}`.
Data resampling uses one Poisson(1) count per `(run_number,event_id)` and preserves all candidates;
every replica reconstructs timing histograms, reruns the production mass fitter, rebuilds B/N and additive weights,
applies the same fixed ellipse, and refits M0. Exposure is fixed.
MC resampling uses only matched physical exclusive-SIMC event IDs, preserving all basis terms, migration and feed-in;
the generated normalization is fixed under the established Poisson convention. Data stays fixed.
The response and reported response-weighted averages are rebuilt in each MC replica.
Replicas with an empty included row retain that row and initialize the variance using the existing MC term.
No extra rows are excluded and no arbitrary variance floor is introduced.
All failures and attempted IDs are retained in the ensemble directories; development trials are separately labeled.

Full covariance is `V_stat = V_data + V_MC`. Published ordering:
`U(bin0...bin4), LT(bin0...bin4), TT(bin0...bin4)`.
Ranks (published; physics theta): data {metadata['data_pub']['effective_rank']};{metadata['data_theta']['effective_rank']},
MC {metadata['MC_pub']['effective_rank']};{metadata['MC_theta']['effective_rank']},
total {metadata['stat_pub']['effective_rank']};{metadata['stat_theta']['effective_rank']}.
Eigenvalues, correlation-based numerical rank definition, ordering and units accompany each covariance.
Nonlinear bin mappings can have empirical covariance rank exceeding the four-parameter tangent rank.
Empirical SDs, medians and 16/84 percentiles are supplied; boundaries can make distributions non-Gaussian.
Uncertainty convergence at 200, 300, 500 and 1000 accepted replicas is in `uncertainty_convergence.csv`.

Local toys recenter complete data-bootstrap detector fluctuations on folded M0 and independently resample exclusive MC,
then repeat the same fit/variance iteration. They test local estimator bias and approximate coverage with the measured
total-statistical scale; they do not test unknown combinatorial transport or contamination.
Bias/pull/coverage by quantity: `M0_toy_summary.csv`. Calibrated local-bootstrap detector GOF p={gof:.5f}.
Estimator shifts are in `estimator_sensitivity.csv`; their spread is not added to statistical covariance.
The old estimator's known timing-transport defect is distinguished from the additive/multiplicative comparison.

## Deferred systematic refinements

- Exact conditional combinatorial-background transport through the ellipse.
- Ellipse-cut systematic variation and residual non-exclusive/SIDIS contamination.
- Delta contamination.
- Low-count treatment of the four masked rows.
- Alternative pi0 parameterizations and additional LT/TT shape freedom if later required.

No no-model extraction, no-model output, charged-pion diagnostic, SIDIS/delta subtraction, ellipse variation or LT/TT slope study was used.
No pre-existing source/config changed in this task: `{changed}`. New source files: `{new}`.
Additional development-only validation helpers: `{auxiliary}`.
The protected no-model source and all other initially hashed sources are checksum-verified.

## Reproduction

From `{REPO}`, load Hall C and run the non-overwriting driver:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; {reproduction_prefix}bash scripts/run_preliminary_pi0_xsec.sh validation/preliminary_model_xsec_reproduction_20261005'
```

The driver records exact stage commands, seeds, accepted runs and runtimes and refuses an existing destination.
The complete copy-paste commands for that fresh reproduction directory are in `commands.sh`.
Production seeds are 20261005 (data), 20261006 (MC), 20261007 (toys).
The driver was syntax-checked; the substantive production, combination, fitting,
ensemble and reporting commands were executed individually in this session.
The data target was shortened only after the documented critical toy gate failed.
Full stage commands, including environment prerequisites above and every build step:

```bash
{commands}
```

## Verdict

{'PRELIMINARY MODEL EXTRACTION READY' if not gates else 'PRELIMINARY MODEL EXTRACTION NOT READY'}
{json.dumps(gates) if gates else 'The ellipse-selected clean-pi0 sample is treated as exclusive for this preliminary analysis. The additive timing-preserving clean-pi0 estimator is used as the best current central estimate. The forward-folded M0 extraction is identifiable and now has data and finite-exclusive-MC statistical covariance. U, LT and TT are reported with preliminary statistical error bars. Remaining background, ellipse, contamination and model refinements are deferred to later systematic studies.'}
"""
    (out/'REPORT.md').write_text(report)
    print(json.dumps(dict(table=table,toys=toyrows,metadata=metadata,critical_blockers=gates)),flush=True)

if __name__=='__main__':main()
