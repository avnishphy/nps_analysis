#!/usr/bin/env python3
"""Render measured objective-stage tables and standalone diagnostic figures."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from study_pi0_objectives import write


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
        ['| '+' | '.join(map(str,r))+' |' for r in rows])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',type=Path,required=True)
    p.add_argument('--validation',type=Path,required=True);a=p.parse_args();base=a.stage
    load=lambda directory:json.loads((base/directory/'summary.json').read_text())
    coverage=load('coverage_v4');boundary=load('boundary_v4')
    screen=load('six_screen_v4')+load('six_low_small_screen_v4')
    onoff=load('onoff_all_screen_v3');fixed=load('fixed_screen')
    ready=json.loads((base/'release_gate/readiness.json').read_text())
    names={0:'Primitive Poisson',1:'Model Gaussian + log V',2:'Observed-variance chi-square'}
    comparison=[]
    for kind,rows in [('fixed shape, on/off',fixed),('free shape, on/off (v3 optimizer)',onoff),('free shape, six categories (v4)',screen)]:
        for r in rows:
            comparison.append(dict(experiment=kind,objective=names[r['mode']],case=r['case'],truth_A=r['truth_A'],truth_Y=r['truth_Y'],
                mean_A=r['mean_A'],A_bias=r['A_bias'],yield_bias=r['yield_bias'],yield_bias_se=r['yield_bias_se'],
                absolute_bias_over_sigma=r['absolute_bias_over_sigma'],successful=r['successful'],failed=r['failed']))
    write(base/'objective_comparison.csv',comparison);write(base/'coverage_table.csv',coverage);write(base/'boundary_table.csv',boundary)
    diagnostic_free=list(csv.DictReader((base/'six_low_small_screen_v4/toys.csv').open()))
    diagnostic_fixed={int(r['experiment']):r for r in csv.DictReader((base/'six_low_small_fixed_v3/toys.csv').open())}
    diff=np.array([float(r['Y'])-float(diagnostic_fixed[int(r['experiment'])]['Y']) for r in diagnostic_free])
    width=np.array([float(r['width']) for r in diagnostic_free]);turn=np.array([float(r['turn']) for r in diagnostic_free])
    independent=load('optimizer_checks_v4')['independent_checks']
    profile=list(csv.DictReader((base/'optimizer_checks_v4/turn_profile.csv').open()))
    pq=np.array([float(r['objective']) for r in profile]);pb=np.array([float(r['background']) for r in profile]);near=pq-pq.min()<.001
    diagnosis=dict(paired_toys=len(diff),mean_free_minus_fixed_yield=float(diff.mean()),paired_difference_se=float(diff.std(ddof=1)/np.sqrt(len(diff))),
        minimum_width_fraction=float(np.mean(width<.00101)),maximum_width_fraction=float(np.mean(width>.09999)),
        minimum_turn_fraction=float(np.mean(turn<.11001)),maximum_turn_fraction=float(np.mean(turn>.21999)),
        independent_optimizer_max_improvement=max(r['objective_improvement'] for r in independent),
        flat_profile_background_range=[float(pb[near].min()),float(pb[near].max())])
    (base/'shape_diagnosis_v4.json').write_text(json.dumps(diagnosis,indent=2)+'\n')
    withheld=dict(reason='Corrected free-shape estimator fails controlled low-signal bias/pull validation.',
        corrected_nominal_yields=None,new_data_covariance=None,new_MC_covariance=None,new_total_covariance=None,
        new_correlations=None,final_cross_sections=None,real_data_rerun=False,
        previous_evidence_preserved='validation/pi0_uncertainty_20261004/',
        old_covariances_reused=False,empirical_bias_correction_applied=False)
    (base/'withheld_outputs.json').write_text(json.dumps(withheld,indent=2)+'\n')
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    truth=np.array([r['truth_A'] for r in boundary]);mean=np.array([r['mean_A'] for r in boundary])
    ax[0].errorbar(truth,mean,yerr=[r['A_mean_se'] for r in boundary],fmt='o-',capsize=3)
    ax[0].plot([0,.4],[0,.4],'k--',lw=1);ax[0].set(xlabel='True amplitude (events/bin)',ylabel='Mean fitted amplitude',title='Six timing categories; free Fermi shape')
    ax[1].errorbar(truth,[100*r['coverage'] for r in boundary],yerr=[100*r['coverage_se'] for r in boundary],fmt='o-',capsize=3)
    ax[1].axhline(68.2689492,color='k',ls='--',lw=1);ax[1].set(xlabel='True amplitude (events/bin)',ylabel='One-sigma coverage (%)',title='600 outer toys; 100 event-Poisson replicas',ylim=(60,76))
    fig.savefig(a.validation/'boundary_scan.pdf');fig.savefig(a.validation/'boundary_scan.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    ax[0].hist(width,bins=np.linspace(.001,.1,31));ax[0].set(xlabel='Fitted width (GeV)',ylabel='Pseudo-experiments',title='A true=0.5, signal=300; 4,000 toys')
    pt=np.array([float(r['turn']) for r in profile]);ax[1].plot(pt,pq-pq.min());ax[1].set(xlabel='Turn (GeV), width fixed at 0.001 GeV',ylabel='Profile deviance - minimum',ylim=(0,3),title='Representative weak-shape toy')
    fig.savefig(a.validation/'weak_shape_diagnosis.pdf');fig.savefig(a.validation/'weak_shape_diagnosis.png',dpi=160);plt.close(fig)
    headers=['Experiment','Objective','A true','Mean A','A bias','Yield bias +/- SE','|bias| / sigma','Accepted']
    rows=[[r['experiment'],r['objective'],f"{r['truth_A']:g}",f"{r['mean_A']:.6f}",f"{r['A_bias']:+.6f}",
        f"{r['yield_bias']:+.4f} +/- {r['yield_bias_se']:.4f}",f"{r['absolute_bias_over_sigma']:.4f}",r['successful']] for r in comparison]
    ctable=table(['Case','A true','Y true','Bias +/- SE','Bootstrap SD','Empirical SD','Pull mean','Pull width','Coverage +/- SE'],
        [[r['case'],r['truth_A'],r['truth_Y'],f"{r['yield_bias']:+.3f} +/- {r['yield_bias_se']:.3f}",f"{r['mean_bootstrap_sd']:.3f}",f"{r['empirical_sd']:.3f}",
          f"{r['pull_mean']:+.4f}",f"{r['pull_width']:.4f}",f"{100*r['coverage']:.2f}% +/- {100*r['coverage_se']:.2f}%"] for r in coverage])
    btable=table(['A true','Mean A +/- SE','Yield bias','Zero fraction','Coverage +/- SE'],
        [[r['truth_A'],f"{r['mean_A']:.5f} +/- {r['A_mean_se']:.5f}",f"{r['yield_bias']:+.3f}",f"{r['zero_fraction']:.4f}",f"{100*r['coverage']:.2f}% +/- {100*r['coverage_se']:.2f}%"] for r in boundary])
    gates=table(['Gate','Result'],list(ready['gates'].items()))
    full='# Objective study tables\n\nAll errors on means/coverage are Monte Carlo standard errors. No truth-fixed fit is a production estimator.\n\n'+table(headers,rows)+'\n\n## Final candidate coverage\n\n'+ctable+'\n\n## Boundary scan\n\n'+btable+'\n\n## Release gates\n\n'+gates+'\n'
    (a.validation/'TABLES.md').write_text(full)
    snippets=dict(coverage=ctable,boundary=btable,gates=gates)
    (a.validation/'report_tables.json').write_text(json.dumps(snippets,indent=2)+'\n')
    print(json.dumps(dict(diagnosis=diagnosis,verdict=ready['verdict']),indent=2))


if __name__=='__main__':main()
