#!/usr/bin/env python3
"""Plot timing-transport closure evidence; no cross-section publication."""
import csv
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    root=Path(sys.argv[1])
    rows=[]
    for campaign in ('toys_10percent','toys'):
        for r in csv.DictReader((root/campaign/'toy_summary.csv').open()):
            r=dict(campaign=campaign,**r)
            sd=float(r['empirical_estimate_spread'])
            r['bias_over_estimate_spread']=float(r['bias'])/sd if sd else 0.
            rows.append(r)
    with (root/'decision_toy_summary.csv.tmp').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (root/'decision_toy_summary.csv.tmp').replace(root/'decision_toy_summary.csv')
    cases=('clean','clean_timing','clean_comb','clean_timing_comb')
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for ax,campaign,title in zip(axes,('toys_10percent','toys'),('Continuum / signal = 10%','Continuum / signal = 300% (stress)')):
        for offset,method in zip((-.2,0,.2),('old','additive','multiplicative')):
            r=[next(r for r in rows if r['campaign']==campaign and r['case']==case and r['method']==method) for case in cases]
            y=[r['bias_over_estimate_spread'] for r in r]
            se=[float(r['empirical_residual_spread'])/np.sqrt(int(r['successful_toys']))/float(r['empirical_estimate_spread']) for r in r]
            ax.errorbar(np.arange(4)+offset,y,yerr=se,fmt='o',label=method,capsize=3)
        ax.axhline(0,color='black',lw=.8);ax.set_xticks(range(4),['clean','+ timing','+ comb.','+ both'])
        ax.set_title(title);ax.set_ylabel('Mean (estimate - injected selected truth) / estimate SD')
        ax.legend(fontsize=9);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Fixed nominal ellipse; 200 successful production mass fits per case\nBars: paired mean-error uncertainty; declared toy rates, not measured data contamination')
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(root/f'timing_transport_closure.{ext}',dpi=160)
    plt.close(fig)
    # Exact populated-bin counterexample: c=1 for 100 signal + 100 continuum;
    # the downstream cut accepts 90 signal and 10 continuum events.
    counterexample=dict(N=200,T=200,B=100,S=100,selected_signal=90,selected_background=10,
                        old_selected=50,additive_selected=50,multiplicative_selected=50)
    result=dict(verdict='PRELIMINARY MODEL EXTRACTION NOT READY',production_estimator=None,
        blocker='Both allowed scalar-bin candidates fail combinatorial transport closure even with oracle background.',
        counterexample=counterexample,corrected_M0_Q=None,corrected_structure_functions=None,
        data_covariance=None,exclusive_MC_covariance=None,total_statistical_errors=None,
        rows_0_5_6_11='Not reassessed: no validated production estimator; no new row exclusions.',
        production_scope='Three runs regenerated and combined for branch validation; full campaign not released.')
    (root/'decision.json.tmp').write_text(json.dumps(result,indent=2)+'\n')
    (root/'decision.json.tmp').replace(root/'decision.json')


if __name__=='__main__':main()
