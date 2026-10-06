#!/usr/bin/env python3
"""Controlled nested event-Poisson coverage tests; does not change production.

Poisson resampling of disjoint observed event counts is exactly equivalent to
assigning independent Poisson(1) multiplicities to those physical events.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import Fitter, write_csv


def estimate(fit, coin, side):
    background=side/6
    # One-run version of the production timing normalization-error term.
    variance=coin+side/36+background**2/max(1.,side.sum())
    code,final,info=fit.fit(coin-background,variance)
    return code,float(final.sum()),info


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--experiments',type=int,default=200)
    p.add_argument('--bootstrap',type=int,default=150)
    p.add_argument('--seed',type=int,default=20261005)
    p.add_argument('--signal-yield',type=float,default=300.)
    p.add_argument('--accidental-per-bin',type=float,default=2.)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    fit=Fitter(a.library);x=(np.arange(200)+.5)*.002
    signal=np.exp(-.5*((x-.135)/.006)**2);signal*=a.signal_yield/signal.sum()
    rows=[];failures=[];summaries=[]
    for case,amplitude in enumerate((0.,.5,12.)):
        continuum=amplitude/(1+np.exp((x-.145)/.025))
        accepted=[];ninner=0;failed_inner=0
        for exp in range(a.experiments):
            rng=np.random.default_rng(np.random.SeedSequence([a.seed,case,exp]))
            c=rng.poisson(signal+continuum+a.accidental_per_bin).astype(float)
            s=rng.poisson(np.full(200,6*a.accidental_per_bin)).astype(float)
            code,y,info=estimate(fit,c,s)
            if code:
                failures.append(dict(case=case,experiment=exp,replica=-1,status=int(info[1]),covariance=int(info[2])))
                continue
            samples=[]
            for b in range(a.bootstrap):
                code,yb,bi=estimate(fit,rng.poisson(c).astype(float),rng.poisson(s).astype(float))
                ninner+=1
                if code:
                    failed_inner+=1
                    failures.append(dict(case=case,experiment=exp,replica=b,status=int(bi[1]),covariance=int(bi[2])))
                else:samples.append(yb)
            if len(samples)<.9*a.bootstrap:continue
            sd=float(np.std(samples,ddof=1));pull=(y-signal.sum())/sd
            row=dict(case=case,amplitude=amplitude,experiment=exp,true_yield=float(signal.sum()),
                     estimate=y,bootstrap_sd=sd,pull=pull,covered=int(abs(pull)<=1),
                     zero_background=int(info[0]),successful_replicas=len(samples))
            rows.append(row);accepted.append(row)
            if exp%20==0:print(json.dumps(dict(case=case,experiment=exp,accepted=len(accepted))),flush=True)
        ys=np.array([r['estimate'] for r in accepted]);sd=np.array([r['bootstrap_sd'] for r in accepted])
        pulls=np.array([r['pull'] for r in accepted]);coverage=float(np.mean(abs(pulls)<=1))
        summaries.append(dict(case=case,amplitude=amplitude,true_background=float(continuum.sum()),
            requested=a.experiments,successful=len(accepted),nominal_failed=sum(r['case']==case and r['replica']==-1 for r in failures),
            inner_failed=failed_inner,inner_requested=ninner,bias=float(ys.mean()-signal.sum()),
            bias_mc_se=float(ys.std(ddof=1)/np.sqrt(len(ys))),empirical_sd=float(ys.std(ddof=1)),
            mean_bootstrap_sd=float(sd.mean()),pull_mean=float(pulls.mean()),pull_width=float(pulls.std(ddof=1)),
            coverage=coverage,coverage_mc_se=float(np.sqrt(coverage*(1-coverage)/len(ys))),
            zero_fraction=float(np.mean([r['zero_background'] for r in accepted]))))
        write_csv(a.output/'pseudoexperiments.csv',rows);write_csv(a.output/'failures.csv',failures)
        (a.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
        print(json.dumps(summaries[-1]),flush=True)
    (a.output/'configuration.json').write_text(json.dumps(vars(a),default=str,indent=2)+'\n')


if __name__=='__main__':main()
