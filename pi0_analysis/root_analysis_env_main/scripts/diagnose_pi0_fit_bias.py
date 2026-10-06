#!/usr/bin/env python3
"""Isolate observed-variance fit bias with the TRUE Fermi shape held fixed.

Oracle variances are a diagnostic control, not available production inputs.
No changes to the production estimator are made.
"""
import argparse,json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import write_csv


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--experiments',type=int,default=10000);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    x=(np.arange(200)+.5)*.002;f=1/(1+np.exp((x-.145)/.025));side=((x>=.01)&(x<=.11))|((x>=.15)&(x<=.4))
    signal=np.exp(-.5*((x-.135)/.006)**2)
    signal[side]=0. # Remove even tiny signal leakage to isolate variance bias.
    signal*=300/signal.sum();out=[]
    for case,A in enumerate((0.,.5,12.)):
        rng=np.random.default_rng(379+case);cmu=signal+A*f+2;bmu=np.full(200,12.)
        c=rng.poisson(cmu,(a.experiments,200));b=rng.poisson(bmu,(a.experiments,200));y=c-b/6
        observed=c+b/36+(b/6)**2/np.maximum(1,b.sum(1))[:,None]
        expected=cmu+bmu/36+(bmu/6)**2/bmu.sum()
        for label,var in [('observed',observed),('oracle_expected',np.broadcast_to(expected,observed.shape))]:
            v=np.where(var[:,side]>0,var[:,side],1.);ff=f[side]
            q=np.sum(y[:,side]*ff/v,axis=1);d=np.sum(ff**2/v,axis=1)
            for constrained in (False,True):
                amp=np.maximum(0,q/d) if constrained else q/d
                est=y.sum(1)-amp*f.sum()
                out.append(dict(amplitude=A,variance=label,nonnegative=int(constrained),experiments=a.experiments,
                    amplitude_mean=float(amp.mean()),yield_bias=float(est.mean()-signal.sum()),
                    bias_mc_se=float(est.std(ddof=1)/np.sqrt(a.experiments)),empirical_sd=float(est.std(ddof=1)),
                    boundary_fraction=float(np.mean(amp==0))))
    write_csv(a.output/'fixed_true_shape.csv',out)
    (a.output/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))


if __name__=='__main__':main()
