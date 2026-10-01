#!/usr/bin/env python3
"""Independent Asimov and compound-Poisson toy check on a saved response.

Usage: python validate_scaled_poisson_toys.py FIT_DIR COMBINED_DATA_ROOT [N_TOYS]
This tests one published U/LT/TT block with other truth blocks held fixed.
It is a diagnostic, not a replacement for a full refit toy coverage study.
"""
import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import uproot
from matplotlib.path import Path as Polygon
with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)
    from scipy.optimize import minimize


def table(path):
    with open(path) as stream:
        return list(csv.DictReader(stream))


def load_problem(fit_dir):
    parameters=table(fit_dir/"migration_parameters.csv")
    p=np.array([float(row["value"]) for row in parameters])
    design=table(fit_dir/"migration_design.csv")
    x=np.zeros((max(int(r["reco_row"]) for r in design)+1,len(p)))
    for r in design:
        x[int(r["reco_row"]),int(r["parameter_index"])]=float(r["response"])
    diagnostics=table(fit_dir/"scaled_poisson_rows.csv")
    keep=np.array([r["included"]=="1" for r in diagnostics])
    return x,p,diagnostics,keep


def empirical_weights(data_root):
    names=["Q2","xB","t","tmin","phi","pi0_weight","scale",
           "charge_uC","run_number","is_exclusive_ellipse_combined"]
    a=uproot.open(data_root)["physics"].arrays(names,library="np")
    q,x,tp=a["Q2"],a["xB"],a["t"]-a["tmin"]
    diamond=Polygon([(.28536635126633153,3.243861779409012),
                     (.33822268816213624,3.6193121389727954),
                     (.43953011480920373,4.698466077676237),
                     (.3721482218990109,4.244414308310078)])
    mask=(a["is_exclusive_ellipse_combined"]>0)&(q>=3.3)&(q<=4.7)&(x>=.29)&(x<=.44)&(tp>=-.75)&(tp<=0)
    mask &= diamond.contains_points(np.column_stack((x,q)))
    runs=np.unique(a["run_number"])
    charges={r:float(a["charge_uC"][np.flatnonzero(a["run_number"]==r)[0]]) for r in runs}
    total_charge=sum(charges.values())
    w=a["pi0_weight"]*a["scale"]*a["charge_uC"]/total_charge/.584
    block=np.searchsorted([-.75,-.5,-.3,-.2,-.1,0],tp,side="right")-1
    pools=[]
    for b in range(5):
        v=w[mask&(block==b)&np.isfinite(w)&(w>0)]
        if not len(v): raise RuntimeError("empty empirical weight block")
        pools.append(v)
    return pools


def minimum_required(lt,tt,eps):
    z=np.array([-1.,1.])
    if tt>0:
        vertex=-np.sqrt(2*eps*(1+eps))*lt/(4*eps*tt)
        if -1<vertex<1: z=np.append(z,vertex)
    return max(0.,-np.min(np.sqrt(2*eps*(1+eps))*lt*z+eps*tt*(2*z*z-1)))


def fit_one(y,sumw2,scales,x,offset,eps,truth,objective):
    unit=1e-8
    def physical(z):
        lt=z[1]*unit;tt=z[2]*unit
        return np.array([minimum_required(lt,tt,eps)+z[0]**2*unit,lt,tt])
    def score(z):
        mu=offset+x@physical(z)
        if np.any(mu<=0): return 1e30
        if objective=="scaled":
            out=np.zeros_like(y)
            pos=y>0
            out[~pos]=2*mu[~pos]/scales[~pos]
            delta=(mu[pos]-y[pos])/y[pos]
            out[pos]=2*y[pos]*(delta-np.log1p(delta))/scales[pos]
            return float(np.sum(np.maximum(out,0)))
        pos=sumw2>0
        return float(np.sum((y[pos]-mu[pos])**2/sumw2[pos]))
    gap=max(1e-5,(truth[0]-minimum_required(truth[1],truth[2],eps))/unit)
    start=np.array([np.sqrt(gap)*1.05,truth[1]/unit*.95,truth[2]/unit*1.05])
    result=minimize(score,start,method="L-BFGS-B",
                    bounds=[(0,None),(None,None),(None,None)],
                    options={"ftol":1e-12,"gtol":1e-8,"maxiter":2000})
    return physical(result.x),float(result.fun),bool(result.success)


def main():
    fit_dir=Path(sys.argv[1])
    data_root=Path(sys.argv[2])
    n_toys=int(sys.argv[3]) if len(sys.argv)>3 else 80
    x,p,diagnostics,keep=load_problem(fit_dir)
    block=1
    eps=float(table(fit_dir/"migration_truth_blocks.csv")[block]["epsilon_max"])
    a=x[keep,3*block:3*block+3]
    nominal=x[keep]@p
    truth=p[3*block:3*block+3].copy()
    truth[0]*=1.35 # interior reference; avoid boundary-coverage effects
    offset=nominal-a@p[3*block:3*block+3]
    expected=offset+a@truth
    if np.any(expected<=0): raise RuntimeError("nonpositive Asimov mean")
    reference=np.array([float(r["s_used"]) for r in diagnostics if r["included"]=="1"])
    asimov,asimov_d,asimov_ok=fit_one(expected,reference*expected,reference,
                                       a,offset,eps,truth,"scaled")
    pools=empirical_weights(data_root)
    means=np.array([pool.mean() for pool in pools])
    rng=np.random.default_rng(20260929)
    estimates={"scaled":[],"gaussian":[]}
    zero_rows=[]
    for _ in range(n_toys):
        y=np.zeros(len(expected)); v=np.zeros(len(expected)); s=reference.copy()
        for i,row in enumerate(np.flatnonzero(keep)):
            weights=pools[row//12]
            count=rng.poisson(expected[i]/means[row//12])
            if count:
                sample=rng.choice(weights,size=count,replace=True)
                y[i]=sample.sum();v[i]=np.dot(sample,sample);s[i]=v[i]/y[i]
        zero_rows.append(int(np.count_nonzero(y==0)))
        for name in estimates:
            result,_,ok=fit_one(y,v,s,a,offset,eps,truth,name)
            if ok: estimates[name].append(result)
    report={"block":block,"toys":n_toys,"truth":truth.tolist(),
            "asimov_fit":asimov.tolist(),"asimov_max_abs_error":float(np.max(np.abs(asimov-truth))),
            "asimov_deviance":asimov_d,"asimov_success":asimov_ok,
            "mean_zero_rows":float(np.mean(zero_rows))}
    for name,points in estimates.items():
        points=np.asarray(points)
        report[name]={"success":len(points),
                      "bias":(points.mean(axis=0)-truth).tolist() if len(points) else None,
                      "spread":points.std(axis=0,ddof=1).tolist() if len(points)>1 else None}
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
