#!/usr/bin/env python3
"""Convergence, same-stream reproducibility and failed-replica sensitivity."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import Fitter,RunEvents,selected_estimate,write_csv


def read(p):return list(csv.DictReader(p.open()))
def corr(x):
    c=np.cov(x,rowvar=False);s=np.sqrt(c.diagonal())
    return np.divide(c,np.outer(s,s),out=np.zeros_like(c),where=np.outer(s,s)>0)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('data','previous','response','campaign','manifest','library','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    z=np.load(a.data/'replicas.npz');old=np.load(a.previous/'replicas.npz')
    for k in old.files:
        if not np.array_equal(old[k],z[k][z['replica']<500]):raise RuntimeError('Prior deterministic prefix differs: '+k)
    params=read(a.response/'migration_parameters.csv');full=z['sigma'].std(0,ddof=1);fc=corr(z['sigma'])
    published=np.array([r['region']=='published' for r in params]);rows=[];summary=[];pairs=[]
    important=sorted([(i,j) for i in range(len(params)) for j in range(i) if published[i] or published[j]],key=lambda ij:abs(fc[ij]),reverse=True)[:10]
    for B in (100,250,500,1000,1500,2000):
        keep=z['replica']<B;sd=z['sigma'][keep].std(0,ddof=1);cc=corr(z['sigma'][keep]);rc=corr(z['reco'][keep])
        for i,r in enumerate(params):
            rows.append(dict(requested=B,successful=int(keep.sum()),parameter_index=i,truth_block=r['truth_block'],component=r['component'],region=r['region'],sigma=sd[i],relative_change=(sd[i]-full[i])/full[i]))
        for i,j in important:pairs.append(dict(requested=B,i=i,j=j,correlation=cc[i,j],reference=fc[i,j],change=cc[i,j]-fc[i,j]))
        summary.append(dict(requested=B,successful=int(keep.sum()),max_published_relative_error_change=float(np.max(abs(sd[published]/full[published]-1))),
            median_published_relative_error_change=float(np.median(abs(sd[published]/full[published]-1))),
            largest_parameter_abs_correlation=float(np.max(abs(cc-np.eye(len(params))))),
            largest_reco_abs_correlation=float(np.max(abs(rc-np.diag(np.diag(rc)))))))
    write_csv(a.output/'parameter_convergence.csv',rows);write_csv(a.output/'convergence.csv',summary);write_csv(a.output/'important_correlations.csv',pairs)
    table=[]
    for i,r in enumerate(params):
        if not published[i]:continue
        ss={B:z['sigma'][z['replica']<B,i].std(ddof=1) for B in (500,1000,2000)}
        table.append(dict(parameter_index=i,truth_block=r['truth_block'],component=r['component'],sigma500=ss[500],sigma1000=ss[1000],sigma2000=ss[2000],change_500_to_2000=ss[2000]/ss[500]-1))
    write_csv(a.output/'published_stability.csv',table)
    accepted=[r for r in read(a.manifest) if r['accepted']=='1'];runs=[RunEvents(a.campaign/'root'/f"signal_events_run{r['run']}.root") for r in accepted]
    charges=np.array([np.float32(float(r['charge_uC'])) for r in accepted],float)
    scales=np.array([np.float32(float(r['prescale'])/(float(r['charge_uC'])/1000*float(r['livetime'])*float(r['efficiency']))) for r in accepted],float)
    normalization=1000/charges.sum()/.584;corrections=scales*charges/1000
    design=np.zeros((z['reco'].shape[1],len(params)))
    for r in read(a.response/'migration_design.csv'):design[int(r['reco_row']),int(r['parameter_index'])]=float(r['response'])
    fixed=np.array([float(r['fixed_prediction']) for r in read(a.response/'migration_reco_rows.csv')])
    fitter=Fitter(a.library);fail=read(a.data/'replica_failures.csv');extra=[];sensitivity=[]
    for b in sorted({int(r['replica']) for r in fail}):
        mult=[np.random.default_rng(np.random.SeedSequence([20261004,b,r.run])).poisson(1,len(r.coefficient)).astype(float) for r in runs]
        y,v,bad=selected_estimate(runs,mult,corrections,fitter)
        code,par,_=fitter.solve(design,y*normalization,v*normalization**2,fixed) if not bad else (1,None,None)
        if not code:extra.append(par)
        sensitivity.append(dict(replica=b,selected_only_valid=int(not code),selected_yield=float(y.sum()*normalization),
            max_parameter_z=float(np.max(abs((par-z['sigma'].mean(0))/full))) if not code else None))
    if extra:
        merged=np.vstack([z['sigma'],extra]);effect=float(np.max(abs(merged.std(0,ddof=1)/full-1)))
    else:effect=None
    for r in fail:r['fit_classification']='unresolved_invalid_fit'
    write_csv(a.output/'failure_classification.csv',fail);write_csv(a.output/'failure_sensitivity.csv',sensitivity)
    result=dict(requested=2000,successful=len(z['replica']),failed=2000-len(z['replica']),failure_fraction=(2000-len(z['replica']))/2000,
        previous_500_bitwise_identical=True,accepted_runs=len(accepted),charge_uC=float(charges.sum()),
        selected_only_sensitivity_successful=len(extra),failure_sensitivity_max_relative_sigma=effect,convergence=summary)
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
