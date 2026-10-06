#!/usr/bin/env python3
"""Published-parameter sensitivity for all rejected data replicas; not production."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import Fitter,RunEvents,selected_estimate,write_csv


def read(p):return list(csv.DictReader(p.open()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('data','response','manifest','campaign','library','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    accepted=[r for r in read(a.manifest) if r['accepted']=='1'];runs=[RunEvents(a.campaign/'root'/f"signal_events_run{r['run']}.root") for r in accepted]
    q=np.array([np.float32(float(r['charge_uC'])) for r in accepted],float)
    scale=np.array([np.float32(float(r['prescale'])/(float(r['charge_uC'])/1000*float(r['livetime'])*float(r['efficiency']))) for r in accepted],float)
    k=q*scale/1000;norm=1000/q.sum()/.584
    params=read(a.response/'migration_parameters.csv');pub=np.array([r['region']=='published' for r in params]);A=np.zeros((60,len(params)))
    for r in read(a.response/'migration_design.csv'):A[int(r['reco_row']),int(r['parameter_index'])]=float(r['response'])
    fixed=np.array([float(r['fixed_prediction']) for r in read(a.response/'migration_reco_rows.csv')])
    fit=Fitter(a.library);rows=[];extra=[]
    for b in sorted({int(r['replica']) for r in read(a.data/'replica_failures.csv')}):
        m=[np.random.default_rng(np.random.SeedSequence([20261004,b,r.run])).poisson(1,len(r.coefficient)).astype(float) for r in runs]
        y,v,bad=selected_estimate(runs,m,k,fit)
        if bad:raise RuntimeError('Selected-only sensitivity cannot evaluate replica '+str(b))
        y*=norm;v*=norm**2;keep=v>0;aa=A[keep]/np.sqrt(v[keep,None]);scale=np.linalg.norm(aa,axis=0);scale=np.where(scale>0,scale,1)
        u,s,vt=np.linalg.svd(aa/scale,full_matrices=False);rank=int(np.sum(s>s[0]*1e-10));null=float(np.linalg.norm(vt[rank:,pub]))
        if null>1e-8:raise RuntimeError('Published coefficients unidentified in failed data replica')
        par=((vt[:rank].T/s[:rank])@(u[:,:rank].T@((y[keep]-fixed[keep])/np.sqrt(v[keep]))))/scale
        extra.append(par[pub]);rows.append(dict(replica=b,rank=rank,published_null_projection=null,zero_variance_rows=';'.join(map(str,np.flatnonzero(~keep)))))
    z=np.load(a.data/'replicas.npz')['sigma'][:,pub];complete=np.vstack([z,extra]);ratio=complete.std(0,ddof=1)/z.std(0,ddof=1)
    write_csv(a.output/'rejected_replica_diagnostics.csv',rows)
    result=dict(rejected_replicas=len(extra),max_published_relative_sigma_change=float(np.max(abs(ratio-1))),
        nominal_ensemble_unchanged=True,method='Selected estimator only, retaining unique published coefficients when nuisance directions are unidentified; sensitivity, not replica acceptance')
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
