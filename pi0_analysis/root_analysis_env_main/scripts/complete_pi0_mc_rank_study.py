#!/usr/bin/env python3
"""Diagnose rank loss and recover ONLY uniquely estimable published MC terms.

The scaled Moore-Penrose solution fixes an arbitrary nuisance coordinate gauge.
It is accepted for physics only when every null vector has zero projection on
published coordinates. Nuisance-coordinate covariance is explicitly not an
identified physical uncertainty. No production central value is changed.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import write_csv,covariance_outputs


def read(p):return list(csv.DictReader(p.open()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('mc','response','config','first-order','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    config=json.loads(a.config.read_text());meta=json.loads((a.mc/'summary.json').read_text())
    z=np.load(a.mc/'event_basis.npz');old=np.load(a.mc/'replicas.npz')
    params=read(a.response/'migration_parameters.csv');reco=read(a.response/'migration_reco_rows.csv')
    npar,nrow=len(params),len(reco);pub=np.array([r['region']=='published' for r in params])
    blocks={int(r['truth_block']):int(r['active_block_index']) for r in params}
    ids=z['unique_event_id'];inverse=np.searchsorted(ids,z['event_id'])
    index=z['row']*npar+np.array([3*blocks[b] for b in z['truth_block']]);basis=z['basis']
    y=np.array([float(r['data']) for r in reco]);var=np.array([float(r['data_variance']) for r in reco]);keep=var>0
    nominal=np.array([float(r['value']) for r in params]);samples=[];ranks=[];bad=[];parity=0.;max_null=0.
    previous={int(b):old['sigma'][i] for i,b in enumerate(old['replica'])}
    guard_counts=[]
    for b in range(meta['requested']):
        m=np.random.default_rng(np.random.SeedSequence([meta['seed'],b])).poisson(1,len(ids))
        response=np.zeros(nrow*npar)
        for term in range(3):response+=np.bincount(index+term,weights=basis[:,term]*m[inverse],minlength=nrow*npar)
        response=response.reshape(nrow,npar);weighted=response[keep]/np.sqrt(var[keep,None])
        scale=np.linalg.norm(weighted,axis=0);scale=np.where(scale>0,scale,1.)
        u,s,vt=np.linalg.svd(weighted/scale,full_matrices=False);rank=int(np.sum(s>s[0]*config['rank_tolerance']))
        null_projection=float(np.linalg.norm(vt[rank:,pub]));max_null=max(max_null,null_projection)
        identified=null_projection<1e-8
        counts={}
        for block in sorted(blocks):
            select=(z['truth_block']==block)&keep[z['row']]
            occupied=len(np.unique(z['row'][select&(m[inverse]>0)]))
            counts[f'block_{block}_occupied_rows']=occupied
        ranks.append(dict(replica=b,rank=rank,published_identifiable=int(identified),published_null_projection=null_projection,**counts))
        if not identified:bad.append(b);continue
        fixed=np.bincount(z['fixed_row'],weights=z['fixed_weight']*m[np.searchsorted(ids,z['fixed_event_id'])],minlength=nrow)
        solution=((vt[:rank].T/s[:rank])@(u[:,:rank].T@((y[keep]-fixed[keep])/np.sqrt(var[keep]))))/scale
        if b in previous:
            delta=float(np.max(abs(solution-previous[b])));parity=max(parity,delta)
            if not np.allclose(solution,previous[b],rtol=1e-8,atol=1e-15):raise RuntimeError('Full-rank production solver mismatch')
        samples.append(solution)
    samples=np.array(samples)
    covariance_outputs(a.output,'mc_canonical_nuisance_gauge',samples)
    covariance_outputs(a.output,'mc_published',samples[:,pub])
    np.savez_compressed(a.output/'replicas.npz',sigma=samples,published_indices=np.flatnonzero(pub))
    write_csv(a.output/'replica_rank.csv',ranks)
    first=np.loadtxt(a.first_order,delimiter=',');sd=samples.std(0,ddof=1);f=np.sqrt(first.diagonal())
    write_csv(a.output/'comparison.csv',[dict(parameter_index=i,truth_block=r['truth_block'],component=r['component'],region=r['region'],
        first_order_sd=f[i],full_mc_sd=sd[i],ratio=sd[i]/f[i],mean_shift=float(samples[:,i].mean()-nominal[i]),
        identified_uncertainty=int(pub[i])) for i,r in enumerate(params)])
    counts=np.unique([r['rank'] for r in ranks],return_counts=True)
    result=dict(requested=meta['requested'],published_successful=len(samples),published_failed=bad,
        rank_counts={str(k):int(v) for k,v in zip(*counts)},max_published_null_projection=max_null,
        full_rank_production_solver_max_difference=parity,
        nuisance_covariance='canonical scaled-minimum-norm gauge only; not identified uncertainty',
        published_covariance='all event-Poisson replicas; unique published least-squares coefficients')
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
