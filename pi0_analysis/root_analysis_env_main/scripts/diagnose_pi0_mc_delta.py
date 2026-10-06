#!/usr/bin/env python3
"""Compare prediction-only delta propagation with full normal-equation derivative."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import write_csv


def read(p):return list(csv.DictReader(p.open()))
def correlation(c):
    s=np.sqrt(c.diagonal());return c/np.outer(s,s)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('mc','rank-aware','response','first-order','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    z=np.load(a.mc/'event_basis.npz');params=read(a.response/'migration_parameters.csv');rows=read(a.response/'migration_reco_rows.csv')
    npar=len(params);design=np.zeros((len(rows),npar))
    for r in read(a.response/'migration_design.csv'):design[int(r['reco_row']),int(r['parameter_index'])]=float(r['response'])
    blocks={int(r['truth_block']):int(r['active_block_index']) for r in params}
    sigma=np.array([float(r['value']) for r in params]);y=np.array([float(r['data']) for r in rows]);v=np.array([float(r['data_variance']) for r in rows]);keep=v>0
    fixed=np.array([float(r['fixed_prediction']) for r in rows])
    aa=design[keep]/np.sqrt(v[keep,None]);scale=np.linalg.norm(aa,axis=0);u,s,vt=np.linalg.svd(aa/scale,full_matrices=False)
    hinv=((vt.T/s**2)@vt)/np.outer(scale,scale)
    op=np.zeros((npar,len(rows)));op[:,keep]=((vt.T/s)@u.T)/scale[:,None]/np.sqrt(v[None,keep])
    residual=y-fixed-design@sigma
    total_events=len(z['event_id'])+len(z['fixed_event_id'])
    grad=np.zeros((total_events,npar));prediction_grad=grad.copy()
    for b,active in blocks.items():
        select=(z['truth_block']==b)&keep[z['row']];row=z['row'][select];basis=z['basis'][select];sl=slice(3*active,3*active+3)
        fitted_index=np.flatnonzero(select)
        prediction_grad[fitted_index]=-op[:,row].T*(basis@sigma[sl])[:,None]
        grad[fitted_index]=prediction_grad[fitted_index]+(basis@hinv[:,sl].T)*(residual[row]/v[row])[:,None]
    offset=len(z['event_id']);fixed_rows=z['fixed_row'];fixed_weights=z['fixed_weight']
    prediction_grad[offset:]=-op[:,fixed_rows].T*fixed_weights[:,None]
    grad[offset:]=prediction_grad[offset:]
    # The current input has no duplicate simulation records; refuse accidental
    # independent treatment if another producer introduces them.
    all_ids=np.concatenate([z['event_id'],z['fixed_event_id']])
    if len(np.unique(all_ids))!=len(all_ids):raise RuntimeError('Group duplicate generator events before covariance construction')
    cov=grad.T@grad;pred=prediction_grad.T@prediction_grad;previous=np.loadtxt(a.first_order,delimiter=',')
    if not np.allclose(pred,previous,rtol=1e-8,atol=1e-25):raise RuntimeError('Prediction-only delta provenance mismatch')
    np.savetxt(a.output/'full_normal_equation_delta_covariance.csv',cov,delimiter=',',fmt='%.17g')
    full=np.loadtxt(a.rank_aware/'mc_canonical_nuisance_gauge_covariance.csv',delimiter=',');cp,cf=correlation(previous),correlation(full)
    table=[]
    for i,r in enumerate(params):
        table.append(dict(parameter_index=i,truth_block=r['truth_block'],component=r['component'],region=r['region'],
            prediction_delta_sd=np.sqrt(pred[i,i]),complete_delta_sd=np.sqrt(cov[i,i]),full_mc_sd=np.sqrt(full[i,i]),
            full_over_complete_delta=np.sqrt(full[i,i]/cov[i,i])))
    write_csv(a.output/'delta_comparison.csv',table)
    pub=np.array([r['region']=='published' for r in params]);pairs=[]
    for i in np.flatnonzero(pub):
        for j in np.flatnonzero(pub):
            if i>j:pairs.append(dict(i=i,j=j,first_order_correlation=cp[i,j],full_mc_correlation=cf[i,j],difference=cf[i,j]-cp[i,j]))
    pairs.sort(key=lambda r:abs(r['difference']),reverse=True);write_csv(a.output/'correlation_comparison.csv',pairs)
    support=[]
    for b,active in blocks.items():
        select=(z['truth_block']==b)&keep[z['row']];weights=z['basis'][select,0]
        support.append(dict(truth_block=b,events=int(select.sum()),occupied_rows=int(len(np.unique(z['row'][select]))),
            effective_U_events=float(weights.sum()**2/(weights@weights))))
    write_csv(a.output/'truth_block_support.csv',support)
    samples=np.load(a.rank_aware/'replicas.npz')['sigma'];ranks=read(a.rank_aware/'replica_rank.csv');groups=[]
    for rank in sorted({int(r['rank']) for r in ranks}):
        selection=np.array([int(r['rank'])==rank for r in ranks])
        for i in np.flatnonzero(pub):groups.append(dict(rank=rank,replicas=int(selection.sum()),parameter_index=i,mean=samples[selection,i].mean(),sd=samples[selection,i].std(ddof=1)))
    write_csv(a.output/'rank_conditioned_means.csv',groups)
    result=dict(max_published_full_vs_prediction_relative_sd=float(np.max(abs(np.sqrt(full.diagonal()[pub]/pred.diagonal()[pub])-1))),
        max_published_full_vs_complete_delta_relative_sd=float(np.max(abs(np.sqrt(full.diagonal()[pub]/cov.diagonal()[pub])-1))),
        maximum_published_correlation_difference=float(max(abs(r['difference']) for r in pairs)),
        expected_rank_loss_for_three_guard_rows_with_1_3_1_events=float(1-(1-np.exp(-1))**2*(1-np.exp(-3))),
        omitted_first_order_term='H^-1 delta(A)^T W residual; evaluated here separately',
        recommendation='Use all-replica identified published MC covariance, not covariance conditioned on full nuisance rank; overall release still requires coverage resolution')
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
