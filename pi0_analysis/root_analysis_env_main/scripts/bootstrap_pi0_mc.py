#!/usr/bin/env python3
"""Full event-Poisson MC response study with fixed real data and calibration.

One multiplier per exclusive raw event_id, shared by all reconstructed records
and all three basis terms. Generated normalization is fixed (Poissonized MC).
Production C++ kinematics/binning are compiled from the supplied configuration.
Nominal response and cell outer-product closure are mandatory before resampling.
"""
import argparse
import csv
import ctypes
import hashlib
import json
from pathlib import Path
import numpy as np
import uproot
from bootstrap_pi0_data import Fitter,write_csv,covariance_outputs


def read(p):return list(csv.DictReader(Path(p).open()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('sim','vertex','response','library','mc-library','config','first-order','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--replicas',type=int,default=2000);p.add_argument('--seed',type=int,default=20261006)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    config=json.loads(a.config.read_text())
    if config['mmiss_select']!='window':raise RuntimeError('Unsupported MC selector')
    params=read(a.response/'migration_parameters.csv');reco=read(a.response/'migration_reco_rows.csv')
    truth=read(a.response/'migration_truth_blocks.csv')
    treatment={int(r['truth_block']):r['fit_treatment'] for r in truth}
    blocks={int(r['truth_block']):int(r['active_block_index']) for r in params}
    nrow,npar=len(reco),len(params)
    if npar!=3*len(blocks) or [r['component'] for r in params]!=['U','LT','TT']*len(blocks):
        raise RuntimeError('Incompatible extraction basis/order')
    reference=np.zeros((nrow,npar))
    for r in read(a.response/'migration_design.csv'):reference[int(r['reco_row']),int(r['parameter_index'])]=float(r['response'])
    lib=ctypes.CDLL(str(a.mc_library.resolve()))
    vd=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    vi=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
    lib.pi0_mc_terms.argtypes=[ctypes.c_int,vd,vi,vi,vd]
    if lib.pi0_mc_rows()!=nrow:raise RuntimeError('Configured/response dimension mismatch')
    names=['Q2','t','tmin','xB','phi','mmiss','full_weight','sigcm']
    raw_names=['Q2i','Wi','ti','phipqi','hsxptari','hsyptari']
    with uproot.open(a.sim) as f:
        events=f[config['simc_tree']].arrays(names+['event_id','is_exclusive'],library='np')
        meta={k.split(';')[0]:f[k].member('fTitle') for k in f.keys() if any(s in k for s in ('normfac','ngen','seed'))}
        sim_uuid=str(f.file.uuid)
    exclusive=events['is_exclusive']!=0
    events={k:v[exclusive] for k,v in events.items()};ids=events['event_id'].astype(np.int64)
    with uproot.open(a.vertex) as f:
        raw=f['h10'].arrays(raw_names+['sigcm'],library='np');raw_count=f['h10'].num_entries;raw_uuid=str(f.file.uuid)
    if np.any(ids<0) or np.any(ids>=raw_count):raise RuntimeError('Unmatched event_id')
    validmass=(events['mmiss']>=config['mmiss_lower_gev'])&(events['mmiss']<=config['mmiss_upper_gev'])
    delta=np.abs(raw['sigcm'][ids].astype(float)-events['sigcm'].astype(float))
    tolerance=1e-6*np.maximum(abs(raw['sigcm'][ids].astype(float)),abs(events['sigcm'].astype(float)))+1e-20
    if np.any(delta[validmass]>tolerance[validmass]):raise RuntimeError('Raw sigcm match failed')
    inputs=np.ascontiguousarray(np.column_stack([events[k].astype(float) for k in names]+[raw[k][ids].astype(float) for k in raw_names]))
    rows=np.full(len(ids),-1,np.int32);origins=rows.copy();basis=np.zeros((len(ids),3))
    if lib.pi0_mc_terms(len(ids),inputs.ravel(),rows,origins,basis.ravel()):raise RuntimeError('C++ response construction failed')
    use=rows>=0;ids,rows,origins,basis=ids[use],rows[use],origins[use],basis[use]
    nominal_weight=events['full_weight'][use].astype(float)
    fitted=np.array([b in blocks for b in origins])
    fixed=np.array([treatment.get(int(b))=='fixed_model_feedin' for b in origins])
    if np.any(~(fitted|fixed)):raise RuntimeError('Unsupported nominal truth-block treatment')
    fitted_ids,fitted_rows,fitted_origins,fitted_basis=ids[fitted],rows[fitted],origins[fitted],basis[fitted]
    fixed_ids,fixed_rows,fixed_weights=ids[fixed],rows[fixed],nominal_weight[fixed]
    columns=np.array([3*blocks[b] for b in fitted_origins]);index=fitted_rows*npar+columns
    unique,inverse=np.unique(ids,return_inverse=True)
    fitted_inverse=inverse[fitted];fixed_inverse=inverse[fixed]
    def design(mult):
        ans=np.zeros(nrow*npar)
        for term in range(3):ans+=np.bincount(index+term,weights=fitted_basis[:,term]*mult[fitted_inverse],minlength=nrow*npar)
        return ans.reshape(nrow,npar)
    def fixed_prediction(mult):
        return np.bincount(fixed_rows,weights=fixed_weights*mult[fixed_inverse],minlength=nrow)
    nominal=design(np.ones(len(unique)))
    nominal_fixed=fixed_prediction(np.ones(len(unique)))
    if not np.allclose(nominal,reference,rtol=2e-12,atol=1e-7):
        np.savetxt(a.output/'response_difference.csv',nominal-reference,delimiter=',')
        raise RuntimeError('Rebuilt response does not match production')
    for r in read(a.response/'migration_response_cells.csv'):
        if int(r['truth_block']) not in blocks:continue
        keep=(fitted_rows==int(r['reco_row']))&(fitted_origins==int(r['truth_block']))
        cov=fitted_basis[keep].T@fitted_basis[keep]
        ref=np.array([float(r[f'cov_{i}_{j}']) for i in ('U','LT','TT') for j in ('U','LT','TT')]).reshape(3,3)
        if int(keep.sum())!=int(r['events']) or not np.allclose(cov,ref,rtol=3e-12,atol=.001):
            raise RuntimeError('Response event-count or outer-product closure failed')
    reference_fixed=np.array([float(r['fixed_prediction']) for r in reco])
    reference_fixed_var=np.array([float(r['fixed_mc_variance']) for r in reco])
    rebuilt_fixed_var=np.bincount(fixed_rows,weights=fixed_weights**2,minlength=nrow)
    if not np.allclose(nominal_fixed,reference_fixed,rtol=2e-12,atol=1e-7) or not np.allclose(rebuilt_fixed_var,reference_fixed_var,rtol=2e-12,atol=.001):
        raise RuntimeError('Fixed feed-in event sum or outer-product closure failed')
    y=np.array([float(r['data']) for r in reco]);v=np.array([float(r['data_variance']) for r in reco])
    fitter=Fitter(a.library);code,center,_=fitter.solve(nominal,y,v,nominal_fixed)
    if code or not np.allclose(center,[float(r['value']) for r in params],rtol=1e-9,atol=1e-17):
        raise RuntimeError('Nominal extraction mismatch')
    np.savez_compressed(a.output/'event_basis.npz',event_id=fitted_ids,row=fitted_rows,truth_block=fitted_origins,basis=fitted_basis,
        unique_event_id=unique,fixed_event_id=fixed_ids,fixed_row=fixed_rows,fixed_weight=fixed_weights)
    samples=[];failures=[];success=[];linear=[]
    keep=v>0;weighted=nominal[keep]/np.sqrt(v[keep,None]);scale=np.linalg.norm(weighted,axis=0)
    u,s,vt=np.linalg.svd(weighted/scale,full_matrices=False)
    op=np.zeros((npar,nrow));op[:,keep]=((vt.T/s)@u.T)/scale[:,None]/np.sqrt(v[None,keep])
    for b in range(a.replicas):
        rng=np.random.default_rng(np.random.SeedSequence([a.seed,b]))
        multiplier=rng.poisson(1,len(unique));response=design(multiplier);replica_fixed=fixed_prediction(multiplier)
        code,sigma,_=fitter.solve(response,y,v,replica_fixed)
        if code:failures.append(dict(replica=b,reason='response_solve_failed',code=code))
        else:
            samples.append(sigma);success.append(b);linear.append(-op@((response-nominal)@center+replica_fixed-nominal_fixed))
        if b%100==0:print(json.dumps(dict(replica=b,successful=len(samples))),flush=True)
    covariance_outputs(a.output,'mc',samples);covariance_outputs(a.output,'paired_linear',linear)
    np.savez_compressed(a.output/'replicas.npz',replica=success,sigma=samples,linear=linear)
    write_csv(a.output/'failures.csv',failures)
    first=np.loadtxt(a.first_order,delimiter=',');full=np.cov(np.array(samples),rowvar=False)
    fd,sd=np.sqrt(first.diagonal()),np.sqrt(full.diagonal())
    write_csv(a.output/'comparison.csv',[dict(parameter_index=i,truth_block=r['truth_block'],component=r['component'],region=r['region'],
        nominal=center[i],first_order_sd=fd[i],full_mc_sd=sd[i],ratio=sd[i]/fd[i],mc_mean_shift=float(np.mean(samples,axis=0)[i]-center[i]),
        paired_linear_sd=float(np.std(linear,axis=0,ddof=1)[i])) for i,r in enumerate(params)])
    summary=dict(requested=a.replicas,successful=len(samples),failed=len(failures),seed=a.seed,raw_events=raw_count,
        response_records=len(ids),fitted_response_records=len(fitted_ids),fixed_response_records=len(fixed_ids),
        independent_response_events=len(unique),duplicate_records=len(ids)-len(unique),
        nominal_response_max_absolute_difference=float(np.max(abs(nominal-reference))),
        nominal_fixed_feedin_max_absolute_difference=float(np.max(abs(nominal_fixed-reference_fixed))),
        sim_uuid=sim_uuid,raw_uuid=raw_uuid,simulation_metadata=meta,smearing_calibration='fixed',
        generated_normalization='fixed; Poissonized event ensemble',basis='production U,LT,TT',
        sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (a.config,a.mc_library,a.library,Path(__file__))})
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
