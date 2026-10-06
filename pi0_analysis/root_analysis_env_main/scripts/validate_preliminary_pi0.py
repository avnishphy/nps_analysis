#!/usr/bin/env python3
"""Small numerical/production checks for the frozen preliminary M0 campaign."""
import argparse
from preliminary_pi0_xsec import *

OLD_DATA = REPO/'validation/xsec_input_recovery_20261004/output/KinC_x36_4/root/combined_branches_LH2.root'

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=args.output
    data=Data(out);model=Model(out);central=json.loads((out/'central_results.json').read_text())['additive']
    assert model.nparams==7 and model.parameter_names==[
        'N_U','DeltaB_U','N_LT','N_TT','feedin_tprime_below_U','feedin_tprime_below_LT','feedin_tprime_below_TT']
    assert len(model.nuisance_blocks)==1 and not set(model.nuisance_blocks)&set(model.fixed_blocks)
    expected=np.loadtxt(out/'accepted_runs.txt',dtype=int)
    np.testing.assert_array_equal(expected,data.manifest['run_number'])
    with uproot.open(OLD_DATA) as f:
        prior=f['physics'].arrays(library='np');prior_runs=f['analysis_runs'].arrays(library='np')
    for k in prior:np.testing.assert_array_equal(data.e[k],prior[k])
    for k in prior_runs:np.testing.assert_array_equal(data.manifest[k],prior_runs[k])
    weight=data.e['pi0_timing_coeff']-(data.e['pi0_timing_bin_mean']-data.e['pi0_weight'])
    regular=(data.e['mpi0_all']>=0)&(data.e['mpi0_all']<.4)
    audit=json.loads((out/'timing_audit/audit.json').read_text())
    save(out/'production_verification.json',dict(runs=len(expected),events=len(weight),exposure_uC=data.charge,
        applicable_events=int(regular.sum()),undefined_applicable=int((~np.isfinite(weight[regular])).sum()),
        undefined_flow_events=int((~np.isfinite(weight[~regular])).sum()),legacy_per_run_branches_exact=41,
        combined_legacy_branches_exact=True,empty_bin_residual=audit['empty_nonzero_residual_sum'],empty_bins=audit['empty_nonzero_residual_bins']))
    y,v=data.central('additive');start=np.array(central['theta']);tests=[]
    for i in range(4):
        initial=start.copy()
        if i==1:initial[:4]=[.7*start[0],start[1]+.5,-start[2],.7*start[3]]
        if i==2:initial[:4]=[1.5*start[0],start[1]-.5,start[2],-start[3]]
        if i==3:initial[:4]=[1,0,1,1];initial[4:]=0;initial[4::3]=1e-7
        r=model.constrained(y,v,initial);tests.append(dict(start=i,**serial(r)))
        save(out/'start_stability.json',tests)
        print('START',i,r['code'],r['info'][0],r['theta'][:4],flush=True)
    old=model.constrained(*data.central('old'))
    save(out/'boundary_solver_old_parity.json',serial(old))
    assert old['code']==0
    np.testing.assert_allclose(old['theta'][:4],model.old[:4],rtol=3e-4,atol=2e-5)
    # Nominal full data-estimator reconstruction must reproduce stored additive
    # weights and the production mass-fit inputs, including all 56 refits.
    yy=np.zeros(60);vv=np.zeros(60);fail=[];residual=0.
    for run in data.runs:
        m=np.ones(len(run.coefficient));t,var,n=run.spectra(m);code,s,info=data.fitter.fit(t,var)
        if code:fail.append(run.run)
        bn=np.divide(t-s,n,out=np.zeros(200),where=n>0);take=run.selected
        w=(run.coefficient[take]-bn[run.mass_bin[take]])*run.factor[take]
        yy+=np.bincount(run.row[take],weights=w,minlength=60)
        vv+=np.bincount(run.row[take],weights=w*w,minlength=60)
        residual+=s[n==0].sum()
    np.testing.assert_allclose(yy,y,atol=1e-11,rtol=1e-10)
    np.testing.assert_allclose(vv,v,atol=1e-12,rtol=1e-10)
    assert not fail
    # Rank uses exactly the fixed preliminary mask, including in MC ensembles.
    j=np.array(central['jacobian']);var=np.array(central['variance'])
    assert j.shape==(model.nrows,7)
    a=j[model.mask]/np.sqrt(var[model.mask,None]);scale=1/np.linalg.norm(a,axis=0);a*=scale
    _,s,_=np.linalg.svd(a,full_matrices=False);un,sn,_=np.linalg.svd(a[:,4:],full_matrices=False)
    nr=int(sum(sn>sn[0]*1e-10));physics=a[:,:4]-un[:,:nr]@(un[:,:nr].T@a[:,:4])
    save(out/'central_rank.json',dict(detector=int(sum(s>s[0]*1e-10)),nuisance=nr,
        physics=int(np.linalg.matrix_rank(physics,tol=1e-10)),included_rows=int(model.mask.sum()),jacobian_rows=int(model.mask.sum()),
        jacobian_columns=j.shape[1],parameter_order=model.parameter_names,singular_values=s.tolist(),condition=float(s[0]/s[-1])))
    save(out/'complete_estimator_parity.json',dict(mass_fit_failures=fail,
        max_row_yield_difference=float(np.max(abs(yy-y))),max_row_variance_difference=float(np.max(abs(vv-v))),
        empty_bin_residual=float(residual),fixed_ellipse=True,exposure_uC=data.charge))
    factors={r.run:float(r.factor[0]) for r in data.runs}
    empty=[r for r in read(out/'timing_audit/timing_weight_inclusive_closure.csv') if 'UNREPRESENTABLE' in r['status']]
    bound=sum(abs(float(r['S_rb']))*factors[int(r['run'])] for r in empty)
    save(out/'empty_support_summary.json',dict(empty_bins=len(empty),signed_fitted_counts=sum(float(r['S_rb']) for r in empty),
        absolute_normalized_bound_per_mC=bound,fraction_of_selected_yield=bound/y.sum()))

if __name__=='__main__':
    import os,sys,traceback
    try:main();code=0
    except BaseException:traceback.print_exc();code=1
    sys.stdout.flush();sys.stderr.flush();os._exit(code)
