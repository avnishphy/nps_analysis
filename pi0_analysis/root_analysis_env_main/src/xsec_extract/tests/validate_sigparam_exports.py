"""Replay cached event predictions/Jacobians/covariance/MC without ROOT.
Usage: python validate_sigparam_exports.py EXTRACTOR_OUTPUT
"""
import csv
from pathlib import Path
import sys
import numpy as np

p=Path(sys.argv[1])
def read(name):
    with (p/name).open() as f:return list(csv.DictReader(f))
pars=read('model_parameters.csv');theta=np.array([float(r['value']) for r in pars])
tau0=float(read('model_context.csv')[0]['tau0_GeV2'])
folded=read('model_reconstructed_yields.csv');n=len(folded);np_=len(theta)
mu=np.zeros(n);mc=np.zeros(n);jac=np.zeros((n,np_));components=np.zeros((n,3))
for e in read('model_event_cache.csv'):
    r=int(e['row']);basis=np.array([float(e['basis_'+k]) for k in ('U','LT','TT')])
    if e['physical']=='1':
        base=np.array([float(e['baseline_'+k]) for k in ('U','LT','TT')]);value=base*theta[[0,2,3]]
        dt=float(e['tau'])-tau0;shape=np.exp(-theta[1]*dt);value[0]*=shape
        jac[r,0]+=basis[0]*base[0]*shape;jac[r,1]+=-dt*basis[0]*value[0]
        jac[r,2:4]+=basis[1:]*base[1:]
    elif e['treatment']=='fitted_tprime_feedin':
        value=theta[4:7];jac[r,4:7]+=basis
    elif e['treatment']=='fixed_model_feedin':
        value=np.array([float(e['baseline_'+k]) for k in ('U','LT','TT')])
    else:raise AssertionError('unknown event treatment '+e['treatment'])
    components[r]+=basis*value;event=basis@value;mu[r]+=event;mc[r]+=event*event
cov=np.zeros((np_,np_))
for r in read('model_covariance.csv'):cov[int(r['i']),int(r['j'])]=float(r['covariance'])
variance=np.einsum('ri,ij,rj->r',jac,cov,jac)
def check(label,a,b):
    np.testing.assert_allclose(a,b,rtol=3e-11,atol=1e-18,equal_nan=True)
    print(label,'max_absolute=',np.nanmax(np.abs(a-b)))
check('folded',mu,np.array([float(r['prediction']) for r in folded]))
check('components',components,np.array([[float(r['component_'+k]) for k in ('U','LT','TT')] for r in folded]))
rows=read('migration_reco_rows.csv')
check('migration_prediction',mu,np.array([float(r['prediction']) for r in rows]))
check('finite_MC',mc,np.array([float(r['mc_variance_at_final']) for r in rows]))
check('parameter_variance',variance,np.array([float(r['parameter_prediction_variance']) for r in rows]))
exportjac=np.zeros_like(jac)
for r in read('model_row_jacobian.csv'):exportjac[int(r['row']),int(r['parameter'])]=float(r['derivative'])
check('Jacobian',jac,exportjac)
print('PASS: event cache reconstructs all authoritative detector exports')
