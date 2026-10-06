#!/usr/bin/env python3
"""Check the primitive dual against direct constrained likelihood and on/off."""
import argparse,ctypes,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from study_pi0_objectives import PrimitiveFit


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();fit=PrimitiveFit(a.library)
    fn=fit.lib.pi0_primitive_deviance;v=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    fn.argtypes=[ctypes.c_int,v,v,v];fn.restype=ctypes.c_double
    c=np.array([1,-1/6,-1/12,-1/12,1/18,1/18]);pair=c[:2].copy()
    rng=np.random.default_rng(18472);par=np.array([2.,.145,.025]);errors=[]
    for _ in range(100):
        counts=rng.poisson(rng.uniform(0,8,size=(2,200))).astype(float)
        six=np.vstack([counts,np.zeros((4,200))]);par[0]=rng.uniform(0,20)
        errors.append(abs(fn(2,counts.ravel(),pair,par)-fn(6,six.ravel(),c,par)))
    code,info=fit.fit(np.zeros((6,200)),c)
    assert code==0 and info[0]==0 and info[4]==1 and info[7]==0
    checks=[];bin_index=40;f=1/(1+np.exp((((bin_index+.5)*.002)-.145)/.025))
    count_cases=[np.zeros(6),np.array([0,5,0,0,0,0]),np.array([3,0,0,0,0,0]),
                 np.array([0,0,0,0,3,2])]+[rng.poisson([2,5,3,3,1,1]).astype(float) for _ in range(20)]
    for counts in count_cases:
        for mu in (0.,.1,3.,20.):
            par=np.array([mu/f,.145,.025]);mass=np.zeros((6,200));mass[:,bin_index]=counts
            dual=(fn(6,mass.ravel(),c,par)-fn(6,np.zeros(1200),c,par))/2+mu
            positive=counts>0
            def nll(lam):
                if np.any(lam[positive]<=0):return 1e100
                return float((lam-counts).sum()+np.sum(counts[positive]*np.log(counts[positive]/lam[positive])))
            start=counts.copy()+.1;delta=mu-c@start
            if delta>=0:start[0]+=delta
            else:start[1]-=6*delta
            res=minimize(nll,start,method='SLSQP',bounds=[(1e-12 if n>0 else 0,None) for n in counts],
                constraints=[dict(type='eq',fun=lambda z:c@z-mu)],options={'ftol':1e-10,'maxiter':1000})
            checks.append(dict(counts=counts.tolist(),mu=mu,dual=dual,direct=float(res.fun),
                difference=float(res.fun-dual),success=bool(res.success),constraint=float(c@res.x-mu)))
    result=dict(zero_branch=True,onoff_max_absolute_difference=max(errors),direct_profile_checks=checks,
        direct_max_absolute_difference=max(abs(r['difference']) for r in checks))
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    assert max(errors)<1e-8,result
    assert all(r['success'] and abs(r['constraint'])<1e-7 and abs(r['difference'])<2e-6 for r in checks),result
    print(json.dumps({k:v for k,v in result.items() if k!='direct_profile_checks'}))


if __name__=='__main__':main()
