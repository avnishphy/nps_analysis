#!/usr/bin/env python3
"""Independent optimizer and paired fixed/free shape diagnostics (toys only)."""
import argparse,csv,ctypes,json
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution,minimize_scalar
from study_pi0_objectives import PrimitiveFit,write


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',type=Path,required=True);p.add_argument('--free',type=Path,required=True)
    p.add_argument('--fixed',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--checks',type=int,default=40)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    free=list(csv.DictReader((a.free/'toys.csv').open()));fixed=list(csv.DictReader((a.fixed/'toys.csv').open()))
    byid={int(r['experiment']):r for r in fixed}
    diff=np.array([float(r['Y'])-float(byid[int(r['experiment'])]['Y']) for r in free])
    width=np.array([float(r['width']) for r in free]);turn=np.array([float(r['turn']) for r in free])
    fit=PrimitiveFit(a.library);fn=fit.lib.pi0_primitive_deviance
    v=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS');fn.argtypes=[ctypes.c_int,v,v,v];fn.restype=ctypes.c_double
    x=(np.arange(200)+.5)*.002;side=((x>=.01)&(x<=.11))|((x>=.15)&(x<=.4))
    f=1/(1+np.exp((x-.145)/.025));signal=np.exp(-.5*((x-.135)/.006)**2);signal[side]=0;signal*=300/signal.sum()
    c=np.array([1,-1/6,-1/12,-1/12,1/18,1/18]);rates=2*np.array([1,6*.4,12*.35,12*.35,18*.05,18*.05])
    means=np.repeat(rates[:,None],200,axis=1);means[0]+=.5*f+signal
    checks=[];profile=[]
    for exp in range(a.checks):
        rng=np.random.default_rng(np.random.SeedSequence([7312026,5,exp]));counts=rng.poisson(means).astype(float)
        code,info=fit.fit(counts,c);upper=max(10,20*max(1,(c@counts)[side].max()))
        flat=counts.ravel();objective=lambda par:fn(6,flat,c,np.ascontiguousarray(par,dtype=float))
        result=differential_evolution(objective,[(0,upper),(.11,.22),(.001,.1)],seed=exp+171,popsize=12,maxiter=250,tol=1e-9,polish=True)
        bg=lambda par:float(np.sum(par[0]/(1+np.exp((x-par[1])/par[2]))))
        checks.append(dict(experiment=exp,original_code=code,original_objective=info[3],independent_objective=result.fun,
            independent_success=int(result.success),objective_improvement=info[3]-result.fun,original_background=info[7],
            independent_background=bg(result.x),background_difference=bg(result.x)-info[7]))
        if exp==0:
            for t in np.linspace(.11,.22,221):
                res=minimize_scalar(lambda A:objective([A,t,.001]),bounds=(0,upper),method='bounded',options={'xatol':1e-10})
                profile.append(dict(turn=t,width=.001,A=res.x,objective=res.fun,background=bg([res.x,t,.001])))
        if exp%10==0:print(json.dumps(checks[-1]),flush=True)
    result=dict(paired_toys=len(diff),mean_free_minus_fixed_yield=float(diff.mean()),paired_difference_se=float(diff.std(ddof=1)/np.sqrt(len(diff))),
        minimum_width_fraction=float(np.mean(width<.00101)),maximum_width_fraction=float(np.mean(width>.09999)),
        minimum_turn_fraction=float(np.mean(turn<.11001)),maximum_turn_fraction=float(np.mean(turn>.21999)),
        independent_checks=checks)
    write(a.output/'independent_optimizer.csv',checks);write(a.output/'turn_profile.csv',profile)
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
