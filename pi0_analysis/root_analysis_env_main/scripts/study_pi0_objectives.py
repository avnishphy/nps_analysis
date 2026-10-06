#!/usr/bin/env python3
"""Primitive-count candidate screening and nested coverage; no production edits."""
import argparse,csv,ctypes,hashlib,json,time
from pathlib import Path
import numpy as np


class PrimitiveFit:
    def __init__(self,library):
        self.lib=ctypes.CDLL(str(Path(library).resolve()))
        v=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        self.lib.pi0_fit_primitive.argtypes=[ctypes.c_int,v,v,ctypes.c_int,ctypes.c_int,v]
    def fit(self,counts,coeff,mode=0,fixed=False):
        out=np.zeros(8)
        code=self.lib.pi0_fit_primitive(len(coeff),np.ascontiguousarray(counts,dtype=float).ravel(),np.ascontiguousarray(coeff,dtype=float),mode,int(fixed),out)
        return code,out


def write(path,rows):
    if not rows:return
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--experiments',type=int,default=1000);p.add_argument('--bootstrap',type=int,default=0)
    p.add_argument('--start',type=int,default=0,help='First independent outer experiment; allows reproducible shards')
    p.add_argument('--amplitudes',help='Optional comma-separated production-scale truth amplitudes for boundary scan')
    p.add_argument('--fixed-shape',action='store_true');p.add_argument('--six-categories',action='store_true')
    p.add_argument('--modes',default='0,1,2');p.add_argument('--cases',default='0,1,2,3,4')
    p.add_argument('--seed',type=int,default=7312026)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);fit=PrimitiveFit(a.library)
    a.library_sha256=hashlib.sha256(a.library.read_bytes()).hexdigest()
    a.driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    x=(np.arange(200)+.5)*.002;shape=1/(1+np.exp((x-.145)/.025))
    signal=np.exp(-.5*((x-.135)/.006)**2)
    side=((x>=.01)&(x<=.11))|((x>=.15)&(x<=.4));signal[side]=0;signal/=signal.sum()
    cases=[(0.,3637.1667,.3766667),(.5,3637.1667,.3766667),(12.,3637.1667,.3766667),(24.,3637.1667,.3766667),(12.,300.,2.),(.5,300.,2.)]
    if a.amplitudes:cases=[(float(v),3637.1667,.3766667) for v in a.amplitudes.split(',')]
    coeff=np.array([1,-1/6,-1/12,-1/12,1/18,1/18] if a.six_categories else [1,-1/6])
    rows=[];summary=[];failures=[]
    for case in map(int,a.cases.split(',')):
        amplitude,ytrue,acc=cases[case]
        if a.six_categories:
            # Same signed contrast; nuisance rates remain independent counts.
            rates=np.array([acc,6*.4*acc,12*.35*acc,12*.35*acc,18*.05*acc,18*.05*acc])
        else:rates=np.array([acc,6*acc])
        means=np.repeat(rates[:,None],200,axis=1);means[0]+=amplitude*shape+ytrue*signal
        for mode in map(int,a.modes.split(',')):
            result=[];failed=0;innerfailed=0;t0=time.time()
            for exp in range(a.start,a.start+a.experiments):
                rng=np.random.default_rng(np.random.SeedSequence([a.seed,case,exp]))
                counts=rng.poisson(means).astype(float);code,info=fit.fit(counts,coeff,mode,a.fixed_shape)
                if code:
                    failed+=1;failures.append(dict(case=case,mode=mode,experiment=exp,replica=-1,code=code,status=info[5]));continue
                yy=float(np.sum(coeff@counts)-info[7]);sd=pull=coverage=None;nsuccess=0
                if a.bootstrap:
                    samples=[]
                    for b in range(a.bootstrap):
                        replica=rng.poisson(counts).astype(float);code,br=fit.fit(replica,coeff,mode,a.fixed_shape)
                        if code:
                            innerfailed+=1;failures.append(dict(case=case,mode=mode,experiment=exp,replica=b,code=code,status=br[5]));continue
                        samples.append(float(np.sum(coeff@replica)-br[7]))
                    nsuccess=len(samples)
                    if nsuccess<a.bootstrap*.9:continue
                    sd=float(np.std(samples,ddof=1));pull=(yy-ytrue)/sd;coverage=int(abs(pull)<=1)
                row=dict(case=case,mode=mode,experiment=exp,truth_A=amplitude,truth_Y=ytrue,A=info[0],turn=info[1],width=info[2],
                    Y=yy,zero=int(info[4]),bootstrap_sd=sd,pull=pull,covered=coverage,successful_bootstrap=nsuccess)
                rows.append(row);result.append(row)
                if exp%100==0:print(json.dumps(dict(case=case,mode=mode,experiment=exp,elapsed=time.time()-t0)),flush=True)
            av=np.array([r['A'] for r in result]);yv=np.array([r['Y'] for r in result]);spread=float(yv.std(ddof=1))
            stats=dict(case=case,mode=mode,fixed_shape=a.fixed_shape,six_categories=a.six_categories,truth_A=amplitude,truth_Y=ytrue,
                requested=a.experiments,successful=len(result),failed=failed,inner_failed=innerfailed,
                mean_A=float(av.mean()),A_bias=float(av.mean()-amplitude),A_mean_se=float(av.std(ddof=1)/np.sqrt(len(av))),
                yield_bias=float(yv.mean()-ytrue),yield_bias_se=spread/np.sqrt(len(yv)),empirical_sd=spread,
                absolute_bias_over_sigma=float(abs(yv.mean()-ytrue)/spread),zero_fraction=float(np.mean([r['zero'] for r in result])))
            if a.bootstrap:
                pulls=np.array([r['pull'] for r in result]);coverage=float(np.mean(abs(pulls)<=1))
                stats.update(mean_bootstrap_sd=float(np.mean([r['bootstrap_sd'] for r in result])),pull_mean=float(pulls.mean()),pull_width=float(pulls.std(ddof=1)),coverage=coverage,coverage_se=float(np.sqrt(coverage*(1-coverage)/len(result))))
            summary.append(stats);write(a.output/'toys.csv',rows);write(a.output/'failures.csv',failures)
            (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(stats),flush=True)
    (a.output/'configuration.json').write_text(json.dumps(vars(a),default=str,indent=2)+'\n')


if __name__=='__main__':main()
