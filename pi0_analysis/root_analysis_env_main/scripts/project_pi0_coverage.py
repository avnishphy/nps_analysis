#!/usr/bin/env python3
"""Forward-fold one toy signal mode through the actual response/production solve.

This is conditional one-direction coverage, not full detector coverage. Other
truth components and response stay fixed. It demonstrates that inversion does
not cure upstream subtraction bias.
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import Fitter,write_csv


def read(p):return list(csv.DictReader(p.open()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('coverage','response','library','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--parameter',type=int,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    pars=read(a.response/'migration_parameters.csv');reco=read(a.response/'migration_reco_rows.csv')
    A=np.zeros((len(reco),len(pars)))
    for r in read(a.response/'migration_design.csv'):A[int(r['reco_row']),int(r['parameter_index'])]=float(r['response'])
    center=np.array([float(r['value']) for r in pars]);v=np.array([float(r['data_variance']) for r in reco]);f=Fitter(a.library)
    fixed=np.array([float(r['fixed_prediction']) for r in reco])
    ytrue=fixed+A@center;out=[];maxdifference=0.
    for r in read(a.coverage/'pseudoexperiments.csv'):
        unit=center[a.parameter]/float(r['true_yield']);y=ytrue+A[:,a.parameter]*unit*(float(r['estimate'])-float(r['true_yield']))
        code,p,_=f.solve(A,y,v,fixed)
        if code:raise RuntimeError('Forward coverage solve failed')
        sd=abs(unit)*float(r['bootstrap_sd']);pull=(p[a.parameter]-center[a.parameter])/sd
        expected=np.sign(unit)*float(r['pull']);maxdifference=max(maxdifference,abs(pull-expected))
        out.append(dict(case=r['case'],experiment=r['experiment'],parameter=a.parameter,true=center[a.parameter],estimate=p[a.parameter],sigma=sd,pull=pull,covered=int(abs(pull)<=1)))
    write_csv(a.output/'pseudoexperiments.csv',out)
    result=[]
    for case in sorted({r['case'] for r in out}):
        rows=[r for r in out if r['case']==case];pull=np.array([r['pull'] for r in rows])
        result.append(dict(case=case,parameter=a.parameter,bias=float(np.mean([r['estimate']-r['true'] for r in rows])),
            pull_width=float(pull.std(ddof=1)),coverage=float(np.mean(abs(pull)<=1)),replicas=len(rows)))
    (a.output/'summary.json').write_text(json.dumps(dict(cases=result,max_yield_vs_cross_section_pull_difference=maxdifference),indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
