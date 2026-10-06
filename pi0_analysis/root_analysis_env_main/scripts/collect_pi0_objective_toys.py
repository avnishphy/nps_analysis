#!/usr/bin/env python3
"""Combine disjoint outer-experiment shards; recompute statistics from rows."""
import argparse,csv,json
from pathlib import Path
import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    rows=[];failures=[];configurations=[];summaries=[];seen=set();requested={}
    for directory in a.inputs:
        config=json.loads((directory/'configuration.json').read_text())
        configurations.append(config)
        for s in json.loads((directory/'summary.json').read_text()):
            key=(s['case'],s['mode']);requested[key]=requested.get(key,0)+s['requested']
        for row in csv.DictReader((directory/'toys.csv').open()):
            key=tuple(int(row[k]) for k in ('case','mode','experiment'))
            if key in seen:raise RuntimeError(f'Duplicate outer experiment {key}')
            seen.add(key);rows.append(row)
        if (directory/'failures.csv').exists():failures.extend(csv.DictReader((directory/'failures.csv').open()))
    for key,nrequested in sorted(requested.items()):
        selected=[r for r in rows if (int(r['case']),int(r['mode']))==key]
        if not selected:raise RuntimeError(f'No accepted outer experiments {key}')
        v=lambda name:np.array([float(r[name]) for r in selected])
        y,A=v('Y'),v('A');n=len(y);sd=y.std(ddof=1);truth=float(selected[0]['truth_Y'])
        failed=[r for r in failures if (int(r['case']),int(r['mode']))==key]
        outer_failed=sum(int(r['replica'])<0 for r in failed)
        out=dict(case=key[0],mode=key[1],truth_A=float(selected[0]['truth_A']),truth_Y=truth,
            requested=nrequested,successful=n,failed=outer_failed,outer_discarded=nrequested-n-outer_failed,
            inner_failed=len(failed)-outer_failed,mean_A=float(A.mean()),A_bias=float(A.mean()-float(selected[0]['truth_A'])),
            A_mean_se=float(A.std(ddof=1)/np.sqrt(n)),yield_bias=float(y.mean()-truth),yield_bias_se=float(sd/np.sqrt(n)),
            empirical_sd=float(sd),absolute_bias_over_sigma=float(abs(y.mean()-truth)/sd),zero_fraction=float(v('zero').mean()))
        if selected[0]['pull']:
            pull=v('pull');coverage=float(np.mean(abs(pull)<=1));boot=v('bootstrap_sd')
            out.update(mean_bootstrap_sd=float(boot.mean()),pull_mean=float(pull.mean()),pull_width=float(pull.std(ddof=1)),
                pull_mean_se=float(pull.std(ddof=1)/np.sqrt(n)),coverage=coverage,coverage_se=float(np.sqrt(coverage*(1-coverage)/n)),
                inner_successful=int(v('successful_bootstrap').sum()))
        summaries.append(out)
    for name,table in [('toys.csv',rows),('failures.csv',failures)]:
        if table:
            with (a.output/name).open('w') as f:
                w=csv.DictWriter(f,fieldnames=table[0]);w.writeheader();w.writerows(table)
    (a.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (a.output/'configuration.json').write_text(json.dumps(dict(inputs=[str(p) for p in a.inputs],shards=configurations),indent=2)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()
