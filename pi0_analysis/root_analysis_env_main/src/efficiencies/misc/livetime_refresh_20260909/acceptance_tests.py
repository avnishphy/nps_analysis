"""Conditional pulser-vs-physics acceptance and correlated block sensitivity.

Null: pulser labels exchangeable among selected trigger inputs; the recorded
subset is fixed. This is a diagnostic model, not a proof of random pulser phase.
"""
from pathlib import Path
import csv,json,math,numpy as np
P=Path(__file__).resolve().parent;out=[]
def exact_label_p(E,N,D,S):
    lo=max(0,N-(S-D));hi=min(N,D)
    # Constants cancel after normalization; avoids optional SciPy ABI/version coupling.
    def logchoose(n,k):return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
    logw=np.array([logchoose(D,k)+logchoose(S-D,N-k) for k in range(lo,hi+1)])
    obs=logw[E-lo];w=np.exp(logw-logw.max())
    return float(w[logw<=obs+1e-7].sum()/w.sum())
assert abs(exact_label_p(8,10,9,16)-0.034965034965034975)<1e-10
for f in sorted((P/'results').glob('run*.json')):
    r=json.loads(f.read_text())
    if r['status']!='CALCULATED' or r['run_type']!='production':continue
    run=r['run'];p=r['prescale_factor'];n=r['nominal']
    rows=[{k:float(v) for k,v in x.items()} for x in csv.DictReader((P/f'results/run{run}_intervals.csv').open()) if x['nominal']=='1']
    for variant in ['tight','raw']:
        E=int(n['E_'+variant]);N=int(n['N']);D=int(n['D']);S=int(n['S']);cells=[[E,N-E],[D-E,S-D-N+E]]
        Lp=p*E/D;Lc=p*(N-E)/(S-D);mu=N*D/S;var=N*(D/S)*(1-D/S)*(S-N)/(S-1)
        valid=min(x for a in cells for x in a)>=0 and var>0
        z=(E-mu)/math.sqrt(var) if valid else None
        pv=exact_label_p(E,N,D,S) if valid else None
        rng=np.random.default_rng(run);blocks={}
        for x in rows:
            key=(int(x['component']),int(x['t_end']//20));blocks.setdefault(key,np.zeros(4));blocks[key]+=np.array([x['N'],x['E_'+variant],x['S'],x['D']])
        a=np.array(list(blocks.values()));s=a[rng.integers(0,len(a),size=(2000,len(a)))].sum(axis=1)
        diff=p*s[:,1]/s[:,3]-p*(s[:,0]-s[:,1])/(s[:,2]-s[:,3]);sigma=float(np.std(diff,ddof=1))
        out.append(dict(run=run,variant=variant,p=p,N=N,E=E,S=S,D=D,EDTM=Lp,physics_CLT=Lc,difference=Lp-Lc,exchangeable_label_z=z,fisher_two_sided_p=pv,block20_difference_sigma=sigma,block20_z=(Lp-Lc)/sigma if sigma else None,counter_suspect_intervals=sum(x['N']-x['A']>4 for x in rows)))
(P/'acceptance_tests.json').write_text(json.dumps(out,indent=2)+'\n')
with (P/'acceptance_tests.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
for r in out:
    if r['variant']=='raw':print(r['run'],'p',r['p'],'diff',round(r['difference'],7),'HG_z',None if r['exchangeable_label_z'] is None else round(r['exchangeable_label_z'],2),'block_z',None if r['block20_z'] is None else round(r['block20_z'],2),'counter flags',r['counter_suspect_intervals'])
