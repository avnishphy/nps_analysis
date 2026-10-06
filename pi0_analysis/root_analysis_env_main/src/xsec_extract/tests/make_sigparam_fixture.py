"""Synthetic varying-vertex closure. Baseline evaluated by the production C++ header.

Usage: python make_sigparam_fixture.py CONFIG NEW_DIRECTORY SIGPARAM_EVALUATE_BINARY
The run-1 manifest describes only synthetic events; no real-data gate is bypassed.
"""
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np
import uproot

config_path, out, evaluator = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
config=json.loads(Path(config_path).read_text())
out.mkdir(parents=True,exist_ok=False)
mp,mpi=.9382720813,.1349768
def forward(q,w):
    e=(w*w-mp*mp-q)/(2*w); ep=(w*w+mpi*mpi-mp*mp)/(2*w)
    return mpi*mpi-q-2*(e*ep-np.sqrt(e*e+q)*np.sqrt(ep*ep-mpi*mpi))
edges=np.array(config['tprime_bin_edges']);centers=(edges[1:]+edges[:-1])/2
nr=len(centers);nphi=config['phi_bins']
it,ip,origin,repeat=np.meshgrid(np.arange(nr),np.arange(nphi),np.arange(nr),np.arange(20),indexing='ij')
it,ip,origin,repeat=(x.ravel() for x in (it,ip,origin,repeat));n=len(it)
idx=np.arange(n);phase=(repeat-9.5)/9.5
q=(4.+.16*phase).astype('float32').astype(float);xb=.36+.003*np.sin(idx*.71)
w=np.sqrt(mp*mp+q*(1/xb-1)).astype('float32').astype(float)
recw=float(np.float32(math.sqrt(mp*mp+4*(1/.36-1))));tf=forward(4.,recw)
sigcm=np.full(n,1e-6,dtype='float32')
raw={'phipqi':(2*np.pi*(ip+.5)/nphi).astype('float32'),'sigcm':sigcm.copy()}
sim={k:np.asarray(v,dtype='float32') for k,v in dict(Q2=np.full(n,4.),xB=np.full(n,.36),W=np.full(n,recw),t=tf+centers[it],tmin=np.full(n,tf),phi=raw['phipqi'],mmiss=np.full(n,.95),full_weight=1e6*np.where(origin==it,1.,.04)*sigcm,sigcm=sigcm,epsilon_i=np.full(n,.53)).items()}
sim.update(event_id=idx.astype('uint64'),is_exclusive=np.ones(n,dtype='int32'))
dit=np.repeat(np.arange(nr),nphi*100);dip=np.tile(np.repeat(np.arange(nphi),100),nr);nd=len(dit)
data={k:np.asarray(v,dtype='float64') for k,v in dict(Q2=np.full(nd,4.),xB=np.full(nd,.36),W=np.full(nd,recw),t=tf+centers[dit],tmin=np.full(nd,tf),phi=2*np.pi*(dip+.5)/nphi,mmiss_all=np.full(nd,.95)).items()}
data.update(scale=np.ones(nd,dtype='float32'),charge_uC=np.ones(nd,dtype='float32'),run_number=np.ones(nd,dtype='int32'))
manifest={'run_number':np.array([1],dtype='int32'),'success':np.array([1],dtype='int32'),'charge_uC':np.array([1],dtype='float32'),'scale':np.array([1],dtype='float32')}
tp=centers[origin]+.009*phase
t=(forward(q,w)+tp)
raw['Q2i']=q.astype('float32');raw['Wi']=w.astype('float32');raw['ti']=(-t).astype('float32')
t=-raw['ti'].astype(float);tp=t-forward(q,w)
raw['hsxptari']=(.003*phase).astype('float32');raw['hsyptari']=(.002*np.sin(idx)).astype('float32')
angle=math.radians(config['hms_theta_deg'])
# Same generated-track geometry as vertex_epsilon_from_exclusive_simc.
xp=raw['hsxptari'].astype(float);yp=raw['hsyptari'].astype(float)
ct=(math.cos(angle)+yp*math.sin(angle))/np.sqrt(1+xp*xp+yp*yp)
nu=(w*w+q-mp*mp)/(2*mp);eps=1/(1+2*(1+nu*nu/q)*(1-ct)/(1+ct))
records=''.join(f'{a:.17g} {b:.17g} {c:.17g} {d:.17g} {e:.17g}\n' for a,b,c,d,e in zip(q,w,t,tp,eps))
baseline=np.loadtxt(subprocess.run([evaluator],input=records,text=True,capture_output=True,check=True).stdout.splitlines())
theta=np.array([1.2,.4,-.8,.5])
response_weight=sim['full_weight'].astype(float)/sim['sigcm'].astype(float)
tau0=float(np.sum(response_weight*(-tp))/np.sum(response_weight))
coeff=baseline*theta[[0,2,3]]
coeff[:,0]*=np.exp(-theta[1]*(-tp-tau0))
phi=raw['phipqi'].astype(float)
basis=np.column_stack([np.ones(n),np.sqrt(2*eps*(1+eps))*np.cos(phi),eps*np.cos(2*phi)])/(2*np.pi)
event_yield=np.sum(basis*coeff,axis=1)*sim['full_weight'].astype(float)/sim['sigcm'].astype(float)
# Original fixture order: reconstructed t,phi,origin,repeat.
nr=len(config['tprime_bin_edges'])-1; nphi=config['phi_bins']
y=event_yield.reshape(nr*nphi,n//(nr*nphi)).sum(axis=1)
assert np.all(y>0)
data['pi0_weight']=np.repeat(y*config['tgt_contam']/100,100)
def write(path,name,arrays,manifest=None):
    with uproot.recreate(path) as f:
        f['fixture_provenance']='SYNTHETIC varying-kinematics SigParam2021 closure; not production data'
        f.mktree(name,{k:v.dtype for k,v in arrays.items()}).extend(arrays)
        if manifest is not None:f.mktree('analysis_runs',{k:v.dtype for k,v in manifest.items()}).extend(manifest)
if len(sys.argv)>4:
    if sys.argv[4]=='matching-failures':
        sim['event_id'][1]=sim['event_id'][0];sim['event_id'][2]=n+100
        sim['sigcm'][3]*=2;raw['hsxptari'][4]=np.nan
    elif sys.argv[4]=='empty-row':data={k:v[100:] for k,v in data.items()}
    else:raise ValueError('Unknown fixture variant')
write(out/'raw.root','h10',raw);write(out/'sim.root','simulation',sim);write(out/'data.root','physics',data,manifest)
(out/'truth.json').write_text(json.dumps({'model':'sigparam2021_pi0','theta':theta.tolist(),'tau0_GeV2':tau0,'synthetic':True},indent=2)+'\n')
print('Varying-event SigParam fixture:',n,'events; injected',theta.tolist())
