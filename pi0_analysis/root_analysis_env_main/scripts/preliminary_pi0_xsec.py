#!/usr/bin/env python3
"""Preliminary additive-weight M0 production and ensembles.

The mass fitter and constrained M0 solver are the existing production C++
implementations. Frozen ellipse, exposure and four-row mask are explicit.
All input ROOT files and previous validation products are read-only.
"""
import argparse
import csv
import ctypes
import json
import multiprocessing as mp
import os
import time
from pathlib import Path
import numpy as np
import uproot
from bootstrap_pi0_data import Fitter, RunEvents, MASS_EDGES, timing_components
from audit_pi0_exclusive_estimator import ellipse
from validate_pi0_timing_transport import geometry

REPO = Path(__file__).resolve().parents[1]


def campaign(out):
    context = json.loads((Path(out)/'campaign_context.json').read_text())
    context = {key: Path(value) if key.endswith(('_file', '_dir')) else value
               for key, value in context.items()}
    context['config'] = json.loads((Path(out)/'config_snapshot.json').read_text())
    return context


def configured_rows(events, config):
    """Reproduce the selected JSON's one-Q2/xB-group reconstructed rows."""
    t_edges = np.asarray(config['tprime_bin_edges'], dtype=float)
    q_edges = np.asarray(config['q2_bin_edges'], dtype=float)
    xb_edges = config['xb_bin_edges_by_q2']
    phi_bins = int(config['phi_bins'])
    if len(q_edges) != 2 or len(xb_edges) != 1 or len(xb_edges[0]) != 2:
        raise ValueError('M0 toy generation currently requires one Q2/xB group')
    q, x = [events[key].astype('float32').astype(float) for key in ('Q2', 'xB')]
    tp = events['t'].astype('float32').astype(float) - events['tmin'].astype('float32').astype(float)
    it = np.searchsorted(t_edges, tp, side='right') - 1
    phi = np.mod(events['phi'], 2*np.pi)
    ip = np.floor(np.where(np.isfinite(phi), phi, 0)*phi_bins/(2*np.pi)).astype(int)
    vertices = np.asarray(config.get('diamond_xb_q2_vertices', []), dtype=float)
    inside = np.ones(len(q), dtype=bool)
    if len(vertices):
        cross = np.array([(b[0]-a[0])*(q-a[1])-(b[1]-a[1])*(x-a[0])
                          for a, b in zip(vertices, np.roll(vertices, -1, axis=0))])
        inside = np.all(cross >= -1e-12, axis=0) | np.all(cross <= 1e-12, axis=0)
    inside &= ((q >= q_edges[0]) & (q <= q_edges[1]) &
               (x >= xb_edges[0][0]) & (x <= xb_edges[0][1]) &
               (it >= 0) & (it < len(t_edges)-1) & np.isfinite(phi))
    return it*phi_bins+ip, inside

def read(p):
    with Path(p).open() as f: return list(csv.DictReader(f))

def save(p, obj):
    p = Path(p); tmp = p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(obj, indent=2)+'\n'); tmp.replace(p)

def csvout(p, records):
    p = Path(p); tmp = p.with_suffix(p.suffix+'.tmp')
    with tmp.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    tmp.replace(p)

def column(records,k): return np.array([float(r[k]) for r in records])

class Model:
    def __init__(self,out):
        self.out=out
        self.context=campaign(out);self.reference=self.context['fit_output_dir']
        self.config=self.context['config']
        self.ev=np.genfromtxt(self.reference/'model_event_cache.csv',delimiter=',',names=True)
        if self.ev.ndim == 0:self.ev=np.array([self.ev],dtype=self.ev.dtype)
        e=self.ev
        a=np.column_stack([e[k] for k in ('row','truth_block','response_weight','tau','epsilon',
            'baseline_U','baseline_LT','baseline_TT','basis_U','basis_LT','basis_TT','tprime')])
        self.lib=ctypes.CDLL(str((out/'build/libprelim.so').resolve()))
        vd=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        vi=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')
        reconstructed=read(self.reference/'model_reconstructed_yields.csv')
        self.nrows=len(reconstructed);self.mask=np.array([int(r['included']) != 0 for r in reconstructed])
        regions={int(r['truth_block']):r['region'] for r in read(self.reference/'migration_truth_blocks.csv')}
        self.physical_blocks=sorted(set(e['truth_block'][e['physical'] != 0].astype(int)))
        exterior=sorted(set(e['truth_block'][e['physical'] == 0].astype(int)))
        self.nuisance_blocks=[block for block in exterior if regions[block]=='tprime_below']
        self.fixed_blocks=[block for block in exterior if regions[block]!='tprime_below']
        if len(self.nuisance_blocks)!=1:
            raise ValueError('standard M0 requires exactly one populated fitted tprime_below feed-in block')
        self.nbins=len(self.physical_blocks);self.nparams=4+3*len(self.nuisance_blocks)
        expected=(len(self.config['tprime_bin_edges'])-1)*int(self.config['phi_bins'])
        if self.nrows != expected or self.nbins != len(self.config['tprime_bin_edges'])-1:
            raise ValueError('fit output dimensions do not match config_snapshot.json')
        self.lib.prelim_init.argtypes=[ctypes.c_int,vd,ctypes.c_double,ctypes.c_int,vi,vi,
                                      ctypes.c_int,ctypes.c_int,vi]
        self.lib.prelim_fit.argtypes=[vd,vd,vi]+[vd]*7
        self.lib.prelim_cone.argtypes=[vd,vd,ctypes.c_double,vd]
        self.pivot=float(read(self.reference/'model_context.csv')[0]['tau0_GeV2'])
        blocks=np.asarray(self.physical_blocks+self.nuisance_blocks+self.fixed_blocks,dtype=np.int32)
        roles=np.asarray([0]*len(self.physical_blocks)+[1]*len(self.nuisance_blocks)+[2]*len(self.fixed_blocks),dtype=np.int32)
        included=np.asarray(self.mask,dtype=np.int32)
        self.lib.prelim_init(len(e),np.ascontiguousarray(a.ravel()),self.pivot,
            len(blocks),blocks,roles,self.nbins,self.nrows,included)
        model_parameters=read(self.reference/'model_parameters.csv')
        physics=[r for r in model_parameters if r['name'] in ('N_U','DeltaB_U','N_LT','N_TT')]
        physics.sort(key=lambda r:('N_U','DeltaB_U','N_LT','N_TT').index(r['name']))
        migration=read(self.reference/'migration_parameters.csv')
        tprime=[r for r in migration if r['region']=='tprime_below']
        tprime.sort(key=lambda r:('U','LT','TT').index(r['component']))
        self.old=np.array([float(r['value']) for r in physics+tprime])
        if len(self.old) != self.nparams:raise ValueError('unexpected semantic model parameter inventory')
        fixed=np.isin(e['truth_block'].astype(int),self.fixed_blocks)
        fixed_baseline=np.column_stack([e['baseline_'+k][fixed] for k in ('U','LT','TT')])
        if np.any(fixed) and (not np.all(np.isfinite(fixed_baseline)) or np.any(np.all(fixed_baseline==0,axis=1))):
            raise ValueError('fixed Q2/xB feed-in lacks nominal event-level model values; regenerate the fit cache')
        self.parameter_names=['N_U','DeltaB_U','N_LT','N_TT']+[f'feedin_tprime_below_{k}' for k in ('U','LT','TT')]
        self.ones=np.ones(len(e),np.int32)

    def constrained(self,y,v,start=None,mult=None):
        """Same M0 objective with an active-physics-boundary SLSQP solver.

        The staged implementation cannot solve a physics KKT boundary. Use
        exact analytical row/constraint gradients, unchanged angular cones,
        and the same data+MC variance iteration for this numerical fallback.
        """
        from scipy.optimize import minimize, nnls
        tic=time.monotonic()
        e=self.ev;rr=e['row'].astype(int);physical=e['physical'].astype(bool)
        mult=self.ones if mult is None else np.asarray(mult)
        delta=e['tau']-self.pivot;unit=np.r_[np.ones(4),np.full(self.nparams-4,1e-6)]
        keep=physical&(mult>0);pr=rr[physical];dt=delta[physical]
        bu=e['baseline_U'][physical];bl=e['baseline_LT'][physical];bt=e['baseline_TT'][physical]
        eu=e['basis_U'][physical]*bu;el=e['basis_LT'][physical]*bl;et=e['basis_TT'][physical]*bt
        pm=mult[physical];eps=e['epsilon'][physical]
        ecoef=np.sqrt(2*eps*(1+eps));support=pm>0
        design=np.zeros((self.nrows,self.nparams))
        for j,k in [(2,'LT'),(3,'TT')]:
            design[:,j]=np.bincount(pr,weights=e['basis_'+k][physical]*e['baseline_'+k][physical]*pm,minlength=self.nrows)
        nuisance=[]
        for b,block in enumerate(self.nuisance_blocks):
            take=e['truth_block']==block
            ep=float(np.max(e['epsilon'][take&(mult>0)])) if np.any(take&(mult>0)) else 0.5
            nuisance.append((take,ep))
            for k,component in enumerate(['U','LT','TT']):
                design[:,4+3*b+k]=np.bincount(rr[take],weights=e['basis_'+component][take]*mult[take],minlength=self.nrows)
        fixed_events=np.isin(e['truth_block'].astype(int),self.fixed_blocks)
        fixed_values=np.column_stack([e['baseline_'+component] for component in ('U','LT','TT')])
        basis=np.column_stack([e['basis_'+component] for component in ('U','LT','TT')])
        fixed_prediction=np.bincount(rr[fixed_events],weights=np.sum(basis[fixed_events]*fixed_values[fixed_events],axis=1)*mult[fixed_events],minlength=self.nrows)
        variance=v.copy();cached_x=None;cached=None
        def obj(x):
            p=x*unit;jac=design.copy();terms=eu*np.exp(-p[1]*dt)*pm
            jac[:,0]=np.bincount(pr,weights=terms,minlength=self.nrows)
            mu=fixed_prediction+jac@p
            jac[:,1]=-p[0]*np.bincount(pr,weights=terms*dt,minlength=self.nrows)
            residual=(mu-y)[self.mask];return np.sum(residual**2/variance[self.mask]),2*((residual/variance[self.mask])@jac[self.mask])*unit
        def minimum(u,l,t,eps):
            bb=np.sqrt(2*eps*(1+eps))*l;aa=2*eps*t
            z=np.where(aa>0,np.clip(np.divide(-bb,2*aa,out=np.zeros_like(bb),where=aa>0),-1,1),np.where(bb>0,-1.,1.))
            return u+bb*z+eps*t*(2*z*z-1),z
        def cons(x):
            nonlocal cached_x,cached
            if cached_x is not None and np.array_equal(x,cached_x):return cached
            p=x*unit;u=p[0]*bu*np.exp(-p[1]*dt);l=p[2]*bl;t=p[3]*bt
            m,z=minimum(u,l,t,eps);m[~support]=np.inf;i=int(np.argmin(m))
            values=np.zeros(1+len(nuisance));cj=np.zeros((1+len(nuisance),self.nparams));values[0]=m[i]*1e8
            cj[0,:4]=[u[i]/p[0],-dt[i]*u[i],ecoef[i]*z[i]*bl[i],eps[i]*(2*z[i]**2-1)*bt[i]]
            cj[0]*=1e8
            for b,(_,ep) in enumerate(nuisance):
                j=4+3*b;u,l,t=p[j:j+3];m,z=minimum(np.array([u]),np.array([l]),np.array([t]),np.array([ep]))
                values[b+1]=m[0]*1e6;cj[b+1,j:j+3]=np.array([1,np.sqrt(2*ep*(1+ep))*z[0],ep*(2*z[0]**2-1)])*1e6
            cached_x=x.copy();cached=(values,cj*unit);return cached
        def event_values(p):
            values=np.zeros((len(e),3));values[physical,0]=p[0]*bu*np.exp(-p[1]*dt)
            values[physical,1]=p[2]*bl;values[physical,2]=p[3]*bt
            for b,(take,_) in enumerate(nuisance):values[take]=p[4+3*b:7+3*b]
            values[fixed_events]=fixed_values[fixed_events]
            return values
        def moments(p):
            values=event_values(p);basis=np.column_stack([e['basis_'+k] for k in ('U','LT','TT')])
            event=np.sum(basis*values,axis=1)
            return np.bincount(rr,weights=event*mult,minlength=self.nrows),np.bincount(rr,weights=event*event*mult,minlength=self.nrows)
        initial=self.old.copy() if start is None else start.copy();initial[0]*=1+1e-7;initial[4::3]+=1e-13
        # A bootstrap can empty an included row. Retain that row and initialize
        # with the existing finite-MC term; no artificial variance is added.
        variance=v+moments(initial)[1]
        if np.any(variance[self.mask]<=0):raise ValueError('Zero total variance in a fixed included row')
        p=initial;success=True;settled=False;history=[]
        initial_margin,_=minimum(p[0]*bu*np.exp(-p[1]*dt),p[2]*bl,p[3]*bt,eps)
        candidates=np.flatnonzero(support)
        physics_ids=set(candidates[np.argsort(initial_margin[candidates])[:32]].tolist())
        for iteration in range(101):
            # Whiten coordinate scales exactly as the production staged
            # solver does. Fixed raw units badly condition sparse nuisances.
            jj=design.copy();terms=eu*np.exp(-p[1]*dt)*pm
            jj[:,0]=np.bincount(pr,weights=terms,minlength=self.nrows)
            jj[:,1]=-p[0]*np.bincount(pr,weights=terms*dt,minlength=self.nrows)
            norm=np.sqrt(np.sum(jj[self.mask]**2/variance[self.mask,None],axis=0))
            unit=np.divide(1.,norm,out=np.r_[np.ones(4),np.full(self.nparams-4,1e-6)],where=norm>0)
            cached_x=None
            nextp=p.copy();stable=0;success=False
            block_matrices=[]
            for b in range(len(nuisance)):
                A=design[self.mask,4+3*b:7+3*b];H=(A.T/variance[self.mask])@A
                block_matrices.append((A,np.ascontiguousarray(H.ravel())))
            for cycle in range(400):
                priorp=nextp.copy();before=obj(nextp/unit)[0]
                # Exact production cone minimization for every nuisance block.
                for b,(_,ep) in enumerate(nuisance):
                    A,H=block_matrices[b];j=4+3*b
                    if H[0]<=0:nextp[j:j+3]=0;continue
                    terms=eu*np.exp(-nextp[1]*dt)*pm
                    mu=fixed_prediction+design@nextp+nextp[0]*np.bincount(pr,weights=terms,minlength=self.nrows)
                    lin=2*(A.T@((mu[self.mask]-y[self.mask]-A@nextp[j:j+3])/variance[self.mask]))
                    solution=np.zeros(3)
                    if self.lib.prelim_cone(H,np.ascontiguousarray(lin),ep,solution):raise ValueError('Nuisance cone solve failed')
                    nextp[j:j+3]=solution
                fixed_scaled=nextp/unit
                def full(x):return np.r_[x,fixed_scaled[4:]]
                def f4(x):
                    f,g=obj(full(x));return f,g[:4]
                x=fixed_scaled[:4]
                for exchange in range(12):
                    ids=np.array(sorted(physics_ids),dtype=int)
                    def physical_constraints(x):
                        pp=x*unit[:4];u=pp[0]*bu[ids]*np.exp(-pp[1]*dt[ids])
                        l=pp[2]*bl[ids];t=pp[3]*bt[ids];ep=eps[ids];B=ecoef[ids]
                        vals=[];jacs=[]
                        # Both endpoints must be explicit at LT=0: the
                        # derivative of min(endpoint+,endpoint-) is set-valued.
                        for z in [-1.,1.]:
                            vals.append(u+B*l*z+ep*t)
                            jacs.append(np.column_stack([u/pp[0],-dt[ids]*u,B*z*bl[ids],ep*bt[ids]]))
                        m,z=minimum(u,l,t,ep);interior=(t>0)&(abs(z)<1)
                        vals.append(np.where(interior,m,1e-8))
                        jacs.append(np.column_stack([u/pp[0],-dt[ids]*u,B*z*bl[ids],ep*(2*z*z-1)*bt[ids]])*interior[:,None])
                        return np.concatenate(vals)*1e8,np.vstack(jacs)*unit[:4]*1e8
                    fit=minimize(f4,x,jac=True,method='SLSQP',
                        bounds=[(1e-12/unit[0],None),(-20/unit[1],20/unit[1]),(None,None),(None,None)],
                        constraints=[dict(type='ineq',fun=lambda x:physical_constraints(x)[0],jac=lambda x:physical_constraints(x)[1])],
                        options=dict(ftol=1e-11,maxiter=100))
                    x=fit.x;pp=x*unit[:4]
                    mm,_=minimum(pp[0]*bu*np.exp(-pp[1]*dt),pp[2]*bl,pp[3]*bt,eps)
                    violating=np.flatnonzero(support&(mm < -1e-15))
                    new=set(violating[np.argsort(mm[violating])[:16]].tolist())-physics_ids
                    if not new:break
                    physics_ids.update(new)
                nextp[:4]=fit.x*unit[:4];c=cons(nextp/unit)[0]
                # Physical constrained stationarity in whitened coordinates.
                grad=obj(nextp/unit)[1][:4]
                residual=grad.copy()
                if c[0]<1e-5:
                    uu=nextp[0]*bu*np.exp(-nextp[1]*dt)
                    mm,zz=minimum(uu,nextp[2]*bl,nextp[3]*bt,eps)
                    active=np.flatnonzero(support&(mm<1e-12))
                    if len(active)>64:active=active[np.argsort(mm[active])[:64]]
                    allnormals=[]
                    for z in [-1.,1.]:
                        endpoint=uu[active]+ecoef[active]*nextp[2]*bl[active]*z+eps[active]*nextp[3]*bt[active]
                        aa=active[endpoint<1e-12]
                        allnormals.append(np.column_stack([uu[aa]/nextp[0],-dt[aa]*uu[aa],ecoef[aa]*z*bl[aa],eps[aa]*bt[aa]]))
                    aa=active[(nextp[3]*bt[active]>0)&(abs(zz[active])<1)]
                    allnormals.append(np.column_stack([uu[aa]/nextp[0],-dt[aa]*uu[aa],ecoef[aa]*zz[aa]*bl[aa],eps[aa]*(2*zz[aa]**2-1)*bt[aa]]))
                    normals=np.vstack(allnormals)*unit[:4]
                    if nextp[0]<=1.00001e-12:normals=np.vstack([normals,[1.,0,0,0]])
                    if nextp[1]>=20-1e-7:normals=np.vstack([normals,[0,-1.,0,0]])
                    if nextp[1]<=-20+1e-7:normals=np.vstack([normals,[0,1.,0,0]])
                    normals=normals.T
                    norms=np.linalg.norm(normals,axis=0);normals=normals[:,norms>0]/norms[norms>0]
                    if normals.shape[1]:
                        try:residual-=normals@nnls(normals,grad,maxiter=1000)[0]
                        except RuntimeError:pass
                kkt=float(np.linalg.norm(residual))
                after=obj(nextp/unit)[0];step=float(np.max(abs(nextp-priorp)/unit))
                stable=stable+1 if abs(after-before)<1e-8 and step<2e-5 and kkt<2e-4 and c[0]>-1e-6 else 0
                if stable>=2:success=True;break
            # Remove only sub-roundoff violations by increasing U. Changes are
            # below 1e-7 nb/GeV2 and are recorded through subsequent refitting.
            if c[0]<0:
                shape=bu*np.exp(-nextp[1]*dt)
                angular,_=minimum(np.zeros(len(bu)),nextp[2]*bl,nextp[3]*bt,eps)
                nextp[0]=np.nextafter(max(nextp[0],float(np.max((-angular/shape)[support]))),np.inf)
            for b in range(len(nuisance)):
                if c[b+1]<0:nextp[4+3*b]+=-c[b+1]/1e6+1e-20
            pred,mcvar=moments(nextp);nv=v+mcvar
            change=max(np.max(abs(nextp-p)/np.maximum(abs(nextp),unit*1e-3)),np.max(abs(nv[self.mask]-variance[self.mask])/nv[self.mask]))
            history.append([iteration,float(fit.fun),float(change),bool(success),cycle,kkt,float(np.min(c))])
            p=nextp
            if not success:break
            if iteration>0 and change<1e-6:settled=True;break
            variance=nv
        values=event_values(p);pub=np.zeros(3*self.nbins);means=np.zeros(self.nbins)
        for b,block in enumerate(self.physical_blocks):
            take=e['truth_block']==block;w=e['response_weight'][take]*mult[take]
            for k in range(3):pub[k*self.nbins+b]=np.average(values[take,k],weights=w)*1e9
            means[b]=np.average(e['tprime'][take],weights=w)
        terms=eu*np.exp(-p[1]*dt)*pm;jac=design.copy()
        jac[:,0]=np.bincount(pr,weights=terms,minlength=self.nrows);jac[:,1]=-p[0]*np.bincount(pr,weights=terms*dt,minlength=self.nrows)
        c=cons(p/unit)[0];info=np.array([float(obj(p/unit)[0]),iteration,float(settled),float(fit.status),c[0]/1e8,min(c[0]/1e8,np.min(c[1:])/1e6),self.pivot])
        return dict(code=0 if success and settled else 1,theta=p,pub=pub,info=info,prediction=pred,
            variance=variance,jacobian=jac,runtime=time.monotonic()-tic,history=history,solver='SLSQP exact M0 boundary fallback',means=means)

    def fit(self,y,v,start=None,mult=None):
        p=np.zeros(self.nparams);z=np.zeros(3*self.nbins);info=np.zeros(7)
        pred=np.zeros(self.nrows);var=np.zeros(self.nrows);jac=np.zeros((self.nrows,self.nparams))
        tic=time.monotonic()
        code=self.lib.prelim_fit(np.ascontiguousarray(y),np.ascontiguousarray(v),
            self.ones if mult is None else np.ascontiguousarray(mult,dtype=np.int32),
            np.ascontiguousarray(self.old if start is None else start),p,z,info,pred,var,jac.ravel())
        return dict(code=code,theta=p,pub=z,info=info,prediction=pred,variance=var,
                    jacobian=jac,runtime=time.monotonic()-tic)

    def prediction_parts(self,parameters,mult=None):
        """Independent event-level audit of physics, fitted tprime, and fixed feed-in."""
        e=self.ev;mult=self.ones if mult is None else np.asarray(mult);rows=e['row'].astype(int)
        physical=e['physical'].astype(bool);parts={name:np.zeros(self.nrows) for name in ('physics','tprime_fitted','q2_xb_fixed')}
        delta=e['tau'][physical]-self.pivot
        values=np.column_stack([parameters[0]*e['baseline_U'][physical]*np.exp(-parameters[1]*delta),
                                parameters[2]*e['baseline_LT'][physical],parameters[3]*e['baseline_TT'][physical]])
        basis=np.column_stack([e['basis_'+component] for component in ('U','LT','TT')])
        parts['physics']=np.bincount(rows[physical],weights=np.sum(basis[physical]*values,axis=1)*mult[physical],minlength=self.nrows)
        for index,block in enumerate(self.nuisance_blocks):
            take=e['truth_block'].astype(int)==block
            parts['tprime_fitted']+=np.bincount(rows[take],weights=(basis[take]@parameters[4+3*index:7+3*index])*mult[take],minlength=self.nrows)
        fixed=np.isin(e['truth_block'].astype(int),self.fixed_blocks)
        nominal=np.column_stack([e['baseline_'+component][fixed] for component in ('U','LT','TT')])
        parts['q2_xb_fixed']=np.bincount(rows[fixed],weights=np.sum(basis[fixed]*nominal,axis=1)*mult[fixed],minlength=self.nrows)
        return parts

    def ids(self):
        """Match cached generated kinematics to exclusive raw physical IDs."""
        dest=self.out/'mc_identity.npz'
        if dest.exists():return np.load(dest)['inverse']
        with uproot.open(self.context['vertex_file']) as f:
            raw=f['h10'].arrays(['Q2i','ti','phipqi'],library='np')
        with uproot.open(self.context['sim_file']) as f:
            s=f['simulation'].arrays(['event_id','is_exclusive'],library='np')
        allowed=set(s['event_id'][s['is_exclusive']!=0].astype(int))
        keys={};duplicates=0
        for i,(q,t,phi) in enumerate(zip(raw['Q2i'],raw['ti'],raw['phipqi'])):
            key=(float(q),-float(t),float(phi))
            if key in keys:duplicates+=1
            keys[key]=i
        if duplicates:raise RuntimeError('Ambiguous raw exclusive identity')
        ids=np.array([keys[(float(e['Q2']),float(e['t']),float(e['phi']))] for e in self.ev])
        assert set(ids)<=allowed
        unique,inverse=np.unique(ids,return_inverse=True)
        np.savez_compressed(dest,event_id=ids,inverse=inverse,unique=unique)
        save(self.out/'mc_identity.json',dict(response_records=len(ids),physical_events=len(unique),
            raw_exclusive_events=len(raw['ti']),all_matched_exclusive=True,duplicate_response_records=len(ids)-len(unique),
            normalization='fixed generated normalization; physical-event Poisson(1)'))
        return inverse

class Data:
    def __init__(self,out):
        context=campaign(out);config=context['config'];path=context['data_file']
        with uproot.open(path) as f:
            e=f['physics'].arrays(library='np');manifest=f['analysis_runs'].arrays(library='np')
        self.e=e; self.manifest=manifest; self.charge=float(np.sum(manifest['charge_uC'].astype(float)))
        self.geometry=geometry(path)
        row,kin=configured_rows(e,config);selected=ellipse(e['mpi0_all'],e['mmiss_all'],self.geometry)&kin
        np.testing.assert_array_equal(ellipse(e['mpi0_all'],e['mmiss_all'],self.geometry),e['is_exclusive_ellipse_combined']!=0)
        self.selected=selected;self.row=row
        self.nrows=(len(config['tprime_bin_edges'])-1)*int(config['phi_bins'])
        self.target_divisor=float(config['tgt_contam'])
        self.factor=e['scale'].astype(float)*e['charge_uC'].astype(float)/self.charge/self.target_divisor
        self.runs=[]
        for run in manifest['run_number']:
            take=e['run_number']==run;idx=np.flatnonzero(take)
            # Construct the validated timing-spectrum helper with a physical
            # identity map that also supports multiple candidates per event.
            r=RunEvents.__new__(RunEvents);r.run=int(run)
            r.events={k:v[take] for k,v in e.items()}
            r.mass_bin=np.searchsorted(MASS_EDGES,r.events['mpi0_all'],side='right')-1
            r.inmass=(r.mass_bin>=0)&(r.mass_bin<200)
            r.components=timing_components(r.events);r.plane_components=timing_components(r.events,True)
            r.template_factors=np.array([0.,1/6,1/12,1/12,-1/18,-1/18])
            r.coefficient=r.events['pi0_timing_coeff']
            unique,r.inverse=np.unique(r.events['event_id'],return_inverse=True);r.nphysical=len(unique)
            r.selected=selected[idx]&r.inmass;r.row=row[idx];r.factor=self.factor[idx]
            self.runs.append(r)
        self.fitter=Fitter(out/'build/libnps_stat.so')

    def central(self,method):
        e=self.e;w=e['pi0_weight'].copy();mean=e['pi0_timing_bin_mean'];c=e['pi0_timing_coeff']
        if method=='additive':w=c-(mean-w)
        if method=='multiplicative':
            b=mean-w
            w=np.divide(c*w,mean,out=np.full(len(c),np.nan),where=mean!=0)
            w[b==0]=c[b==0]
        if not np.all(np.isfinite(w[self.selected])):raise RuntimeError('Undefined selected '+method+' weights')
        x=w[self.selected]*self.factor[self.selected];rr=self.row[self.selected]
        return np.bincount(rr,weights=x,minlength=self.nrows),np.bincount(rr,weights=x*x,minlength=self.nrows)

    def replica(self,rng):
        y=np.zeros(self.nrows);v=np.zeros(self.nrows);failures=[];empty=0.
        for r in self.runs:
            m=rng.poisson(1,r.nphysical)[r.inverse].astype(float)
            t,var,n=r.spectra(m);code,s,info=self.fitter.fit(t,var)
            if code:failures.append(dict(run=r.run,code=int(code),fit_status=float(info[1]),covariance_status=float(info[2])))
            bn=np.divide(t-s,n,out=np.zeros(200),where=n>0)
            take=r.selected;mb=r.mass_bin[take];w=(r.coefficient[take]-bn[mb])*r.factor[take]
            y+=np.bincount(r.row[take],weights=m[take]*w,minlength=self.nrows)
            v+=np.bincount(r.row[take],weights=m[take]*w*w,minlength=self.nrows)
            empty+=float(s[n==0].sum())
        return y,v,failures,empty

def serial(result):
    return {k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in result.items()}

_REPLICA_STATE = None

def available_cpus():
    """Return CPUs granted to this process, respecting batch affinity."""
    try:return len(os.sched_getaffinity(0))
    except (AttributeError,OSError):return os.cpu_count() or 1

def run_replica(index):
    """Generate and fit one deterministic replica in a worker process."""
    model,central,start,y0,v0,data,inverse,mode,seed=_REPLICA_STATE
    rng=np.random.default_rng(np.random.SeedSequence([seed,index]))
    y,v=y0,v0;empty=0.
    if data:
        y,v,mass_fail,empty=data.replica(rng)
        if mass_fail:return index,'failure',dict(replica=index,kind='mass_fit',details=mass_fail)
    mult=None
    if inverse is not None:mult=rng.poisson(1,int(inverse.max())+1)[inverse]
    if mode=='toys':
        # Center full-estimator event-bootstrap fluctuations on folded M0.
        # This is a calibrated local M0 bootstrap pseudoexperiment, not
        # a new combinatorial-background or contamination toy model.
        y=y-y0+np.array(central['prediction'])
    try:result=model.constrained(y,v,start,mult)
    except Exception as exc:
        return index,'failure',dict(replica=index,kind='M0',details=str(exc))
    boundary_hit=bool(abs(result['theta'][1])>=20-1e-7 or result['theta'][0]<=1.00001e-12)
    if boundary_hit and result['code']:
        return index,'failure',dict(replica=index,kind='boundary',details=serial(result))
    if result['code']:
        return index,'failure',dict(replica=index,kind='M0',details=serial(result))
    if result['info'][4]<-1e-15 or abs(result['theta'][1])>20+1e-7 or result['theta'][0]<.99999e-12:
        return index,'failure',dict(replica=index,kind='boundary',details=serial(result))
    record=dict(index=index,theta=result['theta'],pub=result['pub'],y=y,v=v,
        Q=result['info'][0],iterations=result['info'][1],empty=empty,boundary_hit=boundary_hit)
    return index,'accepted',record

def central(args):
    model=Model(args.output)
    ref=read(model.reference/'model_reconstructed_yields.csv')
    check=model.fit(column(ref,'data'),column(ref,'data_sumw2'))
    save(args.output/'old_solver_parity.json',serial(check))
    assert check['code']==0
    np.testing.assert_allclose(check['theta'][:4],model.old[:4],rtol=2e-5,atol=2e-6)
    print('native parity',check['info'],check['runtime'],flush=True)
    data=Data(args.output);results={}
    for method in ('additive','old','multiplicative'):
        y,v=data.central(method);result=model.fit(y,v) if method=='old' else model.constrained(y,v)
        result.update(y=y,data_variance=v);results[method]=serial(result)
        save(args.output/'central_results.json',results)
        print(method,json.dumps({k:results[method][k] for k in ('code','theta','pub','info','runtime')}),flush=True)
        if result['code']:raise RuntimeError('Central fit failed: '+method)
    additive=results['additive'];parts=model.prediction_parts(np.asarray(additive['theta']))
    reconstructed=sum(parts.values());prediction=np.asarray(additive['prediction'])
    difference=reconstructed-prediction
    np.testing.assert_allclose(reconstructed,prediction,rtol=3e-13,atol=3e-15)
    jacobian=np.asarray(additive['jacobian'])[model.mask]/np.sqrt(np.asarray(additive['variance'])[model.mask,None])
    scale=np.linalg.norm(jacobian,axis=0);scaled=jacobian/scale;singular=np.linalg.svd(scaled,compute_uv=False)
    save(args.output/'feedin_refactor_validation.json',dict(
        parameter_order=model.parameter_names,physics_parameters=4,nuisance_parameters=3,total_parameters=model.nparams,
        fitted_nuisance_regions=['tprime_below'],fixed_feedin_blocks=list(map(int,model.fixed_blocks)),
        fixed_feedin_events=int(np.isin(model.ev['truth_block'].astype(int),model.fixed_blocks).sum()),
        prediction_equivalence_max_abs=float(np.max(abs(difference))),prediction_equivalence_max_rel=float(np.max(abs(difference)/np.maximum(abs(prediction),1e-300))),
        jacobian_rows=int(jacobian.shape[0]),jacobian_columns=int(jacobian.shape[1]),jacobian_rank=int(np.linalg.matrix_rank(scaled,tol=singular[0]*1e-10)),
        scaled_jacobian_singular_values=singular.tolist(),scaled_jacobian_condition=float(singular[0]/singular[-1]),
        fixed_prediction_by_row=parts['q2_xb_fixed'].tolist()))
    model.ids()
    save(args.output/'data_definition.json',dict(runs=list(map(int,data.manifest['run_number'])),
        exposure_uC=data.charge,events=len(data.e['event_id']),selected_events=int(data.selected.sum()),
        ellipse=data.geometry,mask=np.flatnonzero(~model.mask).tolist(),
        weight='pi0_timing_coeff - (pi0_timing_bin_mean - pi0_weight)',
        target_divisor=data.target_divisor,
        central_objective='Gaussian diagonal data Sumw2 + iterated finite-exclusive-MC Sumw2'))

def ensemble(args):
    out=args.output;dest=out/args.mode;dest.mkdir(exist_ok=args.resume)
    model=Model(out);central=json.loads((out/'central_results.json').read_text())['additive']
    start=np.array(central['theta']);y0=np.array(central['y']);v0=np.array(central['data_variance'])
    data=Data(out) if args.mode in ('data','toys') else None
    inverse=model.ids() if args.mode in ('mc','toys') else None
    cached={};elapsed=0.
    if args.resume and (dest/'replicas.npz').exists():
        prior=np.load(dest/'replicas.npz')
        for i,index in enumerate(prior['index']):
            record={k:prior[k][i] for k in prior.files}
            record['boundary_hit']=bool(abs(record['theta'][1])>=20-1e-7 or record['theta'][0]<=1.00001e-12)
            cached[int(index)]=record
        elapsed=json.loads((dest/'summary.json').read_text())['runtime_seconds']
        import shutil
        for name in ['replicas.npz','failures.json','summary.json']:
            if (dest/name).exists():shutil.copy2(dest/name,dest/(name+'.before_resume'))
    accepted=[];failures=[];attempt=0;tic=time.monotonic()
    target=min(args.replicas,args.stop_after or args.replicas)
    jobs=available_cpus() if args.jobs==0 else args.jobs
    jobs=max(1,min(jobs,args.replicas*3))
    global _REPLICA_STATE
    _REPLICA_STATE=(model,central,start,y0,v0,data,inverse,args.mode,args.seed)
    missing=(index for index in range(args.replicas*3) if index not in cached)
    pool=None
    if jobs>1:
        # Linux fork shares the large read-only event arrays copy-on-write.
        # Each process owns its fitter state; BLAS/OpenMP stay at one thread.
        pool=mp.get_context('fork').Pool(jobs)
        results=pool.imap(run_replica,missing,chunksize=1)
    else:
        results=map(run_replica,missing)
    try:
        for index in range(args.replicas*3):
            attempt=index+1
            if index in cached:
                accepted.append(cached[index])
                if len(accepted)>=target:break
                continue
            save(dest/'progress.json',dict(attempted=attempt,accepted=len(accepted),replica=index,jobs=jobs))
            replica_index,status,record=next(results)
            if replica_index != index:raise RuntimeError('Replica results arrived out of order')
            if status=='accepted':accepted.append(record)
            else:
                failures.append(record);save(dest/'failures.json',failures)
            if status=='accepted' and (len(accepted)%20==0 or len(accepted)==target):
                np.savez_compressed(dest/'replicas.npz',**{k:np.array([r[k] for r in accepted]) for k in accepted[0]})
                summary=dict(requested=args.replicas,attempted=attempt,accepted=len(accepted),
                    fit_failures=len(failures),mass_fit_failures=sum(r['kind']=='mass_fit' for r in failures),
                    M0_failures=sum(r['kind']=='M0' for r in failures),boundary_failures=sum(r['kind']=='boundary' for r in failures),
                    coordinate_boundary_hits=sum(bool(r['boundary_hit']) for r in accepted),
                    seed=args.seed,jobs=jobs,available_cpus=available_cpus(),runtime_seconds=elapsed+time.monotonic()-tic,
                    critical_stop_reason=args.stop_reason if target<args.replicas and len(accepted)==target else None)
                save(dest/'summary.json',summary);save(dest/'failures.json',failures)
                print(json.dumps(summary),flush=True)
            if len(accepted)>=target:break
    finally:
        if pool is not None:
            pool.terminate();pool.join()
    if len(accepted)<target:raise RuntimeError('Accepted target not reached; inspect all recorded failures')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['central','data','mc','toys'])
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--replicas',type=int,default=1000);p.add_argument('--seed',type=int,default=20261005)
    p.add_argument('--jobs',type=int,default=1,help='Replica worker processes; 0 uses every CPU in the affinity mask')
    p.add_argument('--resume',action='store_true',help='Reuse checkpointed fits; replay missing or previously failed IDs')
    p.add_argument('--stop-after',type=int,help='Shorten only for an explicitly recorded critical blocker')
    p.add_argument('--stop-reason',default='')
    args=p.parse_args()
    if args.jobs<0:p.error('--jobs must be zero or a positive integer')
    if args.stop_after and args.stop_after<args.replicas and (args.stop_after<300 or not args.stop_reason):p.error('Critical early stop requires >=300 accepted and an explicit reason')
    (central if args.mode=='central' else ensemble)(args)

if __name__=='__main__':
    # ROOT dictionaries loaded through two ctypes bridges can crash during
    # interpreter teardown on this installation. All products are explicitly
    # closed before this exit; exceptions remain nonzero and retain traceback.
    import sys,traceback
    try:main();code=0
    except BaseException:traceback.print_exc();code=1
    sys.stdout.flush();sys.stderr.flush();os._exit(code)
