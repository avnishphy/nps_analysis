#!/usr/bin/env python3
"""Event-level Poisson statistics for the validated KinC_x36_4 selected-sample estimator.

Uses production pre-weight caches and the production C++ fitter through a small
C ABI. Inputs and original analysis artifacts are read-only. One Poisson draw
per run/event is shared by every timing count and fitted-weight contribution.
Legacy mass redistribution and independent-bin fits are diagnostic alternatives.
"""
import argparse
import csv
import ctypes
import json
import os
import shutil
import hashlib
from pathlib import Path

import numpy as np
import uproot

MASS_EDGES = np.linspace(0., .4, 201)
WINDOWS = ((139, 141), (141, 143), (143, 145), (155, 157), (157, 159), (159, 161))
TP_EDGES = np.array([-.75, -.55, -.40, -.25, -.13, 0.])
DIAMOND = np.array([[.28536635126633153,3.243861779409012], [.33822268816213624,3.6193121389727954],
                    [.43953011480920373,4.698466077676237], [.3721482218990109,4.244414308310078]])


class Fitter:
    def __init__(self, library):
        self.lib = ctypes.CDLL(str(Path(library).resolve()))
        vector = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
        self.lib.nps_fit_spectrum.argtypes = [vector] * 4
        self.lib.nps_fit_spectrum.restype = ctypes.c_int
        self.lib.nps_solve_response.argtypes = [ctypes.c_int,ctypes.c_int]+[vector]*5
        self.lib.nps_solve_response.restype = ctypes.c_int

    def fit(self, y, variance):
        final, info = np.zeros(200), np.zeros(8)
        code = self.lib.nps_fit_spectrum(np.ascontiguousarray(y), np.ascontiguousarray(variance), final, info)
        return code, final, info

    def solve(self, design, y, variance, fixed=None):
        fixed=np.zeros(len(y)) if fixed is None else np.asarray(fixed,dtype=float)
        if fixed.shape != np.asarray(y).shape or np.any(~np.isfinite(fixed)):
            raise ValueError('Malformed fixed feed-in prediction')
        adjusted=np.asarray(y)-fixed
        supported=np.any(design!=0,axis=1)|(fixed!=0)
        if np.any(~supported & ((y!=0)|(variance>0))):
            return 2,None,None
        mask=supported & (variance>0)
        a=np.ascontiguousarray(design[mask])
        parameters=np.zeros(design.shape[1]); covariance=np.zeros((len(parameters),len(parameters)))
        code=self.lib.nps_solve_response(a.shape[0],a.shape[1],a.ravel(),
            np.ascontiguousarray(adjusted[mask]),np.ascontiguousarray(variance[mask]),parameters,covariance.ravel())
        return code,parameters,covariance


def timing_components(events, binned=False):
    t1, t2 = events['photon_time_1'], events['photon_time_2']
    if binned:
        # The historical normalization-error estimate integrates TH2 cells;
        # ROOT bin centers preserve its boundary convention exactly.
        t1 = 139 + (np.floor(220*(t1 - 139)/22) + .5) * .1
        t2 = 139 + (np.floor(220*(t2 - 139)/22) + .5) * .1
    inside = lambda t, lo, hi: (t > lo) & (t < hi)
    c1, c2 = inside(t1, 149, 151), inside(t2, 149, 151)
    diag = np.zeros(len(t1), bool)
    hor, ver = diag.copy(), diag.copy()
    for lo, hi in WINDOWS:
        s1, s2 = inside(t1, lo, hi), inside(t2, lo, hi)
        diag |= s1 & s2
        hor |= c1 & s2
        ver |= c2 & s1
    full1 = inside(t1, 155, 161) & inside(t2, 139, 145)
    full2 = inside(t2, 155, 161) & inside(t1, 139, 145)
    # Rows remain independent; covariance of overlapping contributions is
    # handled by using the same multiplicity, not independent histogram toys.
    return np.array([c1 & c2, diag, hor, ver, full1, full2], dtype=float)


class RunEvents:
    def __init__(self, path):
        self.run = int(path.stem.split('run')[-1])
        with uproot.open(path) as source:
            self.events = source['physics'].arrays(library='np')
            self.denominator = source['mass_denominator'].values()
        e = self.events
        if len(np.unique(e['event_id'])) != len(e['event_id']):
            raise RuntimeError(f'Run {self.run}: duplicate physical events')
        self.mass_bin = np.searchsorted(MASS_EDGES, e['mpi0_all'], side='right') - 1
        self.inmass = (self.mass_bin >= 0) & (self.mass_bin < 200)
        if not np.array_equal(np.bincount(self.mass_bin[self.inmass],minlength=200),self.denominator):
            raise RuntimeError(f'Run {self.run}: physical event cache differs from the mass denominator')
        self.components = timing_components(e)
        self.plane_components = timing_components(e, True)
        self.template_factors = np.array([0., 1/6, 1/12, 1/12, -1/18, -1/18])
        self.coefficient = self.components[0] - self.template_factors @ self.components
        self.mmiss = (e['mmiss_all'] >= .8) & (e['mmiss_all'] <= 1.1) & self.inmass
        # C++ fill_data_event takes float kinematics; match those conversions.
        q2, xb = e['Q2'].astype(np.float32).astype(float), e['xB'].astype(np.float32).astype(float)
        tp = e['t'].astype(np.float32).astype(float) - e['tmin'].astype(np.float32).astype(float)
        phi = np.mod(e['phi'], 2*np.pi)
        crosses = []
        for a, b in zip(DIAMOND, np.roll(DIAMOND, -1, axis=0)):
            crosses.append((b[0]-a[0])*(q2-a[1])-(b[1]-a[1])*(xb-a[0]))
        crosses = np.array(crosses)
        diamond = np.all(crosses >= -1e-12, axis=0) | np.all(crosses <= 1e-12, axis=0)
        it = np.searchsorted(TP_EDGES, tp, side='right') - 1
        ip = np.floor(np.where(np.isfinite(phi),phi,0)/(2*np.pi)*12).astype(int)
        self.selected = self.mmiss & diamond & np.isfinite(phi) & (q2 >= 3.3) & (q2 <= 4.7) & (xb >= .29) & (xb <= .44) & (it >= 0) & (it < 5)
        self.reco = it*12+ip

    def hist(self, weight):
        return np.bincount(self.mass_bin[self.inmass], weights=weight[self.inmass], minlength=200).astype(float)

    def spectra(self, multiplicity, selection=None):
        mult = multiplicity if selection is None else multiplicity*selection
        counts = np.array([self.hist(mult*c) for c in self.components])
        plane = self.plane_components @ mult
        background = self.template_factors @ counts
        normalization = self.template_factors @ plane
        normvar = self.template_factors**2 @ plane
        fracvar = normvar / normalization**2 if normalization > 0 else 0.
        y = counts[0] - background
        # Match the historical fitter's error model. Multiplicity counts
        # repeated observations (m), not a weighted fill whose Sumw2 is m*m.
        variance = counts[0] + self.template_factors**2 @ counts + background**2*fracvar
        denom = self.hist(mult)
        return y, variance, denom

    def estimate(self, multiplicity, fitter, selection=None):
        y, variance, denom = self.spectra(multiplicity, selection)
        code, final, info = fitter.fit(y, variance)
        p = np.divide(final, denom, out=np.zeros(200), where=denom > 0)
        event_p = np.zeros(len(multiplicity))
        event_p[self.inmass] = p[self.mass_bin[self.inmass]]
        return dict(code=code, final=final, info=info, y=y, variance=variance,
                    denom=denom, p=event_p, empty_residual=float(final[denom == 0].sum()))


def pooled_estimate(runs, multiplicities, corrections, fitter, spectrum_output=None):
    # Corrected event counts are pooled within the fixed reconstructed bins.
    # Exposure and the target factor are common scalars, applied after fitting.
    ys, variances, frozen = np.zeros((60,200)), np.zeros((60,200)), np.zeros(60)
    for run, mult, correction in zip(runs,multiplicities,corrections):
        for cell in np.unique(run.reco[run.selected]):
            select=run.selected & (run.reco==cell)
            y, variance, _=run.spectra(mult,select)
            ys[cell]+=correction*y
            variances[cell]+=correction**2*variance
            frozen[cell]+=correction**2*np.sum(mult[select]*run.coefficient[select]**2)
    if spectrum_output is not None:
        np.savez_compressed(spectrum_output,y=ys,variance=variances)
    yields, classes, backgrounds, failures = np.zeros(60), [], np.zeros(60), []
    for cell in range(60):
        code, final, info=fitter.fit(ys[cell],variances[cell])
        classes.append('zero_background' if info[0] else 'interior_valid')
        yields[cell]=final.sum()
        backgrounds[cell]=info[6]
        if code:
            failures.append(dict(cell=cell,fit_status=int(info[1]),covariance_status=int(info[2]),reason='unresolved_pooled_fit'))
    return yields, frozen, backgrounds, classes, failures


def pooling_check(campaign, output, fitter, manifest):
    quality=list(csv.DictReader(manifest.open()))
    accepted=[r for r in quality if r['accepted']=='1']
    runs=[RunEvents(campaign/'root'/f"signal_events_run{r['run']}.root") for r in accepted]
    corrections=np.array([float(r['prescale'])/(float(r['livetime'])*float(r['efficiency'])) for r in accepted])
    y,v,bg,classes,failures=pooled_estimate(runs,[np.ones(len(r.coefficient)) for r in runs],corrections,fitter,output/'pooled_spectra.npz')
    rows=[dict(reco_row=i,direct=y[i],conditional_variance=v[i],background=bg[i],classification=classes[i]) for i in range(60)]
    write_csv(output/'pooled_selected_rows.csv',rows)
    write_csv(output/'pooled_selected_failures.csv',failures)
    print(json.dumps({'pooled_yield':float(y.sum()),'pooled_background':float(bg.sum()),'pooled_fit_failures':failures}),flush=True)


def selected_estimate(runs,multiplicities,corrections,fitter):
    """Direct timing subtraction; one selected-sample shape and total BG.

    Sharing the selected-sample shape avoids fitting unidentified shapes in
    sparse per-run/reco cells. Nonnegative category amplitudes are profiled
    with their sum constrained to the inclusive background estimate. No model
    contribution is lost merely because an observed mass bin is empty.
    """
    inclusive_y,inclusive_v=np.zeros(200),np.zeros(200)
    ys,vs=np.zeros((60,200)),np.zeros((60,200))
    frozen=np.zeros(60)
    for run,mult,k in zip(runs,multiplicities,corrections):
        y,v,_=run.spectra(mult,run.selected)
        inclusive_y+=k*y;inclusive_v+=k*k*v
        for cell in np.unique(run.reco[run.selected]):
            s=run.selected & (run.reco==cell)
            y,v,_=run.spectra(mult,s)
            ys[cell]+=k*y;vs[cell]+=k*k*v
            frozen[cell]+=k*k*np.sum(mult[s]*run.coefficient[s]**2)
    code,final,info=fitter.fit(inclusive_y,inclusive_v)
    if code:
        return np.zeros(60),frozen,[dict(cell=-1,fit_status=int(info[1]),covariance_status=int(info[2]),reason='unresolved_selected_fit')]
    background=inclusive_y-final
    total=background.sum()
    allocation=np.zeros(60)
    if total>0:
        template=background/total
        x=(MASS_EDGES[:-1]+MASS_EDGES[1:])/2
        side=((x>=.01)&(x<=.11))|((x>=.15)&(x<=.4))
        variance=np.where(vs[:,side]>0,vs[:,side],1.)
        q=np.sum(ys[:,side]*template[side]/variance,axis=1)
        d=np.sum(template[side]**2/variance,axis=1)
        if np.any(d<=0):
            return np.zeros(60),frozen,[dict(cell=-1,fit_status=-1,covariance_status=-1,reason='unidentified_category_amplitudes')]
        lo,hi=np.min(q-total*d),np.max(q)
        for _ in range(100):
            mid=(lo+hi)/2
            if np.maximum(0,(q-mid)/d).sum()>total:lo=mid
            else:hi=mid
        allocation=np.maximum(0,(q-(lo+hi)/2)/d)
        allocation*=total/allocation.sum()
    return ys.sum(axis=1)-allocation,frozen,[]


def covariance_outputs(output, prefix, samples):
    if len(samples)<2:
        return None
    covariance=np.cov(np.array(samples),rowvar=False,ddof=1)
    sd=np.sqrt(np.maximum(0,np.diag(covariance)))
    correlation=np.divide(covariance,np.outer(sd,sd),out=np.zeros_like(covariance),where=np.outer(sd,sd)>0)
    np.savetxt(output/f'{prefix}_covariance.csv',covariance,delimiter=',',fmt='%.17g')
    np.savetxt(output/f'{prefix}_correlation.csv',correlation,delimiter=',',fmt='%.17g')
    return sd


def export_nominal_input(args,runs,y,variance,charge):
    if args.estimator!='selected':
        raise RuntimeError('Only the selected-sample estimator can publish a corrected input')
    destination=args.write_nominal_input
    destination.parent.mkdir(parents=True,exist_ok=True)
    if destination.exists():
        raise RuntimeError(f'Refusing to overwrite existing nominal input: {destination}')
    stage=destination.with_name(destination.name+'.tmp')
    shutil.copy2(args.combined_data,stage)
    with uproot.open(args.combined_data) as source:
        events=source['physics'].arrays(library='np')
        manifest=source['analysis_runs'].arrays(library='np')
    if set(manifest['run_number'])!={r.run for r in runs} or float(manifest['charge_uC'].astype(float).sum())!=charge:
        raise RuntimeError('Combined input and bootstrap run/exposure sets differ')
    events['pi0_weight_legacy_mass']=events['pi0_weight'].copy()
    events['source_event_id']=np.zeros(len(events['pi0_weight']),dtype=np.int64)
    for run in runs:
        select=events['run_number']==run.run
        # Event IDs are omitted by the legacy combiner; verify the complete
        # ordered mass/kinematic columns before associating cached events.
        for column in ('mpi0_all','mmiss_all','Q2','t','tmin','phi','xB','photon_time_1','photon_time_2'):
            if not np.array_equal(events[column][select],run.events[column]):
                raise RuntimeError(f'Cache/combined event identity mismatch: {run.run} {column}')
        events['pi0_weight'][select]=np.where(run.inmass,run.coefficient,0.)
        events['source_event_id'][select]=run.events['event_id']
    angles=np.linspace(0,2*np.pi,13)
    config=json.loads((Path(__file__).resolve().parents[1]/'src/xsec_extract/xsec_config/xsec_config_x36_4.json').read_text())
    vertices=np.array(config['diamond_xb_q2_vertices'])
    center=vertices.mean(axis=0)
    vertices=vertices[np.argsort(np.arctan2(vertices[:,1]-center[1],vertices[:,0]-center[0]))]
    table={'reco_row':np.arange(60,dtype=np.int32),'data':y*.584,'data_variance':variance*.584**2,
           'total_charge_uC':np.full(60,charge),'run_count':np.full(60,len(runs),dtype=np.int32),
           'tprime_lo':np.repeat(TP_EDGES[:-1],12),'tprime_hi':np.repeat(TP_EDGES[1:],12),
           'phi_lo':np.tile(angles[:-1],5),'phi_hi':np.tile(angles[1:],5),
           'q2_lo':np.full(60,3.3),'q2_hi':np.full(60,4.7),'xb_lo':np.full(60,.29),'xb_hi':np.full(60,.44),
           'mmiss_lo':np.full(60,.8),'mmiss_hi':np.full(60,1.1),'diamond':np.tile(vertices.ravel(),(60,1))}
    with uproot.update(stage) as out:
        out.mktree('physics',{k:v.dtype for k,v in events.items()}).extend(events)
        types={k:np.dtype((v.dtype,(8,))) if k=='diamond' else v.dtype for k,v in table.items()}
        out.mktree('analysis_reco_yields',types).extend(table)
        out['analysis_signal_estimator']='selected_common_shape_v1; pi0_weight=timing_event_coefficient; analysis_reco_yields=authoritative_unpolarized_signal'
        if args.covariance_from:
            summaries=list(csv.DictReader((args.covariance_from/'sigma_summary.csv').open()))
            mapping=list(csv.DictReader((args.response_dir/'migration_parameters.csv').open()))
            cov=np.loadtxt(args.covariance_from/'sigma_covariance.csv',delimiter=',')
            mc=np.loadtxt(args.mc_covariance,delimiter=',') if args.mc_covariance else np.zeros_like(cov)
            central=np.array([float(r['nominal']) for r in summaries])
            n=len(central);ii,jj=np.indices((n,n))
            if cov.shape!=(n,n) or mc.shape!=cov.shape or not np.isfinite(cov+mc).all():
                raise RuntimeError('Malformed statistical covariance')
            ctable={'i':ii.ravel().astype(np.int32),'j':jj.ravel().astype(np.int32),
                    'data_covariance':cov.ravel(),'mc_covariance':mc.ravel(),
                    'nominal_i':central[ii.ravel()],'nominal_j':central[jj.ravel()],
                    'truth_i':np.array([int(mapping[i]['truth_block']) for i in ii.ravel()],dtype=np.int32),
                    'truth_j':np.array([int(mapping[j]['truth_block']) for j in jj.ravel()],dtype=np.int32)}
            out.mktree('analysis_sigma_covariance',{k:v.dtype for k,v in ctable.items()}).extend(ctable)
            summary=json.loads((args.covariance_from/'bootstrap_summary.json').read_text())
            summary['MC_error_method']='response_cell_delta' if args.mc_covariance else 'not_included'
            out['analysis_statistical_provenance']=json.dumps(summary,sort_keys=True)
    os.replace(stage,destination)


def bootstrap(args, fitter):
    quality=list(csv.DictReader(args.manifest.open()))
    accepted=[r for r in quality if r['accepted']=='1']
    sys_path=Path(__file__).resolve().parents[1]/'src/analysis'
    import sys
    sys.path.insert(0,str(sys_path))
    from combine_analysis_branches import require_run_success
    for row in accepted:require_run_success(args.campaign/'root',int(row['run']))
    runs=[RunEvents(args.campaign/'root'/f"signal_events_run{r['run']}.root") for r in accepted]
    charges=np.array([np.float32(float(r['charge_uC'])) for r in accepted],dtype=float)
    scales=np.array([np.float32(float(r['prescale'])/(float(r['charge_uC'])/1000*float(r['livetime'])*float(r['efficiency']))) for r in accepted],dtype=float)
    normalization=1000/charges.sum()/.584
    corrections=scales*charges/1000
    design_rows=list(csv.DictReader((args.response_dir/'migration_design.csv').open()))
    npar=max(int(r['parameter_index']) for r in design_rows)+1
    design=np.zeros((60,npar))
    for r in design_rows:
        design[int(r['reco_row']),int(r['parameter_index'])]=float(r['response'])
    reco_rows=list(csv.DictReader((args.response_dir/'migration_reco_rows.csv').open()))
    fixed=np.array([float(r['fixed_prediction']) for r in reco_rows])
    samples, sigmas, differences, statuses, failures, successful_ids = [], [], [], [], [], []
    nominal_y=nominal_v=nominal_sigma=None
    for replica in range(-1,args.replicas):
        # Independent child stream per replica/run: restartable and invariant
        # to processing order; nominal exposure stays fixed.
        multiplicities=[np.ones(len(r.coefficient)) if replica<0 else
            np.random.default_rng(np.random.SeedSequence([args.seed,replica,r.run])).poisson(1,len(r.coefficient)).astype(float) for r in runs]
        legacy_y,legacy_v=np.zeros(60),np.zeros(60)
        direct_timing=np.zeros(60)
        invalid=[]
        for run,mult,correction in zip(runs,multiplicities,corrections):
            estimate=run.estimate(mult,fitter)
            if estimate['code']:
                invalid.append(dict(replica=replica,run=run.run,cell=-1,fit_status=int(estimate['info'][1]),
                                    covariance_status=int(estimate['info'][2]),reason='unresolved_run_fit'))
                continue
            s=run.selected
            legacy_y+=np.bincount(run.reco[s],weights=mult[s]*estimate['p'][s]*correction,minlength=60)
            legacy_v+=np.bincount(run.reco[s],weights=mult[s]*(estimate['p'][s]*correction)**2,minlength=60)
            direct_timing+=np.bincount(run.reco[s],weights=mult[s]*run.coefficient[s]*correction,minlength=60)
        if args.estimator=='selected':
            y,v,bad=selected_estimate(runs,multiplicities,corrections,fitter)
            invalid.extend(dict(replica=replica,run=-1,**r) for r in bad)
        elif args.estimator=='pooled':
            y,v,_,_,bad=pooled_estimate(runs,multiplicities,corrections,fitter)
            invalid.extend(dict(replica=replica,run=-1,**r) for r in bad)
        else:
            y,v=legacy_y,legacy_v
        y=y*normalization; v=v*normalization**2
        sigma_code,sigma,_=fitter.solve(design,y,v,fixed) if not invalid else (1,None,None)
        if sigma_code and not invalid:
            invalid.append(dict(replica=replica,run=-1,cell=-1,fit_status=-1,covariance_status=-1,reason='response_solve_failed'))
        failures.extend(invalid)
        statuses.append(dict(replica=replica,seed=args.seed,success=int(not invalid),failure_count=len(invalid)))
        if replica<0:
            nominal_y,nominal_v,nominal_sigma=y,v,sigma
            write_csv(args.output/'nominal_reco_rows.csv',[dict(reco_row=i,data=y[i],conditional_variance=v[i],legacy=legacy_y[i]*normalization,
                timing_direct=direct_timing[i]*normalization) for i in range(60)])
            if args.write_nominal_input and not invalid:
                export_nominal_input(args,runs,y,v,float(charges.sum()))
        elif not invalid:
            samples.append(y);sigmas.append(sigma)
            successful_ids.append(replica)
            differences.append((direct_timing-legacy_y)*normalization)
        if replica%10==0 or replica==args.replicas-1:
            print(json.dumps({'replica':replica,'successful':len(samples),'failed':sum(not s['success'] for s in statuses if s['replica']>=0)}),flush=True)
    write_csv(args.output/'replica_status.csv',statuses)
    write_csv(args.output/'replica_failures.csv',failures)
    sd=covariance_outputs(args.output,'reco',samples)
    sigma_sd=covariance_outputs(args.output,'sigma',sigmas)
    covariance_outputs(args.output,'timing_minus_legacy',differences)
    if sd is not None:
        mean=np.mean(samples,axis=0)
        write_csv(args.output/'frozen_vs_bootstrap.csv',[dict(reco_row=i,nominal=nominal_y[i],bootstrap_mean=mean[i],
            frozen_sd=np.sqrt(nominal_v[i]),bootstrap_sd=sd[i],ratio=sd[i]/np.sqrt(nominal_v[i]) if nominal_v[i]>0 else None) for i in range(60)])
    if sigma_sd is not None:
        write_csv(args.output/'sigma_summary.csv',[dict(parameter_index=i,nominal=nominal_sigma[i] if nominal_sigma is not None else None,
            bootstrap_mean=np.mean(sigmas,axis=0)[i],bootstrap_sd=sigma_sd[i]) for i in range(npar)])
    np.savez_compressed(args.output/'replicas.npz',replica=np.array(successful_ids),reco=np.array(samples),sigma=np.array(sigmas),timing_minus_legacy=np.array(differences))
    stability=[]
    for count in (50,100,250,500):
        if args.replicas>=count:
            keep=np.array(successful_ids)<count
            if keep.sum()<2:
                continue
            vs=np.std(np.array(samples)[keep],axis=0,ddof=1)
            ss=np.std(np.array(sigmas)[keep],axis=0,ddof=1)
            np.savez(args.output/f'stability_{count}.npz',reco_sd=vs,sigma_sd=ss)
            stability.append(dict(requested_replicas=count,successful_replicas=int(keep.sum()),reco_sd_norm=float(np.linalg.norm(vs)),sigma_sd_norm=float(np.linalg.norm(ss))))
    write_csv(args.output/'stability.csv',stability)
    summary=dict(estimator=args.estimator,nominal_valid=bool(statuses[0]['success']),seed=args.seed,requested=args.replicas,successful=len(samples),
                 failed=args.replicas-len(samples),failure_fraction=1-len(samples)/args.replicas,accepted_runs=[r.run for r in runs],
                 charge_uC=float(charges.sum()),efficiencies_fixed=True,livetimes_fixed=True,MC_response_fixed=True,
                 fitted_weights_recomputed=True,production_solver=True,
                 production_status='data_bootstrap_complete_fixed_response' if args.estimator=='selected' else 'diagnostic_alternative_estimator')
    summary['input_paths']={'campaign':str(args.campaign.resolve()),'response':str(args.response_dir.resolve()),'manifest':str(args.manifest.resolve())}
    identity_paths=[args.library,args.manifest,args.response_dir/'migration_design.csv',Path(__file__),
                    Path(__file__).resolve().parents[1]/'src/analysis/nps_comb_bg_pepsi.h']
    summary['sha256']={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in identity_paths}
    (args.output/'bootstrap_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def audit(campaign, output, fitter):
    rows, closure, differential = [], [], []
    for path in sorted((campaign/'root').glob('signal_events_run*.root')):
        run = RunEvents(path)
        ones = np.ones(len(run.coefficient))
        full = run.estimate(ones, fitter)
        info = full['info']
        rows.append(dict(run=run.run, classification='bad_fit' if full['code'] else 'zero_background' if info[0] else 'interior_valid',
                         amplitude=info[5], fit_status=int(info[1]), covariance_status=int(info[2]), chi2=info[3], ndf=int(info[4]),
                         background_integral=info[6], boundary_score_upper=info[7], selected_events=len(ones)))
        rootpath = campaign/'root'/f'diagnostics_run{run.run}.root'
        with uproot.open(rootpath) as f:
            h = f['h_pi0_coin_bgsub']
            nominal_p = f['physics']['pi0_weight'].array(library='np')
            closure.append(dict(run=run.run, denominator_max_difference=float(np.max(np.abs(full['denom']-run.denominator))),
                                timing_max_difference=float(np.max(np.abs(full['y']-h.values()))),
                                timing_variance_max_difference=float(np.max(np.abs(full['variance']-h.variances()))),
                                weight_max_difference=float(np.max(np.abs(full['p']-nominal_p))),
                                mass_closure_max_difference=float(np.max(np.abs(run.hist(full['p'])-full['final']))),
                                empty_mass_residual=full['empty_residual']))
        for label, selection in [('mmiss',run.mmiss), ('reco_all',run.selected)]:
            selected = run.estimate(ones, fitter, selection)
            differential.append(dict(run=run.run, selection=label, weighted=float(full['p'][selection].sum()),
                                     timing_direct=float(run.coefficient[selection].sum()), direct=float(selected['final'].sum()),
                                     direct_fit_code=selected['code'], direct_class='zero_background' if selected['info'][0] else 'interior_valid',
                                     direct_background=selected['info'][6],
                                     conditional_difference_sd=float(np.linalg.norm((run.coefficient-full['p'])[selection]))))
        for cell in range(60):
            selection=run.selected & (run.reco==cell)
            if not np.any(selection):
                continue
            selected=run.estimate(ones,fitter,selection)
            differential.append(dict(run=run.run,selection=f'reco_{cell}',weighted=float(full['p'][selection].sum()),
                                     timing_direct=float(run.coefficient[selection].sum()),direct=float(selected['final'].sum()),
                                     direct_fit_code=selected['code'],direct_class='zero_background' if selected['info'][0] else 'interior_valid',
                                     direct_background=selected['info'][6],
                                     conditional_difference_sd=float(np.linalg.norm((run.coefficient-full['p'])[selection]))))
    write_csv(output/'run_classification.csv',rows)
    write_csv(output/'mass_closure.csv',closure)
    write_csv(output/'differential_closure.csv',differential)
    print(json.dumps({'runs':len(rows), 'classes':{c:sum(r['classification']==c for r in rows) for c in ['zero_background','interior_valid','bad_fit']},
                      'maximum_mass_closure_residual':max(r['mass_closure_max_difference'] for r in closure)}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--library',type=Path,required=True)
    parser.add_argument('--pooling-check',type=Path,help='Accepted run quality manifest for pooled final-bin subtraction check')
    parser.add_argument('--manifest',type=Path)
    parser.add_argument('--response-dir',type=Path)
    parser.add_argument('--replicas',type=int,default=0)
    parser.add_argument('--seed',type=int,default=20261004)
    parser.add_argument('--estimator',choices=('legacy','pooled','selected'),default='selected')
    parser.add_argument('--combined-data',type=Path)
    parser.add_argument('--write-nominal-input',type=Path)
    parser.add_argument('--covariance-from',type=Path,help='Completed bootstrap output whose covariance is embedded in the exported input')
    parser.add_argument('--mc-covariance',type=Path,help='Optional independently estimated finite-MC sigma covariance')
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if args.write_nominal_input and args.combined_data is None:
        parser.error('--write-nominal-input requires --combined-data')
    fitter=Fitter(args.library)
    if args.replicas:
        if args.replicas<2 or args.manifest is None or args.response_dir is None:
            parser.error('Bootstrap requires replicas>=2, manifest, and response-dir')
        bootstrap(args,fitter)
    elif args.pooling_check:
        pooling_check(args.campaign,args.output,fitter,args.pooling_check)
    else:
        audit(args.campaign,args.output,fitter)


if __name__ == '__main__':
    main()
