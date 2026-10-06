#!/usr/bin/env python3
"""Diagnostic shared-shape/profile studies; never modify nominal estimator.

Differential tests use a conditional fixed-variance Gaussian parametric
calibration of the implemented chi-square model, including its A>=0 boundary.
They are not unconditional Poisson coverage tests (see validate_pi0_coverage).
"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import Fitter,RunEvents,write_csv,pooled_estimate,selected_estimate,TP_EDGES,DIAMOND


def read(p):return list(csv.DictReader(p.open()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('campaign','manifest','config','response','data','library','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--toys',type=int,default=500);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    cfg=json.loads(a.config.read_text())
    if cfg['configured_kinematic']!='KinC_x36_4' or cfg['phi_bins']!=12 or not np.array_equal(cfg['tprime_bin_edges'],TP_EDGES) or not np.array_equal(cfg['diamond_xb_q2_vertices'],DIAMOND):
        raise RuntimeError('Background cache selector requires configuration-specific validation; refusing a copied setting')
    accepted=[r for r in read(a.manifest) if r['accepted']=='1']
    runs=[RunEvents(a.campaign/'root'/f"signal_events_run{r['run']}.root") for r in accepted]
    q=np.array([np.float32(float(r['charge_uC'])) for r in accepted],float)
    scale=np.array([np.float32(float(r['prescale'])/(float(r['charge_uC'])/1000*float(r['livetime'])*float(r['efficiency']))) for r in accepted],float)
    k=q*scale/1000;norm=1000/q.sum()/cfg['tgt_contam'];fit=Fitter(a.library)
    mult=[np.ones(len(r.coefficient)) for r in runs];nominal,var,bad=selected_estimate(runs,mult,k,fit)
    if bad:raise RuntimeError('Nominal selected fit failed')
    params=read(a.response/'migration_parameters.csv');design=np.zeros((len(nominal),len(params)))
    for r in read(a.response/'migration_design.csv'):design[int(r['reco_row']),int(r['parameter_index'])]=float(r['response'])
    fixed=np.array([float(r['fixed_prediction']) for r in read(a.response/'migration_reco_rows.csv')])
    code,center,_=fit.solve(design,nominal*norm,var*norm**2,fixed)
    if code:raise RuntimeError('Nominal response solve failed')
    bs=np.sqrt(np.diag(np.loadtxt(a.data/'sigma_covariance.csv',delimiter=',')))
    ys,vs=np.zeros((60,200)),np.zeros((60,200))
    for run,m,weight in zip(runs,mult,k):
        for cell in np.unique(run.reco[run.selected]):
            y,v,_=run.spectra(m,run.selected&(run.reco==cell));ys[cell]+=weight*y;vs[cell]+=weight**2*v
    np.savez_compressed(a.output/'selected_mass_spectra.npz',y=ys,variance=vs)
    independent=nominal.copy();cellrows=[];invalid=[]
    for cell in range(len(nominal)):
        code,final,info=fit.fit(ys[cell],vs[cell])
        if not code:independent[cell]=final.sum()
        else:invalid.append(cell)
        cellrows.append(dict(cell=cell,valid=int(not code),classification='invalid' if code else ('zero' if info[0] else 'positive'),
            minimizer=int(info[1]),covariance=int(info[2]),background=info[6],nominal=nominal[cell],
            independent_yield=final.sum() if not code else None))
    write_csv(a.output/'independent_bin_fits.csv',cellrows)
    code,alternative,_=fit.solve(design,independent*norm,var*norm**2,fixed)
    if code:raise RuntimeError('Independent-bin diagnostic solve failed')
    variations=[('A_nominal',center,'central'),('B_valid_independent_bins',alternative,'statistical_fit_flexibility; invalid cells retain nominal explicitly')]
    x=(np.arange(200)+.5)*.002;side=((x>=.01)&(x<=.11))|((x>=.15)&(x<=.4))
    turns,widths=np.meshgrid(np.linspace(.11,.22,33),np.geomspace(.001,.10,33),indexing='ij')
    templates=1/(1+np.exp(np.clip((x[None,:]-turns.ravel()[:,None])/widths.ravel()[:,None],-700,700)))
    f=templates[:,side];profiles=[];positive_profiles=[]
    for run in runs:
        code,_,info=fit.fit(*run.spectra(np.ones(len(run.coefficient)))[:2])
        if code or info[0]:continue
        y,v,_=run.spectra(np.ones(len(run.coefficient)));vv=np.where(v[side]>0,v[side],1.)
        qg=f@(y[side]/vv);d=f*f@(1/vv);aa=np.maximum(0,qg/d);chi=np.sum(y[side]**2/vv)-np.maximum(qg,0)**2/d
        # Standard two-parameter profile contour, diagnostic only near a boundary.
        allowed=chi-chi.min()<=2.30
        positive_profiles.append(dict(run=run.run,amplitude=info[5],fitted_background=info[6],
            allowed_shape_fraction=float(allowed.mean()),turn_min=float(turns.ravel()[allowed].min()),turn_max=float(turns.ravel()[allowed].max()),
            width_min=float(widths.ravel()[allowed].min()),width_max=float(widths.ravel()[allowed].max())))
        for idx in np.flatnonzero(allowed):profiles.append(dict(run=run.run,turn=turns.ravel()[idx],width=widths.ravel()[idx],amplitude=aa[idx],delta_chi2=chi[idx]-chi.min()))
    write_csv(a.output/'positive_run_profiles.csv',positive_profiles);write_csv(a.output/'profile_shapes.csv',profiles)
    inclusive_y=ys.sum(0);inclusive_v=vs.sum(0);v=np.where(inclusive_v[side]>0,inclusive_v[side],1.)
    qg=f@(inclusive_y[side]/v);d=f*f@(1/v);amplitudes=np.maximum(0,qg/d)
    # Profiled shape changes cannot alter a certified A=0 optimum.
    if amplitudes.max()>1e-10:raise RuntimeError('Selected nominal zero certificate disagrees with diagnostic profile')
    variations.append(('C_positive_run_profile_shapes',center.copy(),'shape profile statistically allowed; selected amplitude remains zero'))
    write_csv(a.output/'parameter_variations.csv',[dict(variant=name,classification=kind,parameter_index=i,truth_block=r['truth_block'],component=r['component'],region=r['region'],
        nominal=center[i],alternative=values[i],shift=values[i]-center[i],shift_over_data_sigma=(values[i]-center[i])/bs[i])
        for name,values,kind in variations for i,r in enumerate(params)])
    # Shared theta with independent nonnegative category amplitudes vs free theta.
    tests=[];observed=[]
    for axis,ngroup in [('tprime',5),('phi',4),('Q2',2),('xB',2)]:
        gy,gv=np.zeros((ngroup,200)),np.zeros((ngroup,200))
        for run,m,weight in zip(runs,mult,k):
            e=run.events
            if axis=='tprime':group=run.reco//12
            elif axis=='phi':group=(run.reco%12)//3
            elif axis=='Q2':group=(e['Q2'].astype(np.float32)>np.mean(cfg['q2_bin_edges'])).astype(int)
            else:group=(e['xB'].astype(np.float32)>np.mean(cfg['xb_bin_edges_by_q2'][0])).astype(int)
            for g in range(ngroup):
                yy,vv,_=run.spectra(m,run.selected&(group==g));gy[g]+=weight*yy;gv[g]+=weight**2*vv
        vv=np.where(gv[:,side]>0,gv[:,side],1.);den=(1/vv)@(f*f).T
        def statistic(y):
            score=(y/vv)@f.T;improvement=np.maximum(score,0)**2/den
            shared=int(np.argmax(improvement.sum(0)))
            return float(improvement.max(1).sum()-improvement[:,shared].sum()),shared,np.maximum(score[:,shared],0)/den[:,shared]
        actual,best,amp=statistic(gy[:,side]);mean=amp[:,None]*f[best]
        rng=np.random.default_rng(9200+ngroup+ord(axis[0]));stats=[]
        for b in range(a.toys):stats.append(statistic(mean+rng.normal(size=mean.shape)*np.sqrt(vv))[0])
        pvalue=(1+np.sum(np.array(stats)>=actual))/(a.toys+1)
        tests.append(dict(axis=axis,categories=ngroup,delta_chi2=actual,calibration_replicas=a.toys,pvalue=float(pvalue),
            shared_turn=turns.ravel()[best],shared_width=widths.ravel()[best],positive_categories=int(np.sum(amp>0))))
        for g in range(ngroup):
            for mass in np.flatnonzero(side):observed.append(dict(axis=axis,category=g,mass=x[mass],residual=gy[g,mass],variance=gv[g,mass],shared_prediction=amp[g]*templates[best,mass]))
        np.savetxt(a.output/f'{axis}_null_delta_chi2.csv',stats,delimiter=',')
    write_csv(a.output/'differential_shape_tests.csv',tests);write_csv(a.output/'differential_mass_residuals.csv',observed)
    result=dict(nominal_selected_background=0,independent_valid_bins=60-len(invalid),independent_invalid_bins=invalid,
        max_published_B_shift_over_sigma=float(np.max(abs((alternative-center)/bs)[[r['region']=='published' for r in params]])),
        C_shape_variation_max_shift=0,systematic_covariance_recommended=False,
        reason='No demonstrated shape-only model shift; profile/independent-fit changes contain statistical freedom. Separate coverage failure prevents unconditional validation.',
        differential_tests=tests,positive_run_profiles=positive_profiles)
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
