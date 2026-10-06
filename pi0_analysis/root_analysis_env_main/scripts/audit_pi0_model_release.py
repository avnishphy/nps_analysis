#!/usr/bin/env python3
"""Diagnostic release audit of the current ellipse model fit; not an error estimator.

Reads production model exports and the actual combined events. Does not fit a
different cross-section estimator or transfer a covariance between selections.
The conditional diagonal metric is explicitly diagnostic throughout.
"""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'scripts'))
from bootstrap_pi0_data import timing_components


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def csvout(path, rows, fields=None):
    stage = path.with_suffix(path.suffix + '.tmp')
    with stage.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    stage.replace(path)


def jsonout(path, obj):
    stage = path.with_suffix('.json.tmp')
    stage.write_text(json.dumps(obj, indent=2) + '\n')
    stage.replace(path)


def floats(rows, key):
    return np.array([float(r[key]) for r in rows])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nominal', type=Path, required=True)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    nominal = args.nominal
    config = json.loads((nominal / 'pipeline_config.json').read_text())
    if config['configured_kinematic'] != 'KinC_x36_4':
        raise ValueError('This audited selection contract is specific to KinC_x36_4')
    pars = read(nominal / 'model_parameters.csv')
    names = [r['name'] for r in pars]
    p = floats(pars, 'value')
    context = read(nominal / 'model_context.csv')[0]
    status = read(nominal / 'model_fit_status.csv')[0]
    rows = read(nominal / 'model_reconstructed_yields.csv')
    blocks = read(nominal / 'migration_truth_blocks.csv')
    ev = np.genfromtxt(nominal / 'model_event_cache.csv', delimiter=',', names=True)
    nrow, npar = len(rows), len(pars)
    included = np.array([r['included'] == '1' for r in rows])
    y, var = floats(rows, 'data'), floats(rows, 'data_sumw2')
    objvar = floats(rows, 'objective_variance')
    pred = floats(rows, 'prediction')
    jac = np.zeros((nrow, npar))
    for r in read(nominal / 'model_row_jacobian.csv'):
        jac[int(r['row']), int(r['parameter'])] = float(r['derivative'])
    # Independent reconstruction of the event-integrated model and its Jacobian.
    tau0 = float(context['tau0_GeV2'])
    erow, eb = ev['row'].astype(int), ev['truth_block'].astype(int)
    physical = ev['physical'].astype(bool)
    basis = np.array([ev['basis_' + k] for k in ('U', 'LT', 'TT')]).T
    d = ev['tau'] - tau0
    values = np.zeros_like(basis)
    derivatives = np.zeros((len(ev), 3, 4))
    values[physical, 0] = p[0] * np.exp(-p[1]*d[physical]) * ev['baseline_U'][physical]
    values[physical, 1] = p[2] * ev['baseline_LT'][physical]
    values[physical, 2] = p[3] * ev['baseline_TT'][physical]
    derivatives[physical, 0, 0] = values[physical, 0] / p[0]
    derivatives[physical, 0, 1] = -d[physical] * values[physical, 0]
    derivatives[physical, 1, 2] = ev['baseline_LT'][physical]
    derivatives[physical, 2, 3] = ev['baseline_TT'][physical]
    event_jac = np.zeros((len(ev), npar))
    event_jac[:, :4] = np.einsum('ek,ekp->ep', basis, derivatives)
    for j, name in enumerate(names[4:], 4):
        _, b, component = name.split('_')
        k = ('U', 'LT', 'TT').index(component)
        select = eb == int(b)
        values[select, k] = p[j]
        event_jac[select, j] = basis[select, k]
    prediction_check = np.bincount(erow, weights=np.sum(basis*values, axis=1), minlength=nrow)
    jac_check = np.array([np.bincount(erow, weights=event_jac[:, j], minlength=nrow) for j in range(npar)]).T
    np.testing.assert_allclose(prediction_check, pred, rtol=1e-11, atol=1e-14)
    np.testing.assert_allclose(jac_check, jac, rtol=1e-11, atol=1e-12)
    # Central response-weighted observables. Component-major ordering is fixed.
    published = read(nominal / 'model_structure_functions.csv')
    nb = len(published)
    z = np.zeros(3*nb)
    zj = np.zeros((3*nb, npar))
    central_table = []
    for i, r in enumerate(published):
        b = int(r['truth_block'])
        take = eb == b
        w = ev['response_weight'][take]
        entry = dict(bin=b, tprime_lo=config['tprime_bin_edges'][b],
                     tprime_hi=config['tprime_bin_edges'][b+1], tprime_mean=float(r['tprime']))
        for k, component in enumerate(('U', 'LT', 'TT')):
            z[k*nb+i] = np.average(values[take, k], weights=w)
            zj[k*nb+i, :4] = np.average(derivatives[take, k], weights=w, axis=0)
            np.testing.assert_allclose(z[k*nb+i], float(r['sigma_'+component]), rtol=1e-11)
            entry.update({f'sigma_{component}': z[k*nb+i]*1e9,
                          f'{component}_data_stat': 'nan', f'{component}_MC_stat': 'nan',
                          f'{component}_total_stat': 'nan'})
        entry.update(unit='nb/GeV2', model='M0', status='DIAGNOSTIC_ONLY_NOT_PUBLISHED')
        central_table.append(entry)
    csvout(out / 'central_values_NOT_FOR_PUBLICATION.csv', central_table)
    jsonout(out / 'ordering.json', dict(theta=names, z_pub=[f'{k}(bin{i})' for k in ('U','LT','TT') for i in range(nb)],
             internal_units='microbarn/MeV2', display_units='nb/GeV2', display_scale=1e9,
             tprime_definition='t - tmin (negative); tau = -tprime',
             weighting='sum(full_weight/sigcm * fitted_component)/sum(full_weight/sigcm), accepted generated events',
             uncertainties='unavailable: complete ellipse data-statistical estimator not validated'))
    # This is NOT V_data: frozen-weight Sumw2 is only a conditional diagnostic.
    spectrum = np.sort(var[included])[::-1]
    csvout(out / 'conditional_data_eigenvalues.csv', [dict(mode=i, eigenvalue=v) for i,v in enumerate(spectrum)])
    a = jac[included] / np.sqrt(var[included, None])
    units = 1 / np.linalg.norm(a, axis=0)
    u, s, vt = np.linalg.svd(a*units, full_matrices=True)
    rank = int(np.sum(s > 1e-10*s[0]))
    nuisance_rank = int(np.linalg.matrix_rank(a[:,4:]*units[4:], tol=1e-10))
    un, sn, _ = np.linalg.svd(a[:,4:]*units[4:], full_matrices=False)
    physics_residual = a[:,:4]*units[:4] - un[:,:nuisance_rank] @ (un[:,:nuisance_rank].T @ (a[:,:4]*units[:4]))
    physics_rank = int(np.linalg.matrix_rank(physics_residual, tol=1e-10))
    null = vt[rank:]
    null_effect = zj @ (units[:,None]*null.T)
    linear_operator = units[:,None] * (vt[:rank].T/s[:rank]) @ u[:,:rank].T
    local_cov = linear_operator @ linear_operator.T
    local_pub_cov = zj @ local_cov @ zj.T
    pub_scale = np.sqrt(np.diag(local_pub_cov))
    rank_summary = dict(metric='conditional fixed-weight data Sumw2; NOT authoritative data covariance',
        coordinate_scaling='theta = theta_nominal + column_norm_inverse * dimensionless_coordinate',
        relative_svd_tolerance=1e-10, conditional_covariance_rank=int(np.sum(spectrum>1e-10*spectrum[0])),
        authoritative_covariance_rank=None, detector_jacobian_rank=rank, nuisance_rank=nuisance_rank,
        profiled_physics_rank=physics_rank, detector_null_dimension=len(null),
        maximum_published_null_effect=float(np.max(np.abs(null_effect), initial=0)),
        published_local_dimension=int(np.linalg.matrix_rank(zj[:,:4]/np.maximum(abs(z[:,None]),1e-30), tol=1e-10)),
        descriptive_unconstrained_dof=int(included.sum()-rank),
        boundary_caveat='The fitted low-tprime nuisance angular cone is active; unconstrained ranks/DOF do not validate coverage.')
    jsonout(out/'rank_summary.json', rank_summary)
    csvout(out/'singular_values.csv', [dict(mode=i, singular_value=v, relative=v/s[0]) for i,v in enumerate(s)])
    for name, mat in [('identifiable_basis',vt[:rank]), ('null_basis',null)]:
        csvout(out/f'{name}.csv', [dict(mode=i, **dict(zip(names,v))) for i,v in enumerate(mat)], ['mode']+names)
    csvout(out/'coordinate_scales.csv', [dict(parameter=name, scale=units[j]) for j,name in enumerate(names)])
    csvout(out/'rank_tolerance_sweep.csv', [dict(relative_tolerance=t,
        conditional_covariance_rank=int(np.sum(spectrum>t*spectrum[0])), jacobian_rank=int(np.sum(s>t*s[0])),
        central_values_refitted=0, note='algebraic diagnostic only; no validated V_data for fit tolerance test') for t in (1e-8,1e-10,1e-12)])
    nuisance = []
    fitted_feedin=next(r for r in blocks if r.get('fit_treatment')=='fitted_tprime_feedin')
    for j,name in enumerate(names[4:],4):
        b = int(fitted_feedin['truth_block']);br=fitted_feedin
        coupling = np.abs(zj @ local_cov[:,j]) / (pub_scale*np.sqrt(local_cov[j,j]))
        maxima = [float(np.max(coupling[k*nb:(k+1)*nb])) for k in range(3)]
        null_projection = float(np.linalg.norm(null[:,j]))
        category = 'A' if max(maxima)>=.1 else ('C' if null_projection>.99 else 'B')
        nuisance.append(dict(nuisance=name, truth_block=b, region=br['region'], events=br['events'],
            response_support=int(np.count_nonzero(jac[included,j])), response_norm=float(np.linalg.norm(jac[included,j])),
            rank_projection=1-null_projection**2, coupling_U=maxima[0], coupling_LT=maxima[1], coupling_TT=maxima[2],
            classification=category, action='retain physical feed-in; no coordinate confidence errors',
            caveat='conditional unconstrained correlation screening; A threshold 0.1; B estimable residual coupling'))
    csvout(out/'nuisance_audit.csv', nuisance)
    angular=[]
    for br in blocks:
        b=int(br['truth_block'])
        if br.get('fit_treatment')!='fitted_tprime_feedin':continue
        eps=float(br['epsilon_max'])
        cu,cl,ct=[p[names.index(f'feedin_tprime_below_{k}')] for k in ('U','LT','TT')]
        linear=np.sqrt(2*eps*(1+eps))*cl;quadratic=eps*ct
        points=[-1.,1.]
        if quadratic>0 and abs(linear/(4*quadratic))<1:points.append(-linear/(4*quadratic))
        margin=min(cu+linear*x+quadratic*(2*x*x-1) for x in points)
        scale=max(abs(cu),abs(cl),abs(ct))
        angular.append(dict(truth_block=b,epsilon=eps,minimum=margin,relative_margin=margin/scale,
                            active=int(abs(margin/scale)<1e-10)))
    csvout(out/'nuisance_angular_constraints.csv',angular)
    csvout(out/'response_support.csv', blocks)
    excluded = []
    for r in rows:
        if r['included']=='1': continue
        index = int(r['row'])
        support = int(np.sum(erow==index))
        excluded.append(dict(row=index,tprime_slice=r['it'],phi_bin=r['ip'],reason='zero_observed_sumw2',
            response_events=support,response_norm=float(np.linalg.norm(jac[index])),predicted_yield=r['prediction'],
            predetermined_from_MC=0,angular_hole=1,release_gate='FAIL_data_dependent_exclusion'))
    csvout(out/'excluded_rows.csv', excluded)
    # Independent run/exposure and yield reconstruction from the actual input.
    with uproot.open(args.data) as source:
        e=source['physics'].arrays(library='np')
        manifest=source['analysis_runs'].arrays(library='np')
        input_stat_keys=[k for k in source.keys() if any(v in k for v in ('covariance','reco_yields','statistical'))]
    charge=float(manifest['charge_uC'].astype(float).sum())
    assert set(manifest['run_number'])==set(e['run_number'])
    runrows=[]
    for i,run in enumerate(manifest['run_number']):
        run=int(run);take=e['run_number']==run
        runstatus=read(args.data.parent/f'analysis_status_run{run}.csv')[0]
        assert int(manifest['success'][i])==1 and runstatus['success']=='1' and runstatus['fit_valid']=='1'
        assert np.all(e['charge_uC'][take]==manifest['charge_uC'][i])
        assert np.all(e['scale'][take]==manifest['scale'][i])
        runrows.append(dict(run=run,charge_uC=float(manifest['charge_uC'][i]),scale=float(manifest['scale'][i]),
                            events=int(take.sum()),classification=runstatus['classification'],accepted=1))
    csvout(out/'accepted_runs.csv',runrows)
    q=e['Q2'].astype('float32').astype(float); x=e['xB'].astype('float32').astype(float)
    vertices=np.array(config['diamond_xb_q2_vertices']);center=vertices.mean(axis=0)
    vertices=vertices[np.argsort(np.arctan2(vertices[:,1]-center[1],vertices[:,0]-center[0]))]
    cross=np.array([(b[0]-a[0])*(q-a[1])-(b[1]-a[1])*(x-a[0]) for a,b in zip(vertices,np.roll(vertices,-1,axis=0))])
    diamond=np.all(cross>=-1e-12,axis=0)|np.all(cross<=1e-12,axis=0)
    tp=e['t'].astype('float32').astype(float)-e['tmin'].astype('float32').astype(float)
    it=np.searchsorted(config['tprime_bin_edges'],tp,side='right')-1
    phi=np.mod(e['phi'],2*np.pi);ip=np.floor(phi/(2*np.pi)*12).astype(int)
    selected=(e['is_exclusive_ellipse_combined']!=0)&diamond&(q>=3.3)&(q<=4.7)&(x>=.29)&(x<=.44)&(it>=0)&(it<5)
    row=12*it+ip
    factor=e['scale'].astype('float32').astype(float)*e['charge_uC'].astype('float32').astype(float)/charge/config['tgt_contam']
    weights=e['pi0_weight']*factor
    reco=np.bincount(row[selected],weights=weights[selected],minlength=nrow)
    reco_var=np.bincount(row[selected],weights=weights[selected]**2,minlength=nrow)
    np.testing.assert_allclose(reco,y,rtol=2e-12,atol=1e-14)
    np.testing.assert_allclose(reco_var,var,rtol=2e-12,atol=1e-14)
    coefficient=np.array([1.,-1/6,-1/12,-1/12,1/18,1/18])@timing_components(e)
    coefficient*=((e['mpi0_all']>=0)&(e['mpi0_all']<.4))
    zero_runs=[r['run'] for r in runrows if r['classification']=='zero_background']
    zero=np.isin(e['run_number'],zero_runs)
    closure=[]
    for r in range(nrow):
        take=selected&(row==r);zr=take&zero
        closure.append(dict(row=r,nominal=reco[r],independent_timing_no_comb_bg=float(np.sum(coefficient[take]*factor[take])),
            zero_background_runs_nominal=float(np.sum(weights[zr])),
            zero_background_runs_timing=float(np.sum(coefficient[zr]*factor[zr]))))
    csvout(out/'selected_weight_transport.csv',closure)
    zero_legacy=float(np.sum(weights[selected&zero]))
    zero_timing=float(np.sum(coefficient[selected&zero]*factor[selected&zero]))
    # A fixed-geometry, zero-background derivative is sufficient to disprove
    # equivalence of frozen-weight Sumw2 and the event-statistical covariance.
    # It is NOT the full covariance: fitted geometry/background remain missing.
    mb=np.searchsorted(np.linspace(0,.4,201),e['mpi0_all'],side='right')-1
    influence_cov=np.zeros((nrow,nrow));frozen_zero=np.zeros(nrow)
    mass_mean_error=0.;inclusive_error=0.;derivative_error=0.
    for run in zero_runs:
        take=(e['run_number']==run)&(mb>=0)&(mb<200)
        idx=np.where(take)[0];m=mb[idx];c=coefficient[idx];sel=selected[idx];rr=row[idx]
        den=np.bincount(m,minlength=200)
        mass_signal=np.bincount(m,weights=c,minlength=200)
        mean=np.divide(mass_signal,den,out=np.zeros(200),where=den>0)
        counts=np.zeros((200,nrow))
        np.add.at(counts,(m[sel],rr[sel]),1)
        fraction=np.divide(counts,den[:,None],out=np.zeros_like(counts),where=den[:,None]>0)
        influence=(c-mean[m])[:,None]*fraction[m]
        influence[np.where(sel)[0],rr[sel]]+=mean[m[sel]]
        influence*=factor[idx, None]
        influence_cov+=influence.T@influence
        frozen_zero+=np.bincount(rr[sel],weights=weights[idx[sel]]**2,minlength=nrow)
        mass_mean_error=max(mass_mean_error,float(np.max(abs(e['pi0_weight'][idx]-mean[m]))))
        inclusive_error=max(inclusive_error,abs(float(np.sum(e['pi0_weight'][idx])-np.sum(c))))
        # Numerical derivative with one event multiplicity, keeping this
        # explicitly conditional diagnostic estimator fixed.
        j=int(np.flatnonzero(sel)[0]) if np.any(sel) else 0
        def perturbed(step):
            dd=den.astype(float).copy();ss=mass_signal.copy();cc=counts.copy()
            dd[m[j]]+=step;ss[m[j]]+=step*c[j]
            if sel[j]:cc[m[j],rr[j]]+=step
            mu=np.divide(ss,dd,out=np.zeros(200),where=dd>0)
            return (mu@cc)*factor[idx[j]]
        num=(perturbed(1e-4)-perturbed(-1e-4))/2e-4
        derivative_error=max(derivative_error,float(np.max(abs(num-influence[j]))))
    assert mass_mean_error<1e-12 and derivative_error<1e-9
    np.savetxt(out/'diagnostic_zero_background_fixed_geometry_covariance.csv',influence_cov,delimiter=',')
    sd=np.sqrt(influence_cov.diagonal())
    corr=np.divide(influence_cov,sd[:,None]*sd[None,:],out=np.zeros_like(influence_cov),where=sd[:,None]*sd[None,:]>0)
    np.fill_diagonal(corr,0)
    conditional_covariance_check=dict(zero_background_runs=len(zero_runs),
        maximum_mass_mean_weight_error=mass_mean_error,maximum_inclusive_closure_error=inclusive_error,
        maximum_finite_difference_influence_error=derivative_error,
        maximum_offdiagonal_correlation=float(np.max(abs(corr))),
        summed_variance_recomputed_mass_weights=float(influence_cov.sum()),
        summed_variance_frozen_weights=float(frozen_zero.sum()),
        scope='Diagnostic linear influence, 54 zero-background runs, fixed geometry and background branch. NOT V_data.')
    jsonout(out/'weight_covariance_non_equivalence.json',conditional_covariance_check)
    exposure=dict(accepted_runs=[r['run'] for r in runrows],run_count=len(runrows),events=len(weights),
        selected_events=int(selected.sum()),charge_uC=charge,charge_mC=charge/1000,target_divisor=config['tgt_contam'],
        omitted_config_run=6569,omission_reason='waveform input unavailable; not in accepted input manifest',
        input_statistical_objects=input_stat_keys,nominal_yield_sum=float(reco.sum()),
        timing_yield_without_comb_background=float(np.sum(coefficient[selected]*factor[selected])),
        zero_background_runs=len(zero_runs),zero_background_nominal=zero_legacy,zero_background_timing=zero_timing,
        zero_background_transport_difference=zero_timing-zero_legacy,
        weight_covariance_diagnostic=conditional_covariance_check,
        max_yield_reconstruction_error=float(np.max(abs(reco-y))),
        caveat='Timing sums are a transport diagnostic, not a replacement signal estimator or measured cross section.')
    jsonout(out/'input_audit.json',exposure)
    # Save every start and every completed outer solve, including published values.
    def observable(theta):
        ans=[]
        for k in range(3):
            for r in published:
                take=eb==int(r['truth_block']);w=ev['response_weight'][take]
                value=(theta[0]*np.exp(-theta[1]*d[take])*ev['baseline_U'][take] if k==0 else
                       theta[k+1]*ev['baseline_'+('LT' if k==1 else 'TT')][take])
                ans.append(np.average(value,weights=w)*1e9)
        return np.array(ans)
    def prediction(theta):
        terms=np.zeros(len(ev))
        terms[physical]=(basis[physical,0]*theta[0]*np.exp(-theta[1]*d[physical])*ev['baseline_U'][physical]
            +basis[physical,1]*theta[2]*ev['baseline_LT'][physical]
            +basis[physical,2]*theta[3]*ev['baseline_TT'][physical])
        for j,name in enumerate(names[4:],4):
            b=int(name.split('_')[1]);k=('U','LT','TT').index(name.split('_')[2]);take=eb==b
            terms[take]+=basis[take,k]*theta[j]
        return np.bincount(erow,weights=terms,minlength=nrow)
    starts=read(nominal/'model_starts.csv')
    start_z=np.array([observable(np.array([float(r[n]) for n in names])) for r in starts])
    start_spread=float(np.max(np.ptp(start_z,axis=0)))
    csvout(out/'starts.csv',[dict(start=r['start'],objective=r['objective'],mc_iterations=r['mc_iterations'],
        **{f'z{j}':value for j,value in enumerate(start_z[i])}) for i,r in enumerate(starts)])
    last={}
    with (nominal/'staged_solver_history.csv').open() as stream:
        for r in csv.DictReader(stream): last[int(r['solve'])]=r
    history=[]
    prev=None
    for solve,r in sorted(last.items()):
        theta=np.array([float(r[n]) for n in names]);obs=observable(theta)
        record=dict(solve=solve,objective=r['objective_after'],KKT_relative=r['KKT_relative'],
            parameter_scaled_change='nan' if prev is None else float(np.linalg.norm((theta-prev)/units)),
            prediction_change_norm='nan' if prev is None else float(np.linalg.norm(prediction(theta)-prediction(prev))),
            **dict(zip(names,theta)),**{f'z{j}':v for j,v in enumerate(obs)})
        history.append(record);prev=theta
    csvout(out/'iteration_history.csv',history)
    jsonout(out/'audit_summary.json',dict(objective=float(status['objective']),rank=rank_summary,
        input=exposure,maximum_start_spread_nb_GeV2=start_spread,
        max_response_prediction_error=float(np.max(abs(prediction_check-pred))),
        data_replicas_requested=0,data_replicas_accepted=0,mc_replicas_requested=0,mc_replicas_accepted=0,toys=0,
        chosen_production_model=None,candidate='M0',verdict='MODEL EXTRACTION NOT READY',
        ensemble_stop_reason='Upstream ellipse signal estimator and data-derived selection lack a validated statistical contract; mass-weight transport does not close.'))
    gates=[
        ('production run/exposure consistency','PASS','56 accepted runs, input manifest and event charge/scale agree'),
        ('response construction','PASS','independent event prediction and Jacobian reproduce model exports'),
        ('reconstructed-data covariance rank','FAIL','complete covariance for the ellipse estimator unavailable'),
        ('published-parameter estimability','FAIL','conditional local rank 4; production metric and ensemble validation unavailable'),
        ('nuisance null-space treatment','PASS','no detector-null coordinate found; all feed-in retained'),
        ('chosen-model toy bias','FAIL','not run: upstream statistical estimator blocks model validation'),
        ('chosen-model pull widths','FAIL','not run: upstream statistical estimator blocks model validation'),
        ('chosen-model coverage','FAIL','not run: upstream statistical estimator blocks model validation'),
        ('data-bootstrap covariance convergence','FAIL','zero production replicas; selection contract mismatch'),
        ('finite-MC covariance convergence','FAIL','zero production replicas; data estimator unresolved'),
        ('data+MC covariance-addition cross-check','FAIL','not run; component production ensembles unavailable'),
        ('numerical iteration stability','PASS','six starts converge in 15 variance updates; identical central objective'),
        ('initialization stability','FAIL','six starts agree numerically; comparison in validated statistical sigma unavailable'),
        ('detector-level goodness of fit','FAIL','toy-calibrated GOF unavailable'),
        ('influence/leverage stability','FAIL','deletion refits deferred; valid total-stat sigma unavailable'),
        ('final statistical error bars available','FAIL','no validated ellipse data/MC ensembles'),
        ('published covariance saved','FAIL','no authoritative covariance; conditional diagnostics distinctly named'),
        ('report metadata consistency','PASS','free slope, fixed pivot, staged solver, and conditional ranks explicit'),
        ('row exclusion independent of data','FAIL','rows 0,5,6,11 supported by MC, omitted for observed zero variance'),
        ('selected-yield statistical estimator','FAIL','mass-only redistribution transport failure and missing geometry/background propagation'),
    ]
    csvout(out/'release_gates.csv',[dict(gate=g,status=s,evidence=e) for g,s,e in gates])
    csvout(out/'model_comparison.csv',[
        dict(model=m,physics_parameters=n,estimable_physics_rank=physics_rank if m=='M0' else 'not_tested',
             Q_observed=status['objective'] if m=='M0' else 'not_tested',toy_GOF='not_tested',
             DeltaQ='not_tested',toy_DeltaQ_probability='not_tested',bias='not_tested',pull_width='not_tested',
             coverage='not_tested',variance_inflation='not_tested',MC_stability='not_tested',
             verdict='diagnostic_candidate_only' if m=='M0' else 'not_considered_upstream_gate_failed')
        for m,n in [('M0',4),('M1',5),('M2',5),('M3',6)]])
    # Model-only diagnostic PDF. No confidence bands or final-physics claims.
    with PdfPages(out/'model_release_audit.pdf') as pdf:
        def save(fig,name):
            fig.tight_layout();pdf.savefig(fig);fig.savefig(out/(name+'.png'),dpi=130)
            if name.startswith('diagnostic_sigma_'):fig.savefig(out/(name+'.pdf'))
            plt.close(fig)
        fig,ax=plt.subplots(figsize=(11,8));ax.axis('off')
        ax.text(.02,.97,'MODEL EXTRACTION NOT READY\n\n'
            f'M0 diagnostic central reproduction: Q={float(status["objective"]):.10f}\n'
            f'{included.sum()} included rows, {npar} coordinates, 15 finite-MC variance updates\n'
            f'U slope free: DeltaB_U={p[1]:.9g} GeV^-2; fixed tau0={tau0:.9g} GeV2\n'
            f'Exposure: {len(runrows)} runs, {charge:.10f} uC\n'
            f'Conditional local ranks: detector {rank}, nuisance {nuisance_rank}, physics {physics_rank}\n\n'
            'Blocker: current ellipse mass-weight estimator lacks validated complete data statistics.\n'
            'Available physical-event bootstrap uses a different fixed selection and signal estimator.\n'
            'Fitted mass weights and data-derived ellipse geometry cannot silently be frozen.\n'
            f'54 zero-background runs: selected mass-weight transport difference {zero_timing-zero_legacy:.6g} /mC.\n'
            'Four MC-supported rows are excluded based on observed zero variance.\n\n'
            'No model chosen for production. No data/MC covariance or toy GOF claimed.\n'
            'Raw Hessian availability is not the stopping reason.\n'
            'All following pages are diagnostics, not preliminary measurements.',va='top',fontsize=11,linespacing=1.7)
        save(fig,'release_summary')
        fig,ax=plt.subplots(figsize=(12,9));ax.axis('off')
        ax.text(.01,.99,'Release gates: FAIL includes unvalidated gates deferred after upstream failure\n\n'+
                '\n'.join(f'{s:4s}  {g}' for g,s,e in gates),va='top',fontsize=11,linespacing=1.6)
        save(fig,'release_gates')
        fig,axs=plt.subplots(1,2,figsize=(11,4))
        axs[0].semilogy(spectrum,'o');axs[0].set(title='Conditional Sumw2 spectrum (not full V_data)',xlabel='mode',ylabel='eigenvalue')
        axs[1].semilogy(s/s[0],'o');axs[1].set(title='Column-scaled detector Jacobian',xlabel='mode',ylabel='relative singular value')
        save(fig,'conditional_rank')
        for k,component in enumerate(('U','LT','TT')):
            fig,ax=plt.subplots(figsize=(7,5))
            ax.plot([r['tprime_mean'] for r in central_table],z[k*nb:(k+1)*nb]*1e9,'o',color='gray')
            ax.set(xlabel="signed t' [GeV2]",ylabel=f'sigma_{component} [nb/GeV2]',
                   title=f'M0 diagnostic central values: {component}\nNOT FOR PUBLICATION; statistical errors unvalidated')
            save(fig,f'diagnostic_sigma_{component}')
        for sl in range(5):
            take=np.arange(nrow)//12==sl
            fig,axs=plt.subplots(2,1,figsize=(9,7),sharex=True)
            xx=np.arange(nrow)[take]%12
            axs[0].errorbar(xx,y[take],yerr=np.sqrt(var[take]),fmt='o',label='data, conditional fixed-weight SD')
            axs[0].plot(xx,pred[take],label='absolute folded M0');axs[0].legend()
            axs[0].set(title=f'Reconstructed slice {sl}: diagnostic only',ylabel='yield / mC')
            axs[1].plot(xx,(y-pred)[take],'o');axs[1].axhline(0,color='gray');axs[1].set(xlabel='phi bin',ylabel='raw residual')
            for j in np.where(take&~included)[0]:
                for ax in axs:ax.axvspan(j%12-.4,j%12+.4,color='red',alpha=.15)
            save(fig,f'detector_slice_{sl}')
        fig,ax=plt.subplots(figsize=(10,4))
        ax.plot((y[included]-pred[included])/np.sqrt(objvar[included]),'o')
        ax.axhline(0,color='gray');ax.set(xlabel='included row index',ylabel='conditional diagonal residual',
            title='Whitening with actual nominal objective metric; full data covariance unavailable')
        save(fig,'conditional_residual')
    print(json.dumps(dict(objective=status['objective'],ranks=rank_summary,exposure=exposure),indent=2))


if __name__ == '__main__':
    main()
