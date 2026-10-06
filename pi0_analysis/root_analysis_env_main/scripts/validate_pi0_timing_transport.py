#!/usr/bin/env python3
"""Diagnostic timing transport decision; never publishes a production weight.

Uses every accepted run, the existing fixed nominal ellipse, and the production
C++ mass fitter. No generated non-exclusive component is loaded or simulated.
Output directories must be new. Load Hall C ROOT before the toys subcommand.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import uproot
from bootstrap_pi0_data import Fitter, RunEvents, MASS_EDGES, timing_components
from audit_pi0_exclusive_estimator import rows, ellipse

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / 'validation/xsec_input_recovery_20261004/output/KinC_x36_4/root/combined_branches_LH2.root'
COEFF = np.array([1., -1/6, -1/12, -1/12, 1/18, 1/18])


def write(path, records):
    tmp = path.with_suffix(path.suffix + '.tmp')
    with tmp.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader(); writer.writerows(records)
    tmp.replace(path)


def geometry(data):
    result = {}
    for line in data.with_name(data.stem + '_combined_2d_mass_cut_debug.txt').read_text().splitlines():
        if '=' in line:
            key, value = line.split('=', 1)
            try: result[key] = float(value)
            except ValueError: pass
    return result


def candidates(c, mb, n, t, s):
    """No clipping or near-zero division regularization. NaN means undefined."""
    b = t-s
    old_bin = np.divide(s, n, out=np.zeros(200), where=n > 0)
    mean = np.divide(t, n, out=np.full(200, np.nan), where=n > 0)
    add = c - (mean-old_bin)[mb]
    factor = np.divide(s, t, out=np.full(200, np.nan), where=t != 0)
    # Explicit continuous extension for exactly B=0, including T=0.
    factor[b == 0] = 1.
    mult = c*factor[mb]
    return old_bin[mb], add, mult


def audit(args):
    with uproot.open(args.data) as source:
        events = source['physics'].arrays(library='np')
        runs = source['analysis_runs'].arrays(library='np')
    d = geometry(args.data)
    selected = ellipse(events['mpi0_all'], events['mmiss_all'], d)
    np.testing.assert_array_equal(selected, events['is_exclusive_ellipse_combined'] != 0)
    reco, kinematics = rows(events)
    charge = runs['charge_uC'].astype(float).sum()
    scale = events['scale'].astype(float)*events['charge_uC'].astype(float)/charge/.584
    closure, transport = [], []
    maximum = dict(timing=0., old=0., additive=0., multiplicative=0.)
    for run in runs['run_number']:
        take = events['run_number'] == run
        e = {k:v[take] for k,v in events.items()}
        mb = np.searchsorted(MASS_EDGES, e['mpi0_all'], side='right')-1
        valid = (mb>=0)&(mb<200)
        c = COEFF @ timing_components(e)
        with uproot.open(args.data.parent/f'diagnostics_run{run}.root') as f:
            n = f[f'h_mpi0_all_run{run}'].values()
            t = f['h_pi0_coin_bgsub'].values()
            s = f['h_pi0_final'].values()
        hist = lambda w: np.bincount(mb[valid], weights=w[valid], minlength=200)
        np.testing.assert_array_equal(hist(np.ones(len(mb))), n)
        np.testing.assert_allclose(hist(c), t, rtol=1e-12, atol=1e-12)
        maximum['timing'] = max(maximum['timing'], float(np.max(abs(hist(c)-t))))
        ws = [np.zeros(len(mb)) for _ in range(3)]
        for dest, src in zip(ws, candidates(c[valid], mb[valid], n, t, s)):
            dest[valid] = src
        np.testing.assert_allclose(ws[0], e['pi0_weight'], rtol=1e-12, atol=1e-14)
        sums = [hist(w) for w in ws]
        b = t-s
        for j in range(200):
            status = []
            if n[j] == 0:
                status.append('empty_zero' if s[j] == 0 else 'empty_nonzero_residual_UNREPRESENTABLE')
            else:
                for method, total in zip(('old','additive'), sums[:2]):
                    np.testing.assert_allclose(total[j], s[j], rtol=1e-10, atol=1e-10)
                    maximum[method] = max(maximum[method], abs(float(total[j]-s[j])))
                if b[j] == 0 and t[j] == 0: status.append('mult_zero_background_identity_extension')
                elif t[j] == 0: status.append('mult_UNDEFINED_T_zero')
                elif abs(t[j]) < 1e-12*max(1., n[j]): status.append('mult_NEAR_ZERO_T_unregularized')
                if np.isfinite(sums[2][j]) and not any('NEAR_ZERO' in x for x in status):
                    np.testing.assert_allclose(sums[2][j], s[j], rtol=1e-10, atol=1e-10)
                    maximum['multiplicative'] = max(maximum['multiplicative'], abs(float(sums[2][j]-s[j])))
            closure.append(dict(run=int(run),mass_bin=j+1,N_rb=n[j],T_rb=t[j],B_rb=b[j],S_rb=s[j],
                old_sum=sums[0][j],additive_sum=sums[1][j],multiplicative_sum=sums[2][j],
                status=';'.join(status) or 'PASS'))
        status = next(csv.DictReader((args.data.parent/f'analysis_status_run{run}.csv').open()))
        for scope, selection in [('ellipse_only',selected[take]),('ellipse_and_reco',selected[take]&kinematics[take])]:
            selection &= valid
            factor = scale[take][selection]
            transport.append(dict(run=int(run),scope=scope,classification=status['classification'],
                old_sum=float(np.sum(ws[0][selection]*factor)),
                additive_sum=float(np.sum(ws[1][selection]*factor)),
                multiplicative_sum=float(np.sum(ws[2][selection]*factor)),
                direct_timing_sum=float(np.sum(c[selection]*factor)),
                selected_events=int(selection.sum())))
    summary = []
    for scope in ('ellipse_only','ellipse_and_reco'):
        r = [x for x in transport if x['scope']==scope and x['classification']=='zero_background']
        totals = {k:sum(x[k] for x in r) for k in ('old_sum','additive_sum','multiplicative_sum','direct_timing_sum')}
        for k in ('additive_sum','multiplicative_sum'):
            np.testing.assert_allclose(totals[k], totals['direct_timing_sum'], atol=1e-12, rtol=1e-12)
        summary.append(dict(scope=scope,runs=len(r),**totals,
            old_discrepancy_percent=100*(totals['direct_timing_sum']/totals['old_sum']-1),units='per_mC'))
    write(args.output/'timing_weight_inclusive_closure.csv',closure)
    write(args.output/'ellipse_transport_by_run.csv',transport)
    write(args.output/'ellipse_transport_summary.csv',summary)
    report=dict(runs=len(runs['run_number']),events=len(selected),bins=len(closure),max_absolute_defined_closure_error=maximum,
        empty_nonzero_residual_bins=sum('UNREPRESENTABLE' in x['status'] for x in closure),
        empty_nonzero_residual_sum=sum(x['S_rb'] for x in closure if 'UNREPRESENTABLE' in x['status']),
        mult_undefined_bins=sum('UNDEFINED' in x['status'] for x in closure),
        mult_near_zero_bins=sum('NEAR_ZERO' in x['status'] for x in closure),geometry=d,transport=summary)
    (args.output/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


def synthetic(rng, d, signal_mean, background_mean, accidental_mean):
    """One synthetic run with six disjoint production timing categories.

    Accidentals have a common joint mass/missing-mass distribution in all
    timing categories, hence cancel exactly in expectation even after the cut.
    The category rates (1,6,6,6,9,9) follow the actual window areas.
    """
    cov=np.array([[d['cov_mpi0_mpi0'],d['cov_mpi0_mmiss']],
                  [d['cov_mpi0_mmiss'],d['cov_mmiss_mmiss']]])
    signal=rng.multivariate_normal([d['mean_mpi0'],d['mean_mmiss']],cov,size=rng.poisson(signal_mean))
    centers=(MASS_EDGES[:-1]+MASS_EDGES[1:])/2
    prob=1/(1+np.exp((centers-.145)/.025));prob/=prob.sum()
    nb=rng.poisson(background_mean)
    bm=rng.choice(centers,nb,p=prob)+rng.uniform(-.001,.001,nb)
    parts=[(signal[:,0],signal[:,1],np.zeros(len(signal),int),np.ones(len(signal),bool)),
           (bm,rng.uniform(.6,1.5,nb),np.zeros(nb,int),np.zeros(nb,bool))]
    for category, ratio in enumerate((1,6,6,6,9,9)):
        count=rng.poisson(accidental_mean*ratio)
        # Mass-peaked accidentals test timing transport separately from B.
        m=rng.normal(d['mean_mpi0'],np.sqrt(cov[0,0]),count)
        parts.append((m,rng.uniform(.6,1.5,count),np.full(count,category),np.zeros(count,bool)))
    m,x,cat,label=(np.concatenate([p[i] for p in parts]) for i in range(4))
    keep=(m>=0)&(m<.4)
    m,x,cat,label=(v[keep] for v in (m,x,cat,label))
    # Representative points strictly inside each timing category.
    times=np.array([[150,150],[140,140],[150,140],[140,150],[158,142],[142,158]])
    e=dict(mpi0_all=m,mmiss_all=x,photon_time_1=times[cat,0],photon_time_2=times[cat,1])
    run=RunEvents.__new__(RunEvents)
    run.events=e;run.mass_bin=np.searchsorted(MASS_EDGES,m,side='right')-1
    run.inmass=np.ones(len(m),bool);run.components=timing_components(e)
    run.plane_components=timing_components(e,True)
    run.template_factors=-COEFF.copy();run.template_factors[0]=0
    run.coefficient=COEFF@run.components
    return run,label,ellipse(m,x,d)


def toys(args):
    fit=Fitter(args.library);d=geometry(args.data)
    rng=np.random.default_rng(args.seed);records=[];failures=[];summary=[]
    # High-statistics stress test isolates downstream transport; the continuum
    # is EXACTLY the fitted Fermi family, with all truth labels retained.
    cases=[('clean',0.,0.),('clean_timing',0.,150.),
           ('clean_comb',args.background,0.),('clean_timing_comb',args.background,150.)]
    for case,background,accidental in cases:
        successful=attempts=0
        while successful<args.toys and attempts<args.max_attempts:
            attempts+=1
            run,label,selection=synthetic(rng,d,2000.,background,accidental)
            t,v,n=run.spectra(np.ones(len(label)))
            code,s,info=fit.fit(t,v)
            if code:
                failures.append(dict(case=case,attempt=attempts,code=code,minimizer=info[1],covariance=info[2]))
                continue
            weights=candidates(run.coefficient,run.mass_bin,n,t,s)
            truth=int(np.sum(label&selection))
            # Oracle B controls separate transport bias from fit bias. It is
            # diagnostic only and never used to choose a real-data weight.
            centers=(MASS_EDGES[:-1]+MASS_EDGES[1:])/2
            b=1/(1+np.exp((centers-.145)/.025));b*=background/b.sum()
            oracle=candidates(run.coefficient,run.mass_bin,n,t,t-b)
            for method,w,ow in zip(('old','additive','multiplicative'),weights,oracle):
                defined=bool(np.all(np.isfinite(w[selection])))
                estimate=float(w[selection].sum()) if defined else np.nan
                records.append(dict(case=case,replica=successful,attempt=attempts,method=method,
                    truth=truth,estimate=estimate,residual=estimate-truth,
                    oracle_estimate=float(ow[selection].sum()),
                    fit_background=float(info[6]),true_mean_background=background,
                    empty_residual=float(s[n==0].sum()),
                    undefined_selected_events=int(np.sum(~np.isfinite(w[selection])))))
            successful+=1
            if successful%25==0: print(f'{case}: {successful} valid mass fits / {attempts} attempts',flush=True)
        for method in ('old','additive','multiplicative'):
            rr=[r for r in records if r['case']==case and r['method']==method]
            valid=[r for r in rr if np.isfinite(r['estimate'])]
            residual=np.array([r['residual'] for r in valid]);spread=float(residual.std(ddof=1)) if len(valid)>1 else np.nan
            bias=float(residual.mean()) if valid else np.nan
            summary.append(dict(case=case,method=method,requested_toys=args.toys,attempted_toys=attempts,
                successful_toys=len(valid),fit_failures=attempts-successful,undefined_toys=len(rr)-len(valid),
                truth_mean=float(np.mean([r['truth'] for r in valid])) if valid else np.nan,
                estimated_mean=float(np.mean([r['estimate'] for r in valid])) if valid else np.nan,
                bias=bias,empirical_residual_spread=spread,
                empirical_estimate_spread=float(np.std([r['estimate'] for r in valid],ddof=1)) if len(valid)>1 else np.nan,
                bias_over_spread=bias/spread if spread else 0.,
                bias_over_mean_error=bias/spread*np.sqrt(len(valid)) if spread else 0.,
                oracle_bias=float(np.mean([r['oracle_estimate']-r['truth'] for r in valid])) if valid else np.nan))
        write(args.output/'toy_replicas.csv',records)
        write(args.output/'toy_summary.csv',summary)
        write(args.output/'toy_failures.csv',failures or [dict(case='none',attempt=0,code=0,minimizer=0,covariance=0)])
    config=dict(seed=args.seed,requested_successful_per_case=args.toys,signal_mean=2000,
        cases=cases,geometry=d,geometry_refitted=False,fitter='production nps_stat_bridge.cpp / nps_comb_bg_pepsi.h',
        toy_scope='one synthetic run per replica; production timing and mass fitting; fixed nominal ellipse',
        note='Stress rates are declared assumptions, not measured continuum rates. Oracle is diagnostic only.')
    (args.output/'toy_configuration.json').write_text(json.dumps(config,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


def verify(args):
    """Production regeneration parity, including timing-fitter input errors."""
    if not args.production_dir: raise ValueError('--production-dir required')
    with uproot.open(args.production_dir/'combined_branches_LH2.root') as f:
        combined=f['physics'].arrays(library='np')
        run_ids=f['analysis_runs']['run_number'].array(library='np')
    evidence=[]
    for run in run_ids:
        with uproot.open(args.production_dir/f'diagnostics_run{run}.root') as f, \
             uproot.open(args.data.parent/f'diagnostics_run{run}.root') as old:
            e=f['physics'].arrays(library='np');prior=old['physics'].arrays(library='np')
            for key in prior: np.testing.assert_array_equal(e[key],prior[key])
            take=combined['run_number']==run
            for key in e: np.testing.assert_array_equal(combined[key][take],e[key])
            c=COEFF@timing_components(e)
            np.testing.assert_allclose(e['pi0_timing_coeff'],c,rtol=1e-13,atol=1e-14)
            mb=np.searchsorted(MASS_EDGES,e['mpi0_all'],side='right')-1
            valid=(mb>=0)&(mb<200)
            t=f['h_pi0_coin_bgsub'].values();n=f[f'h_mpi0_all_run{run}'].values()
            expected=t[mb[valid]]/n[mb[valid]]
            np.testing.assert_allclose(e['pi0_timing_bin_mean'][valid],expected,rtol=1e-13,atol=1e-14)
            assert np.all(np.isnan(e['pi0_timing_bin_mean'][~valid]))
            raw=RunEvents(args.production_dir/f'signal_events_run{run}.root')
            ty,tv,tn=raw.spectra(np.ones(len(raw.coefficient)))
            np.testing.assert_array_equal(tn,n)
            np.testing.assert_allclose(ty,t,rtol=1e-12,atol=1e-12)
            np.testing.assert_allclose(tv,f['h_pi0_coin_bgsub'].variances(),rtol=1e-12,atol=1e-12)
            evidence.append(dict(run=int(run),events=len(c),legacy_branches_exact=len(prior),
                new_branches=sorted(set(e)-set(prior)),all_events_retained_in_combination=True,
                timing_fitter_contents_and_variances_match=True))
    (args.output/'production_parity.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(evidence),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['audit','toys','verify']);p.add_argument('--data',type=Path,default=DATA)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--library',type=Path)
    p.add_argument('--production-dir',type=Path)
    p.add_argument('--seed',type=int,default=20261004);p.add_argument('--toys',type=int,default=200)
    p.add_argument('--max-attempts',type=int,default=1000)
    p.add_argument('--background',type=float,default=6000.,help='Injected continuum mean; signal mean is 2000')
    args=p.parse_args()
    if args.mode=='toys' and not args.library:p.error('--library required')
    args.output.mkdir(parents=True,exist_ok=False)
    dict(audit=audit,toys=toys,verify=verify)[args.mode](args)


if __name__=='__main__':main()
