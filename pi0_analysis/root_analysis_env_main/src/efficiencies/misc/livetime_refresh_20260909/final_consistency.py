"""Cross-check independently assembled counts and preserve production sources."""
from pathlib import Path
import json,csv,hashlib,numpy as np
P=Path(__file__).resolve().parent
R=[r for r in json.loads((P/'audited_results.json').read_text()) if r['status']=='CALCULATED' and r['run_type']=='production']
assert len(R)==43 and sum(r['segments'] for r in R)==183
rows=[];checks=[];phase=[]
for r in R:
    run=r['run'];n=r['nominal'];p=r['prescale_factor'];q=json.loads((P/f'timestamp/run{run}.json').read_text())['segments'];phase.extend(q)
    assert sum(x['raw']['total'] for x in q)==n['E_raw']
    assert sum(x['core']['total'] for x in q)==n['E_tight']
    lp=p*n['E_raw']/n['D'];lc=p*(n['N']-n['E_raw'])/(n['S']-n['D'])
    residual=n['CLT_all']-((1-n['D']/n['S'])*lc+(n['D']/n['S'])*lp)
    assert abs(residual)<1e-12
    checks.append(dict(run=run,phase_counts_equal_interval_counts=True,shared_count_identity_residual=residual))
    rows.append(dict(run=run,segments=r['segments'],events=r['events'],trigger=r['trigger'],p=p,N=n['N'],E_raw=n['E_raw'],E_tight=n['E_tight'],D=n['D'],S=n['S'],A=n['A'],saved_NewGen=float(r['saved']['NewGen_EDTM_livetime']),same_cache_NewGen=r['old_ratio'],matched_raw_EDTM=lp,matched_tight_EDTM=n['EDTM_tight'],all_trigger_ratio=n['CLT_all'],raw_subtracted_physics_ratio=lc,counter_suspect_intervals=r['closure_exclusion_diagnostic']['removed_intervals'],same_saved_files=r['same_saved_files']))
with (P/'production_summary.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
ph={k:{**{v:sum(x[k][v] for x in phase) for v in ['total','phase125','side500_1500','model_missing']},'phase_windows':{w:sum(x[k]['phase_windows'][w] for x in phase) for w in ['10','25','50','125','250','500']}} for k in ['core','raw','raw_not_core','positive_not_core','zero_raw','all']}
(P/'phase_summary.json').write_text(json.dumps(ph,indent=2)+'\n')
acceptance=json.loads((P/'acceptance_tests.json').read_text());combined=[]
for variant in ['raw','tight']:
    rr=[r for r in acceptance if r['variant']==variant and r['p']==1 and not r['counter_suspect_intervals']]
    d=np.array([r['difference'] for r in rr]);s=np.array([r['block20_difference_sigma'] for r in rr]);w=1/s**2;mean=float(np.sum(w*d)/sum(w));sigma=float(1/np.sqrt(sum(w)))
    combined.append(dict(variant=variant,runs=len(rr),inverse_variance_mean_difference=mean,independent_run_sigma=sigma,z=mean/sigma,chi2_about_mean=float(sum((d-mean)**2/s**2)),positive_differences=int(sum(d>0)),assumption='Common mean and independent run-level block errors; no shared routing, calibration or sampling systematic included'))
(P/'combined_residual_diagnostic.json').write_text(json.dumps(combined,indent=2)+'\n')
source=[]
for r in json.loads((P/'source_provenance.json').read_text()):
    digest=hashlib.sha256(Path(r['source']).read_bytes()).hexdigest();source.append(dict(path=r['source'],sha256=digest,unchanged=digest==r['sha256']))
assert all(x['unchanged'] for x in source)
base=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc');archives=[]
for name in ['livetime_KinC_x60_4b_lh2_20260909','livetime_beamer_20260909']:
    folder=base/name;manifest=folder/'SHA256SUMS'
    if not manifest.exists():archives.append(dict(archive=name,manifest_absent=True));continue
    n=0
    for line in manifest.read_text().splitlines():
        digest,rel=line.split(None,1);rel=rel.lstrip('*');actual=hashlib.sha256((folder/rel).read_bytes()).hexdigest()
        assert actual==digest,(name,rel,'immutable archive mismatch');n+=1
    archives.append(dict(archive=name,files_verified=n))
result=dict(production_runs=len(R),production_segments=sum(r['segments'] for r in R),production_events=sum(r['events'] for r in R),count_checks=checks,production_sources=source,previous_archives=archives)
(P/'refresh_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print('Validated',len(R),'production runs; phase counts and mixture identity agree.')
print('Production sources unchanged:',len(source));print('Previous archives:',archives)
