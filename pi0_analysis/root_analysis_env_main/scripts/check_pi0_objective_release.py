#!/usr/bin/env python3
"""Evaluate objective validation and matching ensembles; fail closed for release.

Diagnostic screening only until every gate passes. Missing artifacts are failures,
never inferred from a prior estimator's covariance or fit convergence.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np

TARGET=.682689492137086
REQUIRED_CASES=(0,1,2,3,4,5)
REQUIRED_GATES=('central_estimator_toy_bias','zero_background_coverage',
    'small_positive_background_coverage','clearly_positive_background_coverage',
    'low_signal_positive_background_coverage','low_signal_small_positive_background_coverage',
    'boundary_transition','data_bootstrap_covariance_available',
    'full_MC_covariance_available','published_parameters_identifiable')


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def bias_ok(r):
    # A bias below 0.1 sigma is negligible at preliminary precision; a larger
    # one must be statistically compatible with zero at the stated toy precision.
    return abs(r['yield_bias'])<=max(.1*r['empirical_sd'],2.58*r['yield_bias_se'])


def coverage_checks(r):
    n=r['successful'];pull_se=r['pull_width']/np.sqrt(n)
    checks=dict(bias=bias_ok(r),
        pull_mean=abs(r['pull_mean'])<=2.58*pull_se,
        pull_width=abs(r['pull_width']-1)<=max(.1,2.58/np.sqrt(2*(n-1))),
        coverage=abs(r['coverage']-TARGET)<=2.58*r['coverage_se'],
        precision=r['coverage_se']<=.021,
        failures=(r.get('failed',0)+r.get('outer_discarded',0))/r['requested']<=.01)
    return {name:bool(value) for name,value in checks.items()}


def covariance(directory,name,expected=None):
    if directory is None:return False,'not regenerated for the candidate estimator',None
    path=directory/name
    if not path.is_file():return False,f'missing {path}',None
    matrix=np.loadtxt(path,delimiter=',')
    valid=matrix.ndim==2 and matrix.shape[0]==matrix.shape[1]
    if expected is not None:valid &= matrix.shape==expected
    valid &= np.all(np.isfinite(matrix)) and np.allclose(matrix,matrix.T,rtol=1e-10,atol=1e-20)
    if valid:
        scale=np.max(np.abs(matrix));valid &= np.linalg.eigvalsh(matrix).min()>=-1e-10*max(scale,1e-300)
    return bool(valid),str(path),matrix


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--screen',type=Path,nargs='+',required=True)
    p.add_argument('--coverage',type=Path,required=True);p.add_argument('--boundary',type=Path,required=True)
    p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--data',type=Path);p.add_argument('--mc',type=Path)
    p.add_argument('--require-ready',action='store_true')
    a=p.parse_args();gates={};evidence={};checksums={};sha=digest(a.library)
    def read(path):
        checksums[str(path.resolve())]=digest(path);return json.loads(path.read_text())
    screen=[]
    for directory in a.screen:screen.extend(read(directory/'summary.json'))
    bycase={r['case']:r for r in screen if r['mode']==0}
    gates['central_estimator_toy_bias']=all(c in bycase and bias_ok(bycase[c]) and bycase[c]['successful']>=1000 for c in REQUIRED_CASES)
    evidence['central_bias']={str(c):dict(pass_bias=bias_ok(r),absolute_bias_over_sigma=r['absolute_bias_over_sigma']) for c,r in bycase.items()}
    rows=read(a.coverage/'summary.json');bycov={r['case']:r for r in rows if r['mode']==0}
    names={0:'zero_background_coverage',1:'small_positive_background_coverage',2:'clearly_positive_background_coverage',
           4:'low_signal_positive_background_coverage',5:'low_signal_small_positive_background_coverage'}
    for c,name in names.items():
        checks=coverage_checks(bycov[c]) if c in bycov else {'present':False}
        gates[name]=all(checks.values());evidence[name]=checks
    boundary=read(a.boundary/'summary.json')
    gates['boundary_transition']=len(boundary)>=5 and all(all(coverage_checks(r).values()) for r in boundary)
    evidence['boundary_checks']=[dict(truth_A=r['truth_A'],checks=coverage_checks(r)) for r in boundary]
    data_ok,data_note,data=covariance(a.data,'sigma_covariance.csv')
    reco_ok,reco_note,_=covariance(a.data,'reco_covariance.csv',(60,60))
    mc_ok,mc_note,mc=covariance(a.mc,'mc_canonical_nuisance_gauge_covariance.csv')
    matching_dimensions=bool(data_ok and mc_ok and data.shape==mc.shape)
    for directory,name in ((a.data,'sigma_covariance.csv'),(a.data,'reco_covariance.csv'),(a.mc,'mc_canonical_nuisance_gauge_covariance.csv')):
        if directory is not None and (directory/name).is_file():checksums[str((directory/name).resolve())]=digest(directory/name)
    # Ensembles must explicitly bind to the corrected estimator/nominal vector.
    provenance=[]
    for directory in (a.data,a.mc):
        if directory is None or not (directory/'objective_provenance.json').is_file():provenance.append(None)
        else:provenance.append(read(directory/'objective_provenance.json'))
    same=(all(provenance) and matching_dimensions and all(r.get('objective_library_sha256')==sha for r in provenance)
          and provenance[0].get('nominal_sha256') is not None
          and provenance[0]['nominal_sha256']==provenance[1].get('nominal_sha256'))
    gates['data_bootstrap_covariance_available']=bool(data_ok and reco_ok and same)
    gates['full_MC_covariance_available']=bool(mc_ok and same and provenance[1].get('method')=='full_event_poisson')
    gates['published_parameters_identifiable']=bool(same and provenance[1].get('all_published_estimable') is True
        and provenance[1].get('requested_replicas',0)>=2000 and provenance[1].get('successful_replicas')==provenance[1].get('requested_replicas'))
    if a.data and provenance[0]:
        gates['data_bootstrap_covariance_available'] &= provenance[0].get('requested_replicas',0)>=2000
    evidence['ensemble_checks']=dict(data=data_note,reco=reco_note,mc=mc_note,matching_dimensions=matching_dimensions,matching_provenance=bool(same))
    ready=set(gates)==set(REQUIRED_GATES) and all(gates.values())
    result=dict(verdict='READY FOR PRELIMINARY' if ready else 'OBJECTIVE STILL INVALID',production_release=ready,
        gates={k:'PASS' if v else 'FAIL' for k,v in gates.items()},evidence=evidence,
        objective_library_sha256=sha,inputs_sha256=checksums,
        policy=dict(practical_bias_sigma=.1,compatibility_z=2.58,target_coverage=TARGET,max_coverage_se=.021,
                    pull_width_minimum_tolerance=.1,max_outer_failure_fraction=.01))
    a.output.mkdir(parents=True,exist_ok=False)
    (a.output/'readiness.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if a.require_ready and not ready:raise SystemExit(2)


if __name__=='__main__':main()
