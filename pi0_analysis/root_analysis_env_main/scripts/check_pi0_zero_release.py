#!/usr/bin/env python3
"""Manifest-aware release gate for a fixed zero-combinatorial-background dataset.

Positive-background toy failures are outside this estimator's scope. Evidence
that the retained data can use B=0 is mandatory; old fit labels are insufficient.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from check_pi0_objective_release import covariance

ESTIMATOR='fixed_zero_background'
ZERO_GATES=('production_run_manifest','zero_background_estimator_bias',
 'low_signal_zero_background_bias','zero_background_coverage',
 'retained_background_adequacy','exposure_consistency','data_covariance',
 'full_MC_covariance','published_rank')

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def evaluate(manifest,coverage,dataset_validation=None,data=None,mc=None):
    checksums={}
    def read(p):
        checksums[str(p.resolve())]=digest(p);return json.loads(p.read_text())
    checksums[str(manifest.resolve())]=digest(manifest)
    with manifest.open() as stream:
        rows=list(csv.DictReader(stream))
    accepted=[r for r in rows if r['accepted']=='1'];runs=[int(r['run']) for r in accepted]
    gates={k:False for k in ZERO_GATES};evidence={};stats=read(coverage)
    bycase={r['case']:r for r in stats}
    valid_manifest=bool(runs) and len(runs)==len(set(runs)) and all(float(r['charge_uC'])>0 and float(r['livetime'])>0 and float(r['efficiency'])>0 and r['input_available']=='1' for r in accepted)
    adequate=read(dataset_validation) if dataset_validation and dataset_validation.is_file() else {}
    linked=adequate.get('manifest_sha256')==digest(manifest)
    # A nominal zero branch is not a bound on residual contamination.
    gates['production_run_manifest']=valid_manifest and linked and adequate.get('frozen_for_production') is True
    gates['retained_background_adequacy']=linked and adequate.get('validated_zero_background') is True and 0<=adequate.get('max_background_bias_over_data_sigma',float('inf'))<=.1 and bool(adequate.get('inputs_sha256'))
    for name,sha in adequate.get('inputs_sha256',{}).items():
        p=Path(name);good=p.is_file() and digest(p)==sha
        gates['retained_background_adequacy'] &= good
        if good:checksums[str(p.resolve())]=sha
    def bias(r):return r.get('experiments',0)>=1000 and abs(r['bias_sigma'])<=.1 and r['yield_bias_se']/r['empirical_sd']<=.04
    gates['zero_background_estimator_bias']='retained' in bycase and bias(bycase['retained'])
    gates['low_signal_zero_background_bias']=all(k in bycase and bias(bycase[k]) for k in ('low','low_half_acc','low_double_acc'))
    covchecks={}
    for k in ('retained','low','low_half_acc','low_double_acc'):
        r=bycase.get(k,{})
        if not r:covchecks[k]=False;continue
        n=r['experiments'];covchecks[k]=bool(bias(r) and abs(r['pull_mean'])<=2.58*r['pull_width']/np.sqrt(n) and abs(r['pull_width']-1)<=max(.05,2.58/np.sqrt(2*(n-1))) and abs(r['coverage']-.682689492137086)<=2.58*r['coverage_se'] and r['coverage_se']<=.021)
    gates['zero_background_coverage']=all(covchecks.values());evidence['coverage_checks']=covchecks
    provenance=[]
    for directory in (data,mc):
        p=directory/'estimator_provenance.json' if directory else None
        provenance.append(read(p) if p and p.is_file() else {})
    same=all(provenance) and all(r.get('estimator')==ESTIMATOR and r.get('manifest_sha256')==digest(manifest) for r in provenance) and bool(provenance[0].get('nominal_sha256')) and provenance[0]['nominal_sha256']==provenance[1].get('nominal_sha256')
    exact=bool(same)
    exposures={str(r['run']):{k:float(r[k]) for k in ('charge_uC','livetime','efficiency','prescale')} for r in accepted}
    for pr in provenance:
        for key in ('yield_runs','charge_runs','efficiency_runs','livetime_runs','bootstrap_runs'):
            values=pr.get(key,[]);exact &= len(values)==len(set(values)) and sorted(values)==sorted(runs)
        exact &= pr.get('run_exposures')==exposures
    gates['exposure_consistency']=bool(exact)
    dok,dnote,dcov=covariance(data,'sigma_covariance.csv');rok,rnote,_=covariance(data,'reco_covariance.csv',(60,60));mok,mnote,mcov=covariance(mc,'mc_canonical_nuisance_gauge_covariance.csv')
    matching_dimensions=bool(dok and mok and dcov.shape==mcov.shape)
    exact &= matching_dimensions
    for directory,name in ((data,'sigma_covariance.csv'),(data,'reco_covariance.csv'),(mc,'mc_canonical_nuisance_gauge_covariance.csv')):
        if directory and (directory/name).is_file():checksums[str((directory/name).resolve())]=digest(directory/name)
    complete=lambda r:r.get('requested_replicas',0)>=2000 and r.get('successful_replicas')==r.get('requested_replicas') and r.get('failed_replicas')==0
    gates['data_covariance']=bool(dok and rok and exact and complete(provenance[0]))
    gates['full_MC_covariance']=bool(mok and exact and complete(provenance[1]) and provenance[1].get('method')=='full_event_poisson')
    gates['published_rank']=bool(exact and provenance[1].get('all_published_estimable') is True and complete(provenance[1]))
    evidence.update(coverage_path=str(coverage.resolve()),accepted_runs=runs,accepted_charge_uC=sum(float(r['charge_uC']) for r in accepted),ensemble_checks=dict(data=dnote,reco=rnote,mc=mnote,matching_dimensions=matching_dimensions,matching_provenance=bool(same)))
    ready=all(gates.values())
    return dict(estimator=ESTIMATOR,production_release=ready,verdict='READY - CONSERVATIVE RUN SELECTION' if ready else 'ZERO-DATASET RELEASE GATES INCOMPLETE',gates={k:'PASS' if v else 'FAIL' for k,v in gates.items()},evidence=evidence,inputs_sha256=checksums)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('manifest','coverage','output'):p.add_argument('--'+k,type=Path,required=True)
    for k in ('dataset-validation','data','mc'):p.add_argument('--'+k,type=Path)
    p.add_argument('--require-ready',action='store_true');a=p.parse_args()
    result=evaluate(a.manifest,a.coverage,a.dataset_validation,a.data,a.mc)
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'readiness.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if a.require_ready and not result['production_release']:raise SystemExit(2)

if __name__=='__main__':main()
