#!/usr/bin/env python3
"""Save uncertainty-study products; refuse production release on failed coverage.

Matrix dimensions follow the exported fitted-coordinate list. Fixed Q2/xB
feed-in is never represented as a covariance coordinate.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import write_csv
from check_pi0_objective_release import REQUIRED_GATES
from check_pi0_zero_release import ESTIMATOR as ZERO_ESTIMATOR, ZERO_GATES


def read(p):return list(csv.DictReader(p.open()))
def save(out,name,cov):
    np.savetxt(out/(name+'_covariance.csv'),cov,delimiter=',',fmt='%.17g')
    sd=np.sqrt(cov.diagonal());corr=np.divide(cov,np.outer(sd,sd),out=np.zeros_like(cov),where=np.outer(sd,sd)>0)
    np.savetxt(out/(name+'_correlation.csv'),corr,delimiter=',',fmt='%.17g')
    if not np.all(np.isfinite(cov)) or not np.allclose(cov,cov.T):raise RuntimeError('Invalid covariance')
    if np.linalg.eigvalsh(corr).min() < -1e-8:raise RuntimeError('Non-positive covariance correlation matrix')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('data','mc','response','coverage','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--release',action='store_true',help='Require validated uncertainty coverage before a production release')
    p.add_argument('--objective-validation',type=Path,help='Generated objective or fixed-zero-dataset readiness.json')
    a=p.parse_args();coverage=json.loads((a.coverage/'summary.json').read_text())
    objective_ready=False
    estimator='primitive_poisson'
    if a.objective_validation:
        validation=json.loads(a.objective_validation.read_text())
        estimator=validation.get('estimator','primitive_poisson')
        required=ZERO_GATES if estimator==ZERO_ESTIMATOR else REQUIRED_GATES
        objective_ready=(validation.get('production_release') is True
            and estimator in ('primitive_poisson',ZERO_ESTIMATOR)
            and set(required).issubset(validation.get('gates',{}))
            and all(v=='PASS' for v in validation['gates'].values())
            and bool(validation.get('inputs_sha256')))
        for path,expected in validation.get('inputs_sha256',{}).items():
            source=Path(path)
            objective_ready &= source.is_file() and hashlib.sha256(source.read_bytes()).hexdigest()==expected
        ensemble=validation.get('evidence',{}).get('ensemble_checks',{})
        objective_ready &= (Path(ensemble.get('data','')).resolve()==(a.data/'sigma_covariance.csv').resolve()
                            and Path(ensemble.get('mc','')).resolve()==(a.mc/'mc_canonical_nuisance_gauge_covariance.csv').resolve())
        if estimator==ZERO_ESTIMATOR:
            objective_ready &= (Path(validation.get('evidence',{}).get('coverage_path','')).resolve()==(a.coverage/'summary.json').resolve())
    if a.release and not objective_ready:
        raise SystemExit('RELEASE GATES INCOMPLETE: preliminary release requires the selected estimator\'s generated bias, coverage, dataset, covariance and published-rank gates for these ensembles. Supply --objective-validation readiness.json after successful validation.')
    if estimator==ZERO_ESTIMATOR:
        coverage=[r for r in coverage if r.get('case') in ('retained','low','low_half_acc','low_double_acc')]
    under=[r for r in coverage if abs(r['coverage']-.682689492)>2.58*r.get('coverage_mc_se',r.get('coverage_se',0.))]
    if a.release and under:
        raise SystemExit('NEW ANALYSIS-LEVEL ISSUE: controlled background-fit toys under-cover; production extraction release withheld. See coverage summary. Central values are preserved.')
    a.output.mkdir(parents=True,exist_ok=False)
    params=read(a.response/'migration_parameters.csv');pub=np.array([r['region']=='published' for r in params])
    data=np.loadtxt(a.data/'sigma_covariance.csv',delimiter=',');mc=np.loadtxt(a.mc/'mc_canonical_nuisance_gauge_covariance.csv',delimiter=',')
    if data.shape!=mc.shape or data.shape!=(len(params),len(params)):raise RuntimeError('Covariance dimension mismatch')
    nominal=read(a.data/'sigma_summary.csv')
    if not np.allclose([float(r['value']) for r in params],[float(r['nominal']) for r in nominal],rtol=1e-10,atol=1e-20):raise RuntimeError('Central value mismatch')
    total=data+mc
    for name,c in [('data',data),('mc_canonical_nuisance_gauge',mc),('total_canonical_nuisance_gauge',total)]:save(a.output,name,c)
    for name,c in [('data_published',data),('mc_published',mc),('total_statistical_published',total)]:save(a.output,name,c[np.ix_(pub,pub)])
    write_csv(a.output/'diagnostic_cross_sections.csv',[dict(parameter_index=i,truth_block=r['truth_block'],component=r['component'],region=r['region'],
        nominal=float(r['value']),data_sd=np.sqrt(data[i,i]),mc_sd=np.sqrt(mc[i,i]),total_statistical_sd=np.sqrt(total[i,i]),
        uncertainty_identifiable=int(pub[i]),validated_for_preliminary=int(not under and objective_ready)) for i,r in enumerate(params)])
    write_csv(a.output/'published_parameter_map.csv',[dict(matrix_index=j,parameter_index=i,truth_block=params[i]['truth_block'],component=params[i]['component']) for j,i in enumerate(np.flatnonzero(pub))])
    result=dict(verdict='READY FOR PRELIMINARY' if not under and objective_ready else 'RELEASE GATES INCOMPLETE',production_release=bool(not under and objective_ready and a.release),
        estimator=estimator,objective_validation=str(a.objective_validation) if a.objective_validation else None,
        issue='Selected-estimator coverage failed' if under else ('Selected-estimator release gates incomplete' if not objective_ready else None),
        undercoverage_cases=under,central_values_changed=False,
        statistical_independence='data and MC event samples independent conditional on fixed smearing calibration, run corrections and generated normalization',
        matrices_full_parameter=f'{len(params)} fitted coordinates; fixed Q2/xB feed-in has no covariance coordinate; any remaining tprime nuisance gauge is explicitly declared',
        matrices_published='15 uniquely estimable published coordinates; release status is determined by the selected-estimator gates')
    (a.output/'readiness.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
