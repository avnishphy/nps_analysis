#!/usr/bin/env python3
"""Report data-bootstrap, correction-factor, and first-order finite-MC scales.

The response estimate uses the production Poissonized response-cell moments;
it is a delta-method diagnostic, not a new response or smearing calibration.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from bootstrap_pi0_data import RunEvents,write_csv


def rows(path):
    return list(csv.DictReader(path.open()))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bootstrap',type=Path,required=True)
    parser.add_argument('--extraction',type=Path,required=True)
    parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    manifest=rows(args.manifest);accepted=[r for r in manifest if r['accepted']=='1']
    params=rows(args.extraction/'migration_parameters.csv')
    cells=rows(args.extraction/'migration_response_cells.csv')
    reco=rows(args.extraction/'migration_reco_rows.csv')
    design=np.zeros((60,len(params)))
    for r in rows(args.extraction/'migration_design.csv'):
        design[int(r['reco_row']),int(r['parameter_index'])]=float(r['response'])
    var=np.array([float(r['data_variance']) for r in reco]);keep=var>0
    a=design[keep]/np.sqrt(var[keep,None]);scale=np.linalg.norm(a,axis=0)
    u,s,vt=np.linalg.svd(a/scale,full_matrices=False)
    left=((vt.T/s)@u.T)/scale[:,None]/np.sqrt(var[None,keep])
    operator=np.zeros((len(params),60));operator[:,keep]=left
    data_cov=np.loadtxt(args.bootstrap/'sigma_covariance.csv',delimiter=',')
    data_sd=np.sqrt(np.diag(data_cov))
    values=np.array([float(r['value']) for r in params])
    blocks={int(r['truth_block']):int(r['active_block_index']) for r in params}
    mc_var=np.zeros(60);mc_u=np.zeros(60);mc_u2=np.zeros(60);mc_count=np.zeros(60,dtype=int)
    cell_output=[]
    names=['U','LT','TT']
    for row in cells:
        j=int(row['reco_row']);b=int(row['truth_block'])
        cov=np.array([float(row[f'cov_{a}_{c}']) for a in names for c in names]).reshape(3,3)
        basis=float(row['basis_U']);mc_u[j]+=basis;mc_u2[j]+=cov[0,0];mc_count[j]+=int(row['events'])
        if b in blocks:
            sigma=values[3*blocks[b]:3*blocks[b]+3]
            mc_var[j]+=sigma@cov@sigma
        cell_output.append(dict(reco_row=j,truth_block=b,events=int(row['events']),
                                effective_U_events=basis*basis/cov[0,0] if cov[0,0]>0 else 0))
    mc_var += np.array([float(r['fixed_mc_variance']) for r in reco])
    mc_cov=(operator*mc_var)@operator.T
    eff_y_cov=np.zeros((60,60));lt_y_cov=np.zeros((60,60));eff_details=[]
    charge=sum(float(np.float32(float(r['charge_uC']))) for r in accepted)
    for row in accepted:
        run=RunEvents(args.campaign/'root'/f"signal_events_run{row['run']}.root")
        q=float(np.float32(float(row['charge_uC'])))
        k=float(np.float32(float(row['prescale'])/(float(row['charge_uC'])/1000*float(row['livetime'])*float(row['efficiency']))))*q/charge/.584
        sel=run.selected
        y=np.bincount(run.reco[sel],weights=run.coefficient[sel]*k,minlength=60)
        er=float(row['efficiency_err'])/float(row['efficiency'])
        lr=float(row['livetime_err'])/float(row['livetime'])
        eff_y_cov+=np.outer(y,y)*er**2;lt_y_cov+=np.outer(y,y)*lr**2
        eff_details.append(dict(run=run.run,efficiency_relative_error=er,livetime_relative_error=lr))
    eff_cov=operator@eff_y_cov@operator.T;lt_cov=operator@lt_y_cov@operator.T
    for name,matrix in [('finite_mc_delta',mc_cov),('efficiency_delta',eff_cov),('livetime_delta',lt_cov)]:
        np.savetxt(args.output/f'{name}_sigma_covariance.csv',matrix,delimiter=',',fmt='%.17g')
    out=[]
    for i,row in enumerate(params):
        out.append(dict(parameter_index=i,truth_block=row['truth_block'],region=row['region'],component=row['component'],
                        nominal=values[i],data_bootstrap_sd=data_sd[i],finite_mc_delta_sd=np.sqrt(max(0,mc_cov[i,i])),
                        efficiency_sd=np.sqrt(max(0,eff_cov[i,i])),livetime_sd=np.sqrt(max(0,lt_cov[i,i])),
                        finite_mc_over_data=np.sqrt(max(0,mc_cov[i,i]))/data_sd[i],
                        efficiency_livetime_over_data=np.sqrt(max(0,eff_cov[i,i]+lt_cov[i,i]))/data_sd[i]))
    write_csv(args.output/'sigma_error_components.csv',out)
    write_csv(args.output/'mc_effective_events_by_cell.csv',cell_output)
    write_csv(args.output/'correction_factor_errors.csv',eff_details)
    write_csv(args.output/'mc_effective_events_by_reco.csv',[dict(reco_row=j,events=mc_count[j],
        effective_U_events=mc_u[j]**2/mc_u2[j] if mc_u2[j]>0 else 0,
        data_bootstrap_sd=np.sqrt(np.loadtxt(args.bootstrap/'reco_covariance.csv',delimiter=',')[j,j]),
        finite_mc_yield_sd=np.sqrt(max(0,mc_var[j]))) for j in range(60)])
    published=[r for r in out if r['region']=='published']
    result=dict(mc_response_events=int(mc_count.sum()),
                published_finite_mc_over_data_min=min(r['finite_mc_over_data'] for r in published),
                published_finite_mc_over_data_max=max(r['finite_mc_over_data'] for r in published),
                published_efficiency_livetime_over_data_max=max(r['efficiency_livetime_over_data'] for r in published),
                finite_mc_method='first_order_response_Poisson_moments_fixed_generated_normalization',
                correction_covariance_assumption='independent_runs_and_efficiency_livetime_measurements',
                selected_nominal_combinatorial_background=0)
    (args.output/'statistical_scale_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
