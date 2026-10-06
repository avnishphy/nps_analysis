#!/usr/bin/env python3
"""Production model diagnostics from extractor CSVs; never participates in fitting.

Curves, folded components, pulls and propagated errors come from the extractor.
No model formula or parameter-name assumptions are used here. NaN confidence
errors are omitted, with an explicit annotation; raw Hessian diagnostics remain.
"""
import argparse
import csv
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

TERMS = ('U', 'LT', 'TT')
COLORS = ('#225ea8', '#d7301f', '#238b45')
DISPLAY_SCALE = 1e9  # microbarn/MeV^2 -> nb/GeV^2; no angular factor.


def display(values):
    return np.asarray(values)*DISPLAY_SCALE


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def numbers(rows, key):
    return np.array([float(r[key]) for r in rows])


def table_csv(path, fields, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def points(ax, x, y, error, label, color='black', marker='o'):
    valid = np.isfinite(error) & (error >= 0)
    ax.plot(x, y, linestyle='none', marker=marker, color=color, ms=4, label=label)
    if valid.any():
        ax.errorbar(x[valid], y[valid], yerr=error[valid], fmt='none', color=color, capsize=2)
    if not valid.all():
        ax.text(.02, .97, 'Some confidence errors unavailable (markers only)',
                transform=ax.transAxes, va='top', fontsize=8, color='#9c2d13',
                bbox={'facecolor':'white','edgecolor':'none','alpha':.8})


class Report:
    def __init__(self, args, status):
        self.args, self.status = args, status
        self.directory = args.out / 'model'
        self.directory.mkdir(exist_ok=True)
        self.pages, self.artifacts = [], []
        self.invalid = status['boundary_or_invalid_covariance'] == '1'
        self.context = read(args.out/'model_context.csv')[0]
        # Match common ROOT pages: white background, sans-serif axes, ticks on
        # all sides, compact unboxed legends, black data and colored prediction.
        plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 11,
                             'axes.labelsize': 12, 'legend.fontsize': 9,
                             'legend.frameon': False, 'xtick.top': True,
                             'ytick.right': True, 'xtick.direction': 'in',
                             'ytick.direction': 'in', 'savefig.dpi': 150})

    def save(self, fig, name, title):
        label = 'SigParam2021-inspired pi0 | PROVISIONAL | U slope + 3 normalizations'
        prefix=f'{self.args.kin} | '
        if self.context['synthetic']=='1':prefix=f'SYNTHETIC VALIDATION | {self.args.kin}\n'
        fig.suptitle(f'{prefix}{label}\n{title}', fontsize=11)
        footer = ('Confidence errors unavailable; raw Hessian diagnostic only'
                  if self.invalid else 'Bands: conditional parameter covariance; not independent data errors')
        fig.text(.5, .012, footer, ha='center', fontsize=8,
                 color='#9c2d13' if self.invalid else '0.3')
        fig.tight_layout(rect=(.015, .045, .985, .86))
        for extension, enabled in (('pdf', not self.args.no_pdf), ('png', not self.args.no_png)):
            if enabled:
                path = self.directory / f'{name}.{extension}'
                fig.savefig(path)
                self.artifacts.append(str(path.relative_to(self.args.out)))
                if extension == 'pdf':
                    self.pages.append(path)
        plt.close(fig)


def sigparam_diagnostics(report):
    path=report.args.out/'model_baseline_diagnostics.csv'
    if not path.exists():
        return
    rows=read(path); x=numbers(rows,'abs_t'); tp=numbers(rows,'tprime')
    context=f"Fixed generated mean: Q2={float(rows[0]['Q2']):.3f}, W2={float(rows[0]['W2']):.3f} GeV2, epsilon={float(rows[0]['epsilon']):.3f}"
    ranges={r['quantity']:(float(r['min']),float(r['max'])) for r in read(report.args.out/'model_kinematic_ranges.csv')}
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for sign,color in [('plus','#2166ac'),('minus','#b2182b')]:
        axes[0].plot(x,numbers(rows,'kernel_charged_'+sign),color=color,label='charged '+sign)
        axes[0].plot(x,numbers(rows,'kernel_neutral_'+sign),color=color,linestyle='--',label='neutralized '+sign)
        axes[1].plot(x,display(numbers(rows,'charged_L_'+sign)),color=color,label='charged pole L '+sign)
        axes[1].plot(x,display(numbers(rows,'L_'+sign)),color=color,linestyle='--',label='neutralized L '+sign)
    axes[1].plot(x,numbers(rows,'old_pi0_L'),color='black',linestyle=':',label='original SIMC pi0: L=0')
    for ax in axes:
        ax.axvspan(*ranges['abs_t'],alpha=.08,color='gray',label='accepted |t| range')
        ax.set_xlabel('|t| [GeV2]');ax.legend(fontsize=8);ax.grid(alpha=.2)
    axes[0].set_ylabel('Kernel (numerical GeV convention)');axes[1].set_ylabel('L [nb/GeV2], model assumed')
    report.save(fig,'neutral_longitudinal_kernel','Charged pole removed; smooth neutral L retained\n'+context)
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax,term in zip(axes.flat,('T','L','LT','TT')):
        for suffix,label,style in [('plus','pi+ ingredient','--'),('minus','pi- ingredient',':'),('baseline','component average','-')]:
            ax.plot(tp,display(numbers(rows,term+'_'+suffix)),style,label=label)
        ax.axvspan(*ranges['tprime'],color='gray',alpha=.08,label='Accepted generated range')
        ax.set(xlabel="signed t' [GeV2]",ylabel=f'{term} [nb/GeV2]');ax.legend();ax.grid(alpha=.2)
    report.save(fig,'charged_average_ingredients','Component-wise average; L ingredients use neutral non-pole kernel\n'+context)
    fig,axes=plt.subplots(1,3,figsize=(13,4.8))
    for ax,term in zip(axes,('U','LT','TT')):
        ax.plot(tp,display(numbers(rows,term+'_baseline')),'--',label='Baseline')
        ax.plot(tp,display(numbers(rows,term+'_fitted')),label='Fitted correction')
        ax.axvspan(*ranges['tprime'],color='gray',alpha=.08,label='Accepted generated range')
        ax.set(xlabel="signed t' [GeV2]",ylabel=f'{term} [nb/GeV2]');ax.legend();ax.grid(alpha=.2)
    report.save(fig,'baseline_vs_fitted','Only U gains a pivoted slope; LT and TT retain fixed shapes\n'+context)
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for key,factor,label in [('T_baseline',1.,'T component - model assumed'),('L_baseline',float(rows[0]['epsilon']),'epsilon L - model assumed'),('U_baseline',1.,'U = T + epsilon L baseline')]:
        axes[0].plot(tp,display(numbers(rows,key)*factor),label=label)
    axes[0].set(xlabel="signed t' [GeV2]",ylabel='Structure function [nb/GeV2]');axes[0].legend()
    axes[1].plot(tp,numbers(rows,'f_L'));axes[1].axhline(.5,color='red',linestyle='--',label='L dominance threshold')
    axes[1].set(xlabel="signed t' [GeV2]",ylabel='epsilon L / U');axes[1].legend()
    for ax in axes:ax.axvspan(*ranges['tprime'],color='gray',alpha=.08)
    report.save(fig,'assumed_T_L_decomposition','MODEL-DEPENDENT T/L DECOMPOSITION: no experimental separation\n'+context)
    positivity=read(report.args.out/'model_event_positivity.csv')
    summary_rows=read(report.args.out/'model_longitudinal_summary.csv')[0]
    fig,ax=plt.subplots(figsize=(12,7));ax.axis('off')
    lines=['PROVISIONAL charged-pion-inspired pi0 starting model; expected to be replaced.',
           'G(Q2) is empirical damping, not a pi0 electromagnetic or transition form factor.',
           'Fitted curves: declared fixed context. Bin markers: response-weighted event averages.',
           'Detector predictions: exact generated-event integration; no-model points are diagnostic only.',
           f"Response-weighted epsilon L / U: {float(summary_rows['response_weighted_L_fraction']):.4f}",
           f"Events with L dominating U: {summary_rows['events_L_fraction_above_half']} / {summary_rows['physical_events']}"]
    for r in positivity:
        lines.append(f"{r['stage']}: angular-negative events {r['negative_points']}/{r['points']}; minimum {float(r['minimum_response'])*DISPLAY_SCALE:.5g} nb/GeV2 (before 1/2pi)")
    lines+=['', 'Accepted generated ranges:']+[f'{name}: {lo:.5g} to {hi:.5g}' for name,(lo,hi) in ranges.items()]
    ax.text(.02,.98,'\n'.join(lines),va='top',fontsize=10,linespacing=1.55)
    report.save(fig,'provisional_model_validity','Physics limitations, longitudinal fraction and full-phi positivity')


def summary(report, parameters, correlations, validity):
    status, args = report.status, report.args
    slope = next(p for p in parameters if p['name'] == 'DeltaB_U')
    slope_state = ('fixed at 0 GeV^-2' if int(report.context['U_slope_fixed'])
                   else f"free; fitted DeltaB_U = {float(slope['value']):.9g} GeV^-2")
    solver = ('staged_feasible' if (args.out / 'model_fit_strategy.txt').exists()
              else 'Minuit')
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.axis('off')
    lines = [f'Objective: {args.objective}; model: {status["model_identifier"]}',
             f'Config: {Path(args.config).name}',
             f"Fixed response pivot tau0: {float(report.context['tau0_GeV2']):.9g} GeV2; U slope {slope_state}",
             f'{solver} status: {status["status"]}; covariance status: {status["covariance_status"]}',
             f'EDM: {float(status["edm"]):.5g}; objective: {float(status["objective"]):.6g}',
             f'Rows: {status["rows"]}; parameters: {status["parameters"]}; nominal DOF: {status["nominal_dof"]}',
             f'Objective/DOF: {float(status["objective_per_dof"]):.6g} (descriptive, not a p-value)',
             f'MC iterations: {status["mc_iterations"]}; MC converged: {status["mc_converged"]}',
             f'Negative folded rows: {validity["negative_reconstructed_rows"]}; angular-negative blocks: {validity["negative_angular_blocks"]}',
             f'Model parameters: {sum(p["role"]=="model" for p in parameters)}; nuisance parameters: {sum(p["role"]!="model" for p in parameters)}',
             f'Coordinate bound hits: {", ".join(p["name"] for p in parameters if p["at_bound"]=="1") or "none"} (angular constraints separate)',
             f'Raw correlation condition: {float(status["correlation_condition"]):.5g}',
             f'Strong raw correlations (|rho| > 0.9): {len(correlations)}',
             'Global QA shape overlays are area normalized; detector pages use absolute folded yields.',
             'Selected migration fractions describe selected MC; they are not acceptance efficiencies.',
             'Full input paths and invocation: pipeline log and pipeline_artifacts.json.']
    ax.text(.02, .97, '\n\n'.join(lines), va='top', fontsize=10)
    report.save(fig, 'run_summary', 'Run and fit quality')


def input_statistics(report, rows, matching):
    fig,ax=plt.subplots(figsize=(11,8)); ax.axis('off')
    lines=['Response events are selected matched exclusive SIMC, not the generated denominator.',
           f'Truth blocks: {len(rows)}; active: {sum(int(r["active_block_index"])>=0 for r in rows)}',
           f'Selected response events: {sum(int(r["events"]) for r in rows)}',
           'Raw/smeared event matching (before rectangular bin cuts):']
    for key,value in matching.items(): lines.append(f'  {key}: {value}')
    ax.text(.02,.95,'\n\n'.join(lines),va='top',fontsize=10)
    report.save(fig,'input_statistics','Input statistics and event matching')


def structure(report, curves, reported, reference):
    for term, old_term, color in zip(TERMS, ('U', 'TL', 'TT'), COLORS):
        rows = [r for r in curves if r['component'] == term]
        x = numbers(rows,'tprime'); y,err=(display(numbers(rows,k)) for k in ('value','error'))
        for compare in (False, True) if reference else (False,):
            fig, ax = plt.subplots(figsize=(9, 6.5))
            points(ax, numbers(reported, 'tprime'), display(numbers(reported, f'sigma_{term}')),
                   display(numbers(reported, f'sigma_{term}_error')), 'Event-averaged fitted model')
            ax.plot(x, y, '--', lw=.9, color=color, label='Model at reference Q2,W2,epsilon - diagnostic only')
            ranges={r['quantity']:(float(r['min']),float(r['max'])) for r in read(report.args.out/'model_kinematic_ranges.csv')}
            ax.axvspan(*ranges['tprime'],color='gray',alpha=.08,label='Accepted generated range; outside is extrapolation')
            if compare:
                rows = [r for r in reference if r['fit_xsec_ok'] == '1']
                points(ax, numbers(rows, 'mean_tprime_vertex_sim'), display(numbers(rows, f'fit_xsec_sigma{old_term}')),
                       display(numbers(rows, f'fit_xsec_sigma{old_term}err')), 'No-model independent extraction', '#8856a7', 's')
            ax.set(xlabel="Signed generated t' [GeV$^2$]",
                   ylabel=rf'$\sigma_{{{term}}}$ [nb/GeV$^2$]')
            ax.axhline(0, color='0.7', lw=.7); ax.legend(); ax.grid(alpha=.15)
            report.save(fig, ('comparison_' if compare else 'sigma_')+term,
                        f'Born-level structure function: {term}'+(' | diagnostic comparison' if compare else ''))


def detector(report, rows):
    groups = [('all_rows', rows)]
    keys = sorted({tuple(int(r[k]) for k in ('it', 'iq', 'ix')) for r in rows})
    groups += [(f'phi_t{t}_q{q}_x{x}', [r for r in rows if tuple(int(r[k]) for k in ('it','iq','ix'))==(t,q,x)])
               for t,q,x in keys]
    for name, group in groups:
        x = numbers(group, 'row') if name=='all_rows' else (numbers(group,'phi_lo')+numbers(group,'phi_hi'))*90/np.pi
        y, mu = numbers(group,'data'), numbers(group,'prediction')
        err = numbers(group,'prediction_error')
        fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, gridspec_kw={'height_ratios': [2,1,1]})
        points(axes[0], x, y, np.sqrt(numbers(group,'data_sumw2')), 'Data')
        axes[0].plot(x, mu, color='red', label='Forward-folded prediction')
        valid = np.isfinite(err)&(err>=0)
        if valid.any() and not report.invalid:
            axes[0].fill_between(x, mu-err, mu+err, where=valid, color='red', alpha=.15)
        axes[0].set_ylabel('Yield / mC'); axes[0].legend()
        axes[1].plot(x, numbers(group,'residual'), 'o', ms=4, color='black')
        axes[1].set_ylabel('Data - prediction\n[yield / mC]')
        pull = numbers(group,'objective_pull'); included = numbers(group,'included')==1
        axes[2].plot(x[included], pull[included], 'o', ms=4, color='black')
        axes[2].set_ylabel('Signed sqrt(deviance)' if report.args.objective=='scaled-poisson' else 'Objective pull')
        axes[2].set_xlabel('Reconstructed row' if name=='all_rows' else r'Reconstructed $\phi$ [degrees]')
        for ax in axes[1:]: ax.axhline(0, color='0.5', lw=.8)
        if (~included).any():
            axes[2].text(.02,.92,f'{sum(~included)} excluded rows: no objective pull',transform=axes[2].transAxes,va='top',fontsize=9)
        report.save(fig, 'detector_'+name, 'Detector-level yield, residual and objective pull | '+name)
        if name=='all_rows': continue
        fig, ax = plt.subplots(figsize=(10,6.5))
        ax.plot(x,mu,color='black',label='Full forward-folded prediction')
        for term,color in zip(TERMS,COLORS):
            ax.plot(x,numbers(group,'component_'+term),'o-',ms=3,color=color,label=term+' contribution')
        ax.axhline(0,color='0.5',lw=.8); ax.legend()
        ax.set(xlabel=r'Reconstructed $\phi$ [degrees]',ylabel='Signed yield / mC')
        report.save(fig,'components_'+name,'Detector-level signed components (including nuisances) | '+name)


def coverage(report, rows):
    fig, axes = plt.subplots(2,1,figsize=(11,8))
    physical_ids={r['truth_block'] for r in read(report.args.out/'model_structure_functions.csv')}
    for role,color,marker in (('physical','#225ea8','o'),('nuisance','#d7301f','s'),('unsupported','0.5','x')):
        group = [r for r in rows if ('unsupported' if int(r['active_block_index'])<0 else
                  'physical' if r['truth_block'] in physical_ids else 'nuisance')==role]
        for ax,key in zip(axes,('events','response_weight')):
            valid = [r for r in group if np.isfinite(float(r['mean_tprime'])) and float(r['response_weight'])>0]
            label='exterior / excluded-group nuisance' if role=='nuisance' else role
            ax.scatter(-numbers(valid,'mean_tprime'),numbers(valid,key),label=label,marker=marker,color=color)
            for r in valid: ax.annotate('b'+r['truth_block'],(-float(r['mean_tprime']),float(r[key])),fontsize=8)
            ax.set_ylabel('Selected events' if key=='events' else 'Integrated response strength')
            ax.legend()
    for ax,key in zip(axes,('events','response_weight')):
        values=numbers(rows,key)
        positive=values[np.isfinite(values)&(values>0)]
        if positive.size and positive.max()/positive.min()>100:
            ax.set_yscale('log')
            ax.set_ylim(positive.min()*.5,positive.max()*3)
        else:
            # Equal fixture weights (and narrow real coverage) must not turn
            # floating-point roundoff into apparent response variation.
            ax.set_ylim(0,positive.max()*1.4 if positive.size else 1)
    absent = [r['truth_block'] for r in rows if int(r['active_block_index'])<0]
    axes[0].text(.02,.96,'Unsupported blocks (no fit response): '+(', '.join(absent) or 'none'),
                 transform=axes[0].transAxes,va='top',fontsize=9)
    axes[1].set_xlabel(r'Generated $\tau=-t\prime$ [GeV$^2$] (response-weighted block mean)')
    report.save(fig,'coverage','Generated coverage constraining the fit; empty blocks have no measured mean')


def parameters_and_correlations(report, parameters, covrows):
    n=len(parameters); matrix=np.full((n,n),np.nan)
    for row in covrows: matrix[int(row['i']),int(row['j'])]=float(row['correlation'])
    model=[i for i,p in enumerate(parameters) if p['role']=='model']
    for name,indices in (('full',list(range(n))),('model',model)):
        fig, ax = plt.subplots(figsize=(max(9,min(16,.45*len(indices)+3)),max(8,min(15,.45*len(indices)+2))))
        im=ax.imshow(matrix[np.ix_(indices,indices)],vmin=-1,vmax=1,cmap='coolwarm')
        labels=[('M: ' if parameters[i]['role']=='model' else 'N: ')+parameters[i]['name'] for i in indices]
        ax.set_xticks(range(len(indices)),labels,rotation=65,ha='right',fontsize=9)
        ax.set_yticks(range(len(indices)),labels,fontsize=9)
        fig.colorbar(im,ax=ax,label='Raw Hessian correlation')
        if name=='model':
            for a,i in enumerate(indices):
                for b,j in enumerate(indices): ax.text(b,a,f'{matrix[i,j]:.2f}',ha='center',va='center',fontsize=9)
        report.save(fig,'correlation_'+name,'Parameter correlation | M: physics model; N: migration nuisance')
    for start in range(0,n,18):
        fig,ax=plt.subplots(figsize=(12,8)); ax.axis('off')
        rows=[]
        for p in parameters[start:start+18]:
            rows.append([p['name'],'model' if p['role']=='model' else 'nuisance',f"{float(p['value']):.5g}",
                         f"{float(p['error']):.3g}" if np.isfinite(float(p['error'])) else 'unavailable',
                         f"{float(p['raw_hessian_error']):.3g}",f"{float(p['lower_bound']):.3g}",f"{float(p['upper_bound']):.3g}",p['at_bound'],
                         f"{float(p['relative_error']):.3g}" if np.isfinite(float(p['relative_error'])) else 'unavailable'])
        table=ax.table(cellText=rows,colLabels=['Parameter','Role','Value','Error','Raw Hesse','Lower','Upper','Bound','Rel. error'],loc='center',cellLoc='left')
        table.auto_set_font_size(False); table.set_fontsize(8); table.scale(1,1.65)
        report.save(fig,f'parameters_{start//18}','Parameter estimates, bounds and identifiability')
    pairs=sorted(((abs(matrix[i,j]),i,j) for i in range(n) for j in range(i+1,n)
                  if np.isfinite(matrix[i,j])),reverse=True)[:12]
    fig,ax=plt.subplots(figsize=(11,8)); ax.axis('off')
    lines=['Largest absolute raw Hessian correlations (diagnostic):','']
    for _,i,j in pairs: lines.append(f'{parameters[i]["name"]} / {parameters[j]["name"]}: {matrix[i,j]:+.5f}')
    lines+=['','Weak parameters (raw Hesse width >= |value|):',
            ', '.join(p['name'] for p in parameters if p['poorly_constrained']=='1') or 'none']
    ax.text(.02,.95,'\n\n'.join(lines),va='top',fontsize=10)
    report.save(fig,'identifiability','Strongest correlations and weak parameter directions')


def diagnostic_ratio(numerator,error,denominator,scale):
    """Denominator held fixed: no claim of independent same-data fit errors."""
    if not np.isfinite([numerator,denominator,scale]).all():
        return np.nan,np.nan,'nonfinite'
    if abs(denominator)<=max(1e-30,1e-3*abs(scale)):
        return np.nan,np.nan,'near_zero_denominator'
    return numerator/denominator,abs(error/denominator) if np.isfinite(error) else np.nan,'ok'


def shape_diagnostics(report,reference,render):
    out=report.args.out
    rows=read(out/'model_reconstructed_yields.csv')
    config=json.loads(Path(report.args.config).read_text())
    edges=config['tprime_bin_edges']
    keys=sorted({tuple(int(r[k]) for k in ('it','iq','ix')) for r in rows})
    table=[]
    for it,iq,ix in keys:
        selected=[r for r in rows if tuple(int(r[k]) for k in ('it','iq','ix'))==(it,iq,ix)]
        included=[r for r in selected if r['included']=='1']
        pull=numbers(included,'objective_pull')
        table.append(dict(it=it,iq=iq,ix=ix,tprime_lo=edges[it],tprime_hi=edges[it+1],tprime=(edges[it]+edges[it+1])/2,
            included_rows=len(included),data=float(sum(numbers(selected,'data'))),prediction=float(sum(numbers(selected,'prediction'))),
            residual=float(sum(numbers(selected,'residual'))),mean_pull=float(np.mean(pull)) if len(pull) else np.nan,
            objective_contribution=float(np.sum(pull**2)) if len(pull) else np.nan))
    table_csv(out/'model_tprime_shape.csv',tuple(table[0]),table);report.artifacts.append('model_tprime_shape.csv')
    if render:
        for iq,ix in sorted({(k[1],k[2]) for k in keys}):
            group=[r for r in table if r['iq']==iq and r['ix']==ix];x=numbers(group,'tprime')
            fig,axes=plt.subplots(3,1,figsize=(10,8),sharex=True)
            axes[0].plot(x,numbers(group,'residual'),'ko-');axes[0].set_ylabel('Summed data - prediction\n[yield / mC]')
            axes[1].plot(x,numbers(group,'mean_pull'),'ko-');axes[1].set_ylabel('Mean objective pull')
            axes[2].plot(x,numbers(group,'objective_contribution'),'ko-');axes[2].set_ylabel('Objective contribution')
            axes[2].set_xlabel("Reconstructed t' bin center [GeV2]")
            for ax in axes:ax.axhline(0,color='.6',lw=.7);ax.grid(alpha=.2)
            report.save(fig,f'shape_tprime_q{iq}_x{ix}','Detector diagnostics | sums by reconstructed t\' bin\nDescriptive contributions; correlated pulls have no assigned binwise significance')
    spread=read(out/'model_charged_spread.csv');summary=[]
    if render:fig,axes=plt.subplots(1,3,figsize=(13,5))
    for i,term in enumerate(TERMS):
        valid=[r for r in spread if r['component']==term and r['status']=='ok'];d=numbers(valid,'D')
        summary.append(dict(component=term,valid_points=len(valid),minimum=float(min(d)) if len(d) else np.nan,
                            maximum=float(max(d)) if len(d) else np.nan,median=float(np.median(d)) if len(d) else np.nan))
        if render:
            axes[i].scatter(numbers(valid,'tprime'),d,s=3,alpha=.15,rasterized=True)
            axes[i].set(xlabel="Generated t' [GeV2]",ylabel=f'D({term})');axes[i].axhline(0,color='.6',lw=.7)
    table_csv(out/'model_charged_spread_summary.csv',tuple(summary[0]),summary);report.artifacts.append('model_charged_spread_summary.csv')
    if render:report.save(fig,'charged_dimensionless_spread','D = (plus-minus) / [0.5 (|plus|+|minus|)]\nDiagnostic of the charged-average assumption; NOT an uncertainty band')
    if not reference:return
    ref={tuple(int(r[k]) for k in ('it','iq','ix')):r for r in reference if r['fit_xsec_ok']=='1'}
    model=read(out/'model_structure_functions.csv');ratios=[];summaries=[]
    for term,old in zip(TERMS,('U','TL','TT')):
        scale=max(abs(float(r['sigma_'+term])) for r in model)
        for r in model:
            key=tuple(int(r[k]) for k in ('it','iq','ix'))
            if key not in ref:continue
            reference_row=ref[key];num=float(reference_row['fit_xsec_sigma'+old]);error=float(reference_row['fit_xsec_sigma'+old+'err'])
            den=float(r['sigma_'+term]);ratio,ratio_error,status=diagnostic_ratio(num,error,den,scale)
            ratios.append(dict(component=term,truth_block=r['truth_block'],it=key[0],iq=key[1],ix=key[2],tprime=float(r['tprime']),
                no_model=num,no_model_error=error,model=den,model_error=float(r['sigma_'+term+'_error']),ratio=ratio,
                ratio_error_numerator_only=ratio_error,status=status,opposite_sign=int(num*den<0)))
        valid=[r for r in ratios if r['component']==term and r['status']=='ok'];v=numbers(valid,'ratio')
        summaries.append(dict(component=term,valid_bins=len(valid),skipped_bins=sum(r['component']==term and r['status']!='ok' for r in ratios),
            ratio_min=float(min(v)) if len(v) else np.nan,ratio_max=float(max(v)) if len(v) else np.nan,
            max_absolute_departure_from_unity=float(max(abs(v-1))) if len(v) else np.nan,
            opposite_sign_bins=sum(r['opposite_sign'] for r in valid)))
    if not ratios:return
    table_csv(out/'model_shape_ratios.csv',tuple(ratios[0]),ratios);table_csv(out/'model_shape_summary.csv',tuple(summaries[0]),summaries)
    report.artifacts.extend(('model_shape_ratios.csv','model_shape_summary.csv'))
    if render:
        for iq,ix in sorted({(r['iq'],r['ix']) for r in ratios}):
            fig,axes=plt.subplots(1,3,figsize=(13,5))
            for ax,term in zip(axes,TERMS):
                group=[r for r in ratios if r['component']==term and r['iq']==iq and r['ix']==ix];valid=[r for r in group if r['status']=='ok']
                points(ax,numbers(valid,'tprime'),numbers(valid,'ratio'),numbers(valid,'ratio_error_numerator_only'),'No-model / normalized baseline')
                ax.axhline(1,color='.5',linestyle='--');ax.set(xlabel="Signed t' [GeV2]",ylabel=f'R({term})')
                ax.text(.03,.03,f'{len(group)-len(valid)} unstable ratios skipped',transform=ax.transAxes,fontsize=8)
            report.save(fig,f'shape_ratios_q{iq}_x{ix}','Shape-only diagnostic: model denominator held fixed; ratios are NOT fitted\nBars: numerator errors only; shared-data/model correlation is not evaluated')


def u_correction(report,parameters):
    pars={p['name']:p for p in parameters}
    slope=float(pars['DeltaB_U']['value']);error=float(pars['DeltaB_U']['error'])
    pivot=float(report.context['tau0_GeV2'])
    interval=next(r for r in read(report.args.out/'model_kinematic_ranges.csv') if r['quantity']=='tau')
    tau=np.linspace(float(interval['min']),float(interval['max']),201)
    correction=np.exp(-slope*(tau-pivot))
    rows=[dict(tau=t,tau0=pivot,correction=c,error=abs((t-pivot)*c*error)) for t,c in zip(tau,correction)]
    table_csv(report.args.out/'model_U_correction.csv',tuple(rows[0]),rows)
    report.artifacts.append('model_U_correction.csv')
    if report.args.no_pdf and report.args.no_png:return
    fig,ax=plt.subplots(figsize=(9,6))
    ax.plot(tau,correction,label=f'DeltaB_U = {slope:.4g} +/- {error:.3g} GeV^-2')
    if np.isfinite(error):
        err=numbers(rows,'error');ax.fill_between(tau,correction-err,correction+err,alpha=.15)
    ax.axhline(1,color='.5',ls='--');ax.axvline(pivot,color='.5',ls=':',label=f'Fixed pivot = {pivot:.4g} GeV2')
    ax.set(xlabel="Accepted tau = -t' [GeV2]",ylabel='C_U(tau) = exp[-DeltaB_U (tau-tau0)]')
    ax.legend();ax.grid(alpha=.2)
    report.save(fig,'U_slope_correction','Empirical U shape correction only; normalization excluded')


def before_after(report):
    before=report.args.before
    if before is None:return
    after=report.args.out
    old_rows=read(before/'model_reconstructed_yields.csv');new_rows=read(after/'model_reconstructed_yields.csv')
    keys=('row','included','data','data_sumw2','it','iq','ix','ip')
    if len(old_rows)!=len(new_rows) or any(any(a[k]!=b[k] for k in keys) for a,b in zip(old_rows,new_rows)):
        raise ValueError('Before/after comparison requires identical data rows and masks')
    for filename in ('model_event_cache.csv','migration_response_cells.csv'):
        if (before/filename).read_bytes()!=(after/filename).read_bytes():
            raise ValueError('Before/after response differs: '+filename)
    key=lambda r:tuple(int(r[k]) for k in ('it','iq','ix'))
    shapes=[{key(r):r for r in read(p/'model_tprime_shape.csv')} for p in (before,after)]
    ratios=[{key(r):r for r in read(p/'model_shape_ratios.csv') if r['component']=='U' and r['status']=='ok'} for p in (before,after)]
    table=[]
    for k in sorted(shapes[0].keys() & shapes[1].keys()):
        r=dict(it=k[0],iq=k[1],ix=k[2],tprime=float(shapes[0][k]['tprime']),tau=-float(shapes[0][k]['tprime']))
        for i,label in enumerate(('before','after')):
            r['residual_'+label]=float(shapes[i][k]['residual']);r['mean_pull_'+label]=float(shapes[i][k]['mean_pull'])
            r['U_ratio_'+label]=float(ratios[i][k]['ratio']) if k in ratios[i] else np.nan
            r['tau_generated_'+label]=-float(ratios[i][k]['tprime']) if k in ratios[i] else np.nan
        table.append(r)
    table_csv(after/'model_shape_comparison.csv',tuple(table[0]),table);report.artifacts.append('model_shape_comparison.csv')
    summary={'same_data_and_response':True,'before':str(before),'groups':[]}
    for label,path,data in zip(('before','after'),(before,after),(old_rows,new_rows)):
        status=read(path/'model_fit_status.csv')[0]
        summary[label+'_fit']=status
        summary[label+'_max_abs_pull']=float(max(abs(numbers(data,'objective_pull'))))
    summary['objective_improvement']=float(summary['before_fit']['objective'])-float(summary['after_fit']['objective'])
    for iq,ix in sorted({(r['iq'],r['ix']) for r in table}):
        group=[r for r in table if (r['iq'],r['ix'])==(iq,ix)]
        item=dict(iq=iq,ix=ix)
        render=not(report.args.no_pdf and report.args.no_png)
        if render:fig,axes=plt.subplots(1,2,figsize=(12,6))
        for label,color in (('before','#777777'),('after','#2166ac')):
            tau=numbers(group,'tau');res=numbers(group,'residual_'+label)
            rt=numbers(group,'tau_generated_'+label);ratio=numbers(group,'U_ratio_'+label)
            good=np.isfinite(rt)&np.isfinite(ratio)
            item['residual_slope_'+label]=float(np.polyfit(tau,res,1)[0]) if len(tau)>1 else None
            item['U_ratio_slope_'+label]=float(np.polyfit(rt[good],ratio[good],1)[0]) if sum(good)>1 else None
            if render:
                caption='3 parameters' if label=='before' else '4 parameters'
                axes[0].plot(tau,res,'o-',color=color,label=caption)
                axes[1].plot(rt[good],ratio[good],'o-',color=color,label=caption)
        summary['groups'].append(item)
        if render:
            axes[0].axhline(0,color='.5',ls='--');axes[1].axhline(1,color='.5',ls='--')
            axes[0].set(xlabel='Reconstructed tau bin center [GeV2]',ylabel='Summed data - prediction [yield/mC]')
            axes[1].set(xlabel='Response-weighted generated tau [GeV2]',ylabel='U no-model / event-averaged model')
            for ax in axes:ax.legend();ax.grid(alpha=.2)
            report.save(fig,f'U_shape_before_after_q{iq}_x{ix}','Before/after U slope correction | descriptive trends\nSame data and response; ratios not fitted; no independent binwise significance')
    (after/'model_before_after.json').write_text(json.dumps(summary,indent=2)+'\n')
    report.artifacts.append('model_before_after.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--combined',type=Path)
    parser.add_argument('--reference',type=Path)
    parser.add_argument('--before',type=Path,help='Prior three-parameter output; diagnostic only')
    parser.add_argument('--kin',required=True)
    parser.add_argument('--config',required=True)
    parser.add_argument('--objective',required=True)
    parser.add_argument('--no-pdf',action='store_true')
    parser.add_argument('--no-png',action='store_true')
    parser.add_argument('--no-diagnostics',action='store_true')
    args=parser.parse_args()
    out=args.out
    parameters=read(out/'model_parameters.csv'); covrows=read(out/'model_covariance.csv')
    strong=[]
    for r in covrows:
        i,j=int(r['i']),int(r['j'])
        if i<j and abs(float(r['correlation']))>.9:
            strong.append({'parameter_i':parameters[i]['name'],'parameter_j':parameters[j]['name'],
                           'correlation':r['correlation'],'interpretation':'raw_hessian_diagnostic'})
    table_csv(out/'model_strong_correlations.csv',('parameter_i','parameter_j','correlation','interpretation'),strong)
    fields=('name','role','value','error','relative_error','at_bound','poorly_constrained')
    table_csv(out/'model_identifiability.csv',fields,[{k:p[k] for k in fields} for p in parameters])
    status=read(out/'model_fit_status.csv')[0]
    report=Report(args,status)
    # Numerical event-model products are required even with plots disabled.
    report.artifacts.extend(('model_baseline_diagnostics.csv','model_kinematic_ranges.csv',
        'model_longitudinal_summary.csv','model_event_positivity.csv','model_event_cache.csv','model_row_jacobian.csv','model_charged_spread.csv','model_context.csv'))
    reference=[]
    if args.reference and args.reference.is_file():
        reference=read(args.reference)
        print(f'[plots] Diagnostic no-model overlay: {args.reference}')
    else:
        print('[plots] No corresponding no-model output; model-only plots retained')
    if not(args.no_pdf and args.no_png):
        summary(report,parameters,strong,read(out/'model_prediction_validity.csv')[0])
        if not args.no_diagnostics:
            truth=read(out/'migration_truth_blocks.csv')
            matching=out/'model_vertex_matching.csv'
            input_statistics(report,truth,read(matching)[0] if matching.is_file() else {'availability':'not in this replay'})
            coverage(report,truth)
            detector(report,read(out/'model_reconstructed_yields.csv'))
        structure(report,read(out/'model_structure_curves.csv'),read(out/'model_structure_functions.csv'),reference)
        if not args.no_diagnostics:
            parameters_and_correlations(report,parameters,covrows)
            sigparam_diagnostics(report)
    shape_diagnostics(report,reference,render=not(args.no_pdf and args.no_png) and not args.no_diagnostics)
    u_correction(report,parameters)
    before_after(report)
    if not args.no_pdf:
        if not args.combined or not args.combined.is_file():
            parser.error('The extractor combined PDF is required for production assembly')
        # Insert run summary first, preserve all common ROOT pages, append the
        # detailed model diagnostics. Atomic replacement protects the PDF.
        temporary=args.combined.with_suffix('.model-merge.pdf')
        subprocess.run(['pdfunite',str(report.pages[0]),str(args.combined),
                        *(str(p) for p in report.pages[1:]),str(temporary)],check=True)
        temporary.replace(args.combined)
    manifest={'schema_version':1,'kin':args.kin,'mode':'simc_model','artifacts':report.artifacts,
              'reference_overlay':str(args.reference) if reference else None,
              'confidence_errors_valid':not report.invalid,'diagnostics_enabled':not args.no_diagnostics}
    (out/'model_plot_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'[plots] {len(report.pages)} model report pages; confidence errors valid={not report.invalid}')


if __name__=='__main__':
    main()
