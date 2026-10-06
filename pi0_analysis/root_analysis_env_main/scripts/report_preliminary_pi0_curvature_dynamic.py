#!/usr/bin/env python3
"""Local Fisher/profile/toy statistical products for a fresh M0 campaign."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.linalg import null_space
from scipy.optimize import minimize

PHYS = ["N_U", "DeltaB_U", "N_LT", "N_TT"]
COMPONENTS = ("U", "LT", "TT")


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    temporary.replace(path)


def write_matrix(path, matrix, names):
    write_rows(path, [{"parameter": row_name, **{name: float(value) for name, value in zip(names, row)}}
                      for row_name, row in zip(names, matrix)])


def write_json(path, payload):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n");temporary.replace(path)


def correlation(matrix):
    scale = np.sqrt(np.maximum(np.diag(matrix), 0.0));denominator = np.outer(scale, scale)
    return np.divide(matrix, denominator, out=np.zeros_like(matrix), where=denominator > 0)


def angular_minimum(u, lt, tt, epsilon):
    linear = np.sqrt(2 * epsilon * (1 + epsilon)) * lt;quadratic = 2 * epsilon * tt
    z = np.where(quadratic > 0,
                 np.clip(np.divide(-linear, 2 * quadratic, out=np.zeros_like(linear), where=quadratic > 0), -1, 1),
                 np.where(linear > 0, -1.0, 1.0))
    return u + linear * z + epsilon * tt * (2 * z * z - 1), z


class Local:
    def __init__(self, source):
        self.source = source
        self.central = json.loads((source / "central_results.json").read_text())["additive"]
        self.parameters = np.asarray(self.central["theta"], dtype=float)
        self.y = np.asarray(self.central["y"], dtype=float);self.variance = np.asarray(self.central["variance"], dtype=float)
        self.pivot = float(self.central["info"][6])
        self.config = json.loads((source / "config_snapshot.json").read_text())
        context = json.loads((source / "campaign_context.json").read_text());self.fit_output = Path(context["fit_output_dir"])
        events = np.genfromtxt(self.fit_output / "model_event_cache.csv", delimiter=",", names=True)
        if events.ndim == 0:events = np.array([events], dtype=events.dtype)
        self.events = events;reconstructed = read_rows(self.fit_output / "model_reconstructed_yields.csv")
        self.mask = np.array([int(row["included"]) != 0 for row in reconstructed]);self.nrows = len(reconstructed)
        self.row = events["row"].astype(int);self.physical = events["physical"].astype(bool)
        regions={int(row["truth_block"]):row["region"] for row in read_rows(self.fit_output / "migration_truth_blocks.csv")}
        self.physical_blocks = sorted(set(events["truth_block"][self.physical].astype(int)))
        exterior=sorted(set(events["truth_block"][~self.physical].astype(int)))
        self.nuisance_blocks = [block for block in exterior if regions[block] == "tprime_below"]
        self.fixed_blocks = [block for block in exterior if regions[block] != "tprime_below"]
        self.nbins = len(self.physical_blocks);self.nparams = 4 + 3 * len(self.nuisance_blocks)
        if len(self.parameters) != self.nparams:raise ValueError("central parameter vector does not match active migration blocks")
        if self.nbins != len(self.config["tprime_bin_edges"]) - 1:raise ValueError("physical blocks do not match config t-prime bins")
        self.names = [f"{component}(bin{index})" for component in COMPONENTS for index in range(self.nbins)]
        self.nuisance_names = [f"feedin_tprime_below_{component}" for _ in self.nuisance_blocks for component in COMPONENTS]
        self.physical_rows = self.row[self.physical];self.delta_tau = events["tau"][self.physical] - self.pivot
        self.baseline_u = events["baseline_U"][self.physical];self.baseline_lt = events["baseline_LT"][self.physical]
        self.baseline_tt = events["baseline_TT"][self.physical];self.epsilon = events["epsilon"][self.physical]
        self.epsilon_factor = np.sqrt(2 * self.epsilon * (1 + self.epsilon));self.design = np.zeros((self.nrows, self.nparams))
        for index, component, baseline in ((2, "LT", self.baseline_lt), (3, "TT", self.baseline_tt)):
            self.design[:, index] = np.bincount(self.physical_rows,
                weights=events[f"basis_{component}"][self.physical] * baseline, minlength=self.nrows)
        self.nuisance_epsilon = []
        for block_index, block in enumerate(self.nuisance_blocks):
            take = events["truth_block"] == block;self.nuisance_epsilon.append(float(events["epsilon"][take].max()))
            for component_index, component in enumerate(COMPONENTS):
                self.design[:, 4 + 3 * block_index + component_index] = np.bincount(
                    self.row[take], weights=events[f"basis_{component}"][take], minlength=self.nrows)
        fixed=np.isin(events["truth_block"].astype(int),self.fixed_blocks)
        fixed_values=np.column_stack([events[f"baseline_{component}"][fixed] for component in COMPONENTS])
        fixed_basis=np.column_stack([events[f"basis_{component}"][fixed] for component in COMPONENTS])
        self.fixed_prediction=np.bincount(self.row[fixed],weights=np.sum(fixed_basis*fixed_values,axis=1),minlength=self.nrows)

    def prediction_jacobian(self, parameters):
        jacobian = self.design.copy();term = (self.events["basis_U"][self.physical] * self.baseline_u *
                np.exp(-parameters[1] * self.delta_tau))
        jacobian[:, 0] = np.bincount(self.physical_rows, weights=term, minlength=self.nrows)
        prediction = self.fixed_prediction + self.design @ parameters + parameters[0] * jacobian[:, 0]
        jacobian[:, 1] = -parameters[0] * np.bincount(self.physical_rows, weights=term*self.delta_tau, minlength=self.nrows)
        return prediction, jacobian

    def objective_gradient(self, parameters):
        prediction, jacobian = self.prediction_jacobian(parameters);residual = (prediction-self.y)[self.mask]
        return float(np.sum(residual*residual/self.variance[self.mask])), 2*(residual/self.variance[self.mask])@jacobian[self.mask]

    def physics_margin(self, parameters):
        u = parameters[0]*self.baseline_u*np.exp(-parameters[1]*self.delta_tau)
        margin,z = angular_minimum(u,parameters[2]*self.baseline_lt,parameters[3]*self.baseline_tt,self.epsilon)
        return margin,z,u

    def nuisance_margins(self, parameters):
        result=[]
        for block_index,epsilon in enumerate(self.nuisance_epsilon):
            index=4+3*block_index;margin,_=angular_minimum(np.array([parameters[index]]),
                np.array([parameters[index+1]]),np.array([parameters[index+2]]),np.array([epsilon]));result.append(float(margin[0]))
        return result

    def published(self, parameters):
        result=np.zeros(3*self.nbins)
        for index,block in enumerate(self.physical_blocks):
            take=self.events["truth_block"]==block;weights=self.events["response_weight"][take]
            delta=self.events["tau"][take]-self.pivot
            result[index]=np.average(parameters[0]*self.events["baseline_U"][take]*np.exp(-parameters[1]*delta),weights=weights)*1e9
            result[self.nbins+index]=np.average(parameters[2]*self.events["baseline_LT"][take],weights=weights)*1e9
            result[2*self.nbins+index]=np.average(parameters[3]*self.events["baseline_TT"][take],weights=weights)*1e9
        return result

    def published_jacobian(self, parameters):
        result=np.zeros((3*self.nbins,4))
        for index,block in enumerate(self.physical_blocks):
            take=self.events["truth_block"]==block;weights=self.events["response_weight"][take]
            delta=self.events["tau"][take]-self.pivot;u=parameters[0]*self.events["baseline_U"][take]*np.exp(-parameters[1]*delta)
            result[index,:2]=[np.average(u/parameters[0],weights=weights)*1e9,np.average(-delta*u,weights=weights)*1e9]
            result[self.nbins+index,2]=np.average(self.events["baseline_LT"][take],weights=weights)*1e9
            result[2*self.nbins+index,3]=np.average(self.events["baseline_TT"][take],weights=weights)*1e9
        return result

    def profile(self, fixed_index, fixed_value, start, scale):
        free=np.array([index for index in range(self.nparams) if index!=fixed_index]);initial=start[free]/scale[free]
        active_ids=set(np.argsort(self.physics_margin(start)[0])[:64].tolist())
        def unpack(values):
            parameters=np.empty(self.nparams);parameters[fixed_index]=fixed_value;parameters[free]=values*scale[free];return parameters
        def objective(values):
            value,gradient=self.objective_gradient(unpack(values));return value,gradient[free]*scale[free]
        def constraints(values):
            parameters=unpack(values);ids=np.array(sorted(active_ids));u=parameters[0]*self.baseline_u[ids]*np.exp(-parameters[1]*self.delta_tau[ids])
            lt,tt=parameters[2]*self.baseline_lt[ids],parameters[3]*self.baseline_tt[ids];values_out=[];jacobians=[]
            for z in (-1.0,1.0):
                values_out.append(u+self.epsilon_factor[ids]*lt*z+self.epsilon[ids]*tt);jacobian=np.zeros((len(ids),self.nparams))
                jacobian[:,:4]=np.c_[u/parameters[0],-self.delta_tau[ids]*u,self.epsilon_factor[ids]*z*self.baseline_lt[ids],self.epsilon[ids]*self.baseline_tt[ids]];jacobians.append(jacobian)
            margin,z=angular_minimum(u,lt,tt,self.epsilon[ids]);interior=(tt>0)&(abs(z)<1);values_out.append(np.where(interior,margin,1e-8))
            jacobian=np.zeros((len(ids),self.nparams));jacobian[:,:4]=np.c_[u/parameters[0],-self.delta_tau[ids]*u,
                self.epsilon_factor[ids]*z*self.baseline_lt[ids],self.epsilon[ids]*(2*z*z-1)*self.baseline_tt[ids]]*interior[:,None];jacobians.append(jacobian)
            for block_index,epsilon in enumerate(self.nuisance_epsilon):
                index=4+3*block_index;margin,z=angular_minimum(np.array([parameters[index]]),np.array([parameters[index+1]]),np.array([parameters[index+2]]),np.array([epsilon]))
                values_out.append(margin);jacobian=np.zeros((1,self.nparams));jacobian[0,index:index+3]=[1,np.sqrt(2*epsilon*(1+epsilon))*z[0],epsilon*(2*z[0]**2-1)];jacobians.append(jacobian)
            values_out=np.concatenate(values_out);jacobians=np.vstack(jacobians);factors=np.r_[np.full(3*len(ids),1e8),np.full(len(self.nuisance_epsilon),1e6)]
            return values_out*factors,jacobians[:,free]*scale[free]*factors[:,None]
        bounds=[(None,None)]*len(free)
        if 0 in free:bounds[int(np.flatnonzero(free==0)[0])]=(1e-12/scale[0],None)
        if 1 in free:bounds[int(np.flatnonzero(free==1)[0])]=(-20/scale[1],20/scale[1])
        for _ in range(12):
            fit=minimize(objective,initial,jac=True,method="SLSQP",bounds=bounds,
                constraints=[dict(type="ineq",fun=lambda x:constraints(x)[0],jac=lambda x:constraints(x)[1])],options=dict(ftol=2e-11,maxiter=1200))
            parameters=unpack(fit.x);margins=self.physics_margin(parameters)[0];bad=np.flatnonzero(margins < -1e-14)
            new=set(bad[np.argsort(margins[bad])[:32]].tolist())-active_ids
            if not new:break
            active_ids.update(new);initial=fit.x
        nuisance_min=min(self.nuisance_margins(parameters));ok=np.min(margins)>=-2e-13 and nuisance_min>=-2e-13 and np.isfinite(fit.fun)
        return dict(value=float(fixed_value),q=float(fit.fun),theta=parameters,success=bool(fit.success and ok),
                    iterations=int(fit.nit),physics_margin=float(np.min(margins)),nuisance_margin=float(nuisance_min))


def profile_scans(model, sigmas, scale, central_q):
    intervals,points=[],[]
    for index,name in enumerate(PHYS):
        roots={}
        for direction,label in ((-1,"lower"),(1,"upper")):
            last=model.profile(index,model.parameters[index],model.parameters,scale);points.append((name,label,last));bracket=None
            for multiple in (.25,.5,.75,1,1.5,2,3,4.5,6,8):
                value=model.parameters[index]+direction*multiple*sigmas[index]
                if index==0:value=max(value,1.0001e-12)
                if index==1:value=float(np.clip(value,-19.999999,19.999999))
                fit=model.profile(index,value,last["theta"],scale);points.append((name,label,fit))
                if fit["success"] and fit["q"]-central_q>=1:bracket=(last,fit);break
                if fit["success"]:last=fit
            if bracket is None:raise RuntimeError(f"profile root not bracketed: {name} {label}")
            lower,upper=bracket
            for _ in range(26):
                value=(lower["value"]+upper["value"])/2;start=lower["theta"] if abs(value-lower["value"])<abs(value-upper["value"]) else upper["theta"]
                fit=model.profile(index,value,start,scale);points.append((name,label,fit))
                if fit["success"] and fit["q"]-central_q<1:lower=fit
                else:upper=fit
                if abs(upper["value"]-lower["value"])<max(1e-8,2e-6*sigmas[index]):break
            roots[label]=min((lower,upper),key=lambda item:abs(item["q"]-central_q-1))["value"]
        intervals.append(dict(parameter=name,central_value=float(model.parameters[index]),curvature_sigma=float(sigmas[index]),
            profile_lower_value=roots["lower"],profile_upper_value=roots["upper"],profile_minus=float(model.parameters[index]-roots["lower"]),
            profile_plus=float(roots["upper"]-model.parameters[index])))
    rows=[dict(parameter=name,direction=direction,value=fit["value"],delta_q=fit["q"]-central_q,success=fit["success"],
        iterations=fit["iterations"],physics_margin=fit["physics_margin"],nuisance_margin=fit["nuisance_margin"]) for name,direction,fit in points]
    return intervals,rows


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--source",type=Path,required=True);parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args();source=args.source.resolve();output=args.output.resolve()
    if output.exists():raise FileExistsError(f"refusing to overwrite {output}")
    output.mkdir(parents=True);model=Local(source);p=model.parameters;mu,jacobian=model.prediction_jacobian(p)
    np.testing.assert_allclose(mu,model.central["prediction"],rtol=3e-13,atol=3e-15);np.testing.assert_allclose(jacobian,model.central["jacobian"],rtol=3e-12,atol=3e-14)
    design=jacobian[model.mask]/np.sqrt(model.variance[model.mask,None]);scale=1/np.linalg.norm(design,axis=0);scaled=design*scale
    singular=np.linalg.svd(scaled,compute_uv=False);fisher_scaled=scaled.T@scaled;fisher_singular=np.linalg.svd(fisher_scaled,compute_uv=False);rank=np.linalg.matrix_rank(scaled)
    covariance=np.diag(scale)@np.linalg.inv(fisher_scaled)@np.diag(scale);fisher=(jacobian[model.mask].T/model.variance[model.mask])@jacobian[model.mask]
    physics_covariance=np.linalg.inv(fisher[:4,:4]-fisher[:4,4:]@np.linalg.solve(fisher[4:,4:],fisher[4:,:4]));block_difference=float(np.max(abs(physics_covariance-covariance[:4,:4])))
    published=model.published(p);published_gradient=model.published_jacobian(p);published_covariance=published_gradient@physics_covariance@published_gradient.T
    curvature_errors=np.sqrt(np.diag(published_covariance));toys=np.load(source/"toys/replicas.npz")["pub"];toy_covariance=np.cov(toys,rowvar=False,ddof=1)
    toy_errors=np.sqrt(np.diag(toy_covariance));bias=toys.mean(0)-published;physics_sigmas=np.sqrt(np.diag(physics_covariance))
    intervals,scan=profile_scans(model,physics_sigmas,scale,model.central["info"][0]);parameter_correlation=correlation(covariance)
    for index,row in enumerate(intervals):
        nuisance_index=4+np.argmax(abs(parameter_correlation[index,4:]));row.update(boundary_status="active coupled physics positivity",
            strongest_nuisance_correlation=model.nuisance_names[nuisance_index-4],correlation=float(parameter_correlation[index,nuisance_index]))
    hessians={}
    for step in (3e-4,1e-4,3e-5):
        hessian=np.zeros((model.nparams,model.nparams))
        for index in range(model.nparams):
            delta=np.zeros(model.nparams);delta[index]=step*scale[index]
            hessian[:,index]=(model.objective_gradient(p+delta)[1]*scale-model.objective_gradient(p-delta)[1]*scale)/(2*step)
        hessians[str(step)]=(hessian+hessian.T)/2
    reference_hessian=hessians["0.0001"];hessian_difference=float(np.linalg.norm(reference_hessian/2-fisher_scaled)/np.linalg.norm(fisher_scaled))
    margins,z,u=model.physics_margin(p);event_index=int(np.argmin(margins));normals=[];active=[];normal=np.zeros(model.nparams)
    normal[:4]=[u[event_index]/p[0],-model.delta_tau[event_index]*u[event_index],model.epsilon_factor[event_index]*z[event_index]*model.baseline_lt[event_index],
                model.epsilon[event_index]*(2*z[event_index]**2-1)*model.baseline_tt[event_index]]
    normals.append(normal*scale);active.append(dict(constraint="physics",event_index=event_index,row=int(model.physical_rows[event_index]),margin=float(margins[event_index]),z=float(z[event_index])))
    for block_index,epsilon in enumerate(model.nuisance_epsilon):
        index=4+3*block_index;margin,zz=angular_minimum(np.array([p[index]]),np.array([p[index+1]]),np.array([p[index+2]]),np.array([epsilon]));normal=np.zeros(model.nparams)
        normal[index:index+3]=[1,np.sqrt(2*epsilon*(1+epsilon))*zz[0],epsilon*(2*zz[0]**2-1)];normals.append(normal*scale)
        active.append(dict(constraint="feedin_tprime_below",margin=float(margin[0]),z=float(zz[0])))
    normals=np.array(normals);tangent=null_space(normals);tangent_fisher=tangent.T@fisher_scaled@tangent;tangent_singular=np.linalg.svd(tangent_fisher,compute_uv=False)
    tangent_covariance=np.diag(scale)@tangent@np.linalg.inv(tangent_fisher)@tangent.T@np.diag(scale)
    write_matrix(output/"physics_parameter_covariance.csv",physics_covariance,PHYS);write_matrix(output/"physics_parameter_correlation.csv",correlation(physics_covariance),PHYS)
    for name,matrix in (("published_covariance_curvature.csv",published_covariance),("published_correlation_curvature.csv",correlation(published_covariance)),
                        ("published_covariance_toy.csv",toy_covariance),("published_correlation_toy.csv",correlation(toy_covariance))):write_matrix(output/name,matrix,model.names)
    write_matrix(output/"tangent_projected_covariance.csv",tangent_covariance,PHYS+model.nuisance_names);write_rows(output/"physics_profile_intervals.csv",intervals);write_rows(output/"profile_scan_points.csv",scan)
    write_rows(output/"fisher_singular_values.csv",[dict(index=index+1,scaled_jacobian_singular_value=singular[index],scaled_fisher_singular_value=fisher_singular[index]) for index in range(model.nparams)])
    comparison=[]
    for index,name in enumerate(model.names):
        comparison.append(dict(quantity=name,central_value=published[index],curvature_sd=curvature_errors[index],toy_empirical_sd=toy_errors[index],toy_over_curvature=toy_errors[index]/curvature_errors[index],
            toy_mean=toys[:,index].mean(),toy_bias=bias[index],bias_over_curvature=bias[index]/curvature_errors[index],bias_over_toy_sd=bias[index]/toy_errors[index],
            coverage_curvature=np.mean(abs(toys[:,index]-published[index])<=curvature_errors[index]),coverage_toy_sd=np.mean(abs(toys[:,index]-published[index])<=toy_errors[index])))
    write_rows(output/"covariance_comparison.csv",comparison);means=np.asarray(model.central["means"]);edges=model.config["tprime_bin_edges"];table=[]
    for index in range(model.nbins):
        table.append(dict(tprime_bin=index,tprime_low=edges[index],tprime_high=edges[index+1],response_weighted_tprime=means[index],
            sigma_U=published[index],U_stat=curvature_errors[index],sigma_LT=published[model.nbins+index],LT_stat=curvature_errors[model.nbins+index],
            sigma_TT=published[2*model.nbins+index],TT_stat=curvature_errors[2*model.nbins+index],units="nb/GeV^2",covariance="published_covariance_curvature.csv"))
    write_rows(output/"preliminary_cross_sections.csv",table);xlimit=max(abs(float(edges[0])),abs(float(edges[-1])))
    for component_index,component in enumerate(COMPONENTS):
        section=slice(model.nbins*component_index,model.nbins*(component_index+1));figure,axis=plt.subplots(figsize=(6.4,4.7))
        axis.errorbar(-means,published[section],yerr=curvature_errors[section],fmt="o",capsize=5,color="#163f6b");axis.axhline(0,color=".55",lw=.8)
        axis.set(xlabel=r"$-t'\ [\mathrm{GeV}^2]$",ylabel=rf"$\sigma_{{{component}}}\ [\mathrm{{nb}}/\mathrm{{GeV}}^2]$",xlim=(0,xlimit))
        axis.set_title(f"PRELIMINARY - {model.config['configured_kinematic']}\nstatistical uncertainties only\nmodel-dependent forward-folded extraction");figure.tight_layout()
        for extension in ("png","pdf"):figure.savefig(output/f"sigma_{component}_preliminary.{extension}",dpi=180)
        plt.close(figure)
    published_eigenvalues=np.linalg.eigvalsh(published_covariance);toy_eigenvalues=np.linalg.eigvalsh(toy_covariance)
    published_rank=int(np.linalg.matrix_rank(published_covariance,tol=published_eigenvalues.max()*1e-10));toy_rank=int(np.linalg.matrix_rank(toy_covariance,tol=toy_eigenvalues.max()*1e-10))
    maximum_bias=float(np.max(abs(bias/curvature_errors)));minimum_coverage=float(min(row["coverage_curvature"] for row in comparison));ready=maximum_bias<=.5 and minimum_coverage>=.5
    write_json(output/"covariance_metadata.json",dict(parameter_order=PHYS+model.nuisance_names,published_order=model.names,
        physics_parameters=4,nuisance_parameters=len(model.nuisance_names),total_parameters=model.nparams,
        fixed_feedin_blocks=list(map(int,model.fixed_blocks)),
        final_variance="data conditional Sumw2 + exclusive-MC Sumw2",mc_covariance_added_again=False,fisher_rank=int(rank),
        scaled_fisher_condition=float(fisher_singular[0]/fisher_singular[-1]),scaled_jacobian_condition=float(singular[0]/singular[-1]),
        inversion="direct scaled-basis inverse",full_inverse_block_max_abs_difference=block_difference,published_curvature_eigenvalues=published_eigenvalues.tolist(),
        published_curvature_rank=published_rank,published_toy_eigenvalues=toy_eigenvalues.tolist(),published_toy_rank=toy_rank))
    write_json(output/"hessian_crosscheck.json",dict(convention="H_chi2/2 compared with F",relative_frobenius_difference=hessian_difference,
        step_convergence={key:float(np.linalg.norm(value-reference_hessian)/np.linalg.norm(reference_hessian)) for key,value in hessians.items()}))
    write_json(output/"boundary_tangent.json",dict(active_constraints=active,active_normal_rank=int(np.linalg.matrix_rank(normals)),tangent_dimension=int(tangent.shape[1]),
        tangent_fisher_singular_values=tangent_singular.tolist(),tangent_fisher_rank=int(np.linalg.matrix_rank(tangent_fisher)),tangent_fisher_condition=float(tangent_singular[0]/tangent_singular[-1]),
        projected_covariance="tangent_projected_covariance.csv",projected_covariance_rank=int(np.linalg.matrix_rank(tangent_covariance))))
    write_json(output/"release_verdict.json",dict(ready=ready,selected_covariance="curvature",selected_reason="primary conventional final-variance covariance; toys diagnose boundary bias",
        blockers=[] if ready else [f"local M0 toy bias reaches {maximum_bias:.3f} curvature sigma (>0.5 gate)"],bootstrap_interpretation="conservative robustness diagnostic only"))
    with PdfPages(output/"final_model_only_diagnostics.pdf") as pdf:
        figure,axis=plt.subplots(figsize=(10,7));axis.axis("off");axis.text(.03,.96,
            f"{model.config['configured_kinematic']} local curvature\n\nQ={model.central['info'][0]:.12f}; rows={model.mask.sum()}\n"
            f"theta={p[:4]}\nranks={rank}/{model.nparams-4}/4\npublished bins={model.nbins}\ncurvature published rank={published_rank}; toy rank={toy_rank}\n"
            f"max |toy bias|/curvature SD={maximum_bias:.4f}\nmin curvature coverage={minimum_coverage:.3f}\n\nVerdict: {'READY' if ready else 'NOT READY'}",va="top",family="monospace")
        pdf.savefig(figure);plt.close(figure);figure,axis=plt.subplots(figsize=(11,5));axis.plot(range(len(model.names)),toy_errors/curvature_errors,"o")
        axis.axhspan(.6,1.4,color=".85");axis.axhline(1,color=".4");axis.set_xticks(range(len(model.names)));axis.set_xticklabels(model.names,rotation=90)
        axis.set_ylabel("toy SD / curvature SD");figure.tight_layout();pdf.savefig(figure);plt.close(figure)
    (output/"REPORT.md").write_text(f"# {model.config['configured_kinematic']} local-curvature extraction\n\nFresh campaign: 500 accepted toys, {model.nbins} t-prime bins, {int(model.mask.sum())} included rows.\n\nSelected t-prime edges: `{edges}`. Published covariance rank: {published_rank}. Maximum toy bias/curvature SD: {maximum_bias:.3f}; minimum symmetric coverage: {minimum_coverage:.3f}.\n")


if __name__ == "__main__":main()
