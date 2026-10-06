"""Boundary-aware publication calibration for the joint event-level M0 fit.

The campaign mirrors the single-setting release contract: additive central
data weights, one Poisson(1) draw per physical data and SIMC event, constrained
finite-MC refits, exactly 500 accepted toys for publication, deterministic
five-fold held-out coverage, explicit bias gates, and immutable provenance.
"""

import csv
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"

import numpy as np
import uproot

from joint_m0_solver import JointM0Problem


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
from audit_pi0_exclusive_estimator import ellipse
from bootstrap_pi0_data import Fitter, MASS_EDGES, RunEvents, timing_components
from validate_pi0_timing_transport import geometry

COMPONENTS = ("U", "LT", "TT")
SCHEMA = "joint_m0_calibration_v1"
PUBLICATION_TOYS = 500


def _need(condition, message):
    if not condition:
        raise ValueError(message)


def _sha256(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path, value):
    path=Path(path);temporary=path.with_suffix(path.suffix+".partial")
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")
    temporary.replace(path)


def _csv(path, records):
    records=list(records)
    _need(records,f"refusing empty CSV: {path}")
    path=Path(path);temporary=path.with_suffix(path.suffix+".partial")
    with temporary.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    temporary.replace(path)


def _correlation(covariance):
    sd=np.sqrt(np.maximum(np.diag(covariance),0));den=np.outer(sd,sd)
    return np.divide(covariance,den,out=np.zeros_like(covariance),where=den>0)


def _matrix(path,names,matrix):
    _csv(path,({"quantity":name,**{other:float(value) for other,value in zip(names,row)}}
               for name,row in zip(names,matrix)))


def _config_rows(events, config):
    """Reproduce extractor reconstructed-row ordering for arbitrary JSON bins."""
    t_edges=np.asarray(config["tprime_bin_edges"],float);q_edges=np.asarray(config["q2_bin_edges"],float)
    xb_rows=[np.asarray(row,float) for row in config["xb_bin_edges_by_q2"]]
    nphi=int(config.get("phi_bins",len(config.get("phi_bin_edges",[]))-1))
    phi_edges=np.linspace(0,2*np.pi,nphi+1) if "phi_bins" in config else np.asarray(config["phi_bin_edges"],float)
    q=events["Q2"].astype("float32").astype(float);x=events["xB"].astype("float32").astype(float)
    tp=events["t"].astype("float32").astype(float)-events["tmin"].astype("float32").astype(float)
    phi=np.mod(events["phi"],2*np.pi)
    it=np.searchsorted(t_edges,tp,side="right")-1;iq=np.searchsorted(q_edges,q,side="right")-1
    ix=np.full(len(q),-1,dtype=int)
    for group,edges in enumerate(xb_rows):
        take=iq==group;ix[take]=np.searchsorted(edges,x[take],side="right")-1
    ip=np.searchsorted(phi_edges,phi,side="right")-1;ip[phi==2*np.pi]=nphi-1
    vertices=np.asarray(config.get("diamond_xb_q2_vertices") or [],float);inside=np.ones(len(q),bool)
    if len(vertices):
        cross=np.array([(b[0]-a[0])*(q-a[1])-(b[1]-a[1])*(x-a[0])
                        for a,b in zip(vertices,np.roll(vertices,-1,axis=0))])
        inside=np.all(cross>=-1e-12,axis=0)|np.all(cross<=1e-12,axis=0)
    nx=len(xb_rows[0])-1
    valid=(inside&(it>=0)&(it<len(t_edges)-1)&(iq>=0)&(iq<len(q_edges)-1)&
           (ix>=0)&(ix<nx)&(ip>=0)&(ip<nphi)&np.isfinite(phi))
    return (((it*len(xb_rows)+iq)*nx+ix)*nphi+ip),valid


class DataSetting:
    """Validated additive-weight data definition and physical-event replicas."""
    def __init__(self, setting, fitter):
        self.label=setting["label"];self.config=setting["config"];self.fitter=fitter
        with uproot.open(setting["data_file"]) as source:
            events=source["physics"].arrays(library="np")
            manifest=source["analysis_runs"].arrays(library="np")
        self.events=events;self.manifest=manifest
        self.charge=float(np.sum(manifest["charge_uC"].astype(float)))
        _need(self.charge>0,f"{self.label}: nonpositive exposure")
        self.geometry=geometry(Path(setting["data_file"]))
        row,kinematic=_config_rows(events,self.config)
        selected=ellipse(events["mpi0_all"],events["mmiss_all"],self.geometry)&kinematic
        stored=events["is_exclusive_ellipse_combined"]!=0
        _need(np.array_equal(ellipse(events["mpi0_all"],events["mmiss_all"],self.geometry),stored),
              f"{self.label}: stored ellipse flags differ from frozen geometry")
        self.selected=selected;self.row=row
        self.nrows=len(setting["reco"]);self.target_divisor=float(self.config["tgt_contam"])
        _need(self.target_divisor>0,f"{self.label}: nonpositive target divisor")
        self.factor=(events["scale"].astype(float)*events["charge_uC"].astype(float)/
                     self.charge/self.target_divisor)
        self.runs=[]
        for run_number in manifest["run_number"]:
            take=events["run_number"]==run_number;indices=np.flatnonzero(take)
            run=RunEvents.__new__(RunEvents);run.run=int(run_number)
            run.events={key:value[take] for key,value in events.items()}
            run.mass_bin=np.searchsorted(MASS_EDGES,run.events["mpi0_all"],side="right")-1
            run.inmass=(run.mass_bin>=0)&(run.mass_bin<200)
            run.components=timing_components(run.events);run.plane_components=timing_components(run.events,True)
            run.template_factors=np.array([0.,1/6,1/12,1/12,-1/18,-1/18])
            run.coefficient=run.events["pi0_timing_coeff"]
            _,run.inverse=np.unique(run.events["event_id"],return_inverse=True)
            run.nphysical=int(run.inverse.max()+1) if len(run.inverse) else 0
            run.selected=selected[indices]&run.inmass;run.row=row[indices];run.factor=self.factor[indices]
            self.runs.append(run)

    def central(self):
        events=self.events;weight=events["pi0_timing_coeff"]-(events["pi0_timing_bin_mean"]-events["pi0_weight"])
        _need(np.all(np.isfinite(weight[self.selected])),f"{self.label}: undefined additive selected weight")
        values=weight[self.selected]*self.factor[self.selected];row=self.row[self.selected]
        return (np.bincount(row,weights=values,minlength=self.nrows),
                np.bincount(row,weights=values*values,minlength=self.nrows))

    def replica(self,rng):
        y=np.zeros(self.nrows);variance=np.zeros(self.nrows);failures=[];empty=0.
        for run in self.runs:
            multiplicity=rng.poisson(1,run.nphysical)[run.inverse].astype(float)
            spectrum,spectrum_variance,denominator=run.spectra(multiplicity)
            code,signal,info=self.fitter.fit(spectrum,spectrum_variance)
            if code:
                failures.append(dict(run=run.run,code=int(code),fit_status=float(info[1]),
                                     covariance_status=float(info[2])))
                continue
            background=np.divide(spectrum-signal,denominator,out=np.zeros(200),where=denominator>0)
            take=run.selected;mass=run.mass_bin[take]
            weight=(run.coefficient[take]-background[mass])*run.factor[take]
            y+=np.bincount(run.row[take],weights=multiplicity[take]*weight,minlength=self.nrows)
            variance+=np.bincount(run.row[take],weights=multiplicity[take]*weight*weight,minlength=self.nrows)
            empty+=float(signal[denominator==0].sum())
        return y,variance,failures,empty


def _mc_inverse(setting, problem_event):
    vertex=Path(setting["metadata"].get("vertex_source", ""))
    _need(vertex.is_file(),f"{setting['label']}: missing raw vertex_source for physical-event toys")
    with uproot.open(vertex) as source:
        raw=source["h10"].arrays(["Q2i","ti","phipqi"],library="np")
    with uproot.open(setting["simc_file"]) as source:
        simulation=source["simulation"].arrays(["event_id","is_exclusive"],library="np")
    allowed=set(simulation["event_id"][simulation["is_exclusive"]!=0].astype(int))
    lookup={};duplicates=0
    for index,(q,t,phi) in enumerate(zip(raw["Q2i"],raw["ti"],raw["phipqi"])):
        key=(float(q),-float(t),float(phi))
        if key in lookup:duplicates+=1
        lookup[key]=index
    _need(duplicates==0,f"{setting['label']}: ambiguous raw exclusive identity")
    cached=problem_event["raw"]
    try:
        ids=np.array([lookup[(float(q),float(t),float(phi))]
                      for q,t,phi in zip(cached["Q2"],cached["t"],cached["phi"])],dtype=int)
    except KeyError as error:
        raise ValueError(f"{setting['label']}: joint cache event is absent from raw exclusive SIMC") from error
    _need(set(ids)<=allowed,f"{setting['label']}: matched event is not marked exclusive in smeared SIMC")
    unique,inverse=np.unique(ids,return_inverse=True)
    return inverse,dict(response_records=len(ids),physical_events=len(unique),raw_exclusive_events=len(raw["ti"]),
                        all_matched_exclusive=True)


def _signature(settings,bins):
    payload=dict(schema=SCHEMA,bins=bins,settings=[dict(
        kinematic=setting["label"],config_sha256=_sha256(setting["config_path"]),
        joint_event_cache_sha256=_sha256(setting["event_path"]),
        data_file=str(Path(setting["data_file"]).resolve()),simc_file=str(Path(setting["simc_file"]).resolve()),
        vertex_file=str(Path(setting["metadata"].get("vertex_source","")).resolve()),
        data_size=Path(setting["data_file"]).stat().st_size,simc_size=Path(setting["simc_file"]).stat().st_size,
        data_mtime_ns=Path(setting["data_file"]).stat().st_mtime_ns,
        simc_mtime_ns=Path(setting["simc_file"]).stat().st_mtime_ns,
        vertex_size=Path(setting["metadata"].get("vertex_source","")).stat().st_size,
        vertex_mtime_ns=Path(setting["metadata"].get("vertex_source","")).stat().st_mtime_ns,
        target_divisor=float(setting["config"]["tgt_contam"]),
        target_divisor_error=float(setting["config"]["tgt_contam_err"])) for setting in settings])
    return json.loads(json.dumps(payload))


def check_campaign(campaign,settings,bins,required=PUBLICATION_TOYS):
    campaign=Path(campaign)
    manifest=json.loads((campaign/"campaign_manifest.json").read_text())
    _need(manifest["signature"]==_signature(settings,bins),"toy campaign inputs/configuration do not match this joint fit")
    summary=json.loads((campaign/"toys/summary.json").read_text())
    toys=np.load(campaign/"toys/replicas.npz")
    _need(summary.get("accepted")==required and toys["pub"].shape[0]==required,
          f"publication requires exactly {required} accepted matching toys")
    _need(len(np.unique(toys["index"]))==required,"toy campaign contains duplicate accepted IDs")
    return manifest,summary,toys


_WORKER=None


def _replica(index):
    problem,data,mc_inverse,central,options,seed=_WORKER
    rng=np.random.default_rng(np.random.SeedSequence([seed,index]))
    ys=[];variances=[];mass_failures=[];empty=0.
    for sample in data:
        y,v,failed,lost=sample.replica(rng);ys.append(y);variances.append(v)
        mass_failures.extend([dict(kinematic=sample.label,**failure) for failure in failed]);empty+=lost
    if mass_failures:
        return index,"failure",dict(replica=index,kind="mass_fit",details=mass_failures)
    multipliers=[]
    for inverse in mc_inverse:
        multipliers.append(rng.poisson(1,int(inverse.max())+1)[inverse].astype(float))
    y=np.concatenate(ys)-central["y"]+central["prediction"];variance=np.concatenate(variances)
    try:
        result=problem.fit(options["variance_mode"],True,1,options["max_iterations"],options["tolerance"],
                           y=y,data_variance=variance,multipliers=multipliers,initial=central["parameters"])
        names,published,means,_,_=problem.published(result["parameters"],multipliers)
    except Exception as error:
        return index,"failure",dict(replica=index,kind="joint_M0",details=str(error))
    if not np.all(np.isfinite(published)):
        return index,"failure",dict(replica=index,kind="published_vector",details="nonfinite derived value")
    return index,"accepted",dict(index=index,theta=result["parameters"],pub=published,means=means,
        q=result["objective"],iterations=result["iterations"],boundary_hit=bool(result["boundary"]),empty=empty)


def _cpus():
    try:return len(os.sched_getaffinity(0))
    except (AttributeError,OSError):return os.cpu_count() or 1


def generate_campaign(campaign,settings,bins,*,replicas=PUBLICATION_TOYS,jobs=0,seed=20261007,
                      variance_mode="finite-mc",starts=6,rank_tolerance=1e-10,
                      max_iterations=30,tolerance=1e-6,environment=None):
    campaign=Path(campaign).resolve()
    _need(replicas==PUBLICATION_TOYS,f"publication campaign requires exactly {PUBLICATION_TOYS} accepted toys")
    _need(not campaign.exists(),f"refusing existing toy campaign: {campaign}")
    campaign.parent.mkdir(parents=True,exist_ok=True);campaign.mkdir();(campaign/"build").mkdir();(campaign/"before").mkdir()
    shutil.copy2(HERE/"joint_m0_solver.py",campaign/"before/joint_m0_solver.py")
    shutil.copy2(Path(__file__),campaign/"before/joint_m0_release.py")
    signature=_signature(settings,bins)
    manifest=dict(schema=SCHEMA,created_epoch=time.time(),signature=signature,seed=seed,
                  accepted_toys_required=PUBLICATION_TOYS,
                  estimator="joint constrained additive-weight event-level SigParam2021 M0",
                  resampling="independent physical-event Poisson(1) data and exclusive SIMC per setting",
                  parameter_contract="per-setting N_U/DeltaB_U; shared N_LT/N_TT; per-setting low-tprime U; shared low-tprime LT/TT")
    _json(campaign/"campaign_manifest.json",manifest)
    library=campaign/"build/libnps_stat.so"
    flags=subprocess.check_output(["root-config","--cflags","--libs"],text=True,env=environment).split()
    root_libdir=subprocess.check_output(["root-config","--libdir"],text=True,env=environment).strip()
    subprocess.run(["g++","-O2","-std=c++17","-shared","-fPIC",str(REPO/"src/analysis/nps_stat_bridge.cpp"),
                    *flags,"-lMinuit2",f"-Wl,-rpath,{root_libdir}","-o",str(library)+".partial"],
                   cwd=REPO,env=environment,check=True)
    Path(str(library)+".partial").replace(library)
    fitter=Fitter(library);problem=JointM0Problem(settings,bins,rank_tolerance)
    data=[DataSetting(setting,fitter) for setting in settings]
    central_y=[];central_v=[]
    for sample in data:
        y,v=sample.central();central_y.append(y);central_v.append(v)
    central_y=np.concatenate(central_y);central_v=np.concatenate(central_v)
    central_result=problem.fit(variance_mode,True,starts,max_iterations,tolerance,
                               y=central_y,data_variance=central_v)
    names,published,means,published_jacobian,metadata=problem.published(central_result["parameters"])
    central=dict(parameters=central_result["parameters"],pub=published,means=means,y=central_y,
                 data_variance=central_v,prediction=central_result["prediction"],variance=central_result["variance"],
                 jacobian=central_result["jacobian"],published_jacobian=published_jacobian,
                 curvature_covariance=central_result["curvature_covariance"],objective=central_result["objective"],
                 rank=central_result["rank"],condition=central_result["condition"],boundary=central_result["boundary"])
    np.savez_compressed(campaign/"central.npz",**central)
    _json(campaign/"central.json",dict(parameter_names=problem.names,published_names=names,
        published_metadata=metadata,parameters=central_result["parameters"].tolist(),published=published.tolist(),
        means=means.tolist(),objective=central_result["objective"],rank=central_result["rank"],
        condition=central_result["condition"],boundary=bool(central_result["boundary"]),
        rows=int(np.count_nonzero(central_result["fit_mask"]))))
    mc_inverse=[];identities=[]
    for setting,event in zip(settings,problem.events):
        inverse,identity=_mc_inverse(setting,event);mc_inverse.append(inverse);identities.append(dict(kinematic=setting["label"],**identity))
    np.savez_compressed(campaign/"mc_identity.npz",**{f"setting_{i}":value for i,value in enumerate(mc_inverse)})
    _json(campaign/"mc_identity.json",identities)
    destination=campaign/"toys";destination.mkdir();accepted=[];failures=[];tic=time.monotonic()
    workers=_cpus() if jobs==0 else jobs;workers=max(1,min(workers,replicas*3))
    options=dict(variance_mode=variance_mode,max_iterations=max_iterations,tolerance=tolerance)
    global _WORKER
    _WORKER=(problem,data,mc_inverse,central,options,seed)
    ids=range(replicas*3);pool=None
    if workers>1:
        pool=mp.get_context("fork").Pool(workers);results=pool.imap(_replica,ids,chunksize=1)
    else:results=map(_replica,ids)
    try:
        for attempted,(index,status,record) in enumerate(results,1):
            _json(destination/"progress.json",dict(attempted=attempted,accepted=len(accepted),replica=index,jobs=workers))
            if status=="accepted":accepted.append(record)
            else:failures.append(record);_json(destination/"failures.json",failures)
            if status=="accepted" and (len(accepted)%20==0 or len(accepted)==replicas):
                np.savez_compressed(destination/"replicas.npz",**{key:np.asarray([item[key] for item in accepted]) for key in accepted[0]})
                summary=dict(requested=replicas,attempted=attempted,accepted=len(accepted),fit_failures=len(failures),
                    mass_fit_failures=sum(item["kind"]=="mass_fit" for item in failures),
                    joint_M0_failures=sum(item["kind"]=="joint_M0" for item in failures),
                    boundary_hits=sum(bool(item["boundary_hit"]) for item in accepted),seed=seed,jobs=workers,
                    available_cpus=_cpus(),runtime_seconds=time.monotonic()-tic)
                _json(destination/"summary.json",summary);_json(destination/"failures.json",failures)
                print(json.dumps(summary),flush=True)
            if len(accepted)>=replicas:break
    finally:
        if pool is not None:pool.terminate();pool.join()
    _need(len(accepted)==replicas,"accepted toy target not reached; inspect campaign failures")
    _json(campaign/"generation_summary.json",dict(status="complete",accepted_toys=replicas,
          finished_epoch=time.time(),freshly_generated=True))
    return campaign


def _calibrate(residuals,toy_ids,names):
    bias=residuals.mean(0);sd=residuals.std(0,ddof=1);rmse=np.sqrt(np.mean(residuals**2,axis=0))
    q16,q50,q84=np.quantile(residuals,[.16,.5,.84],axis=0,method="linear")
    delta=np.quantile(np.abs(residuals),.68,axis=0,method="linear")
    folds=toy_ids%5;decision=np.zeros_like(residuals,dtype=bool);rows=[]
    for quantity,name in enumerate(names):
        for fold in range(5):
            validation=folds==fold;training=~validation
            radius=float(np.quantile(np.abs(residuals[training,quantity]),.68,method="linear"))
            covered=np.abs(residuals[validation,quantity])<=radius;decision[validation,quantity]=covered
            coverage=float(np.mean(covered));rows.append(dict(quantity=name,fold=fold,training_n=int(training.sum()),
                validation_n=int(validation.sum()),delta68_train=radius,covered_n=int(covered.sum()),coverage68=coverage,
                binomial_standard_error=float(np.sqrt(coverage*(1-coverage)/covered.size))))
    coverage=decision.mean(0);se=np.sqrt(coverage*(1-coverage)/len(residuals))
    rows=[dict(quantity=name,fold="aggregate",training_n="see folds",validation_n=len(residuals),
               delta68_train="fold-specific",covered_n=int(decision[:,i].sum()),coverage68=float(coverage[i]),
               binomial_standard_error=float(se[i])) for i,name in enumerate(names)]+rows
    return dict(bias=bias,sd=sd,rmse=rmse,q16=q16,q50=q50,q84=q84,delta=delta,
                coverage=coverage,coverage_se=se,coverage_rows=rows)


def _plot_release(stage,names,metadata,truth,means,stats,curvature,toy_covariance,summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pages=[]
    settings=sorted({item["setting_index"] for item in metadata})
    for component in COMPONENTS:
        figure,axis=plt.subplots(figsize=(8,5.8))
        for setting in settings:
            ids=[i for i,item in enumerate(metadata) if item["setting_index"]==setting and item["component"]==component]
            order=sorted(ids,key=lambda i:metadata[i]["it"]);label=metadata[order[0]]["kinematic"]
            x=-means[order];axis.errorbar(x,truth[order],yerr=stats["delta"][order],fmt="o-",capsize=4,label=label)
        axis.axhline(0,color=".5",lw=.8);axis.set_xlabel(r"$-t'$ [GeV$^2$]")
        axis.set_ylabel(rf"$\sigma_{{{component}}}$ [nb/GeV$^2$]")
        axis.set_title("PRELIMINARY joint M0\ntoy-calibrated local 68% statistical uncertainties")
        axis.legend();figure.tight_layout();pdf=stage/f"joint_sigma_{component}_preliminary_calibrated.pdf"
        figure.savefig(pdf);figure.savefig(pdf.with_suffix(".png"),dpi=180);plt.close(figure);pages.append(pdf)
    with PdfPages(stage/"joint_calibration_diagnostics.pdf") as book:
        figure,axis=plt.subplots(figsize=(11,8.5));axis.axis("off")
        axis.text(.04,.96,summary,va="top",family="monospace",fontsize=10);book.savefig(figure);plt.close(figure)
        figure,axes=plt.subplots(2,1,figsize=(11,8.5),sharex=True);x=np.arange(len(names))
        curvature_sd=np.sqrt(np.maximum(np.diag(curvature),0))
        axes[0].plot(x,stats["sd"]/curvature_sd,"o",label="toy SD / curvature sigma")
        axes[0].plot(x,stats["delta"]/curvature_sd,"s",label="delta68 / curvature sigma")
        axes[0].axhline(1,color=".5");axes[0].legend();axes[0].set_ylabel("ratio")
        axes[1].errorbar(x,stats["coverage"],yerr=stats["coverage_se"],fmt="o",capsize=3)
        axes[1].axhline(.68,color=".5");axes[1].axhspan(.60,.76,color="#d7ecd9",alpha=.7)
        axes[1].set_ylabel("held-out coverage");axes[1].set_xticks(x,names,rotation=90,fontsize=6)
        figure.tight_layout();book.savefig(figure);plt.close(figure)
        figure,axes=plt.subplots(1,2,figsize=(12,5))
        for axis,matrix,title in ((axes[0],_correlation(curvature),"curvature correlation"),
                                  (axes[1],_correlation(toy_covariance),"toy correlation")):
            image=axis.imshow(matrix,vmin=-1,vmax=1,cmap="coolwarm");axis.set_title(title)
            axis.set_xticks(range(len(names)),names,rotation=90,fontsize=5);axis.set_yticks(range(len(names)),names,fontsize=5)
            figure.colorbar(image,ax=axis,fraction=.046)
        figure.tight_layout();book.savefig(figure);plt.close(figure)
    pages.insert(0,stage/"joint_calibration_diagnostics.pdf")
    return pages


def publish_release(campaign,output,settings,bins,*,fit_output=None,final_pdf=None):
    campaign=Path(campaign).resolve();output=Path(output).resolve()
    manifest,toy_summary,toys=check_campaign(campaign,settings,bins)
    _need(not output.exists(),f"refusing existing calibrated release: {output}")
    central_np=np.load(campaign/"central.npz");central_json=json.loads((campaign/"central.json").read_text())
    names=central_json["published_names"];metadata=central_json["published_metadata"]
    truth=np.asarray(central_np["pub"],float);means=np.asarray(central_np["means"],float)
    samples=np.asarray(toys["pub"],float);_need(samples.shape==(PUBLICATION_TOYS,len(names)),"invalid toy published-vector shape")
    toy_ids=np.asarray(toys["index"],int);residuals=samples-truth;stats=_calibrate(residuals,toy_ids,names)
    toy_cov=np.cov(residuals,rowvar=False,ddof=1);toy_corr=_correlation(toy_cov)
    calibrated=np.outer(stats["delta"],stats["delta"])*toy_corr
    parameter_cov=np.asarray(central_np["curvature_covariance"],float)
    parameter_names=central_json["parameter_names"]
    parameter_residuals=np.asarray(toys["theta"],float)-np.asarray(central_np["parameters"],float)
    parameter_toy_cov=np.cov(parameter_residuals,rowvar=False,ddof=1)
    published_jac=np.asarray(central_np["published_jacobian"],float)
    curvature=published_jac@parameter_cov@published_jac.T;curvature_sd=np.sqrt(np.maximum(np.diag(curvature),0))
    minimum_coverage=float(np.min(stats["coverage"]));maximum_bias=float(np.max(np.abs(stats["bias"])/stats["delta"]))
    finite=bool(np.all(np.isfinite(stats["delta"])) and np.all(stats["delta"]>0))
    _need(finite and minimum_coverage>=.55 and maximum_bias<1,
          f"calibrated release gate failed: finite={finite}, min_cv_coverage={minimum_coverage}, max_bias_over_delta68={maximum_bias}")
    stage=Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.",dir=output.parent))
    try:
        detail=[];table=[]
        fractions={i:float(setting["config"]["tgt_contam_err"])/float(setting["config"]["tgt_contam"])
                   for i,setting in enumerate(settings)}
        systematic=np.zeros(len(names))
        for index,(name,item) in enumerate(zip(names,metadata)):
            use_basic=not (stats["q16"][index]<=0<=stats["q84"][index])
            if not use_basic:
                widths=(stats["q84"][index],-stats["q16"][index])
                use_basic=min(widths)<=0 or max(widths)/min(widths)>=3
            systematic[index]=truth[index]*fractions[item["setting_index"]]
            detail.append(dict(quantity=name,central=truth[index],curvature_sigma=curvature_sd[index],
                toy_mean_bias=stats["bias"][index],toy_sd=stats["sd"][index],toy_rmse=stats["rmse"][index],
                delta68=stats["delta"][index],residual_q16=stats["q16"][index],residual_q50=stats["q50"][index],
                residual_q84=stats["q84"][index],basic_interval_lower=truth[index]-stats["q84"][index],
                basic_interval_upper=truth[index]-stats["q16"][index],use_basic_interval=int(use_basic),
                cross_validated_coverage68=stats["coverage"][index],coverage_binomial_standard_error=stats["coverage_se"][index],
                target_scale_fraction=fractions[item["setting_index"]],target_systematic_absolute=abs(systematic[index]),units="nb/GeV^2"))
            table.append(dict(**item,response_weighted_tprime=means[index],value=truth[index],
                delta68_stat=stats["delta"][index],basic_interval_lower=truth[index]-stats["q84"][index],
                basic_interval_upper=truth[index]-stats["q16"][index],target_systematic=abs(systematic[index]),
                units="nb/GeV^2",status="PRELIMINARY; joint local-M0 toy-calibrated 68% statistical uncertainty"))
        _csv(stage/"toy_residual_statistics.csv",detail);_csv(stage/"cross_validated_coverage.csv",stats["coverage_rows"])
        _csv(stage/"joint_preliminary_cross_sections_calibrated.csv",table)
        _matrix(stage/"published_covariance_curvature.csv",names,curvature)
        _matrix(stage/"published_covariance_toy.csv",names,toy_cov)
        _matrix(stage/"published_correlation_toy.csv",names,toy_corr)
        _matrix(stage/"published_covariance_calibrated68.csv",names,calibrated)
        _matrix(stage/"target_systematic_covariance.csv",names,np.outer(systematic,systematic))
        _matrix(stage/"joint_parameter_covariance_curvature.csv",parameter_names,parameter_cov)
        _matrix(stage/"joint_parameter_covariance_toy.csv",parameter_names,parameter_toy_cov)
        _matrix(stage/"joint_parameter_correlation_toy.csv",parameter_names,_correlation(parameter_toy_cov))
        _csv(stage/"joint_parameter_toy_statistics.csv",(
            dict(parameter=name,central=float(central_np["parameters"][index]),
                 toy_mean_bias=float(np.mean(parameter_residuals[:,index])),
                 toy_sd=float(np.std(parameter_residuals[:,index],ddof=1)),
                 residual_q16=float(np.quantile(parameter_residuals[:,index],.16)),
                 residual_q50=float(np.quantile(parameter_residuals[:,index],.50)),
                 residual_q84=float(np.quantile(parameter_residuals[:,index],.84)))
            for index,name in enumerate(parameter_names)))
        summary_text=(f"Joint calibrated preliminary M0 release\n\naccepted toys = {PUBLICATION_TOYS}\n"
            f"fit parameters = {len(central_json['parameter_names'])}; published derived quantities = {len(names)}\n"
            f"central Q = {central_json['objective']:.12g}; rank = {central_json['rank']}\n"
            f"scaled condition = {central_json['condition']:.6g}; boundary active = {central_json['boundary']}\n"
            f"held-out coverage range = {stats['coverage'].min():.3f} to {stats['coverage'].max():.3f}\n"
            f"max |bias| / delta68 = {maximum_bias:.6f}\n\nVerdict: PRELIMINARY JOINT MODEL EXTRACTION READY")
        pages=_plot_release(stage,names,metadata,truth,means,stats,curvature,toy_cov,summary_text)
        _json(stage/"calibration_summary.json",dict(schema=SCHEMA,publication_ready=True,accepted_toys=PUBLICATION_TOYS,
            minimum_cross_validated_coverage=minimum_coverage,maximum_abs_bias_over_delta68=maximum_bias,
            central_values_bias_corrected=False,statistical_interval="delta68=Q0.68(abs(toy fit - generating truth))",
            covariance_representation="diag(delta68) R_toy diag(delta68)",target_systematic="separate fully correlated divisor scale",
            parameter_count=len(central_json["parameter_names"]),published_quantity_count=len(names),toy_summary=toy_summary))
        _json(stage/"input_manifest.json",dict(campaign=str(campaign),campaign_manifest_sha256=_sha256(campaign/"campaign_manifest.json"),
            central_sha256=_sha256(campaign/"central.npz"),replicas_sha256=_sha256(campaign/"toys/replicas.npz"),signature=manifest["signature"]))
        (stage/"REPORT.md").write_text(
            "# Preliminary calibrated joint M0 extraction\n\n"
            f"Exactly {PUBLICATION_TOYS} accepted physical-event data+SIMC toys were used. "
            f"Held-out coverage spans {minimum_coverage:.3f}-{float(np.max(stats['coverage'])):.3f}; "
            f"maximum |bias|/delta68 is {maximum_bias:.3f}.\n\n"
            "Central values are not bias corrected. Statistical radii are local-M0 toy-calibrated delta68 values; "
            "the target-divisor scale covariance is stored separately and is not added in quadrature. "
            "All per-setting U normalizations/slopes and shared LT/TT parameters are refitted in every toy.\n")
        fit_output=Path(fit_output).resolve() if fit_output else None
        final_pdf=Path(final_pdf).resolve() if final_pdf else output/"joint_preliminary_cross_section_report.pdf"
        _need(final_pdf.parent==output,"final report PDF must be directly inside the calibrated output directory")
        sources=pages+([fit_output/"all_joint_xsec_plots.pdf"] if fit_output and (fit_output/"all_joint_xsec_plots.pdf").is_file() else [])
        subprocess.run(["pdfunite",*(str(path) for path in sources),str(stage/final_pdf.name)],check=True)
        _json(stage/"pipeline_linkage.json",dict(schema_version=1,release_estimator=manifest["estimator"],
            uncertainties_apply_to=str(output/"joint_preliminary_cross_sections_calibrated.csv"),
            generic_fit_output=str(fit_output) if fit_output else None,toy_campaign=str(campaign),final_report_pdf=str(final_pdf)))
        output.parent.mkdir(parents=True,exist_ok=True);stage.replace(output)
    except Exception:
        shutil.rmtree(stage,ignore_errors=True);raise
    return output/final_pdf.name
