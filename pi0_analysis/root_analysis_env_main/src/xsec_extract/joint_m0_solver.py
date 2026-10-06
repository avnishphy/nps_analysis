"""Event-level joint SigParam M0 fit for multiple epsilon settings.

The fitted model has an independent U normalization and U slope for every
setting, shared LT/TT normalizations, an independent low-tprime U coefficient
for every setting, and shared low-tprime LT/TT coefficients.  Q2/xB exterior
events are evaluated with their cached nominal model values and never create
fit coordinates.
"""

import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


COMPONENTS = ("U", "LT", "TT")


def _write(path, fields, records):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def _minimum(u, lt, tt, epsilon):
    u, lt, tt, epsilon = np.broadcast_arrays(u, lt, tt, epsilon)
    linear = np.sqrt(2 * epsilon * (1 + epsilon)) * lt
    quadratic = 2 * epsilon * tt
    endpoint = np.where(linear > 0, -1.0, 1.0)
    z = np.where(quadratic > 0,
                 np.clip(np.divide(-linear, 2 * quadratic,
                                   out=np.zeros_like(linear), where=quadratic > 0), -1, 1),
                 endpoint)
    return u + linear*z + epsilon*tt*(2*z*z-1), z


class JointM0Problem:
    def __init__(self, settings, bins, rank_tolerance=1e-10):
        self.settings = settings
        self.nsettings = len(settings)
        self.nparameters = 3*self.nsettings + 4
        self.rank_tolerance = rank_tolerance
        self.nt=len(bins["tprime_bin_edges"])-1
        self.nq=len(bins["q2_bin_edges"])-1
        self.nx=len(bins["xb_bin_edges_by_q2"][0])-1
        self.nrows_by_setting = [len(setting["reco"]) for setting in settings]
        self.offsets = np.cumsum([0] + self.nrows_by_setting[:-1]).astype(int)
        self.nrows = sum(self.nrows_by_setting)
        self.y = np.concatenate([[float(row["data"]) for row in setting["reco"]]
                                 for setting in settings])
        self.data_variance = np.concatenate([
            [float(row["data_variance"]) for row in setting["reco"]]
            for setting in settings])
        self.fit_mask = self.data_variance > 0
        self.events = []
        self.pivots = []
        self.physical_blocks = None
        self.tprime_block = None
        for setting_index, setting in enumerate(settings):
            path = Path(setting["event_path"])
            event = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
            if event.ndim == 0:
                event = np.array([event], dtype=event.dtype)
            treatment = np.asarray(event["treatment"], dtype=str)
            physical = treatment == "physics_model"
            nuisance = treatment == "fitted_tprime_feedin"
            fixed = treatment == "fixed_model_feedin"
            if not np.all(physical | nuisance | fixed):
                raise ValueError(f"{setting['label']}: unknown joint event treatment")
            blocks = sorted(set(event["truth_block"][physical].astype(int)))
            nuisance_blocks = sorted(set(event["truth_block"][nuisance].astype(int)))
            if len(nuisance_blocks) != 1:
                raise ValueError(f"{setting['label']}: expected one populated low-tprime block")
            if self.physical_blocks is None:
                self.physical_blocks = blocks
                self.tprime_block = nuisance_blocks[0]
            elif blocks != self.physical_blocks or nuisance_blocks[0] != self.tprime_block:
                raise ValueError("joint settings differ in populated physical/low-tprime blocks")
            response_weight = np.asarray(event["response_weight"], dtype=float)
            if not np.all(np.isfinite(response_weight)) or np.any(response_weight < 0):
                raise ValueError(f"{setting['label']}: invalid event response weights")
            pivot_weight = response_weight[physical]
            if not np.sum(pivot_weight) > 0:
                raise ValueError(f"{setting['label']}: no physical response for U pivot")
            pivot = float(np.average(event["tau"][physical], weights=pivot_weight))
            self.pivots.append(pivot)
            rows = event["reco_row"].astype(int) + self.offsets[setting_index]
            basis = np.column_stack([event["basis_"+name] for name in COMPONENTS]).astype(float)
            baseline = np.column_stack([event["baseline_"+name] for name in COMPONENTS]).astype(float)
            if not (np.all(np.isfinite(basis)) and np.all(np.isfinite(baseline))):
                raise ValueError(f"{setting['label']}: nonfinite event basis/model value")
            self.events.append(dict(raw=event, row=rows, local_row=event["reco_row"].astype(int),
                block=event["truth_block"].astype(int), physical=physical, nuisance=nuisance,
                fixed=fixed, basis=basis, baseline=baseline, tau=np.asarray(event["tau"],float),
                epsilon=np.asarray(event["epsilon"],float), weight=response_weight,
                delta=np.asarray(event["tau"],float)-pivot))
        self.names = []
        for setting in settings:
            suffix = setting["label"]
            self.names += [f"N_U_{suffix}", f"DeltaB_U_{suffix}"]
        self.names += ["N_LT_shared", "N_TT_shared"]
        for setting in settings:
            self.names.append(f"feedin_tprime_below_U_{setting['label']}")
        self.names += ["feedin_tprime_below_LT_shared", "feedin_tprime_below_TT_shared"]
        self.fixed_prediction = np.zeros(self.nrows)
        self.fixed_mc_variance = np.zeros(self.nrows)
        for event in self.events:
            value = np.sum(event["basis"][event["fixed"]] * event["baseline"][event["fixed"]], axis=1)
            self.fixed_prediction += np.bincount(event["row"][event["fixed"]], weights=value,
                                                 minlength=self.nrows)
            self.fixed_mc_variance += np.bincount(event["row"][event["fixed"]], weights=value*value,
                                                  minlength=self.nrows)
        expected_fixed = np.concatenate([[float(row["fixed_feedin_prediction"])
                                          for row in setting["reco"]] for setting in settings])
        expected_fixed_var = np.concatenate([[float(row["fixed_feedin_mc_variance"])
                                              for row in setting["reco"]] for setting in settings])
        if not np.allclose(self.fixed_prediction, expected_fixed, rtol=2e-12, atol=1e-18):
            raise ValueError("joint fixed-feed-in event sum does not reconstruct row export")
        if not np.allclose(self.fixed_mc_variance, expected_fixed_var, rtol=2e-12, atol=1e-24):
            raise ValueError("joint fixed-feed-in MC variance does not reconstruct row export")

    @property
    def shared_lt(self):
        return 2*self.nsettings

    @property
    def shared_tt(self):
        return 2*self.nsettings + 1

    @property
    def nuisance_u(self):
        return 2*self.nsettings + 2

    @property
    def nuisance_lt(self):
        return 3*self.nsettings + 2

    @property
    def nuisance_tt(self):
        return 3*self.nsettings + 3

    def initial(self):
        p = np.zeros(self.nparameters)
        for setting in range(self.nsettings):
            p[2*setting] = 1.0
        p[self.shared_lt:self.shared_tt+1] = 1.0
        for setting, event in enumerate(self.events):
            take = event["nuisance"]
            response = np.sum(np.abs(event["basis"][take, 0]))
            observed = np.sum(self.y[self.offsets[setting]:self.offsets[setting]+self.nrows_by_setting[setting]])
            p[self.nuisance_u+setting] = max(1e-20, observed/response if response > 0 else 1e-6)
        return self.make_feasible(p)

    def make_feasible(self, parameters):
        p = np.asarray(parameters, float).copy()
        for setting, event in enumerate(self.events):
            take = event["physical"]
            shape = event["baseline"][take, 0] * np.exp(-p[2*setting+1]*event["delta"][take])
            lt = p[self.shared_lt]*event["baseline"][take, 1]
            tt = p[self.shared_tt]*event["baseline"][take, 2]
            angular, _ = _minimum(np.zeros(len(shape)), lt, tt, event["epsilon"][take])
            valid = shape > 0
            if np.any(valid):
                required = np.max(np.maximum(0, -angular[valid]/shape[valid]))
                p[2*setting] = max(p[2*setting], np.nextafter(required, math.inf), 1e-12)
            take = event["nuisance"]
            epsilon = float(np.max(event["epsilon"][take]))
            angular, _ = _minimum(np.array([0.]), np.array([p[self.nuisance_lt]]),
                                  np.array([p[self.nuisance_tt]]), np.array([epsilon]))
            p[self.nuisance_u+setting] = max(p[self.nuisance_u+setting],
                                              np.nextafter(-angular[0], math.inf), 0.)
        return p

    def evaluate(self, parameters, need_jacobian=True, multipliers=None):
        p = np.asarray(parameters, float)
        prediction = self.fixed_prediction.copy() if multipliers is None else np.zeros(self.nrows)
        mc_variance = self.fixed_mc_variance.copy() if multipliers is None else np.zeros(self.nrows)
        jacobian = np.zeros((self.nrows, self.nparameters)) if need_jacobian else None
        for setting, event in enumerate(self.events):
            multiplier = np.ones(len(event["row"])) if multipliers is None else np.asarray(
                multipliers[setting], dtype=float)
            if multiplier.shape != event["row"].shape or np.any(~np.isfinite(multiplier)) or np.any(multiplier < 0):
                raise ValueError(f"{self.settings[setting]['label']}: invalid MC event multipliers")
            physical = event["physical"]
            row = event["row"]
            value = np.zeros((len(row), 3))
            u = p[2*setting]*np.exp(-p[2*setting+1]*event["delta"][physical])*event["baseline"][physical,0]
            value[physical,0] = u
            value[physical,1] = p[self.shared_lt]*event["baseline"][physical,1]
            value[physical,2] = p[self.shared_tt]*event["baseline"][physical,2]
            nuisance = event["nuisance"]
            value[nuisance,0] = p[self.nuisance_u+setting]
            value[nuisance,1] = p[self.nuisance_lt]
            value[nuisance,2] = p[self.nuisance_tt]
            contribution = np.sum(event["basis"]*value, axis=1)
            if multipliers is not None:
                fixed = event["fixed"]
                fixed_value = np.sum(event["basis"][fixed]*event["baseline"][fixed], axis=1)
                prediction += np.bincount(row[fixed], weights=fixed_value*multiplier[fixed],
                                          minlength=self.nrows)
                mc_variance += np.bincount(row[fixed], weights=fixed_value**2*multiplier[fixed],
                                           minlength=self.nrows)
            # Fixed events are accumulated above (or in the cached nominal sums).
            fitted = physical | nuisance
            prediction += np.bincount(row[fitted], weights=contribution[fitted]*multiplier[fitted],
                                       minlength=self.nrows)
            mc_variance += np.bincount(row[fitted], weights=contribution[fitted]**2*multiplier[fitted],
                                       minlength=self.nrows)
            if need_jacobian:
                pr = row[physical]
                factor = event["basis"][physical,0]*u*multiplier[physical]
                jacobian[:,2*setting] += np.bincount(pr, weights=factor/p[2*setting], minlength=self.nrows)
                jacobian[:,2*setting+1] += np.bincount(pr, weights=-event["delta"][physical]*factor,
                                                       minlength=self.nrows)
                jacobian[:,self.shared_lt] += np.bincount(pr,
                    weights=event["basis"][physical,1]*event["baseline"][physical,1]*multiplier[physical],
                    minlength=self.nrows)
                jacobian[:,self.shared_tt] += np.bincount(pr,
                    weights=event["basis"][physical,2]*event["baseline"][physical,2]*multiplier[physical],
                    minlength=self.nrows)
                nr = row[nuisance]
                for index, column in enumerate((self.nuisance_u+setting,
                                                self.nuisance_lt,self.nuisance_tt)):
                    jacobian[:,column] += np.bincount(nr,
                        weights=event["basis"][nuisance,index]*multiplier[nuisance],minlength=self.nrows)
        return prediction, mc_variance, jacobian

    def margins(self, parameters):
        p = np.asarray(parameters, float)
        physical, guard = [], []
        for setting, event in enumerate(self.events):
            take = event["physical"]
            u = p[2*setting]*np.exp(-p[2*setting+1]*event["delta"][take])*event["baseline"][take,0]
            lt = p[self.shared_lt]*event["baseline"][take,1]
            tt = p[self.shared_tt]*event["baseline"][take,2]
            physical.append(_minimum(u,lt,tt,event["epsilon"][take])[0])
            take = event["nuisance"]
            epsilon = float(np.max(event["epsilon"][take]))
            guard.append(float(_minimum(np.array([p[self.nuisance_u+setting]]),
                np.array([p[self.nuisance_lt]]),np.array([p[self.nuisance_tt]]),
                np.array([epsilon]))[0][0]))
        return physical, np.asarray(guard)

    def _constraints(self, active, unit):
        scale = 1e9
        def function(x):
            p=x*unit; values=[]
            for setting,event in enumerate(self.events):
                ids=np.asarray(sorted(active[setting]),dtype=int)
                take=np.flatnonzero(event["physical"])[ids]
                u=p[2*setting]*np.exp(-p[2*setting+1]*event["delta"][take])*event["baseline"][take,0]
                lt=p[self.shared_lt]*event["baseline"][take,1]
                tt=p[self.shared_tt]*event["baseline"][take,2]
                epsilon=event["epsilon"][take];B=np.sqrt(2*epsilon*(1+epsilon))
                values.extend(u-B*lt+epsilon*tt)
                values.extend(u+B*lt+epsilon*tt)
                minimum,z=_minimum(u,lt,tt,epsilon)
                interior=(tt>0)&(np.abs(z)<1)
                values.extend(np.where(interior,minimum,1e-8))
                ntake=event["nuisance"]
                ep=float(np.max(event["epsilon"][ntake]));gu=p[self.nuisance_u+setting]
                gl=p[self.nuisance_lt];gt=p[self.nuisance_tt];gb=math.sqrt(2*ep*(1+ep))
                values.extend((gu-gb*gl+ep*gt,gu+gb*gl+ep*gt))
                gm,gz=_minimum(np.array([gu]),np.array([gl]),np.array([gt]),np.array([ep]))
                values.append(gm[0] if gt>0 and abs(gz[0])<1 else 1e-8)
            return np.asarray(values)*scale
        def jacobian(x):
            p=x*unit; rows=[]
            for setting,event in enumerate(self.events):
                ids=np.asarray(sorted(active[setting]),dtype=int)
                take=np.flatnonzero(event["physical"])[ids]
                u=p[2*setting]*np.exp(-p[2*setting+1]*event["delta"][take])*event["baseline"][take,0]
                lt=p[self.shared_lt]*event["baseline"][take,1]
                tt=p[self.shared_tt]*event["baseline"][take,2]
                epsilon=event["epsilon"][take];B=np.sqrt(2*epsilon*(1+epsilon))
                for z in (-1.,1.):
                    block=np.zeros((len(take),self.nparameters))
                    block[:,2*setting]=u/p[2*setting]
                    block[:,2*setting+1]=-event["delta"][take]*u
                    block[:,self.shared_lt]=B*z*event["baseline"][take,1]
                    block[:,self.shared_tt]=epsilon*event["baseline"][take,2]
                    rows.extend(block)
                minimum,z=_minimum(u,lt,tt,epsilon);interior=(tt>0)&(np.abs(z)<1)
                block=np.zeros((len(take),self.nparameters))
                block[:,2*setting]=u/p[2*setting]
                block[:,2*setting+1]=-event["delta"][take]*u
                block[:,self.shared_lt]=B*z*event["baseline"][take,1]
                block[:,self.shared_tt]=epsilon*(2*z*z-1)*event["baseline"][take,2]
                block[~interior]=0;rows.extend(block)
                ntake=event["nuisance"]
                ep=float(np.max(event["epsilon"][ntake]));gb=math.sqrt(2*ep*(1+ep))
                for z in (-1.,1.):
                    block=np.zeros(self.nparameters);block[self.nuisance_u+setting]=1
                    block[self.nuisance_lt]=gb*z;block[self.nuisance_tt]=ep;rows.append(block)
                gu=p[self.nuisance_u+setting];gl=p[self.nuisance_lt];gt=p[self.nuisance_tt]
                gm,gz=_minimum(np.array([gu]),np.array([gl]),np.array([gt]),np.array([ep]))
                block=np.zeros(self.nparameters)
                if gt>0 and abs(gz[0])<1:
                    block[self.nuisance_u+setting]=1;block[self.nuisance_lt]=gb*gz[0]
                    block[self.nuisance_tt]=ep*(2*gz[0]*gz[0]-1)
                rows.append(block)
            return np.asarray(rows)*unit*scale
        return function,jacobian

    def solve_once(self, start, variance, positive=True, *, y=None, fit_mask=None,
                   multipliers=None):
        y=self.y if y is None else np.asarray(y,float)
        fit_mask=self.fit_mask if fit_mask is None else np.asarray(fit_mask,bool)
        p=self.make_feasible(start);_,_,j=self.evaluate(p,multipliers=multipliers)
        norm=np.sqrt(np.sum(j[fit_mask]**2/variance[fit_mask,None],axis=0))
        unit=np.divide(1.,norm,out=np.ones(self.nparameters),where=norm>0)
        # Keep slopes in useful GeV^-2 units and guard against extreme scales.
        for setting in range(self.nsettings):unit[2*setting+1]=max(unit[2*setting+1],1e-3)
        unit=np.clip(unit,1e-20,1e20)
        active=[]
        physical,_=self.margins(p)
        for margins in physical:
            count=min(48,len(margins));active.append(set(np.argsort(margins)[:count].tolist()))
        def objective(x):
            q=x*unit;mu,_,jac=self.evaluate(q,multipliers=multipliers);residual=mu-y
            value=float(np.sum(residual[fit_mask]**2/variance[fit_mask]))
            gradient=2*(residual[fit_mask]/variance[fit_mask])@jac[fit_mask]
            return value,gradient*unit
        bounds=[]
        for index in range(self.nparameters):
            if index < 2*self.nsettings and index%2==0:bounds.append((1e-12/unit[index],None))
            elif index < 2*self.nsettings:bounds.append((-20/unit[index],20/unit[index]))
            else:bounds.append((None,None))
        result=None
        for exchange in range(20):
            constraints=[]
            if positive:
                function,jacobian=self._constraints(active,unit)
                constraints=[dict(type="ineq",fun=function,jac=jacobian)]
            result=minimize(objective,p/unit,jac=True,method="SLSQP",bounds=bounds,
                constraints=constraints,options=dict(ftol=2e-11,maxiter=1600,disp=False))
            p=result.x*unit
            if not positive:break
            physical,guard=self.margins(p);added=False
            for setting,margins in enumerate(physical):
                violating=np.flatnonzero(margins < -1e-14)
                new=set(violating[np.argsort(margins[violating])[:32]].tolist())-active[setting]
                if new:active[setting].update(new);added=True
            if np.any(guard < -1e-14):added=True
            if not added:break
            p=self.make_feasible(p)
        physical,guard=self.margins(p)
        margin=min([float(np.min(v)) for v in physical]+guard.tolist())
        success=bool(result.success and (not positive or margin>=-2e-13))
        return dict(parameters=p,objective=float(objective(p/unit)[0]),success=success,
                    status=int(result.status),message=str(result.message),calls=int(result.nfev),
                    unit=unit,margin=margin,active=sum(map(len,active)))

    def fit(self, variance_mode="finite-mc", positive=True, starts=6,
            max_iterations=30, tolerance=1e-6, *, y=None, data_variance=None,
            multipliers=None, initial=None):
        y=self.y if y is None else np.asarray(y,float)
        data_variance=self.data_variance if data_variance is None else np.asarray(data_variance,float)
        if y.shape != (self.nrows,) or data_variance.shape != (self.nrows,):
            raise ValueError("joint fit data and variance must contain one value per reconstructed row")
        if np.any(~np.isfinite(y)) or np.any(~np.isfinite(data_variance)) or np.any(data_variance<0):
            raise ValueError("joint fit data or variance is invalid")
        fit_mask=data_variance>0
        if np.count_nonzero(fit_mask)<=self.nparameters:
            raise ValueError("joint fit has too few positive-variance rows")
        seed=self.initial() if initial is None else self.make_feasible(initial);candidates=[]
        for attempt in range(starts):
            start=seed.copy()
            if attempt:
                start[self.shared_lt]*=-1 if attempt%2 else 1
                start[self.shared_tt]*=-1 if (attempt//2)%2 else 1
                for setting in range(self.nsettings):
                    start[2*setting]*=1.5 if attempt%2 else .7
                    start[2*setting+1]+=-.5 if attempt%2 else .5
                start=self.make_feasible(start)
            variance=data_variance.copy();result=None;settled=variance_mode=="data"
            for iteration in range(max_iterations+1):
                result=self.solve_once(start,variance,positive,y=y,fit_mask=fit_mask,
                                       multipliers=multipliers)
                if not result["success"]:break
                mu,mc,jac=self.evaluate(result["parameters"],multipliers=multipliers)
                next_variance=data_variance+(mc if variance_mode=="finite-mc" else 0)
                if variance_mode=="data":settled=True;break
                parameter_change=np.max(np.abs(result["parameters"]-start)/
                    np.maximum(np.abs(result["parameters"]),result["unit"]*1e-3))
                variance_change=np.max(np.abs(next_variance[fit_mask]-variance[fit_mask])/
                                       next_variance[fit_mask])
                start=result["parameters"];variance=next_variance
                if iteration>0 and max(parameter_change,variance_change)<tolerance:
                    settled=True;break
            result.update(attempt=attempt,iterations=iteration,settled=settled,
                          variance=variance)
            candidates.append(result)
        accepted=[r for r in candidates if r["success"] and r["settled"]]
        if not accepted:
            details="; ".join(f"start {item['attempt']}: status={item['status']} "
                f"success={int(item['success'])} settled={int(item['settled'])} "
                f"iterations={item['iterations']} message={item['message']}" for item in candidates)
            raise RuntimeError("no converged joint M0 start; "+details)
        best=min(accepted,key=lambda r:r["objective"])
        best["starts"]=candidates
        mu,mc,jac=self.evaluate(best["parameters"],multipliers=multipliers)
        best.update(prediction=mu,mc_variance=mc,jacobian=jac,y=y,
                    data_variance=data_variance,fit_mask=fit_mask)
        whitened=jac[fit_mask]/np.sqrt(best["variance"][fit_mask,None])
        column_norm=np.linalg.norm(whitened,axis=0)
        scaled=np.divide(whitened,column_norm,out=np.zeros_like(whitened),where=column_norm>0)
        singular=np.linalg.svd(scaled,compute_uv=False)
        best["singular_values"]=singular
        best["rank"]=int(np.sum(singular>self.rank_tolerance*singular[0]))
        best["condition"]=float(singular[0]/singular[-1]) if singular[-1]>0 else math.inf
        curvature=whitened.T@whitened
        best["curvature_covariance"]=np.linalg.inv(curvature)
        physical,guard=self.margins(best["parameters"]);boundary=False
        for setting,event in enumerate(self.events):
            take=event["physical"]
            u=best["parameters"][2*setting]*np.exp(-best["parameters"][2*setting+1]*event["delta"][take])*event["baseline"][take,0]
            lt=best["parameters"][self.shared_lt]*event["baseline"][take,1]
            tt=best["parameters"][self.shared_tt]*event["baseline"][take,2]
            scale=np.maximum.reduce((np.abs(u),np.abs(lt),np.abs(tt),np.full(len(u),1e-30)))
            boundary=boundary or bool(np.any(physical[setting]<=2e-8*scale))
            gscale=max(abs(best["parameters"][self.nuisance_u+setting]),
                       abs(best["parameters"][self.nuisance_lt]),
                       abs(best["parameters"][self.nuisance_tt]),1e-30)
            boundary=boundary or guard[setting]<=2e-8*gscale
        best["boundary"]=boundary
        best["covariance"]=np.full_like(best["curvature_covariance"],np.nan) if best["boundary"] else best["curvature_covariance"]
        return best

    def reporting(self, result):
        p=result["parameters"];records=[];jacobians=[]
        for setting,event in enumerate(self.events):
            for block in self.physical_blocks:
                take=event["physical"]&(event["block"]==block);weight=event["weight"][take]
                if not np.sum(weight)>0:continue
                dt=event["delta"][take]
                u=p[2*setting]*np.exp(-p[2*setting+1]*dt)*event["baseline"][take,0]
                values=[np.average(u,weights=weight),
                        p[self.shared_lt]*np.average(event["baseline"][take,1],weights=weight),
                        p[self.shared_tt]*np.average(event["baseline"][take,2],weights=weight)]
                derivatives=[]
                for component in range(3):
                    derivative=np.zeros(self.nparameters)
                    if component==0:
                        derivative[2*setting]=values[0]/p[2*setting]
                        derivative[2*setting+1]=np.average(-dt*u,weights=weight)
                    elif component==1:
                        derivative[self.shared_lt]=np.average(event["baseline"][take,1],weights=weight)
                    else:
                        derivative[self.shared_tt]=np.average(event["baseline"][take,2],weights=weight)
                    derivatives.append(derivative)
                it=block//(self.nq*self.nx);iq=(block//self.nx)%self.nq;ix=block%self.nx
                for component,value,derivative in zip(COMPONENTS,values,derivatives):
                    records.append(dict(setting_index=setting,truth_block=block,region="published",
                        it=it,iq=iq,ix=ix,component=component,value=value,error_stat_plus_mc=math.nan))
                    jacobians.append(derivative)
            block=self.tprime_block
            for component,index in zip(COMPONENTS,(self.nuisance_u+setting,self.nuisance_lt,self.nuisance_tt)):
                derivative=np.zeros(self.nparameters);derivative[index]=1
                records.append(dict(setting_index=setting,truth_block=block,region="tprime_below",
                    it=-1,iq=-1,ix=-1,component=component,value=p[index],error_stat_plus_mc=math.nan))
                jacobians.append(derivative)
        jacobians=np.asarray(jacobians);cov=jacobians@result["covariance"]@jacobians.T
        for index,record in enumerate(records):
            if np.isfinite(cov[index,index]) and cov[index,index]>=0:
                record["error_stat_plus_mc"]=math.sqrt(cov[index,index])
        return records,cov,jacobians

    def published(self, parameters, multipliers=None, scale=1e9):
        """Return the ordered per-setting physical structure vector.

        Replica MC multiplicities enter both the reported response-weighted
        means and the t-prime abscissae, matching the finite-MC refit toy.
        Q2/xB and low-tprime feed-in blocks are never published coordinates.
        """
        p=np.asarray(parameters,float);names=[];values=[];means=[];jacobians=[];metadata=[]
        for setting,event in enumerate(self.events):
            multiplier=np.ones(len(event["row"])) if multipliers is None else np.asarray(
                multipliers[setting],float)
            for component in COMPONENTS:
                for block in self.physical_blocks:
                    take=event["physical"]&(event["block"]==block)&(multiplier>0)
                    weight=event["weight"][take]*multiplier[take]
                    if not np.sum(weight)>0:
                        raise ValueError(f"{self.settings[setting]['label']}: empty replica support in truth block {block}")
                    derivative=np.zeros(self.nparameters)
                    if component=="U":
                        dt=event["delta"][take]
                        term=p[2*setting]*np.exp(-p[2*setting+1]*dt)*event["baseline"][take,0]
                        value=np.average(term,weights=weight)
                        derivative[2*setting]=value/p[2*setting]
                        derivative[2*setting+1]=np.average(-dt*term,weights=weight)
                    elif component=="LT":
                        base=np.average(event["baseline"][take,1],weights=weight)
                        value=p[self.shared_lt]*base;derivative[self.shared_lt]=base
                    else:
                        base=np.average(event["baseline"][take,2],weights=weight)
                        value=p[self.shared_tt]*base;derivative[self.shared_tt]=base
                    label=self.settings[setting]["label"]
                    names.append(f"{label}:{component}(bin{block})")
                    values.append(value*scale);jacobians.append(derivative*scale)
                    means.append(float(np.average(event["raw"]["tprime"][take],weights=weight)))
                    it=block//(self.nq*self.nx);iq=(block//self.nx)%self.nq;ix=block%self.nx
                    metadata.append(dict(setting_index=int(setting),kinematic=label,component=component,
                                         truth_block=int(block),it=int(it),iq=int(iq),ix=int(ix)))
        return names,np.asarray(values),np.asarray(means),np.asarray(jacobians),metadata


def fit_and_write(settings, bins, out_dir, *, variance_mode="finite-mc", positive=True,
                  rank_tolerance=1e-10, max_iterations=30, tolerance=1e-6, starts=6):
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    problem=JointM0Problem(settings,bins,rank_tolerance)
    result=problem.fit(variance_mode,positive,starts,max_iterations,tolerance)
    p=result["parameters"];cov=result["covariance"];raw=result["curvature_covariance"]
    parameter_records=[]
    for index,name in enumerate(problem.names):
        role=("U_model_per_setting" if index<2*problem.nsettings else
              "shared_interference_model" if index<=problem.shared_tt else
              "tprime_feedin_per_setting" if index<problem.nuisance_lt else "shared_tprime_feedin")
        parameter_records.append(dict(parameter_index=index,name=name,role=role,value=p[index],
            error_stat_plus_mc=math.sqrt(cov[index,index]) if np.isfinite(cov[index,index]) and cov[index,index]>=0 else math.nan))
    _write(out_dir/"joint_model_parameters.csv",parameter_records[0].keys(),parameter_records)
    covariance_records=[dict(parameter_i=i,parameter_j=j,stat_plus_mc_covariance=cov[i,j],
                             unconstrained_curvature_inverse=raw[i,j])
                        for i in range(problem.nparameters) for j in range(problem.nparameters)]
    _write(out_dir/"joint_model_covariance.csv",covariance_records[0].keys(),covariance_records)
    starts=[]
    for item in result["starts"]:
        record=dict(start=item["attempt"],converged=int(item["success"]),settled=int(item["settled"]),
                    status=item["status"],objective=item["objective"],mc_iterations=item["iterations"],
                    positivity_margin=item["margin"],message=item["message"])
        record.update({name:value for name,value in zip(problem.names,item["parameters"])})
        starts.append(record)
    _write(out_dir/"joint_model_starts.csv",starts[0].keys(),starts)
    singular=[dict(index=i,scaled_singular_value=value) for i,value in enumerate(result["singular_values"])]
    _write(out_dir/"joint_singular_values.csv",singular[0].keys(),singular)
    row_records=[];chi2_by_setting=[]
    for setting,(offset,count) in enumerate(zip(problem.offsets,problem.nrows_by_setting)):
        chi2=0
        for local in range(count):
            row=offset+local;included=problem.fit_mask[row]
            residual=problem.y[row]-result["prediction"][row]
            pull=residual/math.sqrt(result["variance"][row]) if included else math.nan
            if included:chi2+=pull*pull
            row_records.append(dict(setting_index=setting,reco_row=local,fit_index=row if included else -1,
                data=problem.y[row],data_variance=problem.data_variance[row],
                fixed_feedin_prediction=problem.fixed_prediction[row],
                fixed_feedin_mc_variance=problem.fixed_mc_variance[row],
                variance_used=result["variance"][row] if included else math.nan,
                mc_variance_final=result["mc_variance"][row],prediction=result["prediction"][row],
                residual=residual,pull=pull))
        chi2_by_setting.append(dict(setting_index=setting,chi2_contribution=chi2))
    _write(out_dir/"joint_rows.csv",row_records[0].keys(),row_records)
    _write(out_dir/"joint_setting_chi2.csv",chi2_by_setting[0].keys(),chi2_by_setting)
    positivity=[];physical,guard=problem.margins(p)
    for setting,event in enumerate(problem.events):
        for block in problem.physical_blocks:
            take=event["physical"]&(event["block"]==block)
            u=p[2*setting]*np.exp(-p[2*setting+1]*event["delta"][take])*event["baseline"][take,0]
            lt=p[problem.shared_lt]*event["baseline"][take,1]
            tt=p[problem.shared_tt]*event["baseline"][take,2]
            margin=float(np.min(_minimum(u,lt,tt,event["epsilon"][take])[0]))
            epsilon=float(np.max(event["epsilon"][take]))
            scale=max(np.max(np.abs(np.r_[u,lt,tt])),1e-30)
            positivity.append(dict(truth_block=block,region="published",setting_index=setting,
                epsilon_max=epsilon,minimum_response_bracket=margin,boundary_tolerance=2e-8*scale,
                feasibility_tolerance=2e-13))
        take=event["nuisance"]
        positivity.append(dict(truth_block=problem.tprime_block,region="tprime_below",setting_index=setting,
            epsilon_max=float(np.max(event["epsilon"][take])),minimum_response_bracket=guard[setting],
            boundary_tolerance=2e-8*max(abs(p[problem.nuisance_u+setting]),abs(p[problem.nuisance_lt]),
                                         abs(p[problem.nuisance_tt]),1e-30),feasibility_tolerance=2e-13))
    _write(out_dir/"joint_positivity.csv",positivity[0].keys(),positivity)
    structures,structure_cov,structure_jac=problem.reporting(result)
    _write(out_dir/"joint_structure_functions.csv",structures[0].keys(),structures)
    structure_cov_records=[dict(parameter_i=i,parameter_j=j,stat_plus_mc_covariance=structure_cov[i,j])
                           for i in range(len(structures)) for j in range(len(structures))]
    _write(out_dir/"joint_structure_covariance.csv",structure_cov_records[0].keys(),structure_cov_records)
    jac_records=[dict(row=row,parameter=column,derivative=result["jacobian"][row,column])
                 for row in range(problem.nrows) for column in range(problem.nparameters)]
    _write(out_dir/"joint_model_row_jacobian.csv",jac_records[0].keys(),jac_records)
    chi2=float(np.sum((problem.y[problem.fit_mask]-result["prediction"][problem.fit_mask])**2/
                      result["variance"][problem.fit_mask]))
    summary=dict(fit_components="M0_U_normalization_and_slope_per_setting_shared_NLT_NTT",
        settings=problem.nsettings,rows_used=int(np.sum(problem.fit_mask)),rows_available=problem.nrows,
        active_truth_blocks=len(problem.physical_blocks)+1,fitted_exterior_blocks="tprime_below_only",
        fixed_feedin="Q2_xB_exterior_nominal_event_model",parameters=problem.nparameters,
        rank=result["rank"],condition=result["condition"],chi2=chi2,
        ndf=int(np.sum(problem.fit_mask))-problem.nparameters,variance_mode=variance_mode,
        mc_iterations=result["iterations"],positive_xsec=int(positive),
        positivity_boundary_active=int(result["boundary"]),
        covariance_status="unavailable_boundary_constrained_estimate" if result["boundary"] else "conditional_inverse_information",
        input_normalization="each_setting_yield_per_mC_and_event_level_one_mC_SIMC_model_response",
        target_factor_uncertainty="not_propagated")
    with (out_dir/"joint_summary.txt").open("w") as stream:
        for key,value in summary.items():stream.write(f"{key}={value}\n")
    import uproot
    with uproot.recreate(out_dir/"joint_xsec_output.root") as root:
        root["joint_model_parameters"]={"parameter_index":np.arange(problem.nparameters,dtype=np.int32),
                                         "value":p}
        root["joint_model_covariance"]={
            "parameter_i":np.repeat(np.arange(problem.nparameters,dtype=np.int32),problem.nparameters),
            "parameter_j":np.tile(np.arange(problem.nparameters,dtype=np.int32),problem.nparameters),
            "stat_plus_mc_covariance":cov.ravel(),
            "unconstrained_curvature_inverse":raw.ravel()}
        root["joint_rows"]={key:np.asarray([record[key] for record in row_records])
                            for key in row_records[0]}
        root["joint_structure_functions"]={key:np.asarray([record[key] for record in structures])
                                            for key in structures[0] if key not in ("region","component")}
        root["analysis_metadata"]=uproot.writing.identify.to_TObjString(
            "method=joint_event_level_sigparam2021_M0\n"+
            "parameter_contract=separate_NU_DeltaBU_per_setting_shared_NLT_NTT\n"+
            "fitted_exterior=tprime_below_only\nfixed_exterior=Q2_xB_nominal_event_model\n")
    return summary,problem,result
