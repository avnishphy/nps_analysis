import csv
import math
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from joint_m0_solver import JointM0Problem, _minimum


FIELDS = ("event_index", "reco_row", "truth_block", "truth_phi", "treatment",
          "response_weight", "Q2", "W2", "t", "tprime", "tau", "theta_cm",
          "epsilon", "phi", "baseline_T", "baseline_L", "baseline_U",
          "baseline_LT", "baseline_TT", "basis_U", "basis_LT", "basis_TT")


class JointM0SolverTest(unittest.TestCase):
    def test_analytic_angular_minimum(self):
        rng=np.random.default_rng(713)
        u=rng.uniform(.1,2,40);lt=rng.normal(0,.4,40);tt=rng.normal(0,.4,40)
        epsilon=rng.uniform(.2,.85,40);minimum,_=_minimum(u,lt,tt,epsilon)
        phi=np.linspace(0,2*np.pi,20001)
        dense=np.min(u[:,None]+np.sqrt(2*epsilon*(1+epsilon))[:,None]*lt[:,None]*np.cos(phi)
                     +epsilon[:,None]*tt[:,None]*np.cos(2*phi),axis=1)
        np.testing.assert_allclose(minimum,dense,rtol=0,atol=2e-8)

    def test_parameter_inventory_fixed_closure_and_jacobian(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);settings=[]
            for setting,epsilon in enumerate((.72,.52)):
                records=[];fixed=np.zeros(12);fixed_var=np.zeros(12)
                for row in range(12):
                    phi=2*np.pi*(row+.5)/12;base=1e-3*(1+.02*row)
                    basis=(base,base*np.sqrt(2*epsilon*(1+epsilon))*np.cos(phi),
                           base*epsilon*np.cos(2*phi))
                    records.append(dict(zip(FIELDS,(len(records),row,0,row,"physics_model",1.,4.,8.,-.3,
                        -.1-.01*row,.1+.01*row,.3,epsilon,phi,1.,.2,1+.03*setting,.08,-.05,*basis))))
                    records.append(dict(zip(FIELDS,(len(records),row,1,row,"fitted_tprime_feedin",1.,4.,8.,-.8,
                        -.72,.72,.7,epsilon,phi,0.,0.,0.,0.,0.,*basis))))
                    nominal=np.array((.4,.03,-.02));value=float(np.dot(basis,nominal))
                    fixed[row]+=value;fixed_var[row]+=value*value
                    records.append(dict(zip(FIELDS,(len(records),row,3,row,"fixed_model_feedin",1.,4.,8.,-.2,
                        -.05,.05,.2,epsilon,phi,.3,.1,*nominal,*basis))))
                path=root/f"events_{setting}.csv"
                with path.open("w",newline="") as stream:
                    writer=csv.DictWriter(stream,fieldnames=FIELDS);writer.writeheader();writer.writerows(records)
                reco=[dict(data=1.,data_variance=.1,fixed_feedin_prediction=fixed[row],
                           fixed_feedin_mc_variance=fixed_var[row]) for row in range(12)]
                settings.append(dict(label=f"s{setting}",event_path=path,reco=reco))
            bins=dict(tprime_bin_edges=[-.7,0.],q2_bin_edges=[3.,5.],xb_bin_edges_by_q2=[[.2,.5]])
            problem=JointM0Problem(settings,bins)
            self.assertEqual(problem.nparameters,10)
            self.assertEqual(problem.names,["N_U_s0","DeltaB_U_s0","N_U_s1","DeltaB_U_s1",
                "N_LT_shared","N_TT_shared","feedin_tprime_below_U_s0","feedin_tprime_below_U_s1",
                "feedin_tprime_below_LT_shared","feedin_tprime_below_TT_shared"])
            self.assertFalse(any("q2" in name.lower() or "xb" in name.lower() for name in problem.names))
            p=problem.initial();prediction,mc_variance,jacobian=problem.evaluate(p)
            self.assertTrue(np.all(np.isfinite(prediction)))
            doubled=[np.full(len(event["row"]),2.) for event in problem.events]
            prediction2,mc_variance2,jacobian2=problem.evaluate(p,multipliers=doubled)
            np.testing.assert_allclose(prediction2,2*prediction,rtol=2e-14,atol=2e-18)
            np.testing.assert_allclose(mc_variance2,2*mc_variance,rtol=2e-14,atol=2e-24)
            np.testing.assert_allclose(jacobian2,2*jacobian,rtol=2e-14,atol=2e-18)
            for column in range(problem.nparameters):
                with self.subTest(parameter=problem.names[column]):
                    step=1e-5*max(abs(p[column]),1.)
                    plus=p.copy();minus=p.copy();plus[column]+=step;minus[column]-=step
                    numeric=(problem.evaluate(plus,False)[0]-problem.evaluate(minus,False)[0])/(2*step)
                    np.testing.assert_allclose(numeric,jacobian[:,column],rtol=2e-6,atol=2e-8)


if __name__ == "__main__":
    unittest.main()
