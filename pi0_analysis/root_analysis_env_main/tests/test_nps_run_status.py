"""Run-success and exposure contract tests; no production files are modified."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd
import uproot

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("combine_status_test", ROOT / "src/analysis/combine_analysis_branches.py")
combine = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = combine
spec.loader.exec_module(combine)


class StatusTests(unittest.TestCase):
    def test_expected_exclusion_requires_complete_manifest(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'quality.csv'
            path.write_text('run,accepted,exclusion_reason\n1,1,\n2,0,missing_input\n')
            self.assertEqual(combine.accepted_manifest_runs(path,{1,2}),{1})
            with self.assertRaises(RuntimeError):
                combine.accepted_manifest_runs(path,{1,2,3})

    def test_certified_zero_does_not_require_shape_covariance(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder)
            (root/'analysis_status_run1.csv').write_text(
                'run,success,stage,fit_valid,minimizer_status,covariance_status,classification\n'
                '1,1,complete,1,1,2,zero_background\n')
            self.assertEqual(combine.require_run_success(root,1)['classification'],'zero_background')

    def test_negative_event_receives_fitted_geometric_selector(self):
        rng=np.random.default_rng(7)
        signal=rng.multivariate_normal([.13498,.938],
            [[.0035**2,-.82*.0035*.065],[-.82*.0035*.065,.065**2]],12000)
        bx=rng.normal(.124,.006,45000)
        by=1.32-18*(bx-.124)+rng.normal(0,.10,45000)
        n=57000
        data=pd.DataFrame({"mpi0_all": np.r_[signal[:,0],bx,.135],
            "mmiss_all": np.r_[signal[:,1],by,.938],
            "pi0_weight": np.r_[np.ones(n),-.25], "scale": np.ones(n+1)})
        debug=combine.add_combined_2d_mass_cut(data)
        self.assertIsNotNone(debug)
        self.assertEqual(data.iloc[-1].is_exclusive_ellipse_combined,1)
        hist=debug["histograms"]["combined_2d_mass_cut_h_mmiss_vs_mpi0_weighted"][0]
        cfg=combine.MASS_CUT_CONFIG
        inside=((data.mpi0_all>=cfg["mpi0_min"]) & (data.mpi0_all<cfg["mpi0_max"]) &
            (data.mmiss_all>=cfg["mmiss_min"]) & (data.mmiss_all<cfg["mmiss_max"]))
        self.assertAlmostEqual(hist.sum(),data.loc[inside,"pi0_weight"].sum())

    def test_failed_or_missing_run_blocks_combine_before_any_exclusion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lookup = {6419: SimpleNamespace(target="LH2")}
            for status in [None, "6419,0,combinatorial_background_fit,0,1,2,invalid_covariance"]:
                if status:
                    (root / "analysis_status_run6419.csv").write_text(
                        "run,success,stage,fit_valid,minimizer_status,covariance_status,reason\n" + status + "\n")
                with self.assertRaises(RuntimeError):
                    combine.combine_branches(lookup,root,"LH2",{})

    def test_exposure_manifest_matches_exact_yield_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = pd.DataFrame({"run_number": np.array([11,11,12],dtype=np.int32),
                "charge_uC": np.array([100,100,250],dtype=np.float32),
                "scale": np.array([2,2,3],dtype=np.float32),
                "pi0_weight": np.array([1,-.25,.5])})
            combine.save_to_root(data,root/"combined.root")
            with uproot.open(root/"combined.root") as f:
                manifest=f["analysis_runs"].arrays(library="np")
                self.assertEqual(set(manifest["run_number"]),{11,12})
                self.assertEqual(manifest["charge_uC"].sum(),350)
                self.assertEqual(f["physics"]["pi0_weight"].array(library="np")[1],-.25)
            quality=root/'quality.csv'
            quality.write_text('run,accepted,exclusion_reason\n11,1,\n12,0,invalid_background_fit\n')
            accepted=combine.accepted_manifest_runs(quality,{11,12})
            combine.save_to_root(data[data.run_number.isin(accepted)],root/"restricted.root")
            with uproot.open(root/"restricted.root") as f:
                self.assertEqual(f["analysis_runs"]["charge_uC"].array(library="np").sum(),100)
                self.assertEqual(set(f["physics"]["run_number"].array(library="np")),{11})


if __name__ == "__main__":
    unittest.main()
