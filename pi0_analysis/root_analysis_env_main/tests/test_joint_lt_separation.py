#!/usr/bin/env python3
"""Numerical post-fit L/T checks; temporary CSVs only, no ROOT required.

Run from the repository root: python3 tests/test_joint_lt_separation.py
"""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src/xsec_extract"))
import run_joint_xsec_fit as driver


class SeparationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="joint_lt_test_")
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    def run_case(self, epsilon, blocks, covariance):
        records = []
        for block, components in enumerate(blocks):
            for setting, component, value in components:
                records.append(dict(parameter_index=len(records), truth_block=block,
                                    region="published" if block == 0 else "guard",
                                    it=block, iq=0, ix=0, setting_index=setting,
                                    component=component, value=value))
        driver.write_csv(self.folder / "joint_parameters.csv", tuple(records[0]), records)
        n = len(records)
        driver.write_csv(self.folder / "joint_covariance.csv",
                         ("parameter_i", "parameter_j", "stat_plus_mc_covariance",
                          "unconstrained_curvature_inverse"),
                         (dict(parameter_i=i, parameter_j=j,
                               stat_plus_mc_covariance=covariance[i, j],
                               unconstrained_curvature_inverse=float(i == j))
                          for i in range(n) for j in range(n)))
        original = [driver.digest(self.folder / name) for name in
                    ("joint_parameters.csv", "joint_covariance.csv")]
        summary = driver.separate_lt(self.folder, [{"epsilon": e} for e in epsilon], 1e-10)
        self.assertEqual(original, [driver.digest(self.folder / name) for name in
                                   ("joint_parameters.csv", "joint_covariance.csv")])
        rows = driver.rows(self.folder / "joint_separated_parameters.csv")
        values = np.array([float(r["value"]) for r in rows])
        cov = np.array([float(r["stat_plus_mc_covariance"]) for r in
                        driver.rows(self.folder / "joint_separated_covariance.csv")])
        return values, cov.reshape(len(rows), len(rows)), summary, rows

    @staticmethod
    def block(epsilon, t=2., longitudinal=5.):
        return [(s, "U", t+e*longitudinal) for s, e in enumerate(epsilon)] + [
            (-1, "LT", .3), (-1, "TT", -.2)]

    def test_two_settings_full_correlations_and_cross_bins(self):
        epsilon = [.25, .75]
        blocks = [self.block(epsilon), self.block(epsilon, 4., -2.)]
        rng = np.random.default_rng(451)
        a = rng.normal(size=(8, 8))
        covariance = a @ a.T + np.eye(8)
        values, cov, summary, rows = self.run_case(epsilon, blocks, covariance)
        # Closed-form T=(e1*U0-e0*U1)/(e1-e0), L=(U1-U0)/(e1-e0).
        jacobian = np.eye(8)
        for start in (0, 4):
            jacobian[start:start+2, start:start+2] = [[1.5, -.5], [-2., 2.]]
        np.testing.assert_allclose(values, [2., 5., .3, -.2, 4., -2., .3, -.2], atol=1e-13)
        np.testing.assert_allclose(cov, jacobian @ covariance @ jacobian.T, atol=1e-12)
        np.testing.assert_allclose([float(r["error_stat_plus_mc"]) for r in rows], np.sqrt(cov.diagonal()))
        self.assertTrue(all(b["status"] == "ok" for b in summary["blocks"]))

    def test_three_settings_correlated_gls(self):
        epsilon = [.2, .5, .85]
        block = self.block(epsilon)
        block[1] = (1, "U", block[1][2]+.4)
        covariance = np.array([[2., .6, .1, .2, 0.], [.6, 1., .3, 0., .1],
                               [.1, .3, 3., .1, .2], [.2, 0., .1, 1., 0.],
                               [0., .1, .2, 0., 1.]])
        values, cov, summary, _ = self.run_case(epsilon, [block], covariance)
        design = np.column_stack((np.ones(3), epsilon))
        precision = np.linalg.inv(covariance[:3, :3])
        weights = np.linalg.inv(design.T @ precision @ design) @ design.T @ precision
        transform = np.zeros((4, 5))
        transform[:2, :3] = weights
        transform[2:, 3:] = np.eye(2)
        np.testing.assert_allclose(values, transform @ np.array([r[2] for r in block]), atol=1e-12)
        np.testing.assert_allclose(cov, transform @ covariance @ transform.T, atol=1e-12)
        residual = np.array([r[2] for r in block[:3]]) - design @ values[:2]
        self.assertAlmostEqual(summary["blocks"][0]["chi2"], residual @ precision @ residual)
        self.assertEqual(summary["blocks"][0]["ndf"], 1)

    def test_degenerate_epsilon_preserves_interference(self):
        for epsilon in ([.5, .5], [.5, .5+1e-12]):
            values, cov, summary, _ = self.run_case(epsilon, [self.block(epsilon)], np.eye(4))
            self.assertTrue(np.isnan(values[:2]).all())
            self.assertTrue(np.isnan(cov[:2]).all())
            np.testing.assert_allclose(values[2:], [.3, -.2])
            np.testing.assert_allclose(cov[2:, 2:], np.eye(2))
            self.assertEqual(summary["blocks"][0]["status"], "degenerate_epsilon")

    def test_one_setting_guard_does_not_contaminate_published(self):
        epsilon = [.25, .75]
        blocks = [self.block(epsilon), [(1, "U", 4.), (-1, "LT", .1), (-1, "TT", .2)]]
        values, cov, summary, _ = self.run_case(epsilon, blocks, np.eye(7))
        np.testing.assert_allclose(values[:4], [2., 5., .3, -.2])
        self.assertTrue(np.isfinite(cov[:4, :4]).all())
        self.assertTrue(np.isnan(values[4:6]).all())
        self.assertEqual(summary["blocks"][1]["status"], "insufficient_settings")

    def test_boundary_never_uses_diagnostic_curvature(self):
        epsilon = [.25, .75]
        values, cov, summary, rows = self.run_case(epsilon, [self.block(epsilon)], np.full((4, 4), np.nan))
        np.testing.assert_allclose(values, [2., 5., .3, -.2])
        self.assertTrue(np.isnan(cov).all())
        self.assertTrue(all(np.isnan(float(r["error_stat_plus_mc"])) for r in rows))
        self.assertEqual(summary["blocks"][0]["status"], "central_only_covariance_unavailable")
        epsilon = [.2, .5, .8]
        values, _, summary, _ = self.run_case(epsilon, [self.block(epsilon)], np.full((5, 5), np.nan))
        self.assertTrue(np.isnan(values[:2]).all())
        self.assertEqual(summary["blocks"][0]["status"], "unavailable_joint_covariance")

    def test_nominal_kinematics_and_override_validation(self):
        configs = [json.loads((driver.HERE / "xsec_config" / name).read_text()) for name in
                   ("xsec_config_x36_5_407.json", "xsec_config_x36_4.json")]
        nominal = driver.nominal_epsilons(configs)
        # Equivalent massless expression in terms of scattered/incident energy.
        for c, entry in zip(configs, nominal):
            ratio = c["hms_p_gev"] / c["ebeam"]
            theta = np.deg2rad(c["hms_theta_deg"])
            expected = 2*ratio*np.cos(theta/2)**2/(1+ratio*ratio+2*ratio*np.sin(theta/2)**2)
            self.assertAlmostEqual(entry["epsilon"], expected, places=14)
        self.assertGreater(nominal[0]["epsilon"], nominal[1]["epsilon"])
        self.assertEqual([r["epsilon"] for r in driver.nominal_epsilons(configs, [.2, .8])], [.2, .8])
        for invalid in ([.2], [.2, 1.], [.2, float("nan")], [.2, -.1]):
            with self.assertRaises(ValueError):
                driver.nominal_epsilons(configs, invalid)
        del configs[0]["hms_p_gev"]
        with self.assertRaisesRegex(ValueError, "hms_p_gev"):
            driver.nominal_epsilons(configs)
        self.assertEqual(driver.nominal_epsilons(configs, [.2, .8])[0]["epsilon"], .2)


if __name__ == "__main__":
    unittest.main()
