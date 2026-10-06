#!/usr/bin/env python3
"""Regression tests for calibrated-release binning compatibility."""

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location(
    "report_preliminary_pi0_calibrated",
    ROOT / "scripts/report_preliminary_pi0_calibrated.py",
)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


def configuration(edges):
    return {
        "configured_kinematic": "KinC_x36_4",
        "phi_bins": 12,
        "tprime_bin_edges": edges,
        "q2_bin_edges": [3.3, 4.7],
        "xb_bin_edges_by_q2": [[0.29, 0.44]],
        "diamond_xb_q2_vertices": [[0.29, 3.3], [0.44, 4.7]],
    }


class CalibratedReleaseBinningTest(unittest.TestCase):
    def test_dimensions_follow_campaign_config(self):
        edges = [-0.7, -0.289767, -0.185069, -0.103884, 0.0]
        result = REPORT.configure_binning(configuration(edges))
        self.assertEqual(result.tolist(), edges)
        self.assertEqual(REPORT.N_TPRIME_BINS, 4)
        self.assertEqual(len(REPORT.NAMES), 12)
        self.assertEqual(REPORT.NAMES[-1], "TT(bin3)")

    def test_selected_binning_mismatch_is_rejected(self):
        selected = configuration([-0.7, -0.289767, -0.185069, -0.103884, 0.0])
        source = configuration([-0.75, -0.55, -0.4, -0.25, -0.13, 0.0])
        with self.assertRaisesRegex(ValueError, "regenerate the central/toy"):
            REPORT.validate_selected_config(selected, source)

    def test_identical_selected_binning_is_accepted(self):
        selected = configuration([-0.7, -0.289767, -0.185069, -0.103884, 0.0])
        REPORT.validate_selected_config(selected, dict(selected))


if __name__ == "__main__":
    unittest.main()
