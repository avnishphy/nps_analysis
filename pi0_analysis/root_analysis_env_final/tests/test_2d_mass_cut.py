#!/usr/bin/env python3
"""Regression tests for the combined Python 2D mass-cut implementation."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO / "src" / "analysis" / "combine_analysis_branches.py"
SPEC = importlib.util.spec_from_file_location("combine_analysis_branches", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class MassCutRegressionTest(unittest.TestCase):
    def test_remote_background_mode_cannot_capture_fit(self) -> None:
        rng = np.random.default_rng(7)
        n_signal = 12_000
        n_background = 45_000
        covariance = np.array([
            [0.0035 ** 2, -0.82 * 0.0035 * 0.065],
            [-0.82 * 0.0035 * 0.065, 0.065 ** 2],
        ])
        signal = rng.multivariate_normal([0.13498, 0.938], covariance, n_signal)
        background_x = rng.normal(0.124, 0.006, n_background)
        background_y = 1.32 - 18.0 * (background_x - 0.124)
        background_y += rng.normal(0.0, 0.10, n_background)

        frame = pd.DataFrame({
            "mpi0_all": np.concatenate((signal[:, 0], background_x)),
            "mmiss_all": np.concatenate((signal[:, 1], background_y)),
            "pi0_weight": np.ones(n_signal + n_background),
            "scale": np.ones(n_signal + n_background),
        })
        debug = MODULE.add_combined_2d_mass_cut(frame)
        self.assertIsNotNone(debug)
        params = debug["params"]

        cfg = MODULE.MASS_CUT_CONFIG
        self.assertLessEqual(
            abs(params["mean_mpi0"] - cfg["seed_mpi0"]),
            cfg["max_model_mpi0_offset"],
        )
        self.assertLessEqual(
            abs(params["mean_mmiss"] - cfg["seed_mmiss"]),
            cfg["max_model_mmiss_offset"],
        )
        self.assertLessEqual(
            params["fit_subset_total_fraction"],
            cfg["auto_max_core_total_fraction"],
        )
        self.assertGreater(params["fit_subset_total_fraction"], 0.10)
        self.assertGreater(params["core_total_fraction"], 0.0)
        self.assertLess(params["core_total_fraction"], params["ellipse_total_fraction"])
        self.assertGreater(params["ellipse_growth_steps"], 0.0)
        self.assertGreater(frame.iloc[:n_signal]["is_exclusive_ellipse_combined"].mean(), 0.95)
        self.assertLess(frame.iloc[n_signal:]["is_exclusive_ellipse_combined"].mean(), 0.15)
        self.assertEqual(params["mcd_valid"], 1.0)
        self.assertLessEqual(
            abs(params["mcd_mean_mmiss"] - cfg["seed_mmiss"]),
            cfg["max_model_mmiss_offset"],
        )

    def test_covariance_regularization_handles_collinear_bins(self) -> None:
        y = np.linspace(0.88, 1.00, 40)
        model = MODULE.compute_cov_model(
            np.full_like(y, 0.135),
            y,
            np.ones_like(y),
            regularization_x=0.000125,
            regularization_y=0.00225,
        )
        self.assertIsNotNone(model)
        self.assertGreater(model["det"], 0.0)


if __name__ == "__main__":
    unittest.main()
