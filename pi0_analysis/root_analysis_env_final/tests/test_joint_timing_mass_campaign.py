#!/usr/bin/env python3
"""Regression for resumable ALG-002B campaign aggregation."""

from __future__ import annotations

import json
import csv
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src" / "background_fit"))

from run_joint_timing_mass_campaign import aggregate_campaign  # noqa: E402


def main() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        campaign = Path(temporary)
        campaign_config = {
            "fit_initial_dir": "/tmp/initializer", "signal_model": "dscb",
            "combinatorial_model": "bernstein3",
            "mass_gradient_backend": "autograd", "seed": 20261001,
            "coordinate_cycles": 8, "mass_maxiter": 300,
            "timing_refit_maxiter": 100, "workers_per_start": 8,
        }
        (campaign / "campaign_config.json").write_text(json.dumps(campaign_config))
        for index in range(20):
            output = campaign / f"start_{index:02d}"
            output.mkdir()
            (output / "provenance.json").write_text(json.dumps({
                "git_head": "0123456789abcdef",
                "run_manifest": {"represented_runs": [6407], "target": "LH2"},
                "fit_initial_dir": campaign_config["fit_initial_dir"],
                "output_directory": str(output.resolve()),
                "config": {
                    "starts": 1, "start_index_offset": index,
                    "signal_model": campaign_config["signal_model"],
                    "combinatorial_model": campaign_config["combinatorial_model"],
                    "mass_gradient_backend": campaign_config["mass_gradient_backend"],
                    "seed": campaign_config["seed"],
                    "coordinate_cycles": campaign_config["coordinate_cycles"],
                    "mass_maxiter": campaign_config["mass_maxiter"],
                    "timing_refit_maxiter": campaign_config["timing_refit_maxiter"],
                    "nproc": campaign_config["workers_per_start"],
                },
                "optimizer_starts": [{
                    "start_index": index,
                    "objective": 1000.0 + index * 1.0e-6,
                    "pi0_yield_sum": 500.0 + index * 1.0e-4,
                    "converged": True,
                }],
            }))
            with (output / "input_manifest.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=("path", "sha256"))
                writer.writeheader(); writer.writerow({"path": "input.root", "sha256": "abc"})
            np.savez_compressed(output / "pi0_yield_covariance.npz",
                                covariance=np.diag([25.0, 25.0]))
        summary = aggregate_campaign(campaign, 20)
        assert summary["status"] == "PASS"
        assert summary["complete_starts"] == 20
        assert summary["identity_consistent"] is True
        (campaign / "start_19" / "provenance.json").unlink()
        incomplete = aggregate_campaign(campaign, 20)
        assert incomplete["status"] == "NOT_PROMOTABLE"
        assert incomplete["complete_starts"] == 19
    print("joint timing/mass campaign aggregation test: PASS")


if __name__ == "__main__":
    main()
