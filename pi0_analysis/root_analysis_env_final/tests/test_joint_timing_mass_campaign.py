#!/usr/bin/env python3
"""Regression for resumable ALG-002B campaign aggregation."""

from __future__ import annotations

import json
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
        for index in range(20):
            output = campaign / f"start_{index:02d}"
            output.mkdir()
            (output / "provenance.json").write_text(json.dumps({
                "optimizer_starts": [{
                    "start_index": index,
                    "objective": 1000.0 + index * 1.0e-6,
                    "pi0_yield_sum": 500.0 + index * 1.0e-4,
                    "converged": True,
                }],
            }))
            np.savez_compressed(output / "pi0_yield_covariance.npz",
                                covariance=np.diag([25.0, 25.0]))
        summary = aggregate_campaign(campaign, 20)
        assert summary["status"] == "PASS"
        assert summary["complete_starts"] == 20
        (campaign / "start_19" / "provenance.json").unlink()
        incomplete = aggregate_campaign(campaign, 20)
        assert incomplete["status"] == "NOT_PROMOTABLE"
        assert incomplete["complete_starts"] == 19
    print("joint timing/mass campaign aggregation test: PASS")


if __name__ == "__main__":
    main()
