#!/usr/bin/env python3
"""Regression for event-level raw-observation bundle comparison."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src" / "background_fit"))

from compare_raw_observation_inputs import compare_bundles  # noqa: E402


def _write(path: Path, run: int, shifted: bool = False) -> None:
    with uproot.recreate(path) as root_file:
        root_file["raw_observation"] = {
            "run_number": np.asarray([run, run], dtype=np.int32),
            "event_id": np.asarray([10, 11], dtype=np.int64),
            "t1_ns": np.asarray([149.5, 150.5 + float(shifted)], dtype=np.float64),
            "mpi0_all": ak.Array([[0.134, 0.14], [0.136]]),
        }


def main() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        reference = root / "reference"; candidate = root / "candidate"
        reference.mkdir(); candidate.mkdir()
        _write(reference / "diagnostics_run1001.root", 1001)
        _write(candidate / "diagnostics_run1001.root", 1001)
        passing = compare_bundles([str(reference)], [str(candidate)])
        assert passing["status"] == "PASS"
        _write(candidate / "diagnostics_run1001.root", 1001, shifted=True)
        failing = compare_bundles([str(reference)], [str(candidate)])
        assert failing["status"] == "FAIL"
        assert failing["failed_runs"] == [1001]
        assert failing["per_run"][0]["mismatched_branches"] == ["t1_ns"]
    print("raw-observation input comparison test: PASS")


if __name__ == "__main__":
    main()
