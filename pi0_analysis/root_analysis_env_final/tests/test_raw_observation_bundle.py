#!/usr/bin/env python3
"""Synthetic ALG-001 schema/ledger test; no experimental data are read."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import uproot


REPO = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO / "src" / "analysis" / "combine_analysis_branches.py"
SPEC = importlib.util.spec_from_file_location("combine_analysis_branches", MODULE_PATH)
combine = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = combine
SPEC.loader.exec_module(combine)


def raw_payload(run: int, event_ids: np.ndarray) -> dict[str, np.ndarray]:
    n = len(event_ids)
    floats = np.arange(n, dtype=np.float64) + 0.5
    return {
        "run_number": np.full(n, run, dtype=np.int32),
        "event_id": event_ids.astype(np.int64),
        "source_tree_number": np.zeros(n, dtype=np.int32),
        "source_entry": event_ids.astype(np.int64),
        "event_number": event_ids.astype(np.float64) + 1000.0,
        "t1_ns": np.full(n, 150.0, dtype=np.float64),
        "t2_ns": np.full(n, 150.0, dtype=np.float64),
        "pair_dt_ns": np.zeros(n, dtype=np.float64),
        "timing_region_mask": np.ones(n, dtype=np.uint32),
        "timing_category": np.ones(n, dtype=np.int32),
        "acquisition_mode": np.ones(n, dtype=np.int32),
        "shifted_sidebands": np.zeros(n, dtype=np.int32),
        "pair_time_diff_max_ns": np.full(n, 10.0, dtype=np.float64),
        "nclust_selected": np.full(n, 2, dtype=np.int32),
        "passes_mmiss_exclusive_cut": np.ones(n, dtype=np.int32),
        "helicity": np.ones(n, dtype=np.int32),
        "mpi0_all": np.full(n, 0.135, dtype=np.float64),
        "mmiss_all": np.full(n, 0.938, dtype=np.float64),
        "mmiss_all_corr": np.full(n, 0.938, dtype=np.float64),
        "Q2": floats + 2.0,
        "W": floats + 2.5,
        "t": -(floats + 0.1),
        "tmin": -(floats + 0.05),
        "phi": floats,
        "xB": np.full(n, 0.3, dtype=np.float64),
    }


def write_tree(root_file: uproot.WritableDirectory, name: str,
               payload: dict[str, np.ndarray]) -> None:
    branch_types = {key: value.dtype for key, value in payload.items()}
    tree = root_file.mktree(name, branch_types)
    tree.extend(payload)


def write_diagnostics(path: Path, run: int, event_ids: np.ndarray) -> None:
    with uproot.recreate(path) as root_file:
        write_tree(
            root_file,
            "physics",
            {
                "event_id": event_ids.astype(np.int32),
                "mpi0_all": np.full(len(event_ids), 0.135, dtype=np.float64),
            },
        )
        write_tree(root_file, "raw_observation", raw_payload(run, event_ids))
        write_tree(
            root_file,
            "raw_observation_segments",
            {
                "run_number": np.asarray([run], dtype=np.int32),
                "source_tree_number": np.asarray([0], dtype=np.int32),
            },
        )


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="nps-raw-observation-") as tmp:
        root_dir = Path(tmp)
        write_diagnostics(root_dir / "diagnostics_run100.root", 100,
                          np.asarray([2, 5], dtype=np.int64))
        write_diagnostics(root_dir / "diagnostics_run101.root", 101,
                          np.asarray([], dtype=np.int64))

        lookup = {
            run: combine.RunConfig("ps2", 2.0, "LH2") for run in (100, 101, 102)
        }
        efficiency = {
            run: combine.RunEfficiencyMeta(
                charge_uC=1000.0,
                tracking_eff=0.98,
                tracking_eff_err=0.01,
                hodo_3of4_eff=0.99,
                hodo_3of4_eff_err=0.01,
                livetime=0.95,
                livetime_err=0.01,
                efficiency=0.98 * 0.99,
                efficiency_err=0.02,
            )
            for run in (100, 101, 102)
        }

        raw, segments, ledger = combine.collect_raw_observation_bundle(
            lookup, root_dir, "LH2", efficiency
        )
        statuses = dict(zip(ledger["run_number"], ledger["status"]))
        assert statuses == {
            100: "ready",
            101: "zero_candidate",
            102: "missing_diagnostics",
        }
        assert list(raw[["run_number", "event_id"]].itertuples(index=False, name=None)) == [
            (100, 2), (100, 5)
        ]
        assert len(segments) == 2

        default_combined = combine.combine_branches(
            {100: lookup[100]}, root_dir, "LH2", {100: efficiency[100]}
        )
        raw_key_combined = combine.combine_branches(
            {100: lookup[100]}, root_dir, "LH2", {100: efficiency[100]},
            preserve_event_id=True,
        )
        assert "event_id" not in default_combined
        assert raw_key_combined["event_id"].tolist() == [2, 5]
        shared_columns = list(default_combined.columns)
        pd.testing.assert_frame_equal(
            default_combined,
            raw_key_combined[shared_columns],
            check_exact=True,
        )

        output = root_dir / "combined.root"
        physics = pd.DataFrame({
            "run_number": np.asarray([100, 100], dtype=np.int32),
            "event_id": np.asarray([2, 5], dtype=np.int32),
        })
        combine.save_to_root(physics, output, raw_observations=raw)
        with uproot.open(output) as root_file:
            assert root_file["physics"].num_entries == 2
            assert root_file["raw_observation"].num_entries == 2

    print("raw observation bundle synthetic test: PASS")


if __name__ == "__main__":
    main()
