#!/usr/bin/env python3
"""Focused mechanics checks for the isolated ALG-002B LH2 shadow model."""

from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from background_fit.joint_timing_mass_model import (  # noqa: E402
    JointMassFitConfig,
    _mass_layout,
    _mass_probabilities,
    build_joint_dataset,
    enforce_lh2_manifest,
    expected_lh2_runs,
    fit_and_write_joint_model,
)
from background_fit.joint_timing_model import (  # noqa: E402
    FitError,
    ObservationBundle,
    _initial_parameters,
)


def write_run_config(path: Path) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=("Kin_old", "run_number", "Type", "target")
        )
        writer.writeheader()
        writer.writerows((
            {"Kin_old": "KinC_test", "run_number": 1001,
             "Type": "production", "target": "LH2"},
            {"Kin_old": "KinC_test", "run_number": 1002,
             "Type": "Production", "target": "LH2"},
            {"Kin_old": "KinC_test", "run_number": 2001,
             "Type": "production", "target": "LD2"},
            {"Kin_old": "KinC_test", "run_number": 2002,
             "Type": "fan test", "target": "LH2"},
            {"Kin_old": "KinC_other", "run_number": 2003,
             "Type": "production", "target": "LH2"},
        ))


def synthetic_bundle(runs: tuple[int, ...] = (1001, 1002)) -> ObservationBundle:
    repeated = np.repeat(np.asarray(runs, dtype=np.int32), 8)
    n = len(repeated)
    offset = np.tile(np.linspace(-0.6, 0.6, 8), len(runs))
    arrays = {
        "run_number": repeated,
        "event_id": np.arange(n, dtype=np.int64),
        "source_tree_number": np.zeros(n, dtype=np.int32),
        "source_entry": np.arange(n, dtype=np.int64),
        "t1_ns": 150.0 + offset,
        "t2_ns": 150.0 - offset,
        "pair_dt_ns": 2.0 * offset,
        "acquisition_mode": np.full(n, 2, dtype=np.int32),
        "pair_time_diff_max_ns": np.full(n, 13.0),
        "nclust_selected": np.full(n, 2, dtype=np.int32),
        "mpi0_all": np.tile(np.linspace(0.11, 0.16, 8), len(runs)),
    }
    manifest = [
        {"path": f"synthetic://run{run}", "sha256": "synthetic",
         "bytes": 0, "entries": 8, "run_number": run}
        for run in runs
    ]
    return ObservationBundle(arrays=arrays, manifest=manifest)


def assert_rejected(bundle: ObservationBundle, expected: list[int]) -> None:
    try:
        enforce_lh2_manifest(bundle, expected, allowed_missing=())
    except FitError:
        return
    raise AssertionError("invalid LH2 run manifest was accepted")


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="nps-alg002b-test-") as temporary:
        config_path = Path(temporary) / "runs.csv"
        write_run_config(config_path)
        expected, metadata = expected_lh2_runs(config_path, "KinC_test")
        assert expected == [1001, 1002]
        assert sorted(metadata) == expected

        bundle = synthetic_bundle()
        manifest = enforce_lh2_manifest(bundle, expected, allowed_missing=())
        assert manifest["target"] == "LH2"
        assert manifest["represented_runs"] == expected

        # LD2, fan-test, wrong-kinematic, and arbitrary runs are never in the
        # expected set. Any such run entering the observations fails closed.
        assert_rejected(synthetic_bundle((1001, 2001)), expected)
        assert_rejected(synthetic_bundle((1001, 2002)), expected)
        assert_rejected(synthetic_bundle((1001, 2003)), expected)
        assert_rejected(synthetic_bundle((1001, 9999)), expected)
        assert_rejected(synthetic_bundle((1001,)), expected)

        try:
            fit_and_write_joint_model(
                bundle, REPO / "output" / ".alg002b-unit-forbidden",
                expected, config=JointMassFitConfig(starts=1),
                allowed_missing=(),
            )
        except FitError:
            pass
        else:
            raise AssertionError("canonical production output was accepted")

        fit_config = JointMassFitConfig(starts=1, coordinate_cycles=1)
        dataset = build_joint_dataset(bundle, fit_config.timing)
        parameters, layout = _mass_layout(dataset, fit_config)
        probabilities, _, penalty = _mass_probabilities(
            parameters, layout, dataset, fit_config
        )
        assert np.isfinite(penalty)
        assert len(probabilities) == 1
        assert probabilities[0].shape == (2, 6, 202)
        assert np.allclose(probabilities[0].sum(axis=2), 1.0, atol=1e-12)
        assert np.all(probabilities[0] >= 0.0)

        timing_dir = Path(temporary) / "timing"
        timing_dir.mkdir()
        for joint in dataset.strata:
            values, names, _, _ = _initial_parameters(joint.data, fit_config.timing)
            np.savez_compressed(
                timing_dir / f"{joint.tag}_profile_covariance.npz",
                parameters=values,
                parameter_names=np.asarray(names),
                run_numbers=joint.data.runs,
            )
        smoke_config = JointMassFitConfig(
            starts=1, coordinate_cycles=1, mass_maxiter=1,
            timing_refit_maxiter=1, yield_em_maxiter=30,
        )
        output = Path(temporary) / "shadow-output"
        result = fit_and_write_joint_model(
            bundle, output, expected, config=smoke_config,
            allowed_missing=(), command=["synthetic-test"],
            timing_initial_dir=timing_dir,
        )
        assert result.closure["absolute_difference"] <= 2.0e-3
        assert (output / "FIT_STATUS.txt").read_text().startswith(
            "ALG002B_SHADOW_NOT_PRODUCTION"
        )
        assert not list(Path(temporary).glob(".shadow-output.partial-*"))
        assert not any("pi0_weight" in path.name for path in output.iterdir())

    print("joint timing/mass LH2 mechanics test: PASS")


if __name__ == "__main__":
    main()
