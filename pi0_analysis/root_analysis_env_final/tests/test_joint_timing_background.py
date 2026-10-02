#!/usr/bin/env python3
"""Synthetic regression checks for the isolated ALG-002A timing pilot."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from background_fit.joint_timing_model import (  # noqa: E402
    ObservationBundle,
    TimingFitConfig,
    accepted_cell_areas,
    build_strata,
    fit_observation_bundle,
    load_raw_observations,
    _zero_sum_contrast,
)


def geometry_checks() -> None:
    config = TimingFitConfig()
    unrestricted = accepted_cell_areas(config, np.inf)
    restricted = accepted_cell_areas(config, 13.0)
    assert np.isclose(unrestricted.sum(), 22.0**2, rtol=0, atol=1e-12)
    assert np.isclose(restricted.sum(), 22.0**2 - (22.0 - 13.0)**2,
                      rtol=0, atol=1e-12)

    edges = config.time_edges
    full1 = np.zeros((22, 22), dtype=bool)
    for i in range(22):
        for j in range(22):
            full1[i, j] = (edges[i] >= 155.0 and edges[i + 1] <= 161.0 and
                           edges[j] >= 139.0 and edges[j + 1] <= 145.0)
    assert np.isclose(unrestricted.reshape(22, 22)[full1].sum(), 36.0,
                      rtol=0, atol=1e-12)
    assert np.isclose(restricted.reshape(22, 22)[full1].sum(), 4.5,
                      rtol=0, atol=1e-12)

    contrast = _zero_sum_contrast(56)
    assert contrast.shape == (56, 55)
    assert np.allclose(contrast.sum(axis=0), 0.0, atol=1e-12)
    assert np.allclose(contrast.T @ contrast, np.eye(55), atol=1e-12)


def synthetic_bundle(seed: int = 20261001) -> ObservationBundle:
    rng = np.random.default_rng(seed)
    rows: dict[str, list[np.ndarray]] = {
        name: [] for name in (
            "run_number", "event_id", "source_tree_number", "source_entry",
            "t1_ns", "t2_ns", "pair_dt_ns", "acquisition_mode",
            "pair_time_diff_max_ns", "nclust_selected", "mpi0_all",
        )
    }
    event_base = 0
    for run, shift, count in ((1001, -0.35, 900), (1002, 0.35, 1100)):
        values: list[tuple[float, float]] = []
        while len(values) < count:
            n = max(256, count - len(values))
            component = rng.choice(5, size=n, p=(0.58, 0.10, 0.10, 0.17, 0.05))
            t1 = rng.uniform(139.0, 161.0, size=n)
            t2 = rng.uniform(139.0, 161.0, size=n)
            central1 = rng.normal(150.0 + shift, 0.72, size=n)
            central2 = (150.0 + shift + 0.25 * (central1 - (150.0 + shift)) +
                        rng.normal(0.0, 0.70, size=n))
            t1[component == 0] = central1[component == 0]
            t2[component == 0] = central2[component == 0]
            t1[component == 1] = central1[component == 1]
            t2[component == 2] = central2[component == 2]
            satellites = rng.choice((140.0, 142.0, 144.0, 156.0, 158.0, 160.0), size=n)
            satellite1 = rng.normal(satellites + shift, 0.65)
            satellite2 = satellites + shift + 0.55 * (satellite1 - satellites - shift) + rng.normal(0, 0.55, size=n)
            t1[component == 4] = satellite1[component == 4]
            t2[component == 4] = satellite2[component == 4]
            accepted = ((t1 >= 139.0) & (t1 < 161.0) &
                        (t2 >= 139.0) & (t2 < 161.0) &
                        (np.abs(t1 - t2) <= 13.0))
            values.extend(zip(t1[accepted], t2[accepted]))
        pair = np.asarray(values[:count])
        event_ids = np.arange(event_base, event_base + count, dtype=np.int64)
        event_base += count
        payload = {
            "run_number": np.full(count, run, dtype=np.int32),
            "event_id": event_ids,
            "source_tree_number": np.zeros(count, dtype=np.int32),
            "source_entry": np.arange(count, dtype=np.int64),
            "t1_ns": pair[:, 0],
            "t2_ns": pair[:, 1],
            "pair_dt_ns": pair[:, 0] - pair[:, 1],
            "acquisition_mode": np.full(count, 2, dtype=np.int32),
            "pair_time_diff_max_ns": np.full(count, 13.0),
            "nclust_selected": np.full(count, 3, dtype=np.int32),
            # Include both flow bins so exact closure is exercised.
            "mpi0_all": np.concatenate((
                rng.normal(0.135, 0.012, size=count - 2),
                np.asarray([-0.01, 0.45]),
            )),
        }
        for name, value in payload.items():
            rows[name].append(np.asarray(value))
    arrays = {name: np.concatenate(parts) for name, parts in rows.items()}
    return ObservationBundle(arrays=arrays, manifest=[{
        "path": "synthetic://joint-timing-test",
        "sha256": "synthetic",
        "bytes": 0,
        "entries": len(arrays["run_number"]),
        "run_number": "1001,1002",
    }])


def write_loader_fixture(path: Path) -> None:
    payload = synthetic_bundle().arrays
    selected = payload["run_number"] == 1001
    with uproot.recreate(path) as root_file:
        data = {name: values[selected] for name, values in payload.items()}
        tree = root_file.mktree("raw_observation", {
            name: values.dtype for name, values in data.items()
        })
        tree.extend(data)


def main() -> None:
    geometry_checks()
    bundle = synthetic_bundle()
    config = TimingFitConfig(optimizer_maxiter=180, optimizer_tolerance=2e-7)
    strata = build_strata(bundle, config)
    assert len(strata) == 1
    assert strata[0].mode == 2 and strata[0].multiplicity == 3
    assert int(strata[0].event_count_by_run.sum()) == 2000

    with tempfile.TemporaryDirectory(prefix="nps-alg002a-test-") as temporary:
        temporary = Path(temporary)
        fixture = temporary / "diagnostics_run1001.root"
        write_loader_fixture(fixture)
        loaded = load_raw_observations([fixture])
        assert len(loaded.arrays["run_number"]) == 900

        output = temporary / "shadow-output"
        fits = fit_observation_bundle(bundle, output, config=config,
                                      command=["synthetic-test"])
        assert len(fits) == 1
        offsets = [row["timing_offset_ns"] for row in fits[0].run_parameter_rows]
        assert np.isclose(sum(offsets), 0.0, atol=1e-10)
        assert offsets[0] < offsets[1]
        assert sum(int(row["observed_count"]) for row in fits[0].component_rows) == 2000
        assert (output / "FIT_STATUS.txt").read_text().startswith(
            "ALG002A_SHADOW_NOT_PRODUCTION"
        )
        assert not any("pi0_weight" in path.name for path in output.iterdir())
        provenance = (output / "provenance.json").read_text()
        assert '"pi0_weight_written": false' in provenance
        assert '"efficiency_inputs_read": false' in provenance

    print("joint timing background synthetic test: PASS")


if __name__ == "__main__":
    main()
