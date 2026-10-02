#!/usr/bin/env python3
"""ALG-002A combined-statistics timing-background pilot.

This module is deliberately isolated from the production analysis.  It reads
integer ALG-001 raw observations, retains run/mode/multiplicity identity, and
fits an extended-Poisson mixture on the accepted two-photon timing lattice.
It does not calculate pi0 purity weights, efficiencies, or cross sections.

The timing fit is performed separately for each acquisition-mode/multiplicity
stratum.  Timing shapes are shared, run centers and widths have constrained
deviations, and the five component yields remain free for every run.  Once the
timing shapes are fitted, every run/mass bin is decomposed independently.  The
true-coincidence mass spectrum is intentionally unparameterized.

Scientific status: shadow validation only; not a publication estimator.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import warnings
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import uproot

with warnings.catch_warnings():
    warnings.simplefilter("ignore", UserWarning)
    from scipy.interpolate import BSpline
    from scipy.optimize import minimize


COMPONENTS = ("true", "horizontal", "vertical", "random", "diagonal")
REQUIRED_BRANCHES = (
    "run_number",
    "event_id",
    "source_tree_number",
    "source_entry",
    "t1_ns",
    "t2_ns",
    "pair_dt_ns",
    "acquisition_mode",
    "pair_time_diff_max_ns",
    "nclust_selected",
    "mpi0_all",
)


class FitError(RuntimeError):
    """Fail-closed input, identifiability, or optimizer error."""


@dataclass(frozen=True)
class TimingFitConfig:
    time_low_ns: float = 139.0
    time_high_ns: float = 161.0
    time_bin_width_ns: float = 1.0
    mass_low_gev: float = 0.0
    mass_high_gev: float = 0.4
    mass_bin_width_gev: float = 0.002
    spline_knots_ns: tuple[float, ...] = (
        139.0, 143.0, 147.0, 149.0, 151.0, 153.0, 157.0, 161.0
    )
    spline_penalty: float = 10.0
    optimizer_maxiter: int = 500
    optimizer_tolerance: float = 1.0e-8
    em_maxiter: int = 300
    em_tolerance: float = 1.0e-8
    minimum_expected: float = 1.0e-12
    run_offset_min_ns: float = -2.0
    run_offset_max_ns: float = 2.0
    run_log_width_min: float = -0.7
    run_log_width_max: float = 0.7
    central_sigma_min_ns: float = 0.15
    central_sigma_max_ns: float = 3.0
    diagonal_sigma_min_ns: float = 0.15
    diagonal_sigma_max_ns: float = 2.5

    @property
    def time_edges(self) -> np.ndarray:
        return np.arange(
            self.time_low_ns,
            self.time_high_ns + 0.5 * self.time_bin_width_ns,
            self.time_bin_width_ns,
            dtype=float,
        )

    @property
    def mass_edges(self) -> np.ndarray:
        return np.arange(
            self.mass_low_gev,
            self.mass_high_gev + 0.5 * self.mass_bin_width_gev,
            self.mass_bin_width_gev,
            dtype=float,
        )


@dataclass
class ObservationBundle:
    arrays: dict[str, np.ndarray]
    manifest: list[dict[str, object]]


@dataclass
class StratumData:
    mode: int
    multiplicity: int
    runs: np.ndarray
    run_index: np.ndarray
    cell_index: np.ndarray
    mass_index: np.ndarray
    counts_by_run_cell: np.ndarray
    pair_cut_by_run: np.ndarray
    event_count_by_run: np.ndarray


@dataclass
class StratumFit:
    mode: int
    multiplicity: int
    runs: np.ndarray
    pair_cuts: np.ndarray
    parameters: np.ndarray
    parameter_names: list[str]
    probabilities: np.ndarray
    integrated_yields: np.ndarray
    nll: float
    optimizer_success: bool
    optimizer_message: str
    optimizer_iterations: int
    profile_covariance: np.ndarray
    covariance_eigenvalues: np.ndarray
    covariance_condition: float
    component_rows: list[dict[str, object]]
    prediction_rows: list[dict[str, object]]
    run_parameter_rows: list[dict[str, object]]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FitError(message)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_raw_observations(paths: Sequence[Path | str]) -> ObservationBundle:
    """Load and validate ALG-001 trees without modifying source files."""
    normalized = sorted({Path(path).resolve() for path in paths})
    _require(bool(normalized), "no input ROOT files supplied")
    chunks: dict[str, list[np.ndarray]] = {name: [] for name in REQUIRED_BRANCHES}
    manifest: list[dict[str, object]] = []
    seen_keys: set[tuple[int, int]] = set()

    for path in normalized:
        _require(path.is_file(), f"missing input ROOT file: {path}")
        with uproot.open(path) as root_file:
            _require("raw_observation" in root_file,
                     f"{path} has no ALG-001 raw_observation tree")
            tree = root_file["raw_observation"]
            missing = sorted(set(REQUIRED_BRANCHES) - set(tree.keys()))
            _require(not missing, f"{path} missing branches: {missing}")
            data = tree.arrays(REQUIRED_BRANCHES, library="np")
            entries = int(tree.num_entries)
            runs = np.unique(data["run_number"].astype(np.int64))
            _require(len(runs) <= 1,
                     f"per-run diagnostics input contains multiple runs: {path}")
            run_value = int(runs[0]) if len(runs) else None
            for event_id in data["event_id"].astype(np.int64):
                key = (run_value if run_value is not None else -1, int(event_id))
                _require(key not in seen_keys, f"duplicate (run,event_id) key {key}")
                seen_keys.add(key)
            for name in REQUIRED_BRANCHES:
                chunks[name].append(np.asarray(data[name]))
            manifest.append({
                "path": str(path),
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
                "entries": entries,
                "run_number": run_value,
            })

    arrays = {
        name: np.concatenate(parts) if parts else np.asarray([], dtype=float)
        for name, parts in chunks.items()
    }
    n = len(arrays["run_number"])
    _require(n > 0, "raw-observation inputs contain zero selected events")
    _require(all(len(value) == n for value in arrays.values()),
             "raw-observation branch lengths differ")
    for name in ("t1_ns", "t2_ns", "pair_dt_ns", "pair_time_diff_max_ns", "mpi0_all"):
        _require(np.all(np.isfinite(arrays[name])), f"nonfinite values in {name}")
    _require(np.all(np.isin(arrays["acquisition_mode"], (1, 2))),
             "acquisition_mode must be 1 (HCANA) or 2 (waveform)")
    _require(np.all(arrays["nclust_selected"] >= 2),
             "raw observation with fewer than two selected clusters")
    return ObservationBundle(arrays=arrays, manifest=manifest)


def _clip_polygon(
    polygon: list[tuple[float, float]], a: float, b: float, c: float
) -> list[tuple[float, float]]:
    """Clip a convex polygon to a*x+b*y <= c."""
    if not polygon:
        return []
    output: list[tuple[float, float]] = []
    previous = polygon[-1]
    previous_value = a * previous[0] + b * previous[1] - c
    for current in polygon:
        current_value = a * current[0] + b * current[1] - c
        previous_inside = previous_value <= 1.0e-12
        current_inside = current_value <= 1.0e-12
        if previous_inside != current_inside:
            fraction = previous_value / (previous_value - current_value)
            output.append((
                previous[0] + fraction * (current[0] - previous[0]),
                previous[1] + fraction * (current[1] - previous[1]),
            ))
        if current_inside:
            output.append(current)
        previous, previous_value = current, current_value
    return output


def _polygon_area(polygon: Sequence[tuple[float, float]]) -> float:
    if len(polygon) < 3:
        return 0.0
    x = np.asarray([point[0] for point in polygon])
    y = np.asarray([point[1] for point in polygon])
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _accepted_polygon(
    x0: float, x1: float, y0: float, y1: float, pair_cut: float
) -> list[tuple[float, float]]:
    polygon = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    if math.isfinite(pair_cut):
        polygon = _clip_polygon(polygon, 1.0, -1.0, pair_cut)
        polygon = _clip_polygon(polygon, -1.0, 1.0, pair_cut)
    return polygon


def accepted_cell_areas(config: TimingFitConfig, pair_cut: float) -> np.ndarray:
    """Exact geometric accepted area in every 1 ns lattice cell."""
    edges = config.time_edges
    nbin = len(edges) - 1
    areas = np.zeros(nbin * nbin, dtype=float)
    for i in range(nbin):
        for j in range(nbin):
            polygon = _accepted_polygon(edges[i], edges[i + 1],
                                        edges[j], edges[j + 1], pair_cut)
            areas[i * nbin + j] = _polygon_area(polygon)
    return areas


def _cell_quadrature(
    config: TimingFitConfig, pair_cut: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Degree-two triangle quadrature over exact accepted cell polygons."""
    edges = config.time_edges
    nbin = len(edges) - 1
    xs: list[float] = []
    ys: list[float] = []
    weights: list[float] = []
    cells: list[int] = []
    barycentric = ((2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
                   (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
                   (1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0))
    for i in range(nbin):
        for j in range(nbin):
            polygon = _accepted_polygon(edges[i], edges[i + 1],
                                        edges[j], edges[j + 1], pair_cut)
            if len(polygon) < 3:
                continue
            anchor = polygon[0]
            for k in range(1, len(polygon) - 1):
                triangle = (anchor, polygon[k], polygon[k + 1])
                area = _polygon_area(triangle)
                if area <= 0.0:
                    continue
                for l0, l1, l2 in barycentric:
                    xs.append(l0 * triangle[0][0] + l1 * triangle[1][0] + l2 * triangle[2][0])
                    ys.append(l0 * triangle[0][1] + l1 * triangle[1][1] + l2 * triangle[2][1])
                    weights.append(area / 3.0)
                    cells.append(i * nbin + j)
    return (np.asarray(xs), np.asarray(ys), np.asarray(weights),
            np.asarray(cells, dtype=np.int64))


def _multiplicity_class(values: np.ndarray) -> np.ndarray:
    return np.where(values == 2, 2, 3).astype(np.int32)


def _mass_indices(values: np.ndarray, config: TimingFitConfig) -> np.ndarray:
    """Return 0=underflow, 1..N regular bins, N+1=overflow."""
    edges = config.mass_edges
    regular = np.searchsorted(edges, values, side="right") - 1
    output = regular + 1
    output[values < edges[0]] = 0
    output[values >= edges[-1]] = len(edges)
    return output.astype(np.int32)


def build_strata(bundle: ObservationBundle, config: TimingFitConfig) -> list[StratumData]:
    arrays = bundle.arrays
    t1 = arrays["t1_ns"].astype(float)
    t2 = arrays["t2_ns"].astype(float)
    edges = config.time_edges
    inside = ((t1 >= edges[0]) & (t1 < edges[-1]) &
              (t2 >= edges[0]) & (t2 < edges[-1]))
    _require(np.all(inside),
             f"{np.count_nonzero(~inside)} observations outside approved timing domain "
             f"[{edges[0]}, {edges[-1]}) ns")
    i = np.searchsorted(edges, t1, side="right") - 1
    j = np.searchsorted(edges, t2, side="right") - 1
    ntime = len(edges) - 1
    cell = (i * ntime + j).astype(np.int32)
    mass = _mass_indices(arrays["mpi0_all"].astype(float), config)
    mult = _multiplicity_class(arrays["nclust_selected"].astype(int))
    strata: list[StratumData] = []

    for mode in sorted(np.unique(arrays["acquisition_mode"].astype(int))):
        for multiplicity in (2, 3):
            selected = ((arrays["acquisition_mode"] == mode) & (mult == multiplicity))
            if not np.any(selected):
                continue
            runs = np.unique(arrays["run_number"][selected].astype(np.int64))
            run_map = {int(run): index for index, run in enumerate(runs)}
            run_index = np.asarray(
                [run_map[int(run)] for run in arrays["run_number"][selected]], dtype=np.int32
            )
            cuts = np.full(len(runs), np.inf, dtype=float)
            event_counts = np.bincount(run_index, minlength=len(runs)).astype(np.int64)
            counts = np.zeros((len(runs), ntime * ntime), dtype=np.int64)
            np.add.at(counts, (run_index, cell[selected]), 1)
            if multiplicity == 3:
                for run_index_value, run in enumerate(runs):
                    q = selected & (arrays["run_number"] == run)
                    unique_cuts = np.unique(np.round(
                        arrays["pair_time_diff_max_ns"][q].astype(float), 9
                    ))
                    _require(len(unique_cuts) == 1,
                             f"run {run} has multiple pair-time cuts in one stratum")
                    cuts[run_index_value] = float(unique_cuts[0])
                    violation = np.abs(arrays["pair_dt_ns"][q]) > cuts[run_index_value] + 1.0e-9
                    _require(not np.any(violation),
                             f"run {run} contains observations outside recorded pair-time cut")
            strata.append(StratumData(
                mode=int(mode),
                multiplicity=multiplicity,
                runs=runs,
                run_index=run_index,
                cell_index=cell[selected],
                mass_index=mass[selected],
                counts_by_run_cell=counts,
                pair_cut_by_run=cuts,
                event_count_by_run=event_counts,
            ))
    _require(sum(int(stratum.event_count_by_run.sum()) for stratum in strata) == len(t1),
             "stratum construction did not preserve every observation")
    return strata


def _bivariate_normal(
    x: np.ndarray,
    y: np.ndarray,
    mean: float,
    sigma: float,
    rho: float,
) -> np.ndarray:
    sx = (x - mean) / sigma
    sy = (y - mean) / sigma
    one_minus = max(1.0e-8, 1.0 - rho * rho)
    exponent = -0.5 * (sx * sx - 2.0 * rho * sx * sy + sy * sy) / one_minus
    return np.exp(exponent) / (2.0 * math.pi * sigma * sigma * math.sqrt(one_minus))


def _bivariate_normal_log(
    x: np.ndarray,
    y: np.ndarray,
    mean: float,
    sigma: float,
    rho: float,
) -> np.ndarray:
    sx = (x - mean) / sigma
    sy = (y - mean) / sigma
    one_minus = max(1.0e-12, 1.0 - rho * rho)
    return (-0.5 * (sx * sx - 2.0 * rho * sx * sy + sy * sy) / one_minus -
            math.log(2.0 * math.pi * sigma * sigma * math.sqrt(one_minus)))


def _normal(x: np.ndarray, mean: float, sigma: float) -> np.ndarray:
    z = (x - mean) / sigma
    return np.exp(-0.5 * z * z) / (math.sqrt(2.0 * math.pi) * sigma)


def _normal_log(x: np.ndarray, mean: float, sigma: float) -> np.ndarray:
    z = (x - mean) / sigma
    return -0.5 * z * z - math.log(math.sqrt(2.0 * math.pi) * sigma)


def _relative_density(log_density: np.ndarray) -> np.ndarray:
    """Exponentiate after removing an irrelevant component-wide log scale."""
    maximum = float(np.max(log_density))
    _require(np.isfinite(maximum), "component log density has no finite support")
    return np.exp(log_density - maximum)


def _spline_basis(values: np.ndarray, config: TimingFitConfig) -> np.ndarray:
    knots = config.spline_knots_ns
    degree = 3
    full_knots = np.asarray(
        [knots[0]] * (degree + 1) + list(knots[1:-1]) + [knots[-1]] * (degree + 1),
        dtype=float,
    )
    # SciPy 1.9 on the Hall C analysis hosts predates the explicit
    # ``extrapolate`` keyword.  Every caller has already validated that values
    # lie inside the clamped boundary knots, so the older API is sufficient.
    return BSpline.design_matrix(values, full_knots, degree).toarray()


@lru_cache(maxsize=32)
def _cached_quadrature_basis(
    config: TimingFitConfig, pair_cut: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Cache immutable geometry and spline design shared by optimizer calls."""
    t1, t2, weights, cells = _cell_quadrature(config, pair_cut)
    return t1, t2, weights, cells, _spline_basis(t1, config), _spline_basis(t2, config)


def _parameter_layout(nrun: int, nbasis: int) -> tuple[list[str], dict[str, slice | int]]:
    names = ["central_mean_ns", "log_central_sigma_ns", "atanh_central_rho",
             "log_diagonal_sigma_ns", "atanh_diagonal_rho"]
    names.extend(f"random_spline_logcoef_{i}" for i in range(1, nbasis))
    names.extend(("log_run_offset_scale_ns", "log_run_width_scale"))
    offset_start = len(names)
    names.extend(f"run_offset_free_{i}" for i in range(max(0, nrun - 1)))
    width_start = len(names)
    names.extend(f"run_log_width_free_{i}" for i in range(max(0, nrun - 1)))
    layout: dict[str, slice | int] = {
        "center": 0,
        "log_sigma": 1,
        "rho": 2,
        "log_diag_sigma": 3,
        "diag_rho": 4,
        "spline": slice(5, 5 + nbasis - 1),
        "log_tau_offset": 5 + nbasis - 1,
        "log_tau_width": 6 + nbasis - 1,
        "offsets": slice(offset_start, offset_start + max(0, nrun - 1)),
        "widths": slice(width_start, width_start + max(0, nrun - 1)),
    }
    return names, layout


@lru_cache(maxsize=32)
def _zero_sum_contrast(nrun: int) -> np.ndarray:
    """Balanced orthonormal basis for the nrun-dimensional zero-sum space."""
    if nrun <= 1:
        return np.zeros((nrun, 0), dtype=float)
    columns: list[np.ndarray] = []

    def split(indices: np.ndarray) -> None:
        if len(indices) <= 1:
            return
        midpoint = len(indices) // 2
        left, right = indices[:midpoint], indices[midpoint:]
        nl, nr = len(left), len(right)
        vector = np.zeros(nrun, dtype=float)
        vector[left] = math.sqrt(nr / (nl * (nl + nr)))
        vector[right] = -math.sqrt(nl / (nr * (nl + nr)))
        columns.append(vector)
        split(left)
        split(right)

    split(np.arange(nrun, dtype=int))
    contrast = np.column_stack(columns)
    _require(contrast.shape == (nrun, nrun - 1), "contrast construction failed")
    _require(np.allclose(contrast.T @ contrast, np.eye(nrun - 1), atol=1.0e-12),
             "run contrast is not orthonormal")
    return contrast


def _zero_sum_deviations(free: np.ndarray, nrun: int) -> np.ndarray:
    if nrun == 1:
        return np.zeros(1, dtype=float)
    _require(len(free) == nrun - 1, "zero-sum contrast dimension mismatch")
    # Balanced orthonormal contrasts span the zero-sum run-deviation subspace.
    # The earlier [free, -sum(free)] basis singled out the final run and could
    # force it to absorb a large compensating offset.  This tree basis has unit
    # norm columns and O(log n) average support, retaining efficient local
    # finite differences without privileging one run.
    deviations = _zero_sum_contrast(nrun) @ np.asarray(free, dtype=float)
    _require(abs(float(np.sum(deviations))) < 1.0e-10,
             "Helmert run deviations do not sum to zero")
    return deviations


def _component_probabilities(
    parameters: np.ndarray,
    stratum: StratumData,
    config: TimingFitConfig,
    layout: dict[str, slice | int],
    run_indices: Sequence[int] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    nrun = len(stratum.runs)
    center = float(parameters[layout["center"]])
    base_sigma = math.exp(float(parameters[layout["log_sigma"]]))
    rho = math.tanh(float(parameters[layout["rho"]]))
    diag_sigma = math.exp(float(parameters[layout["log_diag_sigma"]]))
    diag_rho = math.tanh(float(parameters[layout["diag_rho"]]))
    spline_free = parameters[layout["spline"]]
    spline_coeff = np.concatenate(([0.0], spline_free))
    spline_coeff -= np.mean(spline_coeff)
    spline_positive = np.exp(np.clip(spline_coeff, -20.0, 20.0))
    offsets = _zero_sum_deviations(parameters[layout["offsets"]], nrun)
    width_deviations = _zero_sum_deviations(parameters[layout["widths"]], nrun)
    tau_offset = math.exp(float(parameters[layout["log_tau_offset"]]))
    tau_width = math.exp(float(parameters[layout["log_tau_width"]]))
    ncell = stratum.counts_by_run_cell.shape[1]
    selected_runs = list(range(nrun)) if run_indices is None else list(run_indices)
    probabilities = np.zeros((len(selected_runs), len(COMPONENTS), ncell), dtype=float)
    satellite_centers = np.asarray((140.0, 142.0, 144.0, 156.0, 158.0, 160.0))

    random_cache: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                                    np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    for output_index, r in enumerate(selected_runs):
        cut = float(stratum.pair_cut_by_run[r])
        cache_key = cut if math.isfinite(cut) else math.inf
        t1, t2, qweight, qcell, basis1, basis2 = _cached_quadrature_basis(
            config, cache_key
        )
        if cache_key not in random_cache:
            random_cache[cache_key] = (
                t1, t2, qweight, qcell, basis1, basis2,
                basis1 @ spline_positive, basis2 @ spline_positive,
            )
        t1, t2, qweight, qcell, basis1, basis2, u1, u2 = random_cache[cache_key]
        run_center = center + offsets[r]
        run_sigma = base_sigma * math.exp(width_deviations[r])
        log_g1 = _normal_log(t1, run_center, run_sigma)
        log_g2 = _normal_log(t2, run_center, run_sigma)
        log_u1 = np.log(np.maximum(u1, np.finfo(float).tiny))
        log_u2 = np.log(np.maximum(u2, np.finfo(float).tiny))
        diagonal_logs = np.asarray([
            _bivariate_normal_log(
                t1, t2, value + offsets[r],
                diag_sigma * math.exp(width_deviations[r]), diag_rho
            )
            for value in satellite_centers
        ])
        diagonal_max = np.max(diagonal_logs, axis=0)
        diagonal_logsum = (diagonal_max +
                           np.log(np.sum(np.exp(diagonal_logs - diagonal_max), axis=0)) -
                           math.log(len(satellite_centers)))
        densities = [
            _relative_density(_bivariate_normal_log(
                t1, t2, run_center, run_sigma, rho
            )),
            _relative_density(log_g1 + log_u2),
            _relative_density(log_u1 + log_g2),
            _relative_density(log_u1 + log_u2),
            _relative_density(diagonal_logsum),
        ]
        for component, density in enumerate(densities):
            integrated = np.bincount(
                qcell, weights=qweight * density, minlength=ncell
            ).astype(float)
            total = float(integrated.sum())
            _require(total > 0.0 and np.isfinite(total),
                     f"nonpositive component normalization for run {stratum.runs[r]}")
            probabilities[output_index, component] = integrated / total
    return probabilities, offsets, width_deviations, tau_offset, tau_width


def _profile_one_run_yields(
    counts: np.ndarray,
    probabilities: np.ndarray,
    config: TimingFitConfig,
    initial: np.ndarray | None = None,
) -> tuple[np.ndarray, float, int]:
    total = float(np.sum(counts))
    if total == 0.0:
        return np.zeros(len(COMPONENTS)), 0.0, 0
    if initial is None:
        yields = total * np.asarray((0.60, 0.10, 0.10, 0.15, 0.05))
    else:
        yields = np.maximum(np.asarray(initial, dtype=float), config.minimum_expected)
        yields *= total / yields.sum()
    occupied = counts > 0
    p = probabilities[:, occupied]
    n = counts[occupied].astype(float)
    converged = False
    for iteration in range(config.em_maxiter):
        terms = yields[:, None] * p
        denominator = np.maximum(np.sum(terms, axis=0), config.minimum_expected)
        updated = np.sum(n[None, :] * terms / denominator[None, :], axis=1)
        delta = np.max(np.abs(updated - yields) / np.maximum(1.0, yields))
        yields = np.maximum(updated, 0.0)
        if delta < config.em_tolerance:
            converged = True
            break
    optimizer_iterations = 0
    if not converged:
        # The five-yield problem is convex for fixed component probabilities.
        # EM can nevertheless converge very slowly when two timing components
        # are nearly collinear.  Use EM only as an initializer, then solve the
        # same extended-Poisson likelihood with its exact analytic gradient.
        def yield_objective(candidate: np.ndarray) -> float:
            mu_candidate = np.maximum(candidate @ p, config.minimum_expected)
            return float(np.sum(candidate) - np.dot(n, np.log(mu_candidate)))

        def yield_gradient(candidate: np.ndarray) -> np.ndarray:
            mu_candidate = np.maximum(candidate @ p, config.minimum_expected)
            return 1.0 - p @ (n / mu_candidate)

        yield_result = minimize(
            yield_objective,
            yields,
            jac=yield_gradient,
            method="L-BFGS-B",
            bounds=[(0.0, None)] * len(COMPONENTS),
            options={"maxiter": 1000, "ftol": 1.0e-12, "gtol": 1.0e-9},
        )
        candidate = np.maximum(np.asarray(yield_result.x), 0.0)

        def kkt_residual(candidate_yields: np.ndarray) -> float:
            gradient = yield_gradient(candidate_yields)
            positive_yields = candidate_yields > 1.0e-10 * max(1.0, total)
            residuals = []
            if np.any(positive_yields):
                residuals.append(float(np.max(np.abs(gradient[positive_yields]))))
            if np.any(~positive_yields):
                residuals.append(float(np.max(np.maximum(-gradient[~positive_yields], 0.0))))
            return max(residuals, default=0.0)

        residual = kkt_residual(candidate)
        optimizer_iterations = int(yield_result.nit)
        if not yield_result.success and residual > 2.0e-6:
            second_result = minimize(
                yield_objective,
                candidate,
                jac=yield_gradient,
                method="SLSQP",
                bounds=[(0.0, None)] * len(COMPONENTS),
                options={"maxiter": 1000, "ftol": 1.0e-12, "disp": False},
            )
            second_candidate = np.maximum(np.asarray(second_result.x), 0.0)
            second_residual = kkt_residual(second_candidate)
            if (yield_objective(second_candidate) <= yield_objective(candidate) or
                    second_residual < residual):
                candidate, residual = second_candidate, second_residual
            optimizer_iterations += int(second_result.nit)
            if not second_result.success and residual > 2.0e-6:
                raise FitError(
                    "component-yield convex fallbacks failed KKT check: "
                    f"lbfgs={yield_result.message}; slsqp={second_result.message}; "
                    f"residual={residual:.6g}"
                )
        yields = candidate
    mu = np.maximum(yields @ probabilities, config.minimum_expected)
    nll = float(np.sum(yields) - np.dot(counts[occupied], np.log(mu[occupied])))
    return yields, nll, iteration + 1 + optimizer_iterations


def _profile_all_yields(
    counts: np.ndarray, probabilities: np.ndarray, config: TimingFitConfig
) -> tuple[np.ndarray, float]:
    yields = np.zeros((len(counts), len(COMPONENTS)), dtype=float)
    nll = 0.0
    for r in range(len(counts)):
        yields[r], run_nll, _ = _profile_one_run_yields(counts[r], probabilities[r], config)
        nll += run_nll
    return yields, nll


def _shape_penalty(
    parameters: np.ndarray,
    layout: dict[str, slice | int],
    nrun: int,
    config: TimingFitConfig,
) -> float:
    spline_free = parameters[layout["spline"]]
    spline = np.concatenate(([0.0], spline_free))
    spline -= np.mean(spline)
    penalty = 0.5 * config.spline_penalty * float(np.dot(np.diff(spline, n=2), np.diff(spline, n=2)))
    offsets = _zero_sum_deviations(parameters[layout["offsets"]], nrun)
    widths = _zero_sum_deviations(parameters[layout["widths"]], nrun)
    tau_offset = math.exp(float(parameters[layout["log_tau_offset"]]))
    tau_width = math.exp(float(parameters[layout["log_tau_width"]]))
    if nrun > 1:
        penalty += 0.5 * float(np.dot(offsets, offsets)) / (tau_offset * tau_offset)
        penalty += (nrun - 1) * math.log(tau_offset)
        penalty += 0.5 * float(np.dot(widths, widths)) / (tau_width * tau_width)
        penalty += (nrun - 1) * math.log(tau_width)
    # Weak, explicit hyper-scale regularization prevents an unidentifiable
    # variance component from running to a numerical boundary.
    penalty += 0.5 * (tau_offset / 2.0) ** 2 + 0.5 * (tau_width / 0.5) ** 2
    return penalty


def _initial_parameters(
    stratum: StratumData, config: TimingFitConfig
) -> tuple[np.ndarray, list[str], dict[str, slice | int], list[tuple[float, float]]]:
    nbasis = _spline_basis(np.asarray([(config.time_low_ns + config.time_high_ns) / 2]), config).shape[1]
    names, layout = _parameter_layout(len(stratum.runs), nbasis)
    values = np.zeros(len(names), dtype=float)
    values[layout["center"]] = 150.0
    values[layout["log_sigma"]] = math.log(0.75)
    values[layout["rho"]] = np.arctanh(0.25)
    values[layout["log_diag_sigma"]] = math.log(0.75)
    values[layout["diag_rho"]] = np.arctanh(0.65)
    values[layout["log_tau_offset"]] = math.log(0.30)
    values[layout["log_tau_width"]] = math.log(0.10)
    bounds: list[tuple[float, float]] = [
        (148.0, 152.0),
        (math.log(config.central_sigma_min_ns), math.log(config.central_sigma_max_ns)),
        (-3.0, 3.0),
        (math.log(config.diagonal_sigma_min_ns), math.log(config.diagonal_sigma_max_ns)),
        (-3.0, 3.0),
    ]
    bounds.extend([(-5.0, 5.0)] * (nbasis - 1))
    bounds.extend(((math.log(0.02), math.log(2.0)),
                   (math.log(0.01), math.log(0.5))))
    bounds.extend([(config.run_offset_min_ns, config.run_offset_max_ns)] * max(0, len(stratum.runs) - 1))
    bounds.extend([(config.run_log_width_min, config.run_log_width_max)] * max(0, len(stratum.runs) - 1))
    _require(len(bounds) == len(values), "internal parameter-layout mismatch")
    return values, names, layout, bounds


def _deviance_terms(observed: np.ndarray, expected: np.ndarray) -> np.ndarray:
    expected = np.maximum(expected, 1.0e-300)
    terms = expected - observed
    positive = observed > 0
    terms[positive] += observed[positive] * np.log(observed[positive] / expected[positive])
    return 2.0 * terms


def _conditional_yield_covariance(
    counts: np.ndarray, probabilities: np.ndarray, yields: np.ndarray,
    config: TimingFitConfig,
) -> tuple[np.ndarray, int, float]:
    occupied = counts > 0
    if not np.any(occupied):
        return np.zeros((len(COMPONENTS), len(COMPONENTS))), 0, math.inf
    mu = np.maximum(yields @ probabilities[:, occupied], config.minimum_expected)
    weighted = probabilities[:, occupied] * (np.sqrt(counts[occupied]) / mu)[None, :]
    information = weighted @ weighted.T
    singular = np.linalg.svd(information, compute_uv=False)
    tolerance = max(information.shape) * np.finfo(float).eps * (singular[0] if len(singular) else 0.0)
    rank = int(np.count_nonzero(singular > tolerance))
    covariance = np.linalg.pinv(information, rcond=1.0e-12)
    condition = float(singular[0] / singular[-1]) if len(singular) and singular[-1] > 0 else math.inf
    return covariance, rank, condition


def fit_stratum(stratum: StratumData, config: TimingFitConfig) -> StratumFit:
    start, names, layout, bounds = _initial_parameters(stratum, config)
    evaluation_failures: list[str] = []

    def evaluate(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        try:
            probabilities, _, _, _, _ = _component_probabilities(
                parameters, stratum, config, layout
            )
            run_nll = np.zeros(len(stratum.runs), dtype=float)
            for run_index in range(len(stratum.runs)):
                _, run_nll[run_index], _ = _profile_one_run_yields(
                    stratum.counts_by_run_cell[run_index], probabilities[run_index], config
                )
            value = float(run_nll.sum()) + _shape_penalty(
                parameters, layout, len(stratum.runs), config
            )
            if not np.isfinite(value):
                return 1.0e100, run_nll
            return float(value), run_nll
        except (FitError, FloatingPointError, ValueError) as error:
            if len(evaluation_failures) < 5:
                evaluation_failures.append(str(error))
            return 1.0e100, np.full(len(stratum.runs), 1.0e100)

    cache_x: np.ndarray | None = None
    cache_value = math.nan
    cache_run_nll: np.ndarray | None = None

    def cached_evaluate(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        nonlocal cache_x, cache_value, cache_run_nll
        if cache_x is None or not np.array_equal(cache_x, parameters):
            cache_value, cache_run_nll = evaluate(parameters)
            cache_x = np.array(parameters, copy=True)
        assert cache_run_nll is not None
        return cache_value, cache_run_nll

    def objective(parameters: np.ndarray) -> float:
        return cached_evaluate(parameters)[0]

    offset_slice = layout["offsets"]
    width_slice = layout["widths"]
    assert isinstance(offset_slice, slice) and isinstance(width_slice, slice)
    local_start = offset_slice.start
    tau_indices = {int(layout["log_tau_offset"]), int(layout["log_tau_width"])}

    def gradient(parameters: np.ndarray) -> np.ndarray:
        base_value, base_run_nll = cached_evaluate(parameters)
        derivative = np.zeros_like(parameters)
        nrun = len(stratum.runs)
        for index in range(len(parameters)):
            step = 2.0e-6 * max(1.0, abs(float(parameters[index])))
            direction = 1.0
            if parameters[index] + step > bounds[index][1]:
                step = -step
                direction = -1.0
            trial = np.array(parameters, copy=True)
            trial[index] += step
            if index < local_start and index not in tau_indices:
                trial_value, _ = evaluate(trial)
            elif index in tau_indices:
                base_penalty = _shape_penalty(parameters, layout, nrun, config)
                trial_penalty = _shape_penalty(trial, layout, nrun, config)
                trial_value = base_value + trial_penalty - base_penalty
            else:
                if index < width_slice.start:
                    free_index = index - offset_slice.start
                else:
                    free_index = index - width_slice.start
                affected = np.flatnonzero(
                    np.abs(_zero_sum_contrast(nrun)[:, free_index]) > 0.0
                ).tolist()
                trial_probabilities, _, _, _, _ = _component_probabilities(
                    trial, stratum, config, layout, run_indices=affected
                )
                trial_run_sum = 0.0
                base_affected_sum = 0.0
                for output_index, run_index in enumerate(affected):
                    _, value, _ = _profile_one_run_yields(
                        stratum.counts_by_run_cell[run_index],
                        trial_probabilities[output_index], config
                    )
                    trial_run_sum += value
                    base_affected_sum += base_run_nll[run_index]
                base_penalty = _shape_penalty(parameters, layout, nrun, config)
                trial_penalty = _shape_penalty(trial, layout, nrun, config)
                trial_value = (base_value - base_affected_sum + trial_run_sum -
                               base_penalty + trial_penalty)
            derivative[index] = (trial_value - base_value) / step
            if direction < 0:
                # ``step`` already carries the negative direction; this branch
                # documents the one-sided boundary derivative explicitly.
                derivative[index] = (trial_value - base_value) / step
        return derivative

    result = minimize(
        objective,
        start,
        jac=gradient,
        method="L-BFGS-B",
        bounds=bounds,
        options={
            "maxiter": config.optimizer_maxiter,
            "ftol": config.optimizer_tolerance,
            "gtol": config.optimizer_tolerance,
            "maxls": 50,
        },
    )
    _require(result.success,
             "timing-shape optimizer failed: " + str(result.message) +
             f"; iterations={result.nit}; objective={result.fun:.12g}; "
             f"max_abs_gradient={float(np.max(np.abs(result.jac))):.6g}" +
             ("; evaluations: " + "; ".join(evaluation_failures) if evaluation_failures else ""))
    probabilities, offsets, widths, tau_offset, tau_width = _component_probabilities(
        result.x, stratum, config, layout
    )
    integrated_yields, data_nll = _profile_all_yields(
        stratum.counts_by_run_cell, probabilities, config
    )
    covariance = np.asarray(result.hess_inv.todense(), dtype=float)
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues = np.linalg.eigvalsh(covariance)
    positive_eigenvalues = eigenvalues[eigenvalues > 0]
    covariance_condition = (
        float(positive_eigenvalues[-1] / positive_eigenvalues[0])
        if len(positive_eigenvalues) == len(eigenvalues) and len(eigenvalues) else math.inf
    )

    mass_count = len(config.mass_edges) + 1
    component_rows: list[dict[str, object]] = []
    for run_index, run in enumerate(stratum.runs):
        event_mask = stratum.run_index == run_index
        for mass_index in range(mass_count):
            group = event_mask & (stratum.mass_index == mass_index)
            counts = np.bincount(
                stratum.cell_index[group], minlength=stratum.counts_by_run_cell.shape[1]
            ).astype(float)
            yields, _, iterations = _profile_one_run_yields(
                counts, probabilities[run_index], config
            )
            conditional_cov, rank, condition = _conditional_yield_covariance(
                counts, probabilities[run_index], yields, config
            )
            if mass_index == 0:
                mass_label, mass_low, mass_high = "underflow", -math.inf, config.mass_low_gev
            elif mass_index == mass_count - 1:
                mass_label, mass_low, mass_high = "overflow", config.mass_high_gev, math.inf
            else:
                mass_label = "regular"
                mass_low = config.mass_edges[mass_index - 1]
                mass_high = config.mass_edges[mass_index]
            row: dict[str, object] = {
                "run_number": int(run),
                "acquisition_mode": stratum.mode,
                "multiplicity_class": stratum.multiplicity,
                "mass_index": mass_index,
                "mass_label": mass_label,
                "mass_low_gev": mass_low,
                "mass_high_gev": mass_high,
                "observed_count": int(np.sum(counts)),
                "em_iterations": iterations,
                "conditional_information_rank": rank,
                "conditional_information_condition": condition,
            }
            for k, component in enumerate(COMPONENTS):
                row[f"yield_{component}"] = float(yields[k])
                row[f"conditional_var_{component}"] = float(conditional_cov[k, k])
            for k, first in enumerate(COMPONENTS):
                for l, second in enumerate(COMPONENTS[k + 1:], start=k + 1):
                    row[f"conditional_cov_{first}_{second}"] = float(conditional_cov[k, l])
            component_rows.append(row)

    prediction_rows: list[dict[str, object]] = []
    ntime = len(config.time_edges) - 1
    for run_index, run in enumerate(stratum.runs):
        expected = integrated_yields[run_index] @ probabilities[run_index]
        observed = stratum.counts_by_run_cell[run_index].astype(float)
        deviance = _deviance_terms(observed, expected)
        pearson = (observed - expected) / np.sqrt(np.maximum(expected, config.minimum_expected))
        areas = accepted_cell_areas(config, stratum.pair_cut_by_run[run_index])
        for cell_index in range(len(expected)):
            i, j = divmod(cell_index, ntime)
            prediction_rows.append({
                "run_number": int(run),
                "acquisition_mode": stratum.mode,
                "multiplicity_class": stratum.multiplicity,
                "t1_low_ns": config.time_edges[i],
                "t1_high_ns": config.time_edges[i + 1],
                "t2_low_ns": config.time_edges[j],
                "t2_high_ns": config.time_edges[j + 1],
                "accepted_area_ns2": float(areas[cell_index]),
                "observed": int(observed[cell_index]),
                "expected": float(expected[cell_index]),
                "pearson_residual": float(pearson[cell_index]),
                "deviance_contribution": float(deviance[cell_index]),
            })

    run_parameter_rows: list[dict[str, object]] = []
    center = float(result.x[layout["center"]])
    sigma = math.exp(float(result.x[layout["log_sigma"]]))
    for r, run in enumerate(stratum.runs):
        run_parameter_rows.append({
            "run_number": int(run),
            "acquisition_mode": stratum.mode,
            "multiplicity_class": stratum.multiplicity,
            "event_count": int(stratum.event_count_by_run[r]),
            "pair_time_diff_max_ns": float(stratum.pair_cut_by_run[r]),
            "timing_offset_ns": float(offsets[r]),
            "timing_center_ns": center + float(offsets[r]),
            "log_width_deviation": float(widths[r]),
            "timing_sigma_ns": sigma * math.exp(float(widths[r])),
            "run_offset_population_scale_ns": tau_offset,
            "run_log_width_population_scale": tau_width,
            **{f"integrated_yield_{component}": float(integrated_yields[r, k])
               for k, component in enumerate(COMPONENTS)},
        })

    return StratumFit(
        mode=stratum.mode,
        multiplicity=stratum.multiplicity,
        runs=stratum.runs,
        pair_cuts=stratum.pair_cut_by_run,
        parameters=np.asarray(result.x),
        parameter_names=names,
        probabilities=probabilities,
        integrated_yields=integrated_yields,
        nll=float(data_nll),
        optimizer_success=bool(result.success),
        optimizer_message=str(result.message),
        optimizer_iterations=int(result.nit),
        profile_covariance=covariance,
        covariance_eigenvalues=eigenvalues,
        covariance_condition=covariance_condition,
        component_rows=component_rows,
        prediction_rows=prediction_rows,
        run_parameter_rows=run_parameter_rows,
    )


def _write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    _require(bool(rows), f"refusing to write empty table: {path.name}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _git_head(repo: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def fit_observation_bundle(
    bundle: ObservationBundle,
    output_dir: Path | str,
    config: TimingFitConfig | None = None,
    command: Sequence[str] | None = None,
) -> list[StratumFit]:
    """Fit all strata and write a self-contained shadow result bundle."""
    config = config or TimingFitConfig()
    output = Path(output_dir).resolve()
    repo = Path(__file__).resolve().parents[2]
    canonical = (repo / "output").resolve()
    _require(output != canonical and canonical not in output.parents,
             "ALG-002A refuses canonical production output paths")
    _require(not output.exists() or not any(output.iterdir()),
             f"output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)

    strata = build_strata(bundle, config)
    fits = [fit_stratum(stratum, config) for stratum in strata]
    _require(sum(sum(int(row["observed_count"]) for row in fit.component_rows)
                 for fit in fits) == len(bundle.arrays["run_number"]),
             "mass-spectrum outputs do not close to raw observation count")

    _write_csv(output / "input_manifest.csv", bundle.manifest)
    all_components = [row for fit in fits for row in fit.component_rows]
    _write_csv(output / "component_yields_by_run_mass.csv", all_components)
    true_rows = [{key: value for key, value in row.items()
                  if key in {
                      "run_number", "acquisition_mode", "multiplicity_class",
                      "mass_index", "mass_label", "mass_low_gev", "mass_high_gev",
                      "observed_count", "yield_true", "conditional_var_true",
                      "conditional_information_rank", "conditional_information_condition",
                  }} for row in all_components]
    _write_csv(output / "true_coincidence_mass_spectrum.csv", true_rows)
    _write_csv(output / "timing_predictions.csv",
               [row for fit in fits for row in fit.prediction_rows])
    _write_csv(output / "run_timing_parameters.csv",
               [row for fit in fits for row in fit.run_parameter_rows])

    summaries: list[dict[str, object]] = []
    for fit in fits:
        tag = f"mode{fit.mode}_mult{fit.multiplicity}"
        np.savez_compressed(
            output / f"{tag}_profile_covariance.npz",
            parameter_names=np.asarray(fit.parameter_names),
            parameters=fit.parameters,
            covariance=fit.profile_covariance,
            covariance_eigenvalues=fit.covariance_eigenvalues,
            component_probabilities=fit.probabilities,
            integrated_yields=fit.integrated_yields,
            run_numbers=fit.runs,
        )
        summary = {
            "tag": tag,
            "mode": fit.mode,
            "multiplicity_class": fit.multiplicity,
            "runs": [int(value) for value in fit.runs],
            "event_count": int(sum(row["event_count"] for row in fit.run_parameter_rows)),
            "optimizer_success": fit.optimizer_success,
            "optimizer_message": fit.optimizer_message,
            "optimizer_iterations": fit.optimizer_iterations,
            "data_negative_log_likelihood_without_constants": fit.nll,
            "profile_covariance_condition": fit.covariance_condition,
            "profile_covariance_status": "LBFGS_PROFILE_APPROXIMATION_NOT_COVERAGE_CALIBRATED",
            "parameters": dict(zip(fit.parameter_names, map(float, fit.parameters))),
        }
        summaries.append(summary)
        with (output / f"{tag}_fit_summary.json").open("w") as stream:
            json.dump(summary, stream, indent=2, sort_keys=True, allow_nan=True)

    provenance = {
        "status": "ALG002A_SHADOW_NOT_PRODUCTION",
        "publication_ready": False,
        "pi0_weight_written": False,
        "efficiency_inputs_read": False,
        "cross_section_formed": False,
        "config": asdict(config),
        "command": list(command) if command is not None else None,
        "working_directory": os.getcwd(),
        "git_head": _git_head(repo),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "uproot": uproot.__version__,
        "input_entries": int(len(bundle.arrays["run_number"])),
        "input_run_numbers": sorted(map(int, np.unique(bundle.arrays["run_number"]))),
        "strata": summaries,
        "covariance_limit": (
            "Saved inverse-LBFGS profile curvature is diagnostic. Conditional mass-bin "
            "yield covariance is saved, but cross-bin shape-nuisance propagation and "
            "toy-calibrated profile coverage remain promotion gates."
        ),
    }
    with (output / "provenance.json").open("w") as stream:
        json.dump(provenance, stream, indent=2, sort_keys=True, allow_nan=True)
    (output / "FIT_STATUS.txt").write_text(
        "ALG002A_SHADOW_NOT_PRODUCTION\n"
        "No pi0_weight, efficiency, or cross section is produced.\n"
        "Inspect provenance.json and the thesis-grade method document before use.\n"
    )
    return fits
