#!/usr/bin/env python3
"""ALG-002B simultaneous timing/mass shadow model.

The implementation uses the site-provided NumPy/SciPy stack. A coordinate
profile optimizer minimizes one six-component extended-Poisson likelihood:
per-run component yields are profiled for every trial; mass parameters are
optimized at the current timing iterate; then all ALG-002A timing parameters
are refitted against the mass-resolved likelihood.

No event-level pi0 weight, efficiency, charge normalization, or cross section
is produced. Canonical output remains refused unless the caller explicitly
authorizes the isolated `output/KinC_x36_4/alg002b/` shadow subtree.
"""

from __future__ import annotations

import csv
import json
import math
import os
import platform
import sys
from dataclasses import asdict, dataclass, field
from multiprocessing import get_all_start_methods, get_context
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.special import betainc, erf

try:
    from .joint_timing_model import (
        FitError, ObservationBundle, StratumData, TimingFitConfig,
        _component_probabilities, _git_head, _initial_parameters,
        _multiplicity_class, _shape_penalty, _write_csv,
        _zero_sum_deviations, build_strata, fit_stratum,
    )
except ImportError:  # Direct execution from src/background_fit.
    from joint_timing_model import (  # type: ignore
        FitError, ObservationBundle, StratumData, TimingFitConfig,
        _component_probabilities, _git_head, _initial_parameters,
        _multiplicity_class, _shape_penalty, _write_csv,
        _zero_sum_deviations, build_strata, fit_stratum,
    )


COMPONENTS = ("pi0", "combinatorial", "horizontal", "vertical", "random", "diagonal")
BACKGROUND_GROUPS = ("combinatorial", "horizontal_vertical", "random", "diagonal")
_PARALLEL_OBJECTIVE: Callable[[np.ndarray], float] | None = None


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FitError(message)


@dataclass(frozen=True)
class JointMassFitConfig:
    timing: TimingFitConfig = field(default_factory=TimingFitConfig)
    signal_model: str = "dscb"
    combinatorial_model: str = "bernstein3"
    starts: int = 3
    start_index_offset: int = 0
    seed: int = 20261001
    coordinate_cycles: int = 2
    nproc: int = 1
    mass_maxiter: int = 120
    timing_refit_maxiter: int = 25
    optimizer_tolerance: float = 2.0e-7
    yield_em_maxiter: int = 500
    yield_em_tolerance: float = 1.0e-9
    background_smoothness: float = 0.5
    closure_tolerance: float = 2.0e-3

    def __post_init__(self) -> None:
        if self.signal_model not in {"dscb", "double_gaussian"}:
            raise ValueError(f"unsupported signal model: {self.signal_model}")
        if self.combinatorial_model not in {"logistic", "bernstein3", "bernstein4"}:
            raise ValueError(f"unsupported combinatorial model: {self.combinatorial_model}")
        if self.starts < 1 or self.coordinate_cycles < 1:
            raise ValueError("starts and coordinate_cycles must be positive")
        if self.start_index_offset < 0:
            raise ValueError("start_index_offset must be nonnegative")
        if self.nproc < 1:
            raise ValueError("nproc must be positive")
        if self.mass_maxiter < 1 or self.timing_refit_maxiter < 1:
            raise ValueError("mass and timing iteration limits must be positive")
        if self.yield_em_maxiter < 1 or self.yield_em_tolerance <= 0.0:
            raise ValueError("yield EM controls must be positive")


@dataclass
class JointStratum:
    data: StratumData
    tag: str
    global_run_index: np.ndarray
    event_run_index: np.ndarray
    event_mass_index: np.ndarray
    event_cell_index: np.ndarray
    event_indices_by_run: list[np.ndarray]


@dataclass
class JointDataset:
    bundle: ObservationBundle
    runs: np.ndarray
    strata: list[JointStratum]


@dataclass
class TimingState:
    parameters: np.ndarray
    names: list[str]
    layout: dict[str, slice | int]
    bounds: list[tuple[float, float]]


@dataclass
class MassLayout:
    names: list[str]
    bounds: list[tuple[float, float]]
    index: dict[str, slice | int]


@dataclass
class Evaluation:
    objective: float
    yields: list[np.ndarray]
    timing_probabilities: list[np.ndarray]
    mass_probabilities: list[np.ndarray]
    run_nll: list[np.ndarray]
    yield_iterations: list[np.ndarray]
    yield_kkt_residuals: list[np.ndarray]
    penalty: float


@dataclass
class JointFitResult:
    dataset: JointDataset
    config: JointMassFitConfig
    mass_parameters: np.ndarray
    mass_layout: MassLayout
    timing_states: list[TimingState]
    evaluation: Evaluation
    start_summaries: list[dict[str, object]]
    cycle_summaries: list[dict[str, object]]
    run_yield_rows: list[dict[str, object]]
    parameter_rows: list[dict[str, object]]
    prediction_rows: list[dict[str, object]]
    covariance: np.ndarray
    covariance_labels: list[str]
    profile_rows: list[dict[str, object]]
    closure: dict[str, float]


def _normalize_row(row: dict[str, object]) -> dict[str, str]:
    return {str(key).strip(): str(value or "").strip() for key, value in row.items()}


def expected_lh2_runs(
    config_csv: Path | str, kin: str,
) -> tuple[list[int], dict[int, dict[str, str]]]:
    """Return the exact configured run set after strict LH2 filtering."""
    path = Path(config_csv).resolve()
    _require(path.is_file(), f"missing run configuration: {path}")
    selected: dict[int, dict[str, str]] = {}
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream, skipinitialspace=True):
            row = _normalize_row(raw)
            if row.get("Kin_old") != kin:
                continue
            if row.get("Type", "").casefold() != "production":
                continue
            if row.get("target", "").upper() != "LH2":
                continue
            text = row.get("run_number", "")
            _require(text.isdigit(), f"invalid run number in {path}: {text!r}")
            run = int(text)
            _require(run not in selected, f"duplicate configured production-LH2 run {run}")
            selected[run] = row
    _require(bool(selected), f"no production-LH2 runs found for {kin} in {path}")
    return sorted(selected), selected


def enforce_lh2_manifest(
    bundle: ObservationBundle, expected_runs: Sequence[int],
    allowed_missing: Sequence[int] = (6569,),
) -> dict[str, object]:
    represented = sorted(map(int, np.unique(bundle.arrays["run_number"])))
    expected = sorted({int(value) for value in expected_runs})
    unexpected = sorted(set(represented) - set(expected))
    missing = sorted(set(expected) - set(represented))
    allowed = sorted({int(value) for value in allowed_missing})
    _require(not unexpected,
             "input contains runs outside configured production-LH2 set: " + str(unexpected))
    _require(not (set(missing) - set(allowed)),
             "configured production-LH2 runs missing without approved exclusion: " +
             str(sorted(set(missing) - set(allowed))))
    manifest_runs = sorted(
        int(row["run_number"]) for row in bundle.manifest
        if row.get("run_number") is not None
    )
    _require(manifest_runs == represented,
             "input manifest does not map one file to every represented LH2 run")
    return {
        "target": "LH2", "expected_runs": expected,
        "represented_runs": represented, "missing_runs": missing,
        "allowed_missing_runs": allowed, "unexpected_runs": unexpected,
    }


def build_joint_dataset(
    bundle: ObservationBundle, timing_config: TimingFitConfig | None = None,
) -> JointDataset:
    config = timing_config or TimingFitConfig()
    base = build_strata(bundle, config)
    runs = np.unique(bundle.arrays["run_number"].astype(np.int64))
    global_map = {int(run): index for index, run in enumerate(runs)}
    mult = _multiplicity_class(bundle.arrays["nclust_selected"].astype(int))
    strata: list[JointStratum] = []
    for data in base:
        selected = ((bundle.arrays["acquisition_mode"] == data.mode) &
                    (mult == data.multiplicity))
        _require(np.count_nonzero(selected) == len(data.run_index),
                 "stratum event ordering mismatch")
        strata.append(JointStratum(
            data=data, tag=f"mode{data.mode}_mult{data.multiplicity}",
            global_run_index=np.asarray([global_map[int(run)] for run in data.runs], dtype=np.int64),
            event_run_index=data.run_index.astype(np.int64),
            event_mass_index=data.mass_index.astype(np.int64),
            event_cell_index=data.cell_index.astype(np.int64),
            event_indices_by_run=[np.flatnonzero(data.run_index == index)
                                  for index in range(len(data.runs))],
        ))
    _require(sum(len(item.event_run_index) for item in strata) == len(bundle.arrays["run_number"]),
             "joint strata do not close to raw-observation count")
    return JointDataset(bundle=bundle, runs=runs, strata=strata)


def _bernstein_integrals(edges: np.ndarray, degree: int) -> np.ndarray:
    x = (edges - edges[0]) / (edges[-1] - edges[0])
    basis = np.asarray([np.diff(betainc(k + 1, degree - k + 1, x))
                        for k in range(degree + 1)])
    basis = np.maximum(basis, 0.0)
    basis /= basis.sum(axis=1, keepdims=True)
    return basis


def _softmax(values: np.ndarray) -> np.ndarray:
    shifted = values - np.max(values)
    result = np.exp(np.clip(shifted, -700.0, 0.0))
    return result / result.sum()


def _dscb_cdf(
    x: np.ndarray, mean: np.ndarray, sigma: np.ndarray,
    alpha_left: float, n_left: float, alpha_right: float, n_right: float,
) -> np.ndarray:
    t = (x - mean) / sigma
    al, nl, ar, nr = alpha_left, n_left, alpha_right, n_right
    acl = (nl / al) ** nl * math.exp(-0.5 * al * al); bl = nl / al - al
    acr = (nr / ar) ** nr * math.exp(-0.5 * ar * ar); br = nr / ar - ar
    left_total = acl / (nl - 1.0) * (bl + al) ** (1.0 - nl)
    core_total = math.sqrt(math.pi / 2.0) * (
        erf(ar / math.sqrt(2.0)) + erf(al / math.sqrt(2.0)))
    right_total = acr / (nr - 1.0) * (br + ar) ** (1.0 - nr)
    total = left_total + core_total + right_total
    # Evaluate only the active branch. np.where evaluates both inputs eagerly,
    # which overflows unused clipped powers for far-tail finite differences.
    result = np.empty_like(t, dtype=float)
    left_mask = t < -al
    core_mask = (t >= -al) & (t <= ar)
    right_mask = t > ar
    result[left_mask] = (
        acl / (nl - 1.0) * (bl - t[left_mask]) ** (1.0 - nl))
    result[core_mask] = left_total + math.sqrt(math.pi / 2.0) * (
        erf(t[core_mask] / math.sqrt(2.0)) + erf(al / math.sqrt(2.0)))
    result[right_mask] = left_total + core_total + acr / (nr - 1.0) * (
        (br + ar) ** (1.0 - nr) -
        (br + t[right_mask]) ** (1.0 - nr))
    return np.clip(result / total, 0.0, 1.0)


def _mass_layout(dataset: JointDataset, config: JointMassFitConfig) -> tuple[np.ndarray, MassLayout]:
    names: list[str] = []; bounds: list[tuple[float, float]] = []
    values: list[float] = []; index: dict[str, slice | int] = {}
    def scalar(name: str, value: float, bound: tuple[float, float]) -> None:
        index[name] = len(values); names.append(name); values.append(value); bounds.append(bound)
    def vector(name: str, value: Sequence[float], bound: tuple[float, float]) -> None:
        start = len(values); values.extend(value)
        names.extend(f"{name}_{i}" for i in range(len(value)))
        bounds.extend([bound] * len(value)); index[name] = slice(start, len(values))
    scalar("log_tau_mass_mean", math.log(0.00045), (math.log(2.0e-5), math.log(0.005)))
    scalar("log_tau_mass_width", math.log(0.08), (math.log(0.003), math.log(0.5)))
    free = max(0, len(dataset.runs) - 1)
    vector("mass_mean_free", np.zeros(free), (-0.005, 0.005))
    vector("mass_width_free", np.zeros(free), (-0.7, 0.7))
    if config.signal_model == "dscb":
        scalar("alpha_left", 1.5, (0.5, 5.0)); scalar("n_left", 3.0, (1.05, 30.0))
        scalar("alpha_right", 2.0, (0.5, 6.0)); scalar("n_right", 4.0, (1.05, 30.0))
    else:
        scalar("wide_fraction", 0.15, (0.0, 0.6)); scalar("wide_scale", 2.0, (1.05, 5.0))
    for joint in dataset.strata:
        tag = joint.tag
        scalar(f"{tag}:signal_mean", 0.1348, (0.125, 0.145))
        scalar(f"{tag}:log_signal_sigma", math.log(0.0045), (math.log(0.0015), math.log(0.020)))
        if config.combinatorial_model == "logistic":
            scalar(f"{tag}:comb_turn", 0.18, (0.08, 0.35))
            scalar(f"{tag}:log_comb_width", math.log(0.06), (math.log(0.003), math.log(0.20)))
        else:
            degree = 4 if config.combinatorial_model == "bernstein4" else 3
            vector(f"{tag}:combinatorial_logits", np.zeros(degree), (-10.0, 10.0))
        for group in ("horizontal_vertical", "random", "diagonal"):
            vector(f"{tag}:{group}_logits", np.zeros(3), (-10.0, 10.0))
        for group in BACKGROUND_GROUPS:
            # The overflow logit is the zero reference.  Removing the common
            # softmax shift makes every background-shape parameter identifiable.
            vector(f"{tag}:{group}_flow_logits",
                   np.log(np.asarray((1.0e-5, 0.95)) / 0.05), (-16.0, 16.0))
    return np.asarray(values), MassLayout(names, bounds, index)


def _mass_run_deviations(
    parameters: np.ndarray, layout: MassLayout, nrun: int,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    mean = _zero_sum_deviations(parameters[layout.index["mass_mean_free"]], nrun)
    width = _zero_sum_deviations(parameters[layout.index["mass_width_free"]], nrun)
    return (mean, width,
            math.exp(float(parameters[layout.index["log_tau_mass_mean"]])),
            math.exp(float(parameters[layout.index["log_tau_mass_width"]])))


def _background_regular(
    parameters: np.ndarray, layout: MassLayout, config: JointMassFitConfig,
    tag: str, group: str,
) -> np.ndarray:
    edges = config.timing.mass_edges
    if group == "combinatorial" and config.combinatorial_model == "logistic":
        turn = float(parameters[layout.index[f"{tag}:comb_turn"]])
        width = math.exp(float(parameters[layout.index[f"{tag}:log_comb_width"]]))
        z = np.clip((edges - turn) / width, -700.0, 700.0)
        anti = edges - width * np.logaddexp(0.0, z)
        result = np.maximum(np.diff(anti), 1.0e-15)
        return result / result.sum()
    degree = 4 if group == "combinatorial" and config.combinatorial_model == "bernstein4" else 3
    logits = np.concatenate((
        parameters[layout.index[f"{tag}:{group}_logits"]], [0.0]
    ))
    result = _softmax(logits) @ _bernstein_integrals(edges, degree)
    return result / result.sum()


def _mass_probabilities(
    parameters: np.ndarray, layout: MassLayout,
    dataset: JointDataset, config: JointMassFitConfig,
) -> tuple[list[np.ndarray], list[dict[str, np.ndarray | float]], float]:
    mean_delta, width_delta, tau_mean, tau_width = _mass_run_deviations(
        parameters, layout, len(dataset.runs))
    penalty = 0.0
    if len(dataset.runs) > 1:
        penalty += 0.5 * float(mean_delta @ mean_delta) / tau_mean**2
        penalty += (len(dataset.runs) - 1) * math.log(tau_mean)
        penalty += 0.5 * float(width_delta @ width_delta) / tau_width**2
        penalty += (len(dataset.runs) - 1) * math.log(tau_width)
    outputs: list[np.ndarray] = []; details: list[dict[str, np.ndarray | float]] = []
    edges = config.timing.mass_edges
    for joint in dataset.strata:
        tag = joint.tag
        base_mean = float(parameters[layout.index[f"{tag}:signal_mean"]])
        base_sigma = math.exp(float(parameters[layout.index[f"{tag}:log_signal_sigma"]]))
        means = base_mean + mean_delta[joint.global_run_index]
        sigmas = base_sigma * np.exp(width_delta[joint.global_run_index])
        x = edges[None, :]
        if config.signal_model == "dscb":
            cdf = _dscb_cdf(
                x, means[:, None], sigmas[:, None],
                float(parameters[layout.index["alpha_left"]]),
                float(parameters[layout.index["n_left"]]),
                float(parameters[layout.index["alpha_right"]]),
                float(parameters[layout.index["n_right"]]))
        else:
            fraction = float(parameters[layout.index["wide_fraction"]])
            scale = float(parameters[layout.index["wide_scale"]])
            z1 = (x - means[:, None]) / (sigmas[:, None] * math.sqrt(2.0))
            z2 = (x - means[:, None]) / (sigmas[:, None] * scale * math.sqrt(2.0))
            cdf = ((1.0 - fraction) * 0.5 * (1.0 + erf(z1)) +
                   fraction * 0.5 * (1.0 + erf(z2)))
        cdf = np.clip(cdf, 0.0, 1.0)
        signal = np.concatenate((cdf[:, :1], np.maximum(np.diff(cdf), 0.0),
                                 1.0 - cdf[:, -1:]), axis=1)
        signal /= signal.sum(axis=1, keepdims=True)
        backgrounds = []
        for group in BACKGROUND_GROUPS:
            regular = _background_regular(parameters, layout, config, tag, group)
            flow = _softmax(np.concatenate((
                parameters[layout.index[f"{tag}:{group}_flow_logits"]], [0.0]
            )))
            backgrounds.append(np.concatenate((flow[:1], flow[1] * regular, flow[2:])))
            if not (group == "combinatorial" and config.combinatorial_model == "logistic"):
                logits = np.concatenate((
                    parameters[layout.index[f"{tag}:{group}_logits"]], [0.0]
                ))
                penalty += 0.5 * config.background_smoothness * float(np.diff(logits, n=2) @ np.diff(logits, n=2))
        mass = np.zeros((len(joint.data.runs), 6, signal.shape[1]))
        mass[:, 0, :] = signal; mass[:, 1, :] = backgrounds[0]
        mass[:, 2, :] = backgrounds[1]; mass[:, 3, :] = backgrounds[1]
        mass[:, 4, :] = backgrounds[2]; mass[:, 5, :] = backgrounds[3]
        outputs.append(mass)
        details.append({"base_mean": base_mean, "base_sigma": base_sigma,
                        "means": means, "sigmas": sigmas})
    return outputs, details, penalty


def _profile_yields(
    features: np.ndarray, config: JointMassFitConfig,
) -> tuple[np.ndarray, float, int, float]:
    count = len(features)
    if count == 0:
        return np.zeros(6), 0.0, 0, 0.0
    yields = count * np.asarray((0.48, 0.14, 0.10, 0.10, 0.13, 0.05))
    converged = False
    for iteration in range(config.yield_em_maxiter):
        terms = yields[None, :] * features
        denominator = np.maximum(terms.sum(axis=1), 1.0e-300)
        updated = np.sum(terms / denominator[:, None], axis=0)
        delta = float(np.max(np.abs(updated - yields) / np.maximum(1.0, yields)))
        yields = np.maximum(updated, 0.0)
        if delta < config.yield_em_tolerance:
            converged = True
            break
    def yield_objective(candidate: np.ndarray) -> float:
        intensity = np.maximum(features @ candidate, 1.0e-300)
        return float(candidate.sum() - np.log(intensity).sum())

    def yield_gradient(candidate: np.ndarray) -> np.ndarray:
        intensity = np.maximum(features @ candidate, 1.0e-300)
        return 1.0 - features.T @ (1.0 / intensity)

    def kkt_residual(candidate: np.ndarray) -> float:
        gradient = yield_gradient(candidate)
        active = candidate > 1.0e-10 * max(1.0, count)
        residuals = []
        if np.any(active):
            residuals.append(float(np.max(np.abs(gradient[active]))))
        if np.any(~active):
            residuals.append(float(np.max(np.maximum(-gradient[~active], 0.0))))
        return max(residuals, default=0.0)

    residual = kkt_residual(yields)
    optimizer_iterations = 0
    if not converged or residual > 1.0e-6:
        result = minimize(
            yield_objective, yields, jac=yield_gradient, method="L-BFGS-B",
            bounds=[(0.0, None)] * 6,
            options={"maxiter": 1000, "ftol": 1.0e-12, "gtol": 1.0e-9},
        )
        candidate = np.maximum(np.asarray(result.x), 0.0)
        candidate_residual = kkt_residual(candidate)
        if (yield_objective(candidate) <= yield_objective(yields) or
                candidate_residual < residual):
            yields, residual = candidate, candidate_residual
        optimizer_iterations += int(result.nit)
        if residual > 2.0e-6:
            second = minimize(
                yield_objective, yields, jac=yield_gradient, method="SLSQP",
                bounds=[(0.0, None)] * 6,
                options={"maxiter": 1000, "ftol": 1.0e-12, "disp": False},
            )
            candidate = np.maximum(np.asarray(second.x), 0.0)
            candidate_residual = kkt_residual(candidate)
            if (yield_objective(candidate) <= yield_objective(yields) or
                    candidate_residual < residual):
                yields, residual = candidate, candidate_residual
            optimizer_iterations += int(second.nit)
            _require(
                residual <= 2.0e-6,
                "six-component yield profile failed KKT check: "
                f"lbfgs={result.message}; slsqp={second.message}; "
                f"residual={residual:.6g}",
            )
    intensity = np.maximum(features @ yields, 1.0e-300)
    return (yields, float(yields.sum() - np.log(intensity).sum()),
            iteration + 1 + optimizer_iterations, residual)


def _evaluate(
    mass_parameters: np.ndarray, mass_layout: MassLayout,
    timing_states: Sequence[TimingState], dataset: JointDataset,
    config: JointMassFitConfig, only_stratum: int | None = None,
) -> Evaluation:
    mass_prob, _, mass_penalty = _mass_probabilities(mass_parameters, mass_layout, dataset, config)
    total_nll = 0.0; total_penalty = mass_penalty if only_stratum is None else 0.0
    all_yields: list[np.ndarray] = []; all_timing: list[np.ndarray] = []
    all_run_nll: list[np.ndarray] = []; all_iterations: list[np.ndarray] = []
    all_kkt: list[np.ndarray] = []
    for s, joint in enumerate(dataset.strata):
        timing5, _, _, _, _ = _component_probabilities(
            timing_states[s].parameters, joint.data, config.timing, timing_states[s].layout)
        timing6 = np.stack((timing5[:, 0], timing5[:, 0], timing5[:, 1],
                            timing5[:, 2], timing5[:, 3], timing5[:, 4]), axis=1)
        yields = np.zeros((len(joint.data.runs), 6)); run_nll = np.zeros(len(joint.data.runs))
        iterations = np.zeros(len(joint.data.runs), dtype=np.int32)
        kkt = np.zeros(len(joint.data.runs))
        for r, event_indices in enumerate(joint.event_indices_by_run):
            features = (timing6[r, :, joint.event_cell_index[event_indices]] *
                        mass_prob[s][r, :, joint.event_mass_index[event_indices]])
            yields[r], run_nll[r], iterations[r], kkt[r] = _profile_yields(features, config)
        all_yields.append(yields); all_timing.append(timing6); all_run_nll.append(run_nll)
        all_iterations.append(iterations); all_kkt.append(kkt)
        if only_stratum is None or only_stratum == s:
            total_nll += float(run_nll.sum())
            total_penalty += _shape_penalty(
                timing_states[s].parameters, timing_states[s].layout,
                len(joint.data.runs), config.timing)
    return Evaluation(total_nll + total_penalty, all_yields, all_timing,
                      mass_prob, all_run_nll, all_iterations, all_kkt,
                      total_penalty)


def _load_timing_states(
    dataset: JointDataset, config: JointMassFitConfig,
    initial_dir: Path | str | None,
) -> list[TimingState]:
    output: list[TimingState] = []
    directory = Path(initial_dir).resolve() if initial_dir is not None else None
    for joint in dataset.strata:
        start, names, layout, bounds = _initial_parameters(joint.data, config.timing)
        if directory is not None:
            path = directory / f"{joint.tag}_profile_covariance.npz"
            _require(path.is_file(), f"missing ALG-002A timing initializer: {path}")
            with np.load(path, allow_pickle=False) as arrays:
                _require(arrays["parameter_names"].astype(str).tolist() == names,
                         f"timing parameter schema mismatch: {joint.tag}")
                _require(np.array_equal(arrays["run_numbers"].astype(np.int64), joint.data.runs),
                         f"timing run order mismatch: {joint.tag}")
                start = arrays["parameters"].astype(float)
        else:
            start = fit_stratum(joint.data, config.timing).parameters
        output.append(TimingState(np.asarray(start), names, layout, bounds))
    return output


def _scipy_supports_parallel_workers() -> bool:
    numbers = scipy.__version__.split(".")
    try:
        return (int(numbers[0]), int(numbers[1])) >= (1, 16)
    except (IndexError, ValueError):
        return False


def _evaluate_parallel_objective(values: np.ndarray) -> float:
    if _PARALLEL_OBJECTIVE is None:
        raise RuntimeError("parallel objective was not initialized")
    return _PARALLEL_OBJECTIVE(values)


def _minimize_lbfgsb(
    objective: Callable[[np.ndarray], float], start: np.ndarray,
    bounds: Sequence[tuple[float, float]], maxiter: int,
    config: JointMassFitConfig,
) -> object:
    options: dict[str, object] = {
        "maxiter": maxiter,
        "ftol": config.optimizer_tolerance,
        "gtol": config.optimizer_tolerance,
        "maxls": 40,
    }
    if config.nproc == 1:
        return minimize(
            objective, start, method="L-BFGS-B", bounds=bounds,
            options=options)
    _require(
        _scipy_supports_parallel_workers(),
        "nproc > 1 requires SciPy >= 1.16; use the managed analysis Python",
    )
    _require(
        "fork" in get_all_start_methods(),
        "nproc > 1 requires the POSIX fork multiprocessing start method",
    )
    # Forked workers share the read-only fit dataset through copy-on-write and
    # bypass Python's global interpreter lock during numerical differentiation.
    global _PARALLEL_OBJECTIVE
    _PARALLEL_OBJECTIVE = objective
    try:
        with get_context("fork").Pool(processes=config.nproc) as pool:
            def parallel_map(
                _scipy_function: Callable[[np.ndarray], float],
                values: Sequence[np.ndarray],
            ) -> list[float]:
                return pool.map(_evaluate_parallel_objective, values)
            options["workers"] = parallel_map
            return minimize(
                objective, start, method="L-BFGS-B", bounds=bounds,
                options=options)
    finally:
        _PARALLEL_OBJECTIVE = None


def _fit_mass_block(
    start: np.ndarray, layout: MassLayout, timing: Sequence[TimingState],
    dataset: JointDataset, config: JointMassFitConfig,
) -> tuple[np.ndarray, object]:
    def objective(values: np.ndarray) -> float:
        try:
            value = _evaluate(values, layout, timing, dataset, config).objective
            return value if np.isfinite(value) else 1.0e100
        except (FitError, FloatingPointError, ValueError, OverflowError):
            return 1.0e100
    result = _minimize_lbfgsb(
        objective, start, layout.bounds, config.mass_maxiter, config)
    _require(np.isfinite(result.fun), "nonfinite mass-block optimizer result")
    return np.asarray(result.x), result


def _fit_timing_block(
    s: int, mass: np.ndarray, layout: MassLayout, timing: list[TimingState],
    dataset: JointDataset, config: JointMassFitConfig,
) -> object:
    state = timing[s]
    def objective(values: np.ndarray) -> float:
        trial_timing = list(timing)
        trial_timing[s] = TimingState(
            np.asarray(values), state.names, state.layout, state.bounds)
        try:
            value = _evaluate(
                mass, layout, trial_timing, dataset, config,
                only_stratum=s).objective
            return value if np.isfinite(value) else 1.0e100
        except (FitError, FloatingPointError, ValueError, OverflowError):
            return 1.0e100
    result = _minimize_lbfgsb(
        objective, state.parameters, state.bounds,
        config.timing_refit_maxiter, config)
    _require(np.isfinite(result.fun), f"nonfinite timing result for {dataset.strata[s].tag}")
    state.parameters = np.asarray(result.x)
    return result


def fit_joint_model(
    dataset: JointDataset, config: JointMassFitConfig,
    timing_initial_dir: Path | str | None = None,
) -> tuple[np.ndarray, MassLayout, list[TimingState], Evaluation,
           list[dict[str, object]], list[dict[str, object]]]:
    initial, layout = _mass_layout(dataset, config)
    base_timing = _load_timing_states(dataset, config, timing_initial_dir)
    best = None
    start_summaries: list[dict[str, object]] = []; all_cycles: list[dict[str, object]] = []
    for start_index in range(config.starts):
        absolute_start_index = config.start_index_offset + start_index
        mass = initial.copy()
        timing = [TimingState(x.parameters.copy(), x.names, x.layout, x.bounds) for x in base_timing]
        if absolute_start_index:
            rng = np.random.default_rng(
                np.random.SeedSequence([config.seed, absolute_start_index]))
            lower = np.asarray([x[0] for x in layout.bounds]); upper = np.asarray([x[1] for x in layout.bounds])
            mass = np.clip(mass + rng.normal(0.0, 0.05, len(mass)) * (upper - lower), lower, upper)
            for state in timing:
                lower = np.asarray([x[0] for x in state.bounds]); upper = np.asarray([x[1] for x in state.bounds])
                state.parameters = np.clip(
                    state.parameters + rng.normal(0.0, 0.03, len(state.parameters)) * (upper - lower),
                    lower, upper)
        previous = math.inf; cycles = []
        for cycle in range(config.coordinate_cycles):
            mass, mass_result = _fit_mass_block(mass, layout, timing, dataset, config)
            timing_records = []
            for s, joint in enumerate(dataset.strata):
                timing_result = _fit_timing_block(s, mass, layout, timing, dataset, config)
                timing_records.append({"tag": joint.tag, "success": bool(timing_result.success),
                                       "iterations": int(timing_result.nit),
                                       "message": str(timing_result.message)})
            evaluation = _evaluate(mass, layout, timing, dataset, config)
            cycles.append({"start_index": absolute_start_index, "cycle": cycle,
                           "objective": evaluation.objective,
                           "mass_success": bool(mass_result.success),
                           "mass_iterations": int(mass_result.nit),
                           "mass_message": str(mass_result.message), "timing": timing_records,
                           "relative_change": ((previous - evaluation.objective) /
                                               max(1.0, abs(previous))) if np.isfinite(previous) else None})
            previous = evaluation.objective
        final = _evaluate(mass, layout, timing, dataset, config)
        converged = all(
            bool(record["mass_success"]) and
            all(bool(item["success"]) for item in record["timing"])
            for record in cycles
        )
        start_summaries.append({"start_index": absolute_start_index, "objective": final.objective,
                                "pi0_yield_sum": float(sum(x[:, 0].sum() for x in final.yields)),
                                "finite": bool(np.isfinite(final.objective)),
                                "converged": converged, "cycles": len(cycles)})
        all_cycles.extend(cycles)
        if best is None or final.objective < best[2].objective:
            best = (mass.copy(), timing, final)
    _require(best is not None, "all joint optimizer starts failed")
    return best[0], layout, best[1], best[2], start_summaries, all_cycles


def _conditional_covariance(
    dataset: JointDataset, evaluation: Evaluation,
) -> tuple[np.ndarray, list[str], list[dict[str, object]]]:
    labels = []; variances = []; rows = []
    for s, joint in enumerate(dataset.strata):
        for r, event_indices in enumerate(joint.event_indices_by_run):
            features = (evaluation.timing_probabilities[s][r, :, joint.event_cell_index[event_indices]] *
                        evaluation.mass_probabilities[s][r, :, joint.event_mass_index[event_indices]])
            yields = evaluation.yields[s][r]; mu = np.maximum(features @ yields, 1.0e-300)
            information = (features / mu[:, None]).T @ (features / mu[:, None])
            covariance = np.linalg.pinv(information, rcond=1.0e-10)
            variance = max(float(covariance[0, 0]), 0.0); error = math.sqrt(variance)
            run = int(joint.data.runs[r]); labels.append(f"run{run}_{joint.tag}_pi0"); variances.append(variance)
            rows.append({"run_number": run, "acquisition_mode": joint.data.mode,
                         "multiplicity_class": joint.data.multiplicity,
                         "pi0_yield": float(yields[0]), "conditional_standard_error": error,
                         "conditional_68_low": max(0.0, float(yields[0]) - error),
                         "conditional_68_high": float(yields[0]) + error,
                         "conditional_95_low": max(0.0, float(yields[0]) - 1.96 * error),
                         "conditional_95_high": float(yields[0]) + 1.96 * error,
                         "status": "SHAPES_FIXED_CONDITIONAL_NOT_COVERAGE_CALIBRATED"})
    return np.diag(variances), labels, rows


def _materialize(
    dataset: JointDataset, config: JointMassFitConfig, mass: np.ndarray,
    layout: MassLayout, timing: list[TimingState], evaluation: Evaluation,
    starts: list[dict[str, object]], cycles: list[dict[str, object]],
) -> JointFitResult:
    mass_prob, mass_details, _ = _mass_probabilities(mass, layout, dataset, config)
    mean_delta, width_delta, tau_mean, tau_width = _mass_run_deviations(mass, layout, len(dataset.runs))
    parameter_rows = [
        {"scope": "optimizer", "tag": "mass", "run_number": "",
         "parameter": name, "value": float(value)}
        for name, value in zip(layout.names, mass)
    ]
    for i, run in enumerate(dataset.runs):
        parameter_rows.extend((
            {"scope": "run", "tag": "all", "run_number": int(run),
             "parameter": "pi0_mass_mean_deviation_gev", "value": mean_delta[i]},
            {"scope": "run", "tag": "all", "run_number": int(run),
             "parameter": "pi0_mass_log_width_deviation", "value": width_delta[i]}))
    parameter_rows.extend((
        {"scope": "global", "tag": "all", "run_number": "",
         "parameter": "pi0_mass_mean_population_scale_gev", "value": tau_mean},
        {"scope": "global", "tag": "all", "run_number": "",
         "parameter": "pi0_mass_log_width_population_scale", "value": tau_width}))
    run_rows = []; predictions = []
    for s, joint in enumerate(dataset.strata):
        state = timing[s]
        parameter_rows.extend(
            {"scope": "optimizer", "tag": joint.tag, "run_number": "",
             "parameter": name, "value": float(value)}
            for name, value in zip(state.names, state.parameters)
        )
        _, offsets, time_widths, tau_offset, tau_time_width = _component_probabilities(
            state.parameters, joint.data, config.timing, state.layout)
        center = float(state.parameters[state.layout["center"]]); sigma = math.exp(float(state.parameters[state.layout["log_sigma"]]))
        details = mass_details[s]
        for r, run in enumerate(joint.data.runs):
            yields = evaluation.yields[s][r]
            row = {"run_number": int(run), "acquisition_mode": joint.data.mode,
                   "multiplicity_class": joint.data.multiplicity,
                   "observed_count": int(len(joint.event_indices_by_run[r])),
                   "expected_count": float(yields.sum()),
                   "pi0_mass_mean_gev": float(details["means"][r]),
                   "pi0_mass_sigma_gev": float(details["sigmas"][r]),
                   "timing_center_ns": center + float(offsets[r]),
                   "timing_sigma_ns": sigma * math.exp(float(time_widths[r])),
                   "timing_offset_population_scale_ns": tau_offset,
                   "timing_log_width_population_scale": tau_time_width}
            row["yield_em_iterations"] = int(evaluation.yield_iterations[s][r])
            row["yield_kkt_residual"] = float(evaluation.yield_kkt_residuals[s][r])
            row.update({f"yield_{name}": float(yields[k]) for k, name in enumerate(COMPONENTS)})
            run_rows.append(row)
            event_indices = joint.event_indices_by_run[r]
            features = (evaluation.timing_probabilities[s][r, :, joint.event_cell_index[event_indices]] *
                        mass_prob[s][r, :, joint.event_mass_index[event_indices]])
            mu = features @ yields
            keys = np.stack((joint.event_mass_index[event_indices], joint.event_cell_index[event_indices]), axis=1)
            unique, inverse = np.unique(keys, axis=0, return_inverse=True)
            observed = np.bincount(inverse); expected = np.bincount(inverse, weights=mu) / observed
            for k, (mass_index, cell_index) in enumerate(unique):
                predictions.append({"run_number": int(run), "acquisition_mode": joint.data.mode,
                                    "multiplicity_class": joint.data.multiplicity,
                                    "mass_index": int(mass_index), "timing_cell_index": int(cell_index),
                                    "observed": int(observed[k]), "expected": float(expected[k])})
        parameter_rows.extend((
            {"scope": "stratum", "tag": joint.tag, "run_number": "",
             "parameter": "signal_base_mean_gev", "value": details["base_mean"]},
            {"scope": "stratum", "tag": joint.tag, "run_number": "",
             "parameter": "signal_base_sigma_gev", "value": details["base_sigma"]},
            {"scope": "stratum", "tag": joint.tag, "run_number": "",
             "parameter": "timing_base_center_ns", "value": center},
            {"scope": "stratum", "tag": joint.tag, "run_number": "",
             "parameter": "timing_base_sigma_ns", "value": sigma}))
    covariance, labels, profiles = _conditional_covariance(dataset, evaluation)
    observed_total = float(len(dataset.bundle.arrays["run_number"]))
    expected_total = float(sum(x.sum() for x in evaluation.yields))
    closure = {"observed_total": observed_total, "expected_total": expected_total,
               "absolute_difference": abs(expected_total - observed_total)}
    return JointFitResult(dataset, config, mass, layout, timing, evaluation,
                          starts, cycles, run_rows, parameter_rows, predictions,
                          covariance, labels, profiles, closure)


def _validate_shadow_output_path(
    destination: Path, allow_canonical_shadow_output: bool,
) -> bool:
    repo = Path(__file__).resolve().parents[2]
    canonical = (repo / "output").resolve()
    if destination != canonical and canonical not in destination.parents:
        return False
    authorized_root = (canonical / "KinC_x36_4" / "alg002b").resolve()
    _require(
        allow_canonical_shadow_output and authorized_root in destination.parents,
        "ALG-002B canonical output requires explicit authorization and a "
        "child path under output/KinC_x36_4/alg002b",
    )
    return True


def write_joint_result(
    result: JointFitResult, output_dir: Path | str,
    run_manifest: dict[str, object], command: Sequence[str] | None = None,
    allow_canonical_shadow_output: bool = False,
) -> None:
    destination = Path(output_dir).resolve(); repo = Path(__file__).resolve().parents[2]
    canonical_shadow = _validate_shadow_output_path(
        destination, allow_canonical_shadow_output)
    _require(not destination.exists(), f"output path already exists: {destination}")
    output = destination.with_name(f".{destination.name}.partial-{os.getpid()}")
    _require(not output.exists(), f"staging output path already exists: {output}")
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "input_manifest.csv", result.dataset.bundle.manifest)
    _write_csv(output / "run_pi0_yields.csv", result.run_yield_rows)
    _write_csv(output / "parameter_estimates.csv", result.parameter_rows)
    _write_csv(output / "conditional_intervals.csv", result.profile_rows)
    _write_csv(output / "mass_timing_predictions.csv", result.prediction_rows)
    represented = set(run_manifest["represented_runs"]); missing = set(run_manifest["missing_runs"])
    _write_csv(output / "run_coverage.csv", [{"run_number": run, "target": "LH2",
        "status": "represented" if run in represented else
                  "excluded_missing_input" if run in missing else "missing_unclassified"}
        for run in run_manifest["expected_runs"]])
    np.savez_compressed(output / "pi0_yield_covariance.npz",
                        labels=np.asarray(result.covariance_labels), covariance=result.covariance,
                        status=np.asarray("SHAPES_FIXED_CONDITIONAL_NOT_FULL_COVARIANCE"))
    parameter_boundaries = []
    for scope, tag, values, names, bounds in [
        ("mass", "all", result.mass_parameters, result.mass_layout.names,
         result.mass_layout.bounds),
        *(("timing", joint.tag, state.parameters, state.names, state.bounds)
          for joint, state in zip(result.dataset.strata, result.timing_states)),
    ]:
        for value, name, (lower, upper) in zip(values, names, bounds):
            tolerance = max(1.0e-9, 1.0e-6 * (upper - lower))
            side = ("lower" if abs(float(value) - lower) <= tolerance else
                    "upper" if abs(float(value) - upper) <= tolerance else None)
            if side is not None:
                parameter_boundaries.append({
                    "scope": scope, "tag": tag, "parameter": name,
                    "value": float(value), "bound": side,
                })
    zero_yields = sum(
        int(np.count_nonzero(values <= 1.0e-10)) for values in result.evaluation.yields
    )
    maximum_yield_kkt = max(
        (float(np.max(values)) for values in result.evaluation.yield_kkt_residuals),
        default=0.0,
    )
    yield_profiles_above_tolerance = sum(
        int(np.count_nonzero(values > 1.0e-6))
        for values in result.evaluation.yield_kkt_residuals
    )
    (output / "model_selection.json").write_text(json.dumps({
        "provisional_signal_candidate": result.config.signal_model,
        "provisional_combinatorial_candidate": result.config.combinatorial_model,
        "candidate_comparison_status": "PENDING_FULL_LOO_AND_NULL_TOY_CALIBRATION",
        "selection_locked_before_weight_construction": True}, indent=2, sort_keys=True) + "\n")
    provenance = {
        "status": "ALG002B_SHADOW_NOT_PRODUCTION", "publication_ready": False,
        "production_ready": False, "target": "LH2", "lh2_only_enforced": True,
        "output_directory": str(destination),
        "canonical_shadow_output_authorized": canonical_shadow,
        "pi0_weight_written": False, "efficiency_inputs_read": False,
        "charge_scaling_used": False, "cross_section_formed": False,
        "command": list(command) if command is not None else None,
        "working_directory": os.getcwd(), "git_head": _git_head(repo),
        "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "config": asdict(result.config), "run_manifest": run_manifest,
        "objective_without_count_constants": result.evaluation.objective,
        "optimizer_starts": result.start_summaries,
        "optimizer_cycles": result.cycle_summaries, "closure": result.closure,
        "boundary_diagnostics": {
            "parameters_at_bounds": parameter_boundaries,
            "component_yields_at_zero": zero_yields,
        },
        "yield_profile_diagnostics": {
            "maximum_kkt_residual": maximum_yield_kkt,
            "profiles_above_1e-6": yield_profiles_above_tolerance,
        },
        "interval_status": "CONDITIONAL_NORMAL_APPROXIMATION_NOT_PROFILE_INTERVALS",
        "covariance_status": "CONDITIONAL_SHAPES_FIXED; FULL_REPLICA_COVARIANCE_PENDING"}
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True, allow_nan=False) + "\n")
    (output / "FIT_STATUS.txt").write_text(
        "ALG002B_SHADOW_NOT_PRODUCTION\n"
        "LH2 production runs only; run 6569 is an explicit missing-input exclusion.\n"
        "No pi0_weight, efficiency, charge scaling, or cross section is produced.\n"
        "Full model selection, replica covariance, and coverage gates remain required.\n")
    os.rename(output, destination)


def fit_and_write_joint_model(
    bundle: ObservationBundle, output_dir: Path | str, expected_runs: Sequence[int],
    config: JointMassFitConfig | None = None, allowed_missing: Sequence[int] = (6569,),
    command: Sequence[str] | None = None, timing_initial_dir: Path | str | None = None,
    allow_canonical_shadow_output: bool = False,
) -> JointFitResult:
    config = config or JointMassFitConfig(); output = Path(output_dir).resolve()
    _validate_shadow_output_path(output, allow_canonical_shadow_output)
    _require(not output.exists(), f"output path already exists: {output}")
    manifest = enforce_lh2_manifest(bundle, expected_runs, allowed_missing)
    dataset = build_joint_dataset(bundle, config.timing)
    mass, layout, timing, evaluation, starts, cycles = fit_joint_model(
        dataset, config, timing_initial_dir=timing_initial_dir)
    result = _materialize(dataset, config, mass, layout, timing, evaluation, starts, cycles)
    write_joint_result(
        result, output, manifest, command=command,
        allow_canonical_shadow_output=allow_canonical_shadow_output)
    return result
