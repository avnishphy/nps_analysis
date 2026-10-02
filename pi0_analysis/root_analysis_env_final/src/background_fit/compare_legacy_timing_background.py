#!/usr/bin/env python3
"""Compare legacy timing accidentals with an ALG-002B shadow fit.

The legacy box estimator is reconstructed independently from the exported
timing-region masks.  Its estimate is never passed back into ALG-002B.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

import numpy as np
import uproot

from joint_timing_mass_model import (
    JointMassFitConfig,
    TimingState,
    _evaluate,
    _initial_parameters,
    _mass_layout,
    _write_csv,
    build_joint_dataset,
    enforce_lh2_manifest,
    expected_lh2_runs,
)
from joint_timing_model import FitError, TimingFitConfig, load_raw_observations


REGION_BITS = {
    "prompt": 1 << 0,
    "diagonal": 1 << 1,
    "horizontal": 1 << 2,
    "vertical": 1 << 3,
    "full1": 1 << 4,
    "full2": 1 << 5,
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FitError(message)


def legacy_estimate_from_region_masks(masks: np.ndarray) -> dict[str, float]:
    """Reproduce the production box-area formula from region membership."""
    counts = {
        name: float(np.count_nonzero(np.asarray(masks, dtype=np.uint32) & bit))
        for name, bit in REGION_BITS.items()
    }
    coin_area = 4.0
    areas = {
        "diagonal": 24.0,
        "horizontal": 24.0,
        "vertical": 24.0,
        "full1": 36.0,
        "full2": 36.0,
    }
    normalized = {
        name: counts[name] * coin_area / area for name, area in areas.items()
    }
    estimate = (normalized["diagonal"] +
                0.5 * (normalized["vertical"] + normalized["horizontal"]) -
                0.5 * (normalized["full1"] + normalized["full2"]))
    variance = (
        counts["diagonal"] * (coin_area / areas["diagonal"]) ** 2 +
        0.25 * counts["vertical"] * (coin_area / areas["vertical"]) ** 2 +
        0.25 * counts["horizontal"] * (coin_area / areas["horizontal"]) ** 2 +
        0.25 * counts["full1"] * (coin_area / areas["full1"]) ** 2 +
        0.25 * counts["full2"] * (coin_area / areas["full2"]) ** 2
    )
    output = {f"raw_mask_{name}": value for name, value in counts.items()}
    output.update({f"raw_mask_normalized_{name}": value
                   for name, value in normalized.items()})
    output["raw_mask_box_formula_estimate"] = estimate
    output["raw_mask_box_formula_standard_error"] = math.sqrt(variance)
    return output


def _load_config(provenance: dict[str, object]) -> JointMassFitConfig:
    values = dict(provenance["config"])
    timing_values = dict(values.pop("timing"))
    timing_values["spline_knots_ns"] = tuple(timing_values["spline_knots_ns"])
    return JointMassFitConfig(timing=TimingFitConfig(**timing_values), **values)


def _parameter_lookup(path: Path) -> dict[tuple[str, str], float]:
    rows = csv.DictReader(path.open())
    return {
        (row["tag"], row["parameter"]): float(row["value"])
        for row in rows if row["scope"] == "optimizer"
    }


def _reconstruct_fit(output: Path, bundle, config: JointMassFitConfig):
    dataset = build_joint_dataset(bundle, config.timing)
    initial_mass, mass_layout = _mass_layout(dataset, config)
    lookup = _parameter_lookup(output / "parameter_estimates.csv")
    mass = np.asarray([lookup[("mass", name)] for name in mass_layout.names])
    _require(len(mass) == len(initial_mass), "mass parameter length mismatch")
    timing_states = []
    for joint in dataset.strata:
        _, names, layout, bounds = _initial_parameters(joint.data, config.timing)
        parameters = np.asarray([lookup[(joint.tag, name)] for name in names])
        timing_states.append(TimingState(parameters, names, layout, bounds))
    evaluation = _evaluate(mass, mass_layout, timing_states, dataset, config)
    return dataset, evaluation


def _read_legacy_regions(manifest: list[dict[str, object]]) -> dict[int, dict[str, object]]:
    output: dict[int, dict[str, object]] = {}
    branches = ("run_number", "timing_region_mask", "shifted_sidebands",
                "acquisition_mode")
    for item in manifest:
        path = Path(str(item["path"]))
        with uproot.open(path) as root_file:
            tree = root_file["raw_observation"]
            missing = sorted(set(branches) - set(tree.keys()))
            _require(not missing, f"{path} lacks legacy-comparison branches: {missing}")
            arrays = tree.arrays(branches, library="np")
            stored_estimate = float(root_file["accidental_est"].member("fVal"))
            stored_error = float(root_file["accidental_err"].member("fVal"))
            stored_prompt = float(root_file["coin_raw"].member("fVal"))
        runs = np.unique(arrays["run_number"].astype(np.int64))
        _require(len(runs) == 1, f"legacy comparison expected one run in {path}")
        run = int(runs[0])
        shifts = np.unique(arrays["shifted_sidebands"].astype(np.int64))
        modes = np.unique(arrays["acquisition_mode"].astype(np.int64))
        _require(len(shifts) == 1 and len(modes) == 1,
                 f"inconsistent timing mode within run {run}")
        output[run] = {
            "shifted_sidebands": int(shifts[0]),
            "acquisition_mode": int(modes[0]),
            "masks": arrays["timing_region_mask"].astype(np.uint32),
            "production_prompt_raw": stored_prompt,
            "production_accidental_estimate": stored_estimate,
            "production_accidental_standard_error": stored_error,
        }
    return output


def compare(output: Path, destination: Path, config_csv: Path, kin: str) -> None:
    _require(output.is_dir(), f"missing ALG-002B output: {output}")
    _require(not destination.exists(), f"comparison output already exists: {destination}")
    provenance = json.loads((output / "provenance.json").read_text())
    _require(provenance.get("target") == "LH2", "fit provenance is not LH2")
    fit_manifest = list(csv.DictReader((output / "input_manifest.csv").open()))
    input_paths = [Path(row["path"]) for row in fit_manifest]
    bundle = load_raw_observations(input_paths)
    expected, _ = expected_lh2_runs(config_csv, kin)
    run_manifest = enforce_lh2_manifest(bundle, expected, (6569,))
    actual_hashes = {str(row["path"]): str(row["sha256"])
                     for row in bundle.manifest}
    recorded_hashes = {str(Path(row["path"]).resolve()): row["sha256"]
                       for row in fit_manifest}
    _require(actual_hashes == recorded_hashes,
             "fit input manifest differs from current raw-observation files")

    config = _load_config(provenance)
    dataset, evaluation = _reconstruct_fit(output, bundle, config)
    recorded_objective = float(provenance["objective_without_count_constants"])
    _require(abs(evaluation.objective - recorded_objective) <= 1.0e-5,
             "reconstructed objective does not match fit provenance")
    legacy = _read_legacy_regions(bundle.manifest)

    edges = config.timing.time_edges
    prompt_axis = (edges[:-1] >= 149.0) & (edges[1:] <= 151.0)
    prompt_cells = np.outer(prompt_axis, prompt_axis).reshape(-1)
    by_run: dict[int, dict[str, float]] = {
        int(run): {
            "joint_prompt_expected": 0.0,
            "joint_prompt_true_coincidence": 0.0,
            "joint_prompt_accidental": 0.0,
            "joint_prompt_accidental_conditional_variance": 0.0,
        }
        for run in dataset.runs
    }
    for s, joint in enumerate(dataset.strata):
        for r, run_value in enumerate(joint.data.runs):
            run = int(run_value)
            timing6 = evaluation.timing_probabilities[s][r]
            prompt_probability = timing6[:, prompt_cells].sum(axis=1)
            yields = evaluation.yields[s][r]
            by_run[run]["joint_prompt_expected"] += float(yields @ prompt_probability)
            by_run[run]["joint_prompt_true_coincidence"] += float(
                yields[:2] @ prompt_probability[:2])
            accidental_coefficients = prompt_probability.copy()
            accidental_coefficients[:2] = 0.0
            by_run[run]["joint_prompt_accidental"] += float(
                yields @ accidental_coefficients)

            event_indices = joint.event_indices_by_run[r]
            features = (
                evaluation.timing_probabilities[s][
                    r, :, joint.event_cell_index[event_indices]
                ] *
                evaluation.mass_probabilities[s][
                    r, :, joint.event_mass_index[event_indices]
                ]
            )
            intensity = np.maximum(features @ yields, 1.0e-300)
            information = (features / intensity[:, None]).T @ (
                features / intensity[:, None])
            covariance = np.linalg.pinv(information, rcond=1.0e-10)
            variance = float(accidental_coefficients @ covariance @
                             accidental_coefficients)
            by_run[run]["joint_prompt_accidental_conditional_variance"] += max(
                variance, 0.0)

    converged = bool(provenance.get("optimizer_starts")) and all(
        bool(item.get("converged")) for item in provenance["optimizer_starts"])
    status = ("CONVERGED_MODEL_DIAGNOSTIC" if converged else
              "UNCONVERGED_MODEL_DIAGNOSTIC_ONLY")
    rows = []
    for run in run_manifest["represented_runs"]:
        legacy_values = legacy[run]
        mask_formula = legacy_estimate_from_region_masks(legacy_values["masks"])
        joint = by_run[run]
        joint_error = math.sqrt(
            joint.pop("joint_prompt_accidental_conditional_variance"))
        legacy_estimate = float(legacy_values["production_accidental_estimate"])
        difference = joint["joint_prompt_accidental"] - legacy_estimate
        row = {
            "run_number": run,
            "target": "LH2",
            "acquisition_mode": legacy_values["acquisition_mode"],
            "shifted_sidebands": legacy_values["shifted_sidebands"],
            "legacy_production_prompt_raw": legacy_values["production_prompt_raw"],
            "legacy_production_accidental_estimate": legacy_estimate,
            "legacy_production_accidental_standard_error":
                legacy_values["production_accidental_standard_error"],
            **mask_formula,
            "raw_mask_minus_stored_legacy_accidental":
                mask_formula["raw_mask_box_formula_estimate"] - legacy_estimate,
            **joint,
            "joint_prompt_accidental_conditional_standard_error": joint_error,
            "joint_minus_legacy_accidental": difference,
            "joint_over_legacy_accidental": (
                joint["joint_prompt_accidental"] /
                legacy_estimate if legacy_estimate != 0.0 else None),
            "status": status,
        }
        rows.append(row)

    legacy_total = sum(float(row["legacy_production_accidental_estimate"])
                       for row in rows)
    raw_mask_total = sum(float(row["raw_mask_box_formula_estimate"])
                         for row in rows)
    joint_total = sum(float(row["joint_prompt_accidental"]) for row in rows)
    legacy_vector = np.asarray(
        [row["legacy_production_accidental_estimate"] for row in rows])
    joint_vector = np.asarray([row["joint_prompt_accidental"] for row in rows])
    correlation = (float(np.corrcoef(legacy_vector, joint_vector)[0, 1])
                   if len(rows) > 1 else None)
    summary = {
        "status": status,
        "publication_ready": False,
        "target": "LH2",
        "kinematic_setting": kin,
        "represented_runs": len(rows),
        "missing_runs": run_manifest["missing_runs"],
        "legacy_timing_estimate_used_by_alg002b": False,
        "legacy_subtracted_histogram_used_by_alg002b": False,
        "shared_inputs": ["t1_ns", "t2_ns", "selected event sample"],
        "shared_conceptual_components": [
            "horizontal", "vertical", "two-random", "diagonal"
        ],
        "legacy_method": (
            "nominal-area box algebra: diagonal + 0.5*(horizontal+vertical) "
            "- 0.5*(full1+full2)"
        ),
        "alg002b_method": (
            "accepted-support-normalized simultaneous timing/mass "
            "extended-Poisson mixture"
        ),
        "legacy_accidental_sum": legacy_total,
        "raw_mask_box_formula_sum": raw_mask_total,
        "raw_mask_minus_stored_legacy_sum": raw_mask_total - legacy_total,
        "joint_prompt_accidental_sum": joint_total,
        "joint_minus_legacy_sum": joint_total - legacy_total,
        "run_level_correlation": correlation,
        "uncertainty_note": (
            "joint errors are shapes-fixed conditional diagnostics; both "
            "estimators use the same events, so their difference uncertainty "
            "requires a joint replica calculation"
        ),
        "support_note": (
            "stored production values are authoritative; the raw-mask formula "
            "also includes exported events in shifted 139--140 and 160--161 ns "
            "sideband slices that lie outside the legacy 140--160 ns histogram"
        ),
    }

    staging = destination.with_name(f".{destination.name}.partial-{os.getpid()}")
    _require(not staging.exists(), f"comparison staging path exists: {staging}")
    staging.mkdir(parents=True)
    _write_csv(staging / "legacy_timing_comparison.csv", rows)
    (staging / "legacy_timing_comparison.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    (staging / "STATUS.txt").write_text(
        status + "\n"
        "The legacy timing estimate is an independent comparator, not an "
        "ALG-002B likelihood input.\n")
    os.rename(staging, destination)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fit_output", type=Path)
    parser.add_argument("comparison_output", type=Path)
    parser.add_argument("--config-csv", type=Path,
                        default=Path("config/nps_dvcs_all_kins_main.csv"))
    parser.add_argument("--kin", choices=("KinC_x36_4",), required=True)
    args = parser.parse_args()
    try:
        compare(args.fit_output.resolve(), args.comparison_output.resolve(),
                args.config_csv.resolve(), args.kin)
    except (FitError, KeyError, OSError, ValueError) as error:
        print(f"legacy timing comparison failed: {error}")
        return 2
    print(f"Legacy timing comparison: {args.comparison_output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
