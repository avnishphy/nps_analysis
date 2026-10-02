#!/usr/bin/env python3
"""Write explicit pass/fail/pending gates for an ALG-002A shadow bundle."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def stripped_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = []
        for source in reader:
            rows.append({str(key).strip(): str(value).strip()
                         for key, value in source.items()})
        return rows


def expected_runs(path: Path, kin: str, target: str, types: set[str]) -> list[int]:
    runs = []
    for row in stripped_rows(path):
        if row.get("Kin_old", "").lower() != kin.lower():
            continue
        if row.get("target", "").lower() != target.lower():
            continue
        if row.get("Type", "").lower() not in types:
            continue
        runs.append(int(float(row["run_number"])))
    if not runs:
        raise ValueError("no expected runs selected from config CSV")
    return sorted(set(runs))


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def gate(status: str, measured, criterion: str, explanation: str) -> dict[str, object]:
    return {
        "status": status,
        "measured": measured,
        "criterion": criterion,
        "explanation": explanation,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--config-csv", type=Path, required=True)
    parser.add_argument("--kin", required=True)
    parser.add_argument("--target", default="LH2")
    parser.add_argument("--types", default="production")
    args = parser.parse_args()

    output = args.output_dir.resolve()
    provenance = json.loads((output / "provenance.json").read_text())
    components = read_table(output / "component_yields_by_run_mass.csv")
    predictions = read_table(output / "timing_predictions.csv")
    run_parameters = read_table(output / "run_timing_parameters.csv")
    expected = expected_runs(
        args.config_csv.resolve(), args.kin, args.target,
        {value.strip().lower() for value in args.types.split(",") if value.strip()},
    )
    represented = sorted(map(int, provenance["input_run_numbers"]))
    missing = sorted(set(expected) - set(represented))
    unexpected = sorted(set(represented) - set(expected))

    coverage_rows = []
    for run in expected:
        coverage_rows.append({
            "run_number": run,
            "status": "represented" if run in represented else "missing_input",
        })
    for run in unexpected:
        coverage_rows.append({"run_number": run, "status": "unexpected_input"})
    with (output / "run_coverage.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("run_number", "status"))
        writer.writeheader()
        writer.writerows(coverage_rows)

    input_entries = int(provenance["input_entries"])
    component_closure = sum(int(row["observed_count"]) for row in components)
    prediction_closure = sum(int(row["observed"]) for row in predictions)
    expected_total = sum(float(row["expected"]) for row in predictions)
    nonempty = [row for row in components if int(row["observed_count"]) > 0]
    full_rank = sum(int(row["conditional_information_rank"]) == 5 for row in nonempty)
    rank_fraction = full_rank / len(nonempty) if nonempty else 0.0
    zero_true = sum(float(row["yield_true"]) <= 1.0e-10 for row in nonempty)
    offsets = np.asarray([float(row["timing_offset_ns"]) for row in run_parameters])
    width_deviations = np.asarray([float(row["log_width_deviation"])
                                   for row in run_parameters])

    covariance_conditions = []
    covariance_positive = True
    for path in sorted(output.glob("*_profile_covariance.npz")):
        with np.load(path) as arrays:
            eigenvalues = np.asarray(arrays["covariance_eigenvalues"])
            covariance_positive &= bool(np.all(np.isfinite(eigenvalues)) and
                                        np.all(eigenvalues > 0.0))
            covariance_conditions.append(float(eigenvalues.max() / eigenvalues.min()))
    max_condition = max(covariance_conditions, default=math.inf)
    tau_boundaries = []
    for summary in provenance["strata"]:
        parameters = summary["parameters"]
        tau_offset = math.exp(float(parameters["log_run_offset_scale_ns"]))
        tau_width = math.exp(float(parameters["log_run_width_scale"]))
        if tau_offset <= 0.0200001:
            tau_boundaries.append(summary["tag"] + ":offset")
        if tau_width <= 0.0100001:
            tau_boundaries.append(summary["tag"] + ":width")

    forbidden_names = [path.name for path in output.iterdir()
                       if "pi0_weight" in path.name.lower()]
    gates = {
        "configured_run_coverage": gate(
            "pass" if not missing and not unexpected else "fail",
            {"expected": len(expected), "represented": len(represented),
             "missing": missing, "unexpected": unexpected},
            "every configured setting/target/type run is represented exactly once",
            "Missing inputs block a complete-setting claim but do not erase the shadow fit.",
        ),
        "integer_count_closure": gate(
            "pass" if component_closure == prediction_closure == input_entries and
            abs(expected_total - input_entries) < 1.0e-5 else "fail",
            {"input": input_entries, "mass_rows": component_closure,
             "timing_rows": prediction_closure, "fitted_expected": expected_total},
            "raw, mass-stratified, timing-cell, and fitted total counts agree",
            "This protects against silent event loss or duplication.",
        ),
        "optimizer_convergence": gate(
            "pass" if all(item["optimizer_success"] for item in provenance["strata"]) else "fail",
            [{"tag": item["tag"], "iterations": item["optimizer_iterations"],
              "message": item["optimizer_message"]} for item in provenance["strata"]],
            "every stratum reports successful outer optimization",
            "Inner yields additionally require nonnegative KKT-checked convergence.",
        ),
        "profile_curvature": gate(
            "pass" if covariance_positive and max_condition < 1.0e8 else "fail",
            {"positive_definite": covariance_positive,
             "maximum_condition_number": max_condition},
            "diagnostic inverse-LBFGS curvature is positive with condition < 1e8",
            "This is a numerical diagnostic, not calibrated interval coverage.",
        ),
        "run_deviation_parameterization": gate(
            "pass" if abs(float(offsets.sum())) < 1.0e-10 and
            abs(float(width_deviations.sum())) < 1.0e-10 and
            float(np.max(np.abs(offsets))) < 2.0 else "fail",
            {"offset_sum_ns": float(offsets.sum()),
             "max_abs_offset_ns": float(np.max(np.abs(offsets))),
             "log_width_deviation_sum": float(width_deviations.sum())},
            "zero-mean contrasts close and no fitted run offset reaches 2 ns",
            "The balanced orthonormal basis prevents a privileged accumulator run.",
        ),
        "run_population_scale_interior": gate(
            "warning" if tau_boundaries else "pass",
            {"lower_bound_parameters": tau_boundaries},
            "run-to-run population scales are reported at any numerical boundary",
            "A lower boundary is consistent with shared parameters but gives nonregular variance-component inference.",
        ),
        "mass_bin_identifiability": gate(
            "pass" if rank_fraction >= 0.90 and zero_true == 0 else "fail",
            {"occupied_bins": len(nonempty), "full_rank_bins": full_rank,
             "full_rank_fraction": rank_fraction,
             "occupied_bins_with_zero_true_mle": zero_true},
            ">=90% of occupied run/multiplicity/mass bins have rank 5 and none force true yield to zero",
            "Failure blocks event purity weights and motivates the separately approved mass model.",
        ),
        "no_pi0_weight_output": gate(
            "pass" if not forbidden_names and not provenance["pi0_weight_written"] else "fail",
            {"forbidden_files": forbidden_names,
             "provenance_pi0_weight_written": provenance["pi0_weight_written"]},
            "ALG-002A writes no pi0_weight",
            "Timing-pilot outputs cannot be mistaken for production purity weights.",
        ),
        "spline_penalty_cross_validation": gate(
            "pending", None, "held-out control-region selection is completed",
            "The current fixed penalty is predeclared but not yet cross-validated.",
        ),
        "central_tail_null_calibration": gate(
            "pending", None, "optional-tail likelihood-ratio calibration passes synthetic null tests",
            "Observed residuals motivate testing the already-proposed tail, not enabling it post hoc.",
        ),
        "synthetic_profile_coverage_2000": gate(
            "pending", None, "2,000-toy bias, pull, coverage, and failure criteria pass",
            "The regression test checks mechanics only; it is not a coverage ensemble.",
        ),
        "toy_calibrated_run_goodness": gate(
            "pending", None, "per-run/multiplicity goodness-of-fit is toy calibrated",
            "Raw Pearson maxima are not interpreted with asymptotic p-values in sparse cells.",
        ),
    }
    blocking = [name for name, value in gates.items()
                if value["status"] in ("fail", "pending")]
    report = {
        "status": "PASS" if not blocking else "NOT_PROMOTABLE",
        "publication_ready": False,
        "production_ready": False,
        "kinematic_setting": args.kin,
        "target": args.target,
        "types": sorted({value.strip().lower() for value in args.types.split(",") if value.strip()}),
        "blocking_gates": blocking,
        "gates": gates,
    }
    (output / "validation_gates.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(json.dumps({"status": report["status"], "blocking_gates": blocking}, indent=2))
    return 0 if not blocking else 3


if __name__ == "__main__":
    raise SystemExit(main())
