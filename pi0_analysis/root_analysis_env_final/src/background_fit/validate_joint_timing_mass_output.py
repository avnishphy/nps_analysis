#!/usr/bin/env python3
"""Write fail-closed scientific gates for one ALG-002B shadow result."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from joint_timing_mass_model import expected_lh2_runs


def gate(status: str, observed: object, requirement: str, interpretation: str) -> dict[str, object]:
    return {
        "status": status,
        "observed": observed,
        "requirement": requirement,
        "interpretation": interpretation,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--config-csv", type=Path,
                        default=Path("config/nps_dvcs_all_kins_main.csv"))
    parser.add_argument("--kin", choices=("KinC_x36_4",), required=True)
    parser.add_argument("--target", choices=("LH2",), default="LH2")
    parser.add_argument("--legacy-comparison", type=Path,
                        help="Directory written by compare_legacy_timing_background.py.")
    args = parser.parse_args()

    output = args.output_dir.resolve()
    provenance = json.loads((output / "provenance.json").read_text())
    model_selection = json.loads((output / "model_selection.json").read_text())
    yields = list(csv.DictReader((output / "run_pi0_yields.csv").open()))
    expected, _ = expected_lh2_runs(args.config_csv, args.kin)
    represented = sorted({int(row["run_number"]) for row in yields})
    missing = sorted(set(expected) - set(represented))
    unexpected = sorted(set(represented) - set(expected))
    closure = provenance["closure"]
    starts = provenance["optimizer_starts"]
    objectives = np.asarray([float(item["objective"]) for item in starts])
    start_yields = np.asarray([float(item.get("pi0_yield_sum", math.nan)) for item in starts])
    finite = bool(np.all(np.isfinite(objectives)))
    converged = bool(starts) and all(bool(item.get("converged")) for item in starts)
    relative_spread = (
        float((objectives.max() - objectives.min()) / max(1.0, abs(objectives.min())))
        if len(objectives) else math.inf
    )
    covariance = np.load(output / "pi0_yield_covariance.npz", allow_pickle=False)
    covariance_status = str(covariance["status"])
    conditional_setting_sigma = math.sqrt(max(0.0, float(covariance["covariance"].sum())))
    yield_spread = (float(start_yields.max() - start_yields.min())
                    if len(start_yields) and np.all(np.isfinite(start_yields)) else math.inf)
    boundaries = provenance.get("boundary_diagnostics", {})
    boundary_count = (len(boundaries.get("parameters_at_bounds", [])) +
                      int(boundaries.get("component_yields_at_zero", 0)))
    yield_profile = provenance.get("yield_profile_diagnostics", {})
    maximum_kkt = float(yield_profile.get("maximum_kkt_residual", math.inf))
    forbidden = sorted(path.name for path in output.iterdir()
                       if "pi0_weight" in path.name.lower())
    legacy_comparison = None
    legacy_comparison_valid = False
    if args.legacy_comparison is not None:
        comparison_path = (args.legacy_comparison.resolve() /
                           "legacy_timing_comparison.json")
        legacy_comparison = json.loads(comparison_path.read_text())
        legacy_comparison_valid = (
            legacy_comparison.get("target") == "LH2" and
            legacy_comparison.get("kinematic_setting") == args.kin and
            legacy_comparison.get("represented_runs") == len(represented) and
            legacy_comparison.get("missing_runs") == [6569] and
            legacy_comparison.get("legacy_timing_estimate_used_by_alg002b") is False and
            legacy_comparison.get("legacy_subtracted_histogram_used_by_alg002b") is False and
            legacy_comparison.get("legacy_timing_histogram_ranges_ns") == [[139.0, 161.0]] and
            legacy_comparison.get("legacy_complete_shifted_support") is True
        )

    gates = {
        "lh2_only_manifest": gate(
            "pass" if provenance.get("target") == "LH2" and
            provenance.get("lh2_only_enforced") is True and not unexpected and
            missing == [6569] else "fail",
            {"expected": len(expected), "represented": len(represented),
             "missing": missing, "unexpected": unexpected,
             "target": provenance.get("target")},
            "only configured KinC_x36_4 production-LH2 runs enter; only run 6569 is excluded",
            "This prevents LD2, dummy, fan-test, or unrelated runs from borrowing shape information.",
        ),
        "integer_count_closure": gate(
            "pass" if float(closure["absolute_difference"]) <= 2.0e-3 else "fail",
            closure,
            "summed fitted component yields close to the raw integer observation count within 0.002",
            "Extended-likelihood closure protects against event loss and duplicated normalization.",
        ),
        "optimizer_finite": gate(
            "pass" if finite else "fail",
            {"starts": len(starts), "objectives": objectives.tolist()},
            "every requested deterministic optimizer start is finite",
            "A finite solution is a software/numerical gate, not evidence of model adequacy.",
        ),
        "optimizer_convergence": gate(
            "pass" if converged else "fail",
            {"starts": len(starts), "converged": [item.get("converged") for item in starts]},
            "every mass and timing coordinate block reports numerical convergence",
            "A finite objective at an iteration limit is retained as evidence but is not a fitted result.",
        ),
        "yield_profile_optimality": gate(
            "pass" if maximum_kkt <= 1.0e-6 else "fail",
            yield_profile,
            "every run/stratum nonnegative-yield profile has KKT residual <= 1e-6",
            "The shared-shape objective is trustworthy only when its inner convex yield solves converge.",
        ),
        "optimizer_reproducibility_20_starts": gate(
            "pass" if len(starts) >= 20 and converged and relative_spread <= 1.0e-6 and
            yield_spread <= 0.1 * conditional_setting_sigma else "pending",
            {"starts": len(starts), "relative_objective_spread": relative_spread,
             "pi0_yield_spread": yield_spread,
             "conditional_setting_pi0_standard_error": conditional_setting_sigma},
            ">=20 dispersed starts agree to relative likelihood 1e-6 and yield checks pass",
            "The initial shadow fit may use fewer starts; promotion may not.",
        ),
        "boundary_calibration": gate(
            "pass" if boundary_count == 0 else "pending",
            boundaries,
            "all parameter/yield boundaries are reported and boundary-dependent intervals pass toy coverage",
            "Boundary solutions require calibrated intervals before the affected yield can be interpreted.",
        ),
        "mass_model_selection": gate(
            "pass" if model_selection.get("candidate_comparison_status") == "PASS" else "pending",
            model_selection,
            "signal/background candidates pass leave-one-run-out prediction and null-toy calibration",
            "The current selected candidate is provisional until this comparison is complete.",
        ),
        "full_cross_run_covariance": gate(
            "pass" if covariance_status == "FULL_REPLICA_COVARIANCE_COVERAGE_VALIDATED" else "pending",
            covariance_status,
            "replica covariance includes shared timing/mass nuisance propagation and calibrated coverage",
            "The written conditional covariance fixes fitted shapes and is not the final uncertainty.",
        ),
        "synthetic_coverage_2000": gate(
            "pending", None,
            "2,000 replicas meet the approved bias, pull-width, and 68/95% coverage tolerances",
            "The mechanics regression does not establish statistical coverage.",
        ),
        "leave_one_run_out_prediction": gate(
            "pending", None,
            "all 56 represented runs pass globally corrected predictive p >= 0.01",
            "This checks whether pooled shapes generalize without absorbing an outlying run.",
        ),
        "mass_timing_factorization": gate(
            "pending", None,
            "predeclared interaction test changes yield by <=max(0.25 sigma,1%)",
            "Pair choice can correlate invariant mass and timing, especially for multiplicity >=3.",
        ),
        "timing_model_calibration": gate(
            "pending", None,
            "spline penalty, central tail, and per-run toy goodness gates pass jointly",
            "ALG-002A timing-shape limitations remain active in the joint model.",
        ),
        "legacy_timing_method_comparison": gate(
            ("pass" if legacy_comparison_valid and
             legacy_comparison.get("status") == "CONVERGED_MODEL_DIAGNOSTIC"
             else "pending" if legacy_comparison_valid or legacy_comparison is None
             else "fail"),
            legacy_comparison,
            "all 56 represented LH2 runs use the 139--161 ns legacy histogram and report its estimate beside a converged ALG-002B prompt-accidental prediction",
            "The original estimator is retained as an independent comparator and never used as an ALG-002B likelihood input.",
        ),
        "no_pi0_weight_output": gate(
            "pass" if not forbidden and provenance.get("pi0_weight_written") is False else "fail",
            {"forbidden_files": forbidden,
             "provenance_pi0_weight_written": provenance.get("pi0_weight_written")},
            "ALG-002B produces yields and covariance but no event pi0_weight",
            "Physics-spectrum weighting remains the separately approved ALG-002C decision.",
        ),
        "no_efficiency_or_cross_section": gate(
            "pass" if provenance.get("efficiency_inputs_read") is False and
            provenance.get("cross_section_formed") is False and
            provenance.get("charge_scaling_used") is False else "fail",
            {key: provenance.get(key) for key in
             ("efficiency_inputs_read", "cross_section_formed", "charge_scaling_used")},
            "no efficiency, charge scaling, or cross section enters ALG-002B",
            "Run yields remain free and efficiencies stay frozen.",
        ),
    }
    blocking = [name for name, value in gates.items()
                if value["status"] in ("fail", "pending")]
    report = {
        "status": "PASS" if not blocking else "NOT_PROMOTABLE",
        "publication_ready": False,
        "production_ready": False,
        "kinematic_setting": args.kin,
        "target": "LH2",
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
