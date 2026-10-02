#!/usr/bin/env python3
"""Run and aggregate resumable, independently dispersed ALG-002B starts."""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Sequence

import numpy as np

from joint_timing_mass_model import FitError, _validate_shadow_output_path


def _resolved_inputs(expressions: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    for expression in expressions:
        matches = sorted(glob.glob(expression))
        paths.extend(Path(value).resolve() for value in matches)
    unique = list(dict.fromkeys(paths))
    if not unique:
        raise ValueError("no input ROOT files resolved")
    return unique


def aggregate_campaign(campaign_dir: Path, expected_starts: int) -> dict[str, object]:
    config_path = campaign_dir / "campaign_config.json"
    if not config_path.exists():
        return {
            "status": "NOT_PROMOTABLE", "expected_starts": expected_starts,
            "complete_starts": 0, "all_converged": False,
            "identity_consistent": False,
            "identity_error": "campaign_config.json is missing", "starts": [],
        }
    campaign_config = json.loads(config_path.read_text())
    config_digest = hashlib.sha256(config_path.read_bytes()).hexdigest()
    reference_identity: dict[str, object] | None = None
    rows: list[dict[str, object]] = []
    for start_index in range(expected_starts):
        start_dir = campaign_dir / f"start_{start_index:02d}"
        provenance_path = start_dir / "provenance.json"
        if not provenance_path.exists():
            rows.append({"start_index": start_index, "status": "missing"})
            continue
        provenance = json.loads(provenance_path.read_text())
        summaries = provenance.get("optimizer_starts", [])
        fit_config = provenance.get("config", {})
        manifest_path = start_dir / "input_manifest.csv"
        identity = {
            "git_head": provenance.get("git_head"),
            "run_manifest": provenance.get("run_manifest"),
            "input_manifest_sha256": (
                hashlib.sha256(manifest_path.read_bytes()).hexdigest()
                if manifest_path.exists() else None
            ),
        }
        expected_fit_initial = campaign_config.get("fit_initial_dir")
        identity_matches = (
            len(summaries) == 1 and
            int(summaries[0]["start_index"]) == start_index and
            fit_config.get("starts") == 1 and
            fit_config.get("start_index_offset") == start_index and
            fit_config.get("signal_model") == campaign_config.get("signal_model") and
            fit_config.get("combinatorial_model") == campaign_config.get("combinatorial_model") and
            fit_config.get("mass_gradient_backend") == campaign_config.get("mass_gradient_backend") and
            fit_config.get("seed") == campaign_config.get("seed") and
            fit_config.get("coordinate_cycles") == campaign_config.get("coordinate_cycles") and
            fit_config.get("mass_maxiter") == campaign_config.get("mass_maxiter") and
            fit_config.get("timing_refit_maxiter") == campaign_config.get("timing_refit_maxiter") and
            fit_config.get("nproc") == campaign_config.get("workers_per_start") and
            provenance.get("fit_initial_dir") == expected_fit_initial and
            provenance.get("output_directory") == str(start_dir.resolve()) and
            identity["input_manifest_sha256"] is not None
        )
        if reference_identity is None and identity_matches:
            reference_identity = identity
        elif identity_matches:
            identity_matches = identity == reference_identity
        if not identity_matches:
            rows.append({"start_index": start_index, "status": "identity_mismatch"})
            continue
        summary = summaries[0]
        rows.append({
            "start_index": start_index,
            "status": "complete",
            "objective": float(summary["objective"]),
            "pi0_yield_sum": float(summary["pi0_yield_sum"]),
            "converged": bool(summary["converged"]),
            "output_directory": str(start_dir.resolve()),
        })
    complete = [row for row in rows if row["status"] == "complete"]
    objectives = np.asarray([float(row["objective"]) for row in complete])
    yields = np.asarray([float(row["pi0_yield_sum"]) for row in complete])
    best = min(complete, key=lambda row: float(row["objective"])) if complete else None
    conditional_sigma = math.nan
    if best is not None:
        covariance = np.load(
            Path(str(best["output_directory"])) / "pi0_yield_covariance.npz",
            allow_pickle=False,
        )
        conditional_sigma = math.sqrt(max(0.0, float(covariance["covariance"].sum())))
    relative_spread = (
        float((objectives.max() - objectives.min()) /
              max(1.0, abs(float(objectives.min()))))
        if len(objectives) else math.inf
    )
    yield_spread = float(yields.max() - yields.min()) if len(yields) else math.inf
    all_complete = len(complete) == expected_starts
    all_converged = all(bool(row["converged"]) for row in complete) and all_complete
    identity_consistent = all_complete and reference_identity is not None
    reproducible = (
        expected_starts >= 20 and identity_consistent and all_converged and
        relative_spread <= 1.0e-6 and
        np.isfinite(conditional_sigma) and yield_spread <= 0.1 * conditional_sigma
    )
    return {
        "status": "PASS" if reproducible else "NOT_PROMOTABLE",
        "expected_starts": expected_starts,
        "complete_starts": len(complete),
        "all_converged": all_converged,
        "identity_consistent": identity_consistent,
        "campaign_config_sha256": config_digest,
        "git_head": (reference_identity or {}).get("git_head"),
        "input_manifest_sha256": (
            (reference_identity or {}).get("input_manifest_sha256")),
        "relative_objective_spread": relative_spread,
        "yield_spread": yield_spread,
        "best_conditional_setting_sigma": conditional_sigma,
        "objective_requirement": 1.0e-6,
        "yield_spread_requirement_sigma_fraction": 0.1,
        "best_start": best,
        "starts": rows,
    }


def _write_summary(campaign_dir: Path, summary: dict[str, object]) -> None:
    (campaign_dir / "campaign_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    rows = list(summary["starts"])
    fieldnames = sorted({key for row in rows for key in row})
    with (campaign_dir / "campaign_starts.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--config-csv", type=Path,
                        default=Path("config/nps_dvcs_all_kins_main.csv"))
    parser.add_argument("--kin", choices=("KinC_x36_4",), required=True)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--timing-initial-dir", type=Path, required=True)
    parser.add_argument(
        "--fit-initial-dir", type=Path,
        help="Compatible prior ALG-002B result used as the common warm start.",
    )
    parser.add_argument("--signal-model", choices=("dscb", "double_gaussian"),
                        default="dscb")
    parser.add_argument("--combinatorial-model",
                        choices=("logistic", "bernstein3", "bernstein4"),
                        default="bernstein3")
    parser.add_argument("--mass-gradient-backend", choices=("autograd", "finite"),
                        default="autograd")
    parser.add_argument("--start-count", type=int, default=20)
    parser.add_argument("--parallel-starts", type=int, default=1)
    parser.add_argument("--workers-per-start", type=int, default=8)
    parser.add_argument("--coordinate-cycles", type=int, default=2)
    parser.add_argument("--mass-maxiter", type=int, default=300)
    parser.add_argument("--timing-refit-maxiter", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--allow-canonical-shadow-output", action="store_true")
    args = parser.parse_args()
    if min(args.start_count, args.parallel_starts, args.workers_per_start) < 1:
        parser.error("start and worker counts must be positive")
    cpu_budget = args.parallel_starts * args.workers_per_start
    available = os.cpu_count() or 1
    if cpu_budget > available:
        parser.error(f"requested {cpu_budget} CPUs but only {available} are visible")
    inputs = _resolved_inputs(args.input)
    campaign_dir = args.campaign_dir.resolve()
    try:
        _validate_shadow_output_path(
            campaign_dir / "start_00", args.allow_canonical_shadow_output)
    except FitError as error:
        parser.error(str(error))
    campaign_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "inputs": [str(path) for path in inputs],
        "config_csv": str(args.config_csv.resolve()),
        "kin": args.kin,
        "timing_initial_dir": str(args.timing_initial_dir.resolve()),
        "fit_initial_dir": (str(args.fit_initial_dir.resolve())
                            if args.fit_initial_dir is not None else None),
        "signal_model": args.signal_model,
        "combinatorial_model": args.combinatorial_model,
        "mass_gradient_backend": args.mass_gradient_backend,
        "start_count": args.start_count,
        "parallel_starts": args.parallel_starts,
        "workers_per_start": args.workers_per_start,
        "coordinate_cycles": args.coordinate_cycles,
        "mass_maxiter": args.mass_maxiter,
        "timing_refit_maxiter": args.timing_refit_maxiter,
        "seed": args.seed,
        "cpu_budget": cpu_budget,
    }
    config_path = campaign_dir / "campaign_config.json"
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        parser.error("existing campaign_config.json does not match requested campaign")
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    launcher = Path(__file__).with_name("run_joint_timing_mass_fit.py")
    environment = os.environ.copy()
    environment.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1",
                        "MKL_NUM_THREADS": "1"})

    def run_start(start_index: int) -> tuple[int, int]:
        destination = campaign_dir / f"start_{start_index:02d}"
        provenance = destination / "provenance.json"
        if provenance.exists():
            data = json.loads(provenance.read_text())
            summaries = data.get("optimizer_starts", [])
            if len(summaries) == 1 and int(summaries[0]["start_index"]) == start_index:
                return start_index, 0
            return start_index, 4
        command = [sys.executable, str(launcher)]
        for path in inputs:
            command.extend(("--input", str(path)))
        command.extend((
            "--config-csv", str(args.config_csv), "--kin", args.kin,
            "--output-dir", str(destination), "--timing-initial-dir",
            str(args.timing_initial_dir), "--signal-model", args.signal_model,
            "--combinatorial-model", args.combinatorial_model,
            "--mass-gradient-backend", args.mass_gradient_backend,
            "--starts", "1", "--start-index-offset", str(start_index),
            "--coordinate-cycles", str(args.coordinate_cycles), "--nproc",
            str(args.workers_per_start), "--mass-maxiter", str(args.mass_maxiter),
            "--timing-refit-maxiter", str(args.timing_refit_maxiter),
            "--seed", str(args.seed),
        ))
        if args.fit_initial_dir is not None:
            command.extend(("--fit-initial-dir", str(args.fit_initial_dir)))
        if args.allow_canonical_shadow_output:
            command.append("--allow-canonical-shadow-output")
        with (campaign_dir / f"start_{start_index:02d}.log").open("w") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                    env=environment, check=False)
        return start_index, result.returncode

    failures: list[tuple[int, int]] = []
    with ThreadPoolExecutor(max_workers=args.parallel_starts) as executor:
        futures = {executor.submit(run_start, index): index
                   for index in range(args.start_count)}
        for future in as_completed(futures):
            result = future.result()
            print(f"start {result[0]:02d}: exit {result[1]}", flush=True)
            if result[1] != 0:
                failures.append(result)
    summary = aggregate_campaign(campaign_dir, args.start_count)
    _write_summary(campaign_dir, summary)
    print(json.dumps({key: summary[key] for key in (
        "status", "complete_starts", "all_converged",
        "relative_objective_spread", "yield_spread", "best_start")},
        indent=2, allow_nan=False))
    if failures:
        return 2
    return 0 if summary["status"] == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
