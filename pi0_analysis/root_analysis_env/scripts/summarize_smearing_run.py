#!/usr/bin/env python3
"""Write compact machine-readable smearing run summaries."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_csv(path: Path):
    if not path.exists():
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream))


def f(row, key, default=0.0):
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--strategy", default="jacobi_global_accept_rollback")
    args = parser.parse_args()

    out_dir = args.out_dir
    sweeps = read_csv(out_dir / "smearing_sweep_history.csv")
    objective = read_csv(out_dir / "smearing_objective_breakdown.csv")
    sections = read_csv(out_dir / "section_map.csv")

    unique_sweeps = {}
    for row in sweeps:
        unique_sweeps.setdefault(row.get("sweep", ""), row)
    sweep_rows = list(unique_sweeps.values())

    runtime_rows = [{
        "strategy": args.strategy,
        "n_sweeps": len(sweep_rows),
        "accepted_sweeps": sum(1 for r in sweep_rows if r.get("candidate_accepted", "0") == "1"),
        "rejected_sweeps": sum(1 for r in sweep_rows if r.get("candidate_accepted", "0") == "0"),
        "repeated_candidates": sum(1 for r in sweep_rows if r.get("repeated_candidate", "0") == "1"),
        "stopped_sweeps": sum(1 for r in sweep_rows if r.get("stop_after_sweep", "0") == "1"),
        "stop_reasons": ";".join(sorted({r.get("stop_reason", "") for r in sweep_rows if r.get("stop_reason", "")})),
        "sweep_runtime_sec": sum(f(r, "runtime_sec") for r in sweep_rows),
    }]
    runtime_path = out_dir / "smearing_runtime_summary.csv"
    with runtime_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(runtime_rows[0].keys()))
        writer.writeheader()
        writer.writerows(runtime_rows)

    best = [r for r in objective if r.get("scope") == "global"]
    best_global = best[-1] if best else {}
    best_params = {
        "strategy": args.strategy,
        "global_objective": f(best_global, "selected_objective"),
        "chi2_mpi0": f(best_global, "chi2_mpi0"),
        "chi2_mmiss": f(best_global, "chi2_mmiss"),
        "chi2_mpgg2": f(best_global, "chi2_mpgg2"),
        "sections": sections,
    }
    with (out_dir / "best_parameters.json").open("w") as stream:
        json.dump(best_params, stream, indent=2)

    config = {
        "strategy": args.strategy,
        "energy_mean_model": "a_plus_bE_plus_clnE",
        "objective_observables": ["mpi0", "mmiss", "mpgg2"],
        "section_map": str(out_dir / "section_map.csv"),
        "sweep_history": str(out_dir / "smearing_sweep_history.csv"),
        "objective_breakdown": str(out_dir / "smearing_objective_breakdown.csv"),
    }
    with (out_dir / "optimization_config.json").open("w") as stream:
        json.dump(config, stream, indent=2)

    print(f"Wrote {runtime_path}")
    print(f"Wrote {out_dir / 'best_parameters.json'}")
    print(f"Wrote {out_dir / 'optimization_config.json'}")


if __name__ == "__main__":
    main()
