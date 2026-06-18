#!/usr/bin/env python3
"""Summarize smearing strategy runs into one comparison CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_csv(path: Path):
    if not path.exists():
        return []
    with path.open() as stream:
        return list(csv.DictReader(stream))


def as_float(row, key, default=0.0):
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def summarize(label: str, out_dir: Path, metrics_subdir: str):
    sweep_rows = read_csv(out_dir / "smearing_sweep_history.csv")
    objective_rows = read_csv(out_dir / "smearing_objective_breakdown.csv")
    section_rows = read_csv(out_dir / "section_map.csv")
    metric_rows = read_csv(out_dir / metrics_subdir / "smearing_comparison_metrics.csv")
    metric2d_rows = read_csv(out_dir / metrics_subdir / "smearing_2d_metrics.csv")

    global_rows = [r for r in objective_rows if r.get("scope") == "global"]
    global_row = global_rows[-1] if global_rows else {}

    sweep_seen = {}
    for row in sweep_rows:
        sweep_seen.setdefault(row.get("sweep", ""), row)
    sweeps = list(sweep_seen.values())

    accepted = sum(1 for r in sweeps if r.get("candidate_accepted", r.get("accepted_best", "0")) == "1")
    rejected = sum(1 for r in sweeps if r.get("candidate_accepted", r.get("accepted_best", "0")) == "0")
    repeated = sum(1 for r in sweeps if r.get("repeated_candidate", "0") == "1")
    stopped = sum(1 for r in sweeps if r.get("stop_after_sweep", "0") == "1")
    stop_reasons = ";".join(sorted({r.get("stop_reason", "") for r in sweeps if r.get("stop_reason", "")}))
    runtime = sum(as_float(r, "runtime_sec") for r in sweeps)

    def range_avg(col):
        vals = [as_float(r, col, None) for r in section_rows if r.get(col, "") != ""]
        vals = [v for v in vals if v is not None]
        if not vals:
            return "", "", ""
        return min(vals), max(vals), sum(vals) / len(vals)

    mu_a_min, mu_a_max, mu_a_avg = range_avg("best_mu_a")
    mu_b_min, mu_b_max, mu_b_avg = range_avg("best_mu_b")
    mu_c_min, mu_c_max, mu_c_avg = range_avg("best_mu_c")
    sig_min, sig_max, sig_avg = range_avg("best_sigma")

    out = {
        "strategy": label,
        "out_dir": str(out_dir),
        "global_objective": global_row.get("selected_objective", ""),
        "chi2_mpi0": global_row.get("chi2_mpi0", ""),
        "chi2_mmiss": global_row.get("chi2_mmiss", ""),
        "chi2_mpgg2": global_row.get("chi2_mpgg2", ""),
        "sweeps": len(sweeps),
        "accepted_sweeps": accepted,
        "rejected_sweeps": rejected,
        "repeated_candidates": repeated,
        "stopped_sweeps": stopped,
        "stop_reasons": stop_reasons,
        "runtime_sec_sweeps": runtime,
        "mu_a_min": mu_a_min,
        "mu_a_max": mu_a_max,
        "mu_a_avg": mu_a_avg,
        "mu_b_min": mu_b_min,
        "mu_b_max": mu_b_max,
        "mu_b_avg": mu_b_avg,
        "mu_c_min": mu_c_min,
        "mu_c_max": mu_c_max,
        "mu_c_avg": mu_c_avg,
        "sigma_min": sig_min,
        "sigma_max": sig_max,
        "sigma_avg": sig_avg,
    }

    for row in metric_rows:
        obs = row.get("observable", "")
        if not obs:
            continue
        out[f"{obs}_chi2_unsmeared"] = row.get("chi2_unsmeared", "")
        out[f"{obs}_chi2_smeared"] = row.get("chi2_smeared", "")
        out[f"{obs}_delta_smeared_minus_unsmeared"] = row.get("delta_chi2_smeared_minus_unsmeared", "")
        out[f"{obs}_mean_data"] = row.get("mean_data", "")
        out[f"{obs}_mean_unsmeared"] = row.get("mean_unsmeared", "")
        out[f"{obs}_mean_smeared"] = row.get("mean_smeared", "")
        out[f"{obs}_mean_abs_delta_unsmeared"] = row.get("mean_abs_delta_unsmeared", "")
        out[f"{obs}_mean_abs_delta_smeared"] = row.get("mean_abs_delta_smeared", "")
        out[f"{obs}_std_data"] = row.get("std_data", "")
        out[f"{obs}_std_unsmeared"] = row.get("std_unsmeared", "")
        out[f"{obs}_std_smeared"] = row.get("std_smeared", "")
        out[f"{obs}_std_abs_delta_unsmeared"] = row.get("std_abs_delta_unsmeared", "")
        out[f"{obs}_std_abs_delta_smeared"] = row.get("std_abs_delta_smeared", "")
        out[f"{obs}_peak_mean_data"] = row.get("peak_mean_data", "")
        out[f"{obs}_peak_mean_unsmeared"] = row.get("peak_mean_unsmeared", "")
        out[f"{obs}_peak_mean_smeared"] = row.get("peak_mean_smeared", "")
        out[f"{obs}_peak_mean_abs_delta_unsmeared"] = row.get("peak_mean_abs_delta_unsmeared", "")
        out[f"{obs}_peak_mean_abs_delta_smeared"] = row.get("peak_mean_abs_delta_smeared", "")
        out[f"{obs}_peak_std_data"] = row.get("peak_std_data", "")
        out[f"{obs}_peak_std_unsmeared"] = row.get("peak_std_unsmeared", "")
        out[f"{obs}_peak_std_smeared"] = row.get("peak_std_smeared", "")
        out[f"{obs}_peak_std_abs_delta_unsmeared"] = row.get("peak_std_abs_delta_unsmeared", "")
        out[f"{obs}_peak_std_abs_delta_smeared"] = row.get("peak_std_abs_delta_smeared", "")
        out[f"{obs}_tail_fraction_data"] = row.get("tail_fraction_data", "")
        out[f"{obs}_tail_fraction_unsmeared"] = row.get("tail_fraction_unsmeared", "")
        out[f"{obs}_tail_fraction_smeared"] = row.get("tail_fraction_smeared", "")
        out[f"{obs}_tail_fraction_abs_delta_unsmeared"] = row.get("tail_fraction_abs_delta_unsmeared", "")
        out[f"{obs}_tail_fraction_abs_delta_smeared"] = row.get("tail_fraction_abs_delta_smeared", "")

    for row in metric2d_rows:
        pair = row.get("pair", "")
        if not pair:
            continue
        out[f"{pair}_chi2_unsmeared"] = row.get("chi2_unsmeared", "")
        out[f"{pair}_chi2_smeared"] = row.get("chi2_smeared", "")
        out[f"{pair}_delta_smeared_minus_unsmeared"] = row.get("delta_chi2_smeared_minus_unsmeared", "")

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", required=True,
                        help="Strategy run as label:path")
    parser.add_argument("--metrics-subdir", default="smearing_comparison",
                        help="Per-run subdirectory containing comparison metric CSVs.")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    rows = []
    for item in args.run:
        if ":" not in item:
            raise ValueError("--run must be label:path")
        label, path = item.split(":", 1)
        rows.append(summarize(label, Path(path), args.metrics_subdir))

    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
