#!/usr/bin/env python3
"""Generate detector-side diagnostics coverage report (HMS and NPS)."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Set

import uproot


HMS_TOKENS = [
    "delta",
    "xptar",
    "yptar",
    "xtar",
    "ytar",
    "xfp",
    "yfp",
    "xpfp",
    "ypfp",
    "s1x",
    "s1y",
    "s2x",
    "s2y",
    "hms",
]

NPS_TOKENS = [
    "clust",
    "ncluster",
    "cluster_x",
    "cluster_y",
    "cluster_e",
    "opening_angle",
    "photon",
    "mpi0",
    "mgg",
    "nps",
]


def contains_token(name: str, tokens: List[str]) -> Set[str]:
    lower = name.lower()
    return {tok for tok in tokens if tok in lower}


def parse_summary_columns(path: Path) -> Set[str]:
    if not path.exists():
        return set()
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        return set(reader.fieldnames or [])


def collect_root_tokens(root_file: Path) -> tuple[Set[str], Set[str]]:
    hms: Set[str] = set()
    nps: Set[str] = set()

    with uproot.open(root_file) as f:
        for key in f.keys():
            obj_name = key.split(";")[0]
            hms |= contains_token(obj_name, HMS_TOKENS)
            nps |= contains_token(obj_name, NPS_TOKENS)

        if "physics" in f:
            for branch in f["physics"].keys():
                hms |= contains_token(branch, HMS_TOKENS)
                nps |= contains_token(branch, NPS_TOKENS)

    return hms, nps


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate HMS/NPS diagnostics coverage CSV from outputs"
    )
    parser.add_argument("--output-base", default="output", help="Output base directory")
    parser.add_argument(
        "--out-file",
        default="",
        help=(
            "Coverage CSV path (default: "
            "<output-base>/diagnostics/diagnostics_coverage_hms_nps.csv)"
        ),
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    output_base = (repo_root / args.output_base).resolve()
    if not output_base.exists():
        raise FileNotFoundError(f"Output base not found: {output_base}")

    out_file = (
        (repo_root / args.out_file).resolve()
        if args.out_file
        else output_base / "diagnostics" / "diagnostics_coverage_hms_nps.csv"
    )
    out_file.parent.mkdir(parents=True, exist_ok=True)

    kin_dirs = [d for d in sorted(output_base.iterdir()) if d.is_dir()]
    rows: List[Dict[str, str]] = []

    for kin_dir in kin_dirs:
        kin = kin_dir.name
        root_dir = kin_dir / "root"
        summary_file = kin_dir / "summary" / "summary_all_runs.csv"

        # Skip non-kin helper directories under output base.
        if not root_dir.exists() and not summary_file.exists():
            continue

        summary_cols = parse_summary_columns(summary_file)

        hms_tokens: Set[str] = set()
        nps_tokens: Set[str] = set()

        if "pass_hms" in summary_cols:
            hms_tokens.add("pass_hms")
        if any(c in summary_cols for c in ("s1x_peak", "s1y_peak", "s2x_peak", "s2y_peak")):
            hms_tokens.add("hms_hodoscope")
        if "pass_hms_nps" in summary_cols:
            nps_tokens.add("pass_hms_nps")
        if any(c in summary_cols for c in ("pi0_mu_MeV", "pi0_sigma_MeV", "pi0_signal_counts")):
            nps_tokens.add("pi0_metrics")

        diagnostics_files = sorted(root_dir.glob("diagnostics_run*.root")) if root_dir.exists() else []
        for root_file in diagnostics_files:
            file_hms, file_nps = collect_root_tokens(root_file)
            hms_tokens |= file_hms
            nps_tokens |= file_nps

        has_hms = len(hms_tokens) >= 2
        has_nps = len(nps_tokens) >= 2

        rows.append(
            {
                "kin": kin,
                "diagnostics_root_files": str(len(diagnostics_files)),
                "summary_exists": "yes" if summary_file.exists() else "no",
                "hms_token_count": str(len(hms_tokens)),
                "nps_token_count": str(len(nps_tokens)),
                "has_hms_diagnostics": "yes" if has_hms else "no",
                "has_nps_diagnostics": "yes" if has_nps else "no",
                "hms_tokens_found": ";".join(sorted(hms_tokens)),
                "nps_tokens_found": ";".join(sorted(nps_tokens)),
            }
        )

    with out_file.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "kin",
                "diagnostics_root_files",
                "summary_exists",
                "hms_token_count",
                "nps_token_count",
                "has_hms_diagnostics",
                "has_nps_diagnostics",
                "hms_tokens_found",
                "nps_tokens_found",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    ok_hms = sum(1 for r in rows if r["has_hms_diagnostics"] == "yes")
    ok_nps = sum(1 for r in rows if r["has_nps_diagnostics"] == "yes")

    print(f"Wrote diagnostics coverage: {out_file}")
    print(f"Kinematics scanned: {len(rows)}")
    print(f"Kinematics with HMS diagnostics: {ok_hms}")
    print(f"Kinematics with NPS diagnostics: {ok_nps}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
