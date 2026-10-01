#!/usr/bin/env python3
"""Generate modular diagnostics reports for HMS, NPS, and experiment-wide coverage."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
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

HMS_SUMMARY_COLUMNS = ["pass_hms", "s1x_peak", "s1y_peak", "s2x_peak", "s2y_peak"]
NPS_SUMMARY_COLUMNS = ["pass_hms_nps", "pi0_mu_MeV", "pi0_sigma_MeV", "pi0_signal_counts"]


@dataclass
class KinCoverage:
    kin: str
    diagnostics_root_files: int
    summary_exists: bool
    hms_tokens: Set[str]
    nps_tokens: Set[str]
    hms_summary_columns: Set[str]
    nps_summary_columns: Set[str]

    @property
    def has_hms_diagnostics(self) -> bool:
        return len(self.hms_tokens) >= 2 or len(self.hms_summary_columns) >= 2

    @property
    def has_nps_diagnostics(self) -> bool:
        return len(self.nps_tokens) >= 2 or len(self.nps_summary_columns) >= 2

    @property
    def has_experiment_diagnostics(self) -> bool:
        return (
            self.diagnostics_root_files > 0
            and self.summary_exists
            and self.has_hms_diagnostics
            and self.has_nps_diagnostics
        )


def contains_tokens(name: str, tokens: List[str]) -> Set[str]:
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
            hms |= contains_tokens(obj_name, HMS_TOKENS)
            nps |= contains_tokens(obj_name, NPS_TOKENS)

        if "physics" in f:
            for branch in f["physics"].keys():
                hms |= contains_tokens(branch, HMS_TOKENS)
                nps |= contains_tokens(branch, NPS_TOKENS)

    return hms, nps


def collect_coverage(output_base: Path) -> List[KinCoverage]:
    kin_dirs = [d for d in sorted(output_base.iterdir()) if d.is_dir()]
    rows: List[KinCoverage] = []

    for kin_dir in kin_dirs:
        kin = kin_dir.name
        root_dir = kin_dir / "root"
        summary_file = kin_dir / "summary" / "summary_all_runs.csv"

        # Skip non-kin helper directories under output base.
        if not root_dir.exists() and not summary_file.exists():
            continue

        summary_cols = parse_summary_columns(summary_file)

        hms_summary = {c for c in HMS_SUMMARY_COLUMNS if c in summary_cols}
        nps_summary = {c for c in NPS_SUMMARY_COLUMNS if c in summary_cols}

        hms_tokens: Set[str] = set()
        nps_tokens: Set[str] = set()

        diagnostics_files = sorted(root_dir.glob("diagnostics_run*.root")) if root_dir.exists() else []
        for root_file in diagnostics_files:
            file_hms, file_nps = collect_root_tokens(root_file)
            hms_tokens |= file_hms
            nps_tokens |= file_nps

        rows.append(
            KinCoverage(
                kin=kin,
                diagnostics_root_files=len(diagnostics_files),
                summary_exists=summary_file.exists(),
                hms_tokens=hms_tokens,
                nps_tokens=nps_tokens,
                hms_summary_columns=hms_summary,
                nps_summary_columns=nps_summary,
            )
        )

    return rows


def write_hms_report(rows: List[KinCoverage], out_file: Path) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "kin",
                "diagnostics_root_files",
                "summary_exists",
                "has_hms_diagnostics",
                "hms_token_count",
                "hms_tokens_found",
                "hms_summary_columns_found",
            ],
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(
                {
                    "kin": r.kin,
                    "diagnostics_root_files": r.diagnostics_root_files,
                    "summary_exists": "yes" if r.summary_exists else "no",
                    "has_hms_diagnostics": "yes" if r.has_hms_diagnostics else "no",
                    "hms_token_count": len(r.hms_tokens),
                    "hms_tokens_found": ";".join(sorted(r.hms_tokens)),
                    "hms_summary_columns_found": ";".join(sorted(r.hms_summary_columns)),
                }
            )


def write_nps_report(rows: List[KinCoverage], out_file: Path) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "kin",
                "diagnostics_root_files",
                "summary_exists",
                "has_nps_diagnostics",
                "nps_token_count",
                "nps_tokens_found",
                "nps_summary_columns_found",
            ],
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(
                {
                    "kin": r.kin,
                    "diagnostics_root_files": r.diagnostics_root_files,
                    "summary_exists": "yes" if r.summary_exists else "no",
                    "has_nps_diagnostics": "yes" if r.has_nps_diagnostics else "no",
                    "nps_token_count": len(r.nps_tokens),
                    "nps_tokens_found": ";".join(sorted(r.nps_tokens)),
                    "nps_summary_columns_found": ";".join(sorted(r.nps_summary_columns)),
                }
            )


def read_diagnostics_index_counts(index_path: Path) -> Dict[str, int]:
    if not index_path.exists():
        return {}
    counts: Dict[str, int] = {}
    with index_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            kin = (row.get("kin") or "").strip()
            if not kin:
                continue
            counts[kin] = counts.get(kin, 0) + 1
    return counts


def write_experiment_report(rows: List[KinCoverage], out_file: Path, output_base: Path) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    index_counts = read_diagnostics_index_counts(output_base / "diagnostics" / "diagnostics_index.csv")

    with out_file.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "kin",
                "diagnostics_root_files",
                "summary_exists",
                "index_artifact_count",
                "has_hms_diagnostics",
                "has_nps_diagnostics",
                "has_experiment_diagnostics",
            ],
        )
        writer.writeheader()
        for r in rows:
            writer.writerow(
                {
                    "kin": r.kin,
                    "diagnostics_root_files": r.diagnostics_root_files,
                    "summary_exists": "yes" if r.summary_exists else "no",
                    "index_artifact_count": index_counts.get(r.kin, 0),
                    "has_hms_diagnostics": "yes" if r.has_hms_diagnostics else "no",
                    "has_nps_diagnostics": "yes" if r.has_nps_diagnostics else "no",
                    "has_experiment_diagnostics": "yes" if r.has_experiment_diagnostics else "no",
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate modular diagnostics reports for HMS, NPS, and experiment"
    )
    parser.add_argument("--output-base", default="output", help="Output base directory")
    parser.add_argument(
        "--report",
        choices=["hms", "nps", "experiment", "all"],
        default="all",
        help="Report type to generate",
    )
    parser.add_argument("--hms-out", default="", help="Override HMS report CSV path")
    parser.add_argument("--nps-out", default="", help="Override NPS report CSV path")
    parser.add_argument(
        "--experiment-out",
        default="",
        help="Override experiment report CSV path",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    output_base = (repo_root / args.output_base).resolve()
    if not output_base.exists():
        raise FileNotFoundError(f"Output base not found: {output_base}")

    diagnostics_dir = output_base / "diagnostics"
    hms_out = (repo_root / args.hms_out).resolve() if args.hms_out else diagnostics_dir / "hms_diagnostics_report.csv"
    nps_out = (repo_root / args.nps_out).resolve() if args.nps_out else diagnostics_dir / "nps_diagnostics_report.csv"
    experiment_out = (
        (repo_root / args.experiment_out).resolve()
        if args.experiment_out
        else diagnostics_dir / "experiment_diagnostics_report.csv"
    )

    rows = collect_coverage(output_base)

    if args.report in {"hms", "all"}:
        write_hms_report(rows, hms_out)
        print(f"Wrote HMS diagnostics report: {hms_out}")

    if args.report in {"nps", "all"}:
        write_nps_report(rows, nps_out)
        print(f"Wrote NPS diagnostics report: {nps_out}")

    if args.report in {"experiment", "all"}:
        write_experiment_report(rows, experiment_out, output_base)
        print(f"Wrote experiment diagnostics report: {experiment_out}")

    print(f"Kinematics scanned: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
