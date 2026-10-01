#!/usr/bin/env python3
"""Build a lightweight diagnostics artifact index from output tree.

This helps phase-7 output hygiene by making run/kin/stage artifacts easy to locate.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List


RUN_RE = re.compile(r"run(\d+)")


def classify(stage_dir: str, file_name: str) -> str:
    suffix = Path(file_name).suffix.lower()
    if stage_dir == "logs":
        return "log"
    if stage_dir == "summary":
        return "csv" if suffix == ".csv" else "summary"
    if suffix == ".root":
        return "root"
    if suffix == ".pdf":
        return "pdf"
    if suffix == ".png":
        return "png"
    return "file"


def find_stage(rel_parts: List[str]) -> str:
    if len(rel_parts) < 2:
        return "unknown"
    if rel_parts[1] in {"root", "logs", "summary", "simulation", "xsec"}:
        return rel_parts[1]
    return "other"


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate output diagnostics index CSV")
    parser.add_argument("--output-base", default="output", help="Output base directory")
    parser.add_argument(
        "--out-file",
        default="",
        help="Output index CSV path (default: <output-base>/diagnostics/diagnostics_index.csv)",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    output_base = (repo_root / args.output_base).resolve()
    if not output_base.exists():
        raise FileNotFoundError(f"Output base not found: {output_base}")

    out_file = (
        (repo_root / args.out_file).resolve()
        if args.out_file
        else output_base / "diagnostics" / "diagnostics_index.csv"
    )
    out_file.parent.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, str]] = []
    for path in output_base.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(output_base)
        parts = rel.parts
        if not parts:
            continue

        kin = parts[0]
        stage = find_stage(list(parts))
        artifact_type = classify(stage, path.name)

        run = ""
        m = RUN_RE.search(path.name)
        if m:
            run = m.group(1)

        rows.append(
            {
                "kin": kin,
                "stage": stage,
                "artifact_type": artifact_type,
                "run": run,
                "relative_path": str(rel),
            }
        )

    rows.sort(key=lambda r: (r["kin"], r["stage"], r["artifact_type"], r["relative_path"]))

    with out_file.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["kin", "stage", "artifact_type", "run", "relative_path"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote diagnostics index: {out_file}")
    print(f"Indexed files: {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
