#!/usr/bin/env python3
"""Compare two diagnostic bundles at the raw_observation branch level."""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import re
from pathlib import Path
from typing import Sequence

import awkward as ak
import numpy as np
import uproot


RUN_PATTERN = re.compile(r"diagnostics_run(\d+)\.root$")


def _expand(expressions: Sequence[str]) -> dict[int, Path]:
    paths: list[Path] = []
    for expression in expressions:
        candidate = Path(expression)
        if candidate.is_dir():
            direct = sorted(candidate.glob("diagnostics_run*.root"))
            paths.extend(direct or sorted(candidate.glob("root/diagnostics_run*.root")))
        else:
            matches = sorted(glob.glob(expression))
            paths.extend(Path(value) for value in matches)
    result: dict[int, Path] = {}
    for path in paths:
        match = RUN_PATTERN.search(path.name)
        if match is None:
            raise ValueError(f"unexpected diagnostic filename: {path}")
        run = int(match.group(1))
        if run in result:
            raise ValueError(f"duplicate run {run}: {result[run]} and {path}")
        result[run] = path.resolve()
    if not result:
        raise ValueError("no diagnostics_run*.root inputs resolved")
    return result


def _array_digest(array: ak.Array) -> str:
    form, length, buffers = ak.to_buffers(array)
    digest = hashlib.sha256()
    digest.update(form.to_json().encode())
    digest.update(str(length).encode())
    for key in sorted(buffers):
        values = np.asarray(buffers[key])
        digest.update(key.encode())
        digest.update(values.dtype.str.encode())
        digest.update(repr(values.shape).encode())
        digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _tree_summary(path: Path, tree_name: str) -> dict[str, object]:
    with uproot.open(path) as root_file:
        if tree_name not in root_file:
            raise ValueError(f"{path} has no {tree_name} tree")
        tree = root_file[tree_name]
        branches = list(tree.keys())
        digests = {
            branch: _array_digest(tree[branch].array(library="ak"))
            for branch in branches
        }
        run_values = np.asarray(tree["run_number"].array(library="np"), dtype=np.int64)
        runs = sorted(set(int(value) for value in run_values))
        return {
            "entries": int(tree.num_entries),
            "branches": branches,
            "branch_digests": digests,
            "run_values": runs,
        }


def compare_bundles(
    reference_inputs: Sequence[str], candidate_inputs: Sequence[str],
    tree_name: str = "raw_observation",
) -> dict[str, object]:
    reference = _expand(reference_inputs)
    candidate = _expand(candidate_inputs)
    reference_runs = set(reference)
    candidate_runs = set(candidate)
    per_run: list[dict[str, object]] = []
    for run in sorted(reference_runs & candidate_runs):
        left = _tree_summary(reference[run], tree_name)
        right = _tree_summary(candidate[run], tree_name)
        left_branches = list(left["branches"])
        right_branches = list(right["branches"])
        common = sorted(set(left_branches) & set(right_branches))
        mismatched = [
            branch for branch in common
            if left["branch_digests"][branch] != right["branch_digests"][branch]
        ]
        missing_candidate = sorted(set(left_branches) - set(right_branches))
        extra_candidate = sorted(set(right_branches) - set(left_branches))
        row_status = (
            "PASS" if left["entries"] == right["entries"] and
            left["run_values"] == [run] and right["run_values"] == [run] and
            not mismatched and not missing_candidate and not extra_candidate
            else "FAIL"
        )
        per_run.append({
            "run_number": run,
            "status": row_status,
            "reference_path": str(reference[run]),
            "candidate_path": str(candidate[run]),
            "reference_entries": left["entries"],
            "candidate_entries": right["entries"],
            "mismatched_branches": mismatched,
            "missing_candidate_branches": missing_candidate,
            "extra_candidate_branches": extra_candidate,
        })
    missing_runs = sorted(reference_runs - candidate_runs)
    extra_runs = sorted(candidate_runs - reference_runs)
    failed_runs = [int(row["run_number"]) for row in per_run if row["status"] != "PASS"]
    status = "PASS" if not missing_runs and not extra_runs and not failed_runs else "FAIL"
    return {
        "status": status,
        "tree_name": tree_name,
        "reference_run_count": len(reference),
        "candidate_run_count": len(candidate),
        "reference_entry_count": sum(int(row["reference_entries"]) for row in per_run),
        "candidate_entry_count": sum(int(row["candidate_entries"]) for row in per_run),
        "missing_candidate_runs": missing_runs,
        "extra_candidate_runs": extra_runs,
        "failed_runs": failed_runs,
        "per_run": per_run,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", action="append", required=True,
                        help="Reference directory, ROOT path, or glob; repeatable.")
    parser.add_argument("--candidate", action="append", required=True,
                        help="Candidate directory, ROOT path, or glob; repeatable.")
    parser.add_argument("--tree", default="raw_observation")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        report = compare_bundles(args.reference, args.candidate, args.tree)
    except (OSError, ValueError, KeyError) as error:
        print(f"raw-observation comparison failed: {error}")
        return 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": report["status"],
        "runs": report["candidate_run_count"],
        "entries": report["candidate_entry_count"],
        "failed_runs": report["failed_runs"],
        "output": str(args.output.resolve()),
    }, indent=2))
    return 0 if report["status"] == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
