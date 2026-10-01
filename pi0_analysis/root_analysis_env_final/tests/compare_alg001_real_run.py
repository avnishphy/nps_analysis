#!/usr/bin/env python3
"""Compare matched default-off and ALG-001 opt-in per-run outputs."""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

import numpy as np
import uproot
from uproot.containers import STLMap, STLVector
from uproot.model import Model


RAW_KEYS = {"raw_observation", "raw_observation_segments"}


def arrays_equal(left: np.ndarray, right: np.ndarray) -> bool:
    return np.array_equal(left, right, equal_nan=True)


def compare_serialized(left, right, location: str) -> None:
    """Recursively compare persisted non-canvas uproot objects."""
    if type(left) is not type(right):
        raise AssertionError(
            f"{location}: type differs ({type(left)} != {type(right)})"
        )
    if isinstance(left, np.ndarray):
        if not arrays_equal(left, right):
            raise AssertionError(f"{location}: array differs")
    elif isinstance(left, STLMap):
        compare_serialized(dict(left), dict(right), location)
    elif isinstance(left, STLVector):
        compare_serialized(list(left), list(right), location)
    elif left.__class__.__name__ == "Model_TString":
        if str(left) != str(right):
            raise AssertionError(f"{location}: TString differs")
    elif isinstance(left, Model):
        if left.classname != right.classname:
            raise AssertionError(f"{location}: ROOT class differs")
        if set(left.all_members) != set(right.all_members):
            raise AssertionError(f"{location}: ROOT member set differs")
        for name in sorted(left.all_members):
            compare_serialized(
                left.all_members[name],
                right.all_members[name],
                f"{location}.{name}",
            )
    elif isinstance(left, dict):
        if set(left) != set(right):
            raise AssertionError(f"{location}: mapping keys differ")
        for name in sorted(left):
            compare_serialized(left[name], right[name], f"{location}.{name}")
    elif isinstance(left, (list, tuple)):
        if len(left) != len(right):
            raise AssertionError(f"{location}: sequence length differs")
        for index, (left_item, right_item) in enumerate(zip(left, right)):
            compare_serialized(
                left_item, right_item, f"{location}[{index}]"
            )
    elif isinstance(left, float):
        if not (left == right or (math.isnan(left) and math.isnan(right))):
            raise AssertionError(f"{location}: {left} != {right}")
    elif isinstance(left, (str, int, bool, type(None), bytes, np.generic)):
        if left != right:
            raise AssertionError(f"{location}: {left!r} != {right!r}")
    elif repr(left) != repr(right):
        raise AssertionError(f"{location}: serialized value differs")


def compare_tree(left, right, location: str) -> int:
    left_branches = set(left.keys())
    right_branches = set(right.keys())
    if left_branches != right_branches:
        raise AssertionError(
            f"{location}: branches differ: "
            f"left-only={sorted(left_branches - right_branches)}, "
            f"right-only={sorted(right_branches - left_branches)}"
        )
    left_arrays = left.arrays(sorted(left_branches), library="np")
    right_arrays = right.arrays(sorted(right_branches), library="np")
    for branch in sorted(left_branches):
        if not arrays_equal(left_arrays[branch], right_arrays[branch]):
            raise AssertionError(f"{location}.{branch}: entries differ")
    return len(left_branches)


def compare_diagnostics(
    default_path: Path,
    optin_path: Path,
    run: int,
    expected_source: str | None,
) -> dict[str, object]:
    default = uproot.open(default_path)
    optin = uproot.open(optin_path)
    default_classes = default.classnames(cycle=False)
    optin_classes = optin.classnames(cycle=False)

    new_keys = set(optin_classes) - set(default_classes)
    if new_keys != RAW_KEYS:
        raise AssertionError(f"unexpected opt-in keys: {sorted(new_keys)}")
    if set(default_classes) != set(optin_classes) - RAW_KEYS:
        raise AssertionError("default/opt-in shared ROOT key sets differ")

    physics_branches = compare_tree(
        default["physics"], optin["physics"], "diagnostics.physics"
    )
    histograms = parameters = named_objects = 0
    for name, classname in default_classes.items():
        left = default[name]
        right = optin[name]
        if classname.startswith(("TH1", "TH2", "TH3")):
            if not arrays_equal(left.values(flow=True), right.values(flow=True)):
                raise AssertionError(f"{name}: bin contents differ")
            left_variances = left.variances(flow=True)
            right_variances = right.variances(flow=True)
            if left_variances is None or right_variances is None:
                if left_variances is not right_variances:
                    raise AssertionError(f"{name}: variance availability differs")
            elif not arrays_equal(left_variances, right_variances):
                raise AssertionError(f"{name}: bin variances differ")
            if left.member("fEntries") != right.member("fEntries"):
                raise AssertionError(f"{name}: fill counts differ")
            for left_axis, right_axis in zip(left.axes, right.axes):
                if not arrays_equal(
                    left_axis.edges(flow=True), right_axis.edges(flow=True)
                ):
                    raise AssertionError(f"{name}: axis edges differ")
            histograms += 1
        elif classname.startswith("TParameter"):
            if left.member("fVal") != right.member("fVal"):
                raise AssertionError(f"{name}: parameter value differs")
            parameters += 1
        elif classname == "TNamed":
            if left.member("fTitle") != right.member("fTitle"):
                raise AssertionError(f"{name}: title differs")
            named_objects += 1

    raw = optin["raw_observation"].arrays(library="np")
    raw_count = raw["event_id"].size
    physics_ids = optin["physics"]["event_id"].array(library="np")
    if raw_count != optin["physics"].num_entries:
        raise AssertionError("raw/physics entry counts differ")
    if not arrays_equal(raw["event_id"], physics_ids):
        raise AssertionError("raw/physics event identifiers or order differ")
    if raw_count and not np.all(np.diff(raw["event_id"]) > 0):
        raise AssertionError("raw event identifiers are not strictly increasing")
    if not np.all(raw["run_number"] == run):
        raise AssertionError("raw run_number differs from requested run")
    if not np.all(raw["source_entry"] >= 0):
        raise AssertionError("raw source_entry contains a negative value")
    if not arrays_equal(raw["pair_dt_ns"], raw["t1_ns"] - raw["t2_ns"]):
        raise AssertionError("stored pair_dt_ns differs from t1_ns - t2_ns")

    shared_raw_physics = {
        "mpi0_all",
        "mmiss_all",
        "mmiss_all_corr",
        "Q2",
        "W",
        "t",
        "tmin",
        "phi",
        "xB",
        "nclust_selected",
        "helicity",
    }
    physics = optin["physics"].arrays(sorted(shared_raw_physics), library="np")
    for branch in sorted(shared_raw_physics):
        if not arrays_equal(raw[branch], physics[branch]):
            raise AssertionError(f"raw/physics {branch} values differ")

    segments = optin["raw_observation_segments"].arrays(library="np")
    segment_count = segments["run_number"].size
    if segment_count == 0:
        raise AssertionError("segment ledger is empty")
    if not np.all(segments["run_number"] == run):
        raise AssertionError("segment ledger run_number differs")
    ordered_segment_numbers = segments["source_tree_number"].tolist()
    if ordered_segment_numbers != list(range(segment_count)):
        raise AssertionError("segment numbers are not unique, ordered, and contiguous")
    segment_numbers = set(ordered_segment_numbers)
    if set(raw["source_tree_number"].tolist()) - segment_numbers:
        raise AssertionError("raw rows reference an absent segment")
    for segment_number, tree_name, source_path in zip(
        ordered_segment_numbers,
        segments["source_tree_name"],
        segments["source_path"],
    ):
        if not os.path.isfile(source_path):
            raise AssertionError(f"segment path does not exist: {source_path}")
        source_tree = uproot.open(source_path)[tree_name]
        row_mask = raw["source_tree_number"] == segment_number
        source_entries = raw["source_entry"][row_mask]
        if source_entries.size and (
            source_entries.min() < 0 or source_entries.max() >= source_tree.num_entries
        ):
            raise AssertionError("raw source_entry is outside its source tree")
        if "g.evnum" in source_tree:
            source_event_numbers = source_tree["g.evnum"].array(library="np")[
                source_entries
            ]
            if not arrays_equal(source_event_numbers, raw["event_number"][row_mask]):
                raise AssertionError("raw event_number differs from source g.evnum")
    if expected_source is not None:
        if segment_count != 1 or segments["source_path"].tolist() != [expected_source]:
            raise AssertionError("segment ledger does not match expected source")
    if segment_count == 1:
        if not arrays_equal(raw["source_entry"], raw["event_id"]):
            raise AssertionError("single-segment source_entry/event_id mismatch")

    categories = raw["timing_category"]
    masks = raw["timing_region_mask"]
    bit_to_category = {1: 1, 2: 2, 4: 3, 8: 4, 16: 5, 32: 6}
    expected_categories = np.empty_like(categories)
    for index, mask in enumerate(masks.tolist()):
        if mask == 0:
            expected_categories[index] = 0
        elif mask in bit_to_category:
            expected_categories[index] = bit_to_category[mask]
        else:
            expected_categories[index] = 7
    if not arrays_equal(categories, expected_categories):
        raise AssertionError("timing category/mask mapping differs")
    if np.count_nonzero(categories == 7):
        raise AssertionError("ambiguous timing categories are present")

    category_values, category_frequencies = np.unique(categories, return_counts=True)
    category_counts = {
        int(category): int(count)
        for category, count in zip(category_values, category_frequencies)
    }
    bit_counts = {
        bit: int(np.count_nonzero(masks & bit)) for bit in (1, 2, 4, 8, 16, 32)
    }
    entry_count = lambda name: int(round(optin[name].member("fEntries")))
    if entry_count(f"h_m_pi0_coin_run{run}") != bit_counts[1]:
        raise AssertionError("prompt category/histogram count differs")
    if entry_count(f"h_m_pi0_acc_run{run}") != raw_count - bit_counts[1]:
        raise AssertionError("accidental histogram count differs")
    for stem, bit in (("h_mgg_diag", 2), ("h_mgg_hor", 4), ("h_mgg_ver", 8)):
        observed = sum(
            entry_count(name) for name in optin_classes if name.startswith(stem)
        )
        if observed != bit_counts[bit]:
            raise AssertionError(f"{stem} timing count differs")
    if entry_count(f"h_mgg_full1_run{run}") != bit_counts[16]:
        raise AssertionError("full-box-1 timing count differs")
    if entry_count(f"h_mgg_full2_run{run}") != bit_counts[32]:
        raise AssertionError("full-box-2 timing count differs")

    return {
        "physics_entries": optin["physics"].num_entries,
        "physics_branches": physics_branches,
        "shared_histograms": histograms,
        "parameters": parameters,
        "named_objects": named_objects,
        "raw_entries": raw_count,
        "segment_rows": segment_count,
        "category_counts": category_counts,
        "mask_bit_counts": bit_counts,
        "event_id_range": (
            (int(raw["event_id"][0]), int(raw["event_id"][-1]))
            if raw_count
            else None
        ),
    }


def compare_auxiliary_root(default_dir: Path, optin_dir: Path, run: int) -> dict[str, int]:
    relative_paths = (
        Path(f"plots/mass_cut_run{run}.root"),
        Path(f"plots/run_{run}/combbg_run{run}_order4_results.root"),
    )
    checked: dict[str, int] = {}
    for relative_path in relative_paths:
        left = uproot.open(default_dir / relative_path)
        right = uproot.open(optin_dir / relative_path)
        left_classes = left.classnames(cycle=False)
        right_classes = right.classnames(cycle=False)
        if left_classes != right_classes:
            raise AssertionError(f"{relative_path}: ROOT key/class sets differ")
        count = 0
        for name, classname in left_classes.items():
            if classname in ("TDirectory", "TCanvas"):
                continue
            compare_serialized(
                left[name], right[name], f"{relative_path}:{name}"
            )
            count += 1
        checked[str(relative_path)] = count
    return checked


def compare_ordinary_files(default_dir: Path, optin_dir: Path) -> int:
    default_files = {
        path.relative_to(default_dir) for path in default_dir.rglob("*") if path.is_file()
    }
    optin_files = {
        path.relative_to(optin_dir) for path in optin_dir.rglob("*") if path.is_file()
    }
    if default_files != optin_files:
        raise AssertionError("default/opt-in output file inventories differ")
    checked = 0
    for relative_path in sorted(default_files):
        if relative_path.suffix.lower() in (".root", ".pdf", ".log"):
            continue
        if (default_dir / relative_path).read_bytes() != (
            optin_dir / relative_path
        ).read_bytes():
            raise AssertionError(f"{relative_path}: file bytes differ")
        checked += 1
    return checked


def compare_log_facts(default_dir: Path, optin_dir: Path, run: int) -> int:
    pattern = re.compile(
        rf"^\[INFO\] Run {run}: (?:"
        r"Added |resolved timing cuts |[0-9]+ entries$|using T\.|"
        r"using helicity good-event charge |dead blocks |"
        r"Weighted pass complete |2D mass cuts done )"
    )

    def selected_lines(directory: Path) -> list[str]:
        log_path = directory / f"logs/analysis_main_run{run}.log"
        return [
            line
            for line in log_path.read_text(errors="replace").splitlines()
            if pattern.match(line)
        ]

    default_lines = selected_lines(default_dir)
    optin_lines = selected_lines(optin_dir)
    if default_lines != optin_lines:
        raise AssertionError(
            "physics-relevant default/opt-in log facts differ:\n"
            f"default={default_lines}\nopt-in={optin_lines}"
        )
    if not default_lines:
        raise AssertionError("no physics-relevant log facts were selected")
    return len(default_lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--default-dir", required=True, type=Path)
    parser.add_argument("--optin-dir", required=True, type=Path)
    parser.add_argument("--run", required=True, type=int)
    parser.add_argument("--expected-source")
    args = parser.parse_args()

    diagnostics = compare_diagnostics(
        args.default_dir / f"root/diagnostics_run{args.run}.root",
        args.optin_dir / f"root/diagnostics_run{args.run}.root",
        args.run,
        args.expected_source,
    )
    auxiliary = compare_auxiliary_root(args.default_dir, args.optin_dir, args.run)
    ordinary_files = compare_ordinary_files(args.default_dir, args.optin_dir)
    log_facts = compare_log_facts(args.default_dir, args.optin_dir, args.run)

    print("ALG001_REALRUN_COMPARE=PASS")
    for name, value in diagnostics.items():
        print(f"{name}={value}")
    print(f"auxiliary_root_objects={auxiliary}")
    print(f"byte_identical_ordinary_files={ordinary_files}")
    print(f"identical_physics_log_facts={log_facts}")


if __name__ == "__main__":
    main()
