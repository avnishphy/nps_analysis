#!/usr/bin/env python3
"""Export exact combined 2D mass geometry from an existing combined data ROOT file.

Read-only for the ROOT input. Refit the producer's deterministic geometry and
require exact agreement with both stored event flags before writing the text.
"""
import argparse
from pathlib import Path

import numpy as np
import uproot

from combine_analysis_branches import (
    COMBINED_MASS_CUT_TAG,
    add_combined_2d_mass_cut,
    write_combined_mass_cut_debug_text,
)

COLUMNS = (
    "mpi0_all", "mmiss_all", "pi0_weight", "scale", "is_exclusive",
    "is_exclusive_ellipse_combined", "is_exclusive_mcd_combined",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_root", type=Path)
    parser.add_argument("--out", type=Path, help="Output debug text; default is beside the ROOT input")
    args = parser.parse_args()
    with uproot.open(args.data_root) as source:
        tree = source["physics"]
        missing = [name for name in COLUMNS if name not in tree]
        if missing:
            parser.error(f"Missing combined-data branches: {', '.join(missing)}")
        frame = tree.arrays(list(COLUMNS), library="pd")
    original = {name: frame[name].to_numpy(copy=True) for name in
                ("is_exclusive_ellipse_combined", "is_exclusive_mcd_combined")}
    debug = add_combined_2d_mass_cut(frame)
    if debug is None or not debug["params"].get("valid"):
        parser.error("Combined mass-cut fit is invalid; regenerate the combined data file.")
    for name, stored in original.items():
        mismatch = np.count_nonzero(stored != frame[name].to_numpy())
        if mismatch:
            parser.error(f"{name}: {mismatch} flags differ from the saved data; "
                         "regenerate the combined data file with its producer version.")
        print(f"[OK] {name}: all {len(stored)} stored flags match")
    out = args.out or args.data_root.with_name(
        f"{args.data_root.stem}_{COMBINED_MASS_CUT_TAG}_debug.txt")
    out.parent.mkdir(parents=True, exist_ok=True)
    write_combined_mass_cut_debug_text(debug, args.data_root, out)


if __name__ == "__main__":
    main()
