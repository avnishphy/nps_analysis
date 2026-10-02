#!/usr/bin/env python3
"""Command-line launcher for the isolated ALG-002B LH2 shadow fit."""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

from joint_timing_mass_model import (
    JointMassFitConfig,
    expected_lh2_runs,
    fit_and_write_joint_model,
)
from joint_timing_model import FitError, TimingFitConfig, load_raw_observations


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True,
                        help="ROOT path or glob; repeatable.")
    parser.add_argument("--config-csv", type=Path,
                        default=Path("config/nps_dvcs_all_kins_main.csv"))
    parser.add_argument("--kin", choices=("KinC_x36_4",), required=True,
                        help="Approved ALG-002B validation setting.")
    parser.add_argument("--target", choices=("LH2",), default="LH2",
                        help="ALG-002B is intentionally restricted to LH2.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--signal-model", choices=("dscb", "double_gaussian"),
                        default="dscb")
    parser.add_argument("--combinatorial-model",
                        choices=("logistic", "bernstein3", "bernstein4"),
                        default="bernstein3")
    parser.add_argument("--mass-gradient-backend", choices=("autograd", "finite"),
                        default="autograd")
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument(
        "--start-index-offset", type=int, default=0,
        help="Absolute first start ID for deterministic resumable campaigns.",
    )
    parser.add_argument("--coordinate-cycles", type=int, default=2)
    parser.add_argument(
        "--nproc", type=int, default=1,
        help=("CPU workers for parallel numerical gradients (requires SciPy "
              ">=1.16; pass the shell value with --nproc \"$nproc\")."),
    )
    parser.add_argument("--mass-maxiter", type=int, default=300)
    parser.add_argument("--timing-refit-maxiter", type=int, default=100)
    parser.add_argument("--timing-initial-dir", type=Path,
                        help="Directory containing ALG-002A profile-covariance NPZ files.")
    parser.add_argument(
        "--fit-initial-dir", type=Path,
        help="Prior compatible ALG-002B result used as a checked warm start.",
    )
    parser.add_argument(
        "--allow-canonical-shadow-output", action="store_true",
        help=("Permit an isolated result below output/KinC_x36_4/alg002b; "
              "canonical output remains refused otherwise."),
    )
    parser.add_argument("--seed", type=int, default=20261001)
    args = parser.parse_args()

    paths: list[Path] = []
    for expression in args.input:
        matches = sorted(glob.glob(expression))
        if matches:
            paths.extend(Path(value) for value in matches)
        else:
            paths.append(Path(expression))
    try:
        expected, _ = expected_lh2_runs(args.config_csv, args.kin)
        bundle = load_raw_observations(paths)
        config = JointMassFitConfig(
            timing=TimingFitConfig(),
            signal_model=args.signal_model,
            combinatorial_model=args.combinatorial_model,
            mass_gradient_backend=args.mass_gradient_backend,
            starts=args.starts,
            start_index_offset=args.start_index_offset,
            coordinate_cycles=args.coordinate_cycles,
            nproc=args.nproc,
            mass_maxiter=args.mass_maxiter,
            timing_refit_maxiter=args.timing_refit_maxiter,
            seed=args.seed,
        )
        result = fit_and_write_joint_model(
            bundle=bundle,
            output_dir=args.output_dir,
            expected_runs=expected,
            config=config,
            allowed_missing=(6569,),
            command=sys.argv,
            timing_initial_dir=args.timing_initial_dir,
            fit_initial_dir=args.fit_initial_dir,
            allow_canonical_shadow_output=args.allow_canonical_shadow_output,
        )
    except (FitError, ValueError, OSError) as error:
        print(f"ALG-002B fit failed: {error}", file=sys.stderr)
        return 2
    print(f"ALG-002B shadow fit complete: objective={result.evaluation.objective:.12g}")
    print(f"Output: {args.output_dir.resolve()}")
    print("Status: shadow only; run the validator before interpreting yields.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
