#!/usr/bin/env python3
"""Command-line launcher for the opt-in ALG-002A timing pilot."""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

from joint_timing_model import FitError, TimingFitConfig, fit_observation_bundle, load_raw_observations


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit the ALG-002A combined-statistics timing model. This shadow tool "
            "does not write pi0 weights, efficiencies, or cross sections."
        )
    )
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="Input diagnostics ROOT file or glob; repeat as needed.",
    )
    parser.add_argument("--output-dir", required=True,
                        help="New/empty non-production output directory.")
    parser.add_argument("--spline-penalty", type=float, default=10.0,
                        help="Second-difference penalty for log B-spline coefficients.")
    parser.add_argument("--maxiter", type=int, default=500,
                        help="Maximum outer L-BFGS iterations per stratum.")
    parser.add_argument("--tolerance", type=float, default=1.0e-8,
                        help="Outer L-BFGS function/gradient tolerance.")
    return parser.parse_args()


def expand_inputs(patterns: list[str]) -> list[Path]:
    paths: list[Path] = []
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if not matches and Path(pattern).is_file():
            matches = [pattern]
        paths.extend(Path(match) for match in matches)
    return sorted(set(paths))


def main() -> int:
    args = parse_arguments()
    try:
        inputs = expand_inputs(args.input)
        bundle = load_raw_observations(inputs)
        config = TimingFitConfig(
            spline_penalty=args.spline_penalty,
            optimizer_maxiter=args.maxiter,
            optimizer_tolerance=args.tolerance,
        )
        fits = fit_observation_bundle(
            bundle, args.output_dir, config=config, command=sys.argv
        )
        print(f"ALG-002A shadow fit complete: {len(fits)} stratum/strata")
        print(f"Output: {Path(args.output_dir).resolve()}")
        print("Status: SHADOW_NOT_PRODUCTION; no pi0_weight was written")
        return 0
    except (FitError, OSError, ValueError) as error:
        print(f"ALG-002A ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
