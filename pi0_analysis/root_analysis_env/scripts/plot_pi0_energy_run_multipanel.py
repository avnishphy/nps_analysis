#!/usr/bin/env python3
"""Plot pi0 fit trends vs run number for several energy settings."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(os.environ.get("TMPDIR", "/tmp")) / "matplotlib"),
)

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(
    "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/"
    "root_analysis_env_main/output"
)
KINEMATIC = "KinC_x36_5"
DEFAULT_OUTDIR = Path("output/plots/pi0_energy_comparison")
METRICS = [
    ("pi0_mu_MeV", r"$m_{\gamma\gamma}$ mean (MeV)"),
    ("pi0_sigma_MeV", r"$m_{\gamma\gamma}$ width (MeV)"),
    ("pi0_signal_counts", "Pi0 signal counts"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make a compact multipanel comparison of pi0 mean, width, and "
            "signal counts vs run number for KinC_x36_5 energy settings."
        )
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=BASE_DIR,
        help=f"Directory containing {KINEMATIC}_yaopeng_*gev outputs.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory for the PNG and PDF.",
    )
    parser.add_argument(
        "--status",
        default="OK",
        help="Run-status filter. Use 'all' to keep every row.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show the plot interactively after writing files.",
    )
    return parser.parse_args()


def energy_value(token: str) -> float:
    """Convert folder tokens like 06, 1, 13, 16 to 0.6, 1.0, 1.3, 1.6."""
    if "." in token:
        return float(token)
    if len(token) == 1:
        return float(token)
    return float(int(token)) / 10.0


def energy_label(token: str) -> str:
    return f"{energy_value(token):g} GeV"


def find_summary_files(base_dir: Path) -> list[tuple[float, str, Path]]:
    pattern = f"{KINEMATIC}_yaopeng_*gev/{KINEMATIC}/summary/summary_all_runs.csv"
    files: list[tuple[float, str, Path]] = []
    for path in base_dir.glob(pattern):
        match = re.search(r"_yaopeng_([^/]+)gev/", str(path))
        if not match:
            continue
        token = match.group(1)
        try:
            value = energy_value(token)
        except ValueError:
            continue
        files.append((value, energy_label(token), path))
    return sorted(files, key=lambda item: item[0])


def load_summary(path: Path, label: str, status: str) -> pd.DataFrame:
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    df["energy"] = label

    for col in ["run", *[name for name, _ in METRICS]]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=["run"])
    if status.lower() != "all" and "run_status" in df.columns:
        df = df[df["run_status"].astype(str).str.strip() == status]
    return df.sort_values("run")


def plot_multipanel(data: list[pd.DataFrame], outdir: Path, status: str) -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    colors = plt.get_cmap("tab10")
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.6), sharex=True)

    for ax, (metric, ylabel) in zip(axes, METRICS):
        for idx, df in enumerate(data):
            good = df[["run", metric]].dropna()
            if good.empty:
                continue
            ax.plot(
                good["run"],
                good[metric],
                marker="o",
                ms=3.8,
                lw=1.15,
                color=colors(idx % 10),
                label=df["energy"].iloc[0],
            )
        ax.set_ylabel(ylabel)
        ax.grid(True, which="major", ls=":", lw=0.7, alpha=0.55)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel("Run number")
    axes[0].legend(
        title="Energy",
        ncols=3,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.34),
        frameon=False,
    )

    title = "KinC_x36_5 pi0 fit stability"
    if status.lower() != "all":
        title += f" ({status} runs)"
    fig.suptitle(title, y=0.995, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.965))

    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / "pi0_energy_run_multipanel.png"
    pdf = outdir / "pi0_energy_run_multipanel.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print(f"Wrote {png}")
    print(f"Wrote {pdf}")


def main() -> None:
    args = parse_args()
    summary_files = find_summary_files(args.base_dir)
    if not summary_files:
        raise SystemExit(f"No summary CSV files found under {args.base_dir}")

    data = [
        load_summary(path, label, args.status)
        for _, label, path in summary_files
    ]
    data = [df for df in data if not df.empty]
    if not data:
        raise SystemExit("No rows left after filtering.")

    plot_multipanel(data, args.outdir, args.status)
    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
