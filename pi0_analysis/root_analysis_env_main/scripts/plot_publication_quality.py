#!/usr/bin/env python3
"""
Publication-quality plots for pi0 analysis with prescale corrections.

Inputs are the current workflow CSVs:
  - summary_all_runs.csv from an analysis summary directory
  - efficiency_<kin>.csv from output/efficiency_stuff

Plots are written to the plots directory beside the summary directory. For
example, a summary CSV at <kin>/summary/summary_all_runs.csv writes to
<kin>/plots.
"""
# python3 scripts/plot_publication_quality.py --summary-csv output/KinC_x36_5_yaopeng_08gev/KinC_x36_5/summary/summary_all_runs.csv --efficiency-csv output/efficiency_stuff/efficiency_KinC_x36_5.csv --run-status OK

from pathlib import Path
from typing import Iterable
import argparse
import sys
import warnings

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


DEFAULT_SUMMARY_CSV = (
    "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/"
    "root_analysis_env_main/output/KinC_x36_5_yaopeng_08gev/KinC_x36_5/"
    "summary/summary_all_runs.csv"
)
DEFAULT_EFFICIENCY_CSV = (
    "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/"
    "root_analysis_env_main/output/efficiency_stuff/efficiency_KinC_x36_5.csv"
)

FIGURE_DPI = 300
FIGURE_SIZE_SINGLE = (10, 6)
FIGURE_SIZE_COMPARISON = (14, 10)
FONT_SIZE_TITLE = 16
FONT_SIZE_LABEL = 13
FONT_SIZE_LEGEND = 11
FONT_SIZE_TICK = 11
COLOR_MAIN = "#1f77b4"
EDGE_COLOR = "black"
EDGE_WIDTH = 0.5
MARKER_SIZE = 70
ALPHA_SCATTER = 0.75


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create publication-quality pi0 analysis plots."
    )
    parser.add_argument(
        "--summary-csv",
        default=DEFAULT_SUMMARY_CSV,
        help="Path to summary_all_runs.csv.",
    )
    parser.add_argument(
        "--efficiency-csv",
        default=DEFAULT_EFFICIENCY_CSV,
        help="Path to efficiency_<kin>.csv.",
    )
    parser.add_argument(
        "--run-status",
        default=None,
        help="Optional exact run_status value to keep from the summary CSV.",
    )
    return parser.parse_args()


def plots_dir_from_summary(summary_csv: Path) -> Path:
    summary_dir = summary_csv.parent
    if summary_dir.name == "summary":
        return summary_dir.parent / "plots"
    return summary_dir / "plots"


def strip_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]
    return df


def require_columns(df: pd.DataFrame, path: Path, columns: Iterable[str]) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"{path} is missing required columns: {joined}")


def to_numeric(df: pd.DataFrame, columns: Iterable[str]) -> None:
    for col in columns:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")


def parse_run_ignore_text(text: str) -> set[int]:
    ignored: set[int] = set()
    for token in text.replace(" ", "").split(","):
        if not token:
            continue
        if "-" in token:
            lo_text, hi_text = token.split("-", 1)
            lo = int(lo_text)
            hi = int(hi_text)
            if hi < lo:
                lo, hi = hi, lo
            ignored.update(range(lo, hi + 1))
        else:
            ignored.add(int(token))
    return ignored


def prompt_for_ignored_runs() -> set[int]:
    prompt = (
        "Runs/ranges to ignore "
        "(example: 3728,3729,4300-4310; blank for none): "
    )
    try:
        text = input(prompt)
    except EOFError:
        text = ""
    text = text.strip()
    if not text:
        return set()
    try:
        return parse_run_ignore_text(text)
    except ValueError as exc:
        print(f"[ERROR] Could not parse ignored runs: {exc}", file=sys.stderr)
        sys.exit(2)


def setup_publication_plot(figsize=FIGURE_SIZE_SINGLE):
    fig, ax = plt.subplots(figsize=figsize, dpi=FIGURE_DPI)
    ax.grid(True, linestyle="--", alpha=0.3, linewidth=0.7)
    ax.set_axisbelow(True)
    return fig, ax


def format_publication_plot(ax, xlabel, ylabel, title, tight=True):
    ax.set_xlabel(xlabel, fontsize=FONT_SIZE_LABEL, weight="bold")
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL, weight="bold")
    ax.set_title(title, fontsize=FONT_SIZE_TITLE, weight="bold", pad=15)
    ax.tick_params(labelsize=FONT_SIZE_TICK)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    if tight:
        plt.tight_layout()


def finite_mask(x, y):
    return ~(np.isnan(x) | np.isnan(y))


def scatter_plot(
    ax,
    x,
    y,
    xlabel,
    ylabel,
    title,
    filename,
    output_dir,
    pdf=None,
    ylim=None,
):
    mask = finite_mask(x, y)
    x_clean = x[mask]
    y_clean = y[mask]

    if len(x_clean) == 0:
        print(f"[WARN] No valid data for plot: {filename}")
        plt.close(ax.figure)
        return 0

    ax.scatter(
        x_clean,
        y_clean,
        s=MARKER_SIZE,
        color=COLOR_MAIN,
        edgecolor=EDGE_COLOR,
        linewidth=EDGE_WIDTH,
        alpha=ALPHA_SCATTER,
        zorder=3,
    )
    format_publication_plot(ax, xlabel, ylabel, title)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.figure.savefig(output_dir / filename, dpi=FIGURE_DPI, bbox_inches="tight")
    if pdf is not None:
        pdf.savefig(ax.figure, bbox_inches="tight")
    plt.close(ax.figure)
    print(f"[SAVE] {filename}")
    return 1


def multi_scatter_plot(
    data_dict,
    xlabel,
    ylabel,
    title,
    filename,
    output_dir,
    figsize=FIGURE_SIZE_COMPARISON,
    pdf=None,
):
    n_plots = len(data_dict)
    cols = min(2, n_plots)
    rows = (n_plots + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=figsize, dpi=FIGURE_DPI)
    if rows == 1 and cols == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    any_data = False
    for idx, (label, (x, y)) in enumerate(data_dict.items()):
        ax = axes[idx]
        mask = finite_mask(x, y)
        x_clean = x[mask]
        y_clean = y[mask]

        if len(x_clean) > 0:
            any_data = True
            ax.scatter(
                x_clean,
                y_clean,
                s=MARKER_SIZE,
                color=COLOR_MAIN,
                edgecolor=EDGE_COLOR,
                linewidth=EDGE_WIDTH,
                alpha=ALPHA_SCATTER,
                zorder=3,
            )
        ax.grid(True, linestyle="--", alpha=0.3, linewidth=0.7)
        ax.set_axisbelow(True)
        ax.set_xlabel(xlabel, fontsize=FONT_SIZE_LABEL, weight="bold")
        ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL - 1, weight="bold")
        ax.set_title(label, fontsize=FONT_SIZE_LABEL, weight="bold")
        ax.tick_params(labelsize=FONT_SIZE_TICK - 1)
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)

    for idx in range(len(data_dict), len(axes)):
        axes[idx].set_visible(False)

    if not any_data:
        print(f"[WARN] No valid data for plot: {filename}")
        plt.close(fig)
        return 0

    fig.suptitle(title, fontsize=FONT_SIZE_TITLE, weight="bold", y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(output_dir / filename, dpi=FIGURE_DPI, bbox_inches="tight")
    if pdf is not None:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVE] {filename}")
    return 1


def comparison_line_plot(
    df,
    x_col,
    y_cols_dict,
    xlabel,
    ylabel,
    title,
    filename,
    output_dir,
    pdf=None,
    ylim=None,
):
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE, dpi=FIGURE_DPI)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
    any_data = False

    for idx, (label, col) in enumerate(y_cols_dict.items()):
        y = df[col].to_numpy(dtype=float)
        mask = ~np.isnan(y)
        x_clean = df[x_col].to_numpy(dtype=float)[mask]
        y_clean = y[mask]

        if len(x_clean) > 0:
            any_data = True
            ax.plot(
                x_clean,
                y_clean,
                marker="o",
                linewidth=2,
                markersize=6,
                label=label,
                color=colors[idx % len(colors)],
                alpha=0.8,
                zorder=2,
            )

    if not any_data:
        print(f"[WARN] No valid data for plot: {filename}")
        plt.close(fig)
        return 0

    ax.grid(True, linestyle="--", alpha=0.3, linewidth=0.7)
    ax.set_axisbelow(True)
    format_publication_plot(ax, xlabel, ylabel, title)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.legend(fontsize=FONT_SIZE_LEGEND, loc="best", framealpha=0.95)
    fig.savefig(output_dir / filename, dpi=FIGURE_DPI, bbox_inches="tight")
    if pdf is not None:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVE] {filename}")
    return 1


def main() -> int:
    args = parse_args()
    summary_csv = Path(args.summary_csv).expanduser().resolve()
    efficiency_csv = Path(args.efficiency_csv).expanduser().resolve()
    output_dir = plots_dir_from_summary(summary_csv)

    if not summary_csv.exists():
        print(f"[ERROR] Summary CSV not found: {summary_csv}", file=sys.stderr)
        return 1
    if not efficiency_csv.exists():
        print(f"[ERROR] Efficiency CSV not found: {efficiency_csv}", file=sys.stderr)
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)

    print("[INFO] Loading input CSV files...")
    summary_df = strip_columns(pd.read_csv(summary_csv, skipinitialspace=True))
    efficiency_df = strip_columns(pd.read_csv(efficiency_csv, skipinitialspace=True))
    print(f"[INFO] Loaded {len(summary_df)} rows from summary CSV: {summary_csv}")
    print(f"[INFO] Loaded {len(efficiency_df)} rows from efficiency CSV: {efficiency_csv}")
    print(f"[INFO] Plot output directory: {output_dir}")

    summary_columns = [
        "run",
        "accumulated_charge(mC)",
        "current_mean_uA",
        "Beam_Time(s)",
        "total_entries",
        "total_coin_entries",
        "chi2_ndf_comb_bg",
        "pi0_mu_MeV",
        "pi0_sigma_MeV",
        "pi0_signal_counts",
        "run_status",
    ]
    efficiency_columns = [
        "run_number",
        "HEL_charge_after_cut_uC",
        "HMS_tracking_eff",
        "NewGen_EDTM_livetime",
        "ps_factor",
        "beam_time",
    ]
    require_columns(summary_df, summary_csv, summary_columns)
    require_columns(efficiency_df, efficiency_csv, efficiency_columns)

    summary_df = summary_df.copy()
    efficiency_df = efficiency_df.copy()

    summary_df["run"] = pd.to_numeric(summary_df["run"], errors="coerce")
    summary_df = summary_df.dropna(subset=["run"])
    summary_df["run"] = summary_df["run"].astype(int)

    efficiency_df["run_number"] = pd.to_numeric(
        efficiency_df["run_number"], errors="coerce"
    )
    efficiency_df = efficiency_df.dropna(subset=["run_number"])
    efficiency_df["run_number"] = efficiency_df["run_number"].astype(int)

    numeric_cols = [
        "accumulated_charge(mC)",
        "current_mean_uA",
        "Beam_Time(s)",
        "total_entries",
        "total_coin_entries",
        "chi2_ndf_comb_bg",
        "pi0_mu_MeV",
        "pi0_sigma_MeV",
        "pi0_signal_counts",
        "HEL_charge_after_cut_uC",
        "HMS_tracking_eff",
        "NewGen_EDTM_livetime",
        "ps_factor",
        "beam_time",
    ]
    to_numeric(summary_df, numeric_cols)
    to_numeric(efficiency_df, numeric_cols)

    df = pd.merge(
        summary_df,
        efficiency_df,
        left_on="run",
        right_on="run_number",
        how="inner",
    )
    print(f"[INFO] {len(df)} runs present in both CSVs")

    suspicious_runs = set()
    suspicious_runs.update(df[df["NewGen_EDTM_livetime"] < 0]["run"].unique())
    for col in df.columns:
        if (
            col.startswith("scaler_") or col.startswith("trig_")
        ) and pd.api.types.is_numeric_dtype(df[col]):
            suspicious_runs.update(df[df[col] > 10_000_000]["run"].unique())

    if suspicious_runs:
        df = df[~df["run"].isin(suspicious_runs)].copy()
        print(f"[INFO] Removed {len(suspicious_runs)} suspicious runs")

    ignored_runs = prompt_for_ignored_runs()
    if ignored_runs:
        before = len(df)
        df = df[~df["run"].isin(ignored_runs)].copy()
        print(f"[INFO] Ignored requested runs: {before} -> {len(df)} rows")

    if args.run_status:
        before = len(df)
        df = df[df["run_status"] == args.run_status].copy()
        print(f"[INFO] Filtered run_status={args.run_status}: {before} -> {len(df)} rows")

    if df.empty:
        print("[ERROR] No rows remain after filtering", file=sys.stderr)
        return 1

    df = df.sort_values("run").copy()
    print("[INFO] Calculating rates and metrics...")

    charge_mC = df["HEL_charge_after_cut_uC"] / 1000.0
    livetime = df["NewGen_EDTM_livetime"]
    hms_eff = df["HMS_tracking_eff"]
    target_boiling_correction = 1 - ((10.2 / 100) * df["current_mean_uA"]) / 100

    df["raw_rate"] = np.where(
        df["beam_time"] > 0,
        df["total_coin_entries"] / df["beam_time"],
        np.nan,
    )
    df["lcf_corrected_rate"] = np.where(
        (df["beam_time"] > 0) & (livetime > 0),
        df["total_coin_entries"] / (df["beam_time"] * livetime),
        np.nan,
    )
    df["current_normalized_rate"] = np.where(
        df["current_mean_uA"] > 0,
        df["lcf_corrected_rate"] / df["current_mean_uA"],
        np.nan,
    )
    df["scaled_current_normalized_rate"] = df["current_normalized_rate"] * (
        df["pi0_signal_counts"] / df["total_coin_entries"]
    )
    df["normalized_yield"] = np.where(
        (charge_mC > 0)
        & (livetime > 0)
        & (hms_eff > 0)
        & (target_boiling_correction > 0),
        (
            df["pi0_signal_counts"]
            / (charge_mC * livetime * hms_eff)
            * target_boiling_correction
        )
        * df["ps_factor"],
        np.nan,
    )
    df["signal_to_coin_ratio"] = np.where(
        df["total_coin_entries"] > 0,
        df["pi0_signal_counts"] / df["total_coin_entries"],
        np.nan,
    )
    df["detection_efficiency"] = np.where(
        df["total_entries"] > 0,
        df["pi0_signal_counts"] / df["total_entries"],
        np.nan,
    )

    for col in ["raw_rate", "lcf_corrected_rate", "normalized_yield"]:
        print(f"  - {col}: min={df[col].min():.2e}, max={df[col].max():.2e}")

    pdf_path = output_dir / "all_plots.pdf"
    plot_count = 0
    print(f"\n[INFO] Creating plots and saving to PDF: {pdf_path}")

    with PdfPages(pdf_path) as pdf:
        print("[INFO] Creating individual rate plots...")
        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["raw_rate"].values,
            "Run Number",
            "Raw Rate [Hz]",
            "Raw Rate vs Run",
            "01_raw_rate_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["lcf_corrected_rate"].values,
            "Run Number",
            "LCF Corrected Rate [Hz]",
            "LCF Corrected Rate vs Run",
            "02_lcf_corrected_rate_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["current_normalized_rate"].values,
            "Run Number",
            "Current Normalized Rate [Hz/uA]",
            "Current Normalized Rate vs Run",
            "03_current_normalized_rate_vs_run.png",
            output_dir,
            pdf=pdf,
            ylim=(0, 0.3),
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["scaled_current_normalized_rate"].values,
            "Run Number",
            "Scaled Current Norm. Rate",
            "Scaled Current Normalized Rate vs Run",
            "04_scaled_current_normalized_rate_vs_run.png",
            output_dir,
            pdf=pdf,
            ylim=(0, 0.3),
        )

        print("[INFO] Creating rate vs current plots...")
        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["raw_rate"].values,
            "Beam Current [uA]",
            "Raw Rate [Hz]",
            "Raw Rate vs Beam Current",
            "05_raw_rate_vs_current.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["lcf_corrected_rate"].values,
            "Beam Current [uA]",
            "LCF Corrected Rate [Hz]",
            "LCF Corrected Rate vs Beam Current",
            "06_lcf_corrected_rate_vs_current.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["current_normalized_rate"].values,
            "Beam Current [uA]",
            "Current Normalized Rate [Hz/uA]",
            "Current Normalized Rate vs Beam Current",
            "07_current_normalized_rate_vs_current.png",
            output_dir,
            pdf=pdf,
            ylim=(0, 0.3),
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["scaled_current_normalized_rate"].values,
            "Beam Current [uA]",
            "Scaled Current Norm. Rate",
            "Scaled Current Normalized Rate vs Beam Current",
            "08_scaled_current_normalized_rate_vs_current.png",
            output_dir,
            pdf=pdf,
            ylim=(0, 0.3),
        )

        print("[INFO] Creating rate comparison plots...")
        plot_count += comparison_line_plot(
            df,
            "run",
            {
                "Raw Rate": "raw_rate",
                "LCF Corrected Rate": "lcf_corrected_rate",
                "Current Normalized Rate": "current_normalized_rate",
                "Scaled Norm. Rate": "scaled_current_normalized_rate",
            },
            "Run Number",
            "Rate [various units]",
            "Rate Metrics Comparison vs Run",
            "09_rate_comparison_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        df_norm_temp = df.copy()
        for col in [
            "raw_rate",
            "lcf_corrected_rate",
            "current_normalized_rate",
            "scaled_current_normalized_rate",
        ]:
            denom = df_norm_temp[col].max() - df_norm_temp[col].min()
            if denom > 0:
                df_norm_temp[col] = (df_norm_temp[col] - df_norm_temp[col].min()) / denom

        plot_count += comparison_line_plot(
            df_norm_temp,
            "run",
            {
                "Raw Rate": "raw_rate",
                "LCF Corrected Rate": "lcf_corrected_rate",
                "Current Normalized Rate": "current_normalized_rate",
                "Scaled Norm. Rate": "scaled_current_normalized_rate",
            },
            "Run Number",
            "Normalized Rate (0-1)",
            "Normalized Rate Metrics Comparison vs Run",
            "10_rate_comparison_normalized_vs_run.png",
            output_dir,
            pdf=pdf,
            ylim=(0, 2),
        )

        print("[INFO] Creating normalized yield plots...")
        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["normalized_yield"].values,
            "Run Number",
            "Normalized Yield (Prescale Corrected)",
            "Normalized Yield vs Run",
            "11_normalized_yield_vs_run.png",
            output_dir,
            pdf=pdf,
            ylim=(45, 110),
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["normalized_yield"].values,
            "Beam Current [uA]",
            "Normalized Yield (Prescale Corrected)",
            "Normalized Yield vs Beam Current",
            "12_normalized_yield_vs_current.png",
            output_dir,
            pdf=pdf,
            ylim=(45, 110),
        )

        print("[INFO] Creating pi0 peak and resolution trend plots...")
        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["pi0_mu_MeV"].values,
            "Run Number",
            "pi0 Peak Position mu [MeV]",
            "pi0 Peak Position vs Run",
            "13_pi0_peak_position_vs_run.png",
            output_dir,
            pdf=pdf,
            ylim=(130, 135),
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["pi0_sigma_MeV"].values,
            "Run Number",
            "pi0 Peak Width sigma [MeV]",
            "pi0 Peak Width vs Run",
            "14_pi0_peak_width_vs_run.png",
            output_dir,
            pdf=pdf,
            ylim=(3.5, 6.5),
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["pi0_mu_MeV"].values,
            "Beam Current [uA]",
            "pi0 Peak Position mu [MeV]",
            "pi0 Peak Position vs Beam Current",
            "15_pi0_peak_position_vs_current.png",
            output_dir,
            pdf=pdf,
            ylim=(130, 135),
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["pi0_sigma_MeV"].values,
            "Beam Current [uA]",
            "pi0 Peak Width sigma [MeV]",
            "pi0 Peak Width vs Beam Current",
            "16_pi0_peak_width_vs_current.png",
            output_dir,
            pdf=pdf,
            ylim=(3.5, 6.5),
        )

        print("[INFO] Creating fit quality plots...")
        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["chi2_ndf_comb_bg"].values,
            "Run Number",
            "chi2 / NDF (Combinatorial BG Fit)",
            "Fit Quality vs Run",
            "17_chi2_ndf_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["chi2_ndf_comb_bg"].values,
            "Beam Current [uA]",
            "chi2 / NDF (Combinatorial BG Fit)",
            "Fit Quality vs Beam Current",
            "18_chi2_ndf_vs_current.png",
            output_dir,
            pdf=pdf,
        )

        print("[INFO] Creating signal and efficiency trend plots...")
        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["pi0_signal_counts"].values,
            "Run Number",
            "pi0 Signal Counts",
            "pi0 Signal Counts vs Run",
            "19_pi0_signal_counts_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["pi0_signal_counts"].values,
            "Beam Current [uA]",
            "pi0 Signal Counts",
            "pi0 Signal Counts vs Beam Current",
            "19b_pi0_signal_counts_vs_current.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["signal_to_coin_ratio"].values,
            "Run Number",
            "Signal to Coincidence Ratio",
            "Signal/Coincidence Ratio vs Run",
            "20_signal_to_coin_ratio_vs_run.png",
            output_dir,
            pdf=pdf,
            ylim=(0, 1),
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["detection_efficiency"].values,
            "Run Number",
            "Detection Efficiency (Signal/Total Entries)",
            "Detection Efficiency vs Run",
            "21_detection_efficiency_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            df["signal_to_coin_ratio"].values,
            "Beam Current [uA]",
            "Signal to Coincidence Ratio",
            "Signal/Coincidence Ratio vs Beam Current",
            "22_signal_to_coin_ratio_vs_current.png",
            output_dir,
            pdf=pdf,
            ylim=(0, 1),
        )

        print("[INFO] Creating beam characteristics plots...")
        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["current_mean_uA"].values,
            "Run Number",
            "Mean Beam Current [uA]",
            "Beam Current vs Run",
            "23_beam_current_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["accumulated_charge(mC)"].values,
            "Run Number",
            "Accumulated Charge [mC]",
            "Accumulated Charge vs Run",
            "24_accumulated_charge_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["Beam_Time(s)"].values,
            "Run Number",
            "Beam Time [s]",
            "Beam Time vs Run",
            "25_beam_time_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        print("[INFO] Creating CPU livetime plot...")
        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["run"].values,
            df["NewGen_EDTM_livetime"].values,
            "Run Number",
            "CPU Live Time Factor (NewGen_EDTM_livetime)",
            "CPU Live Time Factor vs Run",
            "26_cpu_lt_factor_vs_run.png",
            output_dir,
            pdf=pdf,
        )

        print("[INFO] Creating multi-panel comparison plots...")
        rate_data = {
            "Raw Rate": (df["run"].values, df["raw_rate"].values),
            "LCF Corrected Rate": (df["run"].values, df["lcf_corrected_rate"].values),
            "Current Normalized Rate": (
                df["run"].values,
                df["current_normalized_rate"].values,
            ),
            "Scaled Current Norm. Rate": (
                df["run"].values,
                df["scaled_current_normalized_rate"].values,
            ),
        }
        plot_count += multi_scatter_plot(
            rate_data,
            "Run Number",
            "Rate [various units]",
            "Rate Metrics Comparison - Individual Panels",
            "27_rate_metrics_multipanel.png",
            output_dir,
            figsize=(14, 10),
            pdf=pdf,
        )

        pi0_data = {
            "Peak Position (mu)": (df["run"].values, df["pi0_mu_MeV"].values),
            "Peak Width (sigma)": (df["run"].values, df["pi0_sigma_MeV"].values),
            "Signal Counts": (df["run"].values, df["pi0_signal_counts"].values),
            "chi2 / NDF": (df["run"].values, df["chi2_ndf_comb_bg"].values),
        }
        plot_count += multi_scatter_plot(
            pi0_data,
            "Run Number",
            "Value",
            "pi0 Properties - Individual Panels",
            "28_pi0_properties_multipanel.png",
            output_dir,
            figsize=(14, 10),
            pdf=pdf,
        )

        print("[INFO] Creating target boiling correction vs current plot...")
        fig, ax = setup_publication_plot()
        plot_count += scatter_plot(
            ax,
            df["current_mean_uA"].values,
            target_boiling_correction.values,
            "Beam Current [uA]",
            "Target Boiling Correction",
            "Target Boiling Correction vs Beam Current",
            "target_boiling_correction_vs_current.png",
            output_dir,
            pdf=pdf,
            ylim=(0.92, 1.0),
        )

    output_csv = output_dir / "summary_with_rates_and_prescale.csv"
    df.to_csv(output_csv, index=False)
    print(f"[INFO] Augmented CSV saved to: {output_csv}")

    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    summary_metrics = [
        ("Raw Rate (Hz)", "raw_rate", ".3e"),
        ("LCF Corrected Rate (Hz)", "lcf_corrected_rate", ".3e"),
        ("Current Normalized Rate (Hz/uA)", "current_normalized_rate", ".3e"),
        ("Normalized Yield (prescale corrected)", "normalized_yield", ".3e"),
        ("pi0 Peak Position mu (MeV)", "pi0_mu_MeV", ".3f"),
        ("pi0 Peak Width sigma (MeV)", "pi0_sigma_MeV", ".3f"),
    ]
    for label, col, fmt in summary_metrics:
        print(f"\n{label}:")
        print(f"  Mean:   {df[col].mean():{fmt}}")
        print(f"  Min:    {df[col].min():{fmt}}")
        print(f"  Max:    {df[col].max():{fmt}}")
        print(f"  Std:    {df[col].std():{fmt}}")

    print(f"\n[SAVE] PDF with all plots saved to: {pdf_path}")
    print("\n" + "=" * 70)
    print(f"All plots saved to: {output_dir}")
    print(f"   Total plots generated: {plot_count}")
    print(f"   PDF compiled: {pdf_path}")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
