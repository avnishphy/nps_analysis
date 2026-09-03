#!/usr/bin/env python3
from __future__ import annotations

"""Publication and post-event-selection plots for the NPS pi0 analysis.

Default mode compares the combined data ``physics`` tree with the three raw
Geant4/SIMC ``nerd`` trees.  It applies the basic simulation event selection,
uses charge-weighted combined-data scaling, reads SIMC ``normfac`` and ``Ngen``
from the matching .hist files, and produces absolute yield-per-mC overlays.

The previous summary/efficiency run-QA workflow remains available through
``--legacy-run-qa``.
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
FIGURE_SIZE_SINGLE = (8.6, 5.4)
FIGURE_SIZE_COMPARISON = (12.0, 8.0)
FONT_SIZE_TITLE = 14
FONT_SIZE_LABEL = 11.5
FONT_SIZE_LEGEND = 9.5
FONT_SIZE_TICK = 9.5
COLOR_MAIN = "#0072B2"
EDGE_COLOR = "#1a1a1a"
EDGE_WIDTH = 0.65
MARKER_SIZE = 48
ALPHA_SCATTER = 0.92

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 0.9,
        "axes.titleweight": "semibold",
        "axes.labelweight": "normal",
        "axes.formatter.use_mathtext": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "legend.frameon": False,
        "savefig.facecolor": "white",
    }
)

REQUIRED_EFFICIENCY_COLUMNS = [
    "run_number",
    "kinematic_setting",
    "HEL_charge_after_cut_uC",
    "HMS_tracking_eff",
    "HMS_hodo_3of4_eff",
    "NewGen_EDTM_livetime",
    "HMS_tracking_eff_err",
    "HMS_hodo_3of4_eff_err",
    "NewGen_EDTM_livetime_err",
]


def parse_args(argv=None) -> argparse.Namespace:
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
    parser.add_argument(
        "--ignore-runs",
        default="",
        help="Comma-separated runs/ranges to omit without an interactive prompt.",
    )
    parser.add_argument(
        "--include-runs",
        default="",
        help="Optional comma-separated runs/ranges to retain (used by post-analysis to match the combined ROOT file).",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional output directory override.",
    )
    return parser.parse_args(argv)


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
    ax.grid(True, axis="y", linestyle="--", color="#b8b8b8", alpha=0.42, linewidth=0.65)
    ax.set_axisbelow(True)
    return fig, ax


def format_publication_plot(ax, xlabel, ylabel, title, tight=True):
    ax.set_xlabel(xlabel, fontsize=FONT_SIZE_LABEL, labelpad=7)
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL, labelpad=8)
    ax.set_title(title, fontsize=FONT_SIZE_TITLE, pad=12)
    ax.minorticks_on()
    ax.tick_params(
        axis="both",
        which="both",
        top=True,
        right=True,
        labelsize=FONT_SIZE_TICK,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
        spine.set_color("#333333")
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
    yerr=None,
    highlight_mask=None,
    highlight_annotations=None,
):
    mask = finite_mask(x, y)
    x_clean = x[mask]
    y_clean = y[mask]

    if len(x_clean) == 0:
        print(f"[WARN] No valid data for plot: {filename}")
        plt.close(ax.figure)
        return 0

    flagged = (
        np.asarray(highlight_mask, dtype=bool)[mask]
        if highlight_mask is not None
        else np.zeros(len(x_clean), dtype=bool)
    )
    ax.scatter(
        x_clean[~flagged],
        y_clean[~flagged],
        s=MARKER_SIZE,
        color=COLOR_MAIN,
        edgecolor=EDGE_COLOR,
        linewidth=EDGE_WIDTH,
        alpha=ALPHA_SCATTER,
        zorder=3,
        label="Nominal" if np.any(flagged) else None,
    )
    if np.any(flagged):
        ax.scatter(
            x_clean[flagged],
            y_clean[flagged],
            s=MARKER_SIZE * 1.25,
            marker="X",
            color="#D55E00",
            edgecolor=EDGE_COLOR,
            linewidth=EDGE_WIDTH,
            zorder=4,
            label="QA flag",
        )
        if highlight_annotations is not None:
            annotations = np.asarray(highlight_annotations)[mask]
            for x_value, y_value, annotation in zip(
                x_clean[flagged], y_clean[flagged], annotations[flagged]
            ):
                ax.annotate(
                    str(annotation),
                    (x_value, y_value),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=FONT_SIZE_TICK - 1,
                    color="#A33B00",
                )
    if yerr is not None:
        err_clean = np.asarray(yerr, dtype=float)[mask]
        err_mask = np.isfinite(err_clean) & (err_clean >= 0)
        if np.any(err_mask):
            ax.errorbar(
                x_clean[err_mask],
                y_clean[err_mask],
                yerr=err_clean[err_mask],
                fmt="none",
                ecolor="#303030",
                elinewidth=1.5,
                capsize=4.0,
                capthick=1.5,
                alpha=0.95,
                zorder=2,
            )
    format_publication_plot(ax, xlabel, ylabel, title)
    if np.any(flagged):
        ax.legend(fontsize=FONT_SIZE_LEGEND, loc="best")
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
        ax.grid(True, axis="y", linestyle="--", color="#b8b8b8", alpha=0.42, linewidth=0.65)
        ax.set_axisbelow(True)
        ax.set_xlabel(xlabel, fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL - 1)
        ax.set_title(label, fontsize=FONT_SIZE_LABEL)
        ax.minorticks_on()
        ax.tick_params(
            axis="both",
            which="both",
            direction="in",
            top=True,
            right=True,
            labelsize=FONT_SIZE_TICK - 1,
        )
        for spine in ax.spines.values():
            spine.set_linewidth(0.9)
            spine.set_color("#333333")

    for idx in range(len(data_dict), len(axes)):
        axes[idx].set_visible(False)

    if not any_data:
        print(f"[WARN] No valid data for plot: {filename}")
        plt.close(fig)
        return 0

    fig.suptitle(title, fontsize=FONT_SIZE_TITLE, weight="semibold", y=0.995)
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
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
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
                linewidth=1.5,
                markersize=5,
                label=label,
                color=colors[idx % len(colors)],
                alpha=0.9,
                zorder=2,
            )

    if not any_data:
        print(f"[WARN] No valid data for plot: {filename}")
        plt.close(fig)
        return 0

    ax.grid(True, axis="y", linestyle="--", color="#b8b8b8", alpha=0.42, linewidth=0.65)
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


def legacy_main(argv=None, prompt_for_runs=True) -> int:
    args = parse_args(argv)
    summary_csv = Path(args.summary_csv).expanduser().resolve()
    efficiency_csv = Path(args.efficiency_csv).expanduser().resolve()
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else plots_dir_from_summary(summary_csv)
    )

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
    efficiency_columns = [*REQUIRED_EFFICIENCY_COLUMNS, "ps_factor"]
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
        "HMS_hodo_3of4_eff",
        "NewGen_EDTM_livetime",
        "HMS_tracking_eff_err",
        "HMS_hodo_3of4_eff_err",
        "NewGen_EDTM_livetime_err",
        "ps_factor",
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

    included_runs = parse_run_ignore_text(args.include_runs) if args.include_runs else set()
    if included_runs:
        before = len(df)
        df = df[df["run"].isin(included_runs)].copy()
        missing_in_csv = included_runs - set(df["run"])
        print(f"[INFO] Restricted stability sample to combined-ROOT runs: {before} -> {len(df)} rows")
        if missing_in_csv:
            print(f"[WARN] {len(missing_in_csv)} combined-ROOT run(s) missing from stability CSV inputs")

    valid_efficiency = (
        (df["HEL_charge_after_cut_uC"] > 0)
        & (df["HMS_tracking_eff"] > 0)
        & (df["HMS_hodo_3of4_eff"] > 0)
        & (df["NewGen_EDTM_livetime"] > 0)
        & (df["HMS_tracking_eff_err"] >= 0)
        & (df["HMS_hodo_3of4_eff_err"] >= 0)
        & (df["NewGen_EDTM_livetime_err"] >= 0)
        & (df["ps_factor"] > 0)
    )
    invalid_efficiency_runs = set(df.loc[~valid_efficiency, "run"].unique())
    if invalid_efficiency_runs:
        df = df[valid_efficiency].copy()
        print(
            f"[INFO] Removed {len(invalid_efficiency_runs)} run(s) with invalid "
            "charge, efficiency, livetime, uncertainty, or prescale values"
        )

    suspicious_runs = set()
    for col in df.columns:
        if (
            col.startswith("scaler_") or col.startswith("trig_")
        ) and pd.api.types.is_numeric_dtype(df[col]):
            suspicious_runs.update(df[df[col] > 10_000_000]["run"].unique())

    if suspicious_runs:
        df = df[~df["run"].isin(suspicious_runs)].copy()
        print(f"[INFO] Removed {len(suspicious_runs)} suspicious runs")

    ignored_runs = (
        parse_run_ignore_text(args.ignore_runs)
        if args.ignore_runs
        else (prompt_for_ignored_runs() if prompt_for_runs else set())
    )
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
    livetime_err = df["NewGen_EDTM_livetime_err"]
    tracking_eff = df["HMS_tracking_eff"]
    tracking_eff_err = df["HMS_tracking_eff_err"]
    hodo_eff = df["HMS_hodo_3of4_eff"]
    hodo_eff_err = df["HMS_hodo_3of4_eff_err"]
    hms_eff = tracking_eff * hodo_eff
    hms_eff_err = np.hypot(hodo_eff * tracking_eff_err, tracking_eff * hodo_eff_err)
    target_boiling_correction = 1 - ((10.2 / 100) * df["current_mean_uA"]) / 100

    if "pi0_signal_counts_err" in df.columns:
        signal_counts_err = pd.to_numeric(df["pi0_signal_counts_err"], errors="coerce")
        invalid_signal_err = ~np.isfinite(signal_counts_err) | (signal_counts_err < 0)
        signal_counts_err = signal_counts_err.mask(
            invalid_signal_err,
            np.sqrt(np.abs(df["pi0_signal_counts"])),
        )
        signal_error_source = "pi0_signal_counts_err (sqrt(|N|) fallback for invalid rows)"
    else:
        signal_counts_err = np.sqrt(np.abs(df["pi0_signal_counts"]))
        signal_error_source = "sqrt(|pi0_signal_counts|) statistical floor"
        print(
            "[WARN] pi0_signal_counts_err is unavailable; normalized-yield errors "
            "use sqrt(|N_pi0|) as a statistical floor."
        )
    df["pi0_signal_counts_err_used"] = signal_counts_err
    df["pi0_signal_counts_err_source"] = signal_error_source

    df["efficiency"] = hms_eff
    df["efficiency_err"] = hms_eff_err
    df["livetime_err"] = livetime_err

    df["raw_rate"] = np.where(
        df["Beam_Time(s)"] > 0,
        df["total_coin_entries"] / df["Beam_Time(s)"],
        np.nan,
    )
    df["lcf_corrected_rate"] = np.where(
        (df["Beam_Time(s)"] > 0) & (livetime > 0),
        df["total_coin_entries"] / (df["Beam_Time(s)"] * livetime),
        np.nan,
    )
    df["lcf_corrected_rate_err"] = np.abs(df["lcf_corrected_rate"]) * (
        livetime_err / livetime
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
            / target_boiling_correction
        )
        * df["ps_factor"],
        np.nan,
    )
    signal_relative_err = np.divide(
        signal_counts_err,
        np.abs(df["pi0_signal_counts"]),
        out=np.full(len(df), np.nan, dtype=float),
        where=np.abs(df["pi0_signal_counts"]) > 0,
    )
    normalization_relative_err = np.sqrt(
        signal_relative_err**2
        + (livetime_err / livetime) ** 2
        + (hms_eff_err / hms_eff) ** 2
    )
    df["normalized_yield_err"] = np.abs(df["normalized_yield"]) * normalization_relative_err

    qa_reasons = pd.Series("", index=df.index, dtype=object)
    livetime_flag = (livetime < 0.8) | (livetime > 1.05)
    qa_reasons.loc[livetime_flag] += "livetime outside [0.8,1.05];"
    current_group = df["current_mean_uA"].round(1)
    for _, group in df.groupby([current_group, "ps_factor"]):
        if len(group) < 5:
            continue
        median_yield = group["normalized_yield"].median()
        mad_yield = (group["normalized_yield"] - median_yield).abs().median()
        outlier_threshold = max(5.0 * mad_yield, 0.20 * abs(median_yield))
        yield_flag = (group["normalized_yield"] - median_yield).abs() > outlier_threshold
        qa_reasons.loc[group.index[yield_flag]] += "yield outlier within current/prescale group;"
    if "which_TRIG" in df.columns:
        valid_triggers = df["which_TRIG"].dropna().astype(str)
        if not valid_triggers.empty:
            nominal_trigger = valid_triggers.mode().iloc[0]
            trigger_flag = df["which_TRIG"].astype(str) != nominal_trigger
            qa_reasons.loc[trigger_flag] += f"trigger differs from {nominal_trigger};"
    df["normalized_yield_qa_reason"] = qa_reasons.str.rstrip(";")
    df["normalized_yield_qa_flag"] = df["normalized_yield_qa_reason"].ne("")
    if df["normalized_yield_qa_flag"].any():
        flagged_text = ", ".join(
            f"{row.run} ({row.normalized_yield_qa_reason})"
            for row in df.loc[df["normalized_yield_qa_flag"], ["run", "normalized_yield_qa_reason"]].itertuples()
        )
        print(f"[WARN] Normalized-yield QA flags: {flagged_text}")
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
            yerr=df["lcf_corrected_rate_err"].values,
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
            yerr=df["lcf_corrected_rate_err"].values,
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
            ylim=(-0.05, 1.05),
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
            yerr=df["normalized_yield_err"].values,
            highlight_mask=df["normalized_yield_qa_flag"].values,
            highlight_annotations=df["run"].values,
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
            yerr=df["normalized_yield_err"].values,
            highlight_mask=df["normalized_yield_qa_flag"].values,
            highlight_annotations=df["run"].values,
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
            yerr=df["NewGen_EDTM_livetime_err"].values,
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

        print("[INFO] Creating stability trends versus corrected rate...")
        rate_stability_plots = (
            ("normalized_yield", "Normalized Yield", "29_normalized_yield_vs_lcf_corrected_rate.png", df["normalized_yield_err"].values),
            ("pi0_mu_MeV", r"pi0 Peak Position $\mu$ [MeV]", "30_pi0_peak_position_vs_lcf_corrected_rate.png", None),
            ("pi0_sigma_MeV", r"pi0 Peak Width $\sigma$ [MeV]", "31_pi0_peak_width_vs_lcf_corrected_rate.png", None),
            ("signal_to_coin_ratio", "Signal/Coincidence Ratio", "32_signal_to_coin_ratio_vs_lcf_corrected_rate.png", None),
            ("detection_efficiency", "Detection Efficiency", "33_detection_efficiency_vs_lcf_corrected_rate.png", None),
        )
        for column, ylabel, filename, yerr in rate_stability_plots:
            fig, ax = setup_publication_plot()
            plot_count += scatter_plot(
                ax,
                df["lcf_corrected_rate"].values,
                df[column].values,
                "LCF Corrected Rate [Hz]",
                ylabel,
                f"{ylabel} vs LCF Corrected Rate",
                filename,
                output_dir,
                pdf=pdf,
                yerr=yerr,
                highlight_mask=(
                    df["normalized_yield_qa_flag"].values
                    if column == "normalized_yield"
                    else None
                ),
                highlight_annotations=(
                    df["run"].values if column == "normalized_yield" else None
                ),
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


# ---------------------------------------------------------------------------
# Combined-data / raw-Geant4 post-event-selection plotting mode
# ---------------------------------------------------------------------------

from dataclasses import dataclass
import json
import re

import uproot


POST_DEFAULT_SIM_BASE = Path("/volatile/hallc/nps/singhav/geant4_simc")
POST_DEFAULT_HIST_DIR = Path(
    "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/"
    "root_analysis_env_main/output/simc/nps_simc_20260824_135058/"
    "outfiles/simc_gfortran_updated/outfiles"
)
POST_DEFAULT_KIN_CONFIG = Path(__file__).resolve().parents[1] / "config/nps_simulation_kinematics.csv"
POST_CHANNELS = ("excl", "sidis", "delta")
POST_CHANNEL_LABELS = {"excl": "Exclusive", "sidis": "SIDIS", "delta": r"$\Delta\pi^0$"}
POST_CHANNEL_COLORS = {"excl": "#009E73", "sidis": "#E69F00", "delta": "#CC79A7"}
POST_DATA_COLOR = "#111111"
POST_SIM_COLOR = "#0072B2"
POST_ELECTRON_MASS_GEV = 0.0005109989461
POST_PROTON_MASS_GEV = 0.9382720813
POST_PI0_MASS_GEV = 0.1349768


@dataclass(frozen=True)
class PostPlotSpec:
    key: str
    title: str
    xlabel: str
    bins: int
    limits: tuple[float, float]
    group: str


@dataclass
class PostHist:
    counts: np.ndarray
    sumw2: np.ndarray
    entries: int = 0
    weight_sum: float = 0.0


POST_SPECS = (
    PostPlotSpec("mpi0", r"Two-photon invariant mass", r"$m_{\gamma\gamma}$ [GeV/$c^2$]", 100, (0.04, 0.30), "mass"),
    # Deliberately broad and never subject to an exclusivity flag.
    PostPlotSpec("mmiss", r"Missing mass (pre-exclusive selection)", r"$M_X(ep\to e'\pi^0X)$ [GeV/$c^2$]", 120, (0.0, 3.0), "mass"),
    PostPlotSpec("Q2", r"Four-momentum transfer", r"$Q^2$ [GeV$^2$]", 90, (0.0, 10.0), "physics"),
    PostPlotSpec("W", r"Hadronic invariant mass", r"$W$ [GeV/$c^2$]", 90, (0.5, 4.5), "physics"),
    PostPlotSpec("nu", r"Energy transfer", r"$\nu$ [GeV]", 90, (0.0, 10.0), "physics"),
    PostPlotSpec("xB", r"Bjorken $x$", r"$x_B$", 90, (0.0, 1.0), "physics"),
    PostPlotSpec("t", r"Mandelstam $t$", r"$t$ [GeV$^2$]", 90, (-5.0, 0.5), "physics"),
    PostPlotSpec("tmin", r"Forward-limit momentum transfer", r"$t_{\min}$ [GeV$^2$]", 90, (-5.0, 0.1), "physics"),
    PostPlotSpec("pt", r"Pion transverse momentum", r"$p_T$ [GeV/$c$]", 90, (0.0, 2.0), "physics"),
    PostPlotSpec("theta", r"Pion angle relative to $q$", r"$\theta_{\pi q}$ [rad]", 90, (0.0, 0.8), "physics"),
    PostPlotSpec("phi", r"Hadron-plane azimuth", r"$\phi$ [rad]", 90, (-np.pi, np.pi), "physics"),
    PostPlotSpec("z", r"Pion energy fraction", r"$z=E_{\pi}/\nu$", 90, (0.0, 1.5), "physics"),
    PostPlotSpec("delta", "HMS momentum acceptance", r"HMS $\delta$ [%]", 80, (-15.0, 15.0), "hms"),
    PostPlotSpec("xptar", "HMS target-plane slope", r"HMS $x'_{tar}$ [rad]", 80, (-0.11, 0.11), "hms"),
    PostPlotSpec("yptar", "HMS target-plane slope", r"HMS $y'_{tar}$ [rad]", 80, (-0.05, 0.05), "hms"),
    PostPlotSpec("ytar", "HMS target-plane position", r"HMS $y_{tar}$ [cm]", 80, (-10.0, 10.0), "hms"),
    PostPlotSpec("xfp", "HMS focal-plane position", r"HMS $x_{fp}$ [cm]", 80, (-60.0, 60.0), "hms"),
    PostPlotSpec("yfp", "HMS focal-plane position", r"HMS $y_{fp}$ [cm]", 80, (-40.0, 40.0), "hms"),
    PostPlotSpec("xpfp", "HMS focal-plane slope", r"HMS $x'_{fp}$ [rad]", 80, (-0.12, 0.12), "hms"),
    PostPlotSpec("ypfp", "HMS focal-plane slope", r"HMS $y'_{fp}$ [rad]", 80, (-0.08, 0.08), "hms"),
    PostPlotSpec("cluster_e_1", "Leading photon-cluster energy", r"$E_{\gamma 1}$ [GeV]", 90, (0.8, 10.0), "nps"),
    PostPlotSpec("cluster_e_2", "Second photon-cluster energy", r"$E_{\gamma 2}$ [GeV]", 90, (0.8, 10.0), "nps"),
    PostPlotSpec("cluster_e_sum", "Two-cluster energy sum", r"$E_{\gamma 1}+E_{\gamma 2}$ [GeV]", 90, (1.6, 12.0), "nps"),
    PostPlotSpec("cluster_e_asym", "Two-cluster energy asymmetry", r"$|E_1-E_2|/(E_1+E_2)$", 80, (0.0, 1.0), "nps"),
    PostPlotSpec("cluster_x_1", "Leading photon-cluster position", r"$x_{\gamma 1}$ [cm]", 80, (-30.0, 30.0), "nps"),
    PostPlotSpec("cluster_y_1", "Leading photon-cluster position", r"$y_{\gamma 1}$ [cm]", 80, (-36.0, 36.0), "nps"),
    PostPlotSpec("cluster_x_2", "Second photon-cluster position", r"$x_{\gamma 2}$ [cm]", 80, (-30.0, 30.0), "nps"),
    PostPlotSpec("cluster_y_2", "Second photon-cluster position", r"$y_{\gamma 2}$ [cm]", 80, (-36.0, 36.0), "nps"),
    PostPlotSpec("cluster_sep", "Photon-cluster separation", r"$\Delta r_{\gamma\gamma}$ [cm]", 90, (0.0, 90.0), "nps"),
    PostPlotSpec("nclust", "Selected NPS cluster multiplicity", "Number of clusters", 12, (0.5, 12.5), "nps"),
)


def post_build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Post-selection combined-data versus Geant4/SIMC plots."
    )
    parser.add_argument("--combined-root", required=True, help="Combined data ROOT file from combine_analysis_branches.py.")
    parser.add_argument("--kin", default=None, help="Kin_old name, e.g. KinC_x25_3; inferred from combined path when possible.")
    parser.add_argument("--simulation-base", default=str(POST_DEFAULT_SIM_BASE), help="Directory containing nps_geant4_* production directories.")
    parser.add_argument("--simulation-production", default=None, help="Explicit nps_geant4_* production directory; disables newest-complete discovery.")
    parser.add_argument("--hist-dir", default=str(POST_DEFAULT_HIST_DIR), help="Directory containing matching SIMC .hist files.")
    parser.add_argument("--kin-config", default=str(POST_DEFAULT_KIN_CONFIG), help="Simulation kinematics CSV.")
    parser.add_argument("--output-dir", default=None, help="Output directory; default: <combined kin>/plots/post_analysis.")
    parser.add_argument("--summary-csv", default=None, help="Run-summary CSV; default: <combined kin>/summary/summary_all_runs.csv.")
    parser.add_argument("--efficiency-csv", default=None, help="Efficiency CSV; default: repository output/efficiency_stuff/efficiency_<kin>.csv.")
    parser.add_argument("--run-status", default=None, help="Optional exact run_status retained in stability plots.")
    parser.add_argument("--ignore-runs", default="", help="Comma-separated runs/ranges omitted from stability plots.")
    parser.add_argument("--skip-run-stability", action="store_true", help="Skip run/current/rate stability and normalized-yield plots.")
    parser.add_argument("--normalization", choices=("absolute", "shape"), default="absolute", help="absolute: physical yield per 1 mC; shape: per-plot SIMC area match after physical channel weighting.")
    parser.add_argument("--step-size", default="100 MB", help="uproot iteration step size, e.g. '100 MB' or '100000 entries'.")
    parser.add_argument("--max-events", type=int, default=None, help="Debug-only per-file entry limit; absolute yield is incomplete when used.")
    parser.add_argument("--formats", default="png,pdf", help="Comma-separated individual figure formats (png,pdf).")
    return parser


def post_infer_kin(path: Path) -> str:
    match = re.search(r"KinC_x[0-9A-Za-z_]+", str(path))
    if not match:
        raise ValueError("Could not infer --kin from --combined-root; pass --kin explicitly.")
    return match.group(0).rstrip("_")


def post_kin_token(kin: str) -> str:
    match = re.fullmatch(r"KinC_(x[0-9A-Za-z_]+)", kin.strip())
    if not match:
        raise ValueError(f"Invalid Kin_old value: {kin}")
    return match.group(1)


def post_fortran_float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def post_read_hist_meta(path: Path) -> dict[str, float]:
    text = path.read_text(errors="replace")
    patterns = {
        "normfac": r"(?mi)^\s*normfac\s*=\s*([+\-0-9.EeDd]+)",
        "ngen": r"(?mi)^\s*Ngen\s*\(request\)\s*=\s*([+\-0-9.EeDd]+)",
        "ebeam_mev": r"(?mi)^\s*Ebeam\s*=\s*([+\-0-9.EeDd]+)",
        "hms_y_offset_cm": r"(?mi)^\s*y offset\s*=\s*([+\-0-9.EeDd]+)",
    }
    out: dict[str, float] = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            out[key] = post_fortran_float(match.group(1))
    missing = [key for key in ("normfac", "ngen", "ebeam_mev", "hms_y_offset_cm") if key not in out]
    if missing:
        raise ValueError(f"{path} lacks SIMC metadata: {', '.join(missing)}")
    if out["normfac"] <= 0 or out["ngen"] <= 0:
        raise ValueError(f"{path} has non-positive normfac or Ngen")
    return out


def post_load_kinematics(path: Path, kin: str) -> dict[str, float]:
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [str(col).strip() for col in df.columns]
    if "kin_old" not in df:
        raise ValueError(f"{path} lacks kin_old column")
    rows = df[df["kin_old"].astype(str).str.strip() == kin]
    if len(rows) != 1:
        raise ValueError(f"Expected one {kin} row in {path}; found {len(rows)}")
    row = rows.iloc[0]
    required = ("ebeam_gev", "nps_theta_deg", "nps_target_distance_cm")
    values = {key: float(row[key]) for key in required}
    if not all(np.isfinite(value) for value in values.values()):
        raise ValueError(f"Non-finite kinematics for {kin} in {path}")
    return values


def post_find_simulation_files(base: Path, explicit: str | None, token: str) -> tuple[Path, dict[str, Path]]:
    roots = [Path(explicit).expanduser().resolve()] if explicit else []
    if not roots:
        base = base.expanduser().resolve()
        roots = [base] + sorted(base.glob("nps_geant4_*"), key=lambda p: p.name, reverse=True)
    names = {
        "excl": f"nps_excl_pi0_{token}_geant4.root",
        "sidis": f"nps_sidis_pi0_{token}_geant4.root",
        "delta": f"nps_delta_pi0_{token}_geant4.root",
    }
    for root in roots:
        files = {channel: root / channel / name for channel, name in names.items()}
        if all(path.is_file() for path in files.values()):
            return root, files
    raise FileNotFoundError(f"No complete Geant4 production found for {token} under {base}")


def post_empty_histograms() -> dict[str, PostHist]:
    return {
        spec.key: PostHist(np.zeros(spec.bins), np.zeros(spec.bins))
        for spec in POST_SPECS
    }


def post_fill_histograms(histograms: dict[str, PostHist], values: dict[str, np.ndarray], weights: np.ndarray) -> None:
    for spec in POST_SPECS:
        value = np.asarray(values[spec.key], dtype=float)
        mask = np.isfinite(value) & np.isfinite(weights)
        if not np.any(mask):
            continue
        selected_value = value[mask]
        selected_weight = weights[mask]
        counts, _ = np.histogram(selected_value, bins=spec.bins, range=spec.limits, weights=selected_weight)
        sumw2, _ = np.histogram(selected_value, bins=spec.bins, range=spec.limits, weights=selected_weight**2)
        hist = histograms[spec.key]
        hist.counts += counts
        hist.sumw2 += sumw2
        hist.entries += len(selected_value)
        hist.weight_sum += float(np.sum(selected_weight))


POST_DATA_BRANCHES = (
    "mpi0_all", "mmiss_all", "pi0_weight", "Q2", "W", "t", "tmin", "pt", "theta", "phi", "xB", "z",
    "delta", "xptar", "yptar", "ytar", "xfp", "yfp", "xpfp", "ypfp", "nclust_selected",
    "cluster_e_1", "cluster_e_2", "cluster_x_1", "cluster_y_1", "cluster_x_2", "cluster_y_2",
    "scale", "charge_uC", "run_number",
)


def post_data_values(arrays: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    e1 = np.asarray(arrays["cluster_e_1"], dtype=float)
    e2 = np.asarray(arrays["cluster_e_2"], dtype=float)
    esum = e1 + e2
    q2 = np.asarray(arrays["Q2"], dtype=float)
    xb = np.asarray(arrays["xB"], dtype=float)
    nu = np.divide(q2, 2.0 * POST_PROTON_MASS_GEV * xb, out=np.full_like(q2, np.nan), where=xb != 0)
    dx = np.asarray(arrays["cluster_x_1"], dtype=float) - np.asarray(arrays["cluster_x_2"], dtype=float)
    dy = np.asarray(arrays["cluster_y_1"], dtype=float) - np.asarray(arrays["cluster_y_2"], dtype=float)
    values = {
        "mpi0": arrays["mpi0_all"], "mmiss": arrays["mmiss_all"], "Q2": q2, "W": arrays["W"],
        "nu": nu, "xB": xb, "t": arrays["t"], "tmin": arrays["tmin"], "pt": arrays["pt"],
        "theta": arrays["theta"], "phi": arrays["phi"], "z": arrays["z"], "delta": arrays["delta"],
        "xptar": arrays["xptar"], "yptar": arrays["yptar"], "ytar": arrays["ytar"], "xfp": arrays["xfp"],
        "yfp": arrays["yfp"], "xpfp": arrays["xpfp"], "ypfp": arrays["ypfp"],
        "cluster_e_1": e1, "cluster_e_2": e2, "cluster_e_sum": esum,
        "cluster_e_asym": np.divide(np.abs(e1 - e2), esum, out=np.full_like(esum, np.nan), where=esum != 0),
        "cluster_x_1": arrays["cluster_x_1"], "cluster_y_1": arrays["cluster_y_1"],
        "cluster_x_2": arrays["cluster_x_2"], "cluster_y_2": arrays["cluster_y_2"],
        "cluster_sep": np.hypot(dx, dy), "nclust": arrays["nclust_selected"],
    }
    return values


def post_total_data_charge(tree, step_size: str, entry_stop: int | None) -> tuple[float, dict[int, float]]:
    run_charge: dict[int, float] = {}
    for arrays in tree.iterate(("run_number", "charge_uC"), step_size=step_size, entry_stop=entry_stop, library="np", how=dict):
        runs = np.asarray(arrays["run_number"], dtype=int)
        charges = np.asarray(arrays["charge_uC"], dtype=float)
        for run in np.unique(runs):
            values = charges[runs == run]
            finite = values[np.isfinite(values) & (values > 0)]
            if finite.size == 0:
                raise ValueError(f"Run {run} has no positive charge_uC")
            charge = float(finite[0])
            if not np.allclose(finite, charge, rtol=2e-6, atol=1e-4):
                raise ValueError(f"Run {run} has inconsistent charge_uC values")
            run_charge[int(run)] = charge
    total = float(sum(run_charge.values()))
    if total <= 0:
        raise ValueError("Combined data has non-positive total unique-run charge")
    return total, run_charge


def post_read_data(path: Path, step_size: str, entry_stop: int | None) -> tuple[dict[str, PostHist], dict[str, object]]:
    histograms = post_empty_histograms()
    with uproot.open(path, handler=uproot.source.file.MultithreadedFileSource) as root_file:
        if "physics" not in root_file:
            raise KeyError(f"{path} lacks physics tree")
        tree = root_file["physics"]
        missing = [name for name in POST_DATA_BRANCHES if name not in tree.keys()]
        if missing:
            raise ValueError(f"{path} lacks required combined branches: {', '.join(missing)}")
        total_charge_uC, run_charge = post_total_data_charge(tree, step_size, entry_stop)
        entries = 0
        selected_weight_sum = 0.0
        for arrays in tree.iterate(POST_DATA_BRANCHES, step_size=step_size, entry_stop=entry_stop, library="np", how=dict):
            charge_fraction = np.asarray(arrays["charge_uC"], dtype=float) / total_charge_uC
            weights = (
                np.asarray(arrays["scale"], dtype=float)
                * np.asarray(arrays["pi0_weight"], dtype=float)
                * charge_fraction
            )
            post_fill_histograms(histograms, post_data_values(arrays), weights)
            entries += len(weights)
            selected_weight_sum += float(np.sum(weights[np.isfinite(weights)]))
    audit = {
        "source": "data", "entries_read": entries, "entries_selected": entries,
        "total_charge_uC": total_charge_uC, "runs": len(run_charge),
        "run_numbers": sorted(run_charge),
        "selected_weight_sum": selected_weight_sum,
        "weight_formula": "scale*pi0_weight*(charge_uC/total_unique_run_charge_uC)",
    }
    return histograms, audit


POST_SIM_BRANCHES = (
    "Weight", "phot1_hit", "phot2_hit", "phot1_Ecal", "phot2_Ecal",
    "phot1_hit_x", "phot1_hit_y", "phot2_hit_x", "phot2_hit_y", "phot1_vx",
    "hsdelta", "hsxptar", "hsyptar", "hsytar", "hsxfp", "hsyfp", "hsxpfp", "hsypfp",
    "sc_e_Px", "sc_e_Py", "sc_e_Pz", "nClusters",
)


def post_sim_selection(arrays: dict[str, np.ndarray], y_mispoint_cm: float) -> tuple[np.ndarray, np.ndarray]:
    theta_c = np.arctan2(arrays["sc_e_Py"], arrays["sc_e_Pz"])
    abs_theta = np.abs(theta_c)
    denom = np.sin(abs_theta) - arrays["hsyptar"] * np.cos(abs_theta)
    numerator = arrays["hsytar"] - y_mispoint_cm - arrays["phot1_vx"] * (
        np.cos(abs_theta) + arrays["hsyptar"] * np.sin(abs_theta)
    )
    hreactz = np.divide(numerator, denom, out=np.full_like(numerator, np.nan, dtype=float), where=np.abs(denom) > 1e-12)
    mask = (
        np.asarray(arrays["phot1_hit"], dtype=bool)
        & np.asarray(arrays["phot2_hit"], dtype=bool)
        & (np.abs(hreactz) <= 8.0)
        & (np.abs(arrays["hsdelta"]) <= 10.0)
        & (np.abs(arrays["hsxptar"]) <= 0.1)
        & (np.abs(arrays["hsyptar"]) <= 0.04)
    )
    for prefix in ("phot1", "phot2"):
        mask &= (
            (arrays[f"{prefix}_Ecal"] > 0.6)
            & (arrays[f"{prefix}_hit_x"] > -26.0) & (arrays[f"{prefix}_hit_x"] < 26.0)
            & (arrays[f"{prefix}_hit_y"] > -32.0) & (arrays[f"{prefix}_hit_y"] < 32.0)
        )
    return mask, hreactz


def post_norm3(vectors: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(0.0, np.einsum("ij,ij->i", vectors, vectors)))


def post_sim_values(arrays: dict[str, np.ndarray], mask: np.ndarray, kin: dict[str, float]) -> dict[str, np.ndarray]:
    def take(name: str) -> np.ndarray:
        return np.asarray(arrays[name], dtype=float)[mask]

    e1, e2 = take("phot1_Ecal"), take("phot2_Ecal")
    x1, y1, x2, y2 = take("phot1_hit_x"), take("phot1_hit_y"), take("phot2_hit_x"), take("phot2_hit_y")
    z_nps = kin["nps_target_distance_cm"]
    theta_nps = np.deg2rad(-kin["nps_theta_deg"])
    ctheta, stheta = np.cos(theta_nps), np.sin(theta_nps)

    def photon(E: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        direction = np.column_stack((x, y, np.full_like(x, z_nps)))
        direction /= post_norm3(direction)[:, None]
        p_nps = E[:, None] * direction
        return np.column_stack((ctheta * p_nps[:, 0] - stheta * p_nps[:, 2], p_nps[:, 1], stheta * p_nps[:, 0] + ctheta * p_nps[:, 2]))

    p1, p2 = photon(e1, x1, y1), photon(e2, x2, y2)
    ppi = p1 + p2
    epi = e1 + e2
    ppi_mag = post_norm3(ppi)
    mpi0 = np.sqrt(np.maximum(0.0, epi**2 - ppi_mag**2))

    px = take("sc_e_Py") / 1000.0
    py = -take("sc_e_Px") / 1000.0
    pz = take("sc_e_Pz") / 1000.0
    pe = np.column_stack((px, py, pz))
    ee = np.sqrt(post_norm3(pe)**2 + POST_ELECTRON_MASS_GEV**2)
    ebeam = kin["ebeam_gev"]
    pbeam = np.sqrt(max(0.0, ebeam**2 - POST_ELECTRON_MASS_GEV**2))
    q0 = ebeam - ee
    q = np.column_stack((-px, -py, pbeam - pz))
    qmag = post_norm3(q)
    q2 = qmag**2 - q0**2
    nu = q0
    xb = np.divide(q2, 2.0 * POST_PROTON_MASS_GEV * nu, out=np.full_like(q2, np.nan), where=nu != 0)
    w2 = POST_PROTON_MASS_GEV**2 + 2.0 * POST_PROTON_MASS_GEV * nu - q2
    w = np.sqrt(np.maximum(0.0, w2))
    z = np.divide(epi, nu, out=np.full_like(epi, np.nan), where=nu != 0)
    q_minus_pi = q - ppi
    t = (q0 - epi)**2 - post_norm3(q_minus_pi)**2
    qhat = np.divide(q, qmag[:, None], out=np.zeros_like(q), where=qmag[:, None] != 0)
    ppar = np.einsum("ij,ij->i", ppi, qhat)
    pt = post_norm3(ppi - ppar[:, None] * qhat)
    cos_theta = np.divide(np.einsum("ij,ij->i", q, ppi), qmag * ppi_mag, out=np.ones_like(qmag), where=(qmag * ppi_mag) != 0)
    pion_theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))

    k = np.column_stack((np.zeros_like(px), np.zeros_like(px), np.full_like(px, pbeam)))
    nlep = np.cross(k, pe)
    nhad = np.cross(q, ppi)
    nlep_norm, nhad_norm = post_norm3(nlep), post_norm3(nhad)
    nlep_hat = np.divide(nlep, nlep_norm[:, None], out=np.zeros_like(nlep), where=nlep_norm[:, None] != 0)
    nhad_hat = np.divide(nhad, nhad_norm[:, None], out=np.zeros_like(nhad), where=nhad_norm[:, None] != 0)
    phi = np.arctan2(np.einsum("ij,ij->i", qhat, np.cross(nlep_hat, nhad_hat)), np.einsum("ij,ij->i", nlep_hat, nhad_hat))

    miss_e = ebeam + POST_PROTON_MASS_GEV - ee - epi
    miss_p = np.column_stack((-px, -py, pbeam - pz)) - ppi
    mmiss = np.sqrt(np.maximum(0.0, miss_e**2 - post_norm3(miss_p)**2))

    threshold = POST_PROTON_MASS_GEV + POST_PI0_MASS_GEV
    tmin = np.zeros_like(w)
    valid_w = w > threshold
    wv, q2v = w[valid_w], q2[valid_w]
    q0_star = (wv**2 - POST_PROTON_MASS_GEV**2 - q2v) / (2.0 * wv)
    q_star = np.sqrt(np.maximum(0.0, q0_star**2 + q2v))
    epi_star = (wv**2 + POST_PI0_MASS_GEV**2 - POST_PROTON_MASS_GEV**2) / (2.0 * wv)
    lam = (wv**2 - (POST_PROTON_MASS_GEV**2 + POST_PI0_MASS_GEV**2))**2 - 4.0 * POST_PROTON_MASS_GEV**2 * POST_PI0_MASS_GEV**2
    ppi_star = np.sqrt(np.maximum(0.0, lam)) / (2.0 * wv)
    tmin[valid_w] = POST_PI0_MASS_GEV**2 - q2v - 2.0 * (q0_star * epi_star - q_star * ppi_star)

    esum = e1 + e2
    return {
        "mpi0": mpi0, "mmiss": mmiss, "Q2": q2, "W": w, "nu": nu, "xB": xb, "t": t,
        "tmin": tmin, "pt": pt, "theta": pion_theta, "phi": phi, "z": z,
        "delta": take("hsdelta"), "xptar": take("hsxptar"), "yptar": take("hsyptar"),
        "ytar": take("hsytar"), "xfp": take("hsxfp"), "yfp": take("hsyfp"),
        "xpfp": take("hsxpfp"), "ypfp": take("hsypfp"),
        "cluster_e_1": e1, "cluster_e_2": e2, "cluster_e_sum": esum,
        "cluster_e_asym": np.divide(np.abs(e1 - e2), esum, out=np.full_like(esum, np.nan), where=esum != 0),
        "cluster_x_1": x1, "cluster_y_1": y1, "cluster_x_2": x2, "cluster_y_2": y2,
        "cluster_sep": np.hypot(x1 - x2, y1 - y2), "nclust": take("nClusters"),
    }


def post_read_simulation(path: Path, meta: dict[str, float], kin: dict[str, float], step_size: str, entry_stop: int | None) -> tuple[dict[str, PostHist], dict[str, object]]:
    histograms = post_empty_histograms()
    entries = 0
    selected = 0
    selected_weight_sum = 0.0
    with uproot.open(path, handler=uproot.source.file.MultithreadedFileSource) as root_file:
        if "nerd" not in root_file:
            raise KeyError(f"{path} lacks nerd tree")
        tree = root_file["nerd"]
        missing = [name for name in POST_SIM_BRANCHES if name not in tree.keys()]
        if missing:
            raise ValueError(f"{path} lacks required simulation branches: {', '.join(missing)}")
        for arrays in tree.iterate(POST_SIM_BRANCHES, step_size=step_size, entry_stop=entry_stop, library="np", how=dict):
            mask, _ = post_sim_selection(arrays, meta["hms_y_offset_cm"])
            weights = meta["normfac"] * np.asarray(arrays["Weight"], dtype=float)[mask] / meta["ngen"]
            post_fill_histograms(histograms, post_sim_values(arrays, mask, kin), weights)
            entries += len(mask)
            selected += int(np.count_nonzero(mask))
            selected_weight_sum += float(np.sum(weights[np.isfinite(weights)]))
        tree_entries = int(tree.num_entries)
    audit = {
        "source": path.name, "tree_entries": tree_entries, "entries_read": entries,
        "entries_selected": selected, "selection_fraction": selected / entries if entries else 0.0,
        "normfac": meta["normfac"], "ngen": meta["ngen"], "hms_y_offset_cm": meta["hms_y_offset_cm"],
        "selected_weight_sum": selected_weight_sum, "weight_formula": "normfac*Weight/Ngen",
    }
    return histograms, audit


def post_style_axis(ax) -> None:
    ax.grid(axis="y", color="#c8c8c8", linestyle=":", linewidth=0.7, alpha=0.7)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(0.9)


def post_plot_comparison(spec: PostPlotSpec, data: PostHist, sim: dict[str, PostHist], kin: str, normalization: str):
    edges = np.linspace(spec.limits[0], spec.limits[1], spec.bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    total_sim = sum((sim[channel].counts for channel in POST_CHANNELS), start=np.zeros(spec.bins))
    total_sim_var = sum((sim[channel].sumw2 for channel in POST_CHANNELS), start=np.zeros(spec.bins))
    sim_factor = 1.0
    if normalization == "shape":
        sim_integral = float(np.sum(total_sim))
        data_integral = float(np.sum(data.counts))
        if sim_integral != 0.0:
            sim_factor = data_integral / sim_integral

    total_sim *= sim_factor
    total_sim_var *= sim_factor**2
    fig, (ax, ratio_ax) = plt.subplots(
        2, 1, figsize=(8.4, 7.0), dpi=FIGURE_DPI, sharex=True,
        gridspec_kw={"height_ratios": (3.2, 1.0), "hspace": 0.05},
    )
    for channel in POST_CHANNELS:
        ax.stairs(
            sim[channel].counts * sim_factor, edges, label=POST_CHANNEL_LABELS[channel],
            color=POST_CHANNEL_COLORS[channel], linewidth=1.15, alpha=0.9,
        )
    sim_err = np.sqrt(np.maximum(0.0, total_sim_var))
    ax.stairs(total_sim, edges, label="SIMC total", color=POST_SIM_COLOR, linewidth=2.0)
    ax.fill_between(centers, total_sim - sim_err, total_sim + sim_err, step="mid", color=POST_SIM_COLOR, alpha=0.18, linewidth=0)
    data_err = np.sqrt(np.maximum(0.0, data.sumw2))
    ax.errorbar(centers, data.counts, yerr=data_err, fmt="o", ms=3.2, color=POST_DATA_COLOR, ecolor=POST_DATA_COLOR, elinewidth=0.8, capsize=0, label="Data", zorder=5)
    ylabel = "Corrected yield / 1 mC / bin" if normalization == "absolute" else "Area-matched weighted yield / bin"
    ax.set_ylabel(ylabel)
    ax.set_title(f"{kin}: {spec.title}", pad=9)
    ax.legend(ncol=2, fontsize=8.5, loc="best")
    post_style_axis(ax)

    valid = np.isfinite(total_sim) & (total_sim != 0.0)
    ratio = np.divide(data.counts, total_sim, out=np.full(spec.bins, np.nan), where=valid)
    ratio_var = np.divide(data.sumw2, total_sim**2, out=np.zeros(spec.bins), where=valid)
    ratio_var += np.divide(data.counts**2 * total_sim_var, total_sim**4, out=np.zeros(spec.bins), where=valid)
    ratio_ax.errorbar(centers[valid], ratio[valid], yerr=np.sqrt(np.maximum(0.0, ratio_var[valid])), fmt="o", ms=3.0, color=POST_DATA_COLOR, elinewidth=0.75)
    ratio_ax.axhline(1.0, color=POST_SIM_COLOR, linestyle="--", linewidth=1.0)
    ratio_ax.set_ylabel("Data/SIMC")
    ratio_ax.set_xlabel(spec.xlabel)
    ratio_ax.set_xlim(spec.limits)
    finite_ratio = ratio[np.isfinite(ratio)]
    if finite_ratio.size:
        upper = min(4.0, max(1.8, float(np.nanpercentile(finite_ratio, 95)) * 1.25))
        ratio_ax.set_ylim(0.0, upper)
    post_style_axis(ratio_ax)
    note = "pre-exclusive selection" if spec.key == "mmiss" else "basic HMS + NPS event selection"
    if normalization == "shape":
        note += f"; SIMC global factor={sim_factor:.4g}"
    ax.text(0.015, 0.975, note, transform=ax.transAxes, va="top", ha="left", fontsize=8.2, color="#444444")
    fig.align_ylabels((ax, ratio_ax))
    return fig, sim_factor


def post_write_outputs(output_dir: Path, data: dict[str, PostHist], sim: dict[str, dict[str, PostHist]], kin: str, normalization: str, formats: set[str]) -> tuple[int, dict[str, float]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / f"post_analysis_{kin}_{normalization}.pdf"
    factors: dict[str, float] = {}
    with PdfPages(pdf_path) as multipage:
        for index, spec in enumerate(POST_SPECS, start=1):
            per_channel = {channel: sim[channel][spec.key] for channel in POST_CHANNELS}
            fig, factor = post_plot_comparison(spec, data[spec.key], per_channel, kin, normalization)
            factors[spec.key] = factor
            stem = f"{index:02d}_{spec.group}_{spec.key}_{normalization}"
            multipage.savefig(fig, bbox_inches="tight")
            for extension in formats:
                fig.savefig(output_dir / f"{stem}.{extension}", dpi=FIGURE_DPI, bbox_inches="tight")
            plt.close(fig)
            print(f"[SAVE] {stem}")
    print(f"[SAVE] {pdf_path}")
    return len(POST_SPECS), factors


def post_main(argv: list[str] | None = None) -> int:
    args = post_build_parser().parse_args(argv)
    combined_root = Path(args.combined_root).expanduser().resolve()
    if not combined_root.is_file():
        raise FileNotFoundError(f"Combined ROOT file not found: {combined_root}")
    kin = args.kin or post_infer_kin(combined_root)
    token = post_kin_token(kin)
    kin_config = Path(args.kin_config).expanduser().resolve()
    hist_dir = Path(args.hist_dir).expanduser().resolve()
    kinematics = post_load_kinematics(kin_config, kin)
    production, sim_files = post_find_simulation_files(Path(args.simulation_base), args.simulation_production, token)
    hist_files = {channel: hist_dir / f"nps_{channel}_pi0_{token}.hist" for channel in POST_CHANNELS}
    missing_hist = [str(path) for path in hist_files.values() if not path.is_file()]
    if missing_hist:
        raise FileNotFoundError("Missing SIMC .hist file(s): " + ", ".join(missing_hist))
    metadata = {channel: post_read_hist_meta(path) for channel, path in hist_files.items()}
    for channel, meta in metadata.items():
        hist_ebeam = meta["ebeam_mev"] / 1000.0
        if not np.isclose(hist_ebeam, kinematics["ebeam_gev"], rtol=0.0, atol=0.005):
            raise ValueError(f"{channel} .hist Ebeam={hist_ebeam} GeV disagrees with config Ebeam={kinematics['ebeam_gev']} GeV")

    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else combined_root.parent.parent / "plots" / "post_analysis"
    stability_dir = output_dir / "run_stability"
    summary_csv = (
        Path(args.summary_csv).expanduser().resolve()
        if args.summary_csv
        else combined_root.parent.parent / "summary" / "summary_all_runs.csv"
    )
    efficiency_csv = (
        Path(args.efficiency_csv).expanduser().resolve()
        if args.efficiency_csv
        else Path(__file__).resolve().parents[1] / "output" / "efficiency_stuff" / f"efficiency_{kin}.csv"
    )
    if not args.skip_run_stability:
        if not summary_csv.is_file():
            raise FileNotFoundError(f"Run stability summary CSV not found: {summary_csv}; pass --summary-csv or --skip-run-stability")
        if not efficiency_csv.is_file():
            raise FileNotFoundError(f"Run stability efficiency CSV not found: {efficiency_csv}; pass --efficiency-csv or --skip-run-stability")
    formats = {item.strip().lower() for item in args.formats.split(",") if item.strip()}
    invalid_formats = formats - {"png", "pdf"}
    if invalid_formats:
        raise ValueError(f"Unsupported --formats: {', '.join(sorted(invalid_formats))}")
    if args.max_events:
        print("[WARN] --max-events makes absolute normalization incomplete; use only for debugging.")

    print(f"[INFO] kin={kin} combined={combined_root}")
    print(f"[INFO] simulation production={production}")
    print(f"[INFO] Ebeam={kinematics['ebeam_gev']:.6g} GeV, NPS theta={kinematics['nps_theta_deg']:.6g} deg, distance={kinematics['nps_target_distance_cm']:.6g} cm")
    data_hist, data_audit = post_read_data(combined_root, args.step_size, args.max_events)
    sim_hist: dict[str, dict[str, PostHist]] = {}
    sim_audit: dict[str, dict[str, object]] = {}
    for channel in POST_CHANNELS:
        print(f"[INFO] reading {channel}: {sim_files[channel]}")
        sim_hist[channel], sim_audit[channel] = post_read_simulation(
            sim_files[channel], metadata[channel], kinematics, args.step_size, args.max_events
        )
        print(f"[INFO] {channel}: selected {sim_audit[channel]['entries_selected']}/{sim_audit[channel]['entries_read']}")

    plot_count, shape_factors = post_write_outputs(output_dir, data_hist, sim_hist, kin, args.normalization, formats)
    stability_enabled = not args.skip_run_stability
    if stability_enabled:
        stability_argv = [
            "--summary-csv", str(summary_csv),
            "--efficiency-csv", str(efficiency_csv),
            "--output-dir", str(stability_dir),
            "--ignore-runs", args.ignore_runs,
            "--include-runs", ",".join(str(run) for run in data_audit["run_numbers"]),
        ]
        if args.run_status:
            stability_argv.extend(("--run-status", args.run_status))
        print(f"[INFO] creating run stability plots in {stability_dir}")
        if legacy_main(stability_argv, prompt_for_runs=False) != 0:
            raise ValueError("Run stability plotting failed")
    audit = {
        "kin": kin, "normalization": args.normalization,
        "data_weight_formula": "scale*pi0_weight*(charge_uC/total_unique_run_charge_uC)",
        "simulation_weight_formula": "normfac*Weight/Ngen",
        "missing_mass_exclusive_cut": False,
        "kinematics": kinematics, "combined_root": str(combined_root),
        "simulation_production": str(production), "hist_dir": str(hist_dir),
        "run_stability": {
            "enabled": stability_enabled,
            "summary_csv": str(summary_csv),
            "efficiency_csv": str(efficiency_csv),
            "output_dir": str(stability_dir),
        },
        "data": data_audit, "simulation": sim_audit,
        "shape_factors": shape_factors if args.normalization == "shape" else {},
        "plots": plot_count,
    }
    audit_path = output_dir / f"post_analysis_{kin}_{args.normalization}_metadata.json"
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(f"[SAVE] {audit_path}")
    print(f"[DONE] {plot_count} comparison plots in {output_dir}")
    if stability_enabled:
        print(f"[DONE] run/current/rate stability and normalized-yield plots in {stability_dir}")
    return 0


def main() -> int:
    if "--legacy-run-qa" in sys.argv:
        sys.argv.remove("--legacy-run-qa")
        return legacy_main()
    try:
        return post_main()
    except (FileNotFoundError, KeyError, ValueError, OSError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
