
from __future__ import annotations
"""
HOW TO USE:
-----------
This script generates per-run scaler/trigger count and livetime plots from a CSV file.

Basic usage (from the command line):
    python plot_per_run_scalers_effs.py --csv <input_csv> --setting <kinematic_setting> --outdir <output_dir>

Arguments:
    --csv      Path to the input CSV file (default: livetime_results_parallel.csv)
    --setting  Kinematic setting label to filter (default: KinC_x60_4b)
    --outdir   Output directory for PNG and PDF plots (default: eff_plots/per_run_scaler_trigger_vs_charge)

Example:
    python plot_per_run_scalers_effs.py --csv livetime_results_parallel.csv --setting KinC_x60_4b --outdir my_plots

The script will:
    - Filter the CSV for the specified kinematic setting (if present)
    - Exclude suspicious runs (negative livetime, abnormally high counts)
    - Generate PNG and PDF plots for scaler/trigger counts and livetimes
    - Output plots to the specified directory, organized by target (LD2/LH2)

Requirements:
    - Python 3
    - pandas, matplotlib
    - Input CSV must have columns: run, charge_uC, scaler_*, trig_*, livetime*, etc.
    - Config CSV path is hardcoded in the script (edit if needed)
"""
#!/usr/bin/env python3
"""Create per-run scaler/trigger count vs charge plots from luminosity CSV results.

Default behavior targets KinC_x60_4b. If a kinematic column exists in the CSV,
rows are filtered to the requested setting. Otherwise the script uses all rows.
"""

import argparse
import math
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make per-run scaler/trigger count vs charge plots."
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("livetime_results_parallel.csv"),
        help="Input CSV path (default: livetime_results_parallel.csv)",
    )
    parser.add_argument(
        "--setting",
        default="KinC_x60_4b",
        help="Kinematic setting label to filter and include in plot titles.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("eff_plots/per_run_scaler_trigger_vs_charge"),
        help="Output directory for PNG plots.",
    )
    return parser.parse_args()


def detect_kinematic_column(df: pd.DataFrame) -> str | None:
    candidates = [
        "kinematic_setting",
        "kin_setting",
        "setting",
        "kinematic",
        "kin",
    ]
    for col in candidates:
        if col in df.columns:
            return col
    return None


def detect_current_column(df: pd.DataFrame) -> str | None:
    candidates = [
        "avg_current_uA",
        "expected_current_uA",
        "current_uA",
        "beam_current_uA",
        "current",
    ]
    for col in candidates:
        if col in df.columns and pd.api.types.is_numeric_dtype(df[col]):
            return col
    return None


def get_count_columns(df: pd.DataFrame) -> List[str]:
    # Keep only scaler/trigger count-like columns requested by the user.
    cols = [
        col
        for col in df.columns
        if (col.startswith("scaler_") or col.startswith("trig_"))
        and pd.api.types.is_numeric_dtype(df[col])
    ]
    return sorted(cols)


def get_livetime_columns(df: pd.DataFrame) -> List[str]:
    cols = [
        col
        for col in df.columns
        if ("livetime" in col.lower()) and pd.api.types.is_numeric_dtype(df[col])
    ]
    return sorted(cols)


def plot_single_quantity(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    setting: str,
    outdir: Path,
    pdf: PdfPages | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))

    sc = ax.scatter(
        df[x_col],
        df[y_col],
        c=df["run"],
        cmap="viridis",
        s=70,
        alpha=0.9,
        edgecolors="black",
        linewidths=0.4,
    )



    ax.set_title(f"{y_col} vs {x_col} ({setting})", fontsize=13)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    if "hms_track_eff" in y_col.lower():
        ax.set_ylim(0.9925, 1.0)
    ax.grid(True, alpha=0.25)

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("run")

    fig.tight_layout()
    out_path = outdir / f"{y_col}_vs_charge.png"
    fig.savefig(out_path, dpi=160)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)

def plot_multipanel_livetime(
    df: pd.DataFrame,
    x_col: str,
    livetime_cols: List[str],
    setting: str,
    outdir: Path,
    pdf: PdfPages | None = None,
) -> None:
    # For each livetime column, make two panels: >0.9 and <=0.9, both vs current_col
    n = len(livetime_cols)
    ncols = 2
    nrows = n
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(13, 4.1 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    for idx, y_col in enumerate(livetime_cols):
        # Panel 1: >0.9
        ax_high = axes[idx][0]
        mask_high = df[y_col] > 0.9
        sc_high = ax_high.scatter(
            df.loc[mask_high, x_col],
            df.loc[mask_high, y_col],
            c=df.loc[mask_high, "run"],
            cmap="viridis",
            s=35,
            alpha=0.85,
            edgecolors="none",
        )
        ax_high.set_title(f"{y_col} > 0.9", fontsize=10)
        ax_high.set_xlabel(x_col)
        ax_high.set_ylabel("livetime")
        if "hms_track_eff" in y_col.lower():
            ax_high.set_ylim(0.9925, 1.0)
        ax_high.grid(True, alpha=0.2)
        # Panel 2: <=0.9
        ax_low = axes[idx][1]
        mask_low = df[y_col] <= 0.9
        sc_low = ax_low.scatter(
            df.loc[mask_low, x_col],
            df.loc[mask_low, y_col],
            c=df.loc[mask_low, "run"],
            cmap="viridis",
            s=35,
            alpha=0.85,
            edgecolors="none",
        )
        ax_low.set_title(f"{y_col} ≤ 0.9", fontsize=10)
        ax_low.set_xlabel(x_col)
        ax_low.set_ylabel("livetime")
        if "hms_track_eff" in y_col.lower():
            ax_low.set_ylim(0.9925, 1.0)
        ax_low.grid(True, alpha=0.2)
    # Hide any unused subplot panes (not needed here since nrows=n, ncols=2)
    # Add colorbars
    for idx, y_col in enumerate(livetime_cols):
        fig.colorbar(axes[idx][0].collections[0], ax=axes[idx][0], fraction=0.018, pad=0.01, label="run")
        fig.colorbar(axes[idx][1].collections[0], ax=axes[idx][1], fraction=0.018, pad=0.01, label="run")
    fig.suptitle(f"Livetime vs {x_col} ({setting})", fontsize=14)
    fig.savefig(outdir / f"all_livetime_vs_{x_col}.png", dpi=170)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)

def plot_multipanel(
    df: pd.DataFrame,
    x_col: str,
    y_cols: List[str],
    setting: str,
    outdir: Path,
    pdf: PdfPages | None = None,
) -> None:
    n = len(y_cols)
    ncols = 2
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(13, 4.1 * nrows),
        squeeze=False,
        constrained_layout=True,
    )

    sorted_df = df.sort_values("run")

    for idx, y_col in enumerate(y_cols):
        r = idx // ncols
        c = idx % ncols
        ax = axes[r][c]

        sc = ax.scatter(
            df[x_col],
            df[y_col],
            c=df["run"],
            cmap="viridis",
            s=35,
            alpha=0.85,
            edgecolors="none",
        )
        ax.set_title(y_col, fontsize=10)
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        if "hms_track_eff" in y_col.lower():
            ax.set_ylim(0.9925, 1.0)
        ax.grid(True, alpha=0.2)

    # Hide any unused subplot panes.
    for idx in range(n, nrows * ncols):
        r = idx // ncols
        c = idx % ncols
        axes[r][c].axis("off")

    # One shared colorbar for run number.
    fig.colorbar(sc, ax=axes, fraction=0.018, pad=0.01, label="run")
    fig.suptitle(f"Scaler/Trigger counts vs charge_uC ({setting})", fontsize=14)
    fig.savefig(outdir / "all_scaler_trigger_counts_vs_charge.png", dpi=170)
    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)


def plot_livetime_overlays(
    df: pd.DataFrame,
    livetime_cols: List[str],
    current_col: str,
    setting: str,
    outdir: Path,
    pdf: PdfPages | None = None,
) -> int:
    if not livetime_cols:
        print("No livetime columns found; skipping livetime overlay plots.")
        return 0

    # Plot only high-livetime points (> 0.9) as a scatter overlay vs current.
    fig_high, ax_high = plt.subplots(figsize=(10.5, 6.2))
    n_high = 0
    for col in livetime_cols:
        mask_high = df[col] > 0.9
        if not mask_high.any():
            continue
        n_high += int(mask_high.sum())
        ax_high.scatter(
            df.loc[mask_high, current_col],
            df.loc[mask_high, col],
            s=26,
            alpha=0.85,
            label=col,
        )
    ax_high.set_title(f"Livetime overlays (>0.9) vs {current_col} ({setting})", fontsize=13)
    ax_high.set_xlabel(current_col)
    ax_high.set_ylabel("livetime")
    ax_high.grid(True, alpha=0.25)
    if n_high > 0:
        ax_high.legend(loc="best", fontsize=8)
    fig_high.tight_layout()
    fig_high.savefig(outdir / "overlay_livetimes_gt_0p9_vs_current.png", dpi=170)
    if pdf is not None:
        pdf.savefig(fig_high)
    plt.close(fig_high)

    # Plot low-livetime points (<= 0.9) as scatter and annotate run at each point.
    fig_low, ax_low = plt.subplots(figsize=(10.5, 6.2))
    n_low = 0
    for col in livetime_cols:
        low_df = df[df[col] <= 0.9][["run", current_col, col]].copy()
        if low_df.empty:
            continue
        n_low += len(low_df)
        ax_low.scatter(
            low_df[current_col],
            low_df[col],
            s=30,
            alpha=0.9,
            label=col,
        )
        for _, row in low_df.iterrows():
            ax_low.annotate(
                str(int(row["run"])),
                (row[current_col], row[col]),
                textcoords="offset points",
                xytext=(4, 3),
                fontsize=7,
                alpha=0.9,
            )

    ax_low.set_title(f"Livetime overlays (<=0.9) vs {current_col} with run labels ({setting})", fontsize=13)
    ax_low.set_xlabel(current_col)
    ax_low.set_ylabel("livetime")
    ax_low.grid(True, alpha=0.25)
    if n_low > 0:
        ax_low.legend(loc="best", fontsize=8)
    fig_low.tight_layout()
    fig_low.savefig(outdir / "overlay_livetimes_le_0p9_vs_current_with_runs.png", dpi=170)
    if pdf is not None:
        pdf.savefig(fig_low)
    plt.close(fig_low)

    return 2


def main() -> None:
    # Dynamically flag suspicious runs (no hard-coded list)
    extra_suspicious_runs = set()

    args = parse_args()

    if not args.csv.exists():
        raise FileNotFoundError(f"CSV not found: {args.csv}")

    # Load main data CSV
    df = pd.read_csv(args.csv)

    required_cols = {"run", "charge_uC"}
    missing = sorted(required_cols - set(df.columns))
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")

    # Force use of 'expected_current_uA' as the current column
    if 'expected_current_uA' not in df.columns or not pd.api.types.is_numeric_dtype(df['expected_current_uA']):
        raise ValueError("CSV missing required column: expected_current_uA")
    current_col = 'expected_current_uA'

    kin_col = detect_kinematic_column(df)
    if kin_col is not None:
        filtered = df[df[kin_col] == args.setting].copy()
        if filtered.empty:
            available = sorted(df[kin_col].dropna().astype(str).unique().tolist())
            raise ValueError(
                f"No rows found for {args.setting} in column {kin_col}. Available values: {available}"
            )
        df = filtered
        print(f"Filtered to {len(df)} rows using {kin_col} == {args.setting}")
    else:
        print(
            "No kinematic-setting column found. Using all rows and labeling plots with "
            f"setting={args.setting}."
        )

    # Load config CSV for run filtering and target info
    config_path = Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/config/nps_dvcs_all_kins_main.csv")
    if not config_path.exists():
        raise FileNotFoundError(f"Config CSV not found: {config_path}")
    config_df = pd.read_csv(config_path)
    # Strip whitespace from all column names (fixes 'target ')
    config_df.columns = config_df.columns.str.strip()

    # Only keep runs present in both, and not junk, and with target LD2 or LH2

    print("[DEBUG] Selection criteria applied to config CSV:")
    print("  - Type != 'junk'")
    print("  - Type == 'production'")
    print("  - target in ['LD2', 'LH2']")
    print("  - run_number convertible to int and not NA")

    print(f"[DEBUG] Config CSV: {len(config_df)} rows before filtering")
    config_df = config_df[config_df["Type"].str.lower() != "junk"]
    print(f"[DEBUG] Config CSV: {len(config_df)} rows after removing 'junk'")
    config_df = config_df[config_df["Type"] == "production"]
    print(f"[DEBUG] Config CSV: {len(config_df)} rows after keeping only 'production' type")
    config_df = config_df[config_df["target"].isin(["LD2", "LH2"])]
    print(f"[DEBUG] Config CSV: {len(config_df)} rows after keeping only LD2/LH2 targets")
    # Ensure run_number is int for merge
    config_df = config_df.copy()
    config_df["run_number"] = pd.to_numeric(config_df["run_number"], errors="coerce")
    config_df = config_df.dropna(subset=["run_number"])
    config_df["run_number"] = config_df["run_number"].astype(int)
    print(f"[DEBUG] Config CSV: {len(config_df)} rows after dropping NA run_number")

    # Ensure df["run"] is int
    df = df.copy()
    df["run"] = pd.to_numeric(df["run"], errors="coerce")
    df = df.dropna(subset=["run"])
    df["run"] = df["run"].astype(int)
    print(f"[DEBUG] Main CSV: {len(df)} runs after cleaning run column")

    # Merge on run/run_number
    merged = pd.merge(df, config_df[["run_number", "target"]], left_on="run", right_on="run_number", how="inner")
    print(f"[DEBUG] Merged: {len(merged)} runs present in both CSVs after filtering")
    print(f"[DEBUG] LH2 runs: {sum(merged['target'] == 'LH2')}, LD2 runs: {sum(merged['target'] == 'LD2')}")

    # --- Suspicious run detection ---
    # 1. Flag runs with any negative livetime
    livetime_cols = get_livetime_columns(merged)
    for col in livetime_cols:
        neg_livetime_runs = merged[merged[col] < 0]["run"].unique()
        if len(neg_livetime_runs) > 0:
            print(f"[WARNING] Runs with negative {col}: {neg_livetime_runs.tolist()}")
            extra_suspicious_runs.update(neg_livetime_runs)

    # 2. Flag runs with any scaler/trigger count > 10,000,000
    count_cols = get_count_columns(merged)
    for col in count_cols:
        high_count_runs = merged[merged[col] > 10_000_000]["run"].unique()
        if len(high_count_runs) > 0:
            print(f"[WARNING] Runs with high {col} (>10,000,000): {high_count_runs.tolist()}")
            extra_suspicious_runs.update(high_count_runs)

    all_suspicious_runs = extra_suspicious_runs
    print(f"[DEBUG] All suspicious runs (flagged): {sorted(all_suspicious_runs)}")

    # Exclude suspicious runs for plotting
    merged_no_suspicious = merged[~merged["run"].isin(all_suspicious_runs)].copy()
    print(f"[DEBUG] After excluding suspicious runs: {len(merged_no_suspicious)} runs remain")
    print(f"[DEBUG] LH2 runs (no suspicious): {sum(merged_no_suspicious['target'] == 'LH2')}, LD2 runs (no suspicious): {sum(merged_no_suspicious['target'] == 'LD2')}")

    # Only make plots for non-suspicious runs
    for target in ["LD2", "LH2"]:
        target_df_no_susp = merged_no_suspicious[merged_no_suspicious["target"] == target].copy()
        print(f"[DEBUG] {target}: {len(target_df_no_susp)} runs to plot (excluding suspicious)")
        if not target_df_no_susp.empty:
            target_outdir_no_susp = args.outdir / f"target_{target}_no_suspicious"
            target_outdir_no_susp.mkdir(parents=True, exist_ok=True)
            pdf_path_no_susp = target_outdir_no_susp / f"all_scaler_trigger_livetime_plots_{target}_no_suspicious.pdf"

            y_cols = get_count_columns(target_df_no_susp)
            livetime_cols = get_livetime_columns(target_df_no_susp)
            hms_eff_cols = [col for col in target_df_no_susp.columns if "hms_track_eff" in col.lower()]
            with PdfPages(pdf_path_no_susp) as pdf:
                # Only save multipanel scaler/trigger plots in PDF
                if y_cols:
                    plot_multipanel(target_df_no_susp, "charge_uC", y_cols, args.setting + f" ({target}, no suspicious)", target_outdir_no_susp, pdf=pdf)
                    plot_multipanel(target_df_no_susp, "run", y_cols, args.setting + f" ({target}, no suspicious)", target_outdir_no_susp, pdf=pdf)
                # Livetime multipanel (vs current and vs run)
                if livetime_cols:
                    plot_multipanel_livetime(target_df_no_susp, current_col, livetime_cols, args.setting + f" ({target}, no suspicious)", target_outdir_no_susp, pdf=pdf)
                    plot_multipanel_livetime(target_df_no_susp, "run", livetime_cols, args.setting + f" ({target}, no suspicious)", target_outdir_no_susp, pdf=pdf)
                    # Also overlay plots as before
                    plot_livetime_overlays(
                        target_df_no_susp,
                        livetime_cols,
                        current_col,
                        args.setting + f" ({target}, no suspicious)",
                        target_outdir_no_susp,
                        pdf=pdf,
                    )
                # HMS track efficiency plots (vs charge and vs run)
                for eff_col in hms_eff_cols:
                    plot_single_quantity(target_df_no_susp, "charge_uC", eff_col, args.setting + f" ({target}, no suspicious)", target_outdir_no_susp, pdf=pdf)
                    plot_single_quantity(target_df_no_susp, "run", eff_col, args.setting + f" ({target}, no suspicious)", target_outdir_no_susp, pdf=pdf)
                    # vs current
                    if "expected_current_uA" in target_df_no_susp.columns:
                        plot_single_quantity(
                            target_df_no_susp,
                            "expected_current_uA",
                            eff_col,
                            args.setting + f" ({target}, no suspicious)",
                            target_outdir_no_susp,
                            pdf=pdf,
                        )
            print(f"Saved plots to: {target_outdir_no_susp}")
            print(f"Saved combined PDF: {pdf_path_no_susp}")
        else:
            print(f"No runs for target {target} after filtering (no suspicious runs).")


if __name__ == "__main__":
    main()
