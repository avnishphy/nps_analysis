#!/usr/bin/env python3
# Standalone use from src/efficiencies:
#   python3 plot_compare_livetimes.py --output-dir ../../output/efficiency_stuff/plots \
#     ../../output/efficiency_stuff/compare_livetimes_KinC_x60_4b.csv
# Pass more merged CSV paths to plot more kinematics; writes one 3-page PDF per target.
# Current and rate pages read selection_report_<kin>.csv and efficiency_<kin>.csv
# beside each comparison CSV. Values above 1.01 get outlier flags, not points.
# No PNG files are produced.
"""Make run, current and HMS S1X-rate multipanel plots of six livetimes."""

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")  # Batch jobs do not have a display.
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D


# Names match compare_livetimes.cxx; colors/markers mirror efficiency overlays.
METRICS = (
    ("total_edtm_lt", "Total EDTM"),
    ("clta_tsh_lt", "Computer all (TSH)"),
    ("clta_tdc_lt", "Computer all (TDC)"),
    ("cltp_tdc_lt", "Computer physics (TDC)"),
    ("newgen_edtm_lt", "NewGen EDTM"),
    ("matched_edtm_lt", "Matched EDTM"),
)
COLORS = ("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#1f78b4")
MARKERS = ("o", "s", "^", "D", "v", "P")
OUTLIER_LIMIT = 1.01
# Light backgrounds distinguish all possible prescale blocks on the run page.
PS_SHADE = {"ps1": "#edf8f1", "ps2": "#fff0f5",
            "ps3": "#fff7e6", "ps4": "#eef7ff",
            "ps5": "#f4efff", "ps6": "#fff1e8"}
# Prescale is encoded redundantly by a large outline shape and its color.
# The inner marker/color continues to identify the livetime definition.
PS_OUTLINE = {"ps1": "#007a4d", "ps2": "#c2185b",
              "ps3": "#7a3e00", "ps4": "#0057b8",
              "ps5": "#6f2da8", "ps6": "#d55e00"}
PS_MARKERS = {"ps1": "^", "ps2": "v", "ps3": "D",
              "ps4": "o", "ps5": "h", "ps6": "s"}


def safe_name(value):
    """Keep target names usable as filenames without changing CSV labels."""
    return "".join(ch if ch.isalnum() else "_" for ch in value).strip("_") or "unknown"


def prescale_group(row):
    """Trigger n in the comparison CSV corresponds to prescale psn."""
    return f"ps{row['trigger'].strip()}" if row.get("trigger", "").strip().isdigit() else "unknown"


def prescale_run_text(rows):
    """List exact contiguous run ranges for each prescale group."""
    grouped = defaultdict(list)
    for row in rows:
        grouped[prescale_group(row)].append(int(row["run"]))
    lines = ["Prescale runs"]
    for group in sorted(grouped):
        ordered = sorted(set(grouped[group]))
        ranges = []
        start = previous = ordered[0]
        for run in ordered[1:] + [None]:
            if run is not None and run == previous + 1:
                previous = run
                continue
            ranges.append(str(start) if start == previous else f"{start}-{previous}")
            if run is not None:
                start = previous = run
        chunks = [", ".join(ranges[i:i + 3]) for i in range(0, len(ranges), 3)]
        lines.append(f"{group} ({len(ordered)}): {chunks[0]}")
        lines.extend(f"  {chunk}" for chunk in chunks[1:])
    return "\n".join(lines)


def read_beam_currents(csv_path):
    """Average measured segment currents per run, as in efficiency overlays."""
    prefix = "compare_livetimes_"
    if not csv_path.name.startswith(prefix):
        return {}
    report = csv_path.with_name("selection_report_" + csv_path.name[len(prefix):])
    if not report.is_file():
        print(f"[plot] No beam-current report: {report}")
        return {}
    segments = defaultdict(list)
    with report.open(newline="") as source:
        for raw in csv.DictReader(source):
            # The selection report may have spaces or quotes around cell values.
            row = {key.strip(): (value or "").strip().strip('"')
                   for key, value in raw.items() if key is not None}
            if row.get("selection_ok") != "1":
                continue
            try:
                run = int(row["run_number"])
                current = float(row["mean_current_uA"])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(current) and current > 0:
                segments[run].append(current)
    return {run: mean(values) for run, values in segments.items()}


def read_s1x_rates(csv_path):
    """Use the per-run time-weighted HMS S1X scaler rate from efficiency CSV."""
    prefix = "compare_livetimes_"
    if not csv_path.name.startswith(prefix):
        return {}
    source_path = csv_path.with_name("efficiency_" + csv_path.name[len(prefix):])
    if not source_path.is_file():
        print(f"[plot] No HMS S1X-rate CSV: {source_path}")
        return {}
    rates = {}
    with source_path.open(newline="") as source:
        for raw in csv.DictReader(source):
            row = {key.strip(): (value or "").strip().strip('"')
                   for key, value in raw.items() if key is not None}
            try:
                run = int(row["run_number"])
                rate = float(row["HMS_S1X_rate_Hz"])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(rate) and rate >= 0:
                rates[run] = rate
    return rates


def metric_series(rows, x_values):
    """Keep ordinary points separate from high livetime outliers."""
    series = []
    for column, _ in METRICS:
        points, outliers = [], []
        for x, row in zip(x_values, rows):
            try:
                value = float(row[column])
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                point = (x, value, int(row["run"]), prescale_group(row))
                (outliers if value > OUTLIER_LIMIT else points).append(point)
        series.append((points, outliers))
    return series


def shade_prescale_blocks(axes, combined, rows):
    """Shade consecutive prescale runs, matching efficiency run plots."""
    states = [prescale_group(row) for row in rows]
    start = 0
    for end in range(1, len(states) + 1):
        if end < len(states) and states[end] == states[start]:
            continue
        group = states[start]
        for ax in axes:
            ax.axvspan(start - 0.5, end - 0.5,
                       color=PS_SHADE.get(group, "#f4f4f4"), alpha=0.45, zorder=-2)
            if end < len(states):
                ax.axvline(end - 0.5, color="#777777", linestyle=":",
                           linewidth=0.8, alpha=0.55)
        combined.text((start + end - 1) / 2, 0.98, group,
                      transform=combined.get_xaxis_transform(), ha="center",
                      va="top", fontsize=8,
                      bbox={"facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.9})
        start = end


def save_multipanel(series, kin, target, xlabel, pdf, width, rows, runs=None):
    """Draw six individual definitions and a shared overlay in one figure."""
    multi = plt.figure(figsize=(width, 13))
    grid = multi.add_gridspec(4, 2, height_ratios=(1, 1, 1, 1.5))
    panels = [multi.add_subplot(grid[row, col]) for row in range(3) for col in range(2)]
    combined = multi.add_subplot(grid[3, :])
    groups = Counter(prescale_group(row) for row in rows)
    outlier_notes = []
    # Connect neighboring runs in run order; current points are independent runs.
    style = "-" if runs is not None else "None"
    for i, ((_, label), (points, outliers)) in enumerate(zip(METRICS, series)):
        panel = panels[i]
        if points:
            # NaN breaks a run-order line at each excluded outlier.
            line_points = [(x, y) for x, y, _, _ in points]
            if runs is not None:
                line_points += [(x, math.nan) for x, _, _, _ in outliers]
                line_points.sort(key=lambda point: point[0])
            xp, yp = zip(*line_points)
            for ax, legend_label in ((panel, None), (combined, label)):
                ax.plot(xp, yp, linestyle=style, color=COLORS[i],
                        marker=MARKERS[i], markersize=4, linewidth=1,
                        label=legend_label)
                # A colored hollow ring identifies the prescale of each run.
                for group in groups:
                    marked = [(x, y) for x, y, _, ps in points if ps == group]
                    if marked:
                        mx, my = zip(*marked)
                        ax.plot(mx, my, linestyle="None",
                                marker=PS_MARKERS.get(group, "h"), markersize=8.5,
                                markerfacecolor="none",
                                markeredgecolor=PS_OUTLINE.get(group, "#777777"),
                                markeredgewidth=1.6, zorder=4)
        else:
            panel.text(0.5, 0.5, "Only outliers" if outliers else "No finite values",
                       transform=panel.transAxes,
                       ha="center", va="center")
        if outliers:
            ox = [x for x, _, _, _ in outliers]
            for ax, height in ((panel, 0.94), (combined, 0.94 - 0.04 * (i % 3))):
                # The X is an axis flag; it does not plot the measured livetime.
                ax.plot(ox, [height] * len(ox), transform=ax.get_xaxis_transform(),
                        color="#b22222", marker="x", markersize=7,
                        linestyle="None", markeredgewidth=1.5)
            outlier_notes.extend(f"{label}: run {run} = {value:.4f}"
                                 for _, value, run, _ in outliers)
        panel.axhline(1.0, color="0.4", linestyle="--", linewidth=0.8)
        panel.set_title(label, fontsize=10)
        panel.set_ylabel("Livetime")
        panel.grid(axis="y", color="0.85", linewidth=0.6)
        panel.tick_params(labelbottom=False)
    combined.axhline(1.0, color="0.4", linestyle="--", linewidth=0.8)
    combined.set_title("All definitions", fontsize=10)
    combined.set_ylabel("Livetime")
    combined.set_xlabel(xlabel)
    combined.grid(axis="y", color="0.85", linewidth=0.6)
    combined.legend(ncol=3, fontsize=8)
    if runs is not None:
        shade_prescale_blocks((*panels, combined), combined, rows)
        # Every run gets a tick. Compact positions avoid large numeric run gaps.
        ticks = list(range(len(runs)))
        for ax in (*panels, combined):
            ax.set_xlim(-0.5, len(runs) - 0.5)
            ax.set_xticks(ticks)
        combined.set_xticklabels([str(run) for run in runs], rotation=90,
                                 ha="center", fontsize=7)
    # Keep the prescale key and exact outlier values beside every page.
    ps_handles = [Line2D([], [], linestyle="None",
                         marker=PS_MARKERS.get(group, "h"), markersize=10,
                         markerfacecolor="none",
                         markeredgecolor=PS_OUTLINE.get(group, "#777777"),
                         markeredgewidth=1.8,
                         label=f"{group}: {groups[group]} runs") for group in sorted(groups)]
    multi.legend(handles=ps_handles, loc="upper left", bbox_to_anchor=(0.85, 0.93),
                 title="Prescale (shape + outline)", fontsize=8, title_fontsize=8)
    multi.text(0.85, 0.83, prescale_run_text(rows), ha="left", va="top", fontsize=7.2,
               bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.94})
    if outlier_notes:
        lines = outlier_notes[:10]
        if len(outlier_notes) > 10:
            lines.append(f"+{len(outlier_notes) - 10} more")
        multi.text(0.85, 0.62, "X at top: outlier >1.01\n" + "\n".join(lines),
                   ha="left", va="top", fontsize=8,
                   bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.94})
    multi.suptitle(f"{kin} | {target} | {len(rows)} runs | {xlabel}")
    multi.tight_layout(rect=(0, 0, 0.84, 0.97))
    pdf.savefig(multi)
    plt.close(multi)


def plot_csv(csv_path, output_dir):
    with csv_path.open(newline="") as source:
        reader = csv.DictReader(source)
        needed = {"run", "kin", "target", *(name for name, _ in METRICS)}
        missing = needed - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{csv_path}: missing columns {sorted(missing)}")
        # Mixed targets share one kinematic CSV; separate them for readable plots.
        by_group = defaultdict(list)
        for row in reader:
            by_group[(row["kin"].strip(), row["target"].strip())].append(row)

    beam_currents = read_beam_currents(csv_path)
    s1x_rates = read_s1x_rates(csv_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    for (kin, target), rows in sorted(by_group.items()):
        rows.sort(key=lambda row: int(row["run"]))
        runs = [int(row["run"]) for row in rows]
        x = list(range(len(rows)))  # Large numeric run gaps do not flatten trends.
        series = metric_series(rows, x)
        if not any(points or outliers for points, outliers in series):
            continue
        suffix = f"{safe_name(kin)}_{safe_name(target).lower()}"
        pdf_path = output_dir / f"compare_livetimes_multipanel_{suffix}.pdf"
        with PdfPages(pdf_path) as pdf:
            # Scale width with run count so adjacent vertical labels stay legible.
            save_multipanel(series, kin, target, "Run number", pdf,
                            max(14, 0.55 * len(runs)), rows, runs)
            for values, xlabel in ((beam_currents, "Mean beam current (uA)"),
                                   (s1x_rates, "HMS S1X rate (Hz)")):
                pairs = [(values[int(row["run"])], row) for row in rows
                         if int(row["run"]) in values]
                if not pairs:
                    print(f"[plot] No {xlabel} values for {kin} | {target}")
                    continue
                print(f"[plot] {kin} | {target} | {xlabel}: {len(pairs)}/{len(rows)} runs")
                pairs.sort(key=lambda pair: pair[0])
                plot_rows = [row for _, row in pairs]
                plot_x = [value for value, _ in pairs]
                save_multipanel(metric_series(plot_rows, plot_x), kin, target,
                                xlabel, pdf, 15, plot_rows)
        print(f"[plot] {pdf_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+", type=Path, help="Merged compare_livetimes_{kin}.csv files")
    parser.add_argument("--output-dir", required=True, type=Path, help="Plot directory")
    args = parser.parse_args()
    for path in args.csv:
        plot_csv(path, args.output_dir)


if __name__ == "__main__":
    main()
