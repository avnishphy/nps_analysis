#!/usr/bin/env python3
"""Make minimalist overlay plots from efficiency CSV files."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_METRICS = [
    ("HMS_pid_eff", "HMS PID"),
    ("HMS_cal_eff_tag_cer", "Cal tag Cer"),
    ("HMS_cer_eff_tag_cal", "Cer tag Cal"),
    ("HMS_tracking_eff", "Tracking"),
    ("HMS_hodo_3of4_eff", "Hodo 3/4"),
    ("NewGen_EDTM_livetime", "NewGen EDTM"),
]

COLORS = [
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#1f78b4",
    "#a6761d",
]

MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
PS_SHADE_COLORS = ["#f4f4f4", "#fff7e6", "#eef7ff", "#f4efff"]
PS_GROUP_COLORS = {
    "ps3": "#fff7e6",
    "ps4": "#eef7ff",
    "ps5": "#f4efff",
    "ps6": "#eef8ee",
}


def clean_cell(value: str | None) -> str:
    if value is None:
        return ""
    return value.strip().strip('"').strip()


def to_float(value: str | None) -> float:
    text = clean_cell(value)
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def to_int(value: str | None) -> int | None:
    x = to_float(value)
    if not math.isfinite(x):
        return None
    return int(round(x))


def finite(value: float) -> bool:
    return math.isfinite(value)


def parse_nominal_current_uA(kinematic: str) -> float:
    match = re.search(r"(?:^|_)x([0-9]+(?:p[0-9]+)?)(?:_|$)", kinematic)
    if not match:
        return math.nan
    return float(match.group(1).replace("p", "."))


def selection_report_path_for(efficiency_csv: Path) -> Path:
    name = efficiency_csv.name
    if name.startswith("efficiency_"):
        return efficiency_csv.with_name("selection_report_" + name[len("efficiency_") :])
    return efficiency_csv.with_name("selection_report_" + name)


def read_selection_mean_currents(path: Path) -> dict[int, float]:
    if not path.exists():
        return {}

    by_run: dict[int, list[float]] = defaultdict(list)
    with path.open("r", newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            return {}

        names = [clean_cell(h) for h in header]
        for raw in reader:
            if not raw or all(not clean_cell(v) for v in raw):
                continue
            row = {name: clean_cell(raw[i]) if i < len(raw) else "" for i, name in enumerate(names)}
            run = to_int(row.get("run_number"))
            current = to_float(row.get("mean_current_uA"))
            if run is not None and finite(current):
                by_run[run].append(current)

    return {run: mean(values) for run, values in by_run.items() if values}


def read_efficiency_csv(path: Path) -> list[dict[str, str]]:
    selection_currents = read_selection_mean_currents(selection_report_path_for(path))
    with path.open("r", newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            return []

        names = [clean_cell(h) for h in header]
        rows: list[dict[str, str]] = []
        for raw in reader:
            if not raw or all(not clean_cell(v) for v in raw):
                continue
            row = {name: clean_cell(raw[i]) if i < len(raw) else "" for i, name in enumerate(names)}
            row["_source_csv"] = str(path)
            run = to_int(row.get("run_number"))
            if run is not None and run in selection_currents:
                row["_current_uA"] = str(selection_currents[run])
            elif "kinematic_setting" in row:
                row["_current_uA"] = str(current_for_row(row))
            rows.append(row)
    return rows


def current_for_row(row: dict[str, str]) -> float:
    explicit = to_float(row.get("mean_current_uA"))
    if finite(explicit):
        return explicit

    nominal = parse_nominal_current_uA(clean_cell(row.get("kinematic_setting")))
    if finite(nominal):
        return nominal

    # Last-resort estimate only. In normal efficiency CSVs, beam_time may be
    # gated differently than charge, so kinematic-name current is preferred.
    charge = to_float(row.get("HEL_charge_after_cut_uC"))
    beam_time = to_float(row.get("beam_time"))
    if finite(charge) and finite(beam_time) and beam_time > 0.0:
        return charge / beam_time
    return math.nan


def discover_csvs(input_dir: Path) -> list[Path]:
    csvs = sorted(input_dir.glob("efficiency_*.csv"))
    return [p for p in csvs if p.name != "efficiency_runs_processed.csv"]


def selected_metrics(rows: Iterable[dict[str, str]], requested: str) -> list[tuple[str, str]]:
    available = set()
    for row in rows:
        available.update(row.keys())

    metrics = DEFAULT_METRICS
    if requested:
        wanted = [m.strip() for m in requested.split(",") if m.strip()]
        label_by_name = dict(DEFAULT_METRICS)
        metrics = [(name, label_by_name.get(name, name)) for name in wanted]

    return [(name, label) for name, label in metrics if name in available]


def prescale_group(row: dict[str, str]) -> str:
    token = clean_cell(row.get("prescale_token"))
    token_match = re.search(r"\b(ps\d+)\b", token, flags=re.IGNORECASE)
    if token_match:
        return token_match.group(1).lower()

    trig = clean_cell(row.get("which_TRIG"))
    trig_match = re.search(r"hTRIG(\d+)", trig, flags=re.IGNORECASE)
    if trig_match:
        return f"ps{trig_match.group(1)}"

    ps_factor = clean_cell(row.get("ps_factor"))
    if ps_factor:
        return f"ps_factor {ps_factor}"

    return "unknown"


def run_positions(rows: list[dict[str, str]]) -> dict[int, int]:
    runs = sorted({run for row in rows if (run := to_int(row.get("run_number"))) is not None})
    return {run: idx for idx, run in enumerate(runs)}


def compact_run_ticks(runs: list[int], max_ticks: int = 18) -> tuple[list[int], list[str]]:
    if not runs:
        return [], []

    if len(runs) <= max_ticks:
        positions = list(range(len(runs)))
    else:
        step = max(1, math.ceil(len(runs) / max_ticks))
        positions = list(range(0, len(runs), step))
        if positions[-1] != len(runs) - 1:
            positions.append(len(runs) - 1)

    return positions, [str(runs[pos]) for pos in positions]


def large_run_gap_boundaries(runs: list[int]) -> list[tuple[float, int, int]]:
    if len(runs) < 2:
        return []

    gaps = [runs[i + 1] - runs[i] for i in range(len(runs) - 1)]
    positive_gaps = sorted(g for g in gaps if g > 0)
    median_gap = positive_gaps[len(positive_gaps) // 2] if positive_gaps else 1
    threshold = max(10, 5 * median_gap)

    boundaries: list[tuple[float, int, int]] = []
    for idx, gap in enumerate(gaps):
        if gap > threshold:
            boundaries.append((idx + 0.5, runs[idx], runs[idx + 1]))
    return boundaries


def run_blocks(runs: list[int]) -> list[tuple[int, int]]:
    if not runs:
        return []

    boundaries = large_run_gap_boundaries(runs)
    if not boundaries:
        return [(runs[0], runs[-1])]

    blocks: list[tuple[int, int]] = []
    start_idx = 0
    for _, left, right in boundaries:
        left_idx = runs.index(left)
        right_idx = runs.index(right)
        blocks.append((runs[start_idx], runs[left_idx]))
        start_idx = right_idx
    blocks.append((runs[start_idx], runs[-1]))
    return blocks


def format_run_blocks(runs: list[int], max_per_line: int = 3) -> str:
    blocks = [f"{lo}-{hi}" if lo != hi else str(lo) for lo, hi in run_blocks(runs)]
    lines = ["Run ranges"]
    for idx in range(0, len(blocks), max_per_line):
        lines.append("; ".join(blocks[idx : idx + max_per_line]))
    return "\n".join(lines)


def setup_compact_run_axis(
    ax: plt.Axes,
    runs: list[int],
    xlabel: bool,
    show_ticklabels: bool = True,
) -> None:
    positions, labels = compact_run_ticks(runs)
    ax.set_xlim(-0.8, len(runs) - 0.2 if runs else 0.2)
    ax.set_xticks(positions)
    if show_ticklabels:
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    else:
        ax.set_xticklabels([])
    ax.set_xlabel("Run number (compact axis; large gaps compressed)" if xlabel else "")

    for x, left, right in large_run_gap_boundaries(runs):
        ax.axvline(x, color="#555555", linestyle="--", linewidth=0.85, alpha=0.5)
        if xlabel:
            ymin, ymax = ax.get_ylim()
            ax.text(
                x,
                ymin + 0.03 * (ymax - ymin),
                f"{left}|{right}",
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=6.5,
                color="#555555",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 0.3},
            )


def percentile(values: list[float], pct: float) -> float:
    vals = sorted(v for v in values if finite(v))
    if not vals:
        return math.nan
    if len(vals) == 1:
        return vals[0]

    pos = (len(vals) - 1) * pct
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return vals[int(pos)]
    frac = pos - lo
    return vals[lo] * (1.0 - frac) + vals[hi] * frac


def robust_y_limits(points: list[tuple[int, float, float]]) -> tuple[float, float, set[int]]:
    values = [y for _, y, _ in points if finite(y)]
    if not values:
        return math.nan, math.nan, set()
    if len(values) < 6:
        ymin = min(values)
        ymax = max(values)
        pad = max(0.005, 0.08 * max(ymax - ymin, 0.01))
        return ymin - pad, ymax + pad, set()

    q1 = percentile(values, 0.25)
    q3 = percentile(values, 0.75)
    iqr = max(q3 - q1, 1.0e-9)
    # Be conservative: most efficiency structure is real run-to-run behavior.
    # The clipping is only meant to keep pathological near-zero failures from
    # flattening otherwise stable panels.
    fence_low = q1 - max(5.0 * iqr, 0.08)

    inliers = [v for v in values if v >= fence_low]
    if len(inliers) < max(4, len(values) // 2):
        inliers = values

    ymin = min(inliers)
    ymax = max(values)
    pad = max(0.002, 0.12 * max(ymax - ymin, 0.01))
    low = ymin - pad
    high = ymax + pad

    outlier_x = {x for x, y, _ in points if finite(y) and y < low}
    return low, high, outlier_x


def apply_robust_y_limits(
    ax: plt.Axes,
    series: list[tuple[str, list[tuple[int, float, float]]]],
    run_by_x: dict[int, int],
    annotate: bool,
) -> None:
    all_points = [point for _, points in series for point in points]
    low, high, _ = robust_y_limits(all_points)
    if not finite(low) or not finite(high) or low >= high:
        return

    ax.set_ylim(low, high)
    if not annotate:
        return

    clipped: list[tuple[int, str, float, str]] = []
    for label, points in series:
        for x, y, _ in points:
            if y < low:
                clipped.append((x, label, y, "low"))

    if not clipped:
        return

    clipped_by_run: dict[float | int, list[tuple[str, float, str]]] = defaultdict(list)
    for x, label, y, side in clipped:
        clipped_by_run[run_by_x.get(x, x)].append((label, y, side))

    for x, _, y, side in clipped:
        edge_y = low if side == "low" else high
        marker = "v" if side == "low" else "^"
        ax.scatter([x], [edge_y], marker=marker, s=28, color="#202020", zorder=6, clip_on=False)

    run_notes = []
    for run in sorted(clipped_by_run):
        values = clipped_by_run[run]
        axis_label = format_axis_value(run)
        if len(values) == 1:
            label, y, _ = values[0]
            run_notes.append(f"{axis_label}: {label}={y:.3g}")
        else:
            run_notes.append(f"{axis_label}: {len(values)} low outliers")

    max_notes = 5
    note_text = "Low outliers clipped\n" + "\n".join(run_notes[:max_notes])
    if len(run_notes) > max_notes:
        note_text += f"\n+{len(run_notes) - max_notes} more"

    ax.text(
        0.012,
        0.02,
        note_text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.1,
        color="#202020",
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.92},
    )


def format_axis_value(value: float | int) -> str:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return str(value)
    if abs(x - round(x)) < 1.0e-6:
        return str(int(round(x)))
    return f"{x:.2f}"


def annotate_run_ranges(ax: plt.Axes, runs: list[int]) -> None:
    ax.text(
        1.015,
        0.70,
        format_run_blocks(runs),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
        linespacing=1.18,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.94},
    )


def annotate_prescale_blocks(
    ax: plt.Axes,
    rows: list[dict[str, str]],
    x_by_run: dict[int, int],
    label_blocks: bool = False,
    show_key: bool = False,
) -> bool:
    run_states = [
        (x_by_run[run], prescale_group(row))
        for row in rows
        if (run := to_int(row.get("run_number"))) is not None and run in x_by_run
    ]
    if not run_states:
        return False

    run_states.sort(key=lambda item: item[0])
    groups = sorted({state for _, state in run_states})

    blocks: list[tuple[int, int, str]] = []
    start_run, last_run, current_state = run_states[0][0], run_states[0][0], run_states[0][1]
    for run, state in run_states[1:]:
        if state != current_state:
            blocks.append((start_run, last_run, current_state))
            ax.axvline((last_run + run) / 2.0, color="#777777", linestyle=":", linewidth=0.8, alpha=0.55)
            start_run = run
            current_state = state
        last_run = run
    blocks.append((start_run, last_run, current_state))

    ymin, ymax = ax.get_ylim()
    span = ymax - ymin if ymax > ymin else 1.0
    if label_blocks:
        ax.set_ylim(ymin, ymax + 0.16 * span)
    label_y0 = ymax + 0.035 * span

    for idx, (lo, hi, state) in enumerate(blocks):
        pad = 0.45 if lo != hi else 0.35
        ax.axvspan(
            lo - pad,
            hi + pad,
            color=PS_GROUP_COLORS.get(state, PS_SHADE_COLORS[idx % len(PS_SHADE_COLORS)]),
            alpha=0.45,
            zorder=-2,
        )
        if not label_blocks:
            continue
        ax.text(
            (lo + hi) / 2.0,
            label_y0 + (idx % 2) * 0.045 * span,
            state,
            ha="center",
            va="bottom",
            fontsize=8.2,
            color="#303030",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "#bbbbbb", "alpha": 0.9},
        )

    if not show_key:
        return True

    key_lines = ["Prescale groups", *groups]
    ax.text(
        1.015,
        0.98,
        "\n".join(key_lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.4,
        linespacing=1.18,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.94},
    )
    return True


def setup_axes(ax: plt.Axes, title: str, xlabel: str) -> None:
    ax.set_title(title, fontsize=13, pad=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Efficiency / livetime")
    ax.grid(True, color="#d8d8d8", linewidth=0.6, alpha=0.8)
    ax.tick_params(direction="out", length=4, width=0.8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def save_figure(fig: plt.Figure, stem: Path, formats: list[str], dpi: int) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        out = stem.with_suffix(f".{fmt}")
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"Wrote {out}")


def metric_points(
    rows: list[dict[str, str]],
    metric_name: str,
    x_by_run: dict[int, int],
) -> list[tuple[int, float, float]]:
    points: list[tuple[int, float, float]] = []
    for row in rows:
        run = to_int(row.get("run_number"))
        if run is None or run not in x_by_run:
            continue
        value = to_float(row.get(metric_name))
        if not finite(value):
            continue
        err = to_float(row.get(f"{metric_name}_err"))
        points.append((x_by_run[run], value, err if finite(err) and err >= 0.0 else 0.0))
    return points


def draw_metric_series(
    ax: plt.Axes,
    points: list[tuple[int, float, float]],
    label: str,
    color: str,
    marker: str,
    use_errorbars: bool,
    linewidth: float = 1.25,
    markersize: float = 4.4,
) -> None:
    if not points:
        return

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    yerr = [p[2] for p in points]
    if use_errorbars and any(e > 0.0 for e in yerr):
        ax.errorbar(
            xs,
            ys,
            yerr=yerr,
            label=label,
            color=color,
            marker=marker,
            markersize=markersize,
            linewidth=linewidth,
            capsize=2.0,
            alpha=0.95,
        )
    else:
        ax.plot(xs, ys, label=label, color=color, marker=marker, markersize=markersize, linewidth=linewidth)


def plot_vs_run(
    rows: list[dict[str, str]],
    metrics: list[tuple[str, str]],
    out_dir: Path,
    formats: list[str],
    dpi: int,
    use_errorbars: bool,
) -> None:
    rows = sorted(rows, key=lambda r: (to_int(r.get("run_number")) or -1))
    if not rows:
        return

    kin = clean_cell(rows[0].get("kinematic_setting")) or Path(rows[0]["_source_csv"]).stem
    x_by_run = run_positions(rows)
    runs = [run for run, _ in sorted(x_by_run.items(), key=lambda item: item[1])]
    if not runs:
        return
    run_by_x = {x: run for run, x in x_by_run.items()}

    ncols = 2
    metric_rows = math.ceil(len(metrics) / ncols)
    fig_height = 4.3 + 2.35 * metric_rows
    fig = plt.figure(figsize=(13.2, fig_height))
    grid = fig.add_gridspec(
        metric_rows + 1,
        ncols,
        height_ratios=[1.45] + [1.0] * metric_rows,
        hspace=0.34,
        wspace=0.18,
    )

    overlay_ax = fig.add_subplot(grid[0, :])
    setup_axes(overlay_ax, f"{kin}: efficiencies vs run", "")

    overlay_series: list[tuple[str, list[tuple[int, float, float]]]] = []
    for idx, (name, label) in enumerate(metrics):
        points = metric_points(rows, name, x_by_run)
        overlay_series.append((label, points))
        color = COLORS[idx % len(COLORS)]
        marker = MARKERS[idx % len(MARKERS)]
        draw_metric_series(overlay_ax, points, label, color, marker, use_errorbars)

    apply_robust_y_limits(overlay_ax, overlay_series, run_by_x, annotate=True)
    has_prescale_key = annotate_prescale_blocks(
        overlay_ax,
        rows,
        x_by_run,
        label_blocks=True,
        show_key=True,
    )
    annotate_run_ranges(overlay_ax, runs)
    setup_compact_run_axis(overlay_ax, runs, xlabel=False, show_ticklabels=True)
    overlay_ax.legend(frameon=False, fontsize=8.8, ncols=3, loc="lower center")

    for idx, (name, label) in enumerate(metrics):
        panel_ax = fig.add_subplot(grid[idx // ncols + 1, idx % ncols])
        setup_axes(panel_ax, label, "")
        points = metric_points(rows, name, x_by_run)
        color = COLORS[idx % len(COLORS)]
        marker = MARKERS[idx % len(MARKERS)]
        draw_metric_series(
            panel_ax,
            points,
            label,
            color,
            marker,
            use_errorbars,
            linewidth=1.15,
            markersize=3.9,
        )
        apply_robust_y_limits(panel_ax, [(label, points)], run_by_x, annotate=True)
        annotate_prescale_blocks(panel_ax, rows, x_by_run)
        show_bottom_axis = idx // ncols == metric_rows - 1
        setup_compact_run_axis(
            panel_ax,
            runs,
            xlabel=show_bottom_axis,
            show_ticklabels=show_bottom_axis,
        )

    total_slots = metric_rows * ncols
    for empty_idx in range(len(metrics), total_slots):
        empty_ax = fig.add_subplot(grid[empty_idx // ncols + 1, empty_idx % ncols])
        empty_ax.axis("off")

    fig.subplots_adjust(
        left=0.07,
        right=0.84 if has_prescale_key else 0.98,
        top=0.96,
        bottom=0.07,
        hspace=0.58,
        wspace=0.20,
    )
    safe_kin = re.sub(r"[^A-Za-z0-9_.-]+", "_", kin)
    save_figure(fig, out_dir / f"efficiency_multipanel_vs_run_{safe_kin}", formats, dpi)
    plt.close(fig)


def sem(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return stdev(values) / math.sqrt(len(values))


def metric_points_by_current(
    rows: list[dict[str, str]],
    metric_name: str,
) -> list[tuple[float, float, float]]:
    by_current: dict[float, list[float]] = defaultdict(list)
    for row in rows:
        current = to_float(row.get("_current_uA"))
        value = to_float(row.get(metric_name))
        if finite(current) and finite(value):
            by_current[current].append(value)

    points: list[tuple[float, float, float]] = []
    for current in sorted(by_current):
        values = by_current[current]
        points.append((current, mean(values), sem(values)))
    return points


def setup_current_axis(ax: plt.Axes, xlabel: bool) -> None:
    ax.set_xlabel("Current (uA)" if xlabel else "")


def plot_vs_current_multipanel(
    rows: list[dict[str, str]],
    metrics: list[tuple[str, str]],
    out_dir: Path,
    formats: list[str],
    dpi: int,
    use_errorbars: bool,
    title: str,
    output_name: str,
) -> None:
    if not rows:
        return

    ncols = 2
    metric_rows = math.ceil(len(metrics) / ncols)
    fig_height = 4.2 + 2.35 * metric_rows
    fig = plt.figure(figsize=(12.4, fig_height))
    grid = fig.add_gridspec(
        metric_rows + 1,
        ncols,
        height_ratios=[1.35] + [1.0] * metric_rows,
        hspace=0.58,
        wspace=0.20,
    )

    overlay_ax = fig.add_subplot(grid[0, :])
    setup_axes(overlay_ax, f"{title}: efficiencies vs current", "")

    overlay_series: list[tuple[str, list[tuple[float, float, float]]]] = []
    current_by_x: dict[float, float] = {}
    for idx, (name, label) in enumerate(metrics):
        points = metric_points_by_current(rows, name)
        overlay_series.append((label, points))
        current_by_x.update({x: x for x, _, _ in points})
        color = COLORS[idx % len(COLORS)]
        marker = MARKERS[idx % len(MARKERS)]
        draw_metric_series(overlay_ax, points, label, color, marker, use_errorbars, markersize=5.0)

    apply_robust_y_limits(overlay_ax, overlay_series, current_by_x, annotate=False)
    setup_current_axis(overlay_ax, xlabel=False)
    overlay_ax.legend(frameon=False, fontsize=8.8, ncols=3, loc="best")

    for idx, (name, label) in enumerate(metrics):
        panel_ax = fig.add_subplot(grid[idx // ncols + 1, idx % ncols])
        setup_axes(panel_ax, label, "")
        points = metric_points_by_current(rows, name)
        current_by_x = {x: x for x, _, _ in points}
        color = COLORS[idx % len(COLORS)]
        marker = MARKERS[idx % len(MARKERS)]
        draw_metric_series(
            panel_ax,
            points,
            label,
            color,
            marker,
            use_errorbars,
            linewidth=1.15,
            markersize=4.5,
        )
        apply_robust_y_limits(panel_ax, [(label, points)], current_by_x, annotate=True)
        setup_current_axis(panel_ax, xlabel=(idx // ncols == metric_rows - 1))

    total_slots = metric_rows * ncols
    for empty_idx in range(len(metrics), total_slots):
        empty_ax = fig.add_subplot(grid[empty_idx // ncols + 1, empty_idx % ncols])
        empty_ax.axis("off")

    fig.subplots_adjust(left=0.07, right=0.98, top=0.96, bottom=0.07, hspace=0.58, wspace=0.20)
    save_figure(fig, out_dir / output_name, formats, dpi)
    plt.close(fig)


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    default_input_dir = repo_root / "output" / "efficiency_stuff"

    parser = argparse.ArgumentParser(
        description="Create high-quality overlay plots from efficiency CSV files."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_input_dir,
        help="Directory containing efficiency_*.csv files.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        nargs="*",
        default=None,
        help="Explicit efficiency CSV files to plot. Overrides --input-dir discovery.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for plots. Default: <input-dir>/plots.",
    )
    parser.add_argument(
        "--formats",
        default="png,pdf",
        help="Comma-separated output formats, e.g. png,pdf.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Raster output DPI.")
    parser.add_argument(
        "--metrics",
        default="",
        help="Comma-separated metric column names. Default: standard HMS/NewGen metrics.",
    )
    parser.add_argument(
        "--no-errorbars",
        action="store_true",
        help="Draw plain line/marker overlays without error bars.",
    )
    args = parser.parse_args()

    csvs = [p.resolve() for p in args.csv] if args.csv else discover_csvs(args.input_dir.resolve())
    if not csvs:
        raise FileNotFoundError(f"No efficiency CSV files found in {args.input_dir}")

    formats = [fmt.strip().lstrip(".") for fmt in args.formats.split(",") if fmt.strip()]
    out_dir = (args.out_dir.resolve() if args.out_dir else args.input_dir.resolve() / "plots")

    rows_by_csv: dict[Path, list[dict[str, str]]] = {}
    all_rows: list[dict[str, str]] = []
    for csv_path in csvs:
        rows = read_efficiency_csv(csv_path)
        if not rows:
            print(f"Skipping empty CSV: {csv_path}")
            continue
        rows_by_csv[csv_path] = rows
        all_rows.extend(rows)

    metrics = selected_metrics(all_rows, args.metrics)
    if not metrics:
        raise RuntimeError("None of the requested/default efficiency metrics were found.")

    for csv_path, rows in rows_by_csv.items():
        plot_vs_run(rows, metrics, out_dir, formats, args.dpi, not args.no_errorbars)
        kin = clean_cell(rows[0].get("kinematic_setting")) or csv_path.stem.replace("efficiency_", "")
        safe_kin = re.sub(r"[^A-Za-z0-9_.-]+", "_", kin)
        plot_vs_current_multipanel(
            rows,
            metrics,
            out_dir,
            formats,
            args.dpi,
            not args.no_errorbars,
            kin,
            f"efficiency_multipanel_vs_current_{safe_kin}",
        )

    plot_vs_current_multipanel(
        all_rows,
        metrics,
        out_dir,
        formats,
        args.dpi,
        not args.no_errorbars,
        "All kinematics",
        "efficiency_multipanel_vs_current_all",
    )
    print(f"Plotted {len(metrics)} metrics from {len(rows_by_csv)} CSV file(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
