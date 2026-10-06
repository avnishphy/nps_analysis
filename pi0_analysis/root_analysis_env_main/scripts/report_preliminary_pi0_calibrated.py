#!/usr/bin/env python3
"""Calibrate M0 statistical intervals from a matching toy campaign.

This is a reporting-only step. It never runs production, generates toys, or
changes the constrained central extraction. The selected JSON binning must
match the toy and curvature campaigns. The output directory must be new.
"""

import argparse
import csv
import hashlib
import json
import shutil
import tempfile
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
COMPONENTS = ("U", "LT", "TT")
UNITS = "nb/GeV^2"
BINNING_KEYS = (
    "configured_kinematic",
    "phi_bins",
    "tprime_bin_edges",
    "q2_bin_edges",
    "xb_bin_edges_by_q2",
    "diamond_xb_q2_vertices",
)
N_TPRIME_BINS = 0
PHI_BINS = 0
NAMES = []


def configure_binning(configuration):
    """Set report dimensions from the toy campaign's snapshotted JSON."""
    global N_TPRIME_BINS, PHI_BINS, NAMES
    edges = np.asarray(configuration.get("tprime_bin_edges", []), dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all():
        raise ValueError("toy source has invalid tprime_bin_edges")
    if not np.all(np.diff(edges) > 0.0):
        raise ValueError("toy source tprime_bin_edges are not strictly increasing")
    PHI_BINS = int(configuration.get("phi_bins", 0))
    if PHI_BINS <= 0:
        raise ValueError("toy source phi_bins must be positive")
    q2_edges = configuration.get("q2_bin_edges", [])
    xb_rows = configuration.get("xb_bin_edges_by_q2", [])
    if len(q2_edges) != 2 or len(xb_rows) != 1 or len(xb_rows[0]) != 2:
        raise ValueError(
            "calibrated M0 reporting currently requires one Q2/xB group; "
            "regenerate and extend the campaign tooling before publishing more groups"
        )
    N_TPRIME_BINS = len(edges) - 1
    NAMES = [
        f"{component}(bin{index})"
        for component in COMPONENTS
        for index in range(N_TPRIME_BINS)
    ]
    return edges


def validate_selected_config(selected, source_configuration):
    mismatches = [
        key
        for key in BINNING_KEYS
        if selected.get(key) != source_configuration.get(key)
    ]
    if mismatches:
        selected_edges = selected.get("tprime_bin_edges")
        source_edges = source_configuration.get("tprime_bin_edges")
        raise ValueError(
            "selected xsec JSON does not match the central/toy campaign "
            f"({', '.join(mismatches)} differ): selected tprime_bin_edges="
            f"{selected_edges}, toy-source tprime_bin_edges={source_edges}. "
            "Do not publish these toys for this fit; regenerate the central/toy "
            "and curvature campaigns with the selected JSON, then pass those "
            "directories via --calibration-source/--curvature-source."
        )


def atomic_text(path, text):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


def write_rows(path, rows):
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def write_matrix(path, matrix):
    write_rows(
        path,
        [
            {"quantity": row_name, **{name: float(value) for name, value in zip(NAMES, row)}}
            for row_name, row in zip(NAMES, matrix)
        ],
    )


def correlation(covariance):
    scale = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denominator = np.outer(scale, scale)
    return np.divide(
        covariance,
        denominator,
        out=np.zeros_like(covariance),
        where=denominator > 0.0,
    )


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def quantile(values, probability, axis=None):
    return np.quantile(values, probability, axis=axis, method="linear")


def read_curvature_table(path):
    rows = list(csv.DictReader(path.open()))
    if len(rows) != N_TPRIME_BINS:
        raise ValueError(
            f"expected {N_TPRIME_BINS} t-prime bins in {path}, found {len(rows)}"
        )
    central = np.array(
        [float(row[column]) for column in ("sigma_U", "sigma_LT", "sigma_TT") for row in rows]
    )
    errors = np.array(
        [float(row[column]) for column in ("U_stat", "LT_stat", "TT_stat") for row in rows]
    )
    return rows, central, errors


def read_matrix(path):
    rows = list(csv.DictReader(path.open()))
    if [row["parameter"] for row in rows] != NAMES:
        raise ValueError(f"unexpected matrix row ordering in {path}")
    matrix = np.array([[float(row[name]) for name in NAMES] for row in rows])
    expected = len(NAMES)
    if matrix.shape != (expected, expected) or not np.isfinite(matrix).all():
        raise ValueError(f"invalid {expected}x{expected} matrix in {path}")
    return matrix


def classify_interval(q16, q84):
    lower_width = q84
    upper_width = -q16
    contains_estimate = q16 <= 0.0 <= q84
    if not contains_estimate:
        return True, "basic interval excludes central estimate because residual bias exceeds one tail"
    if lower_width <= 0.0 or upper_width <= 0.0:
        return True, "basic interval is one-sided relative to central estimate"
    ratio = max(lower_width, upper_width) / min(lower_width, upper_width)
    if ratio >= 3.0:
        return True, f"basic lower/upper half-width ratio is {ratio:.3f}"
    return False, "symmetric delta68 used; basic interval is not extremely asymmetric"


def validate_inputs(
    source,
    source_configuration,
    tprime_edges,
    toys,
    toy_summary,
    central,
    curve_rows,
    curve_central,
    curve_errors,
    c_curvature,
):
    expected_quantities = len(NAMES)
    if toy_summary.get("accepted") != 500 or toys["pub"].shape != (500, expected_quantities):
        raise ValueError(
            "the accepted local M0 ensemble is not exactly "
            f"500 x {expected_quantities} for the selected binning"
        )
    toy_ids = np.asarray(toys["index"], dtype=int)
    if len(np.unique(toy_ids)) != 500:
        raise ValueError("toy IDs are not unique")
    if not np.isfinite(toys["pub"]).all():
        raise ValueError("toy published values contain non-finite entries")
    if central.get("code") != 0:
        raise ValueError("frozen additive M0 central fit is not marked converged")
    if np.asarray(central.get("pub", [])).shape != (expected_quantities,):
        raise ValueError("central published-vector size does not match selected binning")
    np.testing.assert_allclose(central["pub"], curve_central, rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(np.sqrt(np.diag(c_curvature)), curve_errors, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(c_curvature, c_curvature.T, rtol=0.0, atol=2e-12)
    if any(row["units"] != UNITS for row in curve_rows):
        raise ValueError("unexpected units in curvature cross-section table")
    configured_kinematic = source_configuration.get("configured_kinematic")
    if not isinstance(configured_kinematic, str) or not configured_kinematic.strip():
        raise ValueError("source campaign has no configured kinematic setting")
    for index, row in enumerate(curve_rows):
        np.testing.assert_allclose(
            [float(row["tprime_low"]), float(row["tprime_high"])],
            tprime_edges[index : index + 2],
            rtol=0.0,
            atol=5e-12,
            err_msg=f"curvature t-prime bin {index} disagrees with toy config snapshot",
        )
    generator = (source / "before/preliminary_pi0_xsec_v7.py").read_text()
    required_fragments = (
        "central=json.loads((out/'central_results.json').read_text())['additive']",
        "y=y-y0+np.array(central['prediction'])",
        "result=model.constrained(y,v,start,mult)",
    )
    if any(fragment not in generator for fragment in required_fragments):
        raise ValueError("archived toy generator does not show the expected frozen-M0 construction")
    return toy_ids


def calibrate(residuals, toy_ids):
    bias = residuals.mean(axis=0)
    toy_sd = residuals.std(axis=0, ddof=1)
    rmse = np.sqrt(np.mean(residuals * residuals, axis=0))
    q16, q50, q84 = quantile(residuals, [0.16, 0.50, 0.84], axis=0)
    delta68 = quantile(np.abs(residuals), 0.68, axis=0)

    folds = toy_ids % 5
    decisions = np.zeros_like(residuals, dtype=bool)
    fold_rows = []
    for quantity_index, name in enumerate(NAMES):
        for fold in range(5):
            validation = folds == fold
            training = ~validation
            radius = float(quantile(np.abs(residuals[training, quantity_index]), 0.68))
            covered = np.abs(residuals[validation, quantity_index]) <= radius
            decisions[validation, quantity_index] = covered
            fold_coverage = float(covered.mean())
            fold_rows.append(
                {
                    "quantity": name,
                    "fold": fold,
                    "fold_rule": "toy_id modulo 5",
                    "training_n": int(training.sum()),
                    "validation_n": int(validation.sum()),
                    "delta68_train": radius,
                    "covered_n": int(covered.sum()),
                    "coverage68": fold_coverage,
                    "binomial_standard_error": float(
                        np.sqrt(fold_coverage * (1.0 - fold_coverage) / covered.size)
                    ),
                }
            )

    coverage = decisions.mean(axis=0)
    coverage_se = np.sqrt(coverage * (1.0 - coverage) / residuals.shape[0])
    aggregate_rows = [
        {
            "quantity": name,
            "fold": "aggregate",
            "fold_rule": "each toy held out exactly once by toy_id modulo 5",
            "training_n": "400 approximately; see fold rows",
            "validation_n": residuals.shape[0],
            "delta68_train": "fold-specific; see fold rows",
            "covered_n": int(decisions[:, index].sum()),
            "coverage68": float(coverage[index]),
            "binomial_standard_error": float(coverage_se[index]),
        }
        for index, name in enumerate(NAMES)
    ]
    return {
        "bias": bias,
        "toy_sd": toy_sd,
        "rmse": rmse,
        "q16": q16,
        "q50": q50,
        "q84": q84,
        "delta68": delta68,
        "coverage": coverage,
        "coverage_se": coverage_se,
        "coverage_rows": aggregate_rows + fold_rows,
    }


def summarize_bin_statistics(central, curve_rows, truth, stats):
    """Summarize the frozen weighted data vector; this is not an event count."""
    weighted_yields = np.asarray(central["y"], dtype=float)
    weighted_variances = np.asarray(central["data_variance"], dtype=float)
    expected_rows = N_TPRIME_BINS * PHI_BINS
    if weighted_yields.shape != (expected_rows,) or weighted_variances.shape != (expected_rows,):
        raise ValueError(
            f"expected {expected_rows}-row frozen data and variance vectors "
            "for the selected t-prime/phi binning"
        )
    if not np.isfinite(weighted_yields).all() or not np.isfinite(weighted_variances).all():
        raise ValueError("non-finite frozen data or variance entry")
    if np.any(weighted_variances < 0.0):
        raise ValueError("negative frozen data variance entry")

    rows = []
    for bin_index, source_row in enumerate(curve_rows):
        selection = slice(PHI_BINS * bin_index, PHI_BINS * (bin_index + 1))
        values = weighted_yields[selection]
        variances = weighted_variances[selection]
        included = variances > 0.0
        total = float(values[included].sum())
        total_variance = float(variances[included].sum())
        if total <= 0.0 or total_variance <= 0.0:
            raise ValueError(f"non-positive weighted statistics in t-prime bin {bin_index}")
        rows.append(
            {
                "tprime_bin": bin_index,
                "tprime_low": float(source_row["tprime_low"]),
                "tprime_high": float(source_row["tprime_high"]),
                "response_weighted_tprime": float(source_row["response_weighted_tprime"]),
                "included_phi_rows": int(included.sum()),
                "weighted_yield": total,
                "weighted_variance": total_variance,
                "effective_entries": total * total / total_variance,
                "relative_weighted_yield_error": np.sqrt(total_variance) / total,
                "U_delta68_stat": float(stats["delta68"][bin_index]),
                "LT_delta68_stat": float(stats["delta68"][N_TPRIME_BINS + bin_index]),
                "TT_delta68_stat": float(stats["delta68"][2 * N_TPRIME_BINS + bin_index]),
                "U_relative_delta68": float(stats["delta68"][bin_index] / abs(truth[bin_index])),
                "LT_relative_delta68": float(
                    stats["delta68"][N_TPRIME_BINS + bin_index]
                    / abs(truth[N_TPRIME_BINS + bin_index])
                ),
                "TT_relative_delta68": float(
                    stats["delta68"][2 * N_TPRIME_BINS + bin_index]
                    / abs(truth[2 * N_TPRIME_BINS + bin_index])
                ),
            }
        )
    return rows


def make_component_plot(stage, component, means, central, stats, details, kinematic):
    component_index = COMPONENTS.index(component)
    offset = N_TPRIME_BINS * component_index
    indices = np.arange(offset, offset + N_TPRIME_BINS)
    x = -means
    values = central[indices]

    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    axis.plot(x, values, "o", color="#12263a", ms=5.5, zorder=5)
    cap_width = 0.009
    used_symmetric = False
    used_basic = False
    for point, index in enumerate(indices):
        row = details[index]
        if row["use_basic_interval"]:
            low = row["basic_interval_lower"]
            high = row["basic_interval_upper"]
            axis.vlines(x[point], values[point] - stats["delta68"][index], values[point] + stats["delta68"][index], color="0.72", lw=1.2, zorder=1)
            axis.hlines(
                [values[point] - stats["delta68"][index], values[point] + stats["delta68"][index]],
                x[point] - 0.6 * cap_width,
                x[point] + 0.6 * cap_width,
                color="0.72",
                lw=1.2,
                zorder=1,
            )
            axis.vlines(x[point], low, high, color="#1261a0", lw=2.2, zorder=3)
            axis.hlines([low, high], x[point] - cap_width, x[point] + cap_width, color="#1261a0", lw=2.2, zorder=3)
            used_basic = True
        else:
            low = values[point] - stats["delta68"][index]
            high = values[point] + stats["delta68"][index]
            axis.vlines(x[point], low, high, color="#1261a0", lw=2.2, zorder=3)
            axis.hlines([low, high], x[point] - cap_width, x[point] + cap_width, color="#1261a0", lw=2.2, zorder=3)
            used_symmetric = True

    axis.axhline(0.0, color="0.55", lw=0.8)
    axis.set_xlim(0.0, max(0.05, 1.05 * float(np.max(-means))))
    axis.set_xlabel(r"$-t'\ [\mathrm{GeV}^2]$")
    axis.set_ylabel(rf"$\sigma_{{{component}}}\ [\mathrm{{nb}}/\mathrm{{GeV}}^2]$")
    axis.set_title(
        f"PRELIMINARY - {kinematic}\n"
        "model-dependent forward-folded extraction\n"
        "toy-calibrated local 68% statistical uncertainties"
    )
    handles = [Line2D([0], [0], marker="o", linestyle="none", color="#12263a", label="frozen constrained M0 central")]
    if used_symmetric:
        handles.append(Line2D([0], [0], color="#1261a0", lw=2.2, label="symmetric calibrated +/- delta68"))
    if used_basic:
        handles.append(Line2D([0], [0], color="#1261a0", lw=2.2, label="bias-aware basic 16-84% interval"))
        handles.append(Line2D([0], [0], color="0.72", lw=1.2, label="symmetric +/- delta68 reference"))
    axis.legend(handles=handles, fontsize=7.8, loc="best")
    axis.text(0.02, 0.02, "Intervals switch to the basic construction only where residual asymmetry is extreme.", transform=axis.transAxes, fontsize=7.2)
    figure.tight_layout()
    for extension in ("png", "pdf"):
        figure.savefig(stage / f"sigma_{component}_preliminary_calibrated.{extension}", dpi=200)
    plt.close(figure)


def make_fit_diagnostics_scope_page(path):
    figure, axis = plt.subplots(figsize=(11, 8.5))
    axis.axis("off")
    axis.text(0.05, 0.91, "Appendix: complete SIMC-model fit diagnostics", fontsize=19, weight="bold", va="top")
    axis.text(
        0.05,
        0.82,
        textwrap.fill(
            "The following pages are the complete diagnostic book produced by the current "
            "generic upstream pi0_weight SIMC-model fit. They document inputs, missing-mass "
            "selection, detector-level comparisons, response coverage, fit quality, model "
            "assumptions, parameter behavior, and reconstructed-yield closure.",
            105,
        ),
        fontsize=12,
        va="top",
        linespacing=1.45,
    )
    warning = (
        "Estimator boundary: the preliminary cross sections and toy-calibrated intervals in "
        "the preceding section belong to the separately frozen constrained additive-weight M0 "
        "estimator. They do not apply to the generic-fit summary or slice CSVs shown in this "
        "appendix. The appendix is included for extraction diagnostics and provenance, not as "
        "a second calibrated result."
    )
    axis.text(
        0.05,
        0.57,
        textwrap.fill(warning, 100),
        fontsize=12,
        va="top",
        linespacing=1.5,
        bbox={"boxstyle": "round,pad=0.8", "facecolor": "#fff2cc", "edgecolor": "#b08b2e"},
    )
    axis.text(
        0.05,
        0.24,
        "Reading order",
        fontsize=13,
        weight="bold",
        va="top",
    )
    axis.text(
        0.07,
        0.18,
        "1. Calibrated release decision and uncertainty audit\n"
        "2. Toy coverage, bias, and covariance diagnostics\n"
        "3. Preliminary U, LT, and TT cross sections\n"
        "4. Complete generic-fit diagnostic appendix (following pages)",
        fontsize=11,
        va="top",
        linespacing=1.5,
    )
    figure.savefig(path)
    plt.close(figure)


def diagnostics_pdf(
    path,
    curvature_errors,
    stats,
    c_curvature,
    c_toy,
    r_toy,
    c68,
    details,
    summary_text,
    bin_statistics,
):
    with PdfPages(path) as pdf:
        figure, axis = plt.subplots(figsize=(11, 8.5))
        axis.axis("off")
        axis.text(0.03, 0.97, summary_text, va="top", family="monospace", fontsize=10)
        pdf.savefig(figure)
        plt.close(figure)

        figure = plt.figure(figsize=(11, 8.5))
        figure.suptitle("How to interpret the t-prime bins and calibrated error bars", fontsize=16, y=0.97)
        first = bin_statistics[0]
        last = bin_statistics[-1]
        least = min(bin_statistics, key=lambda row: row["effective_entries"])
        largest = max(bin_statistics, key=lambda row: row["effective_entries"])
        explanation = (
            "Bin numbering runs from large -t' toward the forward limit: bin 0 is "
            f"{first['tprime_low']:.6g} < t' < {first['tprime_high']:.6g}, while bin "
            f"{last['tprime_bin']} is {last['tprime_low']:.6g} < t' < "
            f"{last['tprime_high']:.6g}. The frozen weighted data vector identifies bin "
            f"{least['tprime_bin']} as the least-statistics bin and bin "
            f"{largest['tprime_bin']} as the largest-statistics bin. "
            "N_eff = (sum w)^2 / sum(w^2) is a weighted-statistics diagnostic, not a raw "
            f"event count. The {len(NAMES)} published values are derived from four global physics "
            "parameters (published covariance rank 4), including a shared U normalization "
            "and slope. Consequently the toy-calibrated delta68 values are correlated "
            "model-fit uncertainties and need not follow independent-bin 1/sqrt(N_eff) scaling."
        )
        figure.text(0.055, 0.89, textwrap.fill(explanation, 125), va="top", fontsize=10.2, linespacing=1.35)

        table_axis = figure.add_axes([0.055, 0.43, 0.89, 0.33])
        table_axis.axis("off")
        cell_text = []
        for row in bin_statistics:
            cell_text.append(
                [
                    f"{row['tprime_bin']}",
                    f"[{row['tprime_low']:.2f}, {row['tprime_high']:.2f}]",
                    f"{-row['response_weighted_tprime']:.3f}",
                    f"{row['included_phi_rows']}",
                    f"{row['weighted_yield']:.4f}",
                    f"{row['effective_entries']:.1f}",
                    f"{100.0 * row['relative_weighted_yield_error']:.2f}%",
                    f"{row['U_delta68_stat']:.3f}",
                ]
            )
        table = table_axis.table(
            cellText=cell_text,
            colLabels=["bin", "t' range", "-<t'>", "rows", "sum w", "N_eff", "sqrt(sum w^2)/sum w", "U delta68"],
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8.7)
        table.scale(1.0, 1.55)
        for column, width in enumerate((0.05, 0.14, 0.08, 0.06, 0.10, 0.09, 0.20, 0.10)):
            for row_index in range(len(cell_text) + 1):
                table[(row_index, column)].set_width(width)
        for cell in table.get_celld().values():
            cell.set_edgecolor("0.78")
        for column in range(8):
            table[(0, column)].set_facecolor("#dce8f2")
            table[(0, column)].set_text_props(weight="bold")

        bins = np.arange(N_TPRIME_BINS)
        effective = np.array([row["effective_entries"] for row in bin_statistics])
        axis_left = figure.add_axes([0.09, 0.10, 0.36, 0.25])
        axis_left.bar(bins, effective, color="#4c72b0")
        axis_left.set_yscale("log")
        axis_left.set_xticks(bins)
        axis_left.set_xlabel("t-prime bin")
        axis_left.set_ylabel("weighted effective entries")
        axis_left.set_title("Frozen-data statistics")
        axis_left.grid(axis="y", alpha=0.25)

        axis_right = figure.add_axes([0.57, 0.10, 0.36, 0.25])
        for key, label, marker in (
            ("U_relative_delta68", "U", "o"),
            ("LT_relative_delta68", "LT", "s"),
            ("TT_relative_delta68", "TT", "^"),
        ):
            axis_right.plot(bins, [100.0 * row[key] for row in bin_statistics], marker=marker, label=label)
        axis_right.set_xticks(bins)
        axis_right.set_xlabel("t-prime bin")
        axis_right.set_ylabel("delta68 / |central value| [%]")
        axis_right.set_title("Correlated calibrated radii")
        axis_right.grid(alpha=0.25)
        axis_right.legend(fontsize=8)
        pdf.savefig(figure)
        plt.close(figure)

        figure, axes = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)
        x = np.arange(len(NAMES))
        axes[0].plot(x, stats["toy_sd"] / curvature_errors, "o", label="toy SD / curvature sigma")
        axes[0].plot(x, stats["delta68"] / curvature_errors, "s", label="delta68 / curvature sigma")
        axes[0].axhline(1.0, color="0.5", lw=0.8)
        axes[0].set_ylabel("ratio")
        axes[0].legend()
        axes[1].errorbar(x, stats["coverage"], yerr=stats["coverage_se"], fmt="o", capsize=3)
        axes[1].axhline(0.68, color="0.4", lw=0.8)
        axes[1].axhspan(0.60, 0.76, color="#d7ecd9", alpha=0.7)
        axes[1].set_ylabel("five-fold held-out coverage")
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(NAMES, rotation=90)
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        figure, axis = plt.subplots(figsize=(11, 6))
        ratio = np.abs(stats["bias"]) / stats["delta68"]
        colors = ["#c44e52" if row["use_basic_interval"] else "#4c72b0" for row in details]
        axis.bar(np.arange(len(NAMES)), ratio, color=colors)
        axis.axhline(1.0, color="0.25", lw=1.0)
        axis.set_ylabel("|toy mean bias| / delta68")
        axis.set_xticks(np.arange(len(NAMES)))
        axis.set_xticklabels(NAMES, rotation=90)
        axis.set_ylim(0.0, max(1.05, 1.1 * ratio.max()))
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

        figure, axes = plt.subplots(2, 2, figsize=(11, 9))
        for axis, matrix, title, limits in (
            (axes[0, 0], correlation(c_curvature), "curvature correlation", (-1, 1)),
            (axes[0, 1], r_toy, "empirical toy correlation", (-1, 1)),
            (axes[1, 0], c_toy, "raw empirical toy covariance", None),
            (axes[1, 1], c68, "68%-calibrated covariance representation", None),
        ):
            if limits:
                image = axis.imshow(matrix, cmap="coolwarm", vmin=limits[0], vmax=limits[1])
            else:
                bound = np.max(np.abs(matrix))
                image = axis.imshow(matrix, cmap="coolwarm", vmin=-bound, vmax=bound)
            axis.set_title(title)
            axis.set_xticks(range(len(NAMES)))
            axis.set_xticklabels(NAMES, rotation=90, fontsize=6)
            axis.set_yticks(range(len(NAMES)))
            axis.set_yticklabels(NAMES, fontsize=6)
            figure.colorbar(image, ax=axis, fraction=0.046, pad=0.04)
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)


def markdown_table(headers, rows, alignments=None):
    if alignments is None:
        alignments = ["---"] * len(headers)
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(alignments) + "|"]
    lines.extend("| " + " | ".join(map(str, row)) + " |" for row in rows)
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "validation/preliminary_model_xsec_20261005")
    parser.add_argument("--curvature", type=Path, default=ROOT / "validation/preliminary_model_xsec_curvature_20261005")
    parser.add_argument("--output", type=Path, default=ROOT / "validation/preliminary_model_xsec_calibrated_20261005")
    parser.add_argument(
        "--config",
        type=Path,
        help="selected xsec JSON; must match the toy campaign config_snapshot.json",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="validate campaign/config compatibility without writing a release",
    )
    args = parser.parse_args()

    source = args.source.resolve()
    curvature = args.curvature.resolve()
    output = args.output.resolve()
    source_config_path = source / "config_snapshot.json"
    selected_config_path = args.config.resolve() if args.config else source_config_path

    source_paths = [
        source / "toys/replicas.npz",
        source / "toys/summary.json",
        source / "central_results.json",
        source_config_path,
        source / "before/preliminary_pi0_xsec_v7.py",
        curvature / "preliminary_cross_sections.csv",
        curvature / "published_covariance_curvature.csv",
        curvature / "physics_profile_intervals.csv",
        curvature / "covariance_metadata.json",
        curvature / "hessian_crosscheck.json",
    ]
    for path in source_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    if not selected_config_path.is_file():
        raise FileNotFoundError(selected_config_path)

    source_configuration = json.loads(source_config_path.read_text())
    selected_configuration = json.loads(selected_config_path.read_text())
    kinematic = source_configuration.get("configured_kinematic")
    generation_summary_path = source / "generation_summary.json"
    freshly_generated = False
    if generation_summary_path.is_file():
        generation_summary = json.loads(generation_summary_path.read_text())
        freshly_generated = generation_summary.get("freshly_generated") is True
    validate_selected_config(selected_configuration, source_configuration)
    tprime_edges = configure_binning(source_configuration)

    toys = np.load(source / "toys/replicas.npz")
    toy_summary = json.loads((source / "toys/summary.json").read_text())
    central = json.loads((source / "central_results.json").read_text())["additive"]
    curvature_metadata = json.loads((curvature / "covariance_metadata.json").read_text())
    if curvature_metadata.get("published_order") != NAMES:
        raise ValueError("curvature published ordering does not match selected binning")
    curve_rows, curve_central, curvature_errors = read_curvature_table(curvature / "preliminary_cross_sections.csv")
    c_curvature = read_matrix(curvature / "published_covariance_curvature.csv")
    toy_ids = validate_inputs(
        source,
        source_configuration,
        tprime_edges,
        toys,
        toy_summary,
        central,
        curve_rows,
        curve_central,
        curvature_errors,
        c_curvature,
    )

    if args.check_only:
        print(
            json.dumps(
                {
                    "status": "compatible",
                    "selected_config": str(selected_config_path),
                    "toy_source": str(source),
                    "curvature_source": str(curvature),
                    "tprime_bin_edges": tprime_edges.tolist(),
                    "published_quantities": len(NAMES),
                    "accepted_toys": int(toy_summary["accepted"]),
                    "freshly_generated": freshly_generated,
                },
                indent=2,
            )
        )
        return

    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    truth = np.asarray(central["pub"], dtype=float)
    residuals = np.asarray(toys["pub"], dtype=float) - truth
    stats = calibrate(residuals, toy_ids)
    c_toy = np.cov(residuals, rowvar=False, ddof=1)
    r_toy = correlation(c_toy)
    c68 = np.outer(stats["delta68"], stats["delta68"]) * r_toy
    bin_statistics = summarize_bin_statistics(central, curve_rows, truth, stats)

    details = []
    for index, name in enumerate(NAMES):
        use_basic, reason = classify_interval(stats["q16"][index], stats["q84"][index])
        details.append(
            {
                "quantity": name,
                "component": COMPONENTS[index // N_TPRIME_BINS],
                "tprime_bin": index % N_TPRIME_BINS,
                "central": truth[index],
                "curvature_sigma": curvature_errors[index],
                "toy_mean_bias": stats["bias"][index],
                "toy_sd": stats["toy_sd"][index],
                "toy_rmse": stats["rmse"][index],
                "delta68": stats["delta68"][index],
                "residual_q16": stats["q16"][index],
                "residual_q50": stats["q50"][index],
                "residual_q84": stats["q84"][index],
                "basic_interval_lower": truth[index] - stats["q84"][index],
                "basic_interval_upper": truth[index] - stats["q16"][index],
                "abs_bias_over_curvature_sigma": abs(stats["bias"][index]) / curvature_errors[index],
                "abs_bias_over_delta68": abs(stats["bias"][index]) / stats["delta68"][index],
                "toy_sd_over_curvature_sigma": stats["toy_sd"][index] / curvature_errors[index],
                "delta68_over_curvature_sigma": stats["delta68"][index] / curvature_errors[index],
                "cross_validated_coverage68": stats["coverage"][index],
                "coverage_binomial_standard_error": stats["coverage_se"][index],
                "plot_interval": "bias-aware basic 16-84%" if use_basic else "symmetric +/- delta68",
                "use_basic_interval": use_basic,
                "interval_reason": reason,
                "units": UNITS,
            }
        )

    min_coverage = float(stats["coverage"].min())
    max_bias_ratio = float(np.max(np.abs(stats["bias"]) / stats["delta68"]))
    finite_positive = bool(np.isfinite(stats["delta68"]).all() and np.all(stats["delta68"] > 0.0))
    catastrophic = (not finite_positive) or min_coverage < 0.55 or max_bias_ratio >= 1.0
    if catastrophic:
        raise RuntimeError(
            f"calibrated release gate failed: finite_positive={finite_positive}, "
            f"min_cv_coverage={min_coverage}, max_abs_bias_over_delta68={max_bias_ratio}"
        )

    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.staging.", dir=output.parent))
    try:
        manifest_paths = source_paths + (
            [selected_config_path] if selected_config_path != source_config_path else []
        )
        input_manifest = [
            {
                "path": str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in manifest_paths
        ]
        atomic_text(stage / "input_manifest.json", json.dumps(input_manifest, indent=2) + "\n")

        write_rows(stage / "toy_residual_statistics.csv", details)
        write_rows(stage / "cross_validated_coverage.csv", stats["coverage_rows"])
        write_rows(stage / "bin_statistics.csv", bin_statistics)
        write_matrix(stage / "published_covariance_curvature.csv", c_curvature)
        write_matrix(stage / "published_covariance_toy.csv", c_toy)
        write_matrix(stage / "published_correlation_toy.csv", r_toy)
        write_matrix(stage / "published_covariance_calibrated68.csv", c68)
        write_rows(
            stage / "curvature_vs_calibrated_errors.csv",
            [
                {
                    "quantity": row["quantity"],
                    "central": row["central"],
                    "curvature_sigma": row["curvature_sigma"],
                    "toy_sd": row["toy_sd"],
                    "toy_mean_bias": row["toy_mean_bias"],
                    "toy_rmse": row["toy_rmse"],
                    "calibrated_delta68": row["delta68"],
                    "abs_bias_over_curvature_sigma": row["abs_bias_over_curvature_sigma"],
                    "abs_bias_over_delta68": row["abs_bias_over_delta68"],
                    "cross_validated_coverage68": row["cross_validated_coverage68"],
                    "coverage_binomial_standard_error": row["coverage_binomial_standard_error"],
                }
                for row in details
            ],
        )

        means = np.array([float(row["response_weighted_tprime"]) for row in curve_rows])
        final_rows = []
        for bin_index, source_row in enumerate(curve_rows):
            final_rows.append(
                {
                    "tprime_bin": bin_index,
                    "tprime_low": source_row["tprime_low"],
                    "tprime_high": source_row["tprime_high"],
                    "response_weighted_tprime": means[bin_index],
                    "sigma_U": truth[bin_index],
                    "U_delta68_stat": stats["delta68"][bin_index],
                    "sigma_LT": truth[N_TPRIME_BINS + bin_index],
                    "LT_delta68_stat": stats["delta68"][N_TPRIME_BINS + bin_index],
                    "sigma_TT": truth[2 * N_TPRIME_BINS + bin_index],
                    "TT_delta68_stat": stats["delta68"][2 * N_TPRIME_BINS + bin_index],
                    "units": UNITS,
                    "status": "PRELIMINARY; toy-calibrated local 68% statistical uncertainties",
                }
            )
        write_rows(stage / "preliminary_cross_sections_calibrated.csv", final_rows)

        for component in COMPONENTS:
            make_component_plot(stage, component, means, truth, stats, details, kinematic)
        make_fit_diagnostics_scope_page(stage / "fit_diagnostics_scope.pdf")

        included_rows = int(np.count_nonzero(np.asarray(central["data_variance"]) > 0.0))
        parameter_count = len(central["theta"])
        nuisance_count = parameter_count - 4
        summary_text = (
            f"{kinematic} calibrated preliminary model extraction\n\n"
            f"Frozen M0 Q = {central['info'][0]:.12f}; included rows = {included_rows}\n"
            f"theta = {np.asarray(central['theta'][:4])}\n"
            f"fit parameters = {parameter_count} ({nuisance_count} nuisance + 4 physics)\n"
            "accepted local M0 toys = 500; folds = toy_id modulo 5\n"
            f"cross-validated coverage range = {stats['coverage'].min():.3f} to {stats['coverage'].max():.3f}\n"
            f"max |bias| / delta68 = {max_bias_ratio:.6f}\n\n"
            "Verdict: PRELIMINARY MODEL EXTRACTION READY"
        )
        diagnostics_pdf(
            stage / "final_preliminary_model_diagnostics.pdf",
            curvature_errors,
            stats,
            c_curvature,
            c_toy,
            r_toy,
            c68,
            details,
            summary_text,
            bin_statistics,
        )

        comparison_table = markdown_table(
            ["quantity", "central", "curv sigma", "toy SD", "bias", "RMSE", "delta68", "|bias|/curv", "|bias|/delta68", "CV cov +/- SE"],
            [
                [
                    row["quantity"],
                    f"{row['central']:.6f}",
                    f"{row['curvature_sigma']:.6f}",
                    f"{row['toy_sd']:.6f}",
                    f"{row['toy_mean_bias']:+.6f}",
                    f"{row['toy_rmse']:.6f}",
                    f"{row['delta68']:.6f}",
                    f"{row['abs_bias_over_curvature_sigma']:.3f}",
                    f"{row['abs_bias_over_delta68']:.3f}",
                    f"{row['cross_validated_coverage68']:.3f} +/- {row['coverage_binomial_standard_error']:.3f}",
                ]
                for row in details
            ],
            ["---", "---:", "---:", "---:", "---:", "---:", "---:", "---:", "---:", "---:"],
        )
        final_table = markdown_table(
            ["bin", "response-weighted <t'>", "sigma_U", "U delta68 stat", "sigma_LT", "LT delta68 stat", "sigma_TT", "TT delta68 stat"],
            [
                [
                    row["tprime_bin"],
                    f"{row['response_weighted_tprime']:.9f}",
                    f"{row['sigma_U']:.6f}",
                    f"{row['U_delta68_stat']:.6f}",
                    f"{row['sigma_LT']:.6f}",
                    f"{row['LT_delta68_stat']:.6f}",
                    f"{row['sigma_TT']:.6f}",
                    f"{row['TT_delta68_stat']:.6f}",
                ]
                for row in final_rows
            ],
            ["---:"] * 8,
        )
        statistics_table = markdown_table(
            ["bin", "t' range", "-<t'>", "included rows", "sum w", "N_eff", "relative weighted error", "U delta68"],
            [
                [
                    row["tprime_bin"],
                    f"[{row['tprime_low']:.2f}, {row['tprime_high']:.2f}]",
                    f"{-row['response_weighted_tprime']:.6f}",
                    row["included_phi_rows"],
                    f"{row['weighted_yield']:.6f}",
                    f"{row['effective_entries']:.1f}",
                    f"{100.0 * row['relative_weighted_yield_error']:.2f}%",
                    f"{row['U_delta68_stat']:.6f}",
                ]
                for row in bin_statistics
            ],
            ["---:", "---", "---:", "---:", "---:", "---:", "---:", "---:"],
        )
        asymmetric = ", ".join(row["quantity"] for row in details if row["use_basic_interval"])
        least_statistics = min(bin_statistics, key=lambda row: row["effective_entries"])
        largest_statistics = max(bin_statistics, key=lambda row: row["effective_entries"])
        first_bin = bin_statistics[0]
        last_bin = bin_statistics[-1]
        published_rank = int(np.linalg.matrix_rank(c_curvature, tol=1e-8))
        fisher_rank = curvature_metadata.get("fisher_rank", "not recorded")
        fisher_condition = curvature_metadata.get("scaled_fisher_condition", "not recorded")
        jacobian_condition = curvature_metadata.get("scaled_jacobian_condition", "not recorded")
        report = f"""# {kinematic} preliminary calibrated model extraction

## Release decision

**PRELIMINARY MODEL EXTRACTION READY.** The 500 accepted local M0 toys are finite and internally consistent. Deterministic five-fold held-out coverage is `{stats['coverage'].min():.3f}`-`{stats['coverage'].max():.3f}` for all {len(NAMES)} published quantities, and the maximum `|bias| / delta68` is `{max_bias_ratio:.3f}`. No narrowly defined catastrophic calibrated-interval failure is present.

The frozen constrained M0 fit is unchanged: `Q = {central['info'][0]:.12f}` for {included_rows} included rows, with `(N_U, DeltaB_U, N_LT, N_TT) = ({central['theta'][0]:.12f}, {central['theta'][1]:.12f} GeV^-2, {central['theta'][2]:.12f}, {central['theta'][3]:.12f})`. It contains {parameter_count} parameters ({nuisance_count} nuisance and 4 physics); the Fisher rank is `{fisher_rank}`. The scaled Jacobian and normal-matrix conditions are `{jacobian_condition}` and `{fisher_condition}`. The selected JSON binning was checked exactly against the toy campaign snapshot and the curvature table before publication. No extraction code, central fit, production, or toy ensemble was changed or rerun by this reporting step.

## Statistical prescription

Central values are the constrained forward-folded M0 minimum-Q estimates. The conventional local covariance is retained from `F = J^T W J` after nuisance profiling, together with the existing DeltaQ=1 profiles. Positivity constraints produce finite-sample bias in the local toys, so each published residual is `d = z_fit - z_truth` and the primary symmetric preliminary statistical error is `delta68 = Q_0.68(|d|)`. This is a toy-calibrated local 68% statistical uncertainty, not one Gaussian sigma. The central values are not bias-corrected, and the bias is not added again in quadrature.

Five folds are assigned deterministically by `toy_id modulo 5`. In each fold, `delta68` is learned from the other four folds using NumPy's linear sample quantile and tested only on held-out toys. Each of the 500 toys contributes exactly one held-out decision per quantity. The quoted binomial uncertainty is `sqrt(p(1-p)/500)`.

The empirical covariance `C_toy = cov(d)` supplies the preliminary statistical correlation matrix. `C68_representation = diag(delta68) R_toy diag(delta68)` is saved only as a **68%-calibrated covariance representation**; its diagonal matches the plotted symmetric radii, but `delta68^2` is not asserted to be a variance.

## Interpreting bin statistics and correlated U errors

Bin numbering runs from large `-t'` toward the forward limit. Bin 0 is `{first_bin['tprime_low']:.6g} < t' < {first_bin['tprime_high']:.6g}`; bin {last_bin['tprime_bin']} is `{last_bin['tprime_low']:.6g} < t' < {last_bin['tprime_high']:.6g}`. In the frozen additive-weight data vector, bin {least_statistics['tprime_bin']} is the least-statistics bin (`N_eff = {least_statistics['effective_entries']:.1f}`), while bin {largest_statistics['tprime_bin']} has the largest weighted effective statistics (`N_eff = {largest_statistics['effective_entries']:.1f}`). Here `N_eff = (sum w)^2 / sum(w^2)` is a weighted-data diagnostic, not a raw event count.

The {len(NAMES)} displayed cross sections are not {len(NAMES)} independently fitted bin contents. They are derived from four global physics parameters, so the published curvature covariance has rank {published_rank}. In particular, all {N_TPRIME_BINS} U values share one normalization and one slope. The plotted `delta68` values are correlated refit-toy residual radii for that constrained model; they therefore need not obey independent-bin `1/sqrt(N)` scaling. The U radius in bin {largest_statistics['tprime_bin']} is interpreted together with its weighted statistics and the global shape constraint, not as a standalone counting error.

{statistics_table}

## Curvature-to-calibration comparison

All values and errors are in `{UNITS}`.

{comparison_table}

## Final preliminary cross sections

The table reports the frozen central values and symmetric calibrated `delta68` radii. Detailed bias-aware basic intervals are in `toy_residual_statistics.csv`.

{final_table}

## Plot interval treatment

The basic interval is `z_obs - q84` to `z_obs - q16`. It is used in the main plot wherever it excludes the central estimate or its two half-widths differ by at least a factor of 3. Quantities meeting that explicit strong-asymmetry rule are: {asymmetric}. For those points, a thin gray symmetric `+/- delta68` reference is also shown. The TT basic intervals lie below the frozen point estimates because the constrained-estimator residual distribution is shifted positive; the central values remain unshifted. All other points use symmetric calibrated `+/- delta68` bars.

## Products

- Final table: `preliminary_cross_sections_calibrated.csv`
- Residual statistics and basic intervals: `toy_residual_statistics.csv`
- Fold and aggregate coverage: `cross_validated_coverage.csv`
- Frozen weighted-statistics audit: `bin_statistics.csv`
- Conventional covariance: `published_covariance_curvature.csv`
- Raw toy covariance/correlation: `published_covariance_toy.csv`, `published_correlation_toy.csv`
- Calibrated representation: `published_covariance_calibrated68.csv`
- Error comparison: `curvature_vs_calibrated_errors.csv`
- Figures: `sigma_U_preliminary_calibrated.pdf/png`, `sigma_LT_preliminary_calibrated.pdf/png`, `sigma_TT_preliminary_calibrated.pdf/png`
- Combined diagnostics: `final_preliminary_model_diagnostics.pdf`
- Unified-report appendix divider: `fit_diagnostics_scope.pdf`
- Exact source hashes: `input_manifest.json`

## Deferred refinements

The physical-event bootstrap remains a conservative robustness diagnostic, not the primary error prescription. Background/combinatorial timing transport, ellipse systematics, SIDIS and delta studies, no-model comparison, charged-pion/model refinements, LT/TT slopes, and other production-scale investigations remain deferred. None was rerun here.

## Reproduction

From the repository root:

```bash
cd {ROOT}
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 scripts/report_preliminary_pi0_calibrated.py \\
  --source {source} \\
  --curvature {curvature} \\
  --config {selected_config_path} \\
  --output {output}_reproduced
```

The command requires NumPy and Matplotlib, performs no ROOT or production work, validates that selected and campaign binning agree, refuses to overwrite an existing destination, stages all files, and atomically publishes the output directory.
"""
        atomic_text(stage / "REPORT.md", report)

        summary = {
            "status": "PRELIMINARY MODEL EXTRACTION READY",
            "configured_kinematic": kinematic,
            "frozen_m0_q": central["info"][0],
            "frozen_parameters": central["theta"][:4],
            "accepted_toys": int(toy_summary["accepted"]),
            "fold_rule": "toy_id modulo 5",
            "minimum_cross_validated_coverage68": min_coverage,
            "maximum_cross_validated_coverage68": float(stats["coverage"].max()),
            "maximum_abs_bias_over_delta68": max_bias_ratio,
            "published_covariance_rank": published_rank,
            "tprime_bin_edges": tprime_edges.tolist(),
            "selected_config": str(selected_config_path),
            "selected_config_sha256": sha256(selected_config_path),
            "toy_config_snapshot_sha256": sha256(source_config_path),
            "least_statistics_tprime_bin": int(least_statistics["tprime_bin"]),
            "largest_statistics_tprime_bin": int(largest_statistics["tprime_bin"]),
            "central_extraction_changed": False,
            "production_rerun": False,
            "toys_rerun": freshly_generated,
        }
        atomic_text(stage / "calibration_summary.json", json.dumps(summary, indent=2) + "\n")
        stage.replace(output)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise

    print(json.dumps({"output": str(output), "summary": summary}, indent=2))


if __name__ == "__main__":
    main()
