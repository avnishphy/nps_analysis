#!/usr/bin/env python3
"""
Publication-quality cross-section plotting utility for the no-SIMC-model extraction.

This script reads:
1) Slice-fit coefficients CSV (per t', Q2, xB bin)
2) Per-phi cross-section CSV (per slice, per phi bin)

It writes a single multi-page PDF containing:
- Overall coefficient plots: dσ/dt vs -t' for each (iq, ix)
- Per-slice plots: dσ/(dt dphi) vs phi with model components and uncertainty bands

All sigma-like values are displayed in nb/GeV^2 by applying a configurable scale.

Current extraction convention:
- SIMC acceptance basis uses full_weight/sigcm.
- Gamma is kept in the CSV as a diagnostic, but is not multiplied into the
  plotted response basis because full_weight/sigcm retains siglab/sigcm =
  davejac * gtpr * fac.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec


SliceKey = Tuple[int, int, int]   # (it, iq, ix)
GroupKey = Tuple[int, int]        # (iq, ix)

COEFF_META = [
    ("fit_xsec_sigmaU", "fit_xsec_sigmaUerr", r"$\sigma_U$", "#1f77b4", "o"),
    ("fit_xsec_sigmaTL", "fit_xsec_sigmaTLerr", r"$\sigma_{LT}$", "#d62728", "s"),
    ("fit_xsec_sigmaTT", "fit_xsec_sigmaTTerr", r"$\sigma_{TT}$", "#2ca02c", "^"),
    ("fit_asym_sigmaTLp", "fit_asym_sigmaTLperr", r"$\sigma_{TL'}$", "#9467bd", "D"),
]


def to_float(value: str, default: float = float("nan")) -> float:
    try:
        return float(value)
    except Exception:
        return default


def to_int(value: str, default: int = -1) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def clean_err(arr: np.ndarray) -> np.ndarray:
    return np.where(np.isfinite(arr) & (arr >= 0.0), arr, 0.0)


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def base_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 300,
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "legend.fontsize": 9.8,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": "--",
            "axes.linewidth": 1.15,
        }
    )


def parse_slice_csv(slice_csv: str) -> Tuple[Dict[GroupKey, List[dict]], Dict[SliceKey, dict]]:
    by_group: Dict[GroupKey, List[dict]] = defaultdict(list)
    by_slice: Dict[SliceKey, dict] = {}

    with open(slice_csv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if to_int(row.get("fit_xsec_ok", "0"), 0) != 1:
                continue

            it = to_int(row.get("it", "-1"), -1)
            iq = to_int(row.get("iq", "-1"), -1)
            ix = to_int(row.get("ix", "-1"), -1)
            if it < 0 or iq < 0 or ix < 0:
                continue

            rec = {
                "it": it,
                "iq": iq,
                "ix": ix,
                "tprime_lo": to_float(row.get("tprime_lo", "nan")),
                "tprime_hi": to_float(row.get("tprime_hi", "nan")),
                "tprime_center": to_float(row.get("tprime_center", "nan")),
                "q2_lo": to_float(row.get("q2_lo", "nan")),
                "q2_hi": to_float(row.get("q2_hi", "nan")),
                "xb_lo": to_float(row.get("xb_lo", "nan")),
                "xb_hi": to_float(row.get("xb_hi", "nan")),
                "epsilon": to_float(row.get("epsilon", "nan")),
                "gamma_flux": to_float(row.get("gamma_flux", "nan")),
                "fit_xsec_chi2": to_float(row.get("fit_xsec_chi2", "nan")),
                "fit_xsec_ndf": to_float(row.get("fit_xsec_ndf", "nan")),
                "fit_asym_ok": to_int(row.get("fit_asym_ok", "0"), 0),
            }

            for val_key, err_key, _, _, _ in COEFF_META:
                rec[val_key] = to_float(row.get(val_key, "nan"))
                rec[err_key] = to_float(row.get(err_key, "nan"))

            if not math.isfinite(rec["tprime_center"]):
                continue

            gk = (iq, ix)
            sk = (it, iq, ix)
            by_group[gk].append(rec)
            by_slice[sk] = rec

    for gk in by_group:
        by_group[gk].sort(key=lambda r: r["tprime_center"])

    return by_group, by_slice


def parse_phi_csv(phi_csv: str) -> Dict[SliceKey, List[dict]]:
    by_slice_phi: Dict[SliceKey, List[dict]] = defaultdict(list)

    with open(phi_csv, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            it = to_int(row.get("it", "-1"), -1)
            iq = to_int(row.get("iq", "-1"), -1)
            ix = to_int(row.get("ix", "-1"), -1)
            ip = to_int(row.get("ip", "-1"), -1)
            if it < 0 or iq < 0 or ix < 0 or ip < 0:
                continue

            phi_center = to_float(row.get("phi_center", "nan"))
            xsec = to_float(row.get("xsec", "nan"))
            xsec_err = to_float(row.get("xsec_err", "nan"))
            data_y = to_float(row.get("data", "nan"))
            data_err = to_float(row.get("data_err", "nan"))
            sim_y = to_float(row.get("sim", "nan"))
            sim_err = to_float(row.get("sim_err", "nan"))
            ratio = to_float(row.get("ratio", "nan"))
            ratio_err = to_float(row.get("ratio_err", "nan"))
            if not (math.isfinite(phi_center) and math.isfinite(xsec)):
                continue

            rec = {
                "ip": ip,
                "phi_lo": to_float(row.get("phi_lo", "nan")),
                "phi_hi": to_float(row.get("phi_hi", "nan")),
                "phi_center": phi_center,
                "data": data_y,
                "data_err": data_err,
                "sim": sim_y,
                "sim_err": sim_err,
                "ratio": ratio,
                "ratio_err": ratio_err,
                "xsec": xsec,
                "xsec_err": xsec_err,
                "epsilon": to_float(row.get("epsilon", "nan")),
                "gamma_flux": to_float(row.get("gamma_flux", "nan")),
            }
            by_slice_phi[(it, iq, ix)].append(rec)

    for sk in by_slice_phi:
        by_slice_phi[sk].sort(key=lambda r: r["phi_center"])

    return by_slice_phi


def valid_tlp_group(rows: List[dict]) -> bool:
    vals = np.array([r.get("fit_asym_sigmaTLp", np.nan) for r in rows], dtype=float)
    ers = np.array([r.get("fit_asym_sigmaTLperr", np.nan) for r in rows], dtype=float)
    return bool(np.any(np.isfinite(vals) & ((np.abs(vals) > 0.0) | (ers > 0.0))))


def compute_offsets(x: np.ndarray, nseries: int) -> np.ndarray:
    if nseries <= 1:
        return np.zeros(1)

    uniq = np.unique(np.round(x, 8))
    if uniq.size > 1:
        spacing = np.diff(np.sort(uniq))
        dx = float(np.nanmedian(spacing))
    else:
        dx = max(float(np.nanmax(np.abs(x))) * 0.03, 0.004)

    if not math.isfinite(dx) or dx <= 0.0:
        dx = 0.004

    jitter = 0.10 * dx
    mid = 0.5 * (nseries - 1)
    return np.array([(i - mid) * jitter for i in range(nseries)], dtype=float)


def draw_title_page(
    pdf: PdfPages,
    slice_csv: str,
    phi_csv: str,
    n_groups: int,
    n_slices_fit: int,
    n_slices_phi: int,
    nb_scale: float,
) -> None:
    fig = plt.figure(figsize=(11.0, 8.5))
    ax = fig.add_subplot(111)
    ax.axis("off")

    lines = [
        "Pi0 Cross-Section Publication Plot Pack",
        "",
        f"Slice-fit CSV: {slice_csv}",
        f"Per-phi CSV:  {phi_csv}",
        f"(Q2, xB) groups with fitted slices: {n_groups}",
        f"Fitted slices in coefficient tables: {n_slices_fit}",
        f"Slices available in per-phi tables:  {n_slices_phi}",
        "",
        "Units and conventions:",
        f"- Sigma values scaled by {nb_scale:g} to display nb/GeV^2",
        "- Overall plots: y-axis is dσ/dt [nb/GeV^2]",
        "- Per-slice contribution/helicity panels: y-axis is dσ/(dt dφ) [nb/GeV^2]",
        "- SIMC basis is full_weight/sigcm, preserving davejac*gtpr*fac through siglab/sigcm",
        "- Gamma is read from CSV for diagnostics only; it is not multiplied into plotted components",
        "",
        "Per-slice multipanel content:",
        "- Data vs SIMC yields",
        "- Data/SIMC ratio",
        "- Contribution decomposition with ±1σ bands",
        "- Helicity placeholder panel (enabled when fit_asym_ok=1)",
    ]

    y = 0.93
    for i, text in enumerate(lines):
        fs = 18 if i == 0 else 12
        fw = "bold" if i == 0 else "normal"
        ax.text(0.05, y, text, fontsize=fs, fontweight=fw, va="top")
        y -= 0.056 if i == 0 else 0.046

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def draw_overall_group_page(pdf: PdfPages, rows: List[dict], iq: int, ix: int, nb_scale: float) -> None:
    if not rows:
        return

    has_tlp = valid_tlp_group(rows)
    specs = COEFF_META if has_tlp else COEFF_META[:3]

    x = -np.array([r["tprime_center"] for r in rows], dtype=float)
    xerr = 0.5 * np.abs(
        np.array([r["tprime_hi"] for r in rows], dtype=float)
        - np.array([r["tprime_lo"] for r in rows], dtype=float)
    )

    q2_lo = rows[0]["q2_lo"]
    q2_hi = rows[0]["q2_hi"]
    xb_lo = rows[0]["xb_lo"]
    xb_hi = rows[0]["xb_hi"]

    fig = plt.figure(figsize=(11.0, 12.0))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.2, 1.0, 1.0], hspace=0.35, wspace=0.24)

    ax_overlay = fig.add_subplot(gs[0, :])
    ax_u = fig.add_subplot(gs[1, 0])
    ax_lt = fig.add_subplot(gs[1, 1])
    ax_tt = fig.add_subplot(gs[2, 0])
    ax_tlp = fig.add_subplot(gs[2, 1])

    offsets = compute_offsets(x, len(specs))

    for idx, (vkey, ekey, label, color, marker) in enumerate(specs):
        y = nb_scale * np.array([r[vkey] for r in rows], dtype=float)
        yerr = nb_scale * clean_err(np.array([r[ekey] for r in rows], dtype=float))

        ax_overlay.errorbar(
            x + offsets[idx],
            y,
            xerr=xerr,
            yerr=yerr,
            fmt=marker,
            linestyle="none",
            color=color,
            markerfacecolor="white",
            markeredgewidth=1.45,
            markersize=6.8,
            capsize=3.0,
            elinewidth=1.25,
            label=label,
            zorder=3,
        )

    ax_overlay.axhline(0.0, color="0.35", linestyle="--", linewidth=1.1)
    ax_overlay.set_title(
        rf"Fitted Coefficients vs $-t'$  |  $Q^2 \in [{q2_lo:.3f}, {q2_hi:.3f}]$, $x_B \in [{xb_lo:.3f}, {xb_hi:.3f}]$"
    )
    ax_overlay.set_xlabel(r"$-t'$ [GeV$^2$]")
    ax_overlay.set_ylabel(r"$d\sigma/dt$ [nb/GeV$^2$]")
    ax_overlay.legend(loc="best", framealpha=0.96)

    indiv_axes = [ax_u, ax_lt, ax_tt, ax_tlp]
    indiv_specs = [COEFF_META[0], COEFF_META[1], COEFF_META[2], COEFF_META[3]]

    for ax, (vkey, ekey, label, color, marker) in zip(indiv_axes, indiv_specs):
        if label == r"$\sigma_{TL'}$" and not has_tlp:
            ax.axis("off")
            ax.text(0.5, 0.5, "No helicity TL' extraction in this group", ha="center", va="center", transform=ax.transAxes)
            continue

        y = nb_scale * np.array([r[vkey] for r in rows], dtype=float)
        yerr = nb_scale * clean_err(np.array([r[ekey] for r in rows], dtype=float))

        ax.errorbar(
            x,
            y,
            xerr=xerr,
            yerr=yerr,
            fmt=marker,
            linestyle="none",
            color=color,
            markerfacecolor="white",
            markeredgewidth=1.35,
            markersize=6.3,
            capsize=2.8,
            elinewidth=1.15,
            zorder=3,
        )
        ax.axhline(0.0, color="0.4", linestyle="--", linewidth=1.0)
        ax.set_title(label)
        ax.set_xlabel(r"$-t'$ [GeV$^2$]")
        ax.set_ylabel(r"$d\sigma/dt$ [nb/GeV$^2$]")

    chi2_vals = np.array([r["fit_xsec_chi2"] for r in rows], dtype=float)
    ndf_vals = np.array([r["fit_xsec_ndf"] for r in rows], dtype=float)
    chi2_ndf = np.divide(
        chi2_vals,
        ndf_vals,
        out=np.full_like(chi2_vals, np.nan),
        where=(ndf_vals > 0.0),
    )
    finite = chi2_ndf[np.isfinite(chi2_ndf)]
    qtxt = rf"Median $\chi^2/ndf$: {np.median(finite):.2f}" if finite.size else r"Median $\chi^2/ndf$: n/a"

    fig.suptitle(rf"Overall Coefficients for Group (iq={iq}, ix={ix})  |  {qtxt}", fontsize=14, y=0.988)
    fig.text(
        0.5,
        0.012,
        r"Scatter-only markers with error bars. Coefficients interpreted as $d\sigma/dt$ terms.",
        ha="center",
        fontsize=10,
    )

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def component_model(phi: np.ndarray, fit_row: dict) -> Dict[str, np.ndarray]:
    eps = clamp(fit_row.get("epsilon", 0.0), 0.0, 1.0)
    inv2pi = 1.0 / (2.0 * math.pi)
    # No explicit Gamma factor here. The extraction now fits against the
    # full_weight/sigcm SIMC basis, which already keeps siglab/sigcm =
    # davejac*gtpr*fac in the event response.
    k_lt = math.sqrt(max(0.0, 2.0 * eps * (1.0 + eps)))
    k_tlp = math.sqrt(max(0.0, 2.0 * eps * (1.0 - eps)))

    sigma_u = fit_row.get("fit_xsec_sigmaU", float("nan"))
    sigma_lt = fit_row.get("fit_xsec_sigmaTL", float("nan"))
    sigma_tt = fit_row.get("fit_xsec_sigmaTT", float("nan"))
    sigma_tlp = fit_row.get("fit_asym_sigmaTLp", float("nan"))

    sigma_u_err = max(0.0, fit_row.get("fit_xsec_sigmaUerr", 0.0))
    sigma_lt_err = max(0.0, fit_row.get("fit_xsec_sigmaTLerr", 0.0))
    sigma_tt_err = max(0.0, fit_row.get("fit_xsec_sigmaTTerr", 0.0))
    sigma_tlp_err = max(0.0, fit_row.get("fit_asym_sigmaTLperr", 0.0))

    c1 = np.cos(phi)
    c2 = np.cos(2.0 * phi)
    s1 = np.sin(phi)

    u = inv2pi * sigma_u * np.ones_like(phi)
    lt = inv2pi * k_lt * sigma_lt * c1
    tt = inv2pi * eps * sigma_tt * c2
    tlp = inv2pi * k_tlp * sigma_tlp * s1

    # Unpolarized total used in xsec(phi) panel excludes TL'.
    total = u + lt + tt

    eu = inv2pi * sigma_u_err * np.ones_like(phi)
    elt = np.abs(inv2pi * k_lt * sigma_lt_err * c1)
    ett = np.abs(inv2pi * eps * sigma_tt_err * c2)
    etlp = np.abs(inv2pi * k_tlp * sigma_tlp_err * s1)

    # Diagonal approximation for total uncertainty band.
    etotal = np.sqrt(eu * eu + elt * elt + ett * ett)

    return {
        "u": u,
        "lt": lt,
        "tt": tt,
        "tlp": tlp,
        "total": total,
        "eu": eu,
        "elt": elt,
        "ett": ett,
        "etlp": etlp,
        "etotal": etotal,
    }


def draw_per_slice_page(
    pdf: PdfPages,
    fit_row: dict,
    phi_rows: List[dict],
    nb_scale: float,
) -> None:
    if not phi_rows:
        return

    it = fit_row["it"]
    iq = fit_row["iq"]
    ix = fit_row["ix"]

    tlo = fit_row["tprime_lo"]
    thi = fit_row["tprime_hi"]
    q2lo = fit_row["q2_lo"]
    q2hi = fit_row["q2_hi"]
    xblo = fit_row["xb_lo"]
    xbhi = fit_row["xb_hi"]
    eps = fit_row.get("epsilon", float("nan"))
    gamma_flux = fit_row.get("gamma_flux", float("nan"))

    phi_data = np.array([r["phi_center"] for r in phi_rows], dtype=float)
    xsec_data = nb_scale * np.array([r["xsec"] for r in phi_rows], dtype=float)
    xsec_err = nb_scale * clean_err(np.array([r["xsec_err"] for r in phi_rows], dtype=float))

    yield_data = np.array([r.get("data", np.nan) for r in phi_rows], dtype=float)
    yield_data_err = clean_err(np.array([r.get("data_err", np.nan) for r in phi_rows], dtype=float))
    yield_sim = np.array([r.get("sim", np.nan) for r in phi_rows], dtype=float)
    yield_sim_err = clean_err(np.array([r.get("sim_err", np.nan) for r in phi_rows], dtype=float))

    ratio_data = np.array([r.get("ratio", np.nan) for r in phi_rows], dtype=float)
    ratio_err = clean_err(np.array([r.get("ratio_err", np.nan) for r in phi_rows], dtype=float))
    x_err = 0.5 * np.abs(
        np.array([r.get("phi_hi", np.nan) for r in phi_rows], dtype=float)
        - np.array([r.get("phi_lo", np.nan) for r in phi_rows], dtype=float)
    )
    x_err = np.where(np.isfinite(x_err), x_err, 0.0)

    phi_min = float(np.nanmin(phi_data)) if phi_data.size else 0.0
    phi_max = float(np.nanmax(phi_data)) if phi_data.size else 2.0 * math.pi
    phi_curve = np.linspace(phi_min, phi_max, 600)

    model = component_model(phi_curve, fit_row)

    u = nb_scale * model["u"]
    lt = nb_scale * model["lt"]
    tt = nb_scale * model["tt"]
    tlp = nb_scale * model["tlp"]
    total = nb_scale * model["total"]

    eu = nb_scale * model["eu"]
    elt = nb_scale * model["elt"]
    ett = nb_scale * model["ett"]
    etlp = nb_scale * model["etlp"]
    etotal = nb_scale * model["etotal"]

    has_tlp = bool(fit_row.get("fit_asym_ok", 0) == 1 and (np.any(np.abs(tlp) > 0.0) or np.any(etlp > 0.0)))

    fig = plt.figure(figsize=(12, 8))
    if has_tlp:
        gs = GridSpec(2, 2, figure=fig, height_ratios=[0.72, 1.28], hspace=0.30, wspace=0.24)
        ax_yield = fig.add_subplot(gs[0, 0])
        ax_ratio = fig.add_subplot(gs[0, 1])
        ax_xsec = fig.add_subplot(gs[1, 0])
        ax_hel = fig.add_subplot(gs[1, 1])
    else:
        gs = GridSpec(2, 2, figure=fig, height_ratios=[0.68, 1.42], hspace=0.30, wspace=0.24)
        ax_yield = fig.add_subplot(gs[0, 0])
        ax_ratio = fig.add_subplot(gs[0, 1])
        ax_xsec = fig.add_subplot(gs[1, :])
        ax_hel = None

    # Panel 1: data and SIMC yields.
    mask_data = np.isfinite(yield_data)
    if np.any(mask_data):
        ax_yield.errorbar(
            phi_data[mask_data],
            yield_data[mask_data],
            xerr=x_err[mask_data],
            yerr=yield_data_err[mask_data],
            fmt="o",
            linestyle="none",
            color="black",
            markerfacecolor="white",
            markeredgewidth=1.35,
            markersize=5.9,
            capsize=2.8,
            elinewidth=1.05,
            label="Data yield",
            zorder=4,
        )

    mask_sim = np.isfinite(yield_sim)
    if np.any(mask_sim):
        ax_yield.errorbar(
            phi_data[mask_sim],
            yield_sim[mask_sim],
            xerr=x_err[mask_sim],
            yerr=yield_sim_err[mask_sim],
            fmt="s",
            linestyle="none",
            color="#b91c1c",
            markerfacecolor="white",
            markeredgewidth=1.25,
            markersize=5.5,
            capsize=2.8,
            elinewidth=1.0,
            label="SIMC yield",
            zorder=3,
        )

    ax_yield.set_title("Per-slice yields", fontsize=11.5)
    ax_yield.set_xlabel(r"$\phi$ [rad]")
    ax_yield.set_ylabel("Weighted yield")
    ax_yield.legend(loc="best", framealpha=0.95)

    # Panel 2: ratio.
    mask_ratio = np.isfinite(ratio_data)
    if np.any(mask_ratio):
        ax_ratio.errorbar(
            phi_data[mask_ratio],
            ratio_data[mask_ratio],
            xerr=x_err[mask_ratio],
            yerr=ratio_err[mask_ratio],
            fmt="o",
            linestyle="none",
            color="#111827",
            markerfacecolor="white",
            markeredgewidth=1.35,
            markersize=5.9,
            capsize=2.8,
            elinewidth=1.05,
            label="Data/SIMC",
            zorder=4,
        )

    ax_ratio.axhline(1.0, color="#4b5563", linestyle="--", linewidth=1.0)
    ax_ratio.set_title("Per-slice ratio", fontsize=11.5)
    ax_ratio.set_xlabel(r"$\phi$ [rad]")
    ax_ratio.set_ylabel("Data/SIMC")
    ax_ratio.legend(loc="best", framealpha=0.95)

    # Panel 3: cross-section contributions.
    ax_xsec.errorbar(
        phi_data,
        xsec_data,
        xerr=x_err,
        yerr=xsec_err,
        fmt="o",
        linestyle="none",
        color="black",
        markerfacecolor="white",
        markeredgewidth=1.45,
        markersize=6.5,
        capsize=3.0,
        elinewidth=1.15,
        label="Extracted slice points",
        zorder=4,
    )

    # Uncertainty bands.
    ax_xsec.fill_between(phi_curve, total - etotal, total + etotal, color="#6b7280", alpha=0.20, linewidth=0.0, label="Total model ±1σ")
    ax_xsec.fill_between(phi_curve, u - eu, u + eu, color="#1f77b4", alpha=0.10, linewidth=0.0)
    ax_xsec.fill_between(phi_curve, lt - elt, lt + elt, color="#d62728", alpha=0.10, linewidth=0.0)
    ax_xsec.fill_between(phi_curve, tt - ett, tt + ett, color="#2ca02c", alpha=0.10, linewidth=0.0)

    # TL' is helicity-odd and not part of unpolarized total; shown as reference if present.
    if has_tlp:
        ax_xsec.fill_between(phi_curve, tlp - etlp, tlp + etlp, color="#9467bd", alpha=0.08, linewidth=0.0)

    # Model lines.
    ax_xsec.plot(phi_curve, total, color="#374151", linewidth=2.8, label="Total model")
    ax_xsec.plot(phi_curve, u, color="#1f77b4", linewidth=2.0, linestyle="--", label=r"$\sigma_U$ contribution")
    ax_xsec.plot(phi_curve, lt, color="#d62728", linewidth=2.0, linestyle="-.", label=r"$\sigma_{LT}\cos\phi$ contribution")
    ax_xsec.plot(phi_curve, tt, color="#2ca02c", linewidth=2.0, linestyle=":", label=r"$\sigma_{TT}\cos2\phi$ contribution")
    if has_tlp:
        ax_xsec.plot(phi_curve, tlp, color="#9467bd", linewidth=1.9, linestyle=(0, (4, 2)), label=r"$\sigma_{TL'}\sin\phi$ reference")

    ax_xsec.axhline(0.0, color="0.4", linestyle="--", linewidth=1.0)

    ax_xsec.set_xlabel(r"$\phi$ [rad]")
    ax_xsec.set_ylabel(r"$d\sigma/(dt\,d\phi)$ [nb/GeV$^2$]")
    ax_xsec.set_title(
        rf"Per-Slice Cross Section  (it={it}, iq={iq}, ix={ix})"
    )

    info = (
        rf"$Q^2\in[{q2lo:.3f},{q2hi:.3f}]$ GeV$^2$, "
        rf"$x_B\in[{xblo:.3f},{xbhi:.3f}]$, "
        rf"$t'\in[{tlo:.3f},{thi:.3f}]$ GeV$^2$, "
        rf"$\epsilon={eps:.3f}$"
    )
    ax_xsec.text(0.012, 0.98, info, transform=ax_xsec.transAxes, va="top", fontsize=10.0)

    if math.isfinite(gamma_flux):
        ax_xsec.text(
            0.012,
            0.94,
            rf"$\Gamma={gamma_flux:.3e}$ stored only; not multiplied in this basis",
            transform=ax_xsec.transAxes,
            va="top",
            fontsize=9.3,
        )

    chi2 = fit_row.get("fit_xsec_chi2", float("nan"))
    ndf = fit_row.get("fit_xsec_ndf", float("nan"))
    if math.isfinite(chi2) and math.isfinite(ndf) and ndf > 0.0:
        qtxt = rf"Fit quality: $\chi^2/ndf={chi2/ndf:.2f}$"
    else:
        qtxt = r"Fit quality: n/a"
    ax_xsec.text(0.012, 0.86 if math.isfinite(gamma_flux) else 0.90, qtxt, transform=ax_xsec.transAxes, va="top", fontsize=10.0)

    ax_xsec.legend(loc="best", framealpha=0.96)

    # Panel 4: helicity placeholder (only when helicity info is available).
    if has_tlp and ax_hel is not None:
        ax_hel.set_title("Helicity contribution placeholder")
        ax_hel.set_xlabel(r"$\phi$ [rad]")
        ax_hel.set_ylabel(r"$d\sigma/(dt\,d\phi)$ [nb/GeV$^2$]")
        ax_hel.axhline(0.0, color="0.4", linestyle="--", linewidth=1.0)
        ax_hel.fill_between(phi_curve, tlp - etlp, tlp + etlp, color="#9467bd", alpha=0.12, linewidth=0.0, label=r"$\sigma_{TL'}\sin\phi \pm 1\sigma$")
        ax_hel.plot(phi_curve, tlp, color="#7c3aed", linewidth=2.2, linestyle=(0, (4, 2)), label=r"$\sigma_{TL'}\sin\phi$")
        ax_hel.text(
            0.03,
            0.96,
            "Helicity data placeholder active\n(use per-helicity yields for final asymmetry overlay)",
            transform=ax_hel.transAxes,
            va="top",
            fontsize=9.8,
        )
        ax_hel.legend(loc="best", framealpha=0.95)
    else:
        ax_xsec.text(
            0.99,
            0.02,
            "Helicity panel omitted (fit_asym_ok=0)",
            transform=ax_xsec.transAxes,
            ha="right",
            va="bottom",
            fontsize=9.2,
            color="#4b5563",
        )

    fig.suptitle(
        rf"Per-slice multipanel view | it={it}, iq={iq}, ix={ix}",
        fontsize=14,
        y=0.985,
    )

    fig.text(
        0.5,
        0.012,
        "Contribution/helicity bands use coefficient-error propagation in diagonal approximation.",
        ha="center",
        fontsize=9.7,
    )

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate publication-quality overall and per-slice cross-section plots in one PDF."
    )

    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    default_slice_csv = os.path.join(root_dir, "excl_xsec_pi0_analysis_no_simc_model_slice_summary.csv")
    default_phi_csv = os.path.join(root_dir, "excl_xsec_pi0_analysis_no_simc_model_summary.csv")
    default_pdf = os.path.join(root_dir, "output_pi0_xsec_no_simc_model", "sigma_coeff_publication_plots.pdf")

    parser.add_argument("--slice-csv", default=default_slice_csv, help="Path to slice summary CSV")
    parser.add_argument("--phi-csv", default=default_phi_csv, help="Path to per-phi summary CSV")
    parser.add_argument("--output-pdf", default=default_pdf, help="Path to output PDF")
    parser.add_argument(
        "--nb-scale",
        type=float,
        default=1e9,
        help="Scale factor to convert internal sigma units to nb/GeV^2 (default: 1e9)",
    )
    parser.add_argument(
        "--max-slices",
        type=int,
        default=0,
        help="If >0, limit per-slice pages for quick checks (default: 0 means all)",
    )

    return parser



# --- New function: multi-panel coefficients vs -t' for all (Q2, xB) bins ---
def draw_coeffs_overlay_grid_page(pdf: PdfPages, by_group: Dict[GroupKey, List[dict]], nb_scale: float) -> None:
    # Determine unique Q2 and xB bins
    group_keys = sorted(by_group.keys())
    q2_bins = sorted(set(iq for iq, ix in group_keys))
    xb_bins = sorted(set(ix for iq, ix in group_keys))
    n_q2 = len(q2_bins)
    n_xb = len(xb_bins)

    # Make the figure larger for clarity
    fig, axes = plt.subplots(n_q2, n_xb, figsize=(5.5*n_xb, 4.5*n_q2), sharex=True, sharey=True, squeeze=False)
    for i, iq in enumerate(q2_bins):
        for j, ix in enumerate(xb_bins):
            ax = axes[i][j]
            rows = by_group.get((iq, ix), [])
            if not rows:
                ax.axis("off")
                continue
            x = -np.array([r["tprime_center"] for r in rows], dtype=float)
            xerr = 0.5 * np.abs(
                np.array([r["tprime_hi"] for r in rows], dtype=float)
                - np.array([r["tprime_lo"] for r in rows], dtype=float)
            )
            # Determine if TL' is present in this group
            has_tlp = valid_tlp_group(rows)
            specs = COEFF_META if has_tlp else COEFF_META[:3]
            handles, labels_ = [], []
            for idx, (vkey, ekey, label, color, marker) in enumerate(specs):
                y = nb_scale * np.array([r[vkey] for r in rows], dtype=float)
                yerr = nb_scale * clean_err(np.array([r[ekey] for r in rows], dtype=float))
                h = ax.errorbar(
                    x, y, xerr=xerr, yerr=yerr,
                    fmt=marker, linestyle="none", color=color,
                    markerfacecolor="white", markeredgewidth=1.35, markersize=6.5, capsize=3.0, elinewidth=1.3,
                    label=label
                )
                handles.append(h)
                labels_.append(label)
            q2lo = rows[0]["q2_lo"]
            q2hi = rows[0]["q2_hi"]
            xblo = rows[0]["xb_lo"]
            xbhi = rows[0]["xb_hi"]
            ax.set_title(rf"$Q^2\in[{q2lo:.3f},{q2hi:.3f}]$, $x_B\in[{xblo:.3f},{xbhi:.3f}]$", fontsize=14)
            ax.axhline(0.0, color="0.7", linestyle="--", linewidth=1.1)
            ax.set_xlabel(r"$-t'$ [GeV$^2$]", fontsize=13)
            ax.set_ylabel(r"$d\sigma/dt$ [nb/GeV$^2$]", fontsize=13)
            ax.legend(loc="best", framealpha=0.97, fontsize=11)
            ax.tick_params(axis='both', which='major', labelsize=11)
    fig.suptitle(r"Fitted coefficients vs $-t'$ for all $(Q^2, x_B)$ bins", fontsize=18, y=0.995)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    base_style()
    args = build_parser().parse_args()

    if not os.path.exists(args.slice_csv):
        raise FileNotFoundError(f"Slice summary CSV not found: {args.slice_csv}")
    if not os.path.exists(args.phi_csv):
        raise FileNotFoundError(f"Per-phi summary CSV not found: {args.phi_csv}")

    by_group, by_slice = parse_slice_csv(args.slice_csv)
    if not by_group:
        raise RuntimeError("No fitted slices were found in the slice CSV (fit_xsec_ok == 1 required).")

    by_slice_phi = parse_phi_csv(args.phi_csv)
    if not by_slice_phi:
        raise RuntimeError("No per-phi slice entries were found in the per-phi CSV.")

    os.makedirs(os.path.dirname(os.path.abspath(args.output_pdf)), exist_ok=True)

    ordered_groups = sorted(by_group.keys())
    ordered_slice_keys = sorted(by_slice.keys(), key=lambda k: (k[1], k[2], k[0]))

    with PdfPages(args.output_pdf) as pdf:
        draw_title_page(
            pdf,
            os.path.abspath(args.slice_csv),
            os.path.abspath(args.phi_csv),
            len(ordered_groups),
            len(by_slice),
            len(by_slice_phi),
            args.nb_scale,
        )

        # --- New: add multi-panel grid page with all coefficients overlaid ---
        draw_coeffs_overlay_grid_page(pdf, by_group, args.nb_scale)

        for iq, ix in ordered_groups:
            draw_overall_group_page(pdf, by_group[(iq, ix)], iq, ix, args.nb_scale)

        n_written = 0
        for sk in ordered_slice_keys:
            if args.max_slices > 0 and n_written >= args.max_slices:
                break
            fit_row = by_slice[sk]
            phi_rows = by_slice_phi.get(sk, [])
            if not phi_rows:
                continue
            draw_per_slice_page(pdf, fit_row, phi_rows, args.nb_scale)
            n_written += 1

    print(f"[INFO] Wrote publication PDF: {args.output_pdf}")
    print(f"[INFO] Overall groups plotted: {len(ordered_groups)}")
    print(f"[INFO] Per-slice pages plotted: {n_written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
