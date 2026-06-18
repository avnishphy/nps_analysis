#!/usr/bin/env python3
"""Reproducible data/SIMC/smeared-SIMC comparison for the smearing pipeline.

The notebook remains useful for interactive inspection, but this script writes
plain CSV/PNG artifacts that can be archived with a production or reduced test
run.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np


M_PROTON = 0.9382720813


OBSERVABLES = {
    "mpi0": (0.128, 0.142, 50, "M_gg [GeV]"),
    "mmiss": (0.8, 1.0, 120, "M_miss [GeV]"),
    "mpgg2": (7.8, 12.2, 50, "(p_target + gamma gamma)^2 [GeV^2]"),
}

PEAK_WINDOWS = {
    "mpi0": (0.130, 0.140),
    "mmiss": (0.86, 0.98),
    "mpgg2": (9.0, 11.5),
}

CORRELATION_PAIRS = [
    ("mpi0", "mmiss"),
    ("mpi0", "mpgg2"),
]


DATA_BRANCH_CANDIDATES = [
    "mpi0", "mpi0_all", "mmiss", "mmiss_all",
    "cluster_e_1", "cluster_e_2", "clust_e_1", "clust_e_2",
    "pi0_weight", "weight", "scale",
    "is_exclusive_ellipse_combined", "is_exclusive_ellipse", "is_exclusive_mcd",
]

SIM_BRANCH_CANDIDATES = [
    "mpi0", "mpi0_all", "mmiss", "mmiss_all",
    "clust_E", "full_weight", "is_exclusive_ellipse",
]


class HistStats:
    def __init__(self, lo: float, hi: float, nbins: int, peak_lo: float, peak_hi: float):
        self.lo = lo
        self.hi = hi
        self.nbins = nbins
        self.peak_lo = peak_lo
        self.peak_hi = peak_hi
        self.counts = np.zeros(nbins, dtype=float)
        self.sumw = 0.0
        self.sumwx = 0.0
        self.sumwx2 = 0.0
        self.peak_sumw = 0.0
        self.peak_sumwx = 0.0
        self.peak_sumwx2 = 0.0

    def fill(self, value: float, weight: float) -> None:
        if not (np.isfinite(value) and np.isfinite(weight) and weight > 0.0):
            return
        self.sumw += weight
        self.sumwx += weight * value
        self.sumwx2 += weight * value * value
        if self.peak_lo <= value < self.peak_hi:
            self.peak_sumw += weight
            self.peak_sumwx += weight * value
            self.peak_sumwx2 += weight * value * value
        if self.lo <= value < self.hi:
            ibin = int((value - self.lo) / (self.hi - self.lo) * self.nbins)
            ibin = min(max(ibin, 0), self.nbins - 1)
            self.counts[ibin] += weight

    def normalized(self) -> np.ndarray:
        total = self.counts.sum()
        return self.counts / total if total > 0 else self.counts.copy()

    def mean_std(self) -> tuple[float, float]:
        return mean_std_from_sums(self.sumw, self.sumwx, self.sumwx2)

    def peak_mean_std(self) -> tuple[float, float]:
        return mean_std_from_sums(self.peak_sumw, self.peak_sumwx, self.peak_sumwx2)

    def tail_fraction(self) -> float:
        total = self.counts.sum()
        if total <= 0.0:
            return np.nan
        edges = np.linspace(self.lo, self.hi, self.nbins + 1)
        centers = 0.5 * (edges[1:] + edges[:-1])
        peak_mask = (centers >= self.peak_lo) & (centers < self.peak_hi)
        return float(self.counts[~peak_mask].sum() / total)


class Hist2DStats:
    def __init__(self, xname: str, yname: str):
        xlo, xhi, xbins, _ = OBSERVABLES[xname]
        ylo, yhi, ybins, _ = OBSERVABLES[yname]
        self.xname = xname
        self.yname = yname
        self.xlo = xlo
        self.xhi = xhi
        self.xbins = xbins
        self.ylo = ylo
        self.yhi = yhi
        self.ybins = ybins
        self.counts = np.zeros((xbins, ybins), dtype=float)

    def fill(self, xvalue: float, yvalue: float, weight: float) -> None:
        if not (np.isfinite(xvalue) and np.isfinite(yvalue) and np.isfinite(weight) and weight > 0.0):
            return
        if not (self.xlo <= xvalue < self.xhi and self.ylo <= yvalue < self.yhi):
            return
        ix = int((xvalue - self.xlo) / (self.xhi - self.xlo) * self.xbins)
        iy = int((yvalue - self.ylo) / (self.yhi - self.ylo) * self.ybins)
        ix = min(max(ix, 0), self.xbins - 1)
        iy = min(max(iy, 0), self.ybins - 1)
        self.counts[ix, iy] += weight

    def normalized(self) -> np.ndarray:
        total = self.counts.sum()
        return self.counts / total if total > 0 else self.counts.copy()


def mean_std_from_sums(sumw: float, sumwx: float, sumwx2: float) -> tuple[float, float]:
    if sumw <= 0.0:
        return np.nan, np.nan
    mean = sumwx / sumw
    var = max(sumwx2 / sumw - mean * mean, 0.0)
    return mean, float(np.sqrt(var))


def tail_fraction_from_hist(counts: np.ndarray, edges: np.ndarray, peak_lo: float, peak_hi: float) -> float:
    total = counts.sum()
    if total <= 0.0:
        return np.nan
    centers = 0.5 * (edges[1:] + edges[:-1])
    peak_mask = (centers >= peak_lo) & (centers < peak_hi)
    return float(counts[~peak_mask].sum() / total)


def root_branch(tree, names: list[str], required: bool = True) -> str | None:
    for name in names:
        if tree.GetBranch(name):
            return name
    if required:
        raise KeyError(f"Missing branch choices in {tree.GetName()}: {names}")
    return None


def root_value(tree, name: str):
    return getattr(tree, name)


def root_scalar(tree, name: str) -> float:
    return float(root_value(tree, name))


def root_array2_sum(tree, name: str) -> float:
    arr = root_value(tree, name)
    return float(arr[0]) + float(arr[1])


def process_root_tree(path: Path, tree_name: str, data_like: bool,
                      max_entries: int | None, use_sim_full_weight: bool):
    import ROOT

    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Cannot open ROOT file: {path}")
    tree = root_file.Get(tree_name)
    if not tree:
        raise KeyError(f"{tree_name!r} not found in {path}")

    if data_like:
        b_mpi0 = root_branch(tree, ["mpi0", "mpi0_all"])
        b_mmiss = root_branch(tree, ["mmiss", "mmiss_all"])
        b_e1 = root_branch(tree, ["cluster_e_1", "clust_e_1"])
        b_e2 = root_branch(tree, ["cluster_e_2", "clust_e_2"])
        b_w = root_branch(tree, ["pi0_weight", "weight"], required=False)
        b_scale = root_branch(tree, ["scale"], required=False)
        b_excl = root_branch(tree, ["is_exclusive_ellipse_combined", "is_exclusive_ellipse", "is_exclusive_mcd"], required=False)
    else:
        b_mpi0 = root_branch(tree, ["mpi0", "mpi0_all"])
        b_mmiss = root_branch(tree, ["mmiss", "mmiss_all"])
        b_clust_e = root_branch(tree, ["clust_E"])
        b_w = root_branch(tree, ["full_weight"], required=False)
        b_excl = root_branch(tree, ["is_exclusive_ellipse"], required=False)

    stats = {name: HistStats(lo, hi, nbins, *PEAK_WINDOWS[name])
             for name, (lo, hi, nbins, _label) in OBSERVABLES.items()}
    stats2d = {pair: Hist2DStats(*pair) for pair in CORRELATION_PAIRS}
    nentries = int(tree.GetEntries())
    nread = min(nentries, max_entries) if max_entries else nentries
    print(f"Reading {path}:{tree_name} entries={nread}/{nentries} backend=ROOT", flush=True)

    for i in range(nread):
        tree.GetEntry(i)
        mpi0 = root_scalar(tree, b_mpi0)
        mmiss = root_scalar(tree, b_mmiss)
        if data_like:
            esum = root_scalar(tree, b_e1) + root_scalar(tree, b_e2)
            weight = 1.0
            if b_w:
                weight *= root_scalar(tree, b_w)
            if b_scale:
                weight *= root_scalar(tree, b_scale)
            if b_excl:
                weight *= root_scalar(tree, b_excl)
        else:
            esum = root_array2_sum(tree, b_clust_e)
            weight = root_scalar(tree, b_w) if (use_sim_full_weight and b_w) else 1.0
            if b_excl:
                weight *= root_scalar(tree, b_excl)
        mpgg2 = M_PROTON * M_PROTON + mpi0 * mpi0 + 2.0 * M_PROTON * esum
        values = {"mpi0": mpi0, "mmiss": mmiss, "mpgg2": mpgg2}
        for obs_name, obs_value in values.items():
            stats[obs_name].fill(obs_value, weight)
        for pair, hist2d in stats2d.items():
            hist2d.fill(values[pair[0]], values[pair[1]], weight)

    root_file.Close()
    return stats, stats2d


def load_tree(path: Path, tree_name: str, branch_candidates: list[str], max_entries: int | None):
    import uproot

    with uproot.open(path) as root_file:
        if tree_name not in root_file:
            raise KeyError(f"{tree_name!r} not found in {path}")
        tree = root_file[tree_name]
        available = set(tree.keys())
        branches = [name for name in branch_candidates if name in available]
        if not branches:
            raise KeyError(f"No candidate comparison branches found in {path}:{tree_name}")
        return tree.arrays(branches, library="np", entry_stop=max_entries)


def first_existing(arrays: dict[str, np.ndarray], names: list[str]) -> np.ndarray:
    for name in names:
        if name in arrays:
            return arrays[name]
    raise KeyError(f"Missing all branch choices: {names}")


def optional_branch(arrays: dict[str, np.ndarray], names: list[str], default):
    for name in names:
        if name in arrays:
            return arrays[name]
    n = len(next(iter(arrays.values())))
    return np.full(n, default)


def first_two_energy_sum(arrays: dict[str, np.ndarray], data_like: bool) -> np.ndarray:
    if "clust_E" in arrays:
        e = arrays["clust_E"]
        return np.asarray(e[:, 0], dtype=float) + np.asarray(e[:, 1], dtype=float)
    if data_like:
        e1 = first_existing(arrays, ["cluster_e_1", "clust_e_1"])
        e2 = first_existing(arrays, ["cluster_e_2", "clust_e_2"])
        return np.asarray(e1, dtype=float) + np.asarray(e2, dtype=float)
    raise KeyError("Cannot construct photon energy sum")


def collect_observables(arrays: dict[str, np.ndarray], data_like: bool) -> dict[str, np.ndarray]:
    out = {}
    out["mpi0"] = np.asarray(first_existing(arrays, ["mpi0", "mpi0_all"]), dtype=float)
    out["mmiss"] = np.asarray(first_existing(arrays, ["mmiss", "mmiss_all"]), dtype=float)
    energy_sum = first_two_energy_sum(arrays, data_like)
    out["mpgg2"] = M_PROTON * M_PROTON + out["mpi0"] * out["mpi0"] + 2.0 * M_PROTON * energy_sum
    return out


def weights(arrays: dict[str, np.ndarray], data_like: bool, use_sim_full_weight: bool) -> np.ndarray:
    n = len(next(iter(arrays.values())))
    if data_like:
        base = optional_branch(arrays, ["pi0_weight", "weight"], 1.0).astype(float)
        scale = optional_branch(arrays, ["scale"], 1.0).astype(float)
        excl = optional_branch(arrays, ["is_exclusive_ellipse_combined", "is_exclusive_ellipse", "is_exclusive_mcd"], 1.0).astype(float)
        return base * scale * excl
    base = optional_branch(arrays, ["full_weight"], 1.0).astype(float) if use_sim_full_weight else np.ones(n)
    excl = optional_branch(arrays, ["is_exclusive_ellipse"], 1.0).astype(float)
    return base * excl


def normalized_hist(values, w, lo, hi, nbins):
    mask = np.isfinite(values) & np.isfinite(w) & (w > 0)
    counts, edges = np.histogram(values[mask], bins=nbins, range=(lo, hi), weights=w[mask])
    integral = counts.sum()
    if integral > 0:
        counts = counts / integral
    return counts, edges


def baker_cousins(sim, data, floor=1e-12):
    chi2 = 0.0
    empty = 0
    for s, d in zip(sim, data):
        if d > 0:
            if s <= 0:
                s = floor
                empty += 1
            chi2 += 2.0 * (s - d + d * np.log(d / s))
        elif s > 0:
            chi2 += 2.0 * s
    return chi2, empty


def weighted_mean_std(values, w):
    mask = np.isfinite(values) & np.isfinite(w) & (w > 0)
    values = values[mask]
    w = w[mask]
    if values.size == 0 or w.sum() <= 0:
        return np.nan, np.nan
    mean = np.average(values, weights=w)
    var = np.average((values - mean) ** 2, weights=w)
    return mean, np.sqrt(max(var, 0.0))


def weighted_mean_std_window(values, w, lo: float, hi: float):
    mask = np.isfinite(values) & np.isfinite(w) & (w > 0) & (values >= lo) & (values < hi)
    return weighted_mean_std(values[mask], w[mask])


def normalized_hist2d(xvalues, yvalues, w, xname: str, yname: str):
    xlo, xhi, xbins, _ = OBSERVABLES[xname]
    ylo, yhi, ybins, _ = OBSERVABLES[yname]
    mask = np.isfinite(xvalues) & np.isfinite(yvalues) & np.isfinite(w) & (w > 0)
    counts, xedges, yedges = np.histogram2d(
        xvalues[mask], yvalues[mask],
        bins=(xbins, ybins),
        range=((xlo, xhi), (ylo, yhi)),
        weights=w[mask],
    )
    total = counts.sum()
    if total > 0:
        counts = counts / total
    return counts, xedges, yedges


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--sim", required=True, type=Path)
    parser.add_argument("--smeared", required=True, type=Path)
    parser.add_argument("--data-tree", default="physics")
    parser.add_argument("--sim-tree", default="simulation")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--use-sim-full-weight", action="store_true")
    parser.add_argument("--max-entries", type=int, default=None,
                        help="Optional entry limit for quick validation plots.")
    parser.add_argument("--backend", choices=["auto", "root", "uproot"], default="auto")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    use_root_backend = args.backend in ("auto", "root")
    if use_root_backend:
        try:
            data_stats, data_stats2d = process_root_tree(args.data, args.data_tree, True, args.max_entries, False)
            sim_stats, sim_stats2d = process_root_tree(args.sim, args.sim_tree, False, args.max_entries, args.use_sim_full_weight)
            smeared_stats, smeared_stats2d = process_root_tree(args.smeared, args.sim_tree, False, args.max_entries, args.use_sim_full_weight)
        except Exception:
            if args.backend == "root":
                raise
            print("ROOT backend failed; falling back to uproot backend", flush=True)
            use_root_backend = False

    if not use_root_backend:
        data_arrays = load_tree(args.data, args.data_tree, DATA_BRANCH_CANDIDATES, args.max_entries)
        sim_arrays = load_tree(args.sim, args.sim_tree, SIM_BRANCH_CANDIDATES, args.max_entries)
        smeared_arrays = load_tree(args.smeared, args.sim_tree, SIM_BRANCH_CANDIDATES, args.max_entries)

        data_obs = collect_observables(data_arrays, data_like=True)
        sim_obs = collect_observables(sim_arrays, data_like=False)
        smeared_obs = collect_observables(smeared_arrays, data_like=False)

        w_data = weights(data_arrays, data_like=True, use_sim_full_weight=False)
        w_sim = weights(sim_arrays, data_like=False, use_sim_full_weight=args.use_sim_full_weight)
        w_smeared = weights(smeared_arrays, data_like=False, use_sim_full_weight=args.use_sim_full_weight)

    rows = []
    for name, (lo, hi, nbins, _label) in OBSERVABLES.items():
        peak_lo, peak_hi = PEAK_WINDOWS[name]
        if use_root_backend:
            h_data = data_stats[name].normalized()
            h_sim = sim_stats[name].normalized()
            h_smeared = smeared_stats[name].normalized()
            edges = np.linspace(lo, hi, nbins + 1)
        else:
            h_data, edges = normalized_hist(data_obs[name], w_data, lo, hi, nbins)
            h_sim, _ = normalized_hist(sim_obs[name], w_sim, lo, hi, nbins)
            h_smeared, _ = normalized_hist(smeared_obs[name], w_smeared, lo, hi, nbins)
        chi2_unsm, empty_unsm = baker_cousins(h_sim, h_data)
        chi2_sm, empty_sm = baker_cousins(h_smeared, h_data)
        if use_root_backend:
            mean_data, std_data = data_stats[name].mean_std()
            mean_sim, std_sim = sim_stats[name].mean_std()
            mean_smeared, std_smeared = smeared_stats[name].mean_std()
            peak_mean_data, peak_std_data = data_stats[name].peak_mean_std()
            peak_mean_sim, peak_std_sim = sim_stats[name].peak_mean_std()
            peak_mean_smeared, peak_std_smeared = smeared_stats[name].peak_mean_std()
            tail_data = data_stats[name].tail_fraction()
            tail_sim = sim_stats[name].tail_fraction()
            tail_smeared = smeared_stats[name].tail_fraction()
        else:
            mean_data, std_data = weighted_mean_std(data_obs[name], w_data)
            mean_sim, std_sim = weighted_mean_std(sim_obs[name], w_sim)
            mean_smeared, std_smeared = weighted_mean_std(smeared_obs[name], w_smeared)
            peak_mean_data, peak_std_data = weighted_mean_std_window(data_obs[name], w_data, peak_lo, peak_hi)
            peak_mean_sim, peak_std_sim = weighted_mean_std_window(sim_obs[name], w_sim, peak_lo, peak_hi)
            peak_mean_smeared, peak_std_smeared = weighted_mean_std_window(smeared_obs[name], w_smeared, peak_lo, peak_hi)
            tail_data = tail_fraction_from_hist(h_data, edges, peak_lo, peak_hi)
            tail_sim = tail_fraction_from_hist(h_sim, edges, peak_lo, peak_hi)
            tail_smeared = tail_fraction_from_hist(h_smeared, edges, peak_lo, peak_hi)
        mean_abs_delta_unsm = abs(mean_sim - mean_data)
        mean_abs_delta_sm = abs(mean_smeared - mean_data)
        std_abs_delta_unsm = abs(std_sim - std_data)
        std_abs_delta_sm = abs(std_smeared - std_data)
        peak_mean_abs_delta_unsm = abs(peak_mean_sim - peak_mean_data)
        peak_mean_abs_delta_sm = abs(peak_mean_smeared - peak_mean_data)
        peak_std_abs_delta_unsm = abs(peak_std_sim - peak_std_data)
        peak_std_abs_delta_sm = abs(peak_std_smeared - peak_std_data)
        tail_abs_delta_unsm = abs(tail_sim - tail_data)
        tail_abs_delta_sm = abs(tail_smeared - tail_data)
        rows.append({
            "observable": name,
            "fit_window_lo": lo,
            "fit_window_hi": hi,
            "peak_window_lo": peak_lo,
            "peak_window_hi": peak_hi,
            "chi2_unsmeared": chi2_unsm,
            "chi2_smeared": chi2_sm,
            "delta_chi2_smeared_minus_unsmeared": chi2_sm - chi2_unsm,
            "shape_improved": int(chi2_sm < chi2_unsm),
            "empty_bins_unsmeared": empty_unsm,
            "empty_bins_smeared": empty_sm,
            "mean_data": mean_data,
            "mean_unsmeared": mean_sim,
            "mean_smeared": mean_smeared,
            "mean_abs_delta_unsmeared": mean_abs_delta_unsm,
            "mean_abs_delta_smeared": mean_abs_delta_sm,
            "mean_position_improved": int(mean_abs_delta_sm < mean_abs_delta_unsm),
            "std_data": std_data,
            "std_unsmeared": std_sim,
            "std_smeared": std_smeared,
            "std_abs_delta_unsmeared": std_abs_delta_unsm,
            "std_abs_delta_smeared": std_abs_delta_sm,
            "width_improved": int(std_abs_delta_sm < std_abs_delta_unsm),
            "peak_mean_data": peak_mean_data,
            "peak_mean_unsmeared": peak_mean_sim,
            "peak_mean_smeared": peak_mean_smeared,
            "peak_mean_abs_delta_unsmeared": peak_mean_abs_delta_unsm,
            "peak_mean_abs_delta_smeared": peak_mean_abs_delta_sm,
            "peak_position_improved": int(peak_mean_abs_delta_sm < peak_mean_abs_delta_unsm),
            "peak_std_data": peak_std_data,
            "peak_std_unsmeared": peak_std_sim,
            "peak_std_smeared": peak_std_smeared,
            "peak_std_abs_delta_unsmeared": peak_std_abs_delta_unsm,
            "peak_std_abs_delta_smeared": peak_std_abs_delta_sm,
            "peak_width_improved": int(peak_std_abs_delta_sm < peak_std_abs_delta_unsm),
            "tail_fraction_data": tail_data,
            "tail_fraction_unsmeared": tail_sim,
            "tail_fraction_smeared": tail_smeared,
            "tail_fraction_abs_delta_unsmeared": tail_abs_delta_unsm,
            "tail_fraction_abs_delta_smeared": tail_abs_delta_sm,
            "tail_fraction_improved": int(tail_abs_delta_sm < tail_abs_delta_unsm),
        })

    csv_path = args.out_dir / "smearing_comparison_metrics.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {csv_path}", flush=True)

    rows2d = []
    h2d_cache = {}
    for xname, yname in CORRELATION_PAIRS:
        if use_root_backend:
            h_data_2d = data_stats2d[(xname, yname)].normalized()
            h_sim_2d = sim_stats2d[(xname, yname)].normalized()
            h_smeared_2d = smeared_stats2d[(xname, yname)].normalized()
            xlo, xhi, xbins, _ = OBSERVABLES[xname]
            ylo, yhi, ybins, _ = OBSERVABLES[yname]
            xedges = np.linspace(xlo, xhi, xbins + 1)
            yedges = np.linspace(ylo, yhi, ybins + 1)
        else:
            h_data_2d, xedges, yedges = normalized_hist2d(data_obs[xname], data_obs[yname], w_data, xname, yname)
            h_sim_2d, _, _ = normalized_hist2d(sim_obs[xname], sim_obs[yname], w_sim, xname, yname)
            h_smeared_2d, _, _ = normalized_hist2d(smeared_obs[xname], smeared_obs[yname], w_smeared, xname, yname)
        chi2_unsm_2d, empty_unsm_2d = baker_cousins(h_sim_2d.ravel(), h_data_2d.ravel())
        chi2_sm_2d, empty_sm_2d = baker_cousins(h_smeared_2d.ravel(), h_data_2d.ravel())
        rows2d.append({
            "pair": f"{xname}_vs_{yname}",
            "chi2_unsmeared": chi2_unsm_2d,
            "chi2_smeared": chi2_sm_2d,
            "delta_chi2_smeared_minus_unsmeared": chi2_sm_2d - chi2_unsm_2d,
            "shape_improved": int(chi2_sm_2d < chi2_unsm_2d),
            "empty_bins_unsmeared": empty_unsm_2d,
            "empty_bins_smeared": empty_sm_2d,
        })
        h2d_cache[(xname, yname)] = (h_data_2d, h_sim_2d, h_smeared_2d, xedges, yedges)

    csv2d_path = args.out_dir / "smearing_2d_metrics.csv"
    with csv2d_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows2d[0].keys()))
        writer.writeheader()
        writer.writerows(rows2d)
    print(f"Wrote {csv2d_path}", flush=True)

    os.environ.setdefault("MPLCONFIGDIR", str(args.out_dir / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(OBSERVABLES), 2, figsize=(11, 10), constrained_layout=True)
    for row, (name, (lo, hi, nbins, label)) in enumerate(OBSERVABLES.items()):
        if use_root_backend:
            h_data = data_stats[name].normalized()
            h_sim = sim_stats[name].normalized()
            h_smeared = smeared_stats[name].normalized()
            edges = np.linspace(lo, hi, nbins + 1)
        else:
            h_data, edges = normalized_hist(data_obs[name], w_data, lo, hi, nbins)
            h_sim, _ = normalized_hist(sim_obs[name], w_sim, lo, hi, nbins)
            h_smeared, _ = normalized_hist(smeared_obs[name], w_smeared, lo, hi, nbins)
        centers = 0.5 * (edges[1:] + edges[:-1])
        ax = axes[row, 0]
        ax.stairs(h_data, edges, color="black", linewidth=1.6, label="Data")
        ax.stairs(h_sim, edges, color="tab:blue", linestyle="--", linewidth=1.4, label="SIMC")
        ax.stairs(h_smeared, edges, color="tab:red", linewidth=1.4, label="Smeared")
        ax.set_xlabel(label)
        ax.set_ylabel("Area normalized")
        ax.legend(fontsize=8)
        ratio_ax = axes[row, 1]
        ratio_sim = np.divide(h_sim, h_data, out=np.full_like(h_data, np.nan), where=h_data > 0)
        ratio_smeared = np.divide(h_smeared, h_data, out=np.full_like(h_data, np.nan), where=h_data > 0)
        ratio_ax.plot(centers, ratio_sim, "o-", color="tab:blue", markersize=2.5, label="SIMC/Data")
        ratio_ax.plot(centers, ratio_smeared, "s-", color="tab:red", markersize=2.5, label="Smeared/Data")
        ratio_ax.axhline(1.0, color="black", linewidth=1.0)
        ratio_ax.set_xlabel(label)
        ratio_ax.set_ylabel("Ratio")
        ratio_ax.set_ylim(0.0, 2.5)
        ratio_ax.legend(fontsize=8)
    plot_path = args.out_dir / "smearing_comparison.png"
    fig.savefig(plot_path, dpi=160)

    print(f"Wrote {plot_path}", flush=True)

    fig_zoom, zoom_axes = plt.subplots(len(OBSERVABLES), 2, figsize=(11, 10), constrained_layout=True)
    for row, (name, (lo, hi, nbins, label)) in enumerate(OBSERVABLES.items()):
        peak_lo, peak_hi = PEAK_WINDOWS[name]
        if use_root_backend:
            h_data = data_stats[name].normalized()
            h_sim = sim_stats[name].normalized()
            h_smeared = smeared_stats[name].normalized()
            edges = np.linspace(lo, hi, nbins + 1)
        else:
            h_data, edges = normalized_hist(data_obs[name], w_data, lo, hi, nbins)
            h_sim, _ = normalized_hist(sim_obs[name], w_sim, lo, hi, nbins)
            h_smeared, _ = normalized_hist(smeared_obs[name], w_smeared, lo, hi, nbins)
        centers = 0.5 * (edges[1:] + edges[:-1])
        ax = zoom_axes[row, 0]
        ax.stairs(h_data, edges, color="black", linewidth=1.6, label="Data")
        ax.stairs(h_sim, edges, color="tab:blue", linestyle="--", linewidth=1.4, label="SIMC")
        ax.stairs(h_smeared, edges, color="tab:red", linewidth=1.4, label="Smeared")
        ax.axvspan(peak_lo, peak_hi, color="0.9", zorder=-1)
        ax.set_xlim(peak_lo, peak_hi)
        ax.set_xlabel(label)
        ax.set_ylabel("Area normalized")
        ax.legend(fontsize=8)
        resid_ax = zoom_axes[row, 1]
        resid_ax.plot(centers, h_sim - h_data, "o-", color="tab:blue", markersize=2.5, label="SIMC - Data")
        resid_ax.plot(centers, h_smeared - h_data, "s-", color="tab:red", markersize=2.5, label="Smeared - Data")
        resid_ax.axhline(0.0, color="black", linewidth=1.0)
        resid_ax.axvspan(peak_lo, peak_hi, color="0.9", zorder=-1)
        resid_ax.set_xlim(peak_lo, peak_hi)
        resid_ax.set_xlabel(label)
        resid_ax.set_ylabel("Residual")
        resid_ax.legend(fontsize=8)
    zoom_path = args.out_dir / "smearing_comparison_peak_zoom.png"
    fig_zoom.savefig(zoom_path, dpi=160)
    print(f"Wrote {zoom_path}", flush=True)

    fig2d, axes2d = plt.subplots(len(CORRELATION_PAIRS), 3, figsize=(12, 7), constrained_layout=True)
    if len(CORRELATION_PAIRS) == 1:
        axes2d = np.asarray([axes2d])
    for row, (xname, yname) in enumerate(CORRELATION_PAIRS):
        h_data_2d, h_sim_2d, h_smeared_2d, xedges, yedges = h2d_cache[(xname, yname)]
        vmax = max(h_data_2d.max(), h_sim_2d.max(), h_smeared_2d.max(), 1e-12)
        labels = ("Data", "SIMC", "Smeared")
        hists = (h_data_2d, h_sim_2d, h_smeared_2d)
        for col, (title, hist) in enumerate(zip(labels, hists)):
            ax = axes2d[row, col]
            mesh = ax.pcolormesh(xedges, yedges, hist.T, shading="auto", vmin=0.0, vmax=vmax)
            ax.set_title(f"{title}: {xname} vs {yname}")
            ax.set_xlabel(OBSERVABLES[xname][3])
            ax.set_ylabel(OBSERVABLES[yname][3])
            fig2d.colorbar(mesh, ax=ax, fraction=0.046, pad=0.04)
    plot2d_path = args.out_dir / "smearing_2d_correlations.png"
    fig2d.savefig(plot2d_path, dpi=160)
    print(f"Wrote {plot2d_path}", flush=True)


if __name__ == "__main__":
    main()
