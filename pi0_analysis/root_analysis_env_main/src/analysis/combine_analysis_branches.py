#!/usr/bin/env python3
"""
combine_analysis_branches.py

Combine per-run diagnostics ROOT trees into one per-kinematics ROOT file.

This script is designed as a downstream stage of the unified analysis chain:
1) per-run diagnostics ROOT files are produced by nps_analysis_main.C
2) this combiner merges those run files for one Kin_old setting
3) optional QA plots are generated from the merged dataset

The script resolves canonical paths from the same workflow inputs used by
the main driver:
- NPS_KIN / --kin
- NPS_CONFIG_CSV / --config
- NPS_OUTPUT_BASE / --output-base
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from matplotlib.backends.backend_pdf import PdfPages


DEFAULT_TARGET = "LH2"
DEFAULT_CONFIG = "nps_dvcs_all_kins_main.csv"
BRANCHES_TO_EXCLUDE = {"event_id"}

# Debug-mode hard overrides requested by workflow integration task.
DEBUG_LIVETIME_OVERRIDE = 1.0
DEBUG_EFFICIENCY_OVERRIDE = 1.0

# Combined-level 2D mass-cut branches. Recomputed from the merged weighted
# mmiss_all:mpi0_all distribution, not copied from per-run flags.
CREATE_COMBINED_2D_MASS_CUT = True
COMBINED_MASS_CUT_TAG = "combined_2d_mass_cut"
WRITE_COMBINED_MASS_CUT_PARAM_TREES = False

MASS_CUT_CONFIG = {
    "n_mpi0_bins": 160,
    "mpi0_min": 0.11,
    "mpi0_max": 0.15,
    "n_mmiss_bins": 200,
    "mmiss_min": 0.6,
    "mmiss_max": 1.5,
    "seed_mmiss": 0.938,
    "seed_mmiss_half_width": 0.03,
    "peak_fraction": -1.0,
    "ellipse_core_quantile": 0.995,
    "ellipse_padding": 1.05,
    "auto_peak_min": 0.05,
    "auto_peak_max": 0.60,
    "auto_peak_step": 0.005,
    "auto_min_core_total_fraction": 0.05,
    "auto_min_core_bins": 30,
    "mcd_keep_total_fraction": 0.30,
    "mcd_iterations": 8,
    "mcd_ellipse_quantile": 0.995,
    "mcd_padding": 1.05,
}


@dataclass(frozen=True)
class WorkflowConfig:
    cfg_path: Path
    output_base: Path
    root_dir: Path
    out_combined_root: Path
    kin_setting: str
    target_to_combine: str
    efficiency_csv: Path
    analysis_plots_pdf: Path
    fp_debug_pdf: Path
    allowed_types: Tuple[str, ...]
    run_filter: Tuple[int, ...]
    create_analysis_plots: bool
    create_fp_debug_plots: bool


@dataclass(frozen=True)
class RunConfig:
    prescale_token: str
    prescale_value: float
    target: str
    charge_uC_fallback: Optional[float]


@dataclass(frozen=True)
class RunEfficiencyMeta:
    ps_value: float
    charge_uC: Optional[float]
    livetime_raw: float
    efficiency_raw: float
    livetime_used: float
    efficiency_used: float


@dataclass(frozen=True)
class Hist1DSpec:
    name: str
    title: str
    xlabel: str
    bins: int
    value_range: Optional[Tuple[float, float]]
    color: str


def sanitize_token(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", value.strip())


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    out = str(value).strip()
    if len(out) >= 2 and out[0] == out[-1] and out[0] in {'"', "'"}:
        out = out[1:-1]
    return out.strip()


def safe_int(value: Any) -> Optional[int]:
    try:
        text = clean_text(value)
        if text == "":
            return None
        return int(float(text))
    except Exception:
        return None


def safe_float(value: Any) -> Optional[float]:
    try:
        text = clean_text(value)
        if text == "":
            return None
        out = float(text)
        if not np.isfinite(out):
            return None
        return out
    except Exception:
        return None


def normalize_header(header: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", header.lower().strip())


def pick_column(columns: Iterable[str],
                preferred: Sequence[str],
                contains: Optional[str] = None) -> Optional[str]:
    cols = list(columns)
    norm_map = {normalize_header(c): c for c in cols}

    for key in preferred:
        found = norm_map.get(normalize_header(key))
        if found is not None:
            return found

    if contains:
        token = normalize_header(contains)
        for c in cols:
            if token in normalize_header(c):
                return c

    return None


def extract_prescale_r(token: str) -> Optional[int]:
    if not isinstance(token, str):
        return None
    m = re.search(r"ps\d*\s*=\s*(\d+)", token.strip(), flags=re.IGNORECASE)
    if m:
        return int(m.group(1))
    if "=" in token:
        rhs = token.split("=")[-1].strip()
        try:
            return int(rhs)
        except Exception:
            return None
    return None


def prescale_from_token(token: str) -> float:
    r = extract_prescale_r(token)
    if r is None or r <= 0:
        return 1.0
    return float((2 ** (r - 1)) + 1)


def parse_types_csv(types_csv: str) -> Tuple[str, ...]:
    out: List[str] = []
    for tok in types_csv.split(","):
        cleaned = clean_text(tok)
        if cleaned:
            out.append(cleaned)
    return tuple(out)


def parse_run_args(run_args: Sequence[Any]) -> Tuple[int, ...]:
    run_filter: List[int] = []
    for group in run_args:
        items = [group] if isinstance(group, str) else group
        for item in items:
            for token in str(item).split(","):
                cleaned = clean_text(token)
                if not cleaned:
                    continue
                parsed = safe_int(cleaned)
                if parsed is None:
                    raise ValueError(f"Invalid --run value: {cleaned}")
                run_filter.append(parsed)
    return tuple(sorted(set(run_filter)))


def infer_plots_dir(args: argparse.Namespace,
                    output_base: Path,
                    kin_safe: str,
                    root_dir: Path,
                    out_combined_root: Path) -> Path:
    if args.root_dir:
        return root_dir.parent / "plots" if root_dir.name.lower() == "root" else root_dir / "plots"
    if args.out_combined_root and out_combined_root.parent.name.lower() == "root":
        return out_combined_root.parent.parent / "plots"
    return output_base / kin_safe / "plots"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Combine per-run diagnostics trees into one per-kinematics ROOT file."
    )
    parser.add_argument("--kin", default=None,
                        help="Kin_old setting to combine (defaults to NPS_KIN env var).")
    parser.add_argument("--config", default=None,
                        help="Path to config CSV (defaults to NPS_CONFIG_CSV or canonical config).")
    parser.add_argument("--output-base", default=None,
                        help="Canonical output base directory (defaults to NPS_OUTPUT_BASE or repo/output).")
    parser.add_argument("--root-dir", default=None,
                        help="Optional explicit ROOT input directory (default: <output-base>/<kin>/root).")
    parser.add_argument("--out-combined-root", default=None,
                        help="Optional explicit output ROOT file path.")
    parser.add_argument("--efficiency-csv", default=None,
                        help="Optional explicit efficiency CSV path.")
    parser.add_argument("--target", default=DEFAULT_TARGET,
                        help="Target to combine (default: LH2).")
    parser.add_argument("--types", default="production,Production",
                        help="Allowed Type values in config CSV (comma-separated).")
    parser.add_argument("--run", action="append", nargs="+", default=[],
                        help="Restrict to one or more run numbers (repeatable; accepts space- or comma-separated values).")
    parser.add_argument("--no-analysis-plots", action="store_true",
                        help="Skip creation of analysis plot PDF.")
    parser.add_argument("--fp-debug-plots", action="store_true",
                        help="Enable focal-plane debug PDF.")
    return parser


def resolve_workflow_config(args: argparse.Namespace) -> WorkflowConfig:
    repo_root = Path(__file__).resolve().parents[2]

    kin_setting = clean_text(args.kin or os.getenv("NPS_KIN") or "")
    if not kin_setting:
        raise ValueError("KIN_SETTING is required (pass --kin or set NPS_KIN).")

    cfg_default = repo_root / "config" / DEFAULT_CONFIG
    output_default = repo_root / "output"

    cfg_path = Path(args.config or os.getenv("NPS_CONFIG_CSV") or cfg_default).expanduser().resolve()
    output_base = Path(args.output_base or os.getenv("NPS_OUTPUT_BASE") or output_default).expanduser().resolve()

    kin_safe = sanitize_token(kin_setting)
    root_dir = Path(args.root_dir).expanduser().resolve() if args.root_dir else (output_base / kin_safe / "root")

    target = clean_text(args.target or DEFAULT_TARGET)
    if not target:
        raise ValueError("target cannot be empty")
    target_safe = sanitize_token(target)

    out_combined_root = (
        Path(args.out_combined_root).expanduser().resolve()
        if args.out_combined_root
        else (root_dir / f"combined_branches_{target_safe}.root")
    )

    eff_candidates: List[Path] = []
    if args.efficiency_csv:
        eff_candidates.append(Path(args.efficiency_csv).expanduser().resolve())
    else:
        eff_base = output_base / "efficiency_stuff"
        eff_candidates.append((eff_base / f"efficiency_{kin_setting}.csv").resolve())
        eff_candidates.append((eff_base / f"efficiency_{kin_safe}.csv").resolve())

    efficiency_csv = next((p for p in eff_candidates if p.exists()), eff_candidates[0])

    plots_dir = infer_plots_dir(args, output_base, kin_safe, root_dir, out_combined_root)
    analysis_plots_pdf = plots_dir / f"{out_combined_root.stem}_plots.pdf"
    fp_debug_pdf = plots_dir / f"{out_combined_root.stem}_focal_plane_debug.pdf"

    allowed_types = parse_types_csv(args.types)
    if not allowed_types:
        raise ValueError("--types resolved to an empty set.")

    run_filter = parse_run_args(args.run)

    if not cfg_path.exists():
        raise FileNotFoundError(f"CFG_PATH does not exist: {cfg_path}")
    if not output_base.exists():
        raise FileNotFoundError(f"OUTPUT_BASE does not exist: {output_base}")
    if not root_dir.exists():
        raise FileNotFoundError(f"ROOT_DIR does not exist: {root_dir}")
    if not efficiency_csv.exists():
        candidates = "\n".join(str(p) for p in eff_candidates)
        raise FileNotFoundError(
            "Per-kin efficiency CSV not found. Tried:\n" + candidates
        )

    out_combined_root.parent.mkdir(parents=True, exist_ok=True)
    analysis_plots_pdf.parent.mkdir(parents=True, exist_ok=True)
    fp_debug_pdf.parent.mkdir(parents=True, exist_ok=True)

    return WorkflowConfig(
        cfg_path=cfg_path,
        output_base=output_base,
        root_dir=root_dir,
        out_combined_root=out_combined_root,
        kin_setting=kin_setting,
        target_to_combine=target,
        efficiency_csv=efficiency_csv,
        analysis_plots_pdf=analysis_plots_pdf,
        fp_debug_pdf=fp_debug_pdf,
        allowed_types=allowed_types,
        run_filter=run_filter,
        create_analysis_plots=(not args.no_analysis_plots),
        create_fp_debug_plots=bool(args.fp_debug_plots),
    )


def load_config(cfg_path: Path) -> pd.DataFrame:
    df = pd.read_csv(cfg_path, dtype=str, keep_default_na=False, skipinitialspace=True)
    df = df.replace({"": pd.NA})
    return df


def build_lookup(df_cfg: pd.DataFrame,
                 kin_setting: str,
                 allowed_types: Sequence[str]) -> Dict[int, RunConfig]:
    kin_col = pick_column(df_cfg.columns, ["Kin_old", "kin_old", "kinematic_setting"], contains="kin")
    if kin_col is None:
        raise ValueError("Config CSV missing Kin_old/kinematic-setting column.")

    run_col = pick_column(df_cfg.columns, ["run_number", "run"], contains="run")
    if run_col is None:
        raise ValueError("Config CSV missing run number column.")

    ps_col = pick_column(df_cfg.columns, ["prescale", "prescale_token"], contains="prescale")
    if ps_col is None:
        raise ValueError("Config CSV missing prescale column.")

    type_col = pick_column(df_cfg.columns, ["Type", "run_type"], contains="type")
    target_col = pick_column(df_cfg.columns, ["target"], contains="target")
    charge_mc_col = pick_column(df_cfg.columns,
                                ["charge(mC)", "accumulated_charge(mC)", "charge_mc"],
                                contains="charge")

    initial_rows = len(df_cfg)
    df_cfg = df_cfg[df_cfg[kin_col].map(clean_text) == kin_setting]
    print(f"[INFO] Filtered config by Kin_old='{kin_setting}': {initial_rows} -> {len(df_cfg)} rows")

    if type_col is not None:
        allowed_lower = {clean_text(t).lower() for t in allowed_types}
        before = len(df_cfg)
        df_cfg = df_cfg[df_cfg[type_col].map(lambda x: clean_text(x).lower() in allowed_lower)]
        print(f"[INFO] Filtered config by Type in {sorted(allowed_lower)}: {before} -> {len(df_cfg)} rows")

    if df_cfg.empty:
        raise ValueError(f"No config rows selected for kin='{kin_setting}' after type filtering.")

    lookup: Dict[int, RunConfig] = {}
    for _, row in df_cfg.iterrows():
        run = safe_int(row.get(run_col))
        if run is None:
            continue

        token = clean_text(row.get(ps_col))
        ps_value = prescale_from_token(token)

        target = clean_text(row.get(target_col)) if target_col else ""

        charge_uC_fallback: Optional[float] = None
        if charge_mc_col is not None:
            charge_mc = safe_float(row.get(charge_mc_col))
            if charge_mc is not None and charge_mc > 0:
                charge_uC_fallback = float(charge_mc * 1000.0)

        lookup[int(run)] = RunConfig(
            prescale_token=token,
            prescale_value=float(ps_value),
            target=target,
            charge_uC_fallback=charge_uC_fallback,
        )

    if not lookup:
        raise ValueError("No valid runs found after config parsing.")

    return lookup


def load_efficiency_metadata(efficiency_csv: Path,
                             kin_setting: str) -> Dict[int, RunEfficiencyMeta]:
    df = pd.read_csv(efficiency_csv, dtype=str, keep_default_na=False, skipinitialspace=True)
    df.columns = [str(c).strip() for c in df.columns]

    run_col = pick_column(df.columns, ["run_number", "run"], contains="run")
    if run_col is None:
        raise ValueError(f"Efficiency CSV missing run column: {efficiency_csv}")

    kin_col = pick_column(df.columns, ["kinematic_setting", "Kin_old"], contains="kin")
    ps_factor_col = pick_column(df.columns, ["ps_factor", "ps_value"], contains="psfactor")
    ps_token_col = pick_column(df.columns, ["prescale_token", "prescale"], contains="prescale")
    charge_after_col = pick_column(df.columns, ["HEL_charge_after_cut_uC", "hel_charge_after_cut_uc"],
                                   contains="chargeafter")
    charge_before_col = pick_column(df.columns, ["HEL_charge_before_cut_uC", "hel_charge_before_cut_uc"],
                                    contains="chargebefore")
    livetime_col = pick_column(df.columns, ["NewGen_EDTM_livetime", "livetime"], contains="livetime")
    efficiency_col = pick_column(df.columns,
                                 ["efficiency", "HMS_tracking_eff", "HMS_hodo_3of4_eff"],
                                 contains="eff")

    if kin_col is not None:
        before = len(df)
        df = df[df[kin_col].map(clean_text) == kin_setting]
        print(f"[INFO] Filtered efficiency CSV by kinematic_setting='{kin_setting}': {before} -> {len(df)} rows")

    out: Dict[int, RunEfficiencyMeta] = {}
    for _, row in df.iterrows():
        run = safe_int(row.get(run_col))
        if run is None:
            continue

        ps_value = safe_float(row.get(ps_factor_col)) if ps_factor_col else None
        if ps_value is None and ps_token_col is not None:
            ps_value = prescale_from_token(clean_text(row.get(ps_token_col)))
        if ps_value is None or ps_value <= 0:
            ps_value = 1.0

        charge_uC = None
        if charge_after_col is not None:
            charge_uC = safe_float(row.get(charge_after_col))
        if (charge_uC is None or charge_uC <= 0) and charge_before_col is not None:
            charge_uC = safe_float(row.get(charge_before_col))
        if charge_uC is not None and charge_uC <= 0:
            charge_uC = None

        livetime_raw = safe_float(row.get(livetime_col)) if livetime_col else None
        efficiency_raw = safe_float(row.get(efficiency_col)) if efficiency_col else None

        out[int(run)] = RunEfficiencyMeta(
            ps_value=float(ps_value),
            charge_uC=charge_uC,
            livetime_raw=float(livetime_raw) if livetime_raw is not None else 1.0,
            efficiency_raw=float(efficiency_raw) if efficiency_raw is not None else 1.0,
            livetime_used=float(DEBUG_LIVETIME_OVERRIDE),
            efficiency_used=float(DEBUG_EFFICIENCY_OVERRIDE),
        )

    if not out:
        raise ValueError(f"No usable rows found in efficiency CSV: {efficiency_csv}")

    print(
        f"[INFO] Loaded efficiency metadata for {len(out)} runs from {efficiency_csv}\n"
        f"       Debug override in effect: livetime={DEBUG_LIVETIME_OVERRIDE}, "
        f"efficiency={DEBUG_EFFICIENCY_OVERRIDE}"
    )
    return out


def fit_gaussian_from_histogram(
    bin_centers: np.ndarray,
    counts: np.ndarray,
    fit_half_window: float = 0.035,
    fit_range: Tuple[float, float] = (0.09, 0.18),
) -> Optional[Tuple[float, float, float]]:
    if len(bin_centers) != len(counts) or len(bin_centers) < 6:
        return None

    mask_range = (bin_centers >= fit_range[0]) & (bin_centers <= fit_range[1])
    x_range = np.asarray(bin_centers[mask_range], dtype=float)
    y_range = np.asarray(counts[mask_range], dtype=float)
    if len(x_range) < 6:
        return None

    y_range = np.nan_to_num(y_range, nan=0.0, posinf=0.0, neginf=0.0)
    kernel = np.array([1.0, 2.0, 1.0], dtype=float)
    kernel /= kernel.sum()
    y_smooth = np.convolve(y_range, kernel, mode="same")
    if np.all(y_smooth <= 0):
        return None

    baseline = float(np.percentile(y_range, 20))
    signal = np.clip(y_range - baseline, a_min=0.0, a_max=None)
    if np.all(signal <= 0):
        return None

    peak_idx_local = int(np.argmax(y_smooth))
    peak_x = x_range[peak_idx_local]

    mask_window = np.abs(x_range - peak_x) <= fit_half_window
    x_fit = x_range[mask_window]
    y_fit = signal[mask_window]
    if len(x_fit) < 5 or np.sum(y_fit) <= 0:
        return None

    weight_sum = np.sum(y_fit)
    mean = float(np.sum(x_fit * y_fit) / weight_sum)
    variance = float(np.sum(y_fit * (x_fit - mean) ** 2) / weight_sum)
    if variance <= 0:
        return None
    sigma = float(np.sqrt(variance))

    if sigma < 0.0015 or sigma > 0.06:
        y_peak = float(np.max(y_fit))
        half = 0.5 * y_peak
        above = np.where(y_fit >= half)[0]
        if len(above) >= 2:
            fwhm = float(x_fit[above[-1]] - x_fit[above[0]])
            if fwhm > 0:
                sigma = fwhm / 2.354820045

    amplitude = float(np.max(y_fit) + baseline)
    if not np.isfinite(amplitude) or not np.isfinite(mean) or not np.isfinite(sigma) or sigma <= 0:
        return None
    if not (fit_range[0] <= mean <= fit_range[1]):
        return None

    return amplitude, mean, sigma


def combine_branches(lookup: Dict[int, RunConfig],
                     root_dir: Path,
                     target_to_combine: str,
                     efficiency_map: Dict[int, RunEfficiencyMeta],
                     run_filter: Optional[Set[int]] = None) -> pd.DataFrame:
    combined_data: List[pd.DataFrame] = []
    seen_runs: List[int] = []
    total_events = 0

    for run in sorted(lookup.keys()):
        if run_filter and run not in run_filter:
            continue
        if run == 4349:
            print(f"[WARN] Run {run} ignored due to known bad focal-plane data")
            continue

        cfg = lookup[run]
        if target_to_combine and cfg.target and cfg.target.lower() != target_to_combine.lower():
            continue

        fpath = root_dir / f"diagnostics_run{run}.root"
        if not fpath.exists():
            print(f"[WARN] Missing diagnostics ROOT for run {run}: {fpath.name}")
            continue

        eff = efficiency_map.get(run)
        ps_value = eff.ps_value if eff else cfg.prescale_value

        charge_uC = eff.charge_uC if eff and eff.charge_uC else cfg.charge_uC_fallback
        if charge_uC is None or charge_uC <= 0:
            charge_uC = 1000.0
            print(f"[WARN] Missing valid charge for run {run}; using fallback charge_uC=1000")

        livetime = eff.livetime_used if eff else DEBUG_LIVETIME_OVERRIDE
        efficiency = eff.efficiency_used if eff else DEBUG_EFFICIENCY_OVERRIDE

        denom = (charge_uC / 1000.0) * livetime * efficiency
        if denom <= 0:
            denom = 1.0
        scale = float(ps_value) / float(denom)

        print(
            f"[INFO] run={run} ps={ps_value:.6g} charge_uC={charge_uC:.6g} "
            f"livetime={livetime:.6g} efficiency={efficiency:.6g} scale={scale:.6g}"
        )

        try:
            with uproot.open(fpath) as uf:
                if "physics" not in uf:
                    print(f"[WARN] Missing 'physics' tree in {fpath.name}")
                    continue

                physics = uf["physics"]
                branch_names = [b.name for b in physics.branches if b.name not in BRANCHES_TO_EXCLUDE]
                if not branch_names:
                    print(f"[WARN] No branches available for run {run}")
                    continue

                branch_data: Dict[str, np.ndarray] = {}
                reference_len: Optional[int] = None
                for branch_name in branch_names:
                    try:
                        arr = physics[branch_name].array(library="np")
                    except Exception as ex:
                        print(f"[WARN] Could not read branch '{branch_name}' for run {run}: {ex}")
                        continue

                    if reference_len is None:
                        reference_len = len(arr)
                    if len(arr) != reference_len:
                        print(
                            f"[WARN] Skipping branch '{branch_name}' for run {run}: "
                            f"length {len(arr)} != {reference_len}"
                        )
                        continue
                    branch_data[branch_name] = arr

                if not branch_data or reference_len is None or reference_len == 0:
                    print(f"[WARN] No usable branch data for run {run}")
                    continue

                n_events = reference_len
                branch_data["scale"] = np.full(n_events, scale, dtype=np.float32)
                branch_data["run_number"] = np.full(n_events, run, dtype=np.int32)
                branch_data["charge_uC"] = np.full(n_events, charge_uC, dtype=np.float32)
                branch_data["ps_value"] = np.full(n_events, ps_value, dtype=np.float32)
                branch_data["livetime"] = np.full(n_events, livetime, dtype=np.float32)
                branch_data["efficiency"] = np.full(n_events, efficiency, dtype=np.float32)

                df_run = pd.DataFrame(branch_data)
                combined_data.append(df_run)
                seen_runs.append(run)
                total_events += n_events
        except Exception as ex:
            print(f"[ERROR] Failed to process run {run}: {ex}")

    if not combined_data:
        return pd.DataFrame()

    df_combined = pd.concat(combined_data, ignore_index=True)
    print(
        f"[INFO] Combined {len(seen_runs)} runs for target={target_to_combine}; "
        f"events={total_events}; shape={df_combined.shape}"
    )
    return df_combined


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    if len(values) == 0:
        return 0.0
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    total = float(np.sum(w))
    if total <= 0.0:
        return float(v[-1])
    target = max(0.0, min(1.0, quantile)) * total
    running = np.cumsum(w)
    idx = int(np.searchsorted(running, target, side="left"))
    idx = min(idx, len(v) - 1)
    return float(v[idx])


def compute_cov_model(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> Optional[Dict[str, float]]:
    weight_sum = float(np.sum(w))
    if weight_sum <= 0.0 or len(x) < 3:
        return None
    mean_x = float(np.sum(w * x) / weight_sum)
    mean_y = float(np.sum(w * y) / weight_sum)
    dx = x - mean_x
    dy = y - mean_y
    cov_xx = float(np.sum(w * dx * dx) / weight_sum)
    cov_xy = float(np.sum(w * dx * dy) / weight_sum)
    cov_yy = float(np.sum(w * dy * dy) / weight_sum)
    det = cov_xx * cov_yy - cov_xy * cov_xy
    if det <= 0.0 or not np.isfinite(det):
        return None
    return {
        "mean_x": mean_x,
        "mean_y": mean_y,
        "cov_xx": cov_xx,
        "cov_xy": cov_xy,
        "cov_yy": cov_yy,
        "det": det,
        "weight": weight_sum,
        "n_bins": int(len(x)),
    }


def covariance_d2(model: Dict[str, float], x: np.ndarray, y: np.ndarray) -> np.ndarray:
    dx = x - model["mean_x"]
    dy = y - model["mean_y"]
    return (
        model["cov_yy"] * dx * dx
        - 2.0 * model["cov_xy"] * dx * dy
        + model["cov_xx"] * dy * dy
    ) / model["det"]


def ellipse_points(model: Dict[str, float], d2_cut: float, n: int = 240) -> Tuple[np.ndarray, np.ndarray]:
    trace = model["cov_xx"] + model["cov_yy"]
    diff = model["cov_xx"] - model["cov_yy"]
    root = np.sqrt(0.25 * diff * diff + model["cov_xy"] * model["cov_xy"])
    lambda1 = max(0.0, 0.5 * trace + root)
    lambda2 = max(0.0, 0.5 * trace - root)
    angle = 0.5 * np.arctan2(2.0 * model["cov_xy"], diff)
    ca = np.cos(angle)
    sa = np.sin(angle)
    axis1 = np.sqrt(max(0.0, d2_cut * lambda1))
    axis2 = np.sqrt(max(0.0, d2_cut * lambda2))
    t = np.linspace(0.0, 2.0 * np.pi, n + 1)
    u = axis1 * np.cos(t)
    v = axis2 * np.sin(t)
    x = model["mean_x"] + u * ca - v * sa
    y = model["mean_y"] + u * sa + v * ca
    return x, y


def add_combined_2d_mass_cut(df: pd.DataFrame) -> Optional[Dict[str, Any]]:
    """Add combined-level ellipse/MCD exclusivity flags and return ROOT debug payload."""
    required = ["mpi0_all", "mmiss_all", "pi0_weight"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"[WARN] Missing columns for combined 2D mass cut: {missing}")
        df["is_exclusive_ellipse_combined"] = np.zeros(len(df), dtype=np.int32)
        df["is_exclusive_mcd_combined"] = np.zeros(len(df), dtype=np.int32)
        return None

    cfg = MASS_CUT_CONFIG
    x_all = df["mpi0_all"].to_numpy(dtype=float)
    y_all = df["mmiss_all"].to_numpy(dtype=float)
    pi0_w = df["pi0_weight"].to_numpy(dtype=float)
    scale = df["scale"].to_numpy(dtype=float) if "scale" in df.columns else np.ones(len(df), dtype=float)
    event_w = pi0_w * scale

    valid = (
        np.isfinite(x_all) & np.isfinite(y_all) & np.isfinite(event_w) &
        (event_w > 0.0) &
        (x_all >= cfg["mpi0_min"]) & (x_all < cfg["mpi0_max"]) &
        (y_all >= cfg["mmiss_min"]) & (y_all < cfg["mmiss_max"])
    )
    if not np.any(valid):
        print("[WARN] Combined 2D mass-cut histogram is empty after cuts.")
        df["is_exclusive_ellipse_combined"] = np.zeros(len(df), dtype=np.int32)
        df["is_exclusive_mcd_combined"] = np.zeros(len(df), dtype=np.int32)
        return None

    x_edges = np.linspace(cfg["mpi0_min"], cfg["mpi0_max"], cfg["n_mpi0_bins"] + 1)
    y_edges = np.linspace(cfg["mmiss_min"], cfg["mmiss_max"], cfg["n_mmiss_bins"] + 1)
    h2, _, _ = np.histogram2d(x_all[valid], y_all[valid], bins=[x_edges, y_edges], weights=event_w[valid])
    total_weight = float(np.sum(h2))
    if total_weight <= 0.0:
        print("[WARN] Combined 2D mass-cut histogram has zero total weight.")
        df["is_exclusive_ellipse_combined"] = np.zeros(len(df), dtype=np.int32)
        df["is_exclusive_mcd_combined"] = np.zeros(len(df), dtype=np.int32)
        return None

    nx, ny = h2.shape
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    occupied = np.argwhere(h2 > 0.0)
    point_x = x_centers[occupied[:, 0]]
    point_y = y_centers[occupied[:, 1]]
    point_w = h2[occupied[:, 0], occupied[:, 1]]

    seed_mmiss = float(cfg["seed_mmiss"])
    seed_mmiss_half_width = float(cfg["seed_mmiss_half_width"])
    seed_y_bins = np.flatnonzero(np.abs(y_centers - seed_mmiss) <= seed_mmiss_half_width)
    if len(seed_y_bins) == 0:
        print(
            f"[WARN] No mmiss bins lie within {seed_mmiss_half_width:.3f} GeV "
            f"of the requested seed at {seed_mmiss:.3f} GeV."
        )
        df["is_exclusive_ellipse_combined"] = np.zeros(len(df), dtype=np.int32)
        df["is_exclusive_mcd_combined"] = np.zeros(len(df), dtype=np.int32)
        return None

    seed_window = h2[:, seed_y_bins]
    peak_flat = int(np.argmax(seed_window))
    peak_ix, seed_iy_local = np.unravel_index(peak_flat, seed_window.shape)
    peak_iy = int(seed_y_bins[seed_iy_local])
    peak_weight = float(h2[peak_ix, peak_iy])
    if peak_weight <= 0.0:
        print(
            f"[WARN] No occupied 2D mass-cut seed bin lies within "
            f"{seed_mmiss_half_width:.3f} GeV of mmiss={seed_mmiss:.3f} GeV."
        )
        df["is_exclusive_ellipse_combined"] = np.zeros(len(df), dtype=np.int32)
        df["is_exclusive_mcd_combined"] = np.zeros(len(df), dtype=np.int32)
        return None

    def build_core(candidate_peak_fraction: float) -> Tuple[np.ndarray, Dict[str, float]]:
        threshold_weight = candidate_peak_fraction * peak_weight
        mask = np.zeros_like(h2, dtype=bool)
        frontier = deque([(peak_ix, peak_iy)])
        mask[peak_ix, peak_iy] = True

        while frontier:
            ix, iy = frontier.popleft()
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    if dx == 0 and dy == 0:
                        continue
                    jx = ix + dx
                    jy = iy + dy
                    if jx < 0 or jx >= nx or jy < 0 or jy >= ny:
                        continue
                    if mask[jx, jy]:
                        continue
                    if h2[jx, jy] + 1e-12 < threshold_weight:
                        continue
                    mask[jx, jy] = True
                    frontier.append((jx, jy))

        weight = float(np.sum(h2[mask]))
        bins = int(np.count_nonzero(mask))
        stats = {
            "peak_fraction": float(candidate_peak_fraction),
            "weight": weight,
            "total_fraction": weight / total_weight if total_weight > 0.0 else 0.0,
            "bins": bins,
        }
        return mask, stats

    scan_rows = []
    peak_fraction = float(cfg["peak_fraction"])
    auto_jump_ratio = 1.0
    if peak_fraction <= 0.0:
        best_jump = 0.0
        chosen_peak_fraction = float(cfg["auto_peak_max"])
        previous = None
        n_scan = int(np.floor((cfg["auto_peak_max"] - cfg["auto_peak_min"]) / cfg["auto_peak_step"] + 0.5)) + 1
        for i in range(n_scan):
            candidate = max(cfg["auto_peak_min"], cfg["auto_peak_max"] - i * cfg["auto_peak_step"])
            _, current = build_core(candidate)
            weight_ratio = 0.0
            bin_ratio = 0.0
            if previous is not None and previous["bins"] > 0 and previous["weight"] > 0.0:
                weight_ratio = current["weight"] / previous["weight"]
                bin_ratio = current["bins"] / previous["bins"]
                jump = max(weight_ratio, bin_ratio)
                previous_core_is_real = (
                    previous["bins"] >= cfg["auto_min_core_bins"] and
                    previous["total_fraction"] >= cfg["auto_min_core_total_fraction"]
                )
                if previous_core_is_real and jump > best_jump:
                    best_jump = jump
                    chosen_peak_fraction = previous["peak_fraction"]
            scan_rows.append({
                **current,
                "weight_ratio_to_previous": float(weight_ratio),
                "bin_ratio_to_previous": float(bin_ratio),
            })
            previous = current
        peak_fraction = chosen_peak_fraction
        auto_jump_ratio = best_jump

    core_mask, core_stats = build_core(peak_fraction)
    core_idx = np.argwhere(core_mask & (h2 > 0.0))
    if len(core_idx) < 3:
        print("[WARN] Combined 2D mass cut found too few core bins.")
        df["is_exclusive_ellipse_combined"] = np.zeros(len(df), dtype=np.int32)
        df["is_exclusive_mcd_combined"] = np.zeros(len(df), dtype=np.int32)
        return None

    core_model = compute_cov_model(
        x_centers[core_idx[:, 0]],
        y_centers[core_idx[:, 1]],
        h2[core_idx[:, 0], core_idx[:, 1]],
    )
    if core_model is None:
        print("[WARN] Combined 2D mass cut covariance is singular.")
        df["is_exclusive_ellipse_combined"] = np.zeros(len(df), dtype=np.int32)
        df["is_exclusive_mcd_combined"] = np.zeros(len(df), dtype=np.int32)
        return None

    core_d2 = covariance_d2(core_model, x_centers[core_idx[:, 0]], y_centers[core_idx[:, 1]])
    core_w = h2[core_idx[:, 0], core_idx[:, 1]]
    ellipse_d2_cut = weighted_quantile(core_d2, core_w, cfg["ellipse_core_quantile"])
    ellipse_d2_cut *= cfg["ellipse_padding"] * cfg["ellipse_padding"]

    mcd_model = core_model
    best_mcd_model = None
    best_mcd_subset = None
    best_mcd_det = np.inf
    mcd_target_weight = cfg["mcd_keep_total_fraction"] * total_weight

    for _ in range(int(cfg["mcd_iterations"])):
        d2 = covariance_d2(mcd_model, point_x, point_y)
        order = np.argsort(d2)
        running = 0.0
        keep_indices = []
        for idx in order:
            keep_indices.append(idx)
            running += point_w[idx]
            if running >= mcd_target_weight:
                break
        keep_indices = np.asarray(keep_indices, dtype=int)
        candidate_model = compute_cov_model(point_x[keep_indices], point_y[keep_indices], point_w[keep_indices])
        if candidate_model is None:
            break
        mcd_model = candidate_model
        if candidate_model["det"] < best_mcd_det:
            best_mcd_det = candidate_model["det"]
            best_mcd_model = candidate_model
            best_mcd_subset = keep_indices

    if best_mcd_model is None or best_mcd_subset is None:
        best_mcd_model = core_model
        best_mcd_subset = np.arange(len(point_x), dtype=int)

    mcd_subset_d2 = covariance_d2(best_mcd_model, point_x[best_mcd_subset], point_y[best_mcd_subset])
    mcd_subset_w = point_w[best_mcd_subset]
    mcd_d2_cut = weighted_quantile(mcd_subset_d2, mcd_subset_w, cfg["mcd_ellipse_quantile"])
    mcd_d2_cut *= cfg["mcd_padding"] * cfg["mcd_padding"]

    ellipse_flags = np.zeros(len(df), dtype=np.int32)
    mcd_flags = np.zeros(len(df), dtype=np.int32)
    valid_idx = np.where(valid)[0]
    event_ellipse_d2 = covariance_d2(core_model, x_all[valid], y_all[valid])
    event_mcd_d2 = covariance_d2(best_mcd_model, x_all[valid], y_all[valid])
    ellipse_flags[valid_idx[event_ellipse_d2 <= ellipse_d2_cut]] = 1
    mcd_flags[valid_idx[event_mcd_d2 <= mcd_d2_cut]] = 1

    df["is_exclusive_ellipse_combined"] = ellipse_flags
    df["is_exclusive_mcd_combined"] = mcd_flags

    xx, yy = np.meshgrid(x_centers, y_centers, indexing="ij")
    bin_ellipse_mask = covariance_d2(core_model, xx, yy) <= ellipse_d2_cut
    bin_mcd_mask = covariance_d2(best_mcd_model, xx, yy) <= mcd_d2_cut

    ellipse_weight = float(np.sum(event_w[ellipse_flags == 1]))
    mcd_weight = float(np.sum(event_w[mcd_flags == 1]))
    legacy_exclusive_weight = 0.0
    h_legacy_exclusive = np.zeros_like(h2)
    if "is_exclusive" in df.columns:
        legacy_valid = valid & (df["is_exclusive"].to_numpy(dtype=float) > 0.5)
        legacy_exclusive_weight = float(np.sum(event_w[legacy_valid]))
        if np.any(legacy_valid):
            h_legacy_exclusive, _, _ = np.histogram2d(
                x_all[legacy_valid],
                y_all[legacy_valid],
                bins=[x_edges, y_edges],
                weights=event_w[legacy_valid],
            )
    params = {
        "valid": 1.0,
        "peak_fraction": float(peak_fraction),
        "auto_jump_ratio": float(auto_jump_ratio),
        "peak_weight": peak_weight,
        "threshold_weight": float(peak_fraction * peak_weight),
        "total_weight": total_weight,
        "peak_ix": float(peak_ix + 1),
        "peak_iy": float(peak_iy + 1),
        "peak_mpi0": float(x_centers[peak_ix]),
        "peak_mmiss": float(y_centers[peak_iy]),
        "seed_mmiss": seed_mmiss,
        "seed_mmiss_half_width": seed_mmiss_half_width,
        "core_bins": float(core_stats["bins"]),
        "core_weight": float(core_stats["weight"]),
        "core_total_fraction": float(core_stats["total_fraction"]),
        "mean_mpi0": core_model["mean_x"],
        "mean_mmiss": core_model["mean_y"],
        "cov_mpi0_mpi0": core_model["cov_xx"],
        "cov_mpi0_mmiss": core_model["cov_xy"],
        "cov_mmiss_mmiss": core_model["cov_yy"],
        "cov_det": core_model["det"],
        "ellipse_d2_cut": float(ellipse_d2_cut),
        "ellipse_weight": ellipse_weight,
        "ellipse_total_fraction": ellipse_weight / float(np.sum(event_w[valid])),
        "mcd_mean_mpi0": best_mcd_model["mean_x"],
        "mcd_mean_mmiss": best_mcd_model["mean_y"],
        "mcd_cov_mpi0_mpi0": best_mcd_model["cov_xx"],
        "mcd_cov_mpi0_mmiss": best_mcd_model["cov_xy"],
        "mcd_cov_mmiss_mmiss": best_mcd_model["cov_yy"],
        "mcd_det": best_mcd_model["det"],
        "mcd_d2_cut": float(mcd_d2_cut),
        "mcd_subset_weight": float(np.sum(mcd_subset_w)),
        "mcd_subset_bins": float(len(best_mcd_subset)),
        "mcd_weight": mcd_weight,
        "mcd_total_fraction": mcd_weight / float(np.sum(event_w[valid])),
        "legacy_is_exclusive_weight": legacy_exclusive_weight,
        "legacy_is_exclusive_total_fraction": legacy_exclusive_weight / float(np.sum(event_w[valid])),
    }

    print("[INFO] Combined 2D mass cut:")
    print(
        f"       seed mmiss={params['peak_mmiss']:.5f} GeV "
        f"(target {params['seed_mmiss']:.3f} +/- {params['seed_mmiss_half_width']:.3f} GeV)"
    )
    print(f"       peak_fraction={params['peak_fraction']:.3f} core={100.0 * params['core_total_fraction']:.2f}%")
    print(f"       ellipse={int(np.sum(ellipse_flags))} events, MCD={int(np.sum(mcd_flags))} events")

    ellipse_x, ellipse_y = ellipse_points(core_model, ellipse_d2_cut)
    mcd_x, mcd_y = ellipse_points(best_mcd_model, mcd_d2_cut)

    return {
        "histograms": {
            f"{COMBINED_MASS_CUT_TAG}_h_mmiss_vs_mpi0_weighted": (h2, x_edges, y_edges),
            f"{COMBINED_MASS_CUT_TAG}_h_core_selected": (np.where(core_mask, h2, 0.0), x_edges, y_edges),
            f"{COMBINED_MASS_CUT_TAG}_h_legacy_is_exclusive_selected": (h_legacy_exclusive, x_edges, y_edges),
            f"{COMBINED_MASS_CUT_TAG}_h_ellipse_selected": (np.where(bin_ellipse_mask, h2, 0.0), x_edges, y_edges),
            f"{COMBINED_MASS_CUT_TAG}_h_mcd_selected": (np.where(bin_mcd_mask, h2, 0.0), x_edges, y_edges),
            f"{COMBINED_MASS_CUT_TAG}_h_core_mask": (core_mask.astype(float), x_edges, y_edges),
            f"{COMBINED_MASS_CUT_TAG}_h_ellipse_mask": (bin_ellipse_mask.astype(float), x_edges, y_edges),
            f"{COMBINED_MASS_CUT_TAG}_h_mcd_mask": (bin_mcd_mask.astype(float), x_edges, y_edges),
        },
        "params": params,
        "scan": scan_rows,
        "ellipse_line": (ellipse_x, ellipse_y),
        "mcd_line": (mcd_x, mcd_y),
        "peak": (float(x_centers[peak_ix]), float(y_centers[peak_iy])),
    }



def save_to_root(df: pd.DataFrame, output_path: Path, mass_cut_debug: Optional[Dict[str, Any]] = None) -> None:
    print(f"[INFO] Writing combined ROOT tree to {output_path}")
    payload: Dict[str, np.ndarray] = {}
    for col in df.columns:
        payload[col] = df[col].to_numpy()

    with uproot.recreate(str(output_path)) as out_file:
        branch_types = {name: arr.dtype for name, arr in payload.items()}
        tree = out_file.mktree("physics", branch_types)
        tree.extend(payload)
        if mass_cut_debug:
            for name, hist_tuple in mass_cut_debug.get("histograms", {}).items():
                out_file[name] = hist_tuple
            if WRITE_COMBINED_MASS_CUT_PARAM_TREES:
                params = mass_cut_debug.get("params", {})
                if params:
                    out_file[f"{COMBINED_MASS_CUT_TAG}_params"] = {
                        key: np.asarray([value], dtype=np.float64)
                        for key, value in params.items()
                    }
                scan = mass_cut_debug.get("scan", [])
                if scan:
                    out_file[f"{COMBINED_MASS_CUT_TAG}_peak_scan"] = {
                        key: np.asarray([row[key] for row in scan], dtype=np.float64)
                        for key in scan[0].keys()
                    }

    print(f"[INFO] Wrote {len(df)} events with {len(payload)} branches")
    if mass_cut_debug:
        print(f"[INFO] Wrote combined 2D mass-cut debug objects with tag '{COMBINED_MASS_CUT_TAG}'")


def apply_1d_style(ax: plt.Axes,
                   title: str,
                   xlabel: str,
                   weighted: bool = True) -> None:
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("Weighted events / bin" if weighted else "Events / bin", fontsize=10)
    ax.tick_params(axis="both", labelsize=9)
    ax.margins(x=0.01)
    ax.grid(True, alpha=0.3, linewidth=0.6)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def draw_hist1d(ax: plt.Axes,
                data: pd.Series,
                weights: Optional[np.ndarray],
                spec: Hist1DSpec,
                label: Optional[str] = None,
                alpha: float = 0.75,
                weighted: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    values = data.dropna().to_numpy()
    if len(values) == 0:
        ax.text(0.5, 0.5, f"No data for {spec.name}", ha="center", va="center")
        apply_1d_style(ax, spec.title, spec.xlabel, weighted=weighted)
        return np.array([]), np.array([])

    if weights is not None:
        weights = np.asarray(weights)

    counts, edges, _ = ax.hist(
        values,
        bins=spec.bins,
        range=spec.value_range,
        weights=weights,
        alpha=alpha,
        edgecolor="black",
        linewidth=0.6,
        color=spec.color,
        label=label,
    )
    apply_1d_style(ax, spec.title, spec.xlabel, weighted=weighted)
    if spec.value_range is not None:
        ax.set_xlim(spec.value_range)

    # Keep headroom so overlays/annotations are not clipped.
    if len(counts) > 0:
        finite_counts = np.asarray(counts[np.isfinite(counts)], dtype=float)
        if finite_counts.size > 0:
            cmin = float(np.min(finite_counts))
            cmax = float(np.max(finite_counts))
            span = cmax - cmin
            if not (span > 0.0):
                span = abs(cmax) if cmax != 0.0 else 1.0

            pad = 0.12 * span
            target_low = (cmin - pad) if cmin < 0.0 else 0.0
            target_top = cmax + pad

            current_low, current_top = ax.get_ylim()
            ax.set_ylim(min(float(current_low), target_low), max(float(current_top), target_top))

    if label:
        ax.legend(frameon=False, fontsize=9)
    return counts, edges


def write_combined_mass_cut_debug_text(mass_cut_debug: Optional[Dict[str, Any]], output_root: Path):
    if not mass_cut_debug:
        return
    params = mass_cut_debug.get("params", {})
    if not params:
        return

    debug_path = output_root.with_name(f"{output_root.stem}_{COMBINED_MASS_CUT_TAG}_debug.txt")
    keys = [
        "peak_fraction", "auto_jump_ratio",
        "peak_mpi0", "peak_mmiss", "seed_mmiss", "seed_mmiss_half_width",
        "peak_weight", "threshold_weight",
        "core_bins", "core_weight", "core_total_fraction",
        "mean_mpi0", "mean_mmiss", "ellipse_d2_cut",
        "ellipse_weight", "ellipse_total_fraction",
        "mcd_mean_mpi0", "mcd_mean_mmiss", "mcd_d2_cut",
        "mcd_subset_bins", "mcd_subset_weight",
        "mcd_weight", "mcd_total_fraction",
        "legacy_is_exclusive_weight", "legacy_is_exclusive_total_fraction",
    ]

    with debug_path.open("w", encoding="utf-8") as fout:
        fout.write("# Combined 2D mass-cut debug summary\n")
        fout.write(f"tag={COMBINED_MASS_CUT_TAG}\n")
        fout.write(f"weight=pi0_weight*scale\n")
        fout.write(f"ellipse_branch=is_exclusive_ellipse_combined\n")
        fout.write(f"mcd_branch=is_exclusive_mcd_combined\n")
        fout.write("\n[parameters]\n")
        for key in keys:
            if key in params:
                fout.write(f"{key}={params[key]:.12g}\n")
        fout.write("\n[scan_notes]\n")
        fout.write("Full scan not stored in ROOT. To inspect scan, temporarily enable WRITE_COMBINED_MASS_CUT_PARAM_TREES.\n")

    print(f"[INFO] Combined 2D mass-cut debug text saved: {debug_path}")


def write_combined_mass_cut_canvas(mass_cut_debug: Optional[Dict[str, Any]], output_root: Path):
    if not mass_cut_debug:
        return
    histograms = mass_cut_debug.get("histograms", {})
    all_tuple = histograms.get(f"{COMBINED_MASS_CUT_TAG}_h_mmiss_vs_mpi0_weighted")
    core_tuple = histograms.get(f"{COMBINED_MASS_CUT_TAG}_h_core_selected")
    legacy_tuple = histograms.get(f"{COMBINED_MASS_CUT_TAG}_h_legacy_is_exclusive_selected")
    ellipse_tuple = histograms.get(f"{COMBINED_MASS_CUT_TAG}_h_ellipse_selected")
    mcd_tuple = histograms.get(f"{COMBINED_MASS_CUT_TAG}_h_mcd_selected")
    if not all_tuple or not ellipse_tuple or not mcd_tuple:
        return

    h_all, x_edges, y_edges = all_tuple
    h_core, _, _ = core_tuple
    h_legacy = legacy_tuple[0] if legacy_tuple else np.zeros_like(h_all)
    h_ellipse, _, _ = ellipse_tuple
    h_mcd, _, _ = mcd_tuple
    ellipse_x, ellipse_y = mass_cut_debug.get("ellipse_line", (None, None))
    mcd_x, mcd_y = mass_cut_debug.get("mcd_line", (None, None))
    peak_x, peak_y = mass_cut_debug.get("peak", (None, None))
    params = mass_cut_debug.get("params", {})

    canvas_pdf = output_root.with_name(f"{output_root.stem}_{COMBINED_MASS_CUT_TAG}_canvas.pdf")
    canvas_png = output_root.with_name(f"{output_root.stem}_{COMBINED_MASS_CUT_TAG}_canvas.png")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax = axes[0, 0]
    mesh = ax.pcolormesh(x_edges, y_edges, h_all.T, shading="auto", cmap="viridis")
    fig.colorbar(mesh, ax=ax, label="Weighted counts")
    if ellipse_x is not None:
        ax.plot(ellipse_x, ellipse_y, color="magenta", linewidth=2.0, label="ellipse")
    if mcd_x is not None:
        ax.plot(mcd_x, mcd_y, color="lime", linewidth=2.0, label="MCD")
    if np.any(h_legacy > 0.0):
        ax.contour(
            0.5 * (x_edges[:-1] + x_edges[1:]),
            0.5 * (y_edges[:-1] + y_edges[1:]),
            (h_legacy > 0.0).T.astype(float),
            levels=[0.5],
            colors=["blue"],
            linewidths=1.5,
        )
        ax.plot([], [], color="blue", linewidth=1.5, label="de-correlation")
    if peak_x is not None:
        ax.axvline(peak_x, color="black", linestyle="--", linewidth=1.0)
        ax.axhline(peak_y, color="black", linestyle="--", linewidth=1.0)
    ax.set_title("All weighted events with ellipse/MCD overlay")
    ax.set_xlabel("mpi0_all [GeV]")
    ax.set_ylabel("mmiss_all [GeV]")
    ax.legend(loc="upper right")
    ax.text(
        0.02, 0.98,
        f"peak_fraction={params.get('peak_fraction', np.nan):.3f}\n"
        f"core={100.0 * params.get('core_total_fraction', 0.0):.2f}%\n"
        f"de-correlation={100.0 * params.get('legacy_is_exclusive_total_fraction', 0.0):.2f}%\n"
        f"ellipse={100.0 * params.get('ellipse_total_fraction', 0.0):.2f}%\n"
        f"MCD={100.0 * params.get('mcd_total_fraction', 0.0):.2f}%",
        transform=ax.transAxes, va="top", ha="left",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
        fontsize=9,
    )

    ax = axes[0, 1]
    mesh = ax.pcolormesh(x_edges, y_edges, h_mcd.T, shading="auto", cmap="viridis")
    fig.colorbar(mesh, ax=ax, label="Weighted counts")
    if ellipse_x is not None:
        ax.plot(ellipse_x, ellipse_y, color="magenta", linewidth=2.0, label="ellipse")
    if mcd_x is not None:
        ax.plot(mcd_x, mcd_y, color="lime", linewidth=2.0, label="MCD")
    core_idx = np.argwhere(h_core > 0.0)
    if len(core_idx) > 0:
        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
        ax.scatter(x_centers[core_idx[:, 0]], y_centers[core_idx[:, 1]],
                   s=6, color="red", alpha=0.6, label="core bins")
    ax.set_title("MCD-selected weighted events")
    ax.set_xlabel("mpi0_all [GeV]")
    ax.set_ylabel("mmiss_all [GeV]")
    ax.legend(loc="upper right")

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    ax = axes[1, 0]
    ax.step(x_centers, np.sum(h_all, axis=1), where="mid", color="black", label="all")
    ax.step(x_centers, np.sum(h_core, axis=1), where="mid", color="red", label="core")
    ax.step(x_centers, np.sum(h_legacy, axis=1), where="mid", color="blue", label="de-correlation")
    ax.step(x_centers, np.sum(h_ellipse, axis=1), where="mid", color="magenta", label="ellipse")
    ax.step(x_centers, np.sum(h_mcd, axis=1), where="mid", color="green", label="MCD")
    ax.set_title("mpi0_all projection")
    ax.set_xlabel("mpi0_all [GeV]")
    ax.set_ylabel("Weighted counts")
    ax.legend()

    ax = axes[1, 1]
    ax.step(y_centers, np.sum(h_all, axis=0), where="mid", color="black", label="all")
    ax.step(y_centers, np.sum(h_core, axis=0), where="mid", color="red", label="core")
    ax.step(y_centers, np.sum(h_legacy, axis=0), where="mid", color="blue", label="de-correlation")
    ax.step(y_centers, np.sum(h_ellipse, axis=0), where="mid", color="magenta", label="ellipse")
    ax.step(y_centers, np.sum(h_mcd, axis=0), where="mid", color="green", label="MCD")
    ax.set_title("mmiss_all projection")
    ax.set_xlabel("mmiss_all [GeV]")
    ax.set_ylabel("Weighted counts")
    ax.legend()

    fig.suptitle("Combined 2D Mass-Cut Debug", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(canvas_pdf, dpi=150)
    fig.savefig(canvas_png, dpi=150)
    plt.close(fig)
    print(f"[INFO] Combined 2D mass-cut canvas saved: {canvas_pdf}")
    print(f"[INFO] Combined 2D mass-cut canvas saved: {canvas_png}")



def create_focal_plane_debug_plots(df: pd.DataFrame, output_path: Path) -> None:
    print(f"[INFO] Creating focal-plane debug PDF: {output_path}")
    if "run_number" not in df.columns:
        print("[WARN] run_number column missing. Skipping focal-plane debug plots.")
        return

    fp_keywords = ["xfp", "yfp", "xpfp", "ypfp"]
    fp_cols = [c for c in df.columns if any(k in c.lower() for k in fp_keywords)]
    if not fp_cols:
        print("[WARN] No focal-plane columns found.")
        return

    runs = sorted(df["run_number"].unique())
    with PdfPages(str(output_path)) as pdf:
        for run in runs:
            df_run = df[df["run_number"] == run]
            n_vars = len(fp_cols)
            n_cols = 2
            n_rows = int(np.ceil(n_vars / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4 * n_rows))
            axes = np.atleast_1d(axes).flatten()

            fig.suptitle(f"Run {run} focal-plane distributions", fontsize=13, fontweight="bold")

            for idx, col in enumerate(fp_cols):
                spec = Hist1DSpec(
                    name=col,
                    title=f"{col} distribution",
                    xlabel=col,
                    bins=80,
                    value_range=None,
                    color="slateblue",
                )
                draw_hist1d(axes[idx], df_run[col], None, spec, weighted=False)

            for idx in range(n_vars, len(axes)):
                axes[idx].axis("off")

            plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
            pdf.savefig(fig, dpi=130)
            plt.close(fig)

    print("[INFO] Focal-plane debug plots complete")


def create_analysis_plots(df: pd.DataFrame, output_path: Path) -> None:
    print(f"[INFO] Creating analysis PDF: {output_path}")
    if "scale" not in df.columns:
        print("[WARN] Missing 'scale' column; cannot make weighted analysis plots")
        return

    pi0_weight_col = "pi0_weight" if "pi0_weight" in df.columns else (
        "pi0_weights" if "pi0_weights" in df.columns else None
    )
    has_is_exclusive = "is_exclusive" in df.columns

    with PdfPages(str(output_path)) as pdf:
        if "mpi0_all" in df.columns:
            fig, ax = plt.subplots(figsize=(10, 7))
            series = df["mpi0_all"].dropna()
            w_all = df.loc[series.index, "scale"].to_numpy()

            spec_all = Hist1DSpec(
                name="mpi0_all",
                title="pi0 invariant-mass overlay",
                xlabel=r"$m_{\gamma\gamma}$ [GeV/$c^2$]",
                bins=120,
                value_range=(0.0, 0.4),
                color="royalblue",
            )
            draw_hist1d(ax, series, w_all, spec_all, label="All candidates", alpha=0.55)

            if pi0_weight_col:
                w_final = w_all * df.loc[series.index, pi0_weight_col].fillna(0.0).to_numpy()
                spec_fin = Hist1DSpec(
                    name="pi0_weighted",
                    title="pi0 invariant-mass overlay",
                    xlabel=r"$m_{\gamma\gamma}$ [GeV/$c^2$]",
                    bins=120,
                    value_range=(0.0, 0.4),
                    color="firebrick",
                )
                counts_w, edges_w = draw_hist1d(
                    ax, series, w_final, spec_fin, label="Weighted candidates", alpha=0.50
                )
                if len(counts_w) > 0:
                    centers = 0.5 * (edges_w[:-1] + edges_w[1:])
                    fit = fit_gaussian_from_histogram(centers, counts_w)
                    if fit is not None:
                        amp, mu, sigma = fit
                        x_fit = np.linspace(max(0.09, mu - 4 * sigma), min(0.18, mu + 4 * sigma), 300)
                        y_fit = amp * np.exp(-0.5 * ((x_fit - mu) / sigma) ** 2)
                        ax.plot(x_fit, y_fit, color="darkred", linestyle="--", linewidth=2.0,
                                label="Gaussian fit")
                        ax.text(
                            0.98,
                            0.95,
                            f"$\\mu$={mu:.5f} GeV/$c^2$\\n$\\sigma$={sigma:.5f} GeV/$c^2$",
                            transform=ax.transAxes,
                            va="top",
                            ha="right",
                            fontsize=10,
                            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="gray"),
                        )

            ax.set_xlim(0.0, 0.4)
            ax.legend(frameon=False, fontsize=9)
            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close(fig)

        if "mmiss_all" in df.columns:
            fig, ax = plt.subplots(figsize=(10, 7))

            s_mmiss = df["mmiss_all"].dropna()
            w_mmiss = df.loc[s_mmiss.index, "scale"].to_numpy()
            draw_hist1d(
                ax,
                s_mmiss,
                w_mmiss,
                Hist1DSpec(
                    name="mmiss_all",
                    title="Missing-mass overlay",
                    xlabel=r"$M_{miss}$ [GeV/$c^2$]",
                    bins=120,
                    value_range=(0.0, 2.5),
                    color="seagreen",
                ),
                label="mmiss_all",
                alpha=0.55,
            )

            if "mmiss_all_corr" in df.columns:
                s_corr = df["mmiss_all_corr"].dropna()
                w_corr = df.loc[s_corr.index, "scale"].to_numpy()
                draw_hist1d(
                    ax,
                    s_corr,
                    w_corr,
                    Hist1DSpec(
                        name="mmiss_all_corr",
                        title="Missing-mass overlay",
                        xlabel=r"$M_{miss}$ [GeV/$c^2$]",
                        bins=120,
                        value_range=(0.0, 2.5),
                        color="darkorange",
                    ),
                    label="mmiss_all_corr",
                    alpha=0.5,
                )

            ax.set_xlim(0.0, 2.5)
            ax.legend(frameon=False, fontsize=9)
            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close(fig)

        df_plot = df
        if has_is_exclusive:
            df_plot = df[df["is_exclusive"] == True]

        def event_weights(index: pd.Index) -> np.ndarray:
            weights = df_plot.loc[index, "scale"].to_numpy()
            if pi0_weight_col:
                weights = weights * df_plot.loc[index, pi0_weight_col].fillna(1.0).to_numpy()
            return weights

        physics_specs: Dict[str, Hist1DSpec] = {
            "Q2": Hist1DSpec("Q2", "Q2 distribution", r"$Q^2$ [GeV$^2$]", 80, (0, 10), "steelblue"),
            "W": Hist1DSpec("W", "W distribution", r"$W$ [GeV]", 80, (0.5, 4.5), "steelblue"),
            "t": Hist1DSpec("t", "t distribution", r"$t$ [GeV$^2$]", 80, (-5.0, 0.5), "steelblue"),
            "tmin": Hist1DSpec("tmin", "tmin distribution", r"$t_{min}$ [GeV$^2$]", 80, (-5.0, 0.0), "steelblue"),
            "pt": Hist1DSpec("pt", "pt distribution", r"$p_T$ [GeV/$c$]", 80, (0.0, 1.0), "steelblue"),
            "theta": Hist1DSpec("theta", "theta distribution", r"$\theta$ [rad]", 80, (0.0, 0.5), "steelblue"),
            "phi": Hist1DSpec("phi", "phi distribution", r"$\phi$ [rad]", 80, (-3.2, 3.2), "steelblue"),
            "xB": Hist1DSpec("xB", "xB distribution", r"$x_B$", 80, (0.0, 1.0), "steelblue"),
            "z": Hist1DSpec("z", "z distribution", r"$z$", 80, (0.0, 1.2), "steelblue"),
        }

        available_physics = [k for k in physics_specs if k in df_plot.columns]
        if available_physics:
            n_cols = 3
            n_rows = int(np.ceil(len(available_physics) / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows))
            axes = np.atleast_1d(axes).flatten()
            fig.suptitle("Physics 1D distributions", fontsize=13, fontweight="bold")

            for idx, key in enumerate(available_physics):
                s = df_plot[key].dropna()
                w = event_weights(s.index)
                draw_hist1d(axes[idx], s, w, physics_specs[key], alpha=0.75)

            for idx in range(len(available_physics), len(axes)):
                axes[idx].axis("off")
            plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
            pdf.savefig(fig, dpi=150)
            plt.close(fig)

        spec_specs: Dict[str, Hist1DSpec] = {
            "delta": Hist1DSpec("delta", "delta distribution", r"$\delta$ [%]", 80, None, "coral"),
            "xptar": Hist1DSpec("xptar", "xptar distribution", r"$x'_{tar}$ [rad]", 80, None, "coral"),
            "yptar": Hist1DSpec("yptar", "yptar distribution", r"$y'_{tar}$ [rad]", 80, None, "coral"),
            "xtar": Hist1DSpec("xtar", "xtar distribution", r"$x_{tar}$ [cm]", 80, None, "coral"),
            "ytar": Hist1DSpec("ytar", "ytar distribution", r"$y_{tar}$ [cm]", 80, None, "coral"),
            "xfp": Hist1DSpec("xfp", "xfp distribution", r"$x_{fp}$ [cm]", 80, None, "coral"),
            "yfp": Hist1DSpec("yfp", "yfp distribution", r"$y_{fp}$ [cm]", 80, None, "coral"),
            "xpfp": Hist1DSpec("xpfp", "xpfp distribution", r"$x'_{fp}$ [rad]", 80, None, "coral"),
            "ypfp": Hist1DSpec("ypfp", "ypfp distribution", r"$y'_{fp}$ [rad]", 80, None, "coral"),
        }

        available_spec = [k for k in spec_specs if k in df_plot.columns]
        if available_spec:
            n_cols = 3
            n_rows = int(np.ceil(len(available_spec) / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows))
            axes = np.atleast_1d(axes).flatten()
            fig.suptitle("Spectrometer 1D distributions", fontsize=13, fontweight="bold")

            for idx, key in enumerate(available_spec):
                s = df_plot[key].dropna()
                w = event_weights(s.index)
                draw_hist1d(axes[idx], s, w, spec_specs[key], alpha=0.75)

            for idx in range(len(available_spec), len(axes)):
                axes[idx].axis("off")
            plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
            pdf.savefig(fig, dpi=150)
            plt.close(fig)

    print("[INFO] Analysis plots complete")


def print_summary_statistics(df: pd.DataFrame) -> None:
    print("\n" + "=" * 72)
    print("COMBINED DATASET SUMMARY")
    print("=" * 72)
    print(f"events={len(df)}  columns={len(df.columns)}")

    if "run_number" in df.columns:
        runs = sorted(df["run_number"].unique())
        print(f"runs={runs}")
        print(df["run_number"].value_counts().sort_index())

    for col in ["scale", "charge_uC", "ps_value", "livetime", "efficiency"]:
        if col in df.columns:
            print(
                f"{col}: min={df[col].min():.6g} max={df[col].max():.6g} "
                f"mean={df[col].mean():.6g}"
            )


def run(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    cfg = resolve_workflow_config(args)

    print("=" * 72)
    print("COMBINE ANALYSIS BRANCHES")
    print("=" * 72)
    print(f"[INFO] CFG_PATH={cfg.cfg_path}")
    print(f"[INFO] ROOT_DIR={cfg.root_dir}")
    print(f"[INFO] OUT_COMBINED_ROOT={cfg.out_combined_root}")
    print(f"[INFO] KIN_SETTING={cfg.kin_setting}")
    print(f"[INFO] EFFICIENCY_CSV={cfg.efficiency_csv}")
    print(f"[INFO] TARGET={cfg.target_to_combine}")

    df_cfg = load_config(cfg.cfg_path)
    lookup = build_lookup(df_cfg, cfg.kin_setting, cfg.allowed_types)
    eff_map = load_efficiency_metadata(cfg.efficiency_csv, cfg.kin_setting)

    run_filter_set = set(cfg.run_filter) if cfg.run_filter else None
    df_combined = combine_branches(
        lookup=lookup,
        root_dir=cfg.root_dir,
        target_to_combine=cfg.target_to_combine,
        efficiency_map=eff_map,
        run_filter=run_filter_set,
    )

    if df_combined.empty:
        print("[ERROR] No events were combined.")
        return 1

    mass_cut_debug = None
    if CREATE_COMBINED_2D_MASS_CUT:
        mass_cut_debug = add_combined_2d_mass_cut(df_combined)

    print_summary_statistics(df_combined)
    write_combined_mass_cut_debug_text(mass_cut_debug, cfg.out_combined_root)
    write_combined_mass_cut_canvas(mass_cut_debug, cfg.out_combined_root)
    save_to_root(df_combined, cfg.out_combined_root, mass_cut_debug)

    if cfg.create_fp_debug_plots:
        create_focal_plane_debug_plots(df_combined, cfg.fp_debug_pdf)
    if cfg.create_analysis_plots:
        create_analysis_plots(df_combined, cfg.analysis_plots_pdf)

    print("[INFO] Combine stage complete")
    return 0


def main() -> None:
    try:
        status = run()
    except Exception as ex:
        print(f"[ERROR] {ex}", file=sys.stderr)
        status = 2
    raise SystemExit(status)


if __name__ == "__main__":
    main()
