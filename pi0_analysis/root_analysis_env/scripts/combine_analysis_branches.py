#!/usr/bin/env python3
"""
combine_analysis_branches.py

Combine π⁰ run event-level branch data into a single dataset with prescale corrections.

This script:
1. Reads the config CSV
2. Reads per-run ROOT diagnostics files (physics tree)
3. Extracts prescale and target information
4. Concatenates branch data from all runs
5. Adds a 'scale' branch (prescale / computer_live_time) to each event
6. Filters by target type (e.g., LH2)
7. Saves combined branch data to ROOT file

Usage:
    python combine_analysis_branches.py
"""

import os
import re
from collections import deque
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List

import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ===== Configuration =====
CFG_PATH = Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/config/nps_dvcs_all_kins_main.csv")
ROOT_DIR = Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0")
OUT_COMBINED_ROOT = ROOT_DIR / "combined_branches_LH2_wfpi0.root"
# ROOT_DIR = Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/")
# OUT_COMBINED_ROOT = ROOT_DIR / "combined_branches_LH2.root"
TARGET_TO_COMBINE = "LH2"
KIN_SETTING = "KinC_x60_4b"

# Branches to exclude from combining (all other branches will be included)
BRANCHES_TO_EXCLUDE = ["event_id"]

# Debug: Create focal plane diagnostic plots
CREATE_FP_DEBUG_PLOTS = False
# FP_DEBUG_PDF = ROOT_DIR / "focal_plane_debug_wfpi0.pdf"
FP_DEBUG_PDF = ROOT_DIR / "focal_plane_debug.pdf"

# Create analysis plots
CREATE_ANALYSIS_PLOTS = True
# ANALYSIS_PLOTS_PDF = ROOT_DIR / "combined_branches_LH2_wfpi0_plots.pdf"
ANALYSIS_PLOTS_PDF = ROOT_DIR / "combined_branches_LH2_plots.pdf"

# Combined-level 2D mass-cut branches. These are recomputed from the combined
# weighted mmiss_all:mpi0_all distribution, not copied from per-run flags.
CREATE_COMBINED_2D_MASS_CUT = True
COMBINED_MASS_CUT_TAG = "combined_2d_mass_cut"
# Keep these off for now; debug text file carries compact important info.
WRITE_COMBINED_MASS_CUT_PARAM_TREES = False

MASS_CUT_CONFIG = {
    "n_mpi0_bins": 160,
    "mpi0_min": 0.11,
    "mpi0_max": 0.15,
    "n_mmiss_bins": 200,
    "mmiss_min": 0.6,
    "mmiss_max": 1.5,
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


def fit_gaussian_from_histogram(
    bin_centers: np.ndarray,
    counts: np.ndarray,
    fit_half_window: float = 0.035,
    fit_range: Tuple[float, float] = (0.09, 0.18),
) -> Optional[Tuple[float, float, float]]:
    """Estimate Gaussian-like peak (amplitude, mean, sigma) from histogram robustly.

    Uses background subtraction and moment/FWHM estimates, which are more stable
    than log-quadratic fits for weighted histograms.

    Returns None if the estimate cannot be performed robustly.
    """
    if len(bin_centers) != len(counts) or len(bin_centers) < 6:
        return None

    mask_range = (bin_centers >= fit_range[0]) & (bin_centers <= fit_range[1])
    x_range = np.asarray(bin_centers[mask_range], dtype=float)
    y_range = np.asarray(counts[mask_range], dtype=float)

    if len(x_range) < 6:
        return None

    y_range = np.nan_to_num(y_range, nan=0.0, posinf=0.0, neginf=0.0)

    # Smooth lightly for peak finding stability.
    kernel = np.array([1.0, 2.0, 1.0], dtype=float)
    kernel /= kernel.sum()
    y_smooth = np.convolve(y_range, kernel, mode='same')

    if np.all(y_smooth <= 0):
        return None

    # Robust baseline estimate for background-subtracted signal.
    baseline = float(np.percentile(y_range, 20))
    signal = np.clip(y_range - baseline, a_min=0.0, a_max=None)

    if np.all(signal <= 0):
        return None

    peak_idx_local = int(np.argmax(y_smooth))
    peak_x = x_range[peak_idx_local]

    # Peak-centered window limits the influence from wide tails/background.
    mask_window = np.abs(x_range - peak_x) <= fit_half_window
    x_fit = x_range[mask_window]
    y_fit = signal[mask_window]

    if len(x_fit) < 5 or np.sum(y_fit) <= 0:
        return None

    # Moment estimate for Gaussian parameters.
    weight_sum = np.sum(y_fit)
    mean = float(np.sum(x_fit * y_fit) / weight_sum)
    variance = float(np.sum(y_fit * (x_fit - mean) ** 2) / weight_sum)

    if variance <= 0:
        return None

    sigma = float(np.sqrt(variance))

    # FWHM-based fallback if moment sigma is pathological.
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

    # Keep physically sensible pi0 window.
    if not (fit_range[0] <= mean <= fit_range[1]):
        return None

    return float(amplitude), float(mean), float(sigma)

# ===== Helper Functions =====
def safe_int(x: Any) -> Optional[int]:
    try:
        if pd.isna(x):
            return None
        return int(float(str(x).strip()))
    except Exception:
        return None


def extract_prescale_r(token: str) -> Optional[int]:
    """
    Parse prescale token and return integer r.
    Accept tokens like: 'ps4=3', 'PS4=2', ' ps4 = 3 ', etc.
    Returns None on parse failure.
    """
    if not isinstance(token, str):
        return None
    token = token.strip()
    m = re.search(r'ps\d*\s*=\s*(\d+)', token, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))
    if '=' in token:
        try:
            rhs = token.split('=')[-1].strip()
            return int(rhs)
        except Exception:
            return None
    return None


def prescale_from_token(token: str) -> int:
    """Convert token like 'ps4=3' -> prescale_value using rule:
      if r <= 0 or parse fails -> 1
      else prescale_value = 2**(r-1) + 1
    
    This matches plot_publication_quality.py prescale formula.
    """
    r = extract_prescale_r(token)
    if r is None:
        return 1
    if r <= 0:
        return 1
    return 2 ** (r - 1) + 1


def load_config(cfg_path: Path) -> pd.DataFrame:
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config CSV not found: {cfg_path}")
    df_cfg = pd.read_csv(cfg_path, dtype=str, keep_default_na=False)
    df_cfg = df_cfg.replace({"": pd.NA})
    return df_cfg


def build_lookup(df_cfg: pd.DataFrame, kin_setting: str = None) -> Tuple[Dict[int, Tuple[str, int, str, float]], str]:
    """Build lookup table mapping run number to (prescale_token, prescale_value, target, computer_live_time)"""
    # Filter by kinematic setting if provided
    if kin_setting:
        kin_col = None
        for cand in ['Kin_old']:
            if cand in df_cfg.columns:
                kin_col = cand
                break
        if kin_col is None:
            print(f"[WARN] Kin_old column not found in config. Skipping kinematic filtering.")
        else:
            initial_rows = len(df_cfg)
            available_kins = df_cfg[kin_col].astype(str).str.strip().unique()
            print(f"[INFO] Available kinematic settings in config (Kin_old column): {available_kins}")
            df_cfg = df_cfg[df_cfg[kin_col].astype(str).str.strip() == kin_setting]
            filtered_rows = len(df_cfg)
            print(f"[INFO] Filtered by kinematic setting '{kin_setting}': {initial_rows} rows → {filtered_rows} rows")
            if filtered_rows == 0:
                print(f"[WARN] No rows found for kinematic setting '{kin_setting}'!")
    
    # Filter out "junk" runs based on Type column
    type_col = None
    for c in df_cfg.columns:
        if c.lower().strip() == 'type':
            type_col = c
            break
    
    if type_col:
        print(f"[INFO] Using Type column '{type_col}' from config CSV")
        initial_rows = len(df_cfg)
        # Filter out rows where Type == "junk" (case-insensitive)
        df_cfg = df_cfg[df_cfg[type_col].astype(str).str.strip().str.lower() != 'junk']
        filtered_rows = len(df_cfg)
        removed_rows = initial_rows - filtered_rows
        if removed_rows > 0:
            print(f"[INFO] Filtered out {removed_rows} 'junk' runs: {initial_rows} rows → {filtered_rows} rows")
        else:
            print(f"[INFO] No 'junk' runs found to filter")
    else:
        print(f"[WARN] Type column not found in config CSV. Skipping junk run filtering.")
    
    # find prescale column
    ps_cols = [c for c in df_cfg.columns if "prescale" in c.lower()]
    if not ps_cols:
        raise ValueError("Config CSV must contain a column with 'prescale' in its header.")
    ps_col = ps_cols[0]
    print(f"[INFO] Using prescale column '{ps_col}' from config CSV")

    # find run column
    run_cols_candidates = ["run_number"]
    cfg_run_col = None
    for c in run_cols_candidates:
        if c in df_cfg.columns:
            cfg_run_col = c
            break
    if cfg_run_col is None:
        for c in df_cfg.columns:
            if 'run' in c.lower():
                cfg_run_col = c
                break
    if cfg_run_col is None:
        raise ValueError("Could not find run number column in config CSV.")
    print(f"[INFO] Using run-number column '{cfg_run_col}' from config CSV")

    # find target column
    target_col = None
    for c in df_cfg.columns:
        if c.lower().strip() == 'target':
            target_col = c
            break
    if target_col:
        print(f"[INFO] Using target column '{target_col}' from config CSV")
    else:
        print(f"[WARN] Target column not found in config CSV. Target filtering will be disabled.")

    # find computer live time column
    cput_col = None
    for cand in ["Computer_live_time"]:
        if cand in df_cfg.columns:
            cput_col = cand
            break
    if cput_col is None:
        for c in df_cfg.columns:
            if 'computer' in c.lower() and 'live' in c.lower():
                cput_col = c
                break
    if cput_col is None:
        for c in df_cfg.columns:
            if 'live' in c.lower():
                cput_col = c
                break
    if cput_col:
        print(f"[INFO] Using computer-live-time column '{cput_col}' from config CSV")
    else:
        print("[WARN] Could not find Computer_live_time column. Defaulting to 1.0 for live time corrections.")

    # Build lookup
    lookup: Dict[int, Tuple[str, int, str, float]] = {}
    cfg_runs = df_cfg[cfg_run_col].apply(safe_int)
    df_cfg['_run_int'] = cfg_runs
    for idx, row in df_cfg.iterrows():
        r = row['_run_int']
        if r is None or pd.isna(r):
            continue
        r = int(r)
        token = str(row[ps_col]) if pd.notna(row[ps_col]) else ""
        ps_val = prescale_from_token(token)
        tgt = str(row[target_col]).strip() if target_col and pd.notna(row[target_col]) else ""
        cput = 1.0
        if cput_col and pd.notna(row[cput_col]):
            try:
                cput = float(str(row[cput_col]).strip())
            except Exception:
                cput = 1.0
        lookup[r] = (token, ps_val, tgt, cput)
    return lookup, cfg_run_col


def load_charge_map() -> Dict[int, float]:
    """Load run charge map as charge_uC from one of known CSV sources."""
    candidates = [
        Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/luminosity_analysis/livetime_results_parallel.csv"),
        Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/scripts/efficiencies/livetime_results_parallel.csv"),
        Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/scripts/efficiencies/livetime_results_parallel_updatedtrig.csv"),
    ]

    def pick_col(columns: List[str], preferred: List[str], contains: Optional[str] = None) -> Optional[str]:
        lower_map = {c.lower().strip(): c for c in columns}
        for p in preferred:
            if p.lower() in lower_map:
                return lower_map[p.lower()]
        if contains is not None:
            for c in columns:
                if contains in c.lower().strip():
                    return c
        return None

    for csv_path in candidates:
        if not csv_path.exists():
            continue
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"[WARN] Failed reading charge CSV {csv_path}: {e}")
            continue

        if df.empty:
            print(f"[WARN] Charge CSV is empty: {csv_path}")
            continue

        cols = list(df.columns)
        run_col = pick_col(cols, ["run", "run_number"], contains="run")
        charge_uC_col = pick_col(cols, ["charge_uC", "charge_uc"]) 
        charge_mC_col = pick_col(cols, ["accumulated_charge(mC)", "accumulated_charge_mc"], contains="charge")

        if run_col is None:
            print(f"[WARN] No run column found in {csv_path}; columns={cols}")
            continue

        if charge_uC_col is None and charge_mC_col is None:
            print(f"[WARN] No recognized charge column found in {csv_path}; columns={cols}")
            continue

        out: Dict[int, float] = {}
        for _, row in df.iterrows():
            run = safe_int(row.get(run_col))
            if run is None:
                continue

            charge_uC = None
            if charge_uC_col is not None:
                try:
                    v = float(row.get(charge_uC_col))
                    if np.isfinite(v) and v > 0:
                        charge_uC = v
                except Exception:
                    charge_uC = None

            if charge_uC is None and charge_mC_col is not None:
                try:
                    v = float(row.get(charge_mC_col))
                    if np.isfinite(v) and v > 0:
                        charge_uC = v * 1000.0
                except Exception:
                    charge_uC = None

            if charge_uC is not None:
                out[int(run)] = float(charge_uC)

        if out:
            print(f"[INFO] Loaded charge map for {len(out)} runs from {csv_path}")
            return out

        print(f"[WARN] No usable charge rows in {csv_path}")

    print("[WARN] No charge CSV source loaded; combiner will use fallback unit charge per run")
    return {}


def combine_branches(lookup: Dict[int, Tuple[str, int, str, float]], root_dir: Path, target_to_combine: str = 'LH2'):
    """Combine event-level branch data from multiple runs.
    
    Returns:
        DataFrame containing all branch data with added 'scale' and 'run_number' columns
    """
    combined_data = []
    seen_runs = []
    total_events = 0

    all_files = {p.name: p for p in root_dir.glob('*.root')}

    # Load charge info from known CSV sources (run -> charge_uC)
    charge_map = load_charge_map()

    for run, (token, ps_val, tgt, cput) in lookup.items():
        # Skip run 4349 due to bad focal plane data
        if run == 4349:
            print(f"[WARN] Run {run} ignored due to bad focal plane data")
            continue
        
        if target_to_combine and tgt and str(tgt).strip().lower() != target_to_combine.strip().lower():
            continue
        
        fname = f"diagnostics_run{run}.root"
        if fname not in all_files:
            print(f"[WARN] root file not found for run {run}: expected {fname}")
            continue
        
        fpath = all_files[fname]
        seen_runs.append(run)

        cput_val = cput if (cput is not None and cput > 0) else 1.0
        charge_uC = charge_map.get(run, None)
        if charge_uC is not None and charge_uC > 0:
            charge_mC = charge_uC / 1000.0
            charge_uC_for_branch = float(charge_uC)
        else:
            charge_mC = 1.0
            charge_uC_for_branch = float(charge_mC * 1000.0)
            print(f"[WARN] No valid charge_uC for run {run}, using charge_mC=1.0")
        scale = float(ps_val) / (float(cput_val) * float(charge_mC))

        print(f"[INFO] Processing run {run}: token={token} ps={ps_val} cput={cput_val:.4f} charge_mC={charge_mC:.4f} scale={scale:.5f}")
        
        try:
            with uproot.open(str(fpath)) as uf:
                # Access physics tree
                if "physics" not in uf:
                    print(f"[WARN] 'physics' tree not found in {fpath}")
                    continue
                
                physics = uf["physics"]
                
                # Get all available branches and exclude specified ones
                available_branches = [b.name for b in physics.branches]
                branches_to_read = [b for b in available_branches if b not in BRANCHES_TO_EXCLUDE]
                excluded_branches = [b for b in available_branches if b in BRANCHES_TO_EXCLUDE]
                
                if excluded_branches:
                    print(f"[INFO] Excluding branches from run {run}: {excluded_branches}")
                
                print(f"[INFO] Reading {len(branches_to_read)} branches from run {run}")
                
                if not branches_to_read:
                    print(f"[WARN] No valid branches found in run {run}")
                    continue
                
                # Read branches as arrays
                branch_data = {}
                for branch_name in branches_to_read:
                    try:
                        branch_data[branch_name] = physics[branch_name].array(library="np")
                    except Exception as e:
                        print(f"[WARN] Could not read branch '{branch_name}' from run {run}: {e}")
                
                if not branch_data:
                    print(f"[WARN] No branch data extracted from run {run}")
                    continue
                
                # Get number of events (use first branch length)
                n_events = len(next(iter(branch_data.values())))
                total_events += n_events
                
                # Add scale, run_number, and charge_uC for each event
                branch_data['scale'] = np.full(n_events, scale, dtype=np.float32)
                branch_data['run_number'] = np.full(n_events, run, dtype=np.int32)
                # Keep charge_uC finite for downstream weighting logic.
                branch_data['charge_uC'] = np.full(n_events, charge_uC_for_branch, dtype=np.float32)
                
                # Convert to DataFrame for this run
                df_run = pd.DataFrame(branch_data)
                combined_data.append(df_run)
                
                print(f"[INFO] Extracted {n_events} events from run {run}")
                
        except Exception as e:
            print(f"[ERROR] Failed to process run {run}: {e}")
            continue
    
    if not combined_data:
        print("[ERROR] No data was successfully combined!")
        return pd.DataFrame()
    
    # Concatenate all runs
    df_combined = pd.concat(combined_data, ignore_index=True)
    
    print(f"\n[INFO] Combined {len(seen_runs)} runs for target {target_to_combine}")
    print(f"[INFO] Total events: {total_events}")
    print(f"[INFO] Final combined dataset shape: {df_combined.shape}")
    print(f"[INFO] Columns: {list(df_combined.columns)}")
    
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
    charge_fraction = np.ones(len(df), dtype=float)
    total_charge_uC = np.nan
    n_charge_runs = 0
    weight_mode = "pi0_weight*scale"
    if "charge_uC" in df.columns and "run_number" in df.columns:
        charge_uC = pd.to_numeric(df["charge_uC"], errors="coerce").to_numpy(dtype=float)
        run_number = pd.to_numeric(df["run_number"], errors="coerce").to_numpy(dtype=float)
        valid_charge = np.isfinite(charge_uC) & (charge_uC > 0.0) & np.isfinite(run_number)
        if np.any(valid_charge):
            run_charge = pd.DataFrame({
                "run_number": run_number[valid_charge].astype(np.int64),
                "charge_uC": charge_uC[valid_charge],
            }).drop_duplicates("run_number")
            total_charge_uC = float(run_charge["charge_uC"].sum())
            n_charge_runs = int(len(run_charge))
            if np.isfinite(total_charge_uC) and total_charge_uC > 0.0:
                charge_fraction = np.zeros(len(df), dtype=float)
                charge_fraction[valid_charge] = charge_uC[valid_charge] / total_charge_uC
                weight_mode = "pi0_weight*scale*charge_fraction"
            else:
                print("[WARN] Invalid total charge for combined 2D mass cut; using pi0_weight*scale")
        else:
            print("[WARN] No valid charge_uC/run_number rows for combined 2D mass cut; using pi0_weight*scale")
    else:
        print("[WARN] Missing charge_uC/run_number for combined 2D mass cut; using pi0_weight*scale")
    event_w = pi0_w * scale * charge_fraction

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

    peak_flat = int(np.argmax(h2))
    peak_ix, peak_iy = np.unravel_index(peak_flat, h2.shape)
    peak_weight = float(h2[peak_ix, peak_iy])
    peak_tolerance = max(1e-12, abs(peak_weight) * 1e-9)

    def build_core(candidate_peak_fraction: float) -> Tuple[np.ndarray, Dict[str, float]]:
        threshold_weight = candidate_peak_fraction * peak_weight
        mask = np.zeros_like(h2, dtype=bool)
        frontier = deque()

        peak_bins = np.argwhere(np.abs(h2 - peak_weight) <= peak_tolerance)
        for ix, iy in peak_bins:
            if h2[ix, iy] + 1e-12 >= threshold_weight:
                mask[ix, iy] = True
                frontier.append((int(ix), int(iy)))

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
        "total_charge_uC": float(total_charge_uC) if np.isfinite(total_charge_uC) else -1.0,
        "n_charge_runs": float(n_charge_runs),
        "peak_ix": float(peak_ix + 1),
        "peak_iy": float(peak_iy + 1),
        "peak_mpi0": float(x_centers[peak_ix]),
        "peak_mmiss": float(y_centers[peak_iy]),
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
    print(f"       weight={weight_mode}")
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
        "weight_mode": weight_mode,
        "scan": scan_rows,
        "ellipse_line": (ellipse_x, ellipse_y),
        "mcd_line": (mcd_x, mcd_y),
        "peak": (float(x_centers[peak_ix]), float(y_centers[peak_iy])),
    }


def save_to_root(df: pd.DataFrame, output_path: Path, mass_cut_debug: Optional[Dict[str, Any]] = None):
    """Save DataFrame to ROOT file using uproot."""
    print(f"\n[INFO] Saving combined data to ROOT file: {output_path}")
    
    # Convert DataFrame to dictionary of arrays
    data_dict = {}
    for col in df.columns:
        data_dict[col] = df[col].values
    
    # Write to ROOT file
    with uproot.recreate(str(output_path)) as f:
        f.mktree("physics", data_dict)
        if mass_cut_debug:
            for name, hist_tuple in mass_cut_debug.get("histograms", {}).items():
                f[name] = hist_tuple
            # Disabled for now: these trees are verbose and less useful than
            # histograms + compact debug text.
            if WRITE_COMBINED_MASS_CUT_PARAM_TREES:
                params = mass_cut_debug.get("params", {})
                if params:
                    f[f"{COMBINED_MASS_CUT_TAG}_params"] = {
                        key: np.asarray([value], dtype=np.float64)
                        for key, value in params.items()
                    }
                scan = mass_cut_debug.get("scan", [])
                if scan:
                    f[f"{COMBINED_MASS_CUT_TAG}_peak_scan"] = {
                        key: np.asarray([row[key] for row in scan], dtype=np.float64)
                        for key in scan[0].keys()
                    }
    
    print(f"[INFO] Successfully saved {len(df)} events to {output_path}")
    print(f"[INFO] Tree name: 'physics'")
    print(f"[INFO] Branches: {list(data_dict.keys())}")
    if mass_cut_debug:
        print(f"[INFO] Wrote combined 2D mass-cut debug objects with tag '{COMBINED_MASS_CUT_TAG}'")


def write_combined_mass_cut_debug_text(mass_cut_debug: Optional[Dict[str, Any]], output_root: Path):
    if not mass_cut_debug:
        return
    params = mass_cut_debug.get("params", {})
    if not params:
        return

    debug_path = output_root.with_name(f"{output_root.stem}_{COMBINED_MASS_CUT_TAG}_debug.txt")
    keys = [
        "peak_fraction", "auto_jump_ratio",
        "peak_mpi0", "peak_mmiss", "peak_weight", "threshold_weight",
        "total_charge_uC", "n_charge_runs",
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
        fout.write(f"weight={mass_cut_debug.get('weight_mode', 'pi0_weight*scale')}\n")
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


def save_to_parquet(df: pd.DataFrame, output_path: Path):
    """Save DataFrame to Parquet format (optional, more efficient for pandas)."""
    parquet_path = output_path.with_suffix('.parquet')
    print(f"\n[INFO] Saving combined data to Parquet file: {parquet_path}")
    df.to_parquet(str(parquet_path), index=False, compression='snappy')
    print(f"[INFO] Successfully saved {len(df)} events to {parquet_path}")


def save_to_csv(df: pd.DataFrame, output_path: Path):
    """Save DataFrame to CSV format (optional, for inspection)."""
    csv_path = output_path.with_suffix('.csv')
    print(f"\n[INFO] Saving combined data to CSV file: {csv_path}")
    df.to_csv(str(csv_path), index=False)
    print(f"[INFO] Successfully saved {len(df)} events to {csv_path}")


def create_focal_plane_debug_plots(df: pd.DataFrame, output_path: Path):
    """Create diagnostic plots for focal plane variables per run and save to PDF."""
    print("\n" + "="*60)
    print("CREATING FOCAL PLANE DEBUG PLOTS (PER RUN)")
    print("="*60)
    
    # Find focal plane related columns
    fp_keywords = ['xfp', 'yfp', 'xpfp', 'ypfp']
    fp_cols = []
    for col in df.columns:
        col_lower = col.lower()
        if any(kw in col_lower for kw in fp_keywords):
            fp_cols.append(col)
    
    if not fp_cols:
        print("[WARN] No focal plane columns found in dataset")
        return
    
    print(f"[INFO] Found {len(fp_cols)} focal plane columns: {fp_cols}")
    
    if 'run_number' not in df.columns:
        print("[WARN] run_number column not found. Cannot create per-run plots.")
        return
    
    runs = sorted(df['run_number'].unique())
    print(f"[INFO] Creating plots for {len(runs)} runs")
    
    with PdfPages(str(output_path)) as pdf:
        # Create plots for each run
        for run in runs:
            df_run = df[df['run_number'] == run]
            n_events = len(df_run)
            
            print(f"[INFO] Plotting run {run} ({n_events} events)")
            
            # Page 1 for this run: Histograms of all focal plane variables
            n_cols = len(fp_cols)
            n_rows = int(np.ceil(n_cols / 2))
            fig, axes = plt.subplots(n_rows, 2, figsize=(14, 4*n_rows))
            if n_cols == 1:
                axes = np.array([axes])
            axes = axes.flatten()
            
            fig.suptitle(f'Run {run} - Focal Plane Distributions ({n_events} events)', 
                        fontsize=14, fontweight='bold', y=0.995)
            
            for idx, col in enumerate(fp_cols):
                ax = axes[idx]
                data = df_run[col].dropna()
                if len(data) > 0:
                    ax.hist(data, bins=80, alpha=0.7, edgecolor='black', linewidth=0.8)
                    ax.set_xlabel(col, fontsize=10)
                    ax.set_ylabel('Counts', fontsize=10)
                    ax.set_title(f'{col}\nEntries: {len(data)}, Mean: {data.mean():.4f}, Std: {data.std():.4f}',
                               fontsize=9)
                    ax.grid(True, alpha=0.3, linewidth=0.5)
                else:
                    ax.text(0.5, 0.5, f'No data for {col}', ha='center', va='center')
                    ax.set_title(col, fontsize=9)
            
            # Hide unused subplots
            for idx in range(len(fp_cols), len(axes)):
                axes[idx].axis('off')
            
            plt.tight_layout()
            pdf.savefig(fig, dpi=100)
            plt.close()
            
            # Page 2 for this run: 2D correlations between focal plane variables
            if len(fp_cols) >= 2:
                # Common pairs to plot
                common_pairs = []
                for col1 in fp_cols:
                    for col2 in fp_cols:
                        if col1 < col2:  # Avoid duplicates
                            # Prioritize x/y and xp/yp pairs
                            if ('xfp' in col1.lower() and 'yfp' in col2.lower() and \
                                'p' not in col1.lower().replace('xfp', '').replace('yfp', '')):
                                common_pairs.insert(0, (col1, col2))
                            elif ('xpfp' in col1.lower() and 'ypfp' in col2.lower()):
                                common_pairs.insert(0, (col1, col2))
                            else:
                                common_pairs.append((col1, col2))
                
                # Plot up to 4 most interesting pairs per run
                n_pairs = min(4, len(common_pairs))
                if n_pairs > 0:
                    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
                    axes = axes.flatten()
                    
                    fig.suptitle(f'Run {run} - Focal Plane Correlations ({n_events} events)',
                               fontsize=14, fontweight='bold')
                    
                    for idx, (col1, col2) in enumerate(common_pairs[:n_pairs]):
                        ax = axes[idx]
                        df_clean = df_run[[col1, col2]].dropna()
                        
                        if len(df_clean) > 0:
                            h = ax.hist2d(df_clean[col1], df_clean[col2], bins=60,
                                         cmap='viridis', cmin=1)
                            plt.colorbar(h[3], ax=ax, label='Counts')
                            ax.set_xlabel(col1, fontsize=10)
                            ax.set_ylabel(col2, fontsize=10)
                            ax.set_title(f'{col1} vs {col2}\nEntries: {len(df_clean)}',
                                       fontsize=9)
                            ax.grid(True, alpha=0.3, linewidth=0.5)
                        else:
                            ax.text(0.5, 0.5, f'No data for\n{col1} vs {col2}',
                                   ha='center', va='center')
                            ax.set_title(f'{col1} vs {col2}', fontsize=9)
                    
                    # Hide unused subplots
                    for idx in range(n_pairs, 4):
                        axes[idx].axis('off')
                    
                    plt.tight_layout()
                    pdf.savefig(fig, dpi=100)
                    plt.close()
        
        # Summary page: Overview of all runs
        print("[INFO] Creating summary overview page")
        
        for col in fp_cols:
            fig, ax = plt.subplots(figsize=(14, 7))
            
            run_means = []
            run_stds = []
            run_counts = []
            
            for run in runs:
                run_data = df[df['run_number'] == run][col].dropna()
                if len(run_data) > 0:
                    run_means.append(run_data.mean())
                    run_stds.append(run_data.std())
                    run_counts.append(len(run_data))
                else:
                    run_means.append(np.nan)
                    run_stds.append(np.nan)
                    run_counts.append(0)
            
            ax.errorbar(runs, run_means, yerr=run_stds, fmt='o-', capsize=5,
                       markersize=8, linewidth=2, alpha=0.7, label=col)
            ax.set_xlabel('Run Number', fontsize=12)
            ax.set_ylabel(f'{col} (mean ± std)', fontsize=12)
            ax.set_title(f'{col} Summary Across All Runs', fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3)
            
            # Add text with counts
            for i, (run, count) in enumerate(zip(runs, run_counts)):
                if count > 0 and not np.isnan(run_means[i]):
                    ax.text(run, run_means[i], f'n={count}',
                           fontsize=7, rotation=45, va='bottom', ha='left')
            
            plt.tight_layout()
            pdf.savefig(fig, dpi=100)
            plt.close()
    
    print(f"[INFO] Focal plane debug plots saved to {output_path}")


def create_analysis_plots(df: pd.DataFrame, output_path: Path):
    """Create comprehensive analysis plots and save to PDF.
    
    Plots include:
    1. Overlay of mpi0_all (scaled) and pi0_final (weighted)
    2. Overlay of mmiss_all (scaled) and mmiss_all_corr  
    3. 1D histograms of physics variables (Q2, W, t, tmin, pt, theta, phi, xB, z)
    4. 1D histograms of spectrometer variables (delta, xptar, yptar, xtar, ytar, xfp, yfp, xpfp, ypfp)
    5. 2D correlation plots for physics variables
    6. 2D correlation plots for spectrometer variables
    """
    print("\n" + "="*60)
    print("CREATING ANALYSIS PLOTS")
    print("="*60)
    
    # Check required columns
    required_cols = ['scale']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"[WARN] Missing required columns: {missing_cols}. Skipping plots.")
        return
    
    # Check for optional weight and filter columns
    pi0_weight_col = 'pi0_weight' if 'pi0_weight' in df.columns else ('pi0_weights' if 'pi0_weights' in df.columns else None)
    has_pi0_weight = pi0_weight_col is not None
    exclusive_col = (
        'is_exclusive_mcd_combined' if 'is_exclusive_mcd_combined' in df.columns else
        'is_exclusive_mcd' if 'is_exclusive_mcd' in df.columns else
        'is_exclusive' if 'is_exclusive' in df.columns else
        None
    )
    has_is_exclusive = exclusive_col is not None
    
    if not has_pi0_weight:
        print("[WARN] Column 'pi0_weight' not found. Using uniform weights.")
    if not has_is_exclusive:
        print("[WARN] No exclusivity column found. Skipping exclusivity cuts.")
    
    with PdfPages(str(output_path)) as pdf:
        
        # ===== Plot 1: Overlay mpi0_all and pi0_final =====
        print("[INFO] Creating Plot 1: mpi0_all overlay")
        if 'mpi0_all' in df.columns:
            fig, ax = plt.subplots(figsize=(10, 7))
            
            # mpi0_all with scale weights
            data_mpi0_all = df['mpi0_all'].dropna()
            weights_all = df.loc[data_mpi0_all.index, 'scale'].values
            
            ax.hist(data_mpi0_all, bins=100, range=(0, 0.4), weights=weights_all, 
                   alpha=0.6, label='mpi0_all (scaled)', 
                   edgecolor='black', linewidth=0.5, color='blue')
            
            # pi0_final if available (using mpi0_all with pi0_weight and scale)
            if has_pi0_weight:
                weights_final = df.loc[data_mpi0_all.index, 'scale'].values * \
                               df.loc[data_mpi0_all.index, pi0_weight_col].fillna(0).values
                counts_w, edges_w, _ = ax.hist(data_mpi0_all, bins=100, range=(0, 0.4), weights=weights_final,
                                              alpha=0.6, label='pi0_final (weighted)', 
                                              edgecolor='black', linewidth=0.5, color='red')

                # Fit weighted pi0 mass peak (mpi0_all with pi0_weight applied)
                centers_w = 0.5 * (edges_w[:-1] + edges_w[1:])
                fit_result = fit_gaussian_from_histogram(centers_w, counts_w)
                if fit_result is not None:
                    amp, mu, sigma = fit_result
                    x_lo = max(0.09, mu - 4.0 * sigma)
                    x_hi = min(0.18, mu + 4.0 * sigma)
                    x_fit_line = np.linspace(x_lo, x_hi, 300)
                    y_fit_line = amp * np.exp(-0.5 * ((x_fit_line - mu) / sigma) ** 2)
                    ax.plot(x_fit_line, y_fit_line, color='darkred', linewidth=2.0,
                            linestyle='--', label='Gaussian fit (weighted)')

                    fit_text = (
                        f"Weighted fit (mpi0_all × {pi0_weight_col}):\\n"
                        f"$\\mu$ = {mu:.5f} GeV/$c^2$\\n"
                        f"$\\sigma$ = {sigma:.5f} GeV/$c^2$"
                    )
                    ax.text(0.98, 0.95, fit_text, transform=ax.transAxes,
                            fontsize=10, va='top', ha='right',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='gray'))
                    print(f"[INFO] Weighted pi0 fit: mu={mu:.6f} GeV/c^2, sigma={sigma:.6f} GeV/c^2")
                else:
                    print("[WARN] Could not fit weighted pi0 peak for mpi0_all.")
            
            ax.set_xlabel(r'$m_{\pi^0}$ [GeV/$c^2$]', fontsize=12)
            ax.set_ylabel('Weighted Counts', fontsize=12)
            ax.set_title(r'$\pi^0$ Mass Distribution', fontsize=14, fontweight='bold')
            ax.set_xlim(0, 0.4)
            ax.legend(fontsize=11)
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close()
        else:
            print("[WARN] Column 'mpi0_all' not found. Skipping Plot 1.")
        
        # ===== Plot 2: Overlay mmiss_all and mmiss_all_corr =====
        print("[INFO] Creating Plot 2: mmiss overlay")
        if 'mmiss_all' in df.columns:
            fig, ax = plt.subplots(figsize=(10, 7))
            
            # mmiss_all with scale weights
            data_mmiss = df['mmiss_all'].dropna()
            weights_mmiss = df.loc[data_mmiss.index, 'scale'].values
            
            ax.hist(data_mmiss, bins=100, range=(0, 2.5), weights=weights_mmiss,
                   alpha=0.6, label='mmiss_all (scaled)',
                   edgecolor='black', linewidth=0.5, color='green')
            
            # mmiss_all_corr if available
            if 'mmiss_all_corr' in df.columns:
                data_mmiss_corr = df['mmiss_all_corr'].dropna()
                weights_corr = df.loc[data_mmiss_corr.index, 'scale'].values
                ax.hist(data_mmiss_corr, bins=100, range=(0, 2.5), weights=weights_corr,
                       alpha=0.6, label='mmiss_all_corr (scaled)',
                       edgecolor='black', linewidth=0.5, color='orange')
            
            ax.set_xlabel(r'Missing Mass [GeV/$c^2$]', fontsize=12)
            ax.set_ylabel('Weighted Counts', fontsize=12)
            ax.set_title('Missing Mass Distribution', fontsize=14, fontweight='bold')
            ax.set_xlim(0, 2.5)
            ax.legend(fontsize=11)
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close()
        else:
            print("[WARN] Column 'mmiss_all' not found. Skipping Plot 2.")
        
        # ===== Plot 3: Physics variables 1D histograms =====
        print("[INFO] Creating Plot 3: Physics variables 1D histograms")
        physics_vars = ['Q2', 'W', 't', 'tmin', 'pt', 'theta', 'phi', 'xB', 'z']
        physics_labels = {
            'Q2': r'$Q^2$ [GeV$^2$]',
            'W': r'$W$ [GeV]',
            't': r'$-t$ [GeV$^2$]',
            'tmin': r'$-t_{min}$ [GeV$^2$]',
            'pt': r'$p_T$ [GeV/$c$]',
            'theta': r'$\theta$ [rad]',
            'phi': r'$\phi$ [rad]',
            'xB': r'$x_B$',
            'z': r'$z$'
        }
        
        # Filter data if an exclusivity selector is available
        df_phys = df.copy()
        if has_is_exclusive:
            df_phys = df_phys[df_phys[exclusive_col] == True]
            print(f"[INFO] Applied exclusivity cut ({exclusive_col}): {len(df_phys)}/{len(df)} events")
        
        available_physics_vars = [v for v in physics_vars if v in df_phys.columns]
        
        if available_physics_vars:
            n_vars = len(available_physics_vars)
            n_cols = 3
            n_rows = int(np.ceil(n_vars / n_cols))
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
            if n_vars == 1:
                axes = np.array([axes])
            axes = axes.flatten()
            
            fig.suptitle('Physics Variables (Exclusive Events)', fontsize=14, fontweight='bold')
            
            for idx, var in enumerate(available_physics_vars):
                ax = axes[idx]
                data = df_phys[var].dropna()
                
                if len(data) > 0:
                    # Calculate weights: pi0_weight * scale
                    weights = df_phys.loc[data.index, 'scale'].values
                    if has_pi0_weight:
                        weights = weights * df_phys.loc[data.index, pi0_weight_col].fillna(1).values
                    
                    ax.hist(data, bins=60, weights=weights,
                           alpha=0.7, edgecolor='black', linewidth=0.5, color='steelblue')
                    ax.set_xlabel(physics_labels.get(var, var), fontsize=10)
                    ax.set_ylabel('Weighted Counts', fontsize=10)
                    ax.set_title(f'{var}\nEntries: {len(data)}', fontsize=9)
                    ax.grid(True, alpha=0.3)
                else:
                    ax.text(0.5, 0.5, f'No data for {var}', ha='center', va='center')
            
            # Hide unused subplots
            for idx in range(n_vars, len(axes)):
                axes[idx].axis('off')
            
            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close()
        else:
            print("[WARN] No physics variables found. Skipping Plot 3.")
        
        # ===== Plot 4: Spectrometer variables 1D histograms =====
        print("[INFO] Creating Plot 4: Spectrometer variables 1D histograms")
        spec_vars = ['delta', 'xptar', 'yptar', 'xtar', 'ytar', 'xfp', 'yfp', 'xpfp', 'ypfp']
        spec_labels = {
            'delta': r'$\delta$ [%]',
            'xptar': r"$x'_{tar}$ [rad]",
            'yptar': r"$y'_{tar}$ [rad]",
            'xtar': r'$x_{tar}$ [cm]',
            'ytar': r'$y_{tar}$ [cm]',
            'xfp': r'$x_{fp}$ [cm]',
            'yfp': r'$y_{fp}$ [cm]',
            'xpfp': r"$x'_{fp}$ [rad]",
            'ypfp': r"$y'_{fp}$ [rad]"
        }
        
        available_spec_vars = [v for v in spec_vars if v in df_phys.columns]
        
        if available_spec_vars:
            n_vars = len(available_spec_vars)
            n_cols = 3
            n_rows = int(np.ceil(n_vars / n_cols))
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
            if n_vars == 1:
                axes = np.array([axes])
            axes = axes.flatten()
            
            fig.suptitle('Spectrometer Variables (Exclusive Events)', fontsize=14, fontweight='bold')
            
            for idx, var in enumerate(available_spec_vars):
                ax = axes[idx]
                data = df_phys[var].dropna()
                
                if len(data) > 0:
                    # Calculate weights: pi0_weight * scale
                    weights = df_phys.loc[data.index, 'scale'].values
                    if has_pi0_weight:
                        weights = weights * df_phys.loc[data.index, 'pi0_weight'].fillna(1).values
                    
                    ax.hist(data, bins=60, weights=weights,
                           alpha=0.7, edgecolor='black', linewidth=0.5, color='coral')
                    ax.set_xlabel(spec_labels.get(var, var), fontsize=10)
                    ax.set_ylabel('Weighted Counts', fontsize=10)
                    ax.set_title(f'{var}\nEntries: {len(data)}', fontsize=9)
                    ax.grid(True, alpha=0.3)
                else:
                    ax.text(0.5, 0.5, f'No data for {var}', ha='center', va='center')
            
            # Hide unused subplots
            for idx in range(n_vars, len(axes)):
                axes[idx].axis('off')
            
            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close()
        else:
            print("[WARN] No spectrometer variables found. Skipping Plot 4.")
        
        # ===== Plot 5: Physics variables 2D correlations =====
        print("[INFO] Creating Plot 5: Physics variables 2D correlations")
        if len(available_physics_vars) >= 2:
            # Create correlation pairs
            corr_pairs = []
            for i, var1 in enumerate(available_physics_vars):
                for var2 in available_physics_vars[i+1:]:
                    corr_pairs.append((var1, var2))
            
            # Plot in batches of 6 per page
            pairs_per_page = 6
            n_pages = int(np.ceil(len(corr_pairs) / pairs_per_page))
            
            for page in range(n_pages):
                start_idx = page * pairs_per_page
                end_idx = min(start_idx + pairs_per_page, len(corr_pairs))
                page_pairs = corr_pairs[start_idx:end_idx]
                
                n_pairs = len(page_pairs)
                n_cols = 3
                n_rows = int(np.ceil(n_pairs / n_cols))
                
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4.5*n_rows))
                if n_pairs == 1:
                    axes = np.array([axes])
                axes = axes.flatten()
                
                page_title = f'Physics Variables Correlations (Page {page+1}/{n_pages})'
                fig.suptitle(page_title, fontsize=14, fontweight='bold')
                
                for idx, (var1, var2) in enumerate(page_pairs):
                    ax = axes[idx]
                    df_clean = df_phys[[var1, var2]].dropna()
                    
                    if len(df_clean) > 0:
                        # Calculate weights
                        weights = df_phys.loc[df_clean.index, 'scale'].values
                        if has_pi0_weight:
                            weights = weights * df_phys.loc[df_clean.index, 'pi0_weight'].fillna(1).values
                        
                        h = ax.hist2d(df_clean[var1], df_clean[var2], bins=50,
                                     weights=weights, cmap='viridis', cmin=1e-6)
                        plt.colorbar(h[3], ax=ax, label='Weighted Counts')
                        ax.set_xlabel(physics_labels.get(var1, var1), fontsize=9)
                        ax.set_ylabel(physics_labels.get(var2, var2), fontsize=9)
                        ax.set_title(f'{var1} vs {var2}', fontsize=9)
                        ax.grid(True, alpha=0.3)
                    else:
                        ax.text(0.5, 0.5, f'No data', ha='center', va='center')
                        ax.set_title(f'{var1} vs {var2}', fontsize=9)
                
                # Hide unused subplots
                for idx in range(n_pairs, len(axes)):
                    axes[idx].axis('off')
                
                plt.tight_layout()
                pdf.savefig(fig, dpi=150)
                plt.close()
        else:
            print("[WARN] Not enough physics variables for correlations. Skipping Plot 5.")
        
        # ===== Plot 6: Spectrometer variables 2D correlations =====
        print("[INFO] Creating Plot 6: Spectrometer variables 2D correlations")
        if len(available_spec_vars) >= 2:
            # Create correlation pairs
            corr_pairs = []
            for i, var1 in enumerate(available_spec_vars):
                for var2 in available_spec_vars[i+1:]:
                    corr_pairs.append((var1, var2))
            
            # Plot in batches of 6 per page
            pairs_per_page = 6
            n_pages = int(np.ceil(len(corr_pairs) / pairs_per_page))
            
            for page in range(n_pages):
                start_idx = page * pairs_per_page
                end_idx = min(start_idx + pairs_per_page, len(corr_pairs))
                page_pairs = corr_pairs[start_idx:end_idx]
                
                n_pairs = len(page_pairs)
                n_cols = 3
                n_rows = int(np.ceil(n_pairs / n_cols))
                
                fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4.5*n_rows))
                if n_pairs == 1:
                    axes = np.array([axes])
                axes = axes.flatten()
                
                page_title = f'Spectrometer Variables Correlations (Page {page+1}/{n_pages})'
                fig.suptitle(page_title, fontsize=14, fontweight='bold')
                
                for idx, (var1, var2) in enumerate(page_pairs):
                    ax = axes[idx]
                    df_clean = df_phys[[var1, var2]].dropna()
                    
                    if len(df_clean) > 0:
                        # Calculate weights
                        weights = df_phys.loc[df_clean.index, 'scale'].values
                        if has_pi0_weight:
                            weights = weights * df_phys.loc[df_clean.index, 'pi0_weight'].fillna(1).values
                        
                        h = ax.hist2d(df_clean[var1], df_clean[var2], bins=50,
                                     weights=weights, cmap='plasma', cmin=1e-6)
                        plt.colorbar(h[3], ax=ax, label='Weighted Counts')
                        ax.set_xlabel(spec_labels.get(var1, var1), fontsize=9)
                        ax.set_ylabel(spec_labels.get(var2, var2), fontsize=9)
                        ax.set_title(f'{var1} vs {var2}', fontsize=9)
                        ax.grid(True, alpha=0.3)
                    else:
                        ax.text(0.5, 0.5, f'No data', ha='center', va='center')
                        ax.set_title(f'{var1} vs {var2}', fontsize=9)
                
                # Hide unused subplots
                for idx in range(n_pairs, len(axes)):
                    axes[idx].axis('off')
                
                plt.tight_layout()
                pdf.savefig(fig, dpi=150)
                plt.close()
        else:
            print("[WARN] Not enough spectrometer variables for correlations. Skipping Plot 6.")
    
    print(f"[INFO] Analysis plots saved to {output_path}")


def print_summary_statistics(df: pd.DataFrame):
    """Print summary statistics of the combined dataset."""
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    print(f"\nTotal events: {len(df)}")
    print(f"Total columns: {len(df.columns)}")
    
    if 'run_number' in df.columns:
        print(f"\nRuns included: {sorted(df['run_number'].unique())}")
        print(f"Number of runs: {df['run_number'].nunique()}")
        print(f"\nEvents per run:")
        print(df['run_number'].value_counts().sort_index())
    
    if 'scale' in df.columns:
        print(f"\nScale factor range: [{df['scale'].min():.5f}, {df['scale'].max():.5f}]")
        print(f"Unique scale values: {df['scale'].nunique()}")
    
    # Numeric column statistics
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) > 0:
        print(f"\nNumeric column statistics:")
        print(df[numeric_cols].describe())


def main():
    print("="*60)
    print("COMBINE ANALYSIS BRANCHES")
    print("="*60)
    
    # Load and filter config
    cfg = load_config(CFG_PATH)
    lookup, run_col = build_lookup(cfg, KIN_SETTING)

    # Combine branch data from all runs
    df_combined = combine_branches(lookup, ROOT_DIR, TARGET_TO_COMBINE)
    
    if df_combined.empty:
        print("[ERROR] No data to save. Exiting.")
        return

    mass_cut_debug = None
    if CREATE_COMBINED_2D_MASS_CUT:
        mass_cut_debug = add_combined_2d_mass_cut(df_combined)
    
    # Print summary statistics
    print_summary_statistics(df_combined)
    
    # Create focal plane debug plots
    if CREATE_FP_DEBUG_PLOTS:
        create_focal_plane_debug_plots(df_combined, FP_DEBUG_PDF)
    
    # Create analysis plots
    if CREATE_ANALYSIS_PLOTS:
        create_analysis_plots(df_combined, ANALYSIS_PLOTS_PDF)

    write_combined_mass_cut_debug_text(mass_cut_debug, OUT_COMBINED_ROOT)
    write_combined_mass_cut_canvas(mass_cut_debug, OUT_COMBINED_ROOT)
    
    # Save to ROOT file
    save_to_root(df_combined, OUT_COMBINED_ROOT, mass_cut_debug)
    
    # Optionally save to other formats
    # save_to_parquet(df_combined, OUT_COMBINED_ROOT)
    # save_to_csv(df_combined, OUT_COMBINED_ROOT)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!")
    print("="*60)


if __name__ == "__main__":
    main()
