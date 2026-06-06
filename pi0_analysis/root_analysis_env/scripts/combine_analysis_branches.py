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
from pathlib import Path
from typing import Optional, Tuple, Dict, Any, List

import numpy as np
import pandas as pd
import uproot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ===== Configuration =====
CFG_PATH = Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/config/nps_dvcs_all_kins_main.csv")
# ROOT_DIR = Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0")
# OUT_COMBINED_ROOT = ROOT_DIR / "combined_branches_LH2_wfpi0.root"
ROOT_DIR = Path("/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/")
OUT_COMBINED_ROOT = ROOT_DIR / "combined_branches_LH2.root"
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


def save_to_root(df: pd.DataFrame, output_path: Path):
    """Save DataFrame to ROOT file using uproot."""
    print(f"\n[INFO] Saving combined data to ROOT file: {output_path}")
    
    # Convert DataFrame to dictionary of arrays
    data_dict = {}
    for col in df.columns:
        data_dict[col] = df[col].values
    
    # Write to ROOT file
    with uproot.recreate(str(output_path)) as f:
        f.mktree("physics", data_dict)
    
    print(f"[INFO] Successfully saved {len(df)} events to {output_path}")
    print(f"[INFO] Tree name: 'physics'")
    print(f"[INFO] Branches: {list(data_dict.keys())}")


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
    has_is_exclusive = 'is_exclusive' in df.columns
    
    if not has_pi0_weight:
        print("[WARN] Column 'pi0_weight' not found. Using uniform weights.")
    if not has_is_exclusive:
        print("[WARN] Column 'is_exclusive' not found. Skipping exclusivity cuts.")
    
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
        
        # Filter data if is_exclusive is available
        df_phys = df.copy()
        if has_is_exclusive:
            df_phys = df_phys[df_phys['is_exclusive'] == True]
            print(f"[INFO] Applied exclusivity cut: {len(df_phys)}/{len(df)} events")
        
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
    
    # Print summary statistics
    print_summary_statistics(df_combined)
    
    # Create focal plane debug plots
    if CREATE_FP_DEBUG_PLOTS:
        create_focal_plane_debug_plots(df_combined, FP_DEBUG_PDF)
    
    # Create analysis plots
    if CREATE_ANALYSIS_PLOTS:
        create_analysis_plots(df_combined, ANALYSIS_PLOTS_PDF)
    
    # Save to ROOT file
    save_to_root(df_combined, OUT_COMBINED_ROOT)
    
    # Optionally save to other formats
    # save_to_parquet(df_combined, OUT_COMBINED_ROOT)
    # save_to_csv(df_combined, OUT_COMBINED_ROOT)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!")
    print("="*60)


if __name__ == "__main__":
    main()
