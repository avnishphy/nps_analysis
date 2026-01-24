#!/usr/bin/env python3
"""
plot_run_trends_with_prescale.py
Plots diagnostic trends across runs for π⁰ background fits and applies prescale
corrections read from config/nps_dvcs_all_kins_main.csv.

Normalized yield (prescaled) = (pi0_signal_counts / accumulated_charge(mC)) * prescale_value
where prescale_value = 1 if r <= 0 else 2**(r-1) + 1 for a token "psX=r".
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import sys
from typing import Optional

# ===== User Configuration =====
summary_csv = "output/plots/x60_4b/summary_all_runs.csv"        # path to your summary CSV
config_csv  = "config/nps_dvcs_all_kins_main.csv"              # path to config containing prescale tokens
out_dir = "output/plots/x60_4b/trend_plots"                   # output folder for plots
os.makedirs(out_dir, exist_ok=True)

# ===== Helpers =====
def extract_prescale_r(token: str) -> Optional[int]:
    """
    Parse prescale token and return integer r (right hand side).
    Accept tokens like: 'ps4=3', 'PS4=2', ' ps4 = 3 ', etc.
    Returns None on parse failure.
    """
    if not isinstance(token, str):
        return None
    token = token.strip()
    # try regex ps<digits>=<digits>
    m = re.search(r'ps\d*\s*=\s*(\d+)', token, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))
    # fallback: take last '=' piece if present and numeric
    if '=' in token:
        try:
            rhs = token.split('=')[-1].strip()
            return int(rhs)
        except Exception:
            return None
    return None

def prescale_from_token(token: str) -> int:
    """
    Convert token like 'ps4=3' -> prescale_value using rule:
      if r <= 0 or parse fails -> 1
      else prescale_value = 2**(r-1) + 1
    """
    r = extract_prescale_r(token)
    if r is None:
        # parse failed -> assume no prescale (1)
        return 1
    if r <= 0:
        return 1
    return 2 ** (r - 1) + 1

def make_plot(df, xcol, ycol, xlabel, ylabel, title, filename):
    plt.figure(figsize=(7, 5))
    plt.scatter(df[xcol], df[ycol], color="dodgerblue", s=60, edgecolor="black", alpha=0.85)
    plt.xlabel(xlabel, fontsize=13)
    plt.ylabel(ylabel, fontsize=13)
    plt.title(title, fontsize=14, weight="bold")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, filename), dpi=300)
    plt.close()

# ===== Load files & sanity checks =====
if not os.path.exists(summary_csv):
    print(f"[ERROR] Summary CSV not found: {summary_csv}", file=sys.stderr)
    sys.exit(1)
if not os.path.exists(config_csv):
    print(f"[ERROR] Config CSV not found: {config_csv}", file=sys.stderr)
    sys.exit(1)

df = pd.read_csv(summary_csv)
df_cfg = pd.read_csv(config_csv)

# Normalize column names for convenience (strip spaces)
df.columns = [c.strip() for c in df.columns]
df_cfg.columns = [c.strip() for c in df_cfg.columns]

# determine prescale column name in config CSV: search for any column containing 'prescale' (case-insensitive)
ps_cols = [c for c in df_cfg.columns if 'prescale' in c.lower()]
if not ps_cols:
    raise ValueError(f"Config CSV {config_csv} must contain a column with 'prescale' in its header.")
ps_col = ps_cols[0]
print(f"[INFO] Using prescale column '{ps_col}' from {config_csv}")

# determine run-number column name in config CSV: search for likely names
run_cols_candidates = ['run', 'run_number', 'runNumber', 'Run', 'Run_Number']
cfg_run_col = None
for c in run_cols_candidates:
    if c in df_cfg.columns:
        cfg_run_col = c
        break
if cfg_run_col is None:
    # fallback: try to find a numeric-looking column name with 'run' in it
    for c in df_cfg.columns:
        if 'run' in c.lower():
            cfg_run_col = c
            break
if cfg_run_col is None:
    raise ValueError(f"Could not find run number column in {config_csv} (expected names like run_number).")
print(f"[INFO] Using run-number column '{cfg_run_col}' from {config_csv}")

# ensure summary df has 'run' column
if 'run' not in df.columns:
    # try to find something
    candidate = None
    for c in df.columns:
        if 'run' in c.lower():
            candidate = c
            break
    if candidate:
        df = df.rename(columns={candidate: 'run'})
        print(f"[INFO] Renamed summary column '{candidate}' -> 'run'")
    else:
        raise ValueError("summary CSV must have a 'run' column")

# ensure necessary columns exist in summary
if "pi0_signal_counts" not in df.columns:
    raise ValueError("summary CSV must contain 'pi0_signal_counts' column")
if "accumulated_charge(mC)" not in df.columns:
    raise ValueError("summary CSV must contain 'accumulated_charge(mC)' column")

# ===== Build prescale lookup (map run -> prescale_value) =====
# Make sure run types are comparable: cast to int where possible
def safe_int(x):
    try:
        return int(x)
    except Exception:
        return None

cfg_runs = df_cfg[cfg_run_col].apply(safe_int)
df_cfg['_run_int'] = cfg_runs
lookup = {}
missing_tokens = 0
for idx, row in df_cfg.iterrows():
    r = row['_run_int']
    if r is None:
        continue
    token = str(row[ps_col]) if pd.notna(row[ps_col]) else ""
    ps_val = prescale_from_token(token)
    lookup[r] = (token, ps_val)

# map prescale into summary df
prescale_tokens = []
prescale_values = []
not_found = []
for idx, row in df.iterrows():
    run = safe_int(row['run'])
    token, psval = ("", 1)
    if run in lookup:
        token, psval = lookup[run]
    else:
        # try to match as string in cfg 'run_number' if ints fail
        matched = df_cfg[df_cfg[cfg_run_col].astype(str).str.strip() == str(row['run']).strip()]
        if not matched.empty:
            stoken = str(matched.iloc[0][ps_col]) if pd.notna(matched.iloc[0][ps_col]) else ""
            token = stoken
            psval = prescale_from_token(token)
        else:
            not_found.append(run)
            token = ""
            psval = 1
    prescale_tokens.append(token)
    prescale_values.append(psval)
    print(f"[PRESCALE] run={row['run']} token='{token}' -> prescale_value={psval}")

if not_found:
    print(f"[WARN] No prescale token found in config CSV for runs: {not_found}. Defaulting prescale=1 for those.")

df['prescale_token'] = prescale_tokens
df['prescale_value'] = prescale_values

# ===== Compute (un)prescaled normalized yields =====
# safe division: avoid divide-by-zero
df['accumulated_charge(mC)'] = pd.to_numeric(df['accumulated_charge(mC)'], errors='coerce')
df['pi0_signal_counts'] = pd.to_numeric(df['pi0_signal_counts'], errors='coerce')

df['normalized_yield'] = df['pi0_signal_counts'] / df['accumulated_charge(mC)']
df['normalized_yield_prescaled'] = df['normalized_yield'] * df['prescale_value']

# Sorting for consistent plotting
df = df.sort_values("run")

# ===== Plotting =====
def scatter_plot(xcol, ycol, xlabel, ylabel, title, fname):
    make_plot(df, xcol, ycol, xlabel, ylabel, title, fname)
    print(f"[SAVE] {fname}")

# basic trends
scatter_plot("run", "run_current_mode_uA", "Run", "Beam Current [μA]", "Run vs Current", "run_vs_current.png")
scatter_plot("run", "pi0_mu_MeV", "Run", "π⁰ Peak μ [MeV]", "Run vs π⁰ Peak Position", "run_vs_pi0_peak.png")
scatter_plot("run_current_mode_uA", "pi0_mu_MeV", "Beam Current [μA]", "π⁰ Peak μ [MeV]", "Current vs π⁰ Peak Position", "current_vs_pi0_peak.png")
scatter_plot("run", "pi0_sigma_MeV", "Run", "π⁰ Peak σ [MeV]", "Run vs π⁰ Peak Width", "run_vs_pi0_sigma.png")
scatter_plot("run_current_mode_uA", "pi0_sigma_MeV", "Beam Current [μA]", "π⁰ Peak σ [MeV]", "Current vs π⁰ Peak Width", "current_vs_pi0_sigma.png")
scatter_plot("run", "chi2_ndf_comb_bg", "Run", "χ² / NDF (Combinatorial BG Fit)", "Run vs χ²/NDF", "run_vs_chi2_ndf_comb_bg.png")
scatter_plot("run_current_mode_uA", "chi2_ndf_comb_bg", "Beam Current [μA]", "χ² / NDF (Combinatorial BG Fit)", "Current vs χ²/NDF", "current_vs_chi2_ndf_comb_bg.png")

# normalized yield plots (both raw and prescaled)
scatter_plot("run", "normalized_yield", "Run", "Normalized Yield (signal/charge)", "Run vs Normalized Yield", "run_vs_normalized_yield.png")
scatter_plot("run_current_mode_uA", "normalized_yield", "Beam Current [μA]", "Normalized Yield (signal/charge)", "Current vs Normalized Yield", "current_vs_normalized_yield.png")

scatter_plot("run", "normalized_yield_prescaled", "Run", "Prescaled Normalized Yield", "Run vs Prescaled Normalized Yield", "run_vs_normalized_yield_prescaled.png")
scatter_plot("run_current_mode_uA", "normalized_yield_prescaled", "Beam Current [μA]", "Prescaled Normalized Yield", "Current vs Prescaled Normalized Yield", "current_vs_normalized_yield_prescaled.png")

# optional: plot prescale_value itself versus run (diagnostic)
scatter_plot("run", "prescale_value", "Run", "prescale_value", "Run vs prescale value", "run_vs_prescale_value.png")

# Save augmented CSV for inspection
augmented_csv = os.path.join(out_dir, "summary_all_runs_with_prescale.csv")
df.to_csv(augmented_csv, index=False)
print(f"[INFO] Augmented CSV with prescale saved to: {augmented_csv}")
print(f"✅ Plots saved in '{out_dir}/'")
