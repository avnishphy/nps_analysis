#!/usr/bin/env python3
"""
plot_publication_quality.py

Publication-quality plots for π⁰ analysis with prescale corrections.

Generates multiple diagnostic plots from summary data with prescale corrections.
Plots include rate calculations, normalized yields, and trend diagnostics.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import re
import sys
from typing import Optional, Tuple
import warnings
warnings.filterwarnings('ignore')

# ===== Command-line Arguments =====
# Parse key=value style arguments for easy readability
# Usage: python script.py run_status=OK target=LH2
run_status_filter = None
target_filter = None

for arg in sys.argv[1:]:
    if '=' in arg:
        key, value = arg.split('=', 1)
        key = key.strip().lower()
        value = value.strip()
        if key == 'run_status':
            run_status_filter = value
        elif key == 'target':
            target_filter = value

if run_status_filter or target_filter:
    filters_str = []
    if run_status_filter:
        filters_str.append(f"run_status={run_status_filter}")
    if target_filter:
        filters_str.append(f"target={target_filter}")
    print(f"[INFO] Filters applied: {', '.join(filters_str)}")

# ===== Configuration =====
SUMMARY_CSV = "output/plots/x60_4b/summary_all_runs.csv"
CONFIG_CSV = "config/nps_dvcs_all_kins_main.csv"
OUTPUT_DIR = "output/plots/x60_4b/publication_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Publication-quality plot settings
FIGURE_DPI = 300
FIGURE_SIZE_SINGLE = (10, 6)
FIGURE_SIZE_COMPARISON = (14, 10)
FONT_SIZE_TITLE = 16
FONT_SIZE_LABEL = 13
FONT_SIZE_LEGEND = 11
FONT_SIZE_TICK = 11
COLOR_MAIN = "#1f77b4"
COLOR_ACCENT = "#ff7f0e"
EDGE_COLOR = "black"
EDGE_WIDTH = 0.5
MARKER_SIZE = 70
ALPHA_SCATTER = 0.75

# ===== Helper Functions =====
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
    """
    Convert token like 'ps4=3' -> prescale_value using rule:
      if r <= 0 or parse fails -> 1
      else prescale_value = 2**(r-1) + 1
    """
    r = extract_prescale_r(token)
    if r is None:
        return 1
    if r <= 0:
        return 1
    return 2 ** (r - 1) + 1


def safe_int(x):
    """Safely convert to int, return None on failure."""
    try:
        return int(x)
    except Exception:
        return None


def setup_publication_plot(figsize=FIGURE_SIZE_SINGLE):
    """Create a figure with publication-quality settings."""
    fig, ax = plt.subplots(figsize=figsize, dpi=FIGURE_DPI)
    ax.grid(True, linestyle='--', alpha=0.3, linewidth=0.7)
    ax.set_axisbelow(True)
    return fig, ax


def format_publication_plot(ax, xlabel, ylabel, title, tight=True):
    """Apply consistent formatting to publication plots."""
    ax.set_xlabel(xlabel, fontsize=FONT_SIZE_LABEL, weight='bold')
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL, weight='bold')
    ax.set_title(title, fontsize=FONT_SIZE_TITLE, weight='bold', pad=15)
    ax.tick_params(labelsize=FONT_SIZE_TICK)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    if tight:
        plt.tight_layout()


def scatter_plot(ax, x, y, xlabel, ylabel, title, filename, pdf=None, ylim=None):
    """Create a scatter plot with error handling."""
    # Remove NaN values
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]
    y_clean = y[mask]
    
    if len(x_clean) == 0:
        print(f"[WARN] No valid data for plot: {filename}")
        return
    
    ax.scatter(x_clean, y_clean, s=MARKER_SIZE, color=COLOR_MAIN, 
               edgecolor=EDGE_COLOR, linewidth=EDGE_WIDTH, alpha=ALPHA_SCATTER, zorder=3)
    format_publication_plot(ax, xlabel, ylabel, title)
    if ylim is not None:
        ax.set_ylim(ylim)
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=FIGURE_DPI, bbox_inches='tight')
    if pdf is not None:
        pdf.savefig(bbox_inches='tight')
    plt.close()
    print(f"[SAVE] {filename}")


def multi_scatter_plot(data_dict, xlabel, ylabel, title, filename, figsize=FIGURE_SIZE_COMPARISON, pdf=None):
    """Create a multi-panel scatter plot for comparison."""
    n_plots = len(data_dict)
    cols = min(2, n_plots)
    rows = (n_plots + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=figsize, dpi=FIGURE_DPI)
    if rows == 1 and cols == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    
    for idx, (label, (x, y)) in enumerate(data_dict.items()):
        ax = axes[idx]
        mask = ~(np.isnan(x) | np.isnan(y))
        x_clean = x[mask]
        y_clean = y[mask]
        
        if len(x_clean) > 0:
            ax.scatter(x_clean, y_clean, s=MARKER_SIZE, color=COLOR_MAIN,
                      edgecolor=EDGE_COLOR, linewidth=EDGE_WIDTH, alpha=ALPHA_SCATTER, zorder=3)
        ax.grid(True, linestyle='--', alpha=0.3, linewidth=0.7)
        ax.set_axisbelow(True)
        ax.set_xlabel(xlabel, fontsize=FONT_SIZE_LABEL, weight='bold')
        ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL - 1, weight='bold')
        ax.set_title(label, fontsize=FONT_SIZE_LABEL, weight='bold')
        ax.tick_params(labelsize=FONT_SIZE_TICK - 1)
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)
    
    # Hide unused subplots
    for idx in range(len(data_dict), len(axes)):
        axes[idx].set_visible(False)
    
    fig.suptitle(title, fontsize=FONT_SIZE_TITLE, weight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=FIGURE_DPI, bbox_inches='tight')
    if pdf is not None:
        pdf.savefig(bbox_inches='tight')
    plt.close()
    print(f"[SAVE] {filename}")


def comparison_line_plot(df, x_col, y_cols_dict, xlabel, ylabel, title, filename, pdf=None, ylim=None):
    """Create a line plot for comparing multiple metrics on the same canvas."""
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE, dpi=FIGURE_DPI)
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    for idx, (label, col) in enumerate(y_cols_dict.items()):
        if col not in df.columns:
            print(f"[WARN] Column '{col}' not found for plot {filename}")
            continue
        
        y = df[col].values
        mask = ~np.isnan(y)
        x_clean = df[x_col].values[mask]
        y_clean = y[mask]
        
        if len(x_clean) > 0:
            ax.plot(x_clean, y_clean, marker='o', linewidth=2, markersize=6,
                   label=label, color=colors[idx % len(colors)], alpha=0.8, zorder=2)
    
    ax.grid(True, linestyle='--', alpha=0.3, linewidth=0.7)
    ax.set_axisbelow(True)
    format_publication_plot(ax, xlabel, ylabel, title)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.legend(fontsize=FONT_SIZE_LEGEND, loc='best', framealpha=0.95)
    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=FIGURE_DPI, bbox_inches='tight')
    if pdf is not None:
        pdf.savefig(bbox_inches='tight')
    plt.close()
    print(f"[SAVE] {filename}")


# ===== Data Loading =====
print("[INFO] Loading CSV files...")
if not os.path.exists(SUMMARY_CSV):
    print(f"[ERROR] Summary CSV not found: {SUMMARY_CSV}", file=sys.stderr)
    sys.exit(1)
if not os.path.exists(CONFIG_CSV):
    print(f"[ERROR] Config CSV not found: {CONFIG_CSV}", file=sys.stderr)
    sys.exit(1)

df = pd.read_csv(SUMMARY_CSV)
df_cfg = pd.read_csv(CONFIG_CSV)

# Normalize column names
df.columns = [c.strip() for c in df.columns]
df_cfg.columns = [c.strip() for c in df_cfg.columns]

print(f"[INFO] Loaded {len(df)} runs from summary CSV")
print(f"[INFO] Loaded {len(df_cfg)} config entries from config CSV")

# Preliminary NPS efficiency for pi0 detection
# calculated from MC simulations as ("nClusters","phot1_hit==1 && phot2_hit==1 && nClusters==2", "")/Total generated
NPS_EFFICIENCY_PI0 = 0.40494

# ===== Prescale Lookup =====
ps_cols = [c for c in df_cfg.columns if 'prescale' in c.lower()]
if not ps_cols:
    raise ValueError(f"Config CSV must contain a column with 'prescale' in its header.")
ps_col = ps_cols[0]
print(f"[INFO] Using prescale column '{ps_col}' from config CSV")

run_cols_candidates = ['run', 'run_number', 'runNumber', 'Run', 'Run_Number']
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
    raise ValueError(f"Could not find run number column in config CSV.")
print(f"[INFO] Using run-number column '{cfg_run_col}' from config CSV")

# Ensure 'run' column in summary
if 'run' not in df.columns:
    for c in df.columns:
        if 'run' in c.lower():
            df = df.rename(columns={c: 'run'})
            print(f"[INFO] Renamed summary column '{c}' -> 'run'")
            break

# Find target column in config
target_col = None
for c in df_cfg.columns:
    if c.lower() == 'target':
        target_col = c
        break
if target_col:
    print(f"[INFO] Using target column '{target_col}' from config CSV")
else:
    print(f"[WARN] Target column not found in config CSV. Target filtering will be disabled.")

# Build prescale and target lookup
cfg_runs = df_cfg[cfg_run_col].apply(safe_int)
df_cfg['_run_int'] = cfg_runs
lookup = {}
for idx, row in df_cfg.iterrows():
    r = row['_run_int']
    if r is None:
        continue
    token = str(row[ps_col]) if pd.notna(row[ps_col]) else ""
    ps_val = prescale_from_token(token)
    tgt = str(row[target_col]).strip() if target_col and pd.notna(row[target_col]) else ""
    lookup[r] = (token, ps_val, tgt)

# Map prescale and target into summary df
prescale_tokens = []
prescale_values = []
targets = []
for idx, row in df.iterrows():
    run = safe_int(row['run'])
    token, psval, tgt = ("", 1, "")
    if run in lookup:
        token, psval, tgt = lookup[run]
    else:
        matched = df_cfg[df_cfg[cfg_run_col].astype(str).str.strip() == str(row['run']).strip()]
        if not matched.empty:
            stoken = str(matched.iloc[0][ps_col]) if pd.notna(matched.iloc[0][ps_col]) else ""
            token = stoken
            psval = prescale_from_token(token)
            stgt = str(matched.iloc[0][target_col]).strip() if target_col and pd.notna(matched.iloc[0][target_col]) else ""
            tgt = stgt
        else:
            token = ""
            psval = 1
            tgt = ""
    prescale_tokens.append(token)
    prescale_values.append(psval)
    targets.append(tgt)

df['prescale_token'] = prescale_tokens
df['prescale_value'] = prescale_values
df['target'] = targets

# ===== Apply Filters =====
# Apply run_status filter if requested
if run_status_filter:
    if 'run_status' in df.columns:
        initial_count = len(df)
        df = df[df['run_status'] == run_status_filter]
        filtered_count = len(df)
        print(f"[INFO] Filtered run_status={run_status_filter}: {initial_count} -> {filtered_count} runs")
    else:
        print(f"[WARN] run_status column not found in summary CSV. Using all runs.")

# Apply target filter if requested
if target_filter:
    initial_count = len(df)
    df = df[df['target'] == target_filter]
    filtered_count = len(df)
    print(f"[INFO] Filtered target={target_filter}: {initial_count} -> {filtered_count} runs")

# ===== Data Conversion to Numeric =====
numeric_cols = [
    'accumulated_charge(mC)', 'current_mean_uA', 'CPUT_LT', 'Beam_Time(s)',
    'total_entries', 'pass_hms', 'pass_hms_nps', 'total_coin_entries',
    'chi2_ndf_comb_bg', 'pi0_mu_MeV', 'pi0_sigma_MeV', 'pi0_signal_counts',
    'mmiss_p_mean_GeV', 'mmiss_p_sigma_GeV', 'hms_track_eff'
]
for col in numeric_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

# ===== Calculate Rates and Metrics =====
print("[INFO] Calculating rates and metrics...")

# Raw rate = total_coin_entries / beam_time
df['raw_rate'] = df['total_coin_entries'] / df['Beam_Time(s)']

# LCF corrected rate = total_coin_entries / (beam_time * cpu_lt)
df['lcf_corrected_rate'] = df['total_coin_entries'] / (df['Beam_Time(s)'] * df['CPUT_LT'])

# Current normalized rate = lcf_corrected_rate / current_mean_uA
df['current_normalized_rate'] = df['lcf_corrected_rate'] / df['current_mean_uA']

# Scaled current normalized rate = current_normalized_rate * (pi0_signal_counts / total_coin_entries)
df['scaled_current_normalized_rate'] = df['current_normalized_rate'] * (
    df['pi0_signal_counts'] / df['total_coin_entries']
)

# Normalized yield = (pi0_signal_counts / (accumulated_charge(mC) * cpu_lt * track_eff)) * prescale_value
# Note: track_eff is hms_track_eff; handle zero/NaN values to avoid inf
df['normalized_yield'] = np.where(
    (df['accumulated_charge(mC)'] > 0) & (df['CPUT_LT'] > 0) & (df['hms_track_eff'] > 0),
    (df['pi0_signal_counts'] / (df['accumulated_charge(mC)'] * df['CPUT_LT'] * df['hms_track_eff'])) * df['prescale_value'],
    np.nan
)

# Additional useful metrics
df['signal_to_coin_ratio'] = df['pi0_signal_counts'] / df['total_coin_entries']
df['detection_efficiency'] = df['pi0_signal_counts'] / df['total_entries']

# Sort by run for consistent plotting
df = df.sort_values("run")

print("[INFO] Calculated metrics:")
print(f"  - raw_rate: min={df['raw_rate'].min():.2e}, max={df['raw_rate'].max():.2e}")
print(f"  - lcf_corrected_rate: min={df['lcf_corrected_rate'].min():.2e}, max={df['lcf_corrected_rate'].max():.2e}")
print(f"  - normalized_yield: min={df['normalized_yield'].min():.2e}, max={df['normalized_yield'].max():.2e}")

# ===== Initialize PDF =====
pdf_path = os.path.join(OUTPUT_DIR, "all_plots.pdf")
pdf = PdfPages(pdf_path)
print(f"\n[INFO] Creating plots and saving to PDF: {pdf_path}")

# ===== Individual Rate Plots =====
print("[INFO] Creating individual rate plots...")

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['raw_rate'].values,
            "Run Number", "Raw Rate [Hz]", "Raw Rate vs Run",
            "01_raw_rate_vs_run.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['lcf_corrected_rate'].values,
            "Run Number", "LCF Corrected Rate [Hz]", "LCF Corrected Rate vs Run",
            "02_lcf_corrected_rate_vs_run.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['current_normalized_rate'].values,
            "Run Number", "Current Normalized Rate [Hz/μA]", "Current Normalized Rate vs Run",
            "03_current_normalized_rate_vs_run.png", pdf=pdf, ylim=(0, 0.3))

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['scaled_current_normalized_rate'].values,
            "Run Number", "Scaled Current Norm. Rate", "Scaled Current Normalized Rate vs Run",
            "04_scaled_current_normalized_rate_vs_run.png", pdf=pdf, ylim=(0, 0.3))

# ===== Rate vs Current Plots =====
print("[INFO] Creating rate vs current plots...")

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['raw_rate'].values,
            "Beam Current [μA]", "Raw Rate [Hz]", "Raw Rate vs Beam Current",
            "05_raw_rate_vs_current.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['lcf_corrected_rate'].values,
            "Beam Current [μA]", "LCF Corrected Rate [Hz]", "LCF Corrected Rate vs Beam Current",
            "06_lcf_corrected_rate_vs_current.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['current_normalized_rate'].values,
            "Beam Current [μA]", "Current Normalized Rate [Hz/μA]", "Current Normalized Rate vs Beam Current",
            "07_current_normalized_rate_vs_current.png", pdf=pdf, ylim=(0, 0.3))

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['scaled_current_normalized_rate'].values,
            "Beam Current [μA]", "Scaled Current Norm. Rate", "Scaled Current Normalized Rate vs Beam Current",
            "08_scaled_current_normalized_rate_vs_current.png", pdf=pdf, ylim=(0, 0.3))

# ===== Rate Comparison on Single Canvas =====
print("[INFO] Creating rate comparison plots...")

comparison_line_plot(df, 'run',
    {
        'Raw Rate': 'raw_rate',
        'LCF Corrected Rate': 'lcf_corrected_rate',
        'Current Normalized Rate': 'current_normalized_rate',
        'Scaled Norm. Rate': 'scaled_current_normalized_rate'
    },
    "Run Number", "Rate [various units]", "Rate Metrics Comparison vs Run",
    "09_rate_comparison_vs_run.png", pdf=pdf)

# Normalized comparison (scaled to 0-1 for visibility)
df_norm_temp = df.copy()
for col in ['raw_rate', 'lcf_corrected_rate', 'current_normalized_rate', 'scaled_current_normalized_rate']:
    if df_norm_temp[col].max() > 0:
        df_norm_temp[col] = (df_norm_temp[col] - df_norm_temp[col].min()) / (df_norm_temp[col].max() - df_norm_temp[col].min())

comparison_line_plot(df_norm_temp, 'run',
    {
        'Raw Rate': 'raw_rate',
        'LCF Corrected Rate': 'lcf_corrected_rate',
        'Current Normalized Rate': 'current_normalized_rate',
        'Scaled Norm. Rate': 'scaled_current_normalized_rate'
    },
    "Run Number", "Normalized Rate (0-1)", "Normalized Rate Metrics Comparison vs Run",
    "10_rate_comparison_normalized_vs_run.png", pdf=pdf, ylim=(0, 2))

# ===== Normalized Yield Plots =====
print("[INFO] Creating normalized yield plots...")

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['normalized_yield'].values,
            "Run Number", "Normalized Yield (Prescale Corrected)", "Normalized Yield vs Run",
            "11_normalized_yield_vs_run.png", pdf=pdf, ylim=(30, 100))

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['normalized_yield'].values,
            "Beam Current [μA]", "Normalized Yield (Prescale Corrected)", "Normalized Yield vs Beam Current",
            "12_normalized_yield_vs_current.png", pdf=pdf, ylim=(30, 100))

# ===== π⁰ Peak and Resolution Trends =====
print("[INFO] Creating π⁰ peak and resolution trend plots...")

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['pi0_mu_MeV'].values,
            "Run Number", "π⁰ Peak Position μ [MeV]", "π⁰ Peak Position vs Run",
            "13_pi0_peak_position_vs_run.png", pdf=pdf, ylim=(125,135))

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['pi0_sigma_MeV'].values,
            "Run Number", "π⁰ Peak Width σ [MeV]", "π⁰ Peak Width vs Run",
            "14_pi0_peak_width_vs_run.png", pdf=pdf, ylim=(3.5,6.5))

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['pi0_mu_MeV'].values,
            "Beam Current [μA]", "π⁰ Peak Position μ [MeV]", "π⁰ Peak Position vs Beam Current",
            "15_pi0_peak_position_vs_current.png", pdf=pdf, ylim=(125,135))

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['pi0_sigma_MeV'].values,
            "Beam Current [μA]", "π⁰ Peak Width σ [MeV]", "π⁰ Peak Width vs Beam Current",
            "16_pi0_peak_width_vs_current.png", pdf=pdf, ylim=(3.5,6.5))

# ===== Fit Quality =====
print("[INFO] Creating fit quality plots...")

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['chi2_ndf_comb_bg'].values,
            "Run Number", "χ² / NDF (Combinatorial BG Fit)", "Fit Quality vs Run",
            "17_chi2_ndf_vs_run.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['chi2_ndf_comb_bg'].values,
            "Beam Current [μA]", "χ² / NDF (Combinatorial BG Fit)", "Fit Quality vs Beam Current",
            "18_chi2_ndf_vs_current.png", pdf=pdf)

# ===== Signal and Detection Efficiency =====
print("[INFO] Creating signal and efficiency trend plots...")

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['pi0_signal_counts'].values,
            "Run Number", "π⁰ Signal Counts", "π⁰ Signal Counts vs Run",
            "19_pi0_signal_counts_vs_run.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['pi0_signal_counts'].values,
            "Beam Current [μA]", "π⁰ Signal Counts", "π⁰ Signal Counts vs Beam Current",
            "19b_pi0_signal_counts_vs_current.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['signal_to_coin_ratio'].values,
            "Run Number", "Signal to Coincidence Ratio", "Signal/Coincidence Ratio vs Run",
            "20_signal_to_coin_ratio_vs_run.png", pdf=pdf, ylim=(0, 1))

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['detection_efficiency'].values,
            "Run Number", "Detection Efficiency (Signal/Total Entries)", "Detection Efficiency vs Run",
            "21_detection_efficiency_vs_run.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['current_mean_uA'].values, df['signal_to_coin_ratio'].values,
            "Beam Current [μA]", "Signal to Coincidence Ratio", "Signal/Coincidence Ratio vs Beam Current",
            "22_signal_to_coin_ratio_vs_current.png", pdf=pdf, ylim=(0, 1))

# ===== Beam Current and Accumulated Charge =====
print("[INFO] Creating beam characteristics plots...")

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['current_mean_uA'].values,
            "Run Number", "Mean Beam Current [μA]", "Beam Current vs Run",
            "23_beam_current_vs_run.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['accumulated_charge(mC)'].values,
            "Run Number", "Accumulated Charge [mC]", "Accumulated Charge vs Run",
            "24_accumulated_charge_vs_run.png", pdf=pdf)

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['Beam_Time(s)'].values,
            "Run Number", "Beam Time [s]", "Beam Time vs Run",
            "25_beam_time_vs_run.png", pdf=pdf)

# ===== LCF Factor =====
print("[INFO] Creating LCF factor plot...")

fig, ax = setup_publication_plot()
scatter_plot(ax, df['run'].values, df['CPUT_LT'].values,
            "Run Number", "CPU Live Time Factor", "CPU Live Time Factor vs Run",
            "26_cpu_lt_factor_vs_run.png", pdf=pdf)

# ===== Multi-panel Comparison Plots =====
print("[INFO] Creating multi-panel comparison plots...")

# Rate comparison multi-panel
rate_data = {
    'Raw Rate': (df['run'].values, df['raw_rate'].values),
    'LCF Corrected Rate': (df['run'].values, df['lcf_corrected_rate'].values),
    'Current Normalized Rate': (df['run'].values, df['current_normalized_rate'].values),
    'Scaled Current Norm. Rate': (df['run'].values, df['scaled_current_normalized_rate'].values),
}
multi_scatter_plot(rate_data, "Run Number", "Rate [various units]",
                  "Rate Metrics Comparison - Individual Panels",
                  "27_rate_metrics_multipanel.png", figsize=(14, 10), pdf=pdf)

# π⁰ property comparison
pi0_data = {
    'Peak Position (μ)': (df['run'].values, df['pi0_mu_MeV'].values),
    'Peak Width (σ)': (df['run'].values, df['pi0_sigma_MeV'].values),
    'Signal Counts': (df['run'].values, df['pi0_signal_counts'].values),
    'χ² / NDF': (df['run'].values, df['chi2_ndf_comb_bg'].values),
}
multi_scatter_plot(pi0_data, "Run Number", "Value",
                  "π⁰ Properties - Individual Panels",
                  "28_pi0_properties_multipanel.png", figsize=(14, 10), pdf=pdf)

# ===== Save Augmented Data =====
print("[INFO] Saving augmented data...")

output_csv = os.path.join(OUTPUT_DIR, "summary_with_rates_and_prescale.csv")
df.to_csv(output_csv, index=False)
print(f"[INFO] Augmented CSV saved to: {output_csv}")

# ===== Summary Statistics =====
print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)

print(f"\nRaw Rate (Hz):")
print(f"  Mean:   {df['raw_rate'].mean():.3e}")
print(f"  Min:    {df['raw_rate'].min():.3e}")
print(f"  Max:    {df['raw_rate'].max():.3e}")
print(f"  Std:    {df['raw_rate'].std():.3e}")

print(f"\nLCF Corrected Rate (Hz):")
print(f"  Mean:   {df['lcf_corrected_rate'].mean():.3e}")
print(f"  Min:    {df['lcf_corrected_rate'].min():.3e}")
print(f"  Max:    {df['lcf_corrected_rate'].max():.3e}")
print(f"  Std:    {df['lcf_corrected_rate'].std():.3e}")

print(f"\nCurrent Normalized Rate (Hz/μA):")
print(f"  Mean:   {df['current_normalized_rate'].mean():.3e}")
print(f"  Min:    {df['current_normalized_rate'].min():.3e}")
print(f"  Max:    {df['current_normalized_rate'].max():.3e}")
print(f"  Std:    {df['current_normalized_rate'].std():.3e}")

print(f"\nNormalized Yield (prescale corrected):")
print(f"  Mean:   {df['normalized_yield'].mean():.3e}")
print(f"  Min:    {df['normalized_yield'].min():.3e}")
print(f"  Max:    {df['normalized_yield'].max():.3e}")
print(f"  Std:    {df['normalized_yield'].std():.3e}")

print(f"\nπ⁰ Peak Position μ (MeV):")
print(f"  Mean:   {df['pi0_mu_MeV'].mean():.3f}")
print(f"  Min:    {df['pi0_mu_MeV'].min():.3f}")
print(f"  Max:    {df['pi0_mu_MeV'].max():.3f}")
print(f"  Std:    {df['pi0_mu_MeV'].std():.3f}")

print(f"\nπ⁰ Peak Width σ (MeV):")
print(f"  Mean:   {df['pi0_sigma_MeV'].mean():.3f}")
print(f"  Min:    {df['pi0_sigma_MeV'].min():.3f}")
print(f"  Max:    {df['pi0_sigma_MeV'].max():.3f}")
print(f"  Std:    {df['pi0_sigma_MeV'].std():.3f}")

# Close PDF
pdf.close()
print(f"\n[SAVE] PDF with all plots saved to: {pdf_path}")

print("\n" + "="*70)
print(f"✅ All plots saved to: {OUTPUT_DIR}")
print(f"   Total plots generated: 28")
print(f"   PDF compiled: {pdf_path}")
print("="*70)
