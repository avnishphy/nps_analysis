#!/usr/bin/env python3
"""
plot_run_trends.py
Plots diagnostic trends across runs for π⁰ background fits.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os

# ===== User Configuration =====
csv_file = "output/plots/x60_4b/summary_all_runs.csv"        # path to your CSV
out_dir = "output/plots/x60_4b/trend_plots"        # output folder for plots
os.makedirs(out_dir, exist_ok=True)

# ===== Load Data =====
df = pd.read_csv(csv_file)

# Compute normalized yield (signal counts / charge)
df["normalized_yield"] = df["pi0_signal_counts"] / df["accumulated_charge(mC)"]

# Sorting for consistent plotting
df = df.sort_values("run")

# ===== Helper plotting function =====
def make_plot(x, y, xlabel, ylabel, title, filename):
    plt.figure(figsize=(7, 5))
    plt.scatter(df[x], df[y], color="dodgerblue", s=60, edgecolor="black", alpha=0.8)
    # plt.plot(df[x], df[y], color="royalblue", linewidth=1.2, alpha=0.6)  # removed line
    plt.xlabel(xlabel, fontsize=13)
    plt.ylabel(ylabel, fontsize=13)
    plt.title(title, fontsize=14, weight="bold")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, filename), dpi=300)
    plt.close()

# ====== Plots ======
make_plot("run", "run_current_mode_uA",
           "Run", "Beam Current [μA]",
           "Run vs Current",
           "run_vs_current.png")

make_plot("run", "pi0_mu_MeV",
           "Run", "π⁰ Peak μ [MeV]",
           "Run vs π⁰ Peak Position",
           "run_vs_pi0_peak.png")

make_plot("run_current_mode_uA", "pi0_mu_MeV",
           "Beam Current [μA]", "π⁰ Peak μ [MeV]",
           "Current vs π⁰ Peak Position",
           "current_vs_pi0_peak.png")

make_plot("run", "pi0_sigma_MeV",
           "Run", "π⁰ Peak σ [MeV]",
           "Run vs π⁰ Peak Width",
           "run_vs_pi0_sigma.png")

make_plot("run_current_mode_uA", "pi0_sigma_MeV",
           "Beam Current [μA]", "π⁰ Peak σ [MeV]",
           "Current vs π⁰ Peak Width",
           "current_vs_pi0_sigma.png")

make_plot("run", "chi2_ndf_comb_bg",
           "Run", "χ² / NDF (Combinatorial BG Fit)",
           "Run vs χ²/NDF",
           "run_vs_chi2_ndf_comb_bg.png")

make_plot("run_current_mode_uA", "chi2_ndf_comb_bg",
           "Beam Current [μA]", "χ² / NDF (Combinatorial BG Fit)",
           "Current vs χ²/NDF",
           "current_vs_chi2_ndf_comb_bg.png")

make_plot("run", "normalized_yield",
           "Run", "Normalized Yield (Signal Counts / Charge)",
           "Run vs Normalized Yield",
           "run_vs_normalized_yield.png")

make_plot("run_current_mode_uA", "normalized_yield",
           "Beam Current [μA]", "Normalized Yield (Signal Counts / Charge)",
           "Current vs Normalized Yield",
           "current_vs_normalized_yield.png")

print(f"✅ Plots saved in '{out_dir}/'")
