#!/usr/bin/env python3
"""
Data vs Simulation Comparison Script
====================================

Compare experimental data with GEANT4/SIMC simulation for π⁰ DVCS analysis.
Applies HMS electron cuts to both datasets and plots key kinematic variables.

Data: /group/nps/singhav/nps_replay/rootfiles_avnish/root_analysis_env_skim/x60_4b/skim_run4398.root
Simulation: /w/hallc-scshelf2102/nps/singhav/geant4_simc/HallC_NPS/DVCS_evt_gen/DVCS/nps_excl_pi0_x60_4b_geant_500k.root

HMS Electron Cuts:
  - h_react_z: -4.0 to 4.0 cm
  - h_delta: -15.0 to 15.0 %
  - h_gtr_th: -0.1 to 0.1 rad
  - h_gtr_ph: -0.04 to 0.04 rad
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot
from pathlib import Path

# Configuration
DATA_FILE = "/group/nps/singhav/nps_replay/rootfiles_avnish/root_analysis_env_skim/x60_4b/skim_run4398.root"
SIM_FILE = "/w/hallc-scshelf2102/nps/singhav/geant4_simc/HallC_NPS/DVCS_evt_gen/DVCS/nps_excl_pi0_x60_4b_geant_500k.root"

# HMS Electron cut parameters
HMS_CUTS = {
    "h_react_z": (-4.0, 4.0),
    "h_delta": (-15.0, 15.0),
    "h_gtr_th": (-0.1, 0.1),
    "h_gtr_ph": (-0.04, 0.04),
}

# Variables to plot (HMS kinematics)
PLOT_VARS = [
    "h_delta",      # Delta (%)
    "h_gtr_th",     # Theta angle (rad)
    "h_gtr_ph",     # Phi angle (rad)
    "h_react_z",    # Reaction z (cm)
    "h_xptar",      # X target angle (rad)
    "h_yptar",      # Y target angle (rad)
    "h_ytar",       # Y target (cm)
    "h_xfp",        # X focal plane (cm)
    "h_yfp",        # Y focal plane (cm)
    "h_xpfp",       # X focal plane angle (rad)
    "h_ypfp",       # Y focal plane angle (rad)
]

# Labels for plotting
PLOT_LABELS = {
    "h_delta": "Delta (%)",
    "h_gtr_th": "Theta (rad)",
    "h_gtr_ph": "Phi (rad)",
    "h_react_z": "Reaction Z (cm)",
    "h_xptar": "X Target Angle (rad)",
    "h_yptar": "Y Target Angle (rad)",
    "h_ytar": "Y Target (cm)",
    "h_xfp": "X Focal Plane (cm)",
    "h_yfp": "Y Focal Plane (cm)",
    "h_xpfp": "X FP Angle (rad)",
    "h_ypfp": "Y FP Angle (rad)",
}


def load_data(filepath, tree_name="ntp"):
    """Load data from ROOT file."""
    print(f"\n{'='*80}")
    print(f"Loading: {Path(filepath).name}")
    print(f"{'='*80}")
    
    if not Path(filepath).exists():
        print(f"✗ File not found: {filepath}")
        return None
    
    try:
        root_file = uproot.open(filepath)
        print(f"✓ File opened successfully")
        print(f"  Available keys: {list(root_file.keys())}")
        
        if tree_name not in root_file:
            # Try to find any tree
            trees = [k for k in root_file.keys() if isinstance(root_file[k], uproot.TTree)]
            if trees:
                tree_name = trees[0]
                print(f"  Tree '{tree_name}' not found, using '{tree_name}' instead")
            else:
                print(f"✗ No tree found in file")
                return None
        
        tree = root_file[tree_name]
        print(f"✓ Tree '{tree_name}' loaded")
        print(f"  Total entries: {tree.num_entries:,}")
        print(f"  Branches: {len(tree.keys())}")
        
        return tree
    except Exception as e:
        print(f"✗ Error loading file: {e}")
        return None


def apply_hms_cuts(df, cuts_dict):
    """Apply HMS electron cuts to DataFrame."""
    mask = np.ones(len(df), dtype=bool)
    
    for var, (min_val, max_val) in cuts_dict.items():
        if var in df.columns:
            var_mask = (df[var] >= min_val) & (df[var] <= max_val)
            n_before = mask.sum()
            mask = mask & var_mask
            n_after = mask.sum()
            n_removed = n_before - n_after
            pct = 100 * n_removed / n_before if n_before > 0 else 0
            print(f"  {var:15} [{min_val:8.2f}, {max_val:8.2f}] → "
                  f"{n_after:>8,} events ({pct:>5.1f}% removed)")
        else:
            print(f"  {var:15} [MISSING IN DATA]")
    
    return mask


def load_and_process_data(filepath, tree_name="ntp"):
    """Load data and apply cuts."""
    tree = load_data(filepath, tree_name)
    if tree is None:
        return None
    
    # Load required branches
    print(f"\nLoading branches...")
    try:
        data_dict = {}
        for var in PLOT_VARS + list(HMS_CUTS.keys()):
            if var in tree.keys():
                data_dict[var] = tree[var].array(library="np")
            else:
                print(f"  ⚠ {var} not found in tree")
        
        df = pd.DataFrame(data_dict)
        print(f"✓ Loaded {len(df):,} events")
        
        # Apply cuts
        print(f"\nApplying HMS electron cuts:")
        cuts_available = {k: v for k, v in HMS_CUTS.items() if k in df.columns}
        mask = apply_hms_cuts(df, cuts_available)
        
        df_cut = df[mask].reset_index(drop=True)
        print(f"\n✓ After cuts: {len(df_cut):,} events ({100*len(df_cut)/len(df):.2f}% pass)")
        
        return df_cut
    
    except Exception as e:
        print(f"✗ Error processing data: {e}")
        return None


def create_comparison_plots(data_df, sim_df, output_dir="."):
    """Create comparison plots for all variables."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    n_vars = len(PLOT_VARS)
    n_cols = 3
    n_rows = (n_vars + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 5*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    else:
        axes = axes.reshape(-1, 3)
    
    fig.suptitle("Data vs Simulation Comparison (HMS Side Kinematics)", 
                 fontsize=16, fontweight='bold', y=0.995)
    
    for idx, var in enumerate(PLOT_VARS):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        # Skip if variable not in datasets
        if var not in data_df.columns or var not in sim_df.columns:
            ax.text(0.5, 0.5, f'{var}\n(Missing Data)', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=11)
            ax.set_title(f"⚠ {var}", fontsize=11, fontweight='bold')
            continue
        
        # Get data
        data_vals = data_df[var].dropna().values
        sim_vals = sim_df[var].dropna().values
        
        if len(data_vals) == 0 or len(sim_vals) == 0:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=11)
            ax.set_title(f"⚠ {var}", fontsize=11, fontweight='bold')
            continue
        
        # Determine bins
        all_vals = np.concatenate([data_vals, sim_vals])
        vmin, vmax = np.percentile(all_vals, [0.5, 99.5])
        bins = np.linspace(vmin, vmax, 50)
        
        # Plot histograms
        ax.hist(data_vals, bins=bins, alpha=0.6, label=f'Data (n={len(data_vals):,})', 
               color='steelblue', edgecolor='black', linewidth=1, density=True)
        ax.hist(sim_vals, bins=bins, alpha=0.6, label=f'Sim (n={len(sim_vals):,})', 
               color='coral', edgecolor='black', linewidth=1, density=True)
        
        # Statistics
        data_mean = np.mean(data_vals)
        data_std = np.std(data_vals)
        sim_mean = np.mean(sim_vals)
        sim_std = np.std(sim_vals)
        
        # Formatting
        ax.set_xlabel(PLOT_LABELS.get(var, var), fontsize=10, fontweight='bold')
        ax.set_ylabel('Normalized Count', fontsize=10, fontweight='bold')
        ax.set_title(f'{var}\nData: μ={data_mean:.4f} σ={data_std:.4f} | '
                    f'Sim: μ={sim_mean:.4f} σ={sim_std:.4f}', 
                    fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9, loc='best')
    
    # Hide unused subplots
    for idx in range(n_vars, len(axes.flat)):
        axes.flat[idx].set_visible(False)
    
    plt.tight_layout()
    output_file = output_dir / "data_sim_comparison.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\n✓ Comparison plots saved: {output_file}")
    plt.close()


def main():
    """Main analysis script."""
    print("\n" + "="*80)
    print("DATA vs SIMULATION COMPARISON - HMS KINEMATICS")
    print("="*80)
    
    # Load data
    data_df = load_and_process_data(DATA_FILE, tree_name="ntp")
    if data_df is None:
        print("✗ Failed to load data")
        return
    
    # Load simulation
    sim_df = load_and_process_data(SIM_FILE, tree_name="nerd")
    if sim_df is None:
        print("✗ Failed to load simulation")
        return
    
    # Create comparison plots
    print("\n" + "="*80)
    print("Creating comparison plots...")
    print("="*80)
    create_comparison_plots(data_df, sim_df, output_dir=".")
    
    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    for var in PLOT_VARS:
        if var in data_df.columns and var in sim_df.columns:
            data_vals = data_df[var].dropna().values
            sim_vals = sim_df[var].dropna().values
            
            if len(data_vals) > 0 and len(sim_vals) > 0:
                print(f"\n{var}:")
                print(f"  Data: n={len(data_vals):>6,} μ={np.mean(data_vals):>10.6f} σ={np.std(data_vals):>10.6f}")
                print(f"  Sim:  n={len(sim_vals):>6,} μ={np.mean(sim_vals):>10.6f} σ={np.std(sim_vals):>10.6f}")
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
