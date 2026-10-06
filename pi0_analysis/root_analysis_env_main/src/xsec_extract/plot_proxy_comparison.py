"""Diagnostic-only overlay; no independent-fit points enter the model fit.

Usage: python3 plot_proxy_comparison.py MODEL_DIR REFERENCE_SLICE_CSV OUTPUT.pdf
Requires numpy and matplotlib. Internal cross sections are ub/MeV^2.
"""
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

directory, reference, output = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
with (directory / "model_structure_functions.csv").open() as f:
    points = list(csv.DictReader(f))
with reference.open() as f:
    independent = list(csv.DictReader(f))
x = np.array([float(r["tprime"]) for r in points])
grid = np.linspace(min(x), max(x), 250)
with (directory / "model_structure_curves.csv").open() as f:
    evaluated = list(csv.DictReader(f))
fig, axes = plt.subplots(1, 3, figsize=(13, 4))
for k, (term, old_term) in enumerate(zip(("U", "LT", "TT"), ("U", "TL", "TT"))):
    ax = axes[k]
    samples=[r for r in evaluated if r['component']==term]
    grid=np.array([float(r['tprime']) for r in samples]);value=np.array([float(r['value']) for r in samples])
    err=np.array([float(r['error']) for r in samples])
    ax.plot(grid, value*1e9, '--', lw=.9, label="Reference kinematics - diagnostic only")
    values = np.array([float(r[f"sigma_{term}"]) for r in points])
    errors = np.array([float(r[f"sigma_{term}_error"]) for r in points])
    ax.errorbar(x, values*1e9, yerr=np.where(np.isfinite(errors), errors, 0)*1e9, fmt="o", label="Event-averaged fitted model")
    rows = [r for r in independent if r["fit_xsec_ok"] == "1"]
    rx = np.array([float(r["mean_tprime_vertex_sim"]) for r in rows])
    ry = np.array([float(r[f"fit_xsec_sigma{old_term}"]) for r in rows])
    re = np.array([float(r[f"fit_xsec_sigma{old_term}err"]) for r in rows])
    ax.errorbar(rx, ry*1e9, yerr=np.where(np.isfinite(re), re, 0)*1e9, fmt="s",
                label="Independent reference" if np.isfinite(re).all() else "Reference (errors unavailable)")
    with (directory/'model_kinematic_ranges.csv').open() as f:
        accepted=next(r for r in csv.DictReader(f) if r['quantity']=='tprime')
    ax.axvspan(float(accepted['min']),float(accepted['max']),color='gray',alpha=.08)
    ax.set(xlabel="Signed t' [GeV^2]", ylabel=f"sigma {term} [nb/GeV^2]")
    ax.grid(alpha=.2)
axes[0].legend(fontsize=7)
with (directory/'model_context.csv').open() as f:
    context=next(csv.DictReader(f))
fig.suptitle(('SYNTHETIC VALIDATION | ' if context['synthetic']=='1' else '')+'Event-averaged model versus no-model extraction')
fig.tight_layout(rect=(0,0,1,.92))
fig.savefig(output)
print(output)
