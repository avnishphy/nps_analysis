#!/usr/bin/env python3
"""Scan xsec bin counts for KinC_x36_5_407. Requires numpy, uproot, matplotlib."""

import argparse
import csv
from datetime import datetime
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[2]
XSEC = REPO / "src/xsec_extract"
DATA = REPO / "output/KinC_x36_5_407/KinC_x36_5/root/combined_branches_LH2.root"
MASS_CUT = DATA.with_name(DATA.stem + "_combined_2d_mass_cut_debug.txt")
SIM = Path("/volatile/hallc/nps/singhav/nps_smearing/smear_x36_5_407/smearing_output/KinC_x36_5_407/root/simc_pi0_analysis_output_smeared.root")
VERTEX = REPO / "output/simc/simc_x36_5_407/worksim"
SETUP = "/group/nps/singhav/setup.csh"
SUMMARY_NAME = "excl_xsec_pi0_analysis_no_simc_model_slice_summary.csv"


def vector(header, name):
    match = re.search(r"\b" + name + r"\s*=\s*\{([^{}]+)\}", header)
    if not match:
        raise ValueError(f"Cannot read {name} from xsec_config.h")
    return [float(value.strip()) for value in match.group(1).split(",")]


def selected_tprime(t_range, q_range, x_range):
    branches = ["t", "tmin", "Q2", "xB", "phi", "pi0_weight", "scale",
                "is_exclusive_ellipse_combined"]
    with uproot.open(DATA) as root:
        data = root["physics"].arrays(branches, library="np")
    tp = data["t"] - data["tmin"]
    keep = (
        (data["is_exclusive_ellipse_combined"] != 0)
        & np.isfinite(tp) & np.isfinite(data["phi"])
        & np.isfinite(data["pi0_weight"] * data["scale"])
        & (tp >= t_range[0]) & (tp <= t_range[-1])
        & (data["Q2"] >= q_range[0]) & (data["Q2"] <= q_range[-1])
        & (data["xB"] >= x_range[0]) & (data["xB"] <= x_range[-1])
    )
    values = tp[keep]
    if len(values) < 6:
        raise ValueError("Too few selected data events for equal-count t-prime bins")
    return values


def replace_vector(header, name, values):
    pattern = r"(\bstd::vector<double>\s+" + name + r"\s*=\s*)\{[^{}]*\};"
    updated, count = re.subn(pattern, lambda match: match.group(1) + "{" + values + "};",
                             header, count=1)
    if count != 1:
        raise ValueError(f"Cannot replace {name} in xsec_config.h")
    return updated


def read_fit(csv_path, n_t):
    with csv_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != n_t or any(row["iq"] != "0" or row["ix"] != "0" for row in rows):
        raise ValueError("Expected one Q2/xB cell and one row per t-prime bin")
    rows.sort(key=lambda row: int(row["it"]))
    valid = [row for row in rows if row["fit_xsec_ok"] == "1"]
    ratios = [float(row["fit_xsec_chi2"]) / float(row["fit_xsec_ndf"])
              for row in valid if float(row["fit_xsec_ndf"]) > 0]
    if ratios and not np.allclose(ratios, ratios[0]):
        raise ValueError("Expected one global chi2/ndf shared by all slices")
    return rows, (ratios[0] if ratios else float("nan")), len(valid)


def coefficient_errors(case):
    """Match plot_coefficient_error: covariance diagonal for each published fit term."""
    indices = {}
    with (case / "migration_parameters.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["region"] == "published" and row["is_nuisance"] == "0":
                indices[(int(row["it"]), row["component"])] = int(row["parameter_index"])

    toy_path = case / "positivity_refit_toy_covariance.csv"
    diagonal = {}
    if toy_path.exists():
        with toy_path.open(newline="") as stream:
            for row in csv.DictReader(stream):
                if row["kind"] == "parameter" and row["index_i"] == row["index_j"]:
                    diagonal[int(row["index_i"])] = float(row["covariance"])
                    successful = int(row["successful_toys"])
        source_label = f"conditional refit-toy SD (N={successful})"
    else:
        with (case / "migration_covariance.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                if row["parameter_i"] == row["parameter_j"]:
                    diagonal[int(row["parameter_i"])] = float(
                        row["stat_plus_mc_covariance"])
        source_label = "marginal statistical + MC error"
        if not any(np.isfinite(value) and value >= 0 for value in diagonal.values()):
            source_label = "unavailable (positivity refit toys did not converge)"

    errors = {}
    for key, index in indices.items():
        variance = diagonal.get(index, float("nan"))
        errors[key] = np.sqrt(variance) if np.isfinite(variance) and variance >= 0 else float("nan")
    return errors, source_label


def plot_structure(rows, n_t, n_phi, ratio, case, pdf):
    errors, source_label = coefficient_errors(case)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.0), constrained_layout=True)
    missing_errors = False
    for ax, (field, component, label, color) in zip(axes, [
        ("sigmaU", "U", r"$\sigma_U$", "black"),
        ("sigmaTT", "TT", r"$\sigma_{TT}$", "tab:blue"),
        ("sigmaTL", "LT", r"$\sigma_{LT}$", "tab:red"),
    ]):
        good = [row for row in rows if row["fit_xsec_ok"] == "1"
                and np.isfinite(float(row["fit_xsec_" + field]))]
        if good:
            # The C++ extractor plots response-weighted generated means, not bin centers.
            x = -np.array([float(row["mean_tprime_vertex_sim"]) for row in good])
            y = 1e9 * np.array([float(row["fit_xsec_" + field]) for row in good])
            yerr = 1e9 * np.array([errors.get((int(row["it"]), component), float("nan"))
                                   for row in good])
            missing_errors |= not np.all(np.isfinite(yerr))
            if component == "U":
                xerr = np.array([
                    [max(0.0, x[i] + float(row["tprime_hi"])) for i, row in enumerate(good)],
                    [max(0.0, -float(row["tprime_lo"]) - x[i]) for i, row in enumerate(good)],
                ])
            else:
                xerr = None
            ax.errorbar(x, y, xerr=xerr, yerr=np.where(np.isfinite(yerr), yerr, 0),
                        fmt="o", color=color, capsize=3)
        ax.axhline(0, color="0.75", linewidth=0.8)
        ax.set(xlabel=r"$-t'_{\mathrm{gen}}$ [GeV$^2$]",
               ylabel=label + r" [nb/GeV$^2$]")
        ax.grid(alpha=0.2)
    fig.suptitle(f"{n_t} t' bins, {n_phi} phi bins; nominal chi2/ndf = {ratio:.3f}")
    note = f"Vertical bars: {source_label}; target systematic separate. U horizontal bars: generated bin spans."
    if missing_errors:
        note += " Unavailable errors have no vertical bar."
    fig.supxlabel(note, fontsize=9)
    pdf.savefig(fig)
    plt.close(fig)

def plot_summary(records, pdf):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    for ax, key, title, limit in [
        (axes[0], "chi2_ndf", "Nominal chi2/ndf", None),
        (axes[1], "fit_fraction", "Successful t' slices / total", (0, 1)),
    ]:
        grid = np.full((4, 11), np.nan)
        for row in records:
            grid[row["n_tprime"] - 3, row["n_phi"] - 10] = row[key]
        image = ax.imshow(np.ma.masked_invalid(grid), origin="lower", aspect="auto",
                          vmin=limit[0] if limit else None,
                          vmax=limit[1] if limit else None)
        ax.set(xlabel="Number of phi bins", ylabel="Number of t' bins", title=title,
               xticks=range(11), yticks=range(4),
               xticklabels=range(10, 21), yticklabels=range(3, 7))
        attempted = {(row["n_tprime"], row["n_phi"]) for row in records}
        for i in range(4):
            for j in range(11):
                label = (f"{grid[i, j]:.2f}" if np.isfinite(grid[i, j]) else
                         "fail" if (i + 3, j + 10) in attempted else "")
                ax.text(j, i, label, ha="center", va="center", fontsize=7)
        fig.colorbar(image, ax=ax)
    pdf.savefig(fig)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tprime-bins", nargs="+", type=int, default=list(range(3, 7)))
    parser.add_argument("--phi-bins", nargs="+", type=int, default=list(range(10, 21)))
    args = parser.parse_args()
    if any(n not in range(3, 7) for n in args.tprime_bins) or any(
            n not in range(10, 21) for n in args.phi_bins):
        parser.error("t-prime bins must be 3-6 and phi bins 10-20")

    original = (XSEC / "xsec_config.h").read_text()
    t_range = vector(original, "tprime_bin_edges")
    q_range = vector(original, "q2_bin_edges")
    xb_match = re.search(r"\bxb_bin_edges_by_q2\s*=\s*\{\s*\{([^{}]+)\}", original)
    if len(q_range) != 2 or not xb_match:
        raise ValueError("This scan expects one Q2 bin and one xB bin")
    x_range = [float(value.strip()) for value in xb_match.group(1).split(",")]
    if len(x_range) != 2:
        raise ValueError("This scan expects one xB bin")
    for required in (DATA, MASS_CUT, SIM, VERTEX):
        if not required.exists():
            raise FileNotFoundError(required)
    values = selected_tprime(t_range, q_range, x_range)
    parent = REPO / "output/KinC_x36_5_407/xsec"
    parent.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix=datetime.now().strftime("scan_%Y%m%d_%H%M%S_"),
                                   dir=parent))
    print(f"{len(values)} selected data events; results: {output}", flush=True)
    records = []
    plot_cases = []

    with tempfile.TemporaryDirectory(prefix="xsec_bin_build_", dir=output) as temporary:
        copied = Path(temporary) / "src/xsec_extract"
        shutil.copytree(XSEC, copied, ignore=shutil.ignore_patterns("*.so", "*.pcm", "*.d"))
        header = copied / "xsec_config.h"
        for n_t in args.tprime_bins:
            edges = np.r_[t_range[0], np.quantile(values, np.arange(1, n_t) / n_t),
                          t_range[-1]]
            if np.any(np.diff(edges) <= 0):
                raise ValueError(f"Collapsed equal-count edges for {n_t} t-prime bins")
            t_cpp = ", ".join(f"{edge:.17g}" for edge in edges)
            for n_phi in args.phi_bins:
                phi_cpp = ", ".join(["0.0"] +
                                    [f"{2 * i}*TMath::Pi()/{n_phi}" for i in range(1, n_phi)] +
                                    ["2*TMath::Pi()"])
                header.write_text(replace_vector(
                    replace_vector(original, "tprime_bin_edges", t_cpp),
                    "phi_bin_edges", phi_cpp))
                case = output / f"t{n_t}_phi{n_phi}"
                case.mkdir()
                command = [
                    str(copied / "run_xsec_pipeline.sh"),
                    "--kin", "KinC_x36_5_407", "--target", "LH2",
                    "--root-dir", str(DATA.parent), "--sim-file", str(SIM),
                    "--vertex_simc_file", str(VERTEX),
                    "--mmiss-lower", "0.6", "--mmiss-upper", "1.1",
                    "--positive-xsec", "--mmiss_select", "ellipse",
                    "--mmiss-cut-file", str(MASS_CUT), "--out-dir", str(case),
                    "--no-diagnostics", "--no-pdf", "--no-png",
                ]
                print(f"Running t'={n_t}, phi={n_phi}", flush=True)
                shell_command = ("source /usr/share/Modules/init/csh; "
                                 f"source {shlex.quote(SETUP)}; exec {shlex.join(command)}")
                with (case / "pipeline.log").open("w") as log:
                    result = subprocess.run(["csh", "-c", shell_command], cwd=REPO,
                                            stdout=log, stderr=subprocess.STDOUT,
                                            check=False)
                row = dict(n_tprime=n_t, n_phi=n_phi, tprime_edges=";".join(map(str, edges)),
                           return_code=result.returncode, chi2_ndf=float("nan"),
                           fit_fraction=0.0)
                summary = case / SUMMARY_NAME
                if result.returncode == 0 and summary.exists():
                    try:
                        fit_rows, row["chi2_ndf"], fitted = read_fit(summary, n_t)
                        row["fit_fraction"] = fitted / n_t
                        plot_cases.append((fit_rows, n_t, n_phi, row["chi2_ndf"], case))
                    except (ValueError, KeyError, OSError) as error:
                        print(f"  Could not read {summary}: {error}", file=sys.stderr)
                else:
                    print(f"  Extraction failed; see {case / 'pipeline.log'}", file=sys.stderr)
                records.append(row)
                with (output / "scan_summary.csv").open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=row.keys())
                    writer.writeheader()
                    writer.writerows(records)

    pdf_path = output / "xsec_bin_numbers.pdf"
    with PdfPages(pdf_path) as pdf:
        plot_summary(records, pdf)
        for fit_rows, n_t, n_phi, ratio, case in plot_cases:
            plot_structure(fit_rows, n_t, n_phi, ratio, case, pdf)
    print(f"Finished. PDF and scan_summary.csv: {output}")
    return int(any(row["return_code"] != 0 or not np.isfinite(row["chi2_ndf"])
                   for row in records))


if __name__ == "__main__":
    sys.exit(main())
