#!/usr/bin/env python3
"""Fit shared pi0 LT/TT and independent U for each exclusive-SIMC setting.

Each setting supplies its own measured yield-per-mC rows and one-mC response.
The rows are stacked; each truth block has one U per setting and shared LT/TT.
U is refitted with LT/TT, not fixed to the earlier single-setting estimate.
The complete Poissonized, within-event MC covariance is carried into the fit.
After fitting, separate U = T + epsilon_nominal*L in each truth block.
"""

# python3 src/xsec_extract/run_joint_xsec_fit.py --prepare-setting xsec_config_x36_5_407.json output/simc/simc_x36_5_407/worksim/ --prepare-setting xsec_config_x36_4.json output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/ --binning-config xsec_config_x36_4.json --mmiss-select ellipse --mmiss-lower 0.6 --mmiss-upper 1.1 --positive-xsec --fit-variance finite-mc --out-dir output/joint_x36_5_407_x36_4_LH2

import argparse
import csv
import hashlib
import json
import math
import re
import shlex
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import uproot

from generate_xsec_config import diamond_vertices, render

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
COMPONENTS = ("U", "LT", "TT")
BASIS = tuple("basis_" + name for name in COMPONENTS)
COV = tuple("cov_" + first + "_" + second
            for first in COMPONENTS for second in COMPONENTS)
VERTEX_EPSILON = "original_exclusive_h10_Q2i_Wi_hsxptari_hsyptari_SIMC_electron_angle"


def need(condition, message):
    if not condition:
        raise ValueError(message)


def same(a, b):
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, abs_tol=1e-12, rel_tol=0)
    return a == b


def finite(record, key, label):
    try:
        value = float(record[key])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{label}: invalid {key}") from error
    need(math.isfinite(value), f"{label}: nonfinite {key}")
    return value


def integer(record, key, label):
    try:
        value = int(record[key])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{label}: invalid {key}") from error
    return value


def rows(path):
    need(path.is_file(), f"missing extraction file: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        need(reader.fieldnames is not None, f"empty CSV: {path}")
        return list(reader)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def config_bins(config):
    phi = ([math.tau * i / config["phi_bins"] for i in range(config["phi_bins"] + 1)]
           if "phi_bins" in config else config["phi_bin_edges"])
    return {"tprime_bin_edges": config["tprime_bin_edges"],
            "t_bin_edges": config["t_bin_edges"],
            "q2_bin_edges": config["q2_bin_edges"],
            "xb_bin_edges_by_q2": config["xb_bin_edges_by_q2"],
            "phi_bin_edges": phi,
            "diamond_xb_q2_vertices": diamond_vertices(config.get("diamond_xb_q2_vertices"))}


def read_config(path):
    path = path.expanduser()
    if not path.is_file() and len(path.parts) == 1:
        path = HERE / "xsec_config" / path.name
    path = path.resolve()
    need(path.is_file(), f"config does not exist: {path}")
    def reject_constant(value):
        raise ValueError(f"nonfinite JSON number {value}")
    config = json.loads(path.read_text(), parse_constant=reject_constant)
    need(isinstance(config, dict), f"config is not an object: {path}")
    render(config)  # Use the same schema and bin validation as single-setting extraction.
    return path, config, config_bins(config)


def nominal_epsilons(configs, overrides=None):
    """Fixed central electron kinematics, independent of event response epsilon."""
    need(overrides is None or len(overrides) == len(configs),
         "repeat --nominal-epsilon once per setting, in setting order")
    result = []
    for index, config in enumerate(configs):
        label = config["configured_kinematic"]
        if overrides is not None:
            value = overrides[index]
            info = {"source": "command_line_override"}
        else:
            need("hms_p_gev" in config,
                 f"{label}: nominal L/T separation needs hms_p_gev in CONFIG "
                 "or --nominal-epsilon for every setting")
            beam = finite(config, "ebeam", label)
            momentum = finite(config, "hms_p_gev", label)
            angle = finite(config, "hms_theta_deg", label)
            need(0 < momentum < beam and 0 < angle < 180,
                 f"{label}: invalid nominal electron kinematics")
            # Ultra-relativistic electron convention: E' = central HMS momentum.
            half_theta = math.radians(angle) / 2
            q2 = 4 * beam * momentum * math.sin(half_theta)**2
            nu = beam - momentum
            value = 1 / (1 + 2 * (1 + nu*nu/q2) * math.tan(half_theta)**2)
            info = {"source": "preset_central_electron_kinematics_massless",
                    "ebeam_gev": beam, "hms_p_gev": momentum,
                    "hms_theta_deg": angle, "q2_gev2": q2, "nu_gev": nu}
        need(math.isfinite(value) and 0 < value < 1,
             f"{label}: nominal epsilon must be finite and in (0,1)")
        result.append(dict(info, setting_index=index, kinematic=label, epsilon=value))
    return result


def write_csv(path, fields, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)


def separate_lt(out_dir, nominal, rank_tolerance):
    """Per-block GLS separation, retaining all cross-block and LT/TT covariance.

    Two settings define an exact linear transformation, even when constrained
    fit covariance is unavailable. More settings require the U covariance for
    GLS weights. No diagnostic curvature is substituted for a covariance.
    """
    parameters = rows(out_dir / "joint_parameters.csv")
    npar = len(parameters)
    need([int(r["parameter_index"]) for r in parameters] == list(range(npar)),
         "joint parameter indices are not contiguous")
    covariance = np.full((npar, npar), np.nan)
    for record in rows(out_dir / "joint_covariance.csv"):
        covariance[int(record["parameter_i"]), int(record["parameter_j"])] = float(
            record["stat_plus_mc_covariance"])
    values = np.array([float(r["value"]) for r in parameters])
    blocks = sorted({int(r["truth_block"]) for r in parameters})
    output, transform, available, diagnostics = [], [], [], []
    for block in blocks:
        indices = [i for i, r in enumerate(parameters) if int(r["truth_block"]) == block]
        u_indices = [i for i in indices if parameters[i]["component"] == "U"]
        settings = [int(parameters[i]["setting_index"]) for i in u_indices]
        epsilon = np.array([nominal[s]["epsilon"] for s in settings])
        design = np.column_stack((np.ones(len(settings)), epsilon))
        status, weights, chi2 = "ok", None, None
        if len(settings) < 2:
            status = "insufficient_settings"
        else:
            singular = np.linalg.svd(design, compute_uv=False)
            if singular[-1] <= rank_tolerance * singular[0]:
                status = "degenerate_epsilon"
            elif len(settings) == 2:
                weights = np.linalg.solve(design, np.eye(2))
                chi2 = 0.0
            else:
                cov_u = covariance[np.ix_(u_indices, u_indices)]
                if not np.isfinite(cov_u).all():
                    status = "unavailable_joint_covariance"
                else:
                    try:
                        chol = np.linalg.cholesky(cov_u)
                        whitening = np.linalg.solve(chol, np.eye(len(settings)))
                        whitened = whitening @ design
                        left, singular, right = np.linalg.svd(whitened, full_matrices=False)
                        if singular[-1] <= rank_tolerance * singular[0]:
                            status = "degenerate_weighted_design"
                        else:
                            weights = (right.T / singular) @ left.T @ whitening
                            residual = whitening @ (values[u_indices] - design @ (weights @ values[u_indices]))
                            chi2 = float(residual @ residual)
                    except np.linalg.LinAlgError:
                        status = "invalid_joint_covariance"
        if weights is not None and not np.isfinite(covariance[np.ix_(u_indices, u_indices)]).all():
            status = "central_only_covariance_unavailable"
        diagnostics.append({"truth_block": block, "setting_indices": settings,
                            "epsilon_span": float(np.ptp(epsilon)) if len(epsilon) else None,
                            "status": status, "chi2": chi2, "ndf": max(0, len(settings)-2)})
        for term, component in enumerate(("T", "L", "LT", "TT")):
            jacobian = np.zeros(npar)
            valid = True
            row_status = status
            if term < 2:
                valid = weights is not None
                if valid:
                    jacobian[u_indices] = weights[term]
            else:
                source = [i for i in indices if parameters[i]["component"] == component]
                need(len(source) == 1, f"missing shared {component} in truth block {block}")
                jacobian[source[0]] = 1
                row_status = "shared_joint_fit"
            labels = {key: parameters[indices[0]][key] for key in ("truth_block", "region", "it", "iq", "ix")}
            output.append(dict(labels, parameter_index=len(output), component=component,
                               value=float(jacobian @ values) if valid else math.nan,
                               error_stat_plus_mc=math.nan, status=row_status))
            transform.append(jacobian)
            available.append(valid)
    transform = np.array(transform)
    separated_cov = np.full((len(output), len(output)), np.nan)
    # Restrict each product to nonzero support: unavailable covariances must
    # not contaminate unrelated coefficients through 0*NaN.
    for i in range(len(output)):
        if not available[i]:
            continue
        ii = np.flatnonzero(transform[i])
        for j in range(i + 1):
            if not available[j]:
                continue
            jj = np.flatnonzero(transform[j])
            subcov = covariance[np.ix_(ii, jj)]
            if np.isfinite(subcov).all():
                value = float(transform[i, ii] @ subcov @ transform[j, jj])
                separated_cov[i, j] = separated_cov[j, i] = value
        variance = separated_cov[i, i]
        if math.isfinite(variance) and variance >= 0:
            output[i]["error_stat_plus_mc"] = math.sqrt(variance)
    write_csv(out_dir / "joint_separated_parameters.csv",
              ("parameter_index", "truth_block", "region", "it", "iq", "ix", "component",
               "value", "error_stat_plus_mc", "status"), output)
    write_csv(out_dir / "joint_separated_covariance.csv",
              ("parameter_i", "parameter_j", "stat_plus_mc_covariance"),
              ({"parameter_i": i, "parameter_j": j, "stat_plus_mc_covariance": separated_cov[i, j]}
               for i in range(len(output)) for j in range(len(output))))
    summary = {"convention": "sigmaU = sigmaT + epsilon_nominal * sigmaL",
               "method": "per_truth_block_GLS_from_joint_U",
               "epsilon_uncertainty": "not_propagated_nominal_values_fixed",
               "positivity": "no_additional_constraints_on_separated_T_or_L",
               "covariance": "full_joint_stat_plus_mc_propagated_including_cross_blocks_and_LT_TT",
               "settings": nominal, "blocks": diagnostics}
    (out_dir / "joint_lt_separation.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def read_metadata(path):
    need(path.is_file(), f"missing input metadata: {path}")
    if path.suffix == ".txt":
        text = path.read_text()
    else:
        with uproot.open(path) as source:
            need("analysis_metadata" in source, f"missing analysis_metadata: {path}")
            text = str(source["analysis_metadata"])
    metadata = dict(line.split("=", 1) for line in text.splitlines()
                    if "=" in line and not line.lstrip().startswith("iq="))
    for line in text.splitlines():
        if ":" in line and not line.lstrip().startswith("iq="):
            key, value = line.split(":", 1)
            metadata[key.strip() + ":"] = value.strip()
    metadata["xb_rows"] = [line.strip() for line in text.splitlines()
                           if line.strip().startswith("iq=")]
    return metadata


def check_saved_bins(metadata, bins, label):
    keys = {"tprime_bin_edges": "tprime_edges:", "q2_bin_edges": "q2_edges:",
            "phi_bin_edges": "phi_edges:"}
    for name, stored_name in keys.items():
        need(stored_name in metadata, f"{label}: missing saved {name}")
        found = [float(x) for x in metadata[stored_name].split()]
        need(same(found, bins[name]), f"{label}: saved {name} differs from config")
    xb_rows = metadata["xb_rows"]
    need(len(xb_rows) == len(bins["xb_bin_edges_by_q2"]),
         f"{label}: saved xB bin row count differs from config")
    for iq, line in enumerate(xb_rows):
        prefix, values = line.split(":", 1)
        need(prefix == f"iq={iq}" and
             same([float(x) for x in values.split()], bins["xb_bin_edges_by_q2"][iq]),
             f"{label}: saved xB bins differ from config in Q2 row {iq}")
    vertices = metadata.get("diamond_xb_q2_vertices")
    need(vertices is not None, f"{label}: missing saved diamond selection")
    numbers = [float(x) for x in re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", vertices)]
    expected = [value for vertex in bins["diamond_xb_q2_vertices"] for value in vertex]
    need(same(numbers, expected), f"{label}: saved diamond selection differs from config")


def load_epsilon_envelope(truth, out_dir, nb, label, prepared=False):
    """Use saved event maxima, including the legacy positivity diagnostics."""
    if prepared:
        maxima = [finite(record, "epsilon_max", label) for record in truth]
        for record, e in zip(truth, maxima):
            populated = integer(record, "events", label) > 0
            need((0 < e < 1) if populated else e == 0,
                 f"{label}: invalid prepared event epsilon maximum")
        return maxima, out_dir / "migration_truth_blocks.csv"
    diagnostics_path = out_dir / "positivity_diagnostics.csv"
    diagnostics = rows(diagnostics_path)
    maxima = [0.0] * nb
    seen = set()
    for index, record in enumerate(diagnostics):
        b = integer(record, "truth_block", label)
        need(0 <= b < nb and b not in seen, f"{label}: invalid positivity block {b}")
        need(integer(record, "active_block_index", label) == index,
             f"{label}: invalid active positivity block order")
        seen.add(b)
        maxima[b] = finite(record, "epsilon_max", label)
        need(0 < maxima[b] < 1, f"{label}: invalid vertex epsilon maximum in block {b}")
    for b, record in enumerate(truth):
        events = integer(record, "events", label)
        need((events > 0) == (b in seen),
             f"{label}: positivity envelope does not cover populated block {b}")
        if "epsilon_max" in record:
            saved = finite(record, "epsilon_max", label)
            need(same(saved, maxima[b]),
                 f"{label}: truth and positivity epsilon maxima differ in block {b}")
    return maxima, diagnostics_path


def checked_input_paths(metadata, config, label):
    """Preserve explicit ROOT provenance; label legacy preset inference."""
    observed = [metadata.get("input_data_file"), metadata.get("input_simc_file")]
    configured = [config["data_file"], config["simc_file"]]
    for kind, raw in zip(("data", "SIMC"), configured):
        need(Path(raw).is_file(), f"{label}: configured {kind} input is missing: {raw}")
    for kind, saved, raw in zip(("data", "SIMC"), observed, configured):
        if saved:
            need(Path(saved).resolve() == Path(raw).resolve(),
                 f"{label}: saved {kind} input differs from current preset")
    status = "recorded_in_root" if all(observed) else "legacy_preset_inference"
    if status != "recorded_in_root":
        print(f"[WARN] {label}: legacy ROOT omits input file paths; "
              "using validated preset paths in joint manifest", file=sys.stderr)
    return tuple(str(Path(value).resolve()) for value in configured), status


def check_setting_flags(settings, fit_variance, positive):
    """Enforce input selection/normalization; earlier fit choices are irrelevant."""
    fields = ("exclusive_selection", "sim_reconstructed_mass_selection",
              "mmiss_lower_gev", "mmiss_upper_gev", "target_contam_factor",
              "target_contam_factor_err", "target_contam_usage",
              "response_absolute_scale", "sim_weight_mode")
    reference = settings[0]["metadata"]
    for setting in settings:
        label, meta = setting["label"], setting["metadata"]
        for field in fields:
            need(field in meta, f"{label}: missing extraction flag {field}")
            left, right = reference[field], meta[field]
            agree = same(float(left), float(right)) if field in (
                "mmiss_lower_gev", "mmiss_upper_gev", "target_contam_factor",
                "target_contam_factor_err") else left == right
            need(agree, f"{label}: extraction flag {field} differs across settings")
        need(meta["target_contam_usage"] == "data_yield_divided_by_factor" and
             meta["sim_weight_mode"] == "full_weight_over_sigcm",
             f"{label}: unsupported data/SIMC normalization")
        need(meta["exclusive_selection"] in ("window", "ellipse", "mcd"),
             f"{label}: unsupported exclusive selection")
        if meta["exclusive_selection"] in ("ellipse", "mcd"):
            need(meta["sim_reconstructed_mass_selection"] == "combined_data_geometry",
                 f"{label}: SIMC mass selection differs from data geometry")
            data = Path(setting["data_file"])
            expected = data.with_name(data.stem + "_combined_2d_mass_cut_debug.txt")
            cut = Path(meta.get("combined_mass_cut_file", ""))
            need(cut.is_file() and cut.resolve() == expected.resolve(),
                 f"{label}: ellipse/MCD geometry is missing or differs from its data file")
    targets = [re.search(r"combined_branches_([^/]+)\.root$", s["data_file"])
               for s in settings]
    need(all(targets) and len({match.group(1) for match in targets}) == 1,
         "joint settings must use the same target")
    return targets[0].group(1)


def load_response_cells(out_dir):
    """Read U/LT/TT moments, accepting either existing response export."""
    cell_path = out_dir / "migration_response_cells.csv"
    if cell_path.is_file():
        cells = rows(cell_path)
        cell_sources = [cell_path]
        response_source = "direct_U_LT_TT_export"
    else:
        # T and U use the same constant angular basis. Select its U/LT/TT
        # covariance submatrix without an epsilon average or T/L separation.
        cell_path = out_dir / "migration_joint_response_cells.csv"
        cells = rows(cell_path)
        for cell in cells:
            cell["basis_U"] = cell["basis_T"]
            for first in COMPONENTS:
                for second in COMPONENTS:
                    a = "T" if first == "U" else first
                    b = "T" if second == "U" else second
                    cell[f"cov_{first}_{second}"] = cell[f"cov_{a}_{b}"]
        cell_sources = [cell_path]
        response_source = "U_LT_TT_submatrix_of_four_component_export"
    return cells, cell_sources, response_source


def load_setting(config_path, config, bins, out_dir):
    label = config["configured_kinematic"]
    out_dir = out_dir.expanduser().resolve()
    need(out_dir.is_dir(), f"{label}: output directory does not exist: {out_dir}")
    preparation = out_dir / "joint_input_metadata.txt"
    prepared = preparation.is_file()
    root_path = preparation if prepared else out_dir / Path(config["out_root"]).name
    metadata = read_metadata(root_path)
    if prepared:
        need(metadata.get("input_stage") == "prepared_joint_inputs_v1" and
             metadata.get("fit_objective") == "not_run", f"{label}: invalid preparation metadata")
    saved_label = metadata.get("configured_kinematic")
    if saved_label is not None:
        need(saved_label == label,
             f"{label}: ROOT output belongs to another kinematic setting")
    else:
        # Older output predates the kinematic-label metadata. Its exclusive
        # vertex file, bin metadata and mass-cut geometry must identify this
        # preset before it can be used as a legacy input.
        token = label.removeprefix("KinC_")
        vertex = Path(metadata.get("vertex_source", ""))
        need(vertex.is_file() and token in vertex.name,
             f"{label}: legacy ROOT cannot be tied to its exclusive SIMC setting")
    need(metadata.get("vertex_epsilon_source") == VERTEX_EPSILON,
         f"{label}: output must use exclusive-SIMC vertex epsilon")
    need(metadata.get("response_absolute_scale") ==
         "physical_SIMC_normalization_only_no_data_area_matching",
         f"{label}: response normalization is incompatible with joint fit")
    need(same(float(metadata.get("vertex_electron_hms_theta_deg", "nan")),
              config["hms_theta_deg"]), f"{label}: saved HMS angle differs from config")
    need(same(float(metadata.get("target_contam_factor", "nan")), config["tgt_contam"]),
         f"{label}: saved target factor differs from config")
    need(same(float(metadata.get("target_contam_factor_err", "nan")),
              config["tgt_contam_err"]),
         f"{label}: saved target-factor error differs from config")
    (data_file, simc_file), provenance_status = checked_input_paths(metadata, config, label)
    if prepared:
        provenance_status = "recorded_in_preparation_metadata"
    check_saved_bins(metadata, bins, label)

    nt = len(bins["tprime_bin_edges"]) - 1
    nq = len(bins["q2_bin_edges"]) - 1
    nx = len(bins["xb_bin_edges_by_q2"][0]) - 1
    nphi = len(bins["phi_bin_edges"]) - 1
    published, nr = nt * nq * nx, nt * nq * nx * nphi
    nb = published + 6
    slice_path = out_dir / ("joint_input_slices.csv" if prepared else Path(config["out_slice_csv"]).name)
    slices = rows(slice_path)
    need(len(slices) == published, f"{label}: wrong slice count")
    for b, record in enumerate(slices):
        it, iq, ix = b // (nq * nx), (b // nx) % nq, b % nx
        need([integer(record, name, label) for name in ("it", "iq", "ix")] ==
             [it, iq, ix], f"{label}: slice order differs at {b}")
        expected = {"tprime_lo": bins["tprime_bin_edges"][it],
                    "tprime_hi": bins["tprime_bin_edges"][it + 1],
                    "q2_lo": bins["q2_bin_edges"][iq],
                    "q2_hi": bins["q2_bin_edges"][iq + 1],
                    "xb_lo": bins["xb_bin_edges_by_q2"][iq][ix],
                    "xb_hi": bins["xb_bin_edges_by_q2"][iq][ix + 1]}
        for name, value in expected.items():
            need(same(finite(record, name, label), value),
                 f"{label}: saved {name} differs from config in slice {b}")

    reco_path = out_dir / "migration_reco_rows.csv"
    reco = rows(reco_path)
    need(len(reco) == nr, f"{label}: expected {nr} reconstructed rows, got {len(reco)}")
    for r, record in enumerate(reco):
        b, ip = divmod(r, nphi)
        it, iq, ix = b // (nq * nx), (b // nx) % nq, b % nx
        need([integer(record, name, label) for name in
              ("reco_row", "it", "iq", "ix", "ip")] == [r, it, iq, ix, ip],
             f"{label}: reconstructed row order differs at {r}")
        for key, expected in (("phi_lo", bins["phi_bin_edges"][ip]),
                              ("phi_hi", bins["phi_bin_edges"][ip + 1])):
            need(same(finite(record, key, label), expected),
                 f"{label}: saved phi edge differs from config in row {r}")
        need(finite(record, "data_variance", label) >= 0,
             f"{label}: negative observed variance in row {r}")
        finite(record, "data", label)

    truth_path = out_dir / "migration_truth_blocks.csv"
    truth = rows(truth_path)
    need(len(truth) == nb, f"{label}: expected {nb} truth blocks, got {len(truth)}")
    for b, record in enumerate(truth):
        need(integer(record, "truth_block", label) == b,
             f"{label}: truth block order differs at {b}")
    epsilon, diagnostics_path = load_epsilon_envelope(truth, out_dir, nb, label, prepared)

    cells, cell_sources, response_source = load_response_cells(out_dir)
    need(len(cells) == nr * nb,
         f"{label}: expected {nr * nb} response cells, got {len(cells)}")
    for i, record in enumerate(cells):
        r, b = divmod(i, nb)
        need(integer(record, "reco_row", label) == r and
             integer(record, "truth_block", label) == b,
             f"{label}: response cell order differs at {i}")
        values = [finite(record, key, label) for key in BASIS + COV]
        for x in range(3):
            need(values[3 + 3*x + x] >= -1e-18,
                 f"{label}: negative response covariance diagonal at ({r},{b})")
            for y in range(x):
                need(math.isclose(values[3 + 3*x + y], values[3 + 3*y + x],
                                  rel_tol=1e-12, abs_tol=1e-24),
                     f"{label}: asymmetric response covariance at ({r},{b})")
    return {"label": label, "config_path": config_path, "config": config,
            "out_dir": out_dir, "root_path": root_path, "metadata": metadata,
            "data_file": data_file, "simc_file": simc_file,
            "provenance_status": provenance_status, "response_source": response_source,
            "reco": reco, "cells": cells, "epsilon": epsilon,
            "sources": [root_path, slice_path, reco_path, truth_path,
                        diagnostics_path] + cell_sources}


def write_problem(path, settings, bins):
    nt = len(bins["tprime_bin_edges"]) - 1
    nq = len(bins["q2_bin_edges"]) - 1
    nx = len(bins["xb_bin_edges_by_q2"][0]) - 1
    nphi = len(bins["phi_bin_edges"]) - 1
    nb = nt * nq * nx + 6
    nr = sum(len(setting["reco"]) for setting in settings)
    with path.open("w") as output:
        output.write(f"joint_xsec_v3 {len(settings)} {nt} {nq} {nx} {nphi} {nr} {nb}\n")
        for setting in settings:
            output.write(" ".join(format(e, ".17g") for e in setting["epsilon"]) + "\n")
        for si, setting in enumerate(settings):
            for r, record in enumerate(setting["reco"]):
                output.write(f"{si} {r} {record['data']} {record['data_variance']}")
                for b in range(nb):
                    cell = setting["cells"][r * nb + b]
                    output.write(" " + " ".join(cell[key] for key in BASIS + COV))
                output.write("\n")


def root_environment():
    """Load Hall C once and execute commands as argv, without shell path expansion."""
    output = subprocess.check_output([
        "csh", "-c", "source /usr/share/Modules/init/csh; "
        "source /group/nps/singhav/setup.csh; /usr/bin/env -0"], cwd=REPO)
    # setup prints a banner before env; the first environment record is
    # recovered from the last line before the first NUL separator.
    records = output.split(b"\0")
    records[0] = records[0].splitlines()[-1]
    return dict(record.decode().split("=", 1) for record in records if b"=" in record)


def prepare_and_fit(args):
    """Prepare each raw setting, persist reusable inputs, then invoke only the joint fit."""
    out_dir = args.out_dir.expanduser().resolve()
    inputs_dir = (args.inputs_dir or out_dir.with_name(out_dir.name + "_inputs")).expanduser().resolve()
    need(not out_dir.exists(), f"joint output directory already exists: {out_dir}")
    need(not inputs_dir.exists(), f"preparation directory already exists: {inputs_dir}")
    need(inputs_dir != out_dir and out_dir not in inputs_dir.parents and
         inputs_dir not in out_dir.parents, "fit and preparation directories must be separate")
    common = read_config(args.binning_config)[2] if args.binning_config else None
    prepared_configs = []
    for raw_config, raw_vertex in args.prepare_setting:
        original, config, bins = read_config(Path(raw_config))
        config = dict(config)
        if common is not None:
            config.pop("phi_bins", None)
            config.update(common)
            # Normalized corners are tuples; the JSON config schema needs lists.
            config["diamond_xb_q2_vertices"] = [list(point) for point in
                                                common["diamond_xb_q2_vertices"]] or None
        for arg, key in ((args.mmiss_select, "mmiss_select"),
                         (args.mmiss_lower, "mmiss_lower_gev"),
                         (args.mmiss_upper, "mmiss_upper_gev")):
            if arg is not None:
                config[key] = arg
        label = config["configured_kinematic"]
        need(re.fullmatch(r"[A-Za-z0-9_-]+", label), f"unsafe kinematic label: {label}")
        for key in ("data_file", "simc_file"):
            config[key] = str(Path(config[key]).expanduser().resolve())
            need(Path(config[key]).is_file(), f"missing {key}: {config[key]}")
        vertex = Path(raw_vertex).expanduser().resolve()
        if vertex.is_dir():
            vertex /= "nps_excl_pi0_" + label.removeprefix("KinC_") + ".root"
        need(vertex.is_file(), f"missing exclusive SIMC vertex file: {vertex}")
        if config["mmiss_select"] != "window":
            data = Path(config["data_file"])
            cut = data.with_name(data.stem + "_combined_2d_mass_cut_debug.txt")
            need(cut.is_file(), f"missing combined-data mass geometry: {cut}")
        render(config)
        bins = config_bins(config)
        if prepared_configs:
            reference = prepared_configs[0][1]
            for key, value in config_bins(reference).items():
                need(same(value, bins[key]), f"bin mismatch in {key}; use --binning-config for common bins")
            for key in ("mmiss_select", "mmiss_lower_gev", "mmiss_upper_gev", "tgt_contam", "tgt_contam_err"):
                need(same(reference[key], config[key]), f"preparation settings differ in {key}")
        prepared_configs.append((original, config, vertex))
    for key in ("configured_kinematic", "data_file"):
        values = [config[key] for _, config, _ in prepared_configs]
        need(len(set(values)) == len(values), f"duplicate {key} in preparation")
    nominal_epsilons([config for _, config, _ in prepared_configs], args.nominal_epsilon)
    targets = [re.search(r"combined_branches_([^/]+)\.root$", c["data_file"])
               for _, c, _ in prepared_configs]
    need(all(targets) and len({m.group(1) for m in targets}) == 1,
         "joint settings must use the same target")
    env = root_environment()
    flags = shlex.split(subprocess.check_output(["root-config", "--cflags", "--libs"], env=env, text=True))
    (inputs_dir / "configs").mkdir(parents=True)
    pairs, provenance = [], []
    for original, config, vertex in prepared_configs:
        label = config["configured_kinematic"]
        saved_config = inputs_dir / "configs" / (label + ".json")
        saved_config.write_text(json.dumps(config, indent=2) + "\n")
        destination = inputs_dir / label
        with tempfile.TemporaryDirectory(prefix="pi0_joint_prepare_") as temporary:
            build = Path(temporary)
            (build / "xsec_config.h").write_text(render(config))
            binary = build / "prepare_inputs"
            subprocess.run(["g++", "-std=c++17", "-O2", "-I" + str(build), "-I" + str(HERE),
                            str(HERE / "excl_xsec_pi0_analysis_no_simc_model.C"),
                            *flags, "-lMinuit2", "-o", str(binary)], env=env, cwd=REPO, check=True)
            command = [str(binary), "--prepare-joint-inputs", "--kin", label,
                       "--data-file", config["data_file"], "--sim-file", config["simc_file"],
                       "--vertex_simc_file", str(vertex), "--out-dir", str(destination)]
            if config["mmiss_select"] != "window":
                command += ["--mmiss_select", config["mmiss_select"]]
            print(f"[PREPARE] {label}: observations and response only", flush=True)
            subprocess.run(command, env=env, cwd=REPO, check=True)
        pairs += ["--setting", str(saved_config), str(destination)]
        provenance.append({"kinematic": label, "original_config": str(original),
                           "original_config_sha256": digest(original), "effective_config": str(saved_config),
                           "vertex_simc_file": str(vertex), "output_dir": str(destination)})
    (inputs_dir / "preparation_manifest.json").write_text(json.dumps({
        "method": "event_accumulation_without_individual_fits", "settings": provenance,
        "common_binning": common}, indent=2) + "\n")
    print(f"[PREPARED] Reusable inputs written to {inputs_dir}", flush=True)
    command = pairs + ["--out-dir", str(out_dir), "--fit-variance", args.fit_variance,
                       "--rank-tolerance", str(args.rank_tolerance),
                       "--mc-max-iterations", str(args.mc_max_iterations),
                       "--mc-fit-tolerance", str(args.mc_fit_tolerance)]
    if args.positive_xsec:
        command.append("--positive-xsec")
    for epsilon in args.nominal_epsilon or []:
        command += ["--nominal-epsilon", str(epsilon)]
    if args.no_plots:
        command.append("--no-plots")
    if args.partons:
        command += ["--partons", "--partons-warmups", str(args.partons_warmups),
                    "--partons-calls", str(args.partons_calls)]
    return main(command)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--setting", nargs=2, metavar=("CONFIG", "OUTPUT_DIR"), action="append",
                        help="reuse prepared inputs or previous extraction exports; repeat at least twice")
    source.add_argument("--prepare-setting", nargs=2, metavar=("CONFIG", "VERTEX_SIMC"), action="append",
                        help="read data/SIMC paths from CONFIG and prepare inputs without individual fits")
    source.add_argument("--plot-only", type=Path, metavar="JOINT_OUTPUT",
                        help="add the complete plot report to an existing joint fit without refitting")
    parser.add_argument("--inputs-dir", type=Path, help="new preparation directory (default: OUT_DIR_inputs)")
    parser.add_argument("--binning-config", type=Path, help="use this preset's bins/diamond for all prepared settings")
    parser.add_argument("--mmiss-select", choices=("window", "ellipse", "mcd"), help="preparation selection override")
    parser.add_argument("--mmiss-lower", type=float, help="preparation missing-mass lower bound")
    parser.add_argument("--mmiss-upper", type=float, help="preparation missing-mass upper bound")
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--no-plots", action="store_true", help="skip the default PDF/PNG plot report")
    parser.add_argument("--partons", action="store_true", help="include native GK06/GPDGK19 comparison plots")
    parser.add_argument("--partons-warmups", type=int, default=10000)
    parser.add_argument("--partons-calls", type=int, default=100000)
    parser.add_argument("--fit-variance", choices=("data", "finite-mc"), default="finite-mc")
    parser.add_argument("--positive-xsec", action="store_true")
    parser.add_argument("--nominal-epsilon", type=float, action="append", metavar="EPSILON",
                        help="override nominal L/T separation epsilon; repeat once per setting in order")
    parser.add_argument("--rank-tolerance", type=float, default=1e-10)
    parser.add_argument("--mc-max-iterations", type=int, default=30)
    parser.add_argument("--mc-fit-tolerance", type=float, default=1e-6)
    args = parser.parse_args(argv)
    try:
        need(args.partons_warmups > 0 and args.partons_calls > 0, "PARTONS call counts must be positive")
        need(not (args.no_plots and args.partons), "--partons requires plots")
        if args.plot_only:
            need(args.out_dir is None and not args.no_plots and args.nominal_epsilon is None and
                 all(v is None for v in (args.inputs_dir, args.binning_config, args.mmiss_select,
                                         args.mmiss_lower, args.mmiss_upper)),
                 "--plot-only uses the saved fit and cannot change inputs, epsilon, or output directory")
            from joint_xsec_plots import render as render_plots
            render_plots(args.plot_only, args.partons, args.partons_warmups, args.partons_calls)
            return 0
        need(args.out_dir is not None, "--out-dir is required for fitting")
        need(len(args.setting or args.prepare_setting) >= 2, "supply at least two setting pairs")
        need(math.isfinite(args.rank_tolerance) and 0 < args.rank_tolerance < 1,
             "rank tolerance must be in (0,1)")
        need(args.mc_max_iterations > 0 and math.isfinite(args.mc_fit_tolerance) and
             0 < args.mc_fit_tolerance < 1, "invalid MC convergence controls")
        if args.prepare_setting:
            return prepare_and_fit(args)
        need(all(value is None for value in (args.inputs_dir, args.binning_config,
                 args.mmiss_select, args.mmiss_lower, args.mmiss_upper)),
             "preparation options require --prepare-setting; saved input bins/cuts cannot be changed")
        out_dir = args.out_dir.expanduser().resolve()
        need(not out_dir.exists(), f"joint output directory already exists: {out_dir}")
        configs = [read_config(Path(pair[0])) for pair in args.setting]
        nominal = nominal_epsilons([config for _, config, _ in configs], args.nominal_epsilon)
        reference = configs[0][2]
        for path, config, bins in configs[1:]:
            for name in reference:
                need(same(reference[name], bins[name]),
                     f"bin mismatch in {name}: {configs[0][0]} versus {path}")
        settings = [load_setting(path, config, bins, Path(pair[1]))
                    for (path, config, bins), pair in zip(configs, args.setting)]
        for key, description in (("label", "kinematic label"),
                                 ("out_dir", "extraction output")):
            values = [str(s[key]) for s in settings]
            need(len(values) == len(set(values)), f"duplicate {description} in joint input")
        data_files = [s["data_file"] for s in settings]
        need(len(set(data_files)) == len(data_files),
             "duplicate data file provenance in joint inputs")
        target = check_setting_flags(settings, args.fit_variance, args.positive_xsec)
        with tempfile.TemporaryDirectory(prefix="pi0_joint_xsec_") as temporary:
            temp = Path(temporary)
            (temp / "xsec_config.h").write_text(render(configs[0][1]))
            problem = temp / "problem.txt"
            write_problem(problem, settings, reference)
            compile_command = ("source /usr/share/Modules/init/csh; "
                               "source /group/nps/singhav/setup.csh; "
                               f"g++ -std=c++17 -O2 -I{temp} -I{HERE} "
                               f"{HERE / 'xsec_joint_solver.C'} "
                               "`root-config --cflags --libs` "
                               f"-o {temp / 'joint_solver'}")
            subprocess.run(["csh", "-c", compile_command], cwd=REPO, check=True)
            result = temp / "result"
            result.mkdir()
            subprocess.run(["csh", "-c", "source /usr/share/Modules/init/csh; "
                            "source /group/nps/singhav/setup.csh; "
                            f"{temp / 'joint_solver'} {problem} {result} "
                            f"{args.fit_variance} {int(args.positive_xsec)} "
                            f"{args.rank_tolerance:.17g} {args.mc_max_iterations} "
                            f"{args.mc_fit_tolerance:.17g}"], cwd=REPO, check=True)
            separation = separate_lt(result, nominal, args.rank_tolerance)
            manifest = {
                "lt_separation": separation,
                "method": "shared_LT_TT_independent_U_forward_response_fit",
                "parameter_contract": {
                    "shared_components": ["LT", "TT"],
                    "per_setting_components": ["U"],
                    "setting_index": "zero-based --setting order; -1 denotes shared",
                    "U": "independently refitted per setting and truth block",
                    "unsupported_guard_U": "omitted when that setting has no fitted-row response",
                    "positivity": "per-setting event epsilon maximum in each supported truth block"},
                "normalization": "each setting yield per mC and one-mC SIMC response as separate rows",
                "target_factor_uncertainty": "not propagated",
                "target": target,
                "selection": {key: settings[0]["metadata"][key] for key in (
                    "exclusive_selection", "sim_reconstructed_mass_selection",
                    "mmiss_lower_gev", "mmiss_upper_gev", "target_contam_factor",
                    "target_contam_factor_err")},
                "joint_partons_projection": "not_computed",
                "bins": reference,
                "settings": [{"kinematic": s["label"], "config": str(s["config_path"]),
                              "config_sha256": digest(s["config_path"]),
                              "output_dir": str(s["out_dir"]),
                              "input_data_file": s["data_file"],
                              "input_simc_file": s["simc_file"],
                              "input_provenance": s["provenance_status"],
                              "response_source": s["response_source"],
                              "input_stage": s["metadata"].get("input_stage", "individual_extraction_export"),
                              "upstream_fit_objective": s["metadata"].get("fit_objective", "gaussian"),
                              "combined_mass_cut_file": s["metadata"].get("combined_mass_cut_file"),
                              "source_sha256": {str(p): digest(p) for p in s["sources"]}}
                             for s in settings],
                "options": {"fit_objective": "gaussian", "fit_variance": args.fit_variance,
                            "positive_xsec": args.positive_xsec,
                            "rank_tolerance": args.rank_tolerance,
                            "mc_max_iterations": args.mc_max_iterations,
                            "mc_fit_tolerance": args.mc_fit_tolerance}}
            (result / "joint_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
            out_dir.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(result, out_dir)
        print(f"Joint fit written to {out_dir}")
        statuses = {}
        for block in separation["blocks"]:
            statuses[block["status"]] = statuses.get(block["status"], 0) + 1
        print("[L/T] Truth-block separation status: " +
              ", ".join(f"{status}={count}" for status, count in statuses.items()))
        if not args.no_plots:
            from joint_xsec_plots import render as render_plots
            render_plots(out_dir, args.partons, args.partons_warmups, args.partons_calls)
    except (ValueError, OSError, subprocess.CalledProcessError, KeyError) as error:
        print(f"[FATAL] joint xsec fit: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
