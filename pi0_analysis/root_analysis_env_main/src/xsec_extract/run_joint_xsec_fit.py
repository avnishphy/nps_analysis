#!/usr/bin/env python3
"""Joint event-level SigParam M0 fit for exclusive-SIMC epsilon settings.

Every setting has its own U normalization and U slope.  LT/TT model
normalizations are shared.  The low-tprime feed-in has one U per setting and
shared LT/TT; Q2/xB feed-in remains a fixed nominal event-model contribution.
"""

# python3 src/xsec_extract/run_joint_xsec_fit.py --prepare-setting xsec_config_x36_5_407.json output/simc/simc_x36_5_407/worksim/ --prepare-setting xsec_config_x36_4.json output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/ --binning-config xsec_config_x36_4.json --mmiss-select ellipse --mmiss-lower 0.6 --mmiss-upper 1.1 --positive-xsec --fit-variance finite-mc --out-dir output/joint_x36_5_407_x36_4_LH2

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shlex
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

# Replica workers are process-parallel.  Keep BLAS/OpenMP single-threaded to
# avoid multiplying each worker by the host library's default thread count.
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"

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
    need(config.get("model_identifier") == "sigparam2021_pi0",
         f"{path}: joint M0 requires model_identifier=sigparam2021_pi0")
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
    need(prepared, f"{label}: regenerate prepared joint M0 inputs")
    root_path = preparation if prepared else out_dir / Path(config["out_root"]).name
    metadata = read_metadata(root_path)
    if prepared:
        need(metadata.get("input_stage") == "prepared_joint_m0_inputs_v3" and
             metadata.get("fit_objective") == "not_run" and
             metadata.get("joint_model") == "sigparam2021_pi0_event_level_M0",
             f"{label}: invalid or obsolete joint M0 preparation metadata")
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
        need(finite(record, "fixed_feedin_mc_variance", label) >= 0,
             f"{label}: negative fixed-feed-in MC variance in row {r}")
        finite(record, "fixed_feedin_prediction", label)

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
    event_path = out_dir / "joint_model_events.csv"
    event_rows = rows(event_path)
    need(event_rows, f"{label}: empty joint M0 event cache")
    required = {"event_index", "reco_row", "truth_block", "treatment", "response_weight",
                "tau", "epsilon", "baseline_U", "baseline_LT", "baseline_TT",
                "basis_U", "basis_LT", "basis_TT"}
    need(required.issubset(event_rows[0]), f"{label}: incomplete joint M0 event schema")
    for index, record in enumerate(event_rows):
        need(integer(record, "event_index", label) == index,
             f"{label}: noncontiguous joint event index at {index}")
        need(0 <= integer(record, "reco_row", label) < nr,
             f"{label}: event outside reconstructed rows")
        need(0 <= integer(record, "truth_block", label) < nb,
             f"{label}: event outside truth blocks")
        need(record["treatment"] in ("physics_model", "fitted_tprime_feedin", "fixed_model_feedin"),
             f"{label}: invalid joint event treatment")
        for key in ("response_weight", "tau", "epsilon", "baseline_U", "baseline_LT",
                    "baseline_TT", "basis_U", "basis_LT", "basis_TT"):
            finite(record, key, label)
    return {"label": label, "config_path": config_path, "config": config,
            "out_dir": out_dir, "root_path": root_path, "metadata": metadata,
            "data_file": data_file, "simc_file": simc_file,
            "provenance_status": provenance_status, "response_source": response_source,
            "reco": reco, "cells": cells, "epsilon": epsilon,
            "event_path": event_path,
            "sources": [root_path, slice_path, reco_path, truth_path,
                        diagnostics_path, event_path] + cell_sources}


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
                            str(HERE / "excl_xsec_pi0_analysis_simc_model.C"),
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
                       "--fit-strategy", args.fit_strategy, "--model-starts", str(args.model_starts),
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
    command += ["--fit-objective", args.fit_objective]
    if args.publish_calibrated_release:
        command += ["--publish-calibrated-release", "--toy-jobs", str(args.toy_jobs),
                    "--toy-seed", str(args.toy_seed)]
        if args.reuse_toys:
            command.append("--reuse-toys")
        for flag,value in (("--calibration-source",args.calibration_source),
                           ("--calibrated-out-dir",args.calibrated_out_dir),
                           ("--final-report-pdf",args.final_report_pdf)):
            if value is not None:command += [flag,str(value)]
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
    parser.add_argument("--fit-objective", choices=("gaussian", "scaled-poisson"), default="gaussian")
    parser.add_argument("--fit-strategy", choices=("staged_feasible",), default="staged_feasible",
                        help="joint M0 active-set constrained solver")
    parser.add_argument("--model-starts", type=int, default=6)
    positivity=parser.add_mutually_exclusive_group()
    positivity.add_argument("--positive-xsec", action="store_true")
    positivity.add_argument("--no-positive-xsec", dest="positive_xsec", action="store_false")
    parser.set_defaults(positive_xsec=False)
    parser.add_argument("--publish-calibrated-release", action="store_true",
                        help="generate/reuse 500 joint physical-event toys and publish calibrated intervals")
    parser.add_argument("--reuse-toys", action="store_true",
                        help="reuse an exactly matching --calibration-source instead of generating fresh toys")
    parser.add_argument("--toy-jobs", type=int, default=0, help="toy worker processes; 0 uses the affinity mask")
    parser.add_argument("--toy-seed", type=int, default=20261007)
    parser.add_argument("--calibration-source", type=Path)
    parser.add_argument("--calibrated-out-dir", type=Path)
    parser.add_argument("--final-report-pdf", type=Path)
    parser.add_argument("--nominal-epsilon", type=float, action="append", metavar="EPSILON",
                        help="override nominal L/T separation epsilon; repeat once per setting in order")
    parser.add_argument("--rank-tolerance", type=float, default=1e-10)
    parser.add_argument("--mc-max-iterations", type=int, default=30)
    parser.add_argument("--mc-fit-tolerance", type=float, default=1e-6)
    args = parser.parse_args(argv)
    try:
        need(args.partons_warmups > 0 and args.partons_calls > 0, "PARTONS call counts must be positive")
        need(not (args.no_plots and args.partons), "--partons requires plots")
        need(args.toy_jobs >= 0, "--toy-jobs must be zero or positive")
        need(args.toy_seed >= 0, "--toy-seed must be nonnegative")
        if args.plot_only:
            need(args.out_dir is None and not args.no_plots and args.nominal_epsilon is None and
                 not args.publish_calibrated_release and not args.reuse_toys and
                 all(value is None for value in (args.calibration_source,args.calibrated_out_dir,
                                                  args.final_report_pdf)) and
                 all(v is None for v in (args.inputs_dir, args.binning_config, args.mmiss_select,
                                         args.mmiss_lower, args.mmiss_upper)),
                 "--plot-only uses the saved fit and cannot change inputs, release, epsilon, or output directory")
            from joint_xsec_plots import render as render_plots
            render_plots(args.plot_only, args.partons, args.partons_warmups, args.partons_calls)
            return 0
        need(not args.partons,
             "--partons is not implemented for the joint event-level M0 report")
        need(args.fit_objective == "gaussian",
             "joint event-level M0 currently supports --fit-objective gaussian only")
        need(not args.reuse_toys or args.publish_calibrated_release,
             "--reuse-toys requires --publish-calibrated-release")
        need(args.publish_calibrated_release or all(value is None for value in (
             args.calibration_source,args.calibrated_out_dir,args.final_report_pdf)),
             "calibrated-release paths require --publish-calibrated-release")
        if args.publish_calibrated_release:
            need(args.positive_xsec and args.fit_strategy == "staged_feasible" and
                 args.fit_variance == "finite-mc" and not args.no_plots,
                 "calibrated release requires staged_feasible, finite-mc, positive-xsec, and plots")
            need(args.reuse_toys == (args.calibration_source is not None),
                 "use --calibration-source together with --reuse-toys; fresh campaigns choose a new source")
        need(args.out_dir is not None, "--out-dir is required for fitting")
        need(len(args.setting or args.prepare_setting) >= 2, "supply at least two setting pairs")
        need(math.isfinite(args.rank_tolerance) and 0 < args.rank_tolerance < 1,
             "rank tolerance must be in (0,1)")
        need(args.mc_max_iterations > 0 and math.isfinite(args.mc_fit_tolerance) and
             0 < args.mc_fit_tolerance < 1, "invalid MC convergence controls")
        need(args.model_starts > 0, "--model-starts must be positive")
        if args.prepare_setting:
            return prepare_and_fit(args)
        need(all(value is None for value in (args.inputs_dir, args.binning_config,
                 args.mmiss_select, args.mmiss_lower, args.mmiss_upper)),
             "preparation options require --prepare-setting; saved input bins/cuts cannot be changed")
        out_dir = args.out_dir.expanduser().resolve()
        need(not out_dir.exists(), f"joint output directory already exists: {out_dir}")
        if args.publish_calibrated_release:
            early_release=(args.calibrated_out_dir or out_dir.with_name(out_dir.name+"_calibrated")).expanduser().resolve()
            early_pdf=(args.final_report_pdf.expanduser().resolve() if args.final_report_pdf else
                       early_release/"joint_preliminary_cross_section_report.pdf")
            need(not early_release.exists(),f"refusing existing calibrated release: {early_release}")
            need(early_pdf.parent==early_release,
                 "--final-report-pdf must be directly inside --calibrated-out-dir")
            if args.reuse_toys:
                early_campaign=args.calibration_source.expanduser().resolve()
                need(early_campaign.is_dir(),f"calibration campaign does not exist: {early_campaign}")
                for relative in ("campaign_manifest.json","central.npz","central.json",
                                 "toys/summary.json","toys/replicas.npz"):
                    need((early_campaign/relative).is_file(),
                         f"calibration campaign is incomplete: {early_campaign/relative}")
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
        if args.publish_calibrated_release:
            need(settings[0]["metadata"]["exclusive_selection"] == "ellipse",
                 "fresh/reused joint calibrated release requires the frozen ellipse selection")
            if args.reuse_toys:
                from joint_m0_release import check_campaign
                check_campaign(args.calibration_source.expanduser().resolve(),settings,reference)
        with tempfile.TemporaryDirectory(prefix="pi0_joint_xsec_") as temporary:
            result = Path(temporary) / "result"
            from joint_m0_solver import fit_and_write
            summary, _, _ = fit_and_write(settings, reference, result,
                variance_mode=args.fit_variance, positive=args.positive_xsec,
                rank_tolerance=args.rank_tolerance, max_iterations=args.mc_max_iterations,
                tolerance=args.mc_fit_tolerance, starts=args.model_starts)
            manifest = {
                "method": "joint_event_level_sigparam2021_M0",
                "parameter_contract": {
                    "per_setting_model_parameters": ["N_U", "DeltaB_U"],
                    "shared_model_parameters": ["N_LT", "N_TT"],
                    "per_setting_tprime_feedin": ["U"],
                    "shared_tprime_feedin": ["LT", "TT"],
                    "fixed_feedin": ["q2_below", "q2_above", "xb_below", "xb_above", "tprime_above"],
                    "positivity": "physical event epsilon plus per-setting low-tprime epsilon envelope"},
                "normalization": "each setting yield per mC and event-level one-mC SIMC model response",
                "nominal_epsilon_context": nominal,
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
                            "fit_strategy": args.fit_strategy, "model_starts": args.model_starts,
                            "positive_xsec": args.positive_xsec,
                            "rank_tolerance": args.rank_tolerance,
                            "mc_max_iterations": args.mc_max_iterations,
                            "mc_fit_tolerance": args.mc_fit_tolerance}}
            (result / "joint_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
            out_dir.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(result, out_dir)
        print(f"Joint fit written to {out_dir}")
        print(f"[JOINT_M0] parameters={summary['parameters']} rank={summary['rank']} "
              f"chi2/ndf={summary['chi2']}/{summary['ndf']} boundary={summary['positivity_boundary_active']}")
        if not args.no_plots:
            from joint_xsec_plots import render as render_plots
            render_plots(out_dir, args.partons, args.partons_warmups, args.partons_calls)
        release_pdf=None
        if args.publish_calibrated_release:
            from joint_m0_release import generate_campaign, publish_release
            release_dir=(args.calibrated_out_dir or out_dir.with_name(out_dir.name+"_calibrated")).expanduser().resolve()
            final_pdf=(args.final_report_pdf.expanduser().resolve() if args.final_report_pdf else
                       release_dir/"joint_preliminary_cross_section_report.pdf")
            if args.reuse_toys:
                campaign=args.calibration_source.expanduser().resolve()
            else:
                tag=time.strftime("%Y%m%dT%H%M%S",time.gmtime())+f"_{os.getpid()}"
                campaign=out_dir/"toy_campaigns"/tag
                environment=root_environment();os.environ.update({key:value for key,value in environment.items()
                                                                   if key in ("PATH","LD_LIBRARY_PATH","ROOTSYS")})
                print(f"[TOYS] generating 500 accepted joint physical-event toys at {campaign}",flush=True)
                generate_campaign(campaign,settings,reference,jobs=args.toy_jobs,seed=args.toy_seed,
                    variance_mode=args.fit_variance,starts=args.model_starts,
                    rank_tolerance=args.rank_tolerance,max_iterations=args.mc_max_iterations,
                    tolerance=args.mc_fit_tolerance,environment=environment)
            print(f"[RELEASE] validating and publishing calibrated joint intervals to {release_dir}",flush=True)
            release_pdf=publish_release(campaign,release_dir,settings,reference,
                                        fit_output=out_dir,final_pdf=final_pdf)
            manifest=json.loads((out_dir/"joint_manifest.json").read_text())
            manifest["calibrated_release"]={"status":"complete","campaign":str(campaign),
                "output":str(release_dir),"final_report_pdf":str(release_pdf),"accepted_toys":500}
            manifest["target_factor_uncertainty"]="separate_fully_correlated_scale_covariance_in_calibrated_release"
            (out_dir/"joint_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
        artifacts=[out_dir/name for name in ("joint_xsec_output.root","joint_model_parameters.csv",
            "joint_model_covariance.csv","joint_structure_functions.csv","joint_rows.csv",
            "joint_positivity.csv","joint_manifest.json","joint_summary.txt")]
        if not args.no_plots:artifacts.append(out_dir/"all_joint_xsec_plots.pdf")
        if release_pdf is not None:
            release_dir=release_pdf.parent
            artifacts += [release_pdf,release_dir/"calibration_summary.json",
                release_dir/"joint_preliminary_cross_sections_calibrated.csv",
                release_dir/"published_covariance_calibrated68.csv",release_dir/"pipeline_linkage.json",
                campaign/"campaign_manifest.json",campaign/"central.npz",campaign/"central.json",
                campaign/"toys/summary.json",campaign/"toys/replicas.npz"]
        records=[]
        for path in artifacts:
            need(path.is_file() and path.stat().st_size>0,f"missing/empty final artifact: {path}")
            if path.suffix==".pdf":need(path.read_bytes()[:5]==b"%PDF-",f"invalid PDF artifact: {path}")
            records.append({"path":str(path),"bytes":path.stat().st_size,"sha256":digest(path)})
        (out_dir/"joint_pipeline_artifacts.json").write_text(json.dumps({
            "schema_version":1,"verified_ns":time.time_ns(),"artifacts":records},indent=2)+"\n")
        print(f"[VERIFY] {len(records)} joint artifacts verified",flush=True)
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError, KeyError) as error:
        print(f"[FATAL] joint xsec fit: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
