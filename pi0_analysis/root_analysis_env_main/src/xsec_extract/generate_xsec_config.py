#!/usr/bin/env python3
"""Render the C++ xsec configuration header from a validated JSON preset."""

import argparse
import json
import math
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / "xsec_config_template.h.in"

STRINGS = {
    "configured_kinematic", "simc_file", "data_file", "simc_tree", "data_tree",
    "out_root", "out_csv", "out_slice_csv", "out_dir", "out_all_plots_pdf",
    "mmiss_select", "fit_variance_mode", "model_xsec_mode", "model_identifier",
    "model_free_parameters",
}
NUMBERS = {
    "ebeam", "hms_theta_deg", "tgt_contam", "tgt_contam_err", "mmiss_lower_gev",
    "mmiss_upper_gev", "rank_tolerance", "mc_fit_tolerance",
    "simc_yield_scale", "model_tolerance",
}
INTEGERS = {
    "mc_max_iterations", "partons_warmups", "partons_calls",
    "model_max_iterations", "model_max_evaluations",
}
BOOLEANS = {"positive_xsec", "normalize_mmiss"}
VECTORS = {"tprime_bin_edges", "t_bin_edges", "q2_bin_edges"}
MATRICES = {"xb_bin_edges_by_q2"}
STRING_VECTORS = {"model_xsec_candidates"}
OPTIONAL = {"phi_bins", "phi_bin_edges", "diamond_xb_q2_vertices", "hms_p_gev"}
FIELDS = STRINGS | NUMBERS | INTEGERS | BOOLEANS | VECTORS | MATRICES | STRING_VECTORS


def numbers(value, name):
    if not isinstance(value, list) or len(value) < 2:
        raise ValueError(f"{name} needs at least two edges")
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
           for x in value):
        raise ValueError(f"{name} edges must be finite numbers")
    if any(a >= b for a, b in zip(value, value[1:])):
        raise ValueError(f"{name} edges must be strictly increasing")
    return value


def cpp_number(value):
    return repr(float(value))


def cpp_vector(value):
    return "{" + ", ".join(map(cpp_number, value)) + "}"


def diamond_vertices(value):
    if value is None:
        return []
    if not isinstance(value, list) or len(value) != 4:
        raise ValueError("diamond_xb_q2_vertices needs four [xB, Q2] corners or null")
    vertices = []
    for point in value:
        if (not isinstance(point, list) or len(point) != 2 or
                any(isinstance(x, bool) or not isinstance(x, (int, float)) or
                    not math.isfinite(x) for x in point)):
            raise ValueError("diamond corners must be finite [xB, Q2] number pairs")
        xb, q2 = point
        if not (0 < xb < 1 and q2 > 0):
            raise ValueError("diamond corners require 0 < xB < 1 and Q2 > 0")
        vertices.append((float(xb), float(q2)))
    center_x = sum(p[0] for p in vertices) / 4
    center_q = sum(p[1] for p in vertices) / 4
    vertices.sort(key=lambda p: math.atan2(p[1] - center_q, p[0] - center_x))
    for i in range(4):
        a, b, c = (vertices[(i + j) % 4] for j in range(3))
        cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0])
        if cross <= 1e-12:
            raise ValueError("diamond needs four distinct convex corners")
    return vertices


def render(config):
    missing = FIELDS - config.keys()
    extra = config.keys() - FIELDS - OPTIONAL
    if missing or extra:
        raise ValueError(f"Missing keys: {sorted(missing)}; unknown keys: {sorted(extra)}")
    if ("phi_bins" in config) == ("phi_bin_edges" in config):
        raise ValueError("Specify exactly one of phi_bins or phi_bin_edges")
    for name in STRINGS:
        if not isinstance(config[name], str):
            raise ValueError(f"{name} must be a string")
    for name in NUMBERS:
        x = config[name]
        if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
            raise ValueError(f"{name} must be a finite number")
    for name in INTEGERS:
        if isinstance(config[name], bool) or not isinstance(config[name], int):
            raise ValueError(f"{name} must be an integer")
    for name in BOOLEANS:
        if not isinstance(config[name], bool):
            raise ValueError(f"{name} must be true or false")
    for name in VECTORS:
        numbers(config[name], name)
    q2_bins = len(config["q2_bin_edges"]) - 1
    rows = config["xb_bin_edges_by_q2"]
    if not isinstance(rows, list) or len(rows) != q2_bins:
        raise ValueError("xb_bin_edges_by_q2 needs one row per Q2 bin")
    for row in rows:
        numbers(row, "xb_bin_edges_by_q2")
    if any(len(row) != len(rows[0]) or row[0] != rows[0][0] or row[-1] != rows[0][-1]
           for row in rows):
        raise ValueError("xB rows need matching counts and outer edges")
    if not isinstance(config["model_xsec_candidates"], list) or not all(
        isinstance(x, str) for x in config["model_xsec_candidates"]
    ):
        raise ValueError("model_xsec_candidates must be a list of strings")
    if config["mmiss_lower_gev"] >= config["mmiss_upper_gev"]:
        raise ValueError("mmiss_lower_gev must be below mmiss_upper_gev")
    if config["tgt_contam"] <= 0 or config["tgt_contam_err"] < 0:
        raise ValueError("target contamination factor must be positive; error nonnegative")
    if config["ebeam"] <= 0:
        raise ValueError("ebeam must be positive")
    if not 0 < config["hms_theta_deg"] < 180:
        raise ValueError("hms_theta_deg must be an electron-arm angle in (0, 180) degrees")
    # Used only by post-fit L/T separation, not event vertex epsilon.
    if "hms_p_gev" in config:
        p = config["hms_p_gev"]
        if (isinstance(p, bool) or not isinstance(p, (int, float)) or
                not math.isfinite(p) or not 0 < p < config["ebeam"]):
            raise ValueError("hms_p_gev must be finite and satisfy 0 < p < ebeam")

    values = {}
    for name in STRINGS:
        values[name] = json.dumps(config[name])
    for name in NUMBERS:
        values[name] = cpp_number(config[name])
    for name in INTEGERS:
        values[name] = str(config[name])
    for name in BOOLEANS:
        values[name] = "true" if config[name] else "false"
    for name in VECTORS:
        values[name] = cpp_vector(config[name])
    values["xb_bin_edges_by_q2"] = "{" + ", ".join(cpp_vector(row) for row in rows) + "}"
    values["diamond_xb_q2_vertices"] = "{" + ", ".join(
        "{" + ", ".join(cpp_number(x) for x in point) + "}"
        for point in diamond_vertices(config.get("diamond_xb_q2_vertices"))
    ) + "}"
    values["model_xsec_candidates"] = (
        "{" + ", ".join(json.dumps(x) for x in config["model_xsec_candidates"]) + "}"
    )
    if "phi_bins" in config:
        n = config["phi_bins"]
        if isinstance(n, bool) or not isinstance(n, int) or n < 3:
            raise ValueError("phi_bins must be an integer >= 3")
        terms = ["0.0"] + [f"{2*i}*TMath::Pi()/{n}" for i in range(1, n)]
        terms += ["2*TMath::Pi()"]
        values["phi_bin_edges"] = "{" + ", ".join(terms) + "}"
    else:
        edges = numbers(config["phi_bin_edges"], "phi_bin_edges")
        values["phi_bin_edges"] = cpp_vector(edges)

    template = TEMPLATE.read_text()
    tokens = set(re.findall(r"@([a-z][a-z0-9_]*)@", template))
    if tokens != set(values):
        raise ValueError(f"Template mismatch: missing {sorted(set(values)-tokens)}, "
                         f"unfilled {sorted(tokens-set(values))}")
    for name, value in values.items():
        template = template.replace(f"@{name}@", value)
    return template


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="xsec_config_*.json preset")
    parser.add_argument("output", type=Path, help="generated xsec_config.h")
    args = parser.parse_args()
    def reject_constant(value):
        raise ValueError(f"Nonfinite JSON number: {value}")
    try:
        config = json.loads(args.config.read_text(), parse_constant=reject_constant)
        if not isinstance(config, dict):
            raise ValueError("Config root must be a JSON object")
        args.output.write_text(render(config))
    except (OSError, ValueError) as error:
        print(f"[ERROR] Invalid xsec config {args.config}: {error}", file=sys.stderr)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
