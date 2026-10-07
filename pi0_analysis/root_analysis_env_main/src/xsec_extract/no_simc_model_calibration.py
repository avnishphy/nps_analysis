#!/usr/bin/env python3
"""Local calibrated release for the direct-bin no-SIMC-model estimator.

The campaign is conditional on the exported weighted-yield variances and fixed
response.  Gaussian pseudo-data include the converged Poissonized finite-MC row
variance, while every replica reruns the production continuous-positivity solve
and (when selected) the complete finite-MC variance iteration.
"""

import argparse
import csv
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import tempfile

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

PUBLICATION_TOYS = 500
SCHEMA = "no_simc_direct_bin_calibration_v1"
COMPONENTS = ("U", "LT", "TT")
DISPLAY_SCALE = 1.0e9  # ub/MeV^2 -> nb/GeV^2


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def write_rows(path, rows, fieldnames=None):
    rows = list(rows)
    need(rows or fieldnames, f"refusing schema-less empty CSV: {path}")
    names = fieldnames or list(rows[0])
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def write_matrix(path, names, matrix, unit):
    write_rows(path, ({"quantity": name, "unit": unit,
                       **{other: float(value) for other, value in zip(names, row)}}
                      for name, row in zip(names, matrix)))


def correlation(covariance):
    sd = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denominator = np.outer(sd, sd)
    return np.divide(covariance, denominator, out=np.zeros_like(covariance),
                     where=denominator > 0)


def staged_directory(destination):
    destination = Path(destination).resolve()
    need(not destination.exists(), f"refusing to overwrite existing output directory: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{destination.name}.stage.", dir=destination.parent))
    return destination, stage


class Problem:
    def __init__(self, args):
        self.fit = Path(args.fit_output).resolve()
        self.config_path = Path(args.config).resolve()
        self.central_root = Path(args.central_root).resolve()
        self.central_summary = Path(args.central_summary).resolve()
        self.central_slices = Path(args.central_slices).resolve()
        self.data_file = Path(args.data_file).resolve()
        self.sim_file = Path(args.sim_file).resolve()
        self.vertex_file = Path(args.vertex_file).resolve()
        self.bridge_path = Path(args.bridge).resolve()
        required = [self.config_path, self.central_root, self.central_summary,
                    self.central_slices, self.bridge_path, self.data_file, self.sim_file,
                    ]
        required += [self.fit / name for name in (
            "pipeline_config.json", "migration_parameters.csv", "migration_design.csv",
            "migration_response_cells.csv", "migration_reco_rows.csv",
            "migration_truth_blocks.csv", "positivity_diagnostics.csv",
            "fit_status.csv", "fit_attempts.csv")]
        for path in required:
            need(path.is_file() and path.stat().st_size > 0, f"missing/empty calibration input: {path}")
        need(self.vertex_file.exists(), f"missing calibration vertex input: {self.vertex_file}")
        self.config = json.loads(self.config_path.read_text())
        need(json.loads((self.fit / "pipeline_config.json").read_text()) == self.config,
             "selected config does not exactly match the central pipeline_config.json")
        self.rank_tolerance = float(args.rank_tolerance)
        self.max_iterations = int(args.mc_max_iterations)
        self.fit_tolerance = float(args.mc_fit_tolerance)
        self.variance_mode = args.fit_variance
        need(self.variance_mode in ("data", "finite-mc"), "invalid variance mode")
        need(args.fit_objective == "gaussian", "no-model calibration requires Gaussian objective")
        self.target_factor = float(args.target_factor if args.target_factor is not None
                                   else self.config["tgt_contam"])
        self.target_error = float(args.target_error if args.target_error is not None
                                  else self.config["tgt_contam_err"])
        need(self.target_factor > 0 and self.target_error >= 0, "invalid target divisor")

        parameters = sorted(read_rows(self.fit / "migration_parameters.csv"),
                            key=lambda row: int(row["parameter_index"]))
        need(parameters and [int(row["parameter_index"]) for row in parameters] == list(range(len(parameters))),
             "migration parameter indices are not contiguous")
        need(len(parameters) % 3 == 0, "direct estimator parameters are not U/LT/TT triplets")
        self.parameter_rows = parameters
        self.truth = np.asarray([float(row["value"]) for row in parameters], dtype=np.float64)
        need(np.isfinite(self.truth).all(), "nonfinite central parameter vector")
        self.np = len(parameters)
        self.nb = self.np // 3
        active = [int(parameters[3 * block]["active_block_index"]) for block in range(self.nb)]
        blocks = [int(parameters[3 * block]["truth_block"]) for block in range(self.nb)]
        need(active == list(range(self.nb)), "active block indices are not contiguous")
        for block in range(self.nb):
            triplet = parameters[3 * block:3 * block + 3]
            need([row["component"] for row in triplet] == list(COMPONENTS),
                 "parameter triplet is not ordered U/LT/TT")
            need(len({row["truth_block"] for row in triplet}) == 1,
                 "parameter triplet mixes truth blocks")
        self.truth_blocks = blocks

        positivity = {int(row["active_block_index"]): row
                      for row in read_rows(self.fit / "positivity_diagnostics.csv")}
        need(set(positivity) == set(range(self.nb)), "positivity diagnostics do not cover all active blocks")
        need(all(int(row["positivity_enabled"]) == 1 for row in positivity.values()),
             "central estimator was not run with --positive-xsec")
        self.epsilon = np.asarray([float(positivity[index]["epsilon_max"])
                                   for index in range(self.nb)], dtype=np.float64)
        need(np.isfinite(self.epsilon).all() and np.all((self.epsilon >= 0) & (self.epsilon <= 1)),
             "invalid event-epsilon envelope")

        reco_all = read_rows(self.fit / "migration_reco_rows.csv")
        reco = sorted((row for row in reco_all if int(row["fit_index"]) >= 0),
                      key=lambda row: int(row["fit_index"]))
        need(reco and [int(row["fit_index"]) for row in reco] == list(range(len(reco))),
             "included reconstructed fit rows are not contiguous")
        self.reco_rows = np.asarray([int(row["reco_row"]) for row in reco], dtype=np.int32)
        self.nr = len(reco)
        self.data = np.asarray([float(row["data"]) for row in reco], dtype=np.float64)
        self.data_variance = np.asarray([float(row["data_variance"]) for row in reco], dtype=np.float64)
        self.fixed = np.asarray([float(row["fixed_prediction"]) for row in reco], dtype=np.float64)
        self.fixed_mc = np.asarray([float(row["fixed_mc_variance"]) for row in reco], dtype=np.float64)
        self.central_variance = np.asarray([float(row["variance_used"]) for row in reco], dtype=np.float64)
        need(np.isfinite(self.data).all() and np.all(np.isfinite(self.data_variance)) and
             np.all(self.data_variance > 0), "invalid measured vector or observed data variance")
        need(np.isfinite(self.fixed).all() and np.all(np.isfinite(self.fixed_mc)) and
             np.all(self.fixed_mc >= 0), "invalid fixed exterior feed-in")
        need(np.isfinite(self.central_variance).all() and np.all(self.central_variance > 0),
             "invalid converged central variance")

        self.design = np.zeros((self.nr, self.np), dtype=np.float64)
        row_lookup = {int(value): index for index, value in enumerate(self.reco_rows)}
        for row in read_rows(self.fit / "migration_design.csv"):
            reco_row = int(row["reco_row"])
            if reco_row in row_lookup:
                self.design[row_lookup[reco_row], int(row["parameter_index"])] = float(row["response"])
        need(np.isfinite(self.design).all() and np.all(np.any(self.design != 0, axis=0)),
             "invalid or unsupported exported response design")

        block_lookup = {truth_block: index for index, truth_block in enumerate(self.truth_blocks)}
        self.response_covariance = np.zeros((self.nr, self.nb, 9), dtype=np.float64)
        seen = np.zeros((self.nr, self.nb), dtype=bool)
        covariance_names = [f"cov_{first}_{second}" for first in COMPONENTS for second in COMPONENTS]
        for row in read_rows(self.fit / "migration_response_cells.csv"):
            reco_row, truth_block = int(row["reco_row"]), int(row["truth_block"])
            if reco_row in row_lookup and truth_block in block_lookup:
                i, j = row_lookup[reco_row], block_lookup[truth_block]
                self.response_covariance[i, j] = [float(row[name]) for name in covariance_names]
                seen[i, j] = True
        need(seen.all() and np.isfinite(self.response_covariance).all(),
             "response-cell covariance is incomplete or nonfinite")

        published_blocks = []
        self.published_parameter_indices = []
        for block in range(self.nb):
            row = parameters[3 * block]
            if row["region"] == "published" and int(row["is_nuisance"]) == 0:
                published_blocks.append(block)
        need(published_blocks, "central fit has no retained published truth blocks")
        # Component-major release order, while the fitted vector remains block-major.
        self.published_parameter_indices = np.asarray(
            [3 * block + component for component in range(3) for block in published_blocks], dtype=int)
        self.published_blocks = published_blocks
        self.names = []
        for component in COMPONENTS:
            for block in published_blocks:
                row = parameters[3 * block]
                self.names.append(f"{component}(it={row['it']},iq={row['iq']},ix={row['ix']})")

        truth_rows = {int(row["truth_block"]): row
                      for row in read_rows(self.fit / "migration_truth_blocks.csv")}
        self.bin_rows = [truth_rows[self.truth_blocks[block]] for block in published_blocks]
        status = read_rows(self.fit / "fit_status.csv")
        self.retained_groups = [{"iq": int(row["iq"]), "ix": int(row["ix"]),
                                 "fit_ok": int(row["fit_ok"]), "fit_scope": row["fit_scope"],
                                 "failure_reason": row["failure_reason"]} for row in status]
        attempts = read_rows(self.fit / "fit_attempts.csv")
        need(attempts and int(attempts[-1]["fit_ok"]) == 1,
             "central whole-group fallback did not end in a successful fit")
        self.fit_attempts = attempts
        self.generating_mean = self.fixed + self.design @ self.truth

        self.lib = ctypes.CDLL(str(self.bridge_path))
        double_array = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
        function = self.lib.nps_no_simc_calibration_fit
        function.restype = ctypes.c_int
        function.argtypes = [ctypes.c_int, ctypes.c_int] + [double_array] * 6 + [
            ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double] + [double_array] * 6 + [
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int), ctypes.c_char_p, ctypes.c_int]
        self.function = function
        self._contiguous = [np.ascontiguousarray(value, dtype=np.float64) for value in (
            self.design.ravel(), self.data_variance, self.response_covariance.ravel(),
            self.fixed_mc, self.epsilon)]

    def solve(self, total_yield):
        y = np.ascontiguousarray(np.asarray(total_yield, dtype=np.float64) - self.fixed)
        design, data_variance, response_covariance, fixed_mc, epsilon = self._contiguous
        parameters = np.empty(self.np, dtype=np.float64)
        covariance = np.empty(self.np * self.np, dtype=np.float64)
        final_variance = np.empty(self.nr, dtype=np.float64)
        minima = np.empty(self.nb, dtype=np.float64)
        boundary_tolerance = np.empty(self.nb, dtype=np.float64)
        feasibility_tolerance = np.empty(self.nb, dtype=np.float64)
        rank, ndf = ctypes.c_int(), ctypes.c_int()
        mc_iterations, positivity_iterations, boundary = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
        condition, chi2 = ctypes.c_double(), ctypes.c_double()
        error = ctypes.create_string_buffer(2048)
        code = self.function(self.nr, self.np, design, y, data_variance,
            response_covariance, fixed_mc, epsilon,
            int(self.variance_mode == "finite-mc"), self.rank_tolerance,
            self.max_iterations, self.fit_tolerance, parameters, covariance,
            final_variance, minima, boundary_tolerance, feasibility_tolerance,
            ctypes.byref(rank), ctypes.byref(condition), ctypes.byref(chi2), ctypes.byref(ndf),
            ctypes.byref(mc_iterations), ctypes.byref(positivity_iterations),
            ctypes.byref(boundary), error, len(error))
        if code:
            raise RuntimeError(error.value.decode("utf-8", errors="replace"))
        return {"parameters": parameters, "covariance": covariance.reshape(self.np, self.np),
                "final_variance": final_variance, "minima": minima,
                "boundary_tolerance": boundary_tolerance,
                "feasibility_tolerance": feasibility_tolerance,
                "rank": rank.value, "condition": condition.value, "chi2": chi2.value,
                "ndf": ndf.value, "mc_iterations": mc_iterations.value,
                "positivity_iterations": positivity_iterations.value,
                "boundary_active": bool(boundary.value)}

    def signature(self, source_paths):
        critical = [("selected_config", self.config_path),
                    ("central_summary", self.central_summary),
                    ("central_slices", self.central_slices),
                    ("pipeline_config", self.fit / "pipeline_config.json")]
        critical += [(name, self.fit / name) for name in (
            "migration_parameters.csv", "migration_design.csv",
            "migration_response_cells.csv", "migration_reco_rows.csv",
            "migration_truth_blocks.csv", "positivity_diagnostics.csv",
            "fit_status.csv", "fit_attempts.csv")]
        sources = [Path(path).resolve() for path in source_paths]
        def identity(path):
            if path.is_file():
                stat = path.stat()
                return {"path": str(path), "kind": "file", "bytes": stat.st_size,
                        "mtime_ns": stat.st_mtime_ns}
            entries = []
            for child in sorted(path.glob("*.root")):
                stat = child.stat()
                entries.append({"name": child.name, "bytes": stat.st_size,
                                "mtime_ns": stat.st_mtime_ns})
            need(entries, f"vertex input directory has no ROOT files: {path}")
            return {"path": str(path), "kind": "directory", "root_files": entries}
        import uproot
        with uproot.open(self.central_root) as root_file:
            root_keys = sorted(root_file.keys(recursive=True, cycle=False))
        return {"schema": SCHEMA,
                "estimator": "direct-bin no-SIMC U/LT/TT constrained Gaussian forward fit",
                "configured_kinematic": self.config["configured_kinematic"],
                "fit_variance": self.variance_mode, "fit_objective": "gaussian",
                "positive_xsec": True, "rank_tolerance": self.rank_tolerance,
                "mc_max_iterations": self.max_iterations,
                "mc_fit_tolerance": self.fit_tolerance,
                "target_factor": self.target_factor, "target_error": self.target_error,
                "truth_blocks": self.truth_blocks,
                "published_parameter_indices": self.published_parameter_indices.tolist(),
                "parameter_rows": self.parameter_rows,
                "retained_groups": self.retained_groups,
                "inputs": {role: {"basename": path.name, "sha256": sha256(path)}
                           for role, path in critical},
                "upstream_sources": {"data": identity(self.data_file),
                                     "smeared_simc": identity(self.sim_file),
                                     "vertex_simc": identity(self.vertex_file)},
                "central_root_contract": {
                    "basename": self.central_root.name,
                    "key_count": len(root_keys),
                    "keys_sha256": hashlib.sha256("\n".join(root_keys).encode()).hexdigest()},
                "estimator_sources": {str(path): sha256(path) for path in sources}}


def common_problem_arguments(parser):
    parser.add_argument("--fit-output", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--central-root", type=Path, required=True)
    parser.add_argument("--central-summary", type=Path, required=True)
    parser.add_argument("--central-slices", type=Path, required=True)
    parser.add_argument("--data-file", type=Path, required=True)
    parser.add_argument("--sim-file", type=Path, required=True)
    parser.add_argument("--vertex-file", type=Path, required=True)
    parser.add_argument("--bridge", type=Path, required=True)
    parser.add_argument("--source", type=Path, action="append", default=[])
    parser.add_argument("--fit-objective", default="gaussian")
    parser.add_argument("--fit-variance", choices=("data", "finite-mc"), required=True)
    parser.add_argument("--rank-tolerance", type=float, required=True)
    parser.add_argument("--mc-max-iterations", type=int, required=True)
    parser.add_argument("--mc-fit-tolerance", type=float, required=True)
    parser.add_argument("--target-factor", type=float)
    parser.add_argument("--target-error", type=float)


def run_campaign(args):
    accepted_target = int(args.accepted)
    need(accepted_target == PUBLICATION_TOYS or args.allow_test_count,
         f"production campaigns require {PUBLICATION_TOYS} accepted toys")
    need(accepted_target > 4, "a calibration campaign needs more than four toys")
    destination, stage = staged_directory(args.output)
    try:
        problem = Problem(args)
        signature = problem.signature(args.source)
        central = problem.solve(problem.data)
        parameter_delta = np.max(np.abs(central["parameters"] - problem.truth) /
                                 np.maximum(np.abs(problem.truth), 1e-30))
        variance_delta = np.max(np.abs(central["final_variance"] - problem.central_variance) /
                                problem.central_variance)
        need(parameter_delta <= 5e-9 and variance_delta <= 5e-8,
             f"toy adapter fails central parity: parameter={parameter_delta:.3e}, variance={variance_delta:.3e}")
        parity = {"passed": True, "max_relative_parameter_difference": float(parameter_delta),
                  "max_relative_final_variance_difference": float(variance_delta),
                  "rank": central["rank"], "condition": central["condition"],
                  "chi2": central["chi2"], "ndf": central["ndf"],
                  "mc_iterations": central["mc_iterations"],
                  "boundary_active": central["boundary_active"]}
        rng = np.random.default_rng(args.seed)
        accepted, accepted_ids, variances, minima, boundary_tolerances = [], [], [], [], []
        diagnostics, attempts = [], []
        maximum_attempts = args.max_attempts or max(accepted_target * 10, accepted_target + 50)
        for toy_id in range(maximum_attempts):
            toy_y = problem.generating_mean + np.sqrt(problem.central_variance) * rng.standard_normal(problem.nr)
            try:
                result = problem.solve(toy_y)
                need(np.isfinite(result["parameters"]).all(), "nonfinite fitted vector")
                accepted.append(result["parameters"])
                accepted_ids.append(toy_id)
                variances.append(result["final_variance"])
                minima.append(result["minima"])
                boundary_tolerances.append(result["boundary_tolerance"])
                diagnostics.append((result["rank"], result["condition"], result["chi2"], result["ndf"],
                                    result["mc_iterations"], result["positivity_iterations"],
                                    int(result["boundary_active"])))
                attempts.append({"toy_id": toy_id, "status": "accepted", "reason": "",
                                 "accepted_index": len(accepted) - 1,
                                 "mc_iterations": result["mc_iterations"],
                                 "boundary_active": int(result["boundary_active"])})
            except (RuntimeError, ValueError) as error:
                attempts.append({"toy_id": toy_id, "status": "rejected", "reason": str(error),
                                 "accepted_index": -1, "mc_iterations": -1, "boundary_active": -1})
            if len(accepted) == accepted_target:
                break
        need(len(accepted) == accepted_target,
             f"only {len(accepted)} of {accepted_target} toys accepted after {maximum_attempts} attempts")
        accepted = np.asarray(accepted)
        residuals = accepted - problem.truth
        (stage / "toys").mkdir()
        with (stage / "toys" / "replicas.npz.partial").open("wb") as stream:
            np.savez_compressed(stream, index=np.asarray(accepted_ids, dtype=int), fit=accepted,
                residual=residuals, truth=problem.truth, final_variance=np.asarray(variances),
                minima=np.asarray(minima), boundary_tolerance=np.asarray(boundary_tolerances),
                diagnostics=np.asarray(diagnostics, dtype=float))
        os.replace(stage / "toys" / "replicas.npz.partial", stage / "toys" / "replicas.npz")
        write_rows(stage / "toys" / "attempts.csv", attempts)
        fit_rows = []
        for accepted_index, toy_id in enumerate(accepted_ids):
            for parameter_index, row in enumerate(problem.parameter_rows):
                fit_rows.append({"toy_id": toy_id, "accepted_index": accepted_index,
                    "parameter_index": parameter_index, "truth_block": row["truth_block"],
                    "component": row["component"], "fit": accepted[accepted_index, parameter_index],
                    "truth": problem.truth[parameter_index],
                    "residual": residuals[accepted_index, parameter_index]})
        write_rows(stage / "toys" / "fitted_vectors_and_residuals.csv", fit_rows)
        boundary_rows = []
        for accepted_index, toy_id in enumerate(accepted_ids):
            for block in range(problem.nb):
                boundary_rows.append({"toy_id": toy_id, "accepted_index": accepted_index,
                    "active_block_index": block, "truth_block": problem.truth_blocks[block],
                    "minimum_response_bracket": minima[accepted_index][block],
                    "boundary_tolerance": boundary_tolerances[accepted_index][block],
                    "boundary_hit": int(minima[accepted_index][block] <= boundary_tolerances[accepted_index][block])})
        write_rows(stage / "toys" / "boundary_hits.csv", boundary_rows)
        rejected = len(attempts) - accepted_target
        reasons = {}
        for row in attempts:
            if row["status"] == "rejected":
                reasons[row["reason"]] = reasons.get(row["reason"], 0) + 1
        summary = {"schema": SCHEMA, "requested": accepted_target, "attempted": len(attempts),
                   "accepted": accepted_target, "rejected": rejected,
                   "rejection_reasons": reasons, "seed": args.seed,
                   "rng": "numpy.random.Generator(PCG64)",
                   "publication_eligible_count": accepted_target == PUBLICATION_TOYS,
                   "central_parity": parity}
        write_json(stage / "toys" / "summary.json", summary)
        write_json(stage / "central_parity.json", parity)
        write_json(stage / "campaign_manifest.json", {"schema": SCHEMA, "signature": signature,
            "seed": args.seed, "requested_accepted_toys": accepted_target,
            "pseudoexperiment": {"data": "independent Gaussian reconstructed-row draws",
                "mean": "central fitted full forward prediction including fixed exterior feed-in",
                "variance": "central converged data plus Poissonized finite-MC row variance",
                "response": "fixed central response/migration matrix",
                "refit": "exact direct-bin continuous-positivity solve with complete selected variance iteration",
                "subset_policy": "central largest-supported whole-Q2/xB subset is frozen; fixed response and observed variances make toy rank/support invariant",
                "limitations": "no event-level upstream-weight, detector, radiative, or fixed-Ngen multinomial variation"}})
        write_json(stage / "generation_summary.json", summary)
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(f"[no-model toys] accepted={accepted_target} output={destination}")


def load_matching_campaign(args, problem):
    campaign = Path(args.campaign).resolve()
    manifest = json.loads((campaign / "campaign_manifest.json").read_text())
    need(manifest.get("schema") == SCHEMA, "incompatible no-model campaign schema")
    need(manifest.get("signature") == problem.signature(args.source),
         "toy campaign configuration/binning/estimator/source checksums do not match this central fit")
    summary = json.loads((campaign / "toys" / "summary.json").read_text())
    replicas = np.load(campaign / "toys" / "replicas.npz")
    required = int(args.accepted)
    need(required == PUBLICATION_TOYS or args.allow_test_count,
         f"production reporting requires {PUBLICATION_TOYS} accepted toys")
    need(summary["requested"] == required and summary["accepted"] == required and
         replicas["fit"].shape == (required, problem.np),
         "campaign accepted count or fitted-vector dimensions do not match")
    ids = np.asarray(replicas["index"], dtype=int)
    need(len(np.unique(ids)) == required and np.isfinite(replicas["fit"]).all(),
         "campaign toy identifiers or fitted vectors are invalid")
    need(summary["central_parity"]["passed"], "campaign did not pass central estimator parity")
    return campaign, manifest, summary, replicas


def calibrate(residuals, toy_ids, names):
    bias = residuals.mean(axis=0)
    sd = residuals.std(axis=0, ddof=1)
    rmse = np.sqrt(np.mean(residuals * residuals, axis=0))
    quantiles = np.quantile(residuals, [0.16, 0.50, 0.84], axis=0)
    delta = np.quantile(np.abs(residuals), 0.68, axis=0)
    folds = toy_ids % 5
    held = np.zeros_like(residuals, dtype=bool)
    coverage_rows = []
    for fold in range(5):
        train, validation = folds != fold, folds == fold
        need(train.sum() > 0 and validation.sum() > 0, "five-fold split has an empty fold")
        radii = np.quantile(np.abs(residuals[train]), 0.68, axis=0)
        held[validation] = np.abs(residuals[validation]) <= radii
        for index, name in enumerate(names):
            observed = float(held[validation, index].mean())
            coverage_rows.append({"quantity": name, "scope": f"fold_{fold}",
                "training_n": int(train.sum()), "validation_n": int(validation.sum()),
                "radius_learned_elsewhere": float(radii[index]), "coverage": observed,
                "binomial_se": math.sqrt(observed * (1 - observed) / validation.sum())})
    coverage = held.mean(axis=0)
    aggregate = []
    for index, name in enumerate(names):
        observed = float(coverage[index])
        aggregate.append({"quantity": name, "scope": "all_held_out", "training_n": -1,
            "validation_n": len(toy_ids), "radius_learned_elsewhere": float("nan"),
            "coverage": observed,
            "binomial_se": math.sqrt(observed * (1 - observed) / len(toy_ids))})
    return {"bias": bias, "sd": sd, "rmse": rmse, "quantiles": quantiles,
            "delta": delta, "coverage": coverage,
            "coverage_rows": aggregate + coverage_rows}


def run_check(args):
    problem = Problem(args)
    load_matching_campaign(args, problem)
    print("[no-model calibration] exact campaign match")


def run_report(args):
    destination, stage = staged_directory(args.output)
    try:
        problem = Problem(args)
        campaign, manifest, summary, replicas = load_matching_campaign(args, problem)
        ids = np.asarray(replicas["index"], dtype=int)
        indices = problem.published_parameter_indices
        truth_raw = problem.truth[indices]
        fitted_raw = np.asarray(replicas["fit"], dtype=float)[:, indices]
        residuals = (fitted_raw - truth_raw) * DISPLAY_SCALE
        central = truth_raw * DISPLAY_SCALE
        stats = calibrate(residuals, ids, problem.names)
        covariance = np.cov(residuals, rowvar=False, ddof=1)
        corr = correlation(covariance)
        covariance68 = np.outer(stats["delta"], stats["delta"]) * corr
        target_fraction = problem.target_error / problem.target_factor
        target_covariance = np.outer(central, central) * target_fraction ** 2
        diagnostics = np.asarray(replicas["diagnostics"])
        minima = np.asarray(replicas["minima"])
        tolerances = np.asarray(replicas["boundary_tolerance"])
        block_hit_rate = np.mean(minima <= tolerances, axis=0)
        rejected_fraction = summary["rejected"] / summary["attempted"]

        details = []
        for index, name in enumerate(problem.names):
            bias_ratio = abs(stats["bias"][index]) / stats["delta"][index] if stats["delta"][index] > 0 else math.inf
            passed = (np.isfinite(stats["delta"][index]) and stats["delta"][index] > 0 and
                      stats["coverage"][index] >= 0.55 and bias_ratio < 1.0)
            parameter_index = int(indices[index])
            block = parameter_index // 3
            details.append({"quantity": name, "parameter_index": parameter_index,
                "truth_block": problem.truth_blocks[block], "component": COMPONENTS[parameter_index % 3],
                "central_value": central[index], "toy_mean_bias": stats["bias"][index],
                "toy_sd": stats["sd"][index], "toy_rmse": stats["rmse"][index],
                "residual_q16": stats["quantiles"][0, index],
                "residual_q50": stats["quantiles"][1, index],
                "residual_q84": stats["quantiles"][2, index],
                "calibrated_delta68": stats["delta"][index],
                "basic_interval_low": central[index] - stats["quantiles"][2, index],
                "basic_interval_high": central[index] - stats["quantiles"][0, index],
                "cross_validated_coverage": stats["coverage"][index],
                "bias_over_delta68": bias_ratio, "boundary_hit_rate": block_hit_rate[block],
                "quantity_gate_passed": int(passed), "unit": "nb/GeV^2"})
        global_gate = (all(row["quantity_gate_passed"] for row in details) and
                       summary["accepted"] == int(args.accepted) and rejected_fraction <= 0.10 and
                       summary["central_parity"]["passed"])
        need(global_gate, "no-model calibrated release gates failed; campaign retained for audit")

        write_rows(stage / "toy_residual_statistics.csv", details)
        write_rows(stage / "cross_validated_coverage.csv", stats["coverage_rows"])
        write_matrix(stage / "published_covariance_toy.csv", problem.names, covariance,
                     "(nb/GeV^2)^2 literal empirical covariance")
        write_matrix(stage / "published_correlation_toy.csv", problem.names, corr, "dimensionless")
        write_matrix(stage / "published_covariance_calibrated68.csv", problem.names, covariance68,
                     "(nb/GeV^2)^2 calibrated-radius representation, not literal variance")
        write_matrix(stage / "target_systematic_covariance.csv", problem.names, target_covariance,
                     "(nb/GeV^2)^2 separate correlated target scale")

        bin_statistics, table = [], []
        nbin = len(problem.published_blocks)
        by_parameter = {int(row["parameter_index"]): row for row in details}
        edges = problem.config["tprime_bin_edges"]
        for local, (block, truth_row) in enumerate(zip(problem.published_blocks, problem.bin_rows)):
            row = problem.parameter_rows[3 * block]
            bin_statistics.append({"truth_block": problem.truth_blocks[block], "it": row["it"],
                "iq": row["iq"], "ix": row["ix"], "tprime_low": edges[int(row["it"])],
                "tprime_high": edges[int(row["it"]) + 1],
                "response_weight": truth_row["response_weight"], "events": truth_row["events"],
                "mean_q2": truth_row["mean_q2"], "mean_xb": truth_row["mean_xb"],
                "mean_tprime": truth_row["mean_tprime"], "mean_epsilon": truth_row["mean_epsilon"],
                "epsilon_max": truth_row["epsilon_max"],
                "boundary_hit_rate": block_hit_rate[block]})
            output = {"truth_block": problem.truth_blocks[block], "it": row["it"],
                      "iq": row["iq"], "ix": row["ix"],
                      "tprime_low": edges[int(row["it"])],
                      "tprime_high": edges[int(row["it"]) + 1],
                      "mean_tprime": truth_row["mean_tprime"], "unit": "nb/GeV^2",
                      "target_relative_scale_uncertainty": target_fraction,
                      "status": "TEST-ONLY" if int(args.accepted) != PUBLICATION_TOYS else "CALIBRATED RELEASE"}
            for component in range(3):
                detail = by_parameter[3 * block + component]
                label = COMPONENTS[component]
                output[f"sigma_{label}"] = detail["central_value"]
                output[f"{label}_calibrated_delta68"] = detail["calibrated_delta68"]
                output[f"{label}_basic_low"] = detail["basic_interval_low"]
                output[f"{label}_basic_high"] = detail["basic_interval_high"]
                output[f"{label}_target_scale_error"] = abs(detail["central_value"]) * target_fraction
            table.append(output)
        write_rows(stage / "bin_statistics.csv", bin_statistics)
        write_rows(stage / "preliminary_cross_sections_calibrated.csv", table)

        for component_index, component in enumerate(COMPONENTS):
            component_details = details[component_index * nbin:(component_index + 1) * nbin]
            x = np.asarray([-float(problem.bin_rows[i]["mean_tprime"]) for i in range(nbin)])
            y = np.asarray([row["central_value"] for row in component_details])
            error = np.asarray([row["calibrated_delta68"] for row in component_details])
            figure, axis = plt.subplots(figsize=(6.8, 4.9))
            axis.errorbar(x, y, yerr=error, fmt="o", capsize=5, color="#174a73")
            axis.axhline(0, color="0.55", linewidth=0.8)
            axis.set(xlabel=r"$-t'\ [\mathrm{GeV}^2]$",
                     ylabel=rf"$\sigma_{{{component}}}\ [\mathrm{{nb}}/\mathrm{{GeV}}^2]$")
            axis.set_title(f"{problem.config['configured_kinematic']} direct-bin no-SIMC estimator\n"
                           "local toy-calibrated 68% radii; target scale separate")
            figure.tight_layout()
            figure.savefig(stage / f"sigma_{component}_no_simc_model_calibrated.pdf")
            figure.savefig(stage / f"sigma_{component}_no_simc_model_calibrated.png", dpi=180)
            plt.close(figure)

        with PdfPages(stage / "fit_coverage_diagnostics.pdf") as pdf:
            figure, axis = plt.subplots(figsize=(10, 7)); axis.axis("off")
            axis.text(.03, .97, f"Direct-bin no-SIMC calibrated release\n\n"
                f"kinematic: {problem.config['configured_kinematic']}\n"
                f"requested/attempted/accepted/rejected: {summary['requested']}/{summary['attempted']}/"
                f"{summary['accepted']}/{summary['rejected']}\nseed: {summary['seed']}\n"
                f"fit variance: {problem.variance_mode}\nparameters: {problem.np}; published: {len(problem.names)}\n"
                f"minimum held-out coverage: {min(stats['coverage']):.3f}\n"
                f"maximum |bias|/delta68: {max(row['bias_over_delta68'] for row in details):.3f}\n"
                f"rejected fraction: {rejected_fraction:.4f}\n\nGLOBAL GATE: PASS",
                va="top", family="monospace")
            pdf.savefig(figure); plt.close(figure)
            figure, axis = plt.subplots(figsize=(11, 5))
            axis.plot(range(len(problem.names)), stats["coverage"], "o")
            axis.axhline(.68, color="0.25"); axis.axhline(.55, color="firebrick", linestyle="--")
            axis.set_xticks(range(len(problem.names))); axis.set_xticklabels(problem.names, rotation=90)
            axis.set_ylabel("five-fold held-out coverage"); figure.tight_layout(); pdf.savefig(figure); plt.close(figure)
            figure, axis = plt.subplots(figsize=(11, 5))
            axis.plot(range(problem.nb), block_hit_rate, "o")
            axis.set(xlabel="active truth block index", ylabel="boundary-hit rate", ylim=(-.02, 1.02))
            figure.tight_layout(); pdf.savefig(figure); plt.close(figure)

        linkage = {"schema_version": 1,
            "calibrated_estimator": "direct-bin no-SIMC U/LT/TT global constrained forward-response estimator",
            "central_csv": str(problem.central_summary), "central_slice_csv": str(problem.central_slices),
            "central_root": str(problem.central_root),
            "intervals_apply_to": "retained published independent U/LT/TT coefficients listed in migration_parameters and the matching central summary/slice fields",
            "intervals_do_not_apply_to": ["SIMC-model/SigParam or additive-M0 outputs",
                "per-phi experimental-point transformations", "the legacy 256-toy positivity diagnostic"],
            "target_normalization_uncertainty": "separate rank-one scale covariance; not included in calibrated radii",
            "fixed_exterior_feedin": "central fixed prediction retained in every toy; its Poissonized MC variance enters finite-MC row noise and iteration",
            "fitted_exterior_feedin": "low-tprime exterior U/LT/TT triplet is refitted and constrained in every toy",
            "finite_mc_assumption": "fixed response with Gaussian row fluctuation from exported Poissonized event outer products; fixed-Ngen correlations not modeled",
            "campaign": str(campaign)}
        write_json(stage / "pipeline_linkage.json", linkage)
        input_manifest = {"schema": SCHEMA, "campaign_manifest": manifest,
            "current_central_root": str(problem.central_root),
            "current_central_root_sha256": sha256(problem.central_root),
            "campaign_manifest_sha256": sha256(campaign / "campaign_manifest.json"),
            "replicas_sha256": sha256(campaign / "toys" / "replicas.npz"),
            "attempts_sha256": sha256(campaign / "toys" / "attempts.csv")}
        write_json(stage / "input_manifest.json", input_manifest)
        calibration_summary = {"schema": SCHEMA, "status": "PASS",
            "production_release": int(args.accepted) == PUBLICATION_TOYS,
            "requested": summary["requested"], "attempted": summary["attempted"],
            "accepted": summary["accepted"], "rejected": summary["rejected"],
            "rejection_reasons": summary["rejection_reasons"], "seed": summary["seed"],
            "estimator_checksum": hashlib.sha256(json.dumps(manifest["signature"]["estimator_sources"],
                sort_keys=True).encode()).hexdigest(),
            "configuration_checksum": sha256(problem.config_path),
            "truth_vector_raw_ub_per_MeV2": problem.truth.tolist(),
            "published_quantity_order": problem.names,
            "minimum_cross_validated_coverage": float(np.min(stats["coverage"])),
            "maximum_bias_over_delta68": float(max(row["bias_over_delta68"] for row in details)),
            "rejected_fraction": rejected_fraction,
            "global_gate": "PASS", "per_quantity_gates": {row["quantity"]: "PASS" for row in details},
            "covariance_interpretation": {"published_covariance_toy.csv": "literal empirical residual covariance",
                "published_covariance_calibrated68.csv": "diag(delta68) R_toy diag(delta68), not literal variances",
                "target_systematic_covariance.csv": "separate fully correlated target-divisor scale"}}
        write_json(stage / "calibration_summary.json", calibration_summary)
        report = f"""# {problem.config['configured_kinematic']} direct-bin no-SIMC calibrated release

**{'PRODUCTION CALIBRATED RELEASE' if int(args.accepted) == PUBLICATION_TOYS else 'TEST-ONLY REDUCED-TOY RELEASE'}: PASS.**

This release calibrates the independent per-truth-bin U, LT, and TT coefficients from the no-SIMC-model global forward-response estimator. It is not a SigParam, shape-slope, or additive-M0 extraction. The central estimator uses SIMC only for the absolute response/migration matrix and normalization, matched generated coordinates, vertex epsilon, Poissonized finite-MC moments, and nominal fixed contributions from exterior faces other than the fitted low-tprime face.

The campaign requested {summary['requested']} accepted toys, attempted {summary['attempted']}, accepted {summary['accepted']}, and rejected {summary['rejected']} with deterministic seed {summary['seed']}. Every accepted toy reran continuous angular positivity and the complete `{problem.variance_mode}` variance path. Failed fits were retained in the campaign attempt log and were never replaced by the central result. Central adapter parity passed with maximum relative parameter difference {summary['central_parity']['max_relative_parameter_difference']:.3e}.

The local pseudoexperiment is conditional: reconstructed weighted yields are Gaussian around the fitted forward mean with the central converged row variance. The response matrix is fixed. Finite-MC uncertainty enters through the exported event-outer-product row variance both in generation and in every refit's feasible-GLS iteration. This is an approximation to response statistics, not an event-level response bootstrap; fixed-Ngen multinomial correlations, upstream signal-weight refits, detector/radiative effects, and other systematics are outside these intervals.

The central extractor's deterministic largest-supported whole-Q2/xB fallback is recorded and frozen. Because this local ensemble fixes the response, observed row variances, and support, its rank/support decision is toy-invariant; each toy still performs the full-rank check and a singular or nonconvergent solve is rejected and recorded.

The reported `delta68` is the 0.68 quantile of absolute fit-minus-truth residuals. Five-fold coverage is held out by `toy_id modulo 5`. Minimum held-out coverage is {np.min(stats['coverage']):.3f}; maximum |bias|/delta68 is {max(row['bias_over_delta68'] for row in details):.3f}. Central values are not clipped or bias corrected. The calibrated covariance representation has diagonal `delta68^2` but is explicitly not asserted to be a literal variance.

Target-divisor uncertainty remains a separate correlated scale covariance. The intervals apply only to the retained published direct-bin U/LT/TT coefficients named in `pipeline_linkage.json`; they do not apply to per-phi experimental points, SIMC-model products, or the legacy 256-toy plot-spread diagnostic.
"""
        (stage / "REPORT.md").write_text(report)
        artifact_paths = sorted(path for path in stage.iterdir() if path.is_file() and path.name != "artifact_manifest.json")
        write_json(stage / "artifact_manifest.json", {"schema": 1,
            "artifacts": [{"name": path.name, "bytes": path.stat().st_size,
                           "sha256": sha256(path)} for path in artifact_paths]})
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(f"[no-model release] output={destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    campaign = sub.add_parser("campaign")
    common_problem_arguments(campaign)
    campaign.add_argument("--output", type=Path, required=True)
    campaign.add_argument("--accepted", type=int, default=PUBLICATION_TOYS)
    campaign.add_argument("--seed", type=int, default=20261008)
    campaign.add_argument("--max-attempts", type=int)
    campaign.add_argument("--allow-test-count", action="store_true")
    check = sub.add_parser("check")
    common_problem_arguments(check)
    check.add_argument("--campaign", type=Path, required=True)
    check.add_argument("--accepted", type=int, default=PUBLICATION_TOYS)
    check.add_argument("--allow-test-count", action="store_true")
    report = sub.add_parser("report")
    common_problem_arguments(report)
    report.add_argument("--campaign", type=Path, required=True)
    report.add_argument("--output", type=Path, required=True)
    report.add_argument("--accepted", type=int, default=PUBLICATION_TOYS)
    report.add_argument("--allow-test-count", action="store_true")
    args = parser.parse_args()
    if args.command == "campaign":
        run_campaign(args)
    elif args.command == "check":
        run_check(args)
    else:
        run_report(args)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, KeyError) as error:
        raise SystemExit(f"[ERROR] {error}")
