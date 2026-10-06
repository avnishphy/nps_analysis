"""Build an event-level, harmonic forward problem with independent domains.

This is extraction infrastructure, not certification of a physics result.
Input ``data_events.csv`` columns: event_id,run_number,q2,xb,tprime,phi,weight.
Input ``mc_events.csv`` columns: event_id,reco_q2,reco_xb,reco_tprime,reco_phi,
truth_q2,truth_xb,truth_tprime,truth_phi,epsilon,base_weight,nominal_weight. Angles are radians;
tprime is the signed t-t_min in GeV^2. Data weight is the corrected yield/mC,
including the target correction. MC base_weight is full_weight/sigcm, with
the original generator normalization; nominal_weight is full_weight. The common selection in the cache must
precede rectangular fit cuts. This module never renormalizes MC columns.

Configuration JSON (or equivalent dict)::

    {"reco": {"q2_edges": [3.3, 4.7], "xb_edges_by_q2": [[0.29, 0.44]],
              "tprime_edges": [-0.75, -0.4, 0],
              "phi_edges": [0, 1.5707963267948966, 3.141592653589793,
                            4.71238898038469, 6.283185307179586]},
     "truth": {"q2_edges": [3.3, 4.7], "xb_edges_by_q2": [[0.29, 0.44]],
               "tprime_edges": [-0.75, -0.4, 0]},
     "publication_truth_blocks": [1]}

Interior IDs iterate tprime, then Q2, then xB. Variable xB bin counts per Q2
are supported. Six exterior IDs follow all interior IDs, with disjoint corner
priority tprime, Q2, xB. Supported interior blocks and tprime_below retain
independent U/LT/TT coefficients. Q2/xB and all other exterior events enter as
a fixed nominal reconstructed-row contribution with no parameter columns.
Only publication selection changes which interior coefficients are reported.
Edges use [low, high), including the final high edge except periodic phi.
Reco phi edges partition a full period with any origin; phi wraps into
[phi_edges[0], phi_edges[0]+2*pi). Empty count rows are preserved.

The output coefficient unit is nb/GeV^2. Its event basis is
base_weight * 1e-9/(2*pi) * [1,sqrt(2*epsilon*(1+epsilon))*cos(phi),
epsilon*cos(2*phi)]. A fitted binwise constant is response-weighted; it is
not automatically a phase-space average or bin-center cross section.
"""

from pathlib import Path
import csv
import json

import numpy as np


GUARD_NAMES = ("tprime_below", "tprime_above", "q2_below", "q2_above",
               "xb_below", "xb_above")
HARMONICS = ("U", "LT", "TT")
COEFFICIENT_UNIT = "nb/GeV^2"


def _edges(values, name):
    try:
        out = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name}: expected finite numerical bin edges") from error
    if out.ndim != 1 or len(out) < 2 or not np.all(np.isfinite(out)):
        raise ValueError(f"{name}: need at least two finite edges")
    if not np.all(np.diff(out) > 0):
        raise ValueError(f"{name}: edges must be strictly increasing")
    return out


def _domain(config, name):
    if name not in config or not isinstance(config[name], dict):
        raise ValueError(f"Missing {name} domain configuration")
    section = config[name]
    q = _edges(section.get("q2_edges"), f"{name}.q2_edges")
    t = _edges(section.get("tprime_edges"), f"{name}.tprime_edges")
    if q[0] <= 0:
        raise ValueError(f"{name}: Q2 edges must be positive")
    xrows = section.get("xb_edges_by_q2")
    if not isinstance(xrows, list) or len(xrows) != len(q) - 1:
        raise ValueError(f"{name}: one xB edge vector is required per Q2 bin")
    x = [_edges(row, f"{name}.xb_edges_by_q2[{i}]") for i, row in enumerate(xrows)]
    if any(row[0] <= 0 or row[-1] >= 1 for row in x):
        raise ValueError(f"{name}: xB edges must lie inside (0, 1)")
    offsets = np.cumsum([0] + [len(row) - 1 for row in x])
    result = dict(q=q, x=x, t=t, offsets=offsets,
                  nblocks=(len(t) - 1) * int(offsets[-1]))
    if name == "reco":
        p = _edges(section.get("phi_edges"), "reco.phi_edges")
        if not np.isclose(p[-1] - p[0], 2*np.pi, atol=1e-12, rtol=0):
            raise ValueError("reco.phi_edges must span one full 2*pi period in radians")
        p[-1] = p[0] + 2*np.pi
        result["p"] = p
    return result


def _read_csv(path, fields, text_fields):
    columns = {field: [] for field in fields}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = set(fields) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for line, row in enumerate(reader, 2):
            for field in fields:
                value = row.get(field)
                try:
                    if value is None or not value.strip():
                        raise ValueError("missing value")
                    columns[field].append(value.strip() if field in text_fields else float(value))
                except (TypeError, ValueError) as error:
                    raise ValueError(f"{path}:{line}: invalid {field}") from error
    for field in fields:
        columns[field] = np.asarray(columns[field], dtype=str if field in text_fields else float)
        if field not in text_fields and not np.all(np.isfinite(columns[field])):
            raise ValueError(f"{path}: nonfinite {field}")
    return columns


def _bin(edges, values):
    index = np.searchsorted(edges, values, side="right") - 1
    index[values == edges[-1]] = len(edges) - 2
    return index


def _assign_blocks(domain, q, x, t, exterior=False):
    """Vectorized disjoint block allocation; -1 denotes outside reco domain."""
    result = np.full(q.shape, -1, dtype=np.int64)
    nq, nx_per_t = len(domain["q"]) - 1, int(domain["offsets"][-1])
    iq, it = _bin(domain["q"], q), _bin(domain["t"], t)
    for qi in range(nq):
        xe = domain["x"][qi]
        inside = ((iq == qi) & (it >= 0) & (it < len(domain["t"]) - 1)
                  & (x >= xe[0]) & (x <= xe[-1]))
        result[inside] = (it[inside] * nx_per_t + domain["offsets"][qi]
                          + _bin(xe, x[inside]))
    if exterior:
        count = domain["nblocks"]
        conditions = [t < domain["t"][0], t > domain["t"][-1],
                      q < domain["q"][0], q > domain["q"][-1]]
        for face, mask in enumerate(conditions):
            result[(result < 0) & mask] = count + face
        for qi in range(nq):
            xe = domain["x"][qi]
            result[(result < 0) & (iq == qi) & (x < xe[0])] = count + 4
            result[(result < 0) & (iq == qi) & (x > xe[-1])] = count + 5
        if np.any(result < 0):
            raise ValueError("Unassigned truth coordinates")
    return result


def _block_metadata(domain):
    result = []
    for ti in range(len(domain["t"]) - 1):
        for qi, xe in enumerate(domain["x"]):
            for xi in range(len(xe) - 1):
                result.append(dict(global_id=len(result), kind="interior", tprime_index=ti,
                                   q2_index=qi, xb_index=xi,
                                   tprime_bounds=domain["t"][ti:ti+2].tolist(),
                                   q2_bounds=domain["q"][qi:qi+2].tolist(),
                                   xb_bounds=xe[xi:xi+2].tolist()))
    return result


def _row_indices(domain, data, prefix=""):
    b = _assign_blocks(domain, data[prefix + "q2"], data[prefix + "xb"],
                       data[prefix + "tprime"])
    phi = domain["p"][0] + np.mod(data[prefix + "phi"] - domain["p"][0], 2*np.pi)
    rows = b * (len(domain["p"]) - 1) + _bin(domain["p"], phi)
    rows[b < 0] = -1
    return rows


def build_problem(cache_dir, config):
    """Read an immutable cache and return arrays needed for fit/resampling.

    ``data_ids`` identify independent events by (run_number,event_id).
    ``mc_ids`` identify independent generator events, allowing repeated cache
    rows with the same ID to share a resampling multiplier. Dataset namespaces
    are separate. ``bootstrap_data_ids`` and ``bootstrap_mc_ids`` contain the
    unique sorted IDs from the full cache, so identical seeds can pair event
    resampling even when reconstructed domains differ. All reconstructed rows
    remain, including empty data rows.
    """
    if isinstance(config, (str, Path)):
        with Path(config).open(encoding="utf-8") as handle:
            config = json.load(handle)
    if not isinstance(config, dict):
        raise ValueError("Configuration must be a JSON object")
    reco, truth = _domain(config, "reco"), _domain(config, "truth")
    publication = config.get("publication_truth_blocks")
    if (not isinstance(publication, list) or not publication
            or any(type(value) is not int for value in publication)
            or len(set(publication)) != len(publication)
            or any(value < 0 or value >= truth["nblocks"] for value in publication)):
        raise ValueError("publication_truth_blocks must list unique, valid interior integer IDs")
    cache_dir = Path(cache_dir)
    manifest_path = cache_dir / "forward_cache_manifest.json"
    manifest = None
    if manifest_path.exists():
        with manifest_path.open(encoding="utf-8") as handle:
            manifest = json.load(handle)
        if manifest.get("complete") is not True:
            raise ValueError("Event cache manifest is incomplete; export must finish before fitting")
        if manifest.get("schema_version") != 2:
            raise ValueError("Regenerate forward cache with nominal_weight support for fixed Q2/xB feed-in")
    data = _read_csv(cache_dir / "data_events.csv",
                     ["event_id", "run_number", "q2", "xb", "tprime", "phi", "weight"],
                     {"event_id", "run_number"})
    mc = _read_csv(cache_dir / "mc_events.csv",
                   ["event_id", "reco_q2", "reco_xb", "reco_tprime", "reco_phi",
                    "truth_q2", "truth_xb", "truth_tprime", "truth_phi", "epsilon", "base_weight", "nominal_weight"],
                   {"event_id"})
    if np.any(data["weight"] < 0) or np.any(mc["base_weight"] < 0) or np.any(mc["nominal_weight"] < 0):
        raise ValueError("This positive-weight likelihood requires nonnegative data and MC weights")
    if np.any((mc["epsilon"] < 0) | (mc["epsilon"] > 1)):
        raise ValueError("Vertex epsilon must lie in [0, 1]")
    data_row_all, mc_row_all = _row_indices(reco, data), _row_indices(reco, mc, "reco_")
    data_keep, mc_keep = data_row_all >= 0, mc_row_all >= 0
    mc_global_all = _assign_blocks(truth, mc["truth_q2"], mc["truth_xb"],
                                    mc["truth_tprime"], exterior=True)
    supported_global = np.unique(mc_global_all[mc_keep & (mc["base_weight"] > 0)])
    unsupported = sorted(set(publication) - set(supported_global.tolist()))
    if unsupported:
        raise ValueError(f"Published truth blocks have no positive MC response: {unsupported}")
    fitted_global = np.asarray([gid for gid in supported_global
        if gid < truth["nblocks"] or GUARD_NAMES[int(gid)-truth["nblocks"]] == "tprime_below"], dtype=int)
    fitted_keep = mc_keep & (mc["base_weight"] > 0) & np.isin(mc_global_all, fitted_global)
    fixed_keep = mc_keep & (mc["nominal_weight"] > 0) & ~np.isin(mc_global_all, fitted_global)
    rows, weights = data_row_all[data_keep], data["weight"][data_keep]
    mc_rows = mc_row_all[fitted_keep]
    mc_blocks = np.searchsorted(fitted_global, mc_global_all[fitted_keep])
    eps = mc["epsilon"][fitted_keep]
    phi = np.mod(mc["truth_phi"][fitted_keep], 2*np.pi)
    base = mc["base_weight"][fitted_keep] * (1e-9 / (2 * np.pi))
    basis = base[:, None] * np.column_stack((np.ones(len(base)),
             np.sqrt(2 * eps * (1 + eps)) * np.cos(phi), eps * np.cos(2 * phi)))
    nrow, nblock = reco["nblocks"] * (len(reco["p"]) - 1), len(fitted_global)
    design = np.zeros((nrow, 3 * nblock))
    for harmonic in range(3):
        np.add.at(design, (mc_rows, 3 * mc_blocks + harmonic), basis[:, harmonic])
    y = np.bincount(rows, weights=weights, minlength=nrow)
    # Repeated cache rows from one event are correlated. Square the event's
    # total row weight, rather than treating each appearance as independent.
    all_data_ids = np.asarray([json.dumps([r, e], separators=(",", ":"))
                               for r, e in zip(data["run_number"], data["event_id"])], dtype=str)
    data_ids = all_data_ids[data_keep]
    group = np.asarray([f"{row}:{event}" for row, event in zip(rows, data_ids)])
    _, inverse = np.unique(group, return_inverse=True)
    group_weights = np.bincount(inverse, weights=weights)
    group_rows = np.zeros(len(group_weights), dtype=np.int64)
    group_rows[inverse] = rows
    sumw2 = np.bincount(group_rows, weights=group_weights**2, minlength=nrow)
    fixed_rows=mc_row_all[fixed_keep];fixed_ids=mc["event_id"][fixed_keep]
    fixed_weights=mc["nominal_weight"][fixed_keep]
    fixed_keys=np.asarray([f"{row}:{event}" for row,event in zip(fixed_rows,fixed_ids)])
    _,fixed_inverse=np.unique(fixed_keys,return_inverse=True)
    fixed_group_weights=np.bincount(fixed_inverse,weights=fixed_weights)
    fixed_group_rows=np.zeros(len(fixed_group_weights),dtype=np.int64);fixed_group_rows[fixed_inverse]=fixed_rows
    fixed_prediction=np.bincount(fixed_group_rows,weights=fixed_group_weights,minlength=nrow)
    fixed_mc_sumw2=np.bincount(fixed_group_rows,weights=fixed_group_weights**2,minlength=nrow)
    supported_rows = np.any(design[:, ::3] > 0, axis=1) | (fixed_prediction > 0)
    bad_rows = np.flatnonzero((y > 0) & ~supported_rows)
    if len(bad_rows):
        raise ValueError(f"Positive data yield outside MC response support in reco rows {bad_rows.tolist()}")
    if not (np.all(np.isfinite(design)) and np.all(np.isfinite(y)) and np.all(np.isfinite(sumw2))):
        raise ValueError("Event accumulation overflow or nonfinite response")
    epsilon_max = np.zeros(nblock)
    np.maximum.at(epsilon_max, mc_blocks, eps)
    all_blocks = _block_metadata(truth)
    all_blocks += [dict(global_id=truth["nblocks"] + i, kind="guard", name=name)
                   for i, name in enumerate(GUARD_NAMES)]
    blocks = [dict(all_blocks[int(gid)], local_id=i, published=int(gid) in publication,
                   fit_treatment="physics" if int(gid)<truth["nblocks"] else "fitted_tprime_feedin")
              for i, gid in enumerate(fitted_global)]
    parameters = [dict(index=3*i+h, local_block=i, global_block=block["global_id"],
                       harmonic=name, unit=COEFFICIENT_UNIT, published=block["published"])
                  for i, block in enumerate(blocks) for h, name in enumerate(HARMONICS)]
    reco_rows = [dict(block, row=block["global_id"]*(len(reco["p"])-1)+pi,
                      phi_index=pi, phi_bounds=reco["p"][pi:pi+2].tolist())
                 for block in _block_metadata(reco) for pi in range(len(reco["p"])-1)]
    return dict(design=design, y=y, sumw2=sumw2, data_rows=rows, data_weights=weights,
                data_ids=data_ids, mc_rows=mc_rows, mc_blocks=mc_blocks, mc_basis=basis,
                mc_ids=mc["event_id"][fitted_keep], epsilon_max=epsilon_max,
                fixed_prediction=fixed_prediction,fixed_mc_sumw2=fixed_mc_sumw2,
                fixed_mc_rows=fixed_rows,fixed_mc_weights=fixed_weights,fixed_mc_ids=fixed_ids,
                bootstrap_data_ids=np.unique(all_data_ids), bootstrap_mc_ids=np.unique(mc["event_id"]),
                published_blocks=np.asarray([i for i, b in enumerate(blocks) if b["published"]], dtype=int),
                active_global_blocks=fitted_global, truth_blocks=blocks, reco_rows=reco_rows,
                parameter_metadata=parameters, coefficient_unit=COEFFICIENT_UNIT,
                config=config, cache_manifest=manifest, cache_dir=str(cache_dir.resolve()),
                selection_counts=dict(data_cached=len(data_row_all), data_selected=int(data_keep.sum()),
                    mc_cached=len(mc_row_all), mc_selected=int(fitted_keep.sum()+fixed_keep.sum()),
                    mc_fitted=int(fitted_keep.sum()),mc_fixed_feedin=int(fixed_keep.sum()),
                    mc_zero_weight_unsupported=int(np.sum((mc_row_all >= 0) & ~fitted_keep & ~fixed_keep)),
                    fixed_truth_blocks=sorted(set(mc_global_all[fixed_keep].tolist())),
                    unsupported_truth_blocks=sorted(set(range(len(all_blocks))) - set(supported_global.tolist()))))


def mc_prediction_variance(problem, parameters):
    """Poissonized generator-event variance including all harmonic products.

    Repeated rows from the same generator event get one multiplier, retaining
    cross-block covariance. This is the row variance only; inter-row covariance
    remains available through joint event resampling, not a diagonal fit term.
    Fixed-total generator sampling instead needs its multinomial covariance.
    """
    parameters = np.asarray(parameters, dtype=float)
    if parameters.shape != (problem["design"].shape[1],) or not np.all(np.isfinite(parameters)):
        raise ValueError("parameters must contain one finite value per design column")
    fitted = np.einsum("ij,ij->i", problem["mc_basis"], parameters.reshape(-1, 3)[problem["mc_blocks"]])
    prediction=np.r_[fitted,np.asarray(problem["fixed_mc_weights"],dtype=float)]
    event_rows=np.r_[problem["mc_rows"],problem["fixed_mc_rows"]]
    event_ids=np.r_[problem["mc_ids"],problem["fixed_mc_ids"]]
    keys = np.asarray([f"{row}:{event}" for row, event in zip(event_rows,event_ids)])
    _, inverse = np.unique(keys, return_inverse=True)
    grouped = np.bincount(inverse, weights=prediction)
    rows = np.zeros(len(grouped), dtype=int)
    rows[inverse] = event_rows
    result = np.bincount(rows, weights=grouped**2, minlength=len(problem["y"]))
    if not np.all(np.isfinite(result)):
        raise ValueError("MC prediction variance overflow")
    return result


def _svd(matrix):
    u, s, vh = np.linalg.svd(matrix, full_matrices=matrix.shape[0] < matrix.shape[1])
    tolerance = np.finfo(float).eps * max(matrix.shape, default=0) * (s[0] if len(s) else 0)
    return u, s, vh, int(np.count_nonzero(s > tolerance)), float(tolerance)


def diagnostics(problem, parameters=None, row_variance=None):
    """JSON-safe response diagnostics; never rank-truncate the fitted model.

    Unit row weighting is the default, explicitly labeled. For a local
    Gaussian/Fisher diagnostic callers may supply a positive finite variance
    for *every* row; no empty row is silently excluded. Profiled information
    uses orthogonal projection onto the nuisance-column complement. It is
    not a boundary confidence covariance or a publication quality threshold.
    """
    design = np.asarray(problem["design"], dtype=float)
    if row_variance is None:
        whitened, weighting = design, "unit row weights; geometric response diagnostic"
    else:
        variance = np.asarray(row_variance, dtype=float)
        if variance.shape != (len(design),) or not np.all(np.isfinite(variance)) or np.any(variance <= 0):
            raise ValueError("row_variance must supply a finite positive value for every row")
        whitened, weighting = design / np.sqrt(variance)[:, None], "supplied row variances; local unconstrained diagnostic"
    scales = np.linalg.norm(whitened, axis=0)
    scaled = whitened / np.where(scales > 0, scales, 1)
    _, singular, vh, rank, tolerance = _svd(scaled)
    npar = design.shape[1]
    values = np.pad(singular, (0, max(0, npar-len(singular))))
    modes = []
    for index in range(max(0, npar-3), npar):
        direction = vh[index]
        modes.append(dict(singular_value=float(values[index]), scaled_parameter_direction=direction.tolist(),
                          coefficient_direction=(direction / np.where(scales > 0, scales, 1)).tolist()))
    publication_columns = np.asarray([3*b+h for b in problem["published_blocks"] for h in range(3)], dtype=int)
    nuisance_columns = np.setdiff1d(np.arange(npar), publication_columns)
    jp, jn = whitened[:, publication_columns], whitened[:, nuisance_columns]
    if jn.shape[1]:
        un, _, _, nuisance_rank, _ = _svd(jn)
        residual = jp - un[:, :nuisance_rank] @ (un[:, :nuisance_rank].T @ jp)
    else:
        nuisance_rank, residual = 0, jp.copy()
    _, sp, vp, publication_rank, _ = _svd(jp)
    _, profiled_s, _, _, _ = _svd(residual)
    # Compare projection residuals against the original information scale,
    # otherwise roundoff left by an exactly confounded direction looks full rank.
    profile_tolerance = np.finfo(float).eps * max(whitened.shape) * (sp[0] if len(sp) else 0)
    profiled_rank = int(np.count_nonzero(profiled_s > profile_tolerance))
    information_fractions = None
    if publication_rank == len(publication_columns):
        # SVD of Jp provides coordinates with known-nuisance information I.
        # No normal matrix inversion is formed, even for weak nuisance modes.
        transform = vp.T / sp
        retained_s = np.linalg.svd(residual @ transform, compute_uv=False)
        information_fractions = (retained_s**2).tolist()
    total_u = float(np.sum(design[:, ::3]))
    block_diagnostics = []
    parameters = None if parameters is None else np.asarray(parameters, dtype=float)
    if parameters is not None and (parameters.shape != (npar,) or not np.all(np.isfinite(parameters))):
        raise ValueError("parameters must contain one finite value per design column")
    fitted = None if parameters is None else np.asarray(problem["fixed_prediction"]) + design @ parameters
    for i, block in enumerate(problem["truth_blocks"]):
        selected = problem["mc_blocks"] == i
        w = problem["mc_basis"][selected, 0]
        _, inverse = np.unique(problem["mc_ids"][selected], return_inverse=True)
        grouped = np.bincount(inverse, weights=w)
        response = float(w.sum())
        neff = float(response**2 / np.dot(grouped, grouped)) if np.any(grouped) else 0.0
        entry = dict(block, mc_cache_rows=int(selected.sum()), mc_independent_events=len(grouped),
                     u_response=response, u_response_effective_events=neff,
                     u_kernel_fraction=response / total_u if total_u > 0 else None)
        if parameters is not None:
            contribution = design[:, 3*i:3*i+3] @ parameters[3*i:3*i+3]
            entry["fitted_yield"] = float(contribution.sum())
            entry["fitted_fraction_by_reco_row"] = [float(c/f) if f > 0 else None for c, f in zip(contribution, fitted)]
        block_diagnostics.append(entry)
    return dict(weighting=weighting, rows=design.shape[0], parameters=npar, rank=rank,
                rank_tolerance=tolerance, condition_number=float(values[0]/values[-1]) if rank == npar and npar else None,
                column_norms=scales.tolist(), scaled_singular_values=values.tolist(), weakest_modes=modes,
                publication_columns=publication_columns.tolist(), nuisance_rank=nuisance_rank,
                publication_rank_known_nuisance=publication_rank, publication_rank_profiled=profiled_rank,
                profiled_singular_values=profiled_s.tolist(), retained_publication_information_fractions=information_fractions,
                blocks=block_diagnostics,fixed_feedin_by_reco_row=np.asarray(
                    problem.get("fixed_prediction", np.zeros(design.shape[0]))).tolist(),
                interpretation="Kernel fractions are not efficiencies; local information is not confidence-interval coverage.")
