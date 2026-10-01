#!/usr/bin/env python3
"""Independent NumPy validation of exported migration fits; never modifies inputs.

Run in the Hall C/NPS environment for optional --sim-file/--vertex-file event
checks. With only an output-directory argument, NumPy and the standard library
suffice. Deterministic y=X*p closure checks algebra, NOT detector/physics bias.
The optional disjoint event split instead tests statistical closure using an
independent response and pseudo-data sample drawn from the existing simulation.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


def records(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def column(rows, name, dtype=float):
    return np.array([r[name] for r in rows], dtype=dtype)


def solve(x, y, variance, tolerance=1e-10):
    """Independent LAPACK SVD; fail rather than truncate unsupported directions."""
    w = x / np.sqrt(variance[:, None])
    scale = np.linalg.norm(w, axis=0)
    if np.any(scale == 0) or len(y) <= x.shape[1]:
        raise ValueError("unsupported response")
    u, s, vh = np.linalg.svd(w / scale, full_matrices=False)
    if np.count_nonzero(s > tolerance * s[0]) != x.shape[1]:
        raise ValueError("rank-deficient response")
    p = (vh.T @ ((u.T @ (y / np.sqrt(variance))) / s)) / scale
    a = vh.T / s
    covariance = (a @ a.T) / np.outer(scale, scale)
    return p, covariance, s


def relative_error(actual, expected):
    return float(np.linalg.norm(actual - expected) / max(np.linalg.norm(expected), 1e-300))


def manufactured(parameters, zero_interference=False):
    """Known binwise constants with slopes and distinct exterior responses."""
    p = np.zeros(len(parameters))
    for i, row in enumerate(parameters):
        b = int(row['truth_block'])
        if row['region'] == 'published':
            it, iq, ix = (int(row[k]) for k in ('it', 'iq', 'ix'))
            values = [25 + 9 * it + 3 * iq + 2 * ix, (-1)**ix * (2 + .4 * it), -3 + .2 * iq]
        else:
            values = [70 + 3 * b, (-1)**b * 5, -7 + .1 * b]
        component = {'U': 0, 'LT': 1, 'TT': 2}[row['component']]
        p[i] = (0 if zero_interference and component else values[component]) * 1e-9
    return p


def event_split(args, parameters, all_x, rows, report):
    """Independently rebuild event response, then split by generated event ID.

    The input simulation supplies one accepted row per generated event. Both
    halves are doubled to the original luminosity; all variances therefore
    receive a factor four. Holding bin edges fixed isolates response statistics
    from the separate uncertainty of data-adaptive binning.
    """
    import ROOT
    sim_names = ['event_id', 'Q2', 'xB', 't', 'tmin', 'phi', 'mmiss', 'sigcm', 'full_weight', 'is_exclusive']
    raw_names = ['Q2i', 'Wi', 'ti', 'phipqi', 'sigcm']
    sim = ROOT.RDataFrame('simulation', str(args.sim_file)).AsNumpy(sim_names)
    raw = ROOT.RDataFrame('h10', str(args.vertex_file)).AsNumpy(raw_names)
    sim = {k: np.asarray(v).astype(float) if k != 'event_id' else np.asarray(v).astype(np.int64) for k, v in sim.items()}
    raw = {k: np.asarray(v).astype(float) for k, v in raw.items()}
    slices = records(args.directory / 'excl_xsec_pi0_analysis_no_simc_model_slice_summary.csv')
    nt, nq, nx = (1 + max(int(r[k]) for r in slices) for k in ('it', 'iq', 'ix'))
    np_phi = 1 + max(int(r['ip']) for r in rows)
    te = np.array(sorted({float(r[k]) for r in slices for k in ('tprime_lo', 'tprime_hi')}))
    qe = np.array(sorted({float(r[k]) for r in slices for k in ('q2_lo', 'q2_hi')}))
    xe = [np.array(sorted({float(r[k]) for r in slices if int(r['iq']) == iq for k in ('xb_lo', 'xb_hi')})) for iq in range(nq)]
    pe = np.array(sorted({float(r[k]) for r in rows for k in ('phi_lo', 'phi_hi')}))
    tp = sim['t'] - sim['tmin']
    selected = ((sim['mmiss'] >= args.mmiss_lower) & (sim['mmiss'] <= args.mmiss_upper) &
                (sim['Q2'] >= qe[0]) & (sim['Q2'] <= qe[-1]) &
                (sim['xB'] >= min(e[0] for e in xe)) & (sim['xB'] <= max(e[-1] for e in xe)) &
                (tp >= te[0]) & (tp <= te[-1]) & (sim['is_exclusive'] != 0) &
                (np.abs(sim['sigcm']) >= float(np.float32(1e-20))) & np.isfinite(sim['phi']) &
                np.isfinite(sim['full_weight']) & np.isfinite(sim['sigcm']))
    sim = {k: v[selected] for k, v in sim.items()}
    tp = tp[selected]
    ids = sim['event_id']
    assert np.all((ids >= 0) & (ids < len(raw['Q2i'])))
    assert np.all(np.abs(raw['sigcm'][ids] - sim['sigcm']) <=
                  1e-6 * np.maximum(np.abs(raw['sigcm'][ids]), np.abs(sim['sigcm'])) + 1e-20)
    # Inclusive final edge and right-open interior edges match histogram bins.
    def bins(edges, values):
        return np.minimum(np.searchsorted(edges, values, side='right') - 1, len(edges) - 2)
    rq, rt = bins(qe, sim['Q2']), bins(te, tp)
    rx = np.zeros(len(ids), dtype=int)
    for q in range(nq):
        rx[rq == q] = bins(xe[q], sim['xB'][rq == q])
    rphi = bins(pe, np.mod(sim['phi'], 2 * np.pi))
    rr = ((rt * nq + rq) * nx + rx) * np_phi + rphi
    mass, pion, energy = .9382720813, .1349768, args.ebeam
    q, w, t, phi = (raw[k][ids] for k in ('Q2i', 'Wi', 'ti', 'phipqi'))
    x = q / (w * w - mass * mass + q)
    q0 = (w * w - mass * mass - q) / (2 * w)
    epi = (w * w + pion * pion - mass * mass) / (2 * w)
    forward_t = pion * pion - q - 2 * (q0 * epi - np.sqrt(q0*q0+q) * np.sqrt(np.maximum(0, epi*epi-pion*pion)))
    generated_tp = -t - forward_t
    tq = bins(qe, q)
    tx = np.zeros(len(ids), dtype=int)
    for j in range(nq):
        tx[tq == j] = bins(xe[j], x[tq == j])
    truth = (bins(te, generated_tp) * nq + tq) * nx + tx
    published = nt * nq * nx
    # Assign guards from lowest to highest priority, so corners end in t',Q2,xB order.
    for j in range(nq):
        truth[(tq == j) & (x < xe[j][0])] = published + 4
        truth[(tq == j) & (x > xe[j][-1])] = published + 5
    truth[q < qe[0]] = published + 2
    truth[q > qe[-1]] = published + 3
    truth[generated_tp < te[0]] = published
    truth[generated_tp > te[-1]] = published + 1
    active = {int(r['truth_block']): int(r['active_block_index']) for r in parameters}
    ab = np.array([active[int(b)] for b in truth])
    y = q / (2 * mass * x * energy)
    eps = (1-y-q/(4*energy*energy))/(1-y+y*y/2+q/(4*energy*energy))
    phi = np.mod(phi, 2*np.pi)
    basis = (sim['full_weight']/sim['sigcm'])[:, None] * np.column_stack(
        [np.ones(len(ids)), np.sqrt(2*eps*(1+eps))*np.cos(phi), eps*np.cos(2*phi)])/(2*np.pi)
    rebuilt = np.zeros_like(all_x)
    for component in range(3):
        np.add.at(rebuilt, (rr, 3*ab+component), basis[:, component])
    response_error = relative_error(rebuilt, all_x)
    assert response_error < 1e-10, ('independent event response differs', response_error)
    fitted_p = column(parameters, 'value').reshape(-1, 3)
    fitted_event_y = np.sum(basis * fitted_p[ab], axis=1)
    rebuilt_mc_variance = np.bincount(rr, weights=fitted_event_y**2, minlength=len(rows))
    mc_variance_error = relative_error(rebuilt_mc_variance, column(rows, 'mc_variance_at_final'))
    assert mc_variance_error < 1e-10, ('independent finite-MC variance differs', mc_variance_error)
    # Split unique event IDs so duplicate accepted records could never enter
    # both halves. Seed and event count are persisted for repeatability.
    unique, inverse = np.unique(ids, return_inverse=True)
    train = (np.random.default_rng(args.seed).random(len(unique)) < .5)[inverse]
    train_x = np.zeros_like(all_x)
    truth_p = manufactured(parameters)
    event_y = np.sum(basis * truth_p.reshape(-1, 3)[ab], axis=1)
    pseudo_y = np.bincount(rr[~train], weights=2*event_y[~train], minlength=len(rows))
    test_variance = np.bincount(rr[~train], weights=(2*event_y[~train])**2, minlength=len(rows))
    for component in range(3):
        np.add.at(train_x, (rr[train], 3*ab[train]+component), 2*basis[train, component])
    # Known-injection train variance isolates sampling closure without sharing
    # pseudo-data and response events. This is not the real-data covariance.
    train_variance = np.bincount(rr[train], weights=(2*event_y[train])**2, minlength=len(rows))
    variance = test_variance + train_variance
    keep = (variance > 0) & np.any(train_x != 0, axis=1)
    split = dict(seed=args.seed, events=len(ids), unique_events=len(unique), train_events=int(train.sum()),
                 test_events=int((~train).sum()), independent_response_relative_error=response_error,
                 independent_mc_variance_relative_error=mc_variance_error,
                 unsupported_test_rows=int(np.count_nonzero((pseudo_y != 0) & ~np.any(train_x != 0, axis=1))))
    try:
        p, cov, _ = solve(train_x[keep], pseudo_y[keep], variance[keep])
        marginal_pulls = (p - truth_p)/np.sqrt(np.diag(cov))
        delta = p - truth_p
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        coefficient_chi2 = float(np.sum((eigenvectors.T@delta)**2/eigenvalues))
        split.update(status='completed_statistical_test', max_abs_marginal_pull=float(np.max(np.abs(marginal_pulls))),
                     rms_marginal_pull=float(np.sqrt(np.mean(marginal_pulls**2))),
                     joint_coefficient_chi2=coefficient_chi2, joint_coefficient_ndf=len(p),
                     residual_chi2=float(np.sum((pseudo_y[keep]-train_x[keep]@p)**2/variance[keep])),
                     residual_ndf=int(keep.sum()-len(p)))
    except ValueError as error:
        split.update(status='insufficient_independent_response_rank', reason=str(error))
    report['disjoint_event_split'] = split


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--json', type=Path, required=True)
    parser.add_argument('--sim-file', type=Path)
    parser.add_argument('--vertex-file', type=Path)
    parser.add_argument('--seed', type=int, default=20260916)
    parser.add_argument('--target-contam', type=float, default=.584)
    parser.add_argument('--target-contam-err', type=float, default=.014)
    parser.add_argument('--mmiss-lower', type=float, default=.6)
    parser.add_argument('--mmiss-upper', type=float, default=1.1)
    parser.add_argument('--ebeam', type=float, default=10.538)
    args = parser.parse_args()
    parameters = records(args.directory/'migration_parameters.csv')
    rows = records(args.directory/'migration_reco_rows.csv')
    nr, np_ = len(rows), len(parameters)
    x = np.zeros((nr, np_))
    for cell in records(args.directory/'migration_design.csv'):
        x[int(cell['reco_row']), int(cell['parameter_index'])] = float(cell['response'])
    order = np.argsort(column(rows, 'fit_index', int))
    order = order[column(rows, 'fit_index', int)[order] >= 0]
    fit_x, y, variance = x[order], column(rows, 'data')[order], column(rows, 'variance_used')[order]
    p, covariance, singular = solve(fit_x, y, variance)
    exported_p = column(parameters, 'value')
    exported_covariance, target_covariance = np.zeros_like(covariance), np.zeros_like(covariance)
    for cell in records(args.directory/'migration_covariance.csv'):
        i, j = int(cell['parameter_i']), int(cell['parameter_j'])
        exported_covariance[i, j] = float(cell['stat_plus_mc_covariance'])
        target_covariance[i, j] = float(cell['target_correlated_covariance'])
    expected_correlation = covariance / np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
    exported_correlation = np.zeros_like(covariance)
    for cell in records(args.directory/'migration_correlation.csv'):
        i, j = int(cell['parameter_i']), int(cell['parameter_j'])
        exported_correlation[i, j] = float(cell['correlation'])
    errors = dict(parameter_relative_error=relative_error(p, exported_p), covariance_relative_error=relative_error(covariance, exported_covariance))
    errors['correlation_relative_error'] = relative_error(expected_correlation, exported_correlation)
    assert max(errors.values()) < 1e-7, errors
    expected_target = np.outer(exported_p, exported_p)*(args.target_contam_err/args.target_contam)**2
    assert relative_error(target_covariance, expected_target) < 1e-13
    offdiagonal = covariance - np.diag(np.diag(covariance))
    assert np.any(np.abs(offdiagonal) > 0)
    closure = {}
    for name, zero in [('sloped_U_signed_LT_TT_and_distinct_guards', False), ('zero_LT_TT', True)]:
        injection = manufactured(parameters, zero)
        recovered, _, _ = solve(fit_x, fit_x@injection, variance)
        error = relative_error(recovered, injection)
        assert error < 1e-7, (name, error)
        closure[name] = error
    injection = manufactured(parameters)
    published_columns = np.array([r['region'] == 'published' for r in parameters])
    reduced, _, _ = solve(fit_x[:, published_columns], fit_x@injection, variance)
    guard_omission_bias = float(np.max(np.abs(reduced-injection[published_columns]))*1e9)
    assert guard_omission_bias > 1e-6, 'fixture does not demonstrate guard feed-in'
    scaled_p, scaled_cov, _ = solve(fit_x, 2*y, 4*variance)
    assert relative_error(scaled_p, 2*p) < 1e-10 and relative_error(scaled_cov, 4*covariance) < 1e-10
    duplicated = fit_x.copy(); duplicated[:, -1] = duplicated[:, 0]
    try:
        solve(duplicated, y, variance)
    except ValueError:
        rank_rejected = True
    else:
        rank_rejected = False
    assert rank_rejected
    report = dict(status='passed_deterministic_checks', reconstructed_rows=nr, fitted_rows=len(order), parameters=np_,
                  numpy_reproduction=errors, deterministic_algebra_closure=closure,
                  omit_guard_max_published_bias_nb_per_GeV2=guard_omission_bias,
                  target_covariance_relative_error=relative_error(target_covariance, expected_target),
                  target_scaling_check=True, rank_deficiency_rejected=True,
                  full_covariance_offdiagonal_norm_fraction=float(np.linalg.norm(offdiagonal)/np.linalg.norm(covariance)),
                  scaled_condition=float(singular[0]/singular[-1]),
                  limitation='No assertion of physical agreement or detector/model adequacy; same-matrix closure tests algebra only.')
    if bool(args.sim_file) != bool(args.vertex_file):
        parser.error('--sim-file and --vertex-file must be supplied together')
    if args.sim_file:
        event_split(args, parameters, x, rows, report)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
