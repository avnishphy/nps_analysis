#!/usr/bin/env python3
"""Constrained forward inference for positive, event-weighted pi0 yields.

This is an explicit interim statistical model: scaled-Poisson row scales are
estimated from event weights and frozen within each fit. Empty rows borrow the
pooled scale. Conditional event bootstrap and fixed-response profile curves
are diagnostics, not calibrated confidence intervals or publication approval.
No background-subtracted or signed-weight input is accepted.
"""

import numpy as np
from scipy.optimize import LinearConstraint, minimize


class FitError(ValueError):
    """Input, identifiability, or optimizer failure; never silently regularized."""


def _require(condition, message):
    if not condition:
        raise FitError(message)


def angular_minimum(coefficients, epsilon):
    """Return exact min and cos(phi) for U+a LT cos(phi)+epsilon TT cos(2phi).

    At fixed epsilon this is quadratic in cos(phi). Its global angular minimum
    is nonincreasing in epsilon >= 0, so testing epsilon_max protects the whole
    sampled epsilon interval. Harmonic coefficients themselves remain signed.
    """
    u, lt, tt = np.asarray(coefficients, dtype=float)
    a = np.sqrt(2.0 * epsilon * (1.0 + epsilon))
    candidates = [-1.0, 1.0]
    if epsilon * tt > 0.0:
        vertex = -a * lt / (4.0 * epsilon * tt)
        if -1.0 < vertex < 1.0:
            candidates.append(float(vertex))
    values = [u + a * lt * x + epsilon * tt * (2.0 * x * x - 1.0)
              for x in candidates]
    index = int(np.argmin(values))
    return float(values[index]), candidates[index]


def _validate(problem):
    design = np.asarray(problem['design'], dtype=float)
    y = np.asarray(problem['y'], dtype=float)
    sumw2 = np.asarray(problem['sumw2'], dtype=float)
    epsilon = np.asarray(problem['epsilon_max'], dtype=float)
    _require(design.ndim == 2 and design.shape[1] > 0 and
             design.shape[1] % 3 == 0, 'design must have three columns per truth block')
    _require(y.shape == sumw2.shape == (design.shape[0],), 'yield array shape mismatch')
    _require(epsilon.shape == (design.shape[1] // 3,), 'epsilon_max shape mismatch')
    _require(all(np.all(np.isfinite(v)) for v in (design, y, sumw2, epsilon)),
             'non-finite fit input')
    _require(np.all(y >= 0) and np.all(sumw2 >= 0),
             'scaled-Poisson inference requires nonnegative, unsubtracted event weights')
    _require(np.all((y > 0) == (sumw2 > 0)), 'yield and sumw2 zero patterns disagree')
    _require(np.all((epsilon >= 0) & (epsilon <= 1)), 'epsilon must lie in [0,1]')
    _require(y.sum() > 0, 'no positive data yield')
    _require(np.all(design[:, ::3] >= 0), 'negative unpolarized response contribution')
    supported = np.any(design != 0, axis=1)
    _require(not np.any((y > 0) & ~supported), 'positive data yield in row without MC support')
    pooled = float(sumw2.sum() / y.sum())
    scales = np.full_like(y, pooled)
    np.divide(sumw2, y, out=scales, where=y > 0)
    return design, y, sumw2, epsilon, scales


def _angular_row(npar, block, epsilon, x):
    row = np.zeros(npar)
    row[3 * block:3 * block + 3] = (
        1.0, np.sqrt(2.0 * epsilon * (1.0 + epsilon)) * x,
        epsilon * (2.0 * x * x - 1.0))
    return row


def _fit(problem, initial=None, maxiter=2000, tolerance=1e-8, equality=None):
    design, y, sumw2, epsilon, scales = _validate(problem)
    nrow, npar = design.shape
    # Row scaling estimates variance for numerical conditioning only. It does
    # not change the scaled-Poisson objective or drop empty observations.
    numerical_sigma = np.sqrt(scales * np.maximum(y, scales))
    weighted = design / numerical_sigma[:, None]
    norms = np.linalg.norm(weighted, axis=0)
    _require(np.all(norms > 0), 'unconstrained response column; merge or measure that truth region')
    parameter_scale = 1.0 / norms
    singular = np.linalg.svd(weighted * parameter_scale, compute_uv=False)
    _require(len(singular) == npar and singular[-1] > 1e-12 * singular[0],
             'rank-deficient response; no SVD truncation or hidden regularization is applied')
    scaled_design = design * parameter_scale
    if initial is None:
        start = np.zeros(npar)
        start[::3] = y.sum() / design[:, ::3].sum()
    else:
        start = np.asarray(initial, dtype=float)
        _require(start.shape == (npar,) and np.all(np.isfinite(start)), 'invalid initial parameters')
    z = start / parameter_scale
    positive = y > 0
    observed_effective = y / scales
    mu_floor = float(y.mean()) * 1e-13
    _require(mu_floor > 0, 'yield scale below floating-point numerical range')

    def objective(v):
        mu = scaled_design @ v
        safe_mu = np.maximum(mu, mu_floor * 1e-6)
        terms = mu / scales - observed_effective
        terms[positive] += observed_effective[positive] * np.log(y[positive] / safe_mu[positive])
        return float(terms.sum())

    def gradient(v):
        mu = np.maximum(scaled_design @ v, mu_floor * 1e-6)
        return scaled_design.T @ ((1.0 - y / mu) / scales)

    angular_rows = [_angular_row(npar, block, eps, x)
                    for block, eps in enumerate(epsilon)
                    for x in (-1.0, -0.5, 0.0, 0.5, 1.0)]
    # Explicit prediction constraints also protect profiles with trial negative
    # yields. Entirely unsupported empty rows contribute zero and stay recorded.
    supported = np.any(design != 0, axis=1)
    prediction = LinearConstraint(scaled_design[supported] / numerical_sigma[supported, None],
                                 np.where(positive[supported], mu_floor, 0.0)
                                 / numerical_sigma[supported], np.inf)
    fixed = []
    if equality is not None:
        functional, value = equality
        row = np.asarray(functional, dtype=float) * parameter_scale
        _require(row.shape == (npar,) and np.any(row != 0) and np.all(np.isfinite(row)),
                 'invalid profile functional')
        _require(np.isfinite(value), 'non-finite profile value')
        fixed = [LinearConstraint(row[None, :], value, value)]
    iterations = 0
    for cutting_iteration in range(50):
        positivity = LinearConstraint(np.asarray(angular_rows) * parameter_scale, 0.0, np.inf)
        result = minimize(objective, z, jac=gradient, method='SLSQP',
                          constraints=[positivity, prediction] + fixed,
                          options={'ftol': tolerance, 'maxiter': maxiter, 'disp': False})
        iterations += int(result.nit)
        _require(result.success, 'optimizer failed: ' + str(result.message))
        z = result.x
        parameters = z * parameter_scale
        minima = [angular_minimum(parameters[3 * block:3 * block + 3], eps)
                  for block, eps in enumerate(epsilon)]
        violation = False
        for block, (minimum, x) in enumerate(minima):
            threshold = tolerance * max(1.0, float(np.max(np.abs(parameters[3 * block:3 * block + 3]))))
            if minimum < -threshold:
                angular_rows.append(_angular_row(npar, block, epsilon[block], x))
                violation = True
        if not violation:
            break
    else:
        raise FitError('angular cutting-plane constraints did not converge')
    predicted = design @ parameters
    _require(np.all(predicted >= -1e-8 * float(y.max())), 'negative fitted prediction')
    _require(np.all(predicted[positive] > 0), 'nonpositive fitted prediction in occupied row')
    if equality is not None:
        attained = float(np.asarray(equality[0]) @ parameters)
        _require(abs(attained - equality[1]) <= 1e-6 * max(1.0, abs(equality[1])),
                 'profile equality constraint not satisfied')
    return {
        'status': 'converged', 'parameters': parameters.tolist(),
        'predicted': predicted.tolist(), 'objective': objective(z),
        'deviance': 2.0 * objective(z), 'iterations': iterations,
        'angular_cutting_iterations': cutting_iteration + 1,
        'angular_minima': [v[0] for v in minima],
        'angular_minimum_cosphi': [v[1] for v in minima],
        'positivity_boundary_blocks': [b for b, v in enumerate(minima)
            if v[0] <= 1e-6 * max(1.0, abs(parameters[3 * b]))],
        'scaled_condition_number': float(singular[0] / singular[-1]),
        'scaled_singular_values': singular.tolist(),
        'row_scales': scales.tolist(), 'pooled_empty_row_scale': float(sumw2.sum() / y.sum()),
        'empty_rows_borrowing_pooled_scale': np.flatnonzero(y == 0).tolist(),
        'empty_rows_without_mc_support': np.flatnonzero((y == 0) & ~supported).tolist(),
        'included_rows': nrow, 'parameter_count': npar,
        'inference_status': 'INTERIM_SCALED_POISSON_FROZEN_WEIGHT_SCALES',
        'publication_ready': False,
    }


def fit_problem(problem, *, initial=None, maxiter=2000, tolerance=1e-8):
    """Fit every truth block, including exterior nuisance blocks, without priors."""
    return _fit(problem, initial=initial, maxiter=maxiter, tolerance=tolerance)


def _event_arrays(problem):
    design, y, sumw2, _, _ = _validate(problem)
    data_rows = np.asarray(problem['data_rows'], dtype=int)
    data_weights = np.asarray(problem['data_weights'], dtype=float)
    mc_rows = np.asarray(problem['mc_rows'], dtype=int)
    mc_blocks = np.asarray(problem['mc_blocks'], dtype=int)
    mc_basis = np.asarray(problem['mc_basis'], dtype=float)
    mc_ids = np.asarray(problem['mc_ids'])
    _require(data_rows.shape == data_weights.shape and data_rows.ndim == 1,
             'data event arrays mismatch')
    _require(np.all(np.isfinite(data_weights)) and np.all(data_weights >= 0),
             'bootstrap requires nonnegative, unsubtracted data event weights')
    _require(np.all((data_rows >= 0) & (data_rows < len(y))), 'invalid data event row')
    data_ids = np.asarray(problem.get('data_ids', np.arange(len(data_rows))))
    _require(data_ids.shape == data_rows.shape, 'data_ids shape mismatch')
    # Multiple records of one original event are correlated. Sum their weights
    # within each row before computing its sumw2 or assigning a multiplier.
    unique_ids, id_inverse = np.unique(data_ids, return_inverse=True)
    group_keys, group_inverse = np.unique(id_inverse * len(y) + data_rows, return_inverse=True)
    data_weights = np.bincount(group_inverse, weights=data_weights)
    data_rows = group_keys % len(y)
    data_ids = unique_ids[group_keys // len(y)]
    _require(mc_rows.shape == mc_blocks.shape == mc_ids.shape and mc_rows.ndim == 1 and
             mc_basis.shape == (len(mc_rows), 3), 'MC event arrays mismatch')
    _require(np.all((mc_rows >= 0) & (mc_rows < len(y))) and
             np.all((mc_blocks >= 0) & (mc_blocks < design.shape[1] // 3)), 'invalid MC event index')
    _require(np.all(np.isfinite(mc_basis)) and np.all(mc_basis[:, 0] >= 0),
             'invalid MC event basis')
    # Every snapshot must reproduce the estimator before resampling. Otherwise
    # a filtered or misnormalized cache could yield plausible, incorrect errors.
    _require(np.allclose(np.bincount(data_rows, weights=data_weights, minlength=len(y)), y,
                         rtol=2e-10, atol=1e-14 * np.max(y)), 'data cache does not reproduce nominal yield')
    _require(np.allclose(np.bincount(data_rows, weights=data_weights**2, minlength=len(y)), sumw2,
                         rtol=2e-10, atol=1e-14 * np.max(sumw2)), 'data cache does not reproduce nominal sumw2')
    rebuilt = np.zeros_like(design)
    for harmonic in range(3):
        np.add.at(rebuilt, (mc_rows, 3 * mc_blocks + harmonic), mc_basis[:, harmonic])
    _require(np.allclose(rebuilt, design, rtol=2e-10, atol=1e-14 * np.max(np.abs(design))),
             'MC cache does not reproduce nominal response normalization')
    return data_rows, data_weights, mc_rows, mc_blocks, mc_basis, mc_ids, data_ids


def _resampling_plan(problem):
    arrays = _event_arrays(problem)
    data_rows, weights, mc_rows, mc_blocks, basis, mc_ids, data_ids = arrays

    def catalog_indices(ids, catalog_key):
        catalog = np.unique(np.asarray(problem.get(catalog_key, ids)))
        _require(len(catalog) > 0, 'empty bootstrap event catalog')
        indices = np.searchsorted(catalog, ids)
        _require(np.all(indices < len(catalog)), 'selected event absent from bootstrap catalog')
        _require(np.array_equal(catalog[indices], ids), 'selected event absent from bootstrap catalog')
        return len(catalog), indices

    return (arrays, catalog_indices(data_ids, 'bootstrap_data_ids'),
            catalog_indices(mc_ids, 'bootstrap_mc_ids'))


def _resample_from_plan(problem, rng, plan, resample_data=True, resample_mc=True):
    arrays, (ndata, data_index), (nmc, mc_index) = plan
    data_rows, weights, mc_rows, mc_blocks, basis, _, _ = arrays
    design = np.asarray(problem['design'], dtype=float)
    result = dict(problem)
    if resample_data:
        multiplier = rng.poisson(1.0, ndata)[data_index]
        result['y'] = np.bincount(data_rows, weights=weights * multiplier, minlength=len(design))
        result['sumw2'] = np.bincount(data_rows, weights=weights**2 * multiplier, minlength=len(design))
    if resample_mc:
        multiplier = rng.poisson(1.0, nmc)[mc_index]
        rebuilt = np.zeros_like(design)
        for harmonic in range(3):
            np.add.at(rebuilt, (mc_rows, 3 * mc_blocks + harmonic), basis[:, harmonic] * multiplier)
        result['design'] = rebuilt
    return result


def resample_problem(problem, rng, *, resample_data=True, resample_mc=True):
    """Independent Poisson multipliers per original event, shared across copies.

    MC normalization stays fixed. This is a Poissonized empirical bootstrap
    conditional on the accepted-event cache, generated phase space, and model.
    It cannot recover absent MC support or uncertainty in normalization inputs.
    Optional complete bootstrap_*_ids catalogs fix random-draw ordering across
    binning configurations, including changes to their selected event domains.
    """
    return _resample_from_plan(problem, rng, _resampling_plan(problem), resample_data, resample_mc)


def bootstrap_problem(problem, replicates, seed, *, fit_options=None):
    """Full refits with fluctuating data and response; outputs are conditional spreads."""
    _require(isinstance(replicates, (int, np.integer)) and replicates >= 2,
             'at least two bootstrap replicates required')
    options = dict(fit_options or {})
    nominal = fit_problem(problem, **options)
    # Validate and group original-event IDs once; sorting multi-million-event
    # catalogs in every replica would dominate the full-fit resampling cost.
    plan = _resampling_plan(problem)
    options.setdefault('initial', nominal['parameters'])
    rng = np.random.default_rng(seed)
    samples, failures, changed_responses, boundary_samples = [], [], 0, 0
    successful_indices = []
    for index in range(replicates):
        try:
            replica = _resample_from_plan(problem, rng, plan)
            changed_responses += int(not np.array_equal(replica['design'], problem['design']))
            fitted = fit_problem(replica, **options)
            samples.append(fitted['parameters'])
            successful_indices.append(index)
            boundary_samples += bool(fitted['positivity_boundary_blocks'])
        except (FitError, np.linalg.LinAlgError) as error:
            failures.append({'replicate': index, 'reason': str(error)})
    result = {'requested': int(replicates), 'successful': len(samples), 'seed': int(seed),
              'failures': failures, 'response_changed_replicates': changed_responses,
              'boundary_replicates': boundary_samples, 'samples': samples,
              'successful_replicate_indices': successful_indices,
              'nominal_parameters': nominal['parameters'],
              'inference_status': 'CONDITIONAL_POISSON_EVENT_BOOTSTRAP_NOT_COVERAGE_CALIBRATED',
              'publication_ready': False,
              'normalization_convention': 'fixed MC normalization; Poisson original-event multipliers',
              'spread_conditioning': 'successful refits only; inspect all failures before use',
              'limitations': [
                  'Originally empty data rows remain empty under empirical event resampling.',
                  'Unsampled MC support cannot be generated by resampling accepted events.',
                  'Frozen-scale weighted-data approximation is not a full generative likelihood.',
                  'Normalization, detector, radiative and background uncertainties are excluded.',
                  'Quantiles are conditional spreads, not coverage-calibrated confidence limits.'],
              'covariance': None, 'quantiles': None,
              'all_replicates_successful': not failures}
    result['status'] = 'complete' if not failures else 'REQUIRES_REEVALUATION_FAILED_REFITS'
    if len(samples) >= 2:
        result['covariance'] = np.cov(np.asarray(samples), rowvar=False, ddof=1).tolist()
        result['quantiles'] = dict(zip(('0.025', '0.16', '0.5', '0.84', '0.975'),
            np.quantile(samples, [0.025, 0.16, 0.5, 0.84, 0.975], axis=0).tolist()))
    return result


def profile_problem(problem, functional, values, *, fit_result=None, fit_options=None):
    """Profile a coefficient or linear bin integral, refitting all nuisance terms.

    The response and observed weight scales are fixed. The curve is nominal;
    Wilks thresholds are not calibrated with boundaries, finite MC, or weights.
    """
    options = dict(fit_options or {})
    nominal = fit_result if fit_result is not None else fit_problem(problem, **options)
    options.setdefault('initial', nominal['parameters'])
    functional = np.asarray(functional, dtype=float)
    _require(functional.shape == (np.asarray(problem['design']).shape[1],),
             'profile functional shape mismatch')
    points = []
    for value in values:
        try:
            fit = _fit(problem, equality=(functional, float(value)), **options)
            delta = 2.0 * (fit['objective'] - nominal['objective'])
            _require(delta >= -1e-5 * max(1.0, nominal['objective']),
                     'profile improved nominal optimum; nominal fit requires investigation')
            points.append({'value': float(value), 'delta_deviance': max(0.0, delta),
                           'parameters': fit['parameters'], 'status': 'converged'})
        except (FitError, np.linalg.LinAlgError) as error:
            points.append({'value': float(value), 'delta_deviance': None,
                           'status': 'failed', 'reason': str(error)})
    return {'functional': functional.tolist(),
            'nominal_value': float(functional @ np.asarray(nominal['parameters'])),
            'points': points, 'inference_status': 'FIXED_RESPONSE_PROFILE_UNCALIBRATED_THRESHOLDS',
            'publication_ready': False}
