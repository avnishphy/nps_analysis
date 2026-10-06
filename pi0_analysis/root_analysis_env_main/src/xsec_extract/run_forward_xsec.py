#!/usr/bin/env python3
"""Forward harmonic fit with independent grids and event-level MC replicas.

This implements conditional weighted-yield inference, not publication
certification: upstream signal-weight, detector and radiative uncertainties
and interval coverage require separate validation. Existing outputs are never
overwritten. See FORWARD_EXTRACTION.md for cache preparation and conventions.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import scipy

from forward_xsec_problem import build_problem, diagnostics, mc_prediction_variance
from forward_xsec_statistics import fit_problem, bootstrap_problem, profile_problem


def json_safe(value):
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path, value):
    temporary = path.with_name(path.name + '.partial')
    with temporary.open('x') as stream:
        json.dump(json_safe(value), stream, indent=2, allow_nan=False)
        stream.write('\n')
    os.replace(temporary, path)


def fingerprint(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return {'path': str(path.resolve()), 'sha256': h.hexdigest(), 'bytes': path.stat().st_size}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bootstrap', type=int, default=200,
                        help='Full data+MC Poisson event replicas; 0 is central-value diagnostics only')
    parser.add_argument('--seed', type=int, default=20261001)
    parser.add_argument('--maxiter', type=int, default=2000)
    parser.add_argument('--profiles', type=Path,
                        help='JSON list of {name,functional,values}; fixed-response profiles, no calibrated coverage claim')
    args = parser.parse_args(argv)
    if args.bootstrap < 0 or args.bootstrap == 1 or args.maxiter <= 0:
        parser.error('bootstrap must be 0 or at least 2 and maxiter positive')
    if args.output.exists():
        parser.error('Output already exists; choose a fresh directory')
    manifest_path = args.cache / 'forward_cache_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('complete') is not True or manifest.get('schema_version') != 2:
        parser.error('Cache is incomplete or has an unsupported schema')
    config = json.loads(args.config.read_text())
    if not (np.isfinite(manifest.get('target_divisor', np.nan)) and manifest['target_divisor'] > 0
            and np.isfinite(manifest.get('target_divisor_error', np.nan)) and manifest['target_divisor_error'] >= 0):
        parser.error('Cache requires a valid target divisor and its uncertainty')
    problem = build_problem(args.cache, config)
    args.output.mkdir(parents=True, exist_ok=False)
    source_dir = Path(__file__).resolve().parent
    provenance = {
        'argv': sys.argv if argv is None else argv,
        'config': config,
        'cache_manifest': manifest,
        'inputs': [fingerprint(p) for p in (args.config, manifest_path,
                                           args.cache/'data_events.csv', args.cache/'mc_events.csv')],
        'sources': [fingerprint(source_dir/name) for name in
                    ('run_forward_xsec.py', 'forward_xsec_problem.py', 'forward_xsec_statistics.py')],
        'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__, 'seed': args.seed,
    }
    for name in ('forward_export_provenance.json', 'forward_export_config.h'):
        if (args.cache/name).exists():
            provenance['inputs'].append(fingerprint(args.cache/name))
    if args.profiles:
        provenance['inputs'].append(fingerprint(args.profiles))
    write_json(args.output/'provenance.json', provenance)
    snapshot = args.output/'source_snapshot'
    snapshot.mkdir()
    for item in provenance['sources']:
        source = Path(item['path'])
        with (snapshot/source.name).open('xb') as stream:
            stream.write(source.read_bytes())
    fit_succeeded = False
    try:
        fit = fit_problem(problem, maxiter=args.maxiter)
        fit_succeeded = True
        write_json(args.output/'fit.json', fit)
        coefficients = np.asarray(fit['parameters'], dtype=float)
        write_json(args.output/'response_diagnostics.json', diagnostics(problem, parameters=coefficients))
        expected_data_variance = np.asarray(fit['row_scales']) * np.maximum(fit['predicted'], 0.)
        expected_mc_variance = mc_prediction_variance(problem, coefficients)
        variance_diagnostics = {
            'interpretation': 'local conditional row-variance decomposition, not the fitted likelihood or confidence covariance',
            'expected_data_variance': expected_data_variance,
            'poissonized_mc_prediction_variance': expected_mc_variance,
        }
        total_variance = expected_data_variance + expected_mc_variance
        if np.all(total_variance > 0):
            variance_diagnostics['local_profiled_information'] = diagnostics(
                problem, parameters=coefficients, row_variance=total_variance)
        write_json(args.output/'variance_diagnostics.json', variance_diagnostics)
        np.savez_compressed(args.output/'response_problem.npz', design=problem['design'],
                            y=problem['y'], sumw2=problem['sumw2'], epsilon_max=problem['epsilon_max'],
                            fixed_prediction=problem['fixed_prediction'],fixed_mc_sumw2=problem['fixed_mc_sumw2'])
        write_json(args.output/'bin_metadata.json', {
            key: problem[key] for key in ('truth_blocks', 'parameter_metadata', 'reco_rows',
                                         'published_blocks', 'active_global_blocks')})
        with (args.output/'coefficients.csv').open('x', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['parameter_index', 'active_truth_block', 'publication', 'component',
                             'coefficient_nb_per_GeV2'])
            published = set(problem['published_blocks'])
            for i, value in enumerate(coefficients):
                writer.writerow([i, i//3, i//3 in published, ('U','LT','TT')[i%3], value])
        target_fraction = manifest['target_divisor_error'] / manifest['target_divisor']
        write_json(args.output/'target_covariance.json', {
            'kind': 'separate common target normalization; first-order rank-one',
            'units': '(nb/GeV^2)^2',
            'covariance': np.outer(coefficients, coefficients)*target_fraction**2})
        replicas = None
        if args.bootstrap:
            replicas = bootstrap_problem(problem, args.bootstrap, args.seed,
                                         fit_options={'maxiter': args.maxiter})
            write_json(args.output/'conditional_bootstrap.json', replicas)
        if args.profiles:
            profiles = []
            for request in json.loads(args.profiles.read_text()):
                profiles.append({'name': request['name'], 'result': profile_problem(
                    problem, request['functional'], request['values'], fit_result=fit,
                    fit_options={'maxiter': args.maxiter})})
            write_json(args.output/'fixed_response_profiles.json', profiles)
        write_json(args.output/'status.json', {
            'complete': True, 'fit_succeeded': True, 'publication_ready': False,
            'inference': 'conditional scaled-Poisson; event bootstrap includes data and finite MC',
            'fit_parameterization': 'independent interior plus tprime_below U/LT/TT; Q2/xB feed-in fixed',
            'missing_publication_validation': [
                'Upstream pi0_weight/background-fit shared uncertainty',
                'Detector/radiative response and normalization validation',
                'Independent full-experiment interval coverage and within-bin/exterior shape stress tests'],
            'bootstrap_requested': args.bootstrap,
            'bootstrap_status': 'not_requested' if replicas is None else replicas['status'],
            'all_bootstrap_refits_succeeded': None if replicas is None else replicas['all_replicates_successful'],
            'profile_interpretation': 'nominal fixed-response likelihood scans, not calibrated confidence intervals',
        })
    except Exception as error:
        write_json(args.output/'status.json', {'complete': False, 'fit_succeeded': fit_succeeded,
                                              'publication_ready': False, 'error': str(error)})
        raise
    print(f'Forward fit complete: {args.output}; publication validation still required.')


if __name__ == '__main__':
    main()
