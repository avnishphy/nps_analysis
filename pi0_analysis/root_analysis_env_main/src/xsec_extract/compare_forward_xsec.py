#!/usr/bin/env python3
"""Compare aligned tprime integrals using paired conditional event bootstraps.

Diagnostic only: conditional bootstrap spreads do not establish confidence
coverage. Both fits must use the same complete event cache and RNG convention.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np


def _need(condition, message):
    if not condition:
        raise ValueError(message)


def _bounds(value, name):
    result = np.asarray(value, dtype=float)
    _need(result.shape == (2,) and np.all(np.isfinite(result)) and result[0] < result[1],
          f'{name} must have two finite increasing bounds')
    return result


def _same(a, b):
    return np.allclose(a, b, rtol=0, atol=1e-12)


def functional(blocks, nparameters, observable):
    """Whole published truth-bin integrals only; no sub-bin interpolation."""
    _need(observable.get('component') in ('U', 'LT', 'TT'), 'Unknown harmonic component')
    q2, xb, tp = [_bounds(observable[k], k) for k in ('q2', 'xb', 'tprime')]
    _need(nparameters == 3 * len(blocks), 'Parameter/truth-block dimension mismatch')
    _need(sorted(b['local_id'] for b in blocks) == list(range(len(blocks))),
          'Truth local_id values must be unique and contiguous')
    result = np.zeros(nparameters)
    selected = []
    harmonic = ('U', 'LT', 'TT').index(observable['component'])
    for block in blocks:
        if block['kind'] != 'interior' or not block['published']:
            continue
        if not (_same(_bounds(block['q2_bounds'], 'q2'), q2) and
                _same(_bounds(block['xb_bounds'], 'xb'), xb)):
            continue
        low, high = _bounds(block['tprime_bounds'], 'tprime')
        overlap = min(high, tp[1]) - max(low, tp[0])
        if overlap <= 1e-12:
            continue
        _need(_same([max(low, tp[0]), min(high, tp[1])], [low, high]),
              'Partial truth-bin overlap; use common aligned tprime integrals')
        selected.append((low, high))
        result[3 * block['local_id'] + harmonic] = overlap
    selected.sort()
    cursor = tp[0]
    for low, high in selected:
        _need(_same(low, cursor), 'Published tprime coverage has a gap or overlapping bins')
        cursor = high
    _need(bool(selected) and _same(cursor, tp[1]),
          'Requested integral lacks full published coverage at exactly matching q2/xb bounds')
    return result


def _hash_record(records, basename):
    hashes = [r.get('sha256') for r in records if Path(r['path']).name == basename]
    _need(len(hashes) == 1 and isinstance(hashes[0], str) and len(hashes[0]) == 64,
          f'Missing or ambiguous SHA256 provenance: {basename}')
    return hashes[0]


def _read_run(path):
    names = ('fit.json', 'bin_metadata.json', 'conditional_bootstrap.json', 'provenance.json')
    contents = {name: (path / name).read_bytes() for name in names}
    fit, bins, boot, provenance = [json.loads(contents[name]) for name in names]
    parameters = np.asarray(fit['parameters'], dtype=float)
    samples = np.asarray(boot['samples'], dtype=float)
    indices = boot['successful_replicate_indices']
    _need(fit['status'] == 'converged' and parameters.ndim == 1 and np.all(np.isfinite(parameters)),
          'Fit must have finite converged parameters')
    _need(samples.ndim == 2 and samples.shape == (len(indices), len(parameters)) and
          np.all(np.isfinite(samples)), 'Invalid bootstrap sample dimensions or values')
    _need(type(boot['requested']) is int and boot['requested'] >= 2 and
          all(type(i) is int and 0 <= i < boot['requested'] for i in indices) and
          len(set(indices)) == len(indices) and boot['successful'] == len(indices),
          'Invalid successful bootstrap replicate indices')
    _need(boot['seed'] == provenance['seed'], 'Bootstrap/provenance seed mismatch')
    _need(np.allclose(boot['nominal_parameters'], parameters, rtol=1e-8, atol=1e-10),
          'Bootstrap nominal parameters disagree with saved fit')
    return dict(parameters=parameters, samples=samples, indices=indices, bootstrap=boot,
                blocks=bins['truth_blocks'], provenance=provenance,
                files={name: hashlib.sha256(contents[name]).hexdigest() for name in names})


def compare(run_a, run_b, observables):
    a, b = [_read_run(Path(path)) for path in (run_a, run_b)]
    _need(isinstance(observables, list) and observables, 'Provide a nonempty observable list')
    names = [o.get('name') for o in observables]
    _need(all(isinstance(n, str) and n for n in names) and len(set(names)) == len(names),
          'Observable names must be unique nonempty strings')
    for name in ('data_events.csv', 'mc_events.csv', 'forward_cache_manifest.json'):
        _need(_hash_record(a['provenance']['inputs'], name) ==
              _hash_record(b['provenance']['inputs'], name),
              f'Cannot pair different cache/exposure inputs: {name}')
    for name in ('forward_xsec_statistics.py', 'forward_xsec_problem.py'):
        _need(_hash_record(a['provenance']['sources'], name) ==
              _hash_record(b['provenance']['sources'], name),
              f'Cannot pair different resampling source: {name}')
    _need(a['bootstrap']['seed'] == b['bootstrap']['seed'] and
          a['provenance']['numpy'] == b['provenance']['numpy'],
          'Pairing requires identical bootstrap seed and NumPy RNG version')
    transforms = [np.vstack([functional(r['blocks'], len(r['parameters']), o) for o in observables])
                  for r in (a, b)]
    common = sorted(set(a['indices']) & set(b['indices']))
    _need(len(common) >= 2, 'At least two shared successful bootstrap replicates are required')
    central, paired = [], []
    for run, transform in zip((a, b), transforms):
        lookup = {replicate: i for i, replicate in enumerate(run['indices'])}
        paired.append(run['samples'][[lookup[i] for i in common]] @ transform.T)
        central.append(transform @ run['parameters'])
    centered = [sample - sample.mean(axis=0) for sample in paired]
    difference = paired[1] - paired[0]
    dc = difference - difference.mean(axis=0)
    denominator = len(common) - 1
    return {
        'schema_version': 1, 'publication_ready': False,
        'observables': observables, 'integral_unit': 'nb', 'covariance_unit': 'nb^2',
        'difference_convention': 'B minus A', 'central_A': central[0].tolist(),
        'central_B': central[1].tolist(), 'central_difference': (central[1]-central[0]).tolist(),
        'paired_difference_covariance': (dc.T @ dc / denominator).tolist(),
        'paired_crosscovariance_A_B': (centered[0].T @ centered[1] / denominator).tolist(),
        'paired_covariance_A': (centered[0].T @ centered[0] / denominator).tolist(),
        'paired_covariance_B': (centered[1].T @ centered[1] / denominator).tolist(),
        'paired_replicate_indices': common, 'paired_count': len(common),
        'replicates': {label: {
            'requested': run['bootstrap']['requested'], 'successful': len(run['indices']),
            'failed_indices': sorted(set(range(run['bootstrap']['requested'])) - set(run['indices'])),
            'successful_but_unpaired_indices': sorted(set(run['indices'])-set(common)),
            'reported_failures': run['bootstrap'].get('failures', []),
        } for label, run in zip(('A', 'B'), (a, b))},
        'input_directories': [str(Path(p).resolve()) for p in (run_a, run_b)],
        'input_sha256': {'A': a['files'], 'B': b['files']},
        'limitations': [
            'Conditional paired bootstrap spreads are not coverage-calibrated confidence intervals.',
            'Only common successful refits enter covariance; excluded failures can bias these spreads.',
            'Shared target, background, detector and radiative uncertainties are not included.',
            'Aligned bin integrals retain each fit\'s within-bin shape assumptions.',
            'No independent-fit assumption or chi-square p-value is used.'],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('A', type=Path)
    parser.add_argument('B', type=Path)
    parser.add_argument('--observables', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    temporary = None
    try:
        _need(not args.output.exists(), 'Output exists; choose a new JSON path')
        result = compare(args.A, args.B, json.loads(args.observables.read_text()))
        payload = json.dumps(result, indent=2, allow_nan=False) + '\n'
        with tempfile.NamedTemporaryFile(mode='w', dir=args.output.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, args.output)  # Atomic publication, refusing existing output.
    except (ValueError, KeyError, OSError) as error:
        parser.error(str(error))
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    print(f'Paired conditional comparison written: {args.output}')


if __name__ == '__main__':
    main()
