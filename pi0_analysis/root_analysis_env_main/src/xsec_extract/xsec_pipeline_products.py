#!/usr/bin/env python3
"""Production preflight and fresh-artifact verification; never repairs inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import time


def preflight(data, sim, raw, kin):
    import uproot
    for path, producer in ((data, 'analysis/combine'), (sim, 'simulation/smearing')):
        if not path.is_file():
            raise ValueError(f'kin={kin}: missing {path}; expected producer: {producer}')
    if not raw or not Path(raw).exists():
        raise ValueError(f'kin={kin}: missing raw SIMC file/directory {raw!r}; '
                         'supply --vertex_simc_file from exclusive SIMC production')
    with uproot.open(data) as source:
        if 'analysis_runs' not in source:
            raise ValueError(f'kin={kin}: {data} lacks analysis_runs; run validated per-run '
                             'analysis then src/analysis/combine_analysis_branches.py. '
                             'No extraction was attempted.')
        if source['analysis_runs'].num_entries == 0:
            raise ValueError(f'kin={kin}: empty analysis_runs in {data}; regenerate analysis/combine')
    # Detailed exposure, status, and matching validation remains authoritative
    # inside the extractor; do not implement a weaker alternative here.


def verify(args):
    import uproot
    required = [Path(x) for x in args.files]
    out = args.out
    if args.prepare:
        required += [out / x for x in ('data_events.csv', 'mc_events.csv',
                                      'forward_cache_manifest.json')]
    else:
        required += [out / x for x in ('migration_covariance.csv', 'fit_status.csv')]
        if args.mode == 'simc-model':
            required += [out / x for x in (
                'model_fit.root', 'model_parameters.csv', 'model_fit_status.csv',
                'model_covariance.csv', 'model_structure_functions.csv',
                'model_structure_covariance.csv', 'model_starts.csv',
                'model_structure_curves.csv', 'model_reconstructed_yields.csv',
                'model_vertex_matching.csv', 'model_vertex_rejections.csv',
                'model_vertex_epsilon.root',
                'model_prediction_validity.csv', 'model_identifiability.csv',
                'model_strong_correlations.csv', 'model_plot_manifest.json')]
            plot_manifest = out / 'model_plot_manifest.json'
            if plot_manifest.is_file():
                required += [out / p for p in json.loads(plot_manifest.read_text())['artifacts']]
    records = []
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'missing/empty required artifact: {path}')
        if path.stat().st_mtime_ns < args.started:
            raise ValueError(f'stale artifact from a previous invocation: {path}')
        if path.suffix == '.root':
            with uproot.open(path) as f:
                if not f.keys():
                    raise ValueError(f'empty ROOT directory: {path}')
        if path.suffix == '.pdf' and path.read_bytes()[:5] != b'%PDF-':
            raise ValueError(f'invalid PDF artifact: {path}')
        records.append({'path': str(path), 'bytes': path.stat().st_size,
                        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    result = {'schema_version': 1, 'mode': args.mode, 'kin': args.kin,
              'config': args.config, 'started_ns': args.started,
              'verified_ns': time.time_ns(), 'artifacts': records}
    tmp = out / 'pipeline_artifacts.json.partial'
    tmp.write_text(json.dumps(result, indent=2) + '\n')
    tmp.replace(out / 'pipeline_artifacts.json')
    print(f'[verify] {len(records)} fresh artifacts verified')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    pre = sub.add_parser('preflight')
    pre.add_argument('--data', type=Path, required=True)
    pre.add_argument('--sim', type=Path, required=True)
    pre.add_argument('--raw', required=True)
    pre.add_argument('--kin', required=True)
    check = sub.add_parser('verify')
    check.add_argument('--out', type=Path, required=True)
    check.add_argument('--started', type=int, required=True)
    check.add_argument('--mode', required=True)
    check.add_argument('--kin', required=True)
    check.add_argument('--config', required=True)
    check.add_argument('--prepare', action='store_true')
    check.add_argument('files', nargs='*')
    args = parser.parse_args()
    try:
        if args.command == 'preflight':
            preflight(args.data, args.sim, args.raw, args.kin)
        else:
            verify(args)
    except (ValueError, OSError, KeyError) as error:
        parser.exit(1, f'[ERROR] {error}\n')


if __name__ == '__main__':
    main()
