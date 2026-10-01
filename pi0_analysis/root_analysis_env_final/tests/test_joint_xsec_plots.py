#!/usr/bin/env python3
"""Compare native plot coverage with joint rendering on synthetic ROOT inputs.

python3 tests/test_joint_xsec_plots.py [--partons] [--keep-dir /tmp/new-directory]
PARTONS uses small smoke-test integration counts, not production precision.
"""
import argparse
import json
import copy
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import uproot

from test_joint_xsec_preparation import fixture, driver, tree
import joint_xsec_plots


def ellipse_fixture(config_path):
    config = json.loads(config_path.read_text())
    config['mmiss_select'] = 'ellipse'
    config_path.write_text(json.dumps(config))
    for path, name, fields in [
            (Path(config['data_file']), 'physics', {'mpi0_all': .135, 'mmiss_all_corr': .95,
                                                  'is_exclusive_ellipse_combined': 1}),
            (Path(config['simc_file']), 'simulation', {'mpi0': .135})]:
        with uproot.open(path) as source:
            columns = source[name].arrays(library='np')
        for key, value in fields.items():
            dtype = 'int32' if key.startswith('is_') else ('float64' if name == 'physics' else 'float32')
            columns[key] = np.full(12, value, dtype=dtype)
        tree(path, name, columns)
    data = Path(config['data_file'])
    metadata = dict(tag='combined_2d_mass_cut', ellipse_valid=1, mean_mpi0=.135,
                    mean_mmiss=.95, cov_mpi0_mpi0=.0001, cov_mpi0_mmiss=0,
                    cov_mmiss_mmiss=.01, ellipse_d2_cut=4., mpi0_min=.1, mpi0_max=.17,
                    mmiss_min=.5, mmiss_max=1.5)
    data.with_name(data.stem+'_combined_2d_mass_cut_debug.txt').write_text(
        ''.join(f'{k}={v}\n' for k, v in metadata.items()))


def check(folder, partons):
    reference_hash = driver.digest(driver.HERE/'run_xsec_pipeline.sh')
    fixtures = [fixture(folder/str(s), s) for s in range(2)]
    for config, _, _ in fixtures:
        ellipse_fixture(config)
    output = folder/'joint'
    command = [sys.executable, str(driver.HERE/'run_joint_xsec_fit.py')]
    for config, vertex, _ in fixtures:
        command += ['--prepare-setting', str(config), str(vertex)]
    command += ['--binning-config', str(fixtures[0][0]), '--out-dir', str(output),
                '--fit-variance', 'finite-mc', '--positive-xsec']
    if partons:
        command += ['--partons', '--partons-warmups', '100', '--partons-calls', '1000']
    with (folder/'joint.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    report = json.loads((output/'plots/plot_manifest.json').read_text())
    assert report['status'] == 'complete' and not report['individual_fits_run']
    assert report['partons'] == partons
    original = driver.rows(output/'joint_parameters.csv')
    for s in range(2):
        setting = output/'plots/settings'/f'KinC_test{s}'
        log = (setting/'render.log').read_text()
        assert 'no individual fit or refit toys were run' in log
        assert 'Global fit attempt' not in log
        fitted = driver.rows(setting/'joint_plot_slices.csv')[0]
        for component, field in [('U', 'sigmaU'), ('LT', 'sigmaTL'), ('TT', 'sigmaTT')]:
            row = next(r for r in original if r['component'] == component and
                       int(r['setting_index']) == (s if component == 'U' else -1))
            np.testing.assert_allclose(float(fitted['fit_xsec_'+field]), float(row['value']), rtol=1e-12)
            np.testing.assert_allclose(float(fitted['fit_xsec_'+field+'err']), float(row['error_stat_plus_mc']), rtol=1e-12)
        assert fitted['fit_scope'] == 'joint_shared_LT_TT_independent_U'
        points = driver.rows(setting/'experimental_points.csv')
        assert all(r['status'] == 'central_only_joint' for r in points)
        assert all(np.isnan(float(r['sigma_exp_error'])) for r in points)
    baseline = folder/'individual_reference'
    # Run the same native extractor as the reference pipeline. Its shell file
    # and presets remain read-only; this baseline uses private test configs.
    config = json.loads(fixtures[0][0].read_text())
    joint_xsec_plots.run_native(config, fixtures[0][1], baseline, driver.root_environment(),
                               partons, 100, 1000)
    reference = {str(p.relative_to(baseline)) for p in baseline.rglob('*.pdf')}
    native = output/'plots/settings/KinC_test0'
    actual = {str(p.relative_to(native)) for p in native.rglob('*.pdf')}
    assert reference <= actual, ('missing native pages', reference-actual)
    assert len(reference) >= 14, reference
    extras = output/'plots/joint'
    assert {p.stem for p in extras.glob('*.pdf')} >= {
        'rosenbluth_truth_0', 'coefficient_ellipses_truth_0', 'separated_terms_q0_x0',
        'per_setting_U_q0_x0', 'joint_correlation', 'separated_correlation',
        'joint_fit_quality', 'positivity_and_epsilon'}
    for relative in report['pages']:
        pdf = output/relative
        assert pdf.is_file() and pdf.stat().st_size > 1000
        assert pdf.with_suffix('.png').is_file()
    info = subprocess.check_output(['pdfinfo', str(output/report['combined_pdf'])], text=True)
    pages = int(next(line.split(':')[1] for line in info.splitlines() if line.startswith('Pages:')))
    assert pages == report['page_count'] == len(report['pages'])
    assert driver.digest(driver.HERE/'run_xsec_pipeline.sh') == reference_hash
    edge_figures(output)
    print(f'PASS {len(reference)} native PDF products covered, {len(list(extras.glob("*.pdf")))} additional joint figures, {pages} combined pages')


def edge_figures(output):
    """Unavailable uncertainties/epsilon separation still produce honest plots."""
    manifest = json.loads((output/'joint_manifest.json').read_text())
    for case in ('boundary', 'degenerate_epsilon'):
        with tempfile.TemporaryDirectory(prefix='joint_plot_edge_') as temporary:
            destination = Path(temporary)
            for name in ('joint_parameters.csv', 'joint_covariance.csv', 'joint_rows.csv', 'joint_positivity.csv'):
                shutil.copy2(output/name, destination/name)
            current = copy.deepcopy(manifest)
            if case == 'boundary':
                covariance = driver.rows(destination/'joint_covariance.csv')
                for row in covariance:
                    row['stat_plus_mc_covariance'] = 'nan'
                driver.write_csv(destination/'joint_covariance.csv', tuple(covariance[0]), covariance)
            else:
                for setting in current['lt_separation']['settings']:
                    setting['epsilon'] = .5
            current['lt_separation'] = driver.separate_lt(
                destination, current['lt_separation']['settings'], 1e-10)
            pages = joint_xsec_plots.joint_figures(destination, current)
            assert len(pages) == 8 and all(p.is_file() for p in pages)
    protected = output/'plots/plot_manifest.json'
    before = driver.digest(protected)
    try:
        joint_xsec_plots.render(output)
    except ValueError as error:
        assert 'already complete' in str(error)
    else:
        raise AssertionError('completed report was overwritten')
    assert driver.digest(protected) == before
    print('PASS boundary/degenerate plots and completed-report overwrite refusal')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--partons', action='store_true')
    parser.add_argument('--keep-dir', type=Path)
    args = parser.parse_args()
    if args.keep_dir:
        folder = args.keep_dir.resolve()
        folder.mkdir(parents=True, exist_ok=False)
        check(folder, args.partons)
    else:
        with tempfile.TemporaryDirectory(prefix='joint_plot_test_') as temporary:
            check(Path(temporary), args.partons)
