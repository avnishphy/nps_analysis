#!/usr/bin/env python3
"""Synthetic joint-fit regression; requires NumPy, uproot and the Hall C ROOT env.

Run from the repository: csh -c 'source /usr/share/Modules/init/csh;
source /group/nps/singhav/setup.csh; python3 tests/test_joint_xsec_fit.py'
Only temporary files are written. Closure verifies fit algebra, not physics bias.
"""
import csv
import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

import numpy as np

REPO = Path(__file__).resolve().parents[1]
SOURCE = REPO / "src/xsec_extract"
sys.path.insert(0, str(SOURCE))
import run_joint_xsec_fit as driver


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def fixture(nsettings=2, guard=False, noisy=False, physical="interior"):
    nphi, nb = 24, 7
    bins = {"tprime_bin_edges": [0, 1], "q2_bin_edges": [1, 2],
            "xb_bin_edges_by_q2": [[.2, .4]],
            "phi_bin_edges": np.linspace(0, 2*np.pi, nphi+1).tolist()}
    expected = {(0, s, "U"): 10. + 4*s for s in range(nsettings)}
    expected.update({(0, -1, "LT"): -.7, (0, -1, "TT"): .4})
    if physical == "envelope":
        expected.update({(0, 0, "U"): .6, (0, 1, "U"): 5.,
                         (0, -1, "LT"): .5, (0, -1, "TT"): .1})
    elif physical == "boundary":
        expected.update({(0, 0, "U"): .05, (0, 1, "U"): 2.,
                         (0, -1, "LT"): .9, (0, -1, "TT"): .8})
    if guard:
        expected.update({(1, 0, "U"): 7., (1, -1, "LT"): .3,
                         (1, -1, "TT"): -.2})
    settings = []
    for s in range(nsettings):
        e = (.25 if s == 0 else .8) if physical != "interior" else .6
        setting = {"reco": [], "cells": [], "epsilon": [e] + [0.] * 6}
        if guard and s == 0:
            setting["epsilon"][1] = e
        for r in range(nphi):
            phi = (r+.5)*2*np.pi/nphi
            data = 0.
            for b in range(nb):
                basis = np.zeros(3)
                if b == 0:
                    acceptance = (1+.2*np.sin(phi)+.1*np.cos(3*phi))*(s+1)
                    basis = acceptance*np.array([1, np.sqrt(2*e*(1+e))*np.cos(phi),
                                                 e*np.cos(2*phi)])
                elif guard and b == 1 and s == 0:
                    basis = np.array([.3+.1*np.cos(3*phi), .1*np.sin(phi), .1*np.sin(2*phi)])
                if np.any(basis):
                    p = np.array([expected[b, s, "U"], expected[b, -1, "LT"],
                                  expected[b, -1, "TT"]])
                    data += basis @ p
                # Poissonized event moments with nonzero cross covariances.
                cov = .002*np.outer(basis, basis)
                cell = dict(zip(driver.BASIS + driver.COV,
                                map(lambda x: format(x, '.17g'), [*basis, *cov.ravel()])))
                setting["cells"].append(cell)
            if noisy:
                data += .2*np.sin(5*phi+.2*s)+.07*np.cos(phi)
            setting["reco"].append({"data": format(data, '.17g'),
                                    "data_variance": str(.5+.03*r+.1*s)})
        settings.append(setting)
    return settings, bins, expected


def check_case(solver, temp, name, *, finite=False, positive=False, **options):
    settings, bins, expected = fixture(**options)
    problem, result = temp / (name + '.txt'), temp / name
    result.mkdir()
    driver.write_problem(problem, settings, bins)
    subprocess.run([str(solver), str(problem), str(result),
                    "finite-mc" if finite else "data", str(int(positive)),
                    "1e-10", "100", "1e-10"], check=True)
    records = read_csv(result / 'joint_parameters.csv')
    keys = [(int(r['truth_block']), int(r['setting_index']), r['component']) for r in records]
    assert set(keys) == set(expected), (name, keys)
    p = np.array([float(r['value']) for r in records])
    n = len(p)
    x, y, var, row_cov = [], [], [], []
    for s, setting in enumerate(settings):
        for r, reco in enumerate(setting['reco']):
            design, cov = np.zeros(n), np.zeros((n, n))
            for b in range(7):
                if (b, s, 'U') not in keys:
                    continue
                index = [keys.index((b, s if c == 'U' else -1, c)) for c in driver.COMPONENTS]
                cell = setting['cells'][7*r+b]
                design[index] = [float(cell[k]) for k in driver.BASIS]
                cov[np.ix_(index, index)] = np.array([float(cell[k]) for k in driver.COV]).reshape(3, 3)
            x.append(design)
            y.append(float(reco['data']))
            var.append(float(reco['data_variance']))
            row_cov.append(cov)
    x, y, var, row_cov = map(np.array, (x, y, var, row_cov))
    exported = read_csv(result / 'joint_rows.csv')
    np.testing.assert_allclose([float(r['prediction']) for r in exported], x @ p, atol=1e-11)
    np.testing.assert_allclose([float(r['mc_variance_final']) for r in exported],
                               np.einsum('i,rij,j->r', p, row_cov, p), atol=1e-11)
    variance = np.array([float(r['variance_used']) for r in exported])
    summary = dict(line.split('=', 1) for line in (result/'joint_summary.txt').read_text().splitlines())
    assert int(summary['ndf']) == len(y)-n
    assert int(summary['rank']) == n
    np.testing.assert_allclose(float(summary['chi2']), np.sum((y-x@p)**2/variance), atol=1e-10)
    boundary = options.get('physical') == 'boundary'
    assert int(summary['positivity_boundary_active']) == int(boundary)
    if not boundary:
        reference_variance = var.copy()
        for _ in range(100):
            w = x / np.sqrt(reference_variance[:, None])
            reference = np.linalg.lstsq(w, y/np.sqrt(reference_variance), rcond=None)[0]
            updated = var + np.einsum('i,rij,j->r', reference, row_cov, reference) if finite else var
            if np.max(abs(updated-reference_variance)/updated) < 1e-12:
                break
            reference_variance = updated
        else:
            raise AssertionError('independent MC iteration failed')
        np.testing.assert_allclose(p, reference, rtol=1e-8, atol=1e-10)
        covariance = np.array([float(r['stat_plus_mc_covariance'])
                               for r in read_csv(result/'joint_covariance.csv')]).reshape(n, n)
        np.testing.assert_allclose(covariance, np.linalg.inv(w.T @ w), rtol=1e-8, atol=1e-10)
        if not options.get('noisy'):
            np.testing.assert_allclose(p, [expected[k] for k in keys], atol=1e-10)
    else:
        assert all(np.isnan(float(r['error_stat_plus_mc'])) for r in records)
        assert float(summary['chi2']) > 0
    if positive:
        diagnostics = read_csv(result/'joint_positivity.csv')
        assert len(diagnostics) == sum(c == 'U' for _, _, c in keys)
        for row in diagnostics:
            b, s = int(row['truth_block']), int(row['setting_index'])
            e = settings[s]['epsilon'][b]
            assert float(row['epsilon_max']) == e
            u, lt, tt = [p[keys.index((b, s if c == 'U' else -1, c))] for c in driver.COMPONENTS]
            z = np.linspace(-1, 1, 10001)
            angular = u+np.sqrt(2*e*(1+e))*lt*z+e*tt*(2*z*z-1)
            assert angular.min() >= -1e-8
            assert float(row['minimum_response_bracket']) >= -float(row['feasibility_tolerance'])
    print('PASS', name)


def main():
    with tempfile.TemporaryDirectory(prefix='test_joint_xsec_') as folder:
        temp = Path(folder)
        # Both old three-component and newer four-component exports must
        # preserve the same U/LT/TT moments, including all cross terms.
        direct = temp/'direct'
        four = temp/'four'
        direct.mkdir()
        four.mkdir()
        basis = np.array([2., 1.1, -.3, .4])
        covariance = np.outer(basis, basis)
        components = ('T', 'L', 'LT', 'TT')
        record = {'reco_row': '0', 'truth_block': '0', 'events': '1'}
        record.update(zip(('basis_'+c for c in components), map(str, basis)))
        record.update(zip(('cov_'+a+'_'+b for a in components for b in components),
                          map(str, covariance.ravel())))
        path = four/'migration_joint_response_cells.csv'
        with path.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(record))
            writer.writeheader()
            writer.writerow(record)
        index = [0, 2, 3]
        expected = dict(zip(driver.BASIS + driver.COV,
                            map(str, [*basis[index], *covariance[np.ix_(index, index)].ravel()])))
        with (direct/'migration_response_cells.csv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(expected))
            writer.writeheader()
            writer.writerow(expected)
        for directory in (direct, four):
            cells, sources, _ = driver.load_response_cells(directory)
            assert all(cells[0][key] == value for key, value in expected.items())
            assert len(sources) == 1 and sources[0].is_file()
        print('PASS direct and four-component response inputs')
        preset = json.loads((SOURCE/'xsec_config/xsec_config_x36_4.json').read_text())
        (temp/'xsec_config.h').write_text(driver.render(preset))
        solver = temp/'solver'
        flags = shlex.split(subprocess.check_output(['root-config', '--cflags', '--libs'], text=True))
        subprocess.run(['g++', '-std=c++17', '-O2', '-I'+str(temp), '-I'+str(SOURCE),
                        str(SOURCE/'xsec_joint_solver.C'), *flags, '-o', str(solver)], check=True)
        check_case(solver, temp, 'two_settings_equal_epsilon')
        check_case(solver, temp, 'three_settings_guard_finite_mc', nsettings=3,
                   guard=True, noisy=True, finite=True, positive=True)
        check_case(solver, temp, 'per_setting_epsilon', physical='envelope', positive=True)
        check_case(solver, temp, 'positivity_boundary', physical='boundary', positive=True)
        check_case(solver, temp, 'positivity_boundary_finite_mc', physical='boundary',
                   positive=True, finite=True)
        settings, bins, _ = fixture()
        for setting in settings:
            for cell in setting['cells']:
                cell['basis_LT'] = '0'
        problem = temp/'singular.txt'
        driver.write_problem(problem, settings, bins)
        failed = subprocess.run([str(solver), str(problem), str(temp), 'data', '0',
                                 '1e-10', '30', '1e-6'], capture_output=True, text=True)
        assert failed.returncode != 0 and 'unsupported truth coefficient' in failed.stderr, failed.stderr
        print('PASS rejects unsupported shared LT')


if __name__ == '__main__':
    main()
