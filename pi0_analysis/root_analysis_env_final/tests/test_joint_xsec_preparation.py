#!/usr/bin/env python3
"""End-to-end raw ROOT -> prepared response -> joint fit, with temporary files.

Run with python3 tests/test_joint_xsec_preparation.py. The driver loads Hall C
itself. Requires NumPy/uproot and access to /group/nps/singhav/setup.csh.
"""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src/xsec_extract'))
import run_joint_xsec_fit as driver


def tree(path, name, columns):
    with uproot.recreate(path) as output:
        output.mktree(name, {key: value.dtype for key, value in columns.items()})
        output[name].extend(columns)


def fixture(folder, index):
    folder.mkdir()
    preset = 'xsec_config_x36_4.json' if index == 0 else 'xsec_config_x36_5_407.json'
    config = json.loads((driver.HERE/'xsec_config'/preset).read_text())
    q, mp, mpi = 5., .9382720813, .1349768
    w = np.float32(np.sqrt(mp*mp+q*(1/.4-1)))
    q0 = (float(w)**2-mp*mp-q)/(2*float(w))
    epi = (float(w)**2+mpi*mpi-mp*mp)/(2*float(w))
    tfwd = mpi*mpi-q-2*(q0*epi-np.sqrt(q0*q0+q)*np.sqrt(epi*epi-mpi*mpi))
    theta = np.deg2rad(config['hms_theta_deg'])
    nu = (float(w)**2+q-mp*mp)/(2*mp)
    e = 1/(1+2*(1+nu*nu/q)*(1-np.cos(theta))/(1+np.cos(theta)))
    phi = np.array([(i+.5)*2*np.pi/12 for i in range(12)], dtype='float32')
    u, lt, tt = 6.+4*index, .3, -.2
    y = (u+np.sqrt(2*e*(1+e))*lt*np.cos(phi.astype(float))+
         e*tt*np.cos(2*phi.astype(float)))/(2*np.pi)
    def const(value, dtype='float64'):
        return np.full(12, value, dtype=dtype)
    data = folder/'combined_branches_LH2.root'
    sim, vertex = folder/'sim.root', folder/'vertex.root'
    tree(data, 'physics', {'Q2': const(q), 't': const(-.4), 'tmin': const(-.3),
         'xB': const(.4), 'phi': phi.astype('float64'), 'W': const(w),
         'mmiss_all': const(.95), 'pi0_weight': y*config['tgt_contam'],
         'scale': const(1, 'float32'), 'charge_uC': const(1000, 'float32'),
         'run_number': const(index+1, 'int32')})
    tree(sim, 'simulation', {'Q2': const(q, 'float32'), 't': const(-.4, 'float32'),
         'tmin': const(-.3, 'float32'), 'xB': const(.4, 'float32'), 'phi': phi,
         'W': const(w, 'float32'), 'mmiss': const(.95, 'float32'),
         'full_weight': const(1, 'float32'), 'sigcm': const(1, 'float32'),
         'is_exclusive': const(1, 'int32'), 'event_id': np.arange(12, dtype='uint64')})
    tree(vertex, 'h10', {'Q2i': const(q, 'float32'), 'Wi': const(w, 'float32'),
         'ti': const(-(tfwd-.1), 'float32'), 'phipqi': phi,
         'sigcm': const(1, 'float32'), 'hsxptari': const(0, 'float32'),
         'hsyptari': const(0, 'float32')})
    config.update(configured_kinematic=f'KinC_test{index}', data_file=str(data), simc_file=str(sim),
                  tprime_bin_edges=[-.5-.1*index, 0.], t_bin_edges=[-1., 0.],
                  q2_bin_edges=[4., 6.], xb_bin_edges_by_q2=[[.3, .5]],
                  diamond_xb_q2_vertices=[[.32, 5.], [.4, 4.1], [.48, 5.], [.4, 5.9]],
                  phi_bins=12, mmiss_select='window')
    config.pop('phi_bin_edges', None)
    path = folder/'config.json'
    path.write_text(json.dumps(config))
    return path, vertex, y


def main():
    with tempfile.TemporaryDirectory(prefix='joint_preparation_test_') as temporary:
        folder = Path(temporary)
        fixtures = [fixture(folder/str(s), s) for s in range(2)]
        original_hashes = [driver.digest(p) for p, _, _ in fixtures]
        output, inputs = folder/'fit', folder/'prepared'
        command = [sys.executable, str(driver.HERE/'run_joint_xsec_fit.py')]
        for config, vertex, _ in fixtures:
            command += ['--prepare-setting', str(config), str(vertex)]
        command += ['--binning-config', str(fixtures[0][0]), '--out-dir', str(output),
                    '--inputs-dir', str(inputs), '--fit-variance', 'data', '--positive-xsec', '--no-plots']
        subprocess.run(command, check=True)
        assert [driver.digest(p) for p, _, _ in fixtures] == original_hashes
        manifest = json.loads((output/'joint_manifest.json').read_text())
        assert manifest['options']['fit_objective'] == 'gaussian'
        assert all(s['input_stage'] == 'prepared_joint_inputs_v1' for s in manifest['settings'])
        expected_diamond = json.loads(fixtures[0][0].read_text())['diamond_xb_q2_vertices']
        assert driver.same(manifest['bins']['diamond_xb_q2_vertices'],
                           driver.diamond_vertices(expected_diamond))
        parameters = driver.rows(output/'joint_parameters.csv')
        actual = {(int(r['setting_index']), r['component']): float(r['value']) for r in parameters}
        for key, expected in {(0, 'U'): 6., (1, 'U'): 10., (-1, 'LT'): .3, (-1, 'TT'): -.2}.items():
            np.testing.assert_allclose(actual[key], expected, rtol=1e-7, atol=1e-8)
        eps = [s['epsilon'] for s in manifest['lt_separation']['settings']]
        separated = driver.rows(output/'joint_separated_parameters.csv')
        expected_l = 4./(eps[1]-eps[0])
        np.testing.assert_allclose([float(r['value']) for r in separated],
                                   [6.-eps[0]*expected_l, expected_l, .3, -.2], rtol=1e-7, atol=1e-8)
        assert all(r['status'] == 'ok' for r in manifest['lt_separation']['blocks'])
        reuse = [sys.executable, str(driver.HERE/'run_joint_xsec_fit.py')]
        for s, (_, _, y) in enumerate(fixtures):
            setting = inputs/f'KinC_test{s}'
            assert not list(setting.glob('*.root'))
            assert not (setting/'migration_parameters.csv').exists()
            assert not (setting/'positivity_diagnostics.csv').exists()
            rows = driver.rows(setting/'migration_reco_rows.csv')
            np.testing.assert_allclose([float(r['data']) for r in rows], y, rtol=1e-12)
            np.testing.assert_allclose([float(r['data_variance']) for r in rows], y*y, rtol=1e-12)
            assert driver.read_metadata(setting/'joint_input_metadata.txt')['fit_objective'] == 'not_run'
            config = inputs/'configs'/f'KinC_test{s}.json'
            reuse += ['--setting', str(config), str(setting)]
        # A second joint fit changes variance/positivity without repeating
        # event preparation or satisfying individual-fit flags.
        reuse += ['--out-dir', str(folder/'refit'), '--fit-variance', 'finite-mc', '--no-plots',
                  '--nominal-epsilon', '.25', '--nominal-epsilon', '.75']
        subprocess.run(reuse, check=True)
        separated = driver.rows(folder/'refit'/'joint_separated_parameters.csv')
        np.testing.assert_allclose([float(r['value']) for r in separated], [4., 8., .3, -.2], rtol=1e-7)
        rejected = subprocess.run(command, capture_output=True, text=True)
        assert rejected.returncode != 0 and 'already exists' in rejected.stderr
        print('PASS raw preparation, common binning, normalization, joint closure, reuse, and overwrite refusal')


if __name__ == '__main__':
    main()
