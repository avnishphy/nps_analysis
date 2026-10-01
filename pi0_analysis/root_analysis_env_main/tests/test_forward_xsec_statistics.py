"""Independent toy checks for the interim forward-fit statistical machinery.

These checks exercise the specified conditional model. They do not demonstrate
coverage for weighted backgrounds, imperfect detector response, or real data.
"""

import importlib.util
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


MODULE = Path(__file__).resolve().parents[1] / 'src/xsec_extract/forward_xsec_statistics.py'
SPEC = importlib.util.spec_from_file_location('forward_xsec_statistics', MODULE)
stats = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stats)


def toy_problem(seed=7103, *, sample=True, weights=False, exposure=4.0):
    """Two truth cells with migration, independently integrated harmonic response."""
    rng = np.random.default_rng(seed)
    nphi, blocks, epsilon = 12, 2, 0.72
    # Analytic bin integrals, independent of the fitting positivity evaluator.
    edges = np.linspace(0.0, 2.0 * np.pi, nphi + 1)
    width = np.diff(edges)
    c1 = (np.sin(edges[1:]) - np.sin(edges[:-1])) / width
    c2 = (np.sin(2 * edges[1:]) - np.sin(2 * edges[:-1])) / (2 * width)
    angular = np.column_stack((np.ones(nphi),
        np.sqrt(2 * epsilon * (1 + epsilon)) * c1, epsilon * c2))
    migration = np.array([[0.8, 0.3], [0.2, 0.7]])
    design = np.zeros((blocks * nphi, blocks * 3))
    mc_rows, mc_blocks, mc_basis, mc_ids = [], [], [], []
    # Response copies share an original ID across reconstructed rows. Testing
    # this grouping catches independently fluctuating fragments of one event.
    events_per_cell = 40
    for block in range(blocks):
        for phi in range(nphi):
            for reco in range(blocks):
                row = reco * nphi + phi
                contribution = exposure * migration[reco, block] * angular[phi]
                design[row, 3 * block:3 * block + 3] += contribution
                for event in range(events_per_cell):
                    mc_rows.append(row)
                    mc_blocks.append(block)
                    mc_basis.append(contribution / events_per_cell)
                    mc_ids.append((block * nphi + phi) * events_per_cell + event)
    truth = np.array([28.0, -4.0, -5.0, 19.0, 3.0, 2.0])
    expected = design @ truth
    event_weights = np.linspace(0.8, 1.2, len(expected)) if weights else np.ones(len(expected))
    counts = rng.poisson(expected / event_weights) if sample else expected / event_weights
    y, sumw2 = counts * event_weights, counts * event_weights**2
    problem = {'design': design, 'y': y, 'sumw2': sumw2,
               'epsilon_max': np.full(blocks, epsilon), 'truth': truth,
               'mc_rows': np.array(mc_rows), 'mc_blocks': np.array(mc_blocks),
               'mc_basis': np.array(mc_basis), 'mc_ids': np.array(mc_ids),
               'published_blocks': [0], 'truth_blocks': [{'kind': 'interior'}, {'kind': 'guard'}]}
    if sample:
        problem['data_rows'] = np.repeat(np.arange(len(y)), counts)
        problem['data_weights'] = np.repeat(event_weights, counts)
        problem['data_ids'] = np.arange(len(problem['data_rows']))
    return problem


class ForwardStatisticsTests(unittest.TestCase):
    def test_representable_migrating_truth_recovery(self):
        problem = toy_problem(sample=False)
        fit = stats.fit_problem(problem)
        np.testing.assert_allclose(fit['parameters'], problem['truth'], atol=2e-4)
        self.assertEqual(fit['parameter_count'], 6)  # Includes guard coefficients.
        self.assertLess(abs(fit['deviance']), 1e-7)
        self.assertFalse(fit['publication_ready'])

    def test_yield_unit_invariance(self):
        problem = toy_problem(weights=True)
        nominal = stats.fit_problem(problem)
        for scale in (1e-6, 1e6):
            other = dict(problem)
            for name in ('design', 'y', 'mc_basis', 'data_weights'):
                other[name] = problem[name] * scale
            other['sumw2'] = problem['sumw2'] * scale**2
            fit = stats.fit_problem(other)
            np.testing.assert_allclose(fit['parameters'], nominal['parameters'], rtol=1e-7, atol=1e-6)
            self.assertAlmostEqual(fit['deviance'], nominal['deviance'], places=7)
            stats.resample_problem(other, np.random.default_rng(521))

    def test_empty_rows_retained_with_reported_weight_scale(self):
        problem = toy_problem(weights=True)
        problem['y'][3] = problem['sumw2'][3] = 0.0
        fit = stats.fit_problem(problem)
        self.assertEqual(fit['included_rows'], 24)
        self.assertIn(3, fit['empty_rows_borrowing_pooled_scale'])
        self.assertAlmostEqual(fit['row_scales'][3], problem['sumw2'].sum() / problem['y'].sum())
        self.assertTrue(all(v >= -1e-6 for v in fit['angular_minima']))

    def test_exact_angular_minimum_between_grid_points(self):
        epsilon = 0.8
        triple = np.array([1.0, -0.95, 0.7])
        minimum, x = stats.angular_minimum(triple, epsilon)
        grid = np.linspace(-1, 1, 100001)
        reference = triple[0] + np.sqrt(2 * epsilon * (1 + epsilon)) * triple[1] * grid
        reference += epsilon * triple[2] * (2 * grid**2 - 1)
        self.assertAlmostEqual(minimum, reference.min(), places=8)
        self.assertTrue(-1 < x < 1)
        for eps in np.linspace(0, epsilon, 31):
            self.assertGreaterEqual(stats.angular_minimum(triple, eps)[0], minimum - 1e-12)

    def test_fit_boundary_and_signed_harmonics(self):
        problem = toy_problem(sample=False)
        # Strong angular contrast forces the optimizer onto a physical boundary.
        y = problem['y'].copy()
        y[[0, 1, 10, 11, 12, 13, 22, 23]] *= 0.015
        problem['y'] = problem['sumw2'] = y
        fit = stats.fit_problem(problem)
        self.assertTrue(fit['positivity_boundary_blocks'])
        self.assertTrue(all(v >= -1e-5 for v in fit['angular_minima']))
        self.assertTrue(any(v < 0 for v in fit['parameters']))

    def test_rank_failure_is_explicit(self):
        problem = toy_problem()
        problem['design'][:, 4] = problem['design'][:, 1]
        with self.assertRaisesRegex(stats.FitError, 'rank-deficient'):
            stats.fit_problem(problem)

    def test_negative_data_fails(self):
        problem = toy_problem()
        problem['y'][0] = -1
        with self.assertRaisesRegex(stats.FitError, 'nonnegative'):
            stats.fit_problem(problem)

    def test_cache_normalization_mismatch_fails(self):
        problem = toy_problem()
        problem['mc_basis'] *= 1.1
        with self.assertRaisesRegex(stats.FitError, 'normalization'):
            stats.resample_problem(problem, np.random.default_rng(8))

    def test_mc_original_event_grouping_and_response_rebuild(self):
        problem = toy_problem()
        replica = stats.resample_problem(problem, np.random.default_rng(19), resample_data=False)
        self.assertFalse(np.array_equal(replica['design'], problem['design']))
        # Shared generator records preserve the two-row migration ratio exactly.
        np.testing.assert_allclose(replica['design'][:12, :3] * 0.25,
                                   replica['design'][12:, :3], atol=1e-12)
        np.testing.assert_array_equal(replica['y'], problem['y'])

    def test_bootstrap_refits_response_and_all_nuisance_coefficients(self):
        problem = toy_problem(weights=True)
        result = stats.bootstrap_problem(problem, replicates=12, seed=739)
        self.assertEqual(result['successful'], 12, result['failures'])
        self.assertEqual(result['response_changed_replicates'], 12)
        self.assertEqual(np.asarray(result['samples']).shape, (12, 6))
        self.assertGreater(np.asarray(result['covariance'])[0, 0], 0)
        self.assertIn('NOT_COVERAGE_CALIBRATED', result['inference_status'])

    def test_complete_catalog_keeps_binning_resamples_paired(self):
        full = toy_problem()
        full['bootstrap_data_ids'] = full['data_ids'].copy()
        full['bootstrap_mc_ids'] = np.unique(full['mc_ids'])
        selected = dict(full)
        # Keep only first reconstructed region; selected arrays must still
        # reproduce that reduced design and yield before resampling.
        selected['design'] = full['design'][:12].copy()
        selected['y'] = full['y'][:12].copy()
        selected['sumw2'] = full['sumw2'][:12].copy()
        data_mask = full['data_rows'] < 12
        mc_mask = full['mc_rows'] < 12
        for name in ('data_rows', 'data_weights', 'data_ids'):
            selected[name] = full[name][data_mask]
        for name in ('mc_rows', 'mc_blocks', 'mc_basis', 'mc_ids'):
            selected[name] = full[name][mc_mask]
        left = stats.resample_problem(full, np.random.default_rng(234))
        right = stats.resample_problem(selected, np.random.default_rng(234))
        np.testing.assert_array_equal(left['y'][:12], right['y'])
        np.testing.assert_array_equal(left['sumw2'][:12], right['sumw2'])
        np.testing.assert_array_equal(left['design'][:12], right['design'])

    def test_repeated_data_records_share_weight_variance_and_multiplier(self):
        full = toy_problem()
        split = dict(full)
        split['data_rows'] = np.repeat(full['data_rows'], 2)
        split['data_weights'] = np.repeat(full['data_weights'] / 2, 2)
        split['data_ids'] = np.repeat(full['data_ids'], 2)
        left = stats.resample_problem(full, np.random.default_rng(932))
        right = stats.resample_problem(split, np.random.default_rng(932))
        np.testing.assert_array_equal(left['y'], right['y'])
        np.testing.assert_array_equal(left['sumw2'], right['sumw2'])
        np.testing.assert_array_equal(left['design'], right['design'])

    def test_zero_weight_data_event_is_permitted(self):
        problem = toy_problem()
        problem['data_rows'] = np.append(problem['data_rows'], 0)
        problem['data_weights'] = np.append(problem['data_weights'], 0.0)
        problem['data_ids'] = np.append(problem['data_ids'], -1)
        replica = stats.resample_problem(problem, np.random.default_rng(42))
        self.assertTrue(np.all(replica['y'] >= 0))

    def test_profile_refits_guard_terms_at_fixed_linear_integral(self):
        problem = toy_problem()
        nominal = stats.fit_problem(problem)
        functional = np.array([0.6, 0, 0, 0.4, 0, 0])
        value = functional @ nominal['parameters']
        profile = stats.profile_problem(problem, functional, [value - 1, value, value + 1],
                                        fit_result=nominal)
        self.assertTrue(all(p['status'] == 'converged' for p in profile['points']), profile)
        self.assertLess(profile['points'][1]['delta_deviance'], 1e-6)
        for point in profile['points']:
            self.assertAlmostEqual(functional @ point['parameters'], point['value'], places=6)
        self.assertIn('UNCALIBRATED', profile['inference_status'])

    def test_independent_conditional_poisson_ensemble(self):
        # A deliberately broad regression guard, not a precision coverage claim:
        # data come from independent Poisson draws with exact fixed response and
        # interior truth. The real-data bootstrap interval is never certified.
        covered, trials = 0, 24
        functional = np.array([1, 0, 0, 0, 0, 0])
        for seed in range(1200, 1200 + trials):
            problem = toy_problem(seed=seed)
            nominal = stats.fit_problem(problem)
            profile = stats.profile_problem(problem, functional, [problem['truth'][0]],
                                            fit_result=nominal)
            point = profile['points'][0]
            self.assertEqual(point['status'], 'converged', point)
            covered += point['delta_deviance'] <= 3.841459
        self.assertGreaterEqual(covered, 19,
            'conditional fixed-response Poisson 95% profile threshold fails coarse toy guard')
        self.assertFalse(nominal['publication_ready'])


class ForwardRunnerIntegrationTests(unittest.TestCase):
    @staticmethod
    def write_cache(directory, *, rank_deficient=False):
        cache = directory / 'cache'
        cache.mkdir()
        nphi, epsilon, exposure = 12, 0.72, 4.0
        phi = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
        harmonic = np.column_stack((np.ones(nphi),
            np.sqrt(2 * epsilon * (1 + epsilon)) * np.cos(phi), epsilon * np.cos(2 * phi)))
        migration = np.array([[0.8, 0.3], [0.2, 0.7]])
        truth = np.array([[28.0, -4.0, -5.0], [19.0, 3.0, 2.0]])
        tprime = [-0.6, -0.2]
        expected = exposure * migration @ (truth @ harmonic.T)
        counts = np.random.default_rng(915).poisson(expected)
        with (cache / 'data_events.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['event_id', 'run_number', 'q2', 'xb', 'tprime', 'phi', 'weight'])
            event = 0
            for row in range(2):
                for angular in range(nphi):
                    for _ in range(counts[row, angular]):
                        writer.writerow([event, 100, 4, 0.35, tprime[row], phi[angular], 1])
                        event += 1
        with (cache / 'mc_events.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['event_id', 'reco_q2', 'reco_xb', 'reco_tprime', 'reco_phi',
                             'truth_q2', 'truth_xb', 'truth_tprime', 'truth_phi', 'epsilon', 'base_weight'])
            for block in range(2):
                for angular in range(nphi):
                    for event in range(40):
                        for row in range(2):
                            base = 1e9 * 2 * np.pi * exposure * migration[row, block] / 40
                            writer.writerow([(block * nphi + angular) * 40 + event,
                                4, 0.35, tprime[row], phi[angular], 4, 0.35, tprime[block],
                                0 if rank_deficient else phi[angular], epsilon, base])
        (cache / 'forward_cache_manifest.json').write_text(json.dumps({
            'complete': True, 'schema_version': 1, 'target_divisor': 1.0,
            'target_divisor_error': 0.02}))
        domain = {'q2_edges': [3, 5], 'xb_edges_by_q2': [[0.2, 0.5]],
                  'tprime_edges': [-0.8, -0.4, 0]}
        config = {'truth': domain, 'reco': dict(domain, phi_edges=np.linspace(0, 2*np.pi, 13).tolist()),
                  'publication_truth_blocks': [1]}
        config_path = directory / 'config.json'
        config_path.write_text(json.dumps(config))
        return cache, config_path

    @staticmethod
    def command(cache, config, output, *extra):
        runner = MODULE.parent / 'run_forward_xsec.py'
        return subprocess.run([sys.executable, str(runner), '--cache', str(cache),
            '--config', str(config), '--output', str(output), *extra],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)

    def test_cli_bootstrap_profiles_and_existing_output_preservation(self):
        with tempfile.TemporaryDirectory(prefix='forward-cli-test-') as name:
            directory = Path(name)
            cache, config = self.write_cache(directory)
            profiles = directory / 'profiles.json'
            profiles.write_text(json.dumps([{'name': 'U0', 'functional': [1, 0, 0, 0, 0, 0],
                                             'values': [28]}]))
            output = directory / 'result'
            result = self.command(cache, config, output, '--bootstrap', '4', '--profiles', str(profiles))
            self.assertEqual(result.returncode, 0, result.stderr)
            status = json.loads((output / 'status.json').read_text())
            self.assertTrue(status['complete'])
            self.assertTrue(status['fit_succeeded'])
            self.assertFalse(status['publication_ready'])
            bootstrap = json.loads((output / 'conditional_bootstrap.json').read_text())
            self.assertEqual(bootstrap['successful'], 4)
            self.assertEqual(bootstrap['successful_replicate_indices'], [0, 1, 2, 3])
            self.assertEqual(bootstrap['response_changed_replicates'], 4)
            self.assertEqual(len(json.loads((output / 'fit.json').read_text())['parameters']), 6)
            profile = json.loads((output / 'fixed_response_profiles.json').read_text())
            self.assertEqual(profile[0]['result']['points'][0]['status'], 'converged')
            before = {str(path.relative_to(output)): path.read_bytes()
                      for path in output.rglob('*') if path.is_file()}
            repeat = self.command(cache, config, output, '--bootstrap', '0')
            self.assertNotEqual(repeat.returncode, 0)
            self.assertIn('Output already exists', repeat.stderr)
            self.assertEqual(before, {str(path.relative_to(output)): path.read_bytes()
                                      for path in output.rglob('*') if path.is_file()})

    def test_cli_rank_failure_writes_explicit_incomplete_status(self):
        with tempfile.TemporaryDirectory(prefix='forward-cli-failure-') as name:
            directory = Path(name)
            cache, config = self.write_cache(directory, rank_deficient=True)
            output = directory / 'result'
            result = self.command(cache, config, output, '--bootstrap', '0')
            self.assertNotEqual(result.returncode, 0)
            status = json.loads((output / 'status.json').read_text())
            self.assertFalse(status['complete'])
            self.assertFalse(status['fit_succeeded'])
            self.assertIn('rank-deficient', status['error'])
            self.assertFalse((output / 'coefficients.csv').exists())


if __name__ == '__main__':
    unittest.main()
