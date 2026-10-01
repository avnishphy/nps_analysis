"""Deterministic checks of common-integral paired covariance and provenance."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

MODULE = Path(__file__).resolve().parents[1] / 'src/xsec_extract/compare_forward_xsec.py'
SPEC = importlib.util.spec_from_file_location('compare_forward_xsec', MODULE)
comparison = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(comparison)


class ComparisonTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.a, self.b = [Path(self.temp.name) / name for name in ('a', 'b')]
        self.observables = [dict(name='U', component='U', q2=[2, 3], xb=[.2, .4], tprime=[-1, 0])]
        for path, edges, central in ((self.a, [-1, -.5, 0], [10, 10]), (self.b, [-1, 0], [11])):
            path.mkdir()
            blocks = [dict(local_id=i, global_id=i, kind='interior', published=True,
                           q2_bounds=[2, 3], xb_bounds=[.2, .4], tprime_bounds=edges[i:i+2])
                      for i in range(len(central))]
            parameters = np.zeros(3*len(central))
            parameters[::3] = central
            samples = np.tile(parameters, (4, 1))
            samples[:, ::3] += np.array([-2, -1, 1, 2])[:, None]
            records = lambda names: [dict(path='/cache/'+n, sha256='a'*64) for n in names]
            payloads = {
                'fit.json': dict(status='converged', parameters=parameters.tolist()),
                'bin_metadata.json': dict(truth_blocks=blocks),
                'conditional_bootstrap.json': dict(samples=samples.tolist(), successful=4,
                    requested=4, seed=17, successful_replicate_indices=[0, 1, 2, 3],
                    nominal_parameters=parameters.tolist(), failures=[]),
                'provenance.json': dict(seed=17, numpy=np.__version__,
                    inputs=records(['data_events.csv', 'mc_events.csv', 'forward_cache_manifest.json']),
                    sources=records(['forward_xsec_statistics.py', 'forward_xsec_problem.py'])),
            }
            for name, value in payloads.items():
                (path/name).write_text(json.dumps(value))

    def test_shared_noise_cancels_and_failures_align(self):
        result = comparison.compare(self.a, self.b, self.observables)
        np.testing.assert_allclose(result['central_difference'], [1])
        np.testing.assert_allclose(result['paired_difference_covariance'], [[0]])
        np.testing.assert_allclose(result['paired_crosscovariance_A_B'], [[10/3]])
        path = self.b/'conditional_bootstrap.json'
        boot = json.loads(path.read_text())
        boot['successful_replicate_indices'] = [0, 2, 3]
        boot['samples'].pop(1)
        boot['successful'] = 3
        path.write_text(json.dumps(boot))
        result = comparison.compare(self.a, self.b, self.observables)
        self.assertEqual(result['paired_replicate_indices'], [0, 2, 3])
        self.assertEqual(result['replicates']['B']['failed_indices'], [1])
        self.assertEqual(result['replicates']['A']['successful_but_unpaired_indices'], [1])
        np.testing.assert_allclose(result['paired_difference_covariance'], [[0]])

    def test_cache_exposure_source_and_seed_mismatch_rejected(self):
        path = self.b/'provenance.json'
        original = json.loads(path.read_text())
        for section in ('inputs', 'sources'):
            for index in range(len(original[section])):
                changed = json.loads(json.dumps(original))
                changed[section][index]['sha256'] = 'b'*64
                path.write_text(json.dumps(changed))
                with self.assertRaises(ValueError):
                    comparison.compare(self.a, self.b, self.observables)
        original['seed'] = 18
        path.write_text(json.dumps(original))
        with self.assertRaises(ValueError):
            comparison.compare(self.a, self.b, self.observables)

    def test_partial_integral_and_wrong_q2_rejected(self):
        for key, value in [('tprime', [-1, -.5]), ('q2', [2, 2.5])]:
            request = [dict(self.observables[0], **{key: value})]
            with self.assertRaises(ValueError):
                comparison.compare(self.a, self.b, request)

    def test_cli_preserves_existing_output(self):
        root = Path(self.temp.name)
        observations, output = root/'observables.json', root/'result.json'
        observations.write_text(json.dumps(self.observables))
        args = [str(self.a), str(self.b), '--observables', str(observations), '--output', str(output)]
        comparison.main(args)
        before = output.read_bytes()
        with self.assertRaises(SystemExit):
            comparison.main(args)
        self.assertEqual(output.read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
