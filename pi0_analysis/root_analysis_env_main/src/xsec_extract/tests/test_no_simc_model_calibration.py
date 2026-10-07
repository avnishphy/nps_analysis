import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location(
    "no_simc_model_calibration", ROOT / "src/xsec_extract/no_simc_model_calibration.py")
CALIBRATION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CALIBRATION)


class NoSimcCalibrationTest(unittest.TestCase):
    def test_five_fold_calibration_is_deterministic(self):
        rng = np.random.default_rng(90210)
        residuals = rng.normal(size=(500, 6))
        ids = np.arange(500)
        first = CALIBRATION.calibrate(residuals, ids, [f"q{i}" for i in range(6)])
        second = CALIBRATION.calibrate(residuals, ids, [f"q{i}" for i in range(6)])
        np.testing.assert_array_equal(first["delta"], second["delta"])
        np.testing.assert_array_equal(first["coverage"], second["coverage"])
        self.assertEqual(len(first["coverage_rows"]), 36)
        self.assertTrue(np.all(first["delta"] > 0))

    def test_staging_refuses_overwrite_and_uses_sibling(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "release"
            destination.mkdir()
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                CALIBRATION.staged_directory(destination)
            destination.rmdir()
            resolved, stage = CALIBRATION.staged_directory(destination)
            self.assertEqual(resolved, destination.resolve())
            self.assertEqual(stage.parent, destination.resolve().parent)
            stage.rmdir()

    def test_covariance_representation_keeps_calibrated_diagonal(self):
        covariance = np.array([[4., -1.], [-1., 9.]])
        radii = np.array([3., 5.])
        representation = np.outer(radii, radii) * CALIBRATION.correlation(covariance)
        np.testing.assert_allclose(np.diag(representation), radii ** 2)
        self.assertFalse(np.array_equal(representation, covariance))


if __name__ == "__main__":
    unittest.main()
