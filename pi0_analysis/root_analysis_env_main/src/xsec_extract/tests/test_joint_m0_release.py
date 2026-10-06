from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from joint_m0_release import _calibrate, PUBLICATION_TOYS


class JointM0ReleaseTest(unittest.TestCase):
    def test_deterministic_five_fold_calibration(self):
        rng=np.random.default_rng(2718)
        residuals=rng.normal(size=(PUBLICATION_TOYS,6))
        ids=np.arange(PUBLICATION_TOYS)
        result=_calibrate(residuals,ids,[f"q{i}" for i in range(6)])
        self.assertEqual(len(result["coverage_rows"]),36)
        self.assertTrue(np.all(result["delta"]>0))
        self.assertTrue(np.all((result["coverage"]>=0)&(result["coverage"]<=1)))
        aggregate=result["coverage_rows"][:6]
        self.assertTrue(all(row["validation_n"]==PUBLICATION_TOYS for row in aggregate))
        fold_rows=result["coverage_rows"][6:]
        self.assertTrue(all(row["training_n"]==400 and row["validation_n"]==100 for row in fold_rows))


if __name__=="__main__":
    unittest.main()
