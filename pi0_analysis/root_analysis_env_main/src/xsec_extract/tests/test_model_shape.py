import importlib.util
from pathlib import Path
import unittest
import math

spec=importlib.util.spec_from_file_location('model_plots',Path(__file__).resolve().parents[1]/'plot_model_diagnostics.py')
plots=importlib.util.module_from_spec(spec);spec.loader.exec_module(plots)

class RatioGuardTest(unittest.TestCase):
    def test_display_units(self):
        # 0.2 microbarn/GeV2 == 200 nb/GeV2 == 2e-7 microbarn/MeV2.
        self.assertAlmostEqual(float(plots.display(2e-7)),200.)
        self.assertAlmostEqual(float(plots.display(2e-8)),20.)
    def test_signed_ratio(self):
        value,error,status=plots.diagnostic_ratio(-2.,.4,-1.,1.)
        self.assertEqual((value,error,status),(2.,.4,'ok'))
    def test_zero_crossing(self):
        value,error,status=plots.diagnostic_ratio(1e-8,1e-9,1e-13,1e-8)
        self.assertTrue(math.isnan(value));self.assertTrue(math.isnan(error))
        self.assertEqual(status,'near_zero_denominator')
    def test_missing_error_is_not_zero(self):
        value,error,status=plots.diagnostic_ratio(2.,math.nan,1.,1.)
        self.assertEqual(value,2.);self.assertTrue(math.isnan(error));self.assertEqual(status,'ok')
