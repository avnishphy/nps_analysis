"""Display/selection boundary regressions; no changes to production events."""
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from PIL import Image

REPO = Path(__file__).resolve().parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, REPO / path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


collector = module("collector", "scripts/collect_kinematic_plots.py")
combine = module("combine_diagnostics", "src/analysis/combine_analysis_branches.py")


class DiagnosticTests(unittest.TestCase):
    def test_collector_needs_only_existing_images_and_leaves_inputs_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            kin = Path(tmp)/"input"/"KinC_test"
            plots = kin/"plots"
            plots.mkdir(parents=True)
            source = plots/"Q2_overlay_run123.png"
            Image.new("RGB", (200, 150), "white").save(source)
            original = source.read_bytes()
            output, runs, pages = collector.make_kinematic_pdf(kin, Path(tmp)/"pdfs")
            self.assertTrue(output.read_bytes().startswith(b"%PDF"))
            self.assertEqual((runs, pages), (1, 1))
            self.assertEqual(source.read_bytes(), original)
            self.assertEqual(list(kin.rglob("*")), [plots, source])

    def test_layout_keeps_every_representation_without_overlap(self):
        plots = [(Path(f"plot{i}.png"), f"Q2_run{i}.png") for i in range(20)]
        plots += [(Path("mass.png"), "mass_cut_run1.png"),
                  (Path("mass.pdf"), "mass_cut_run1.pdf [1/1]")]
        pages = collector.paginate(plots)
        seen = []
        for page in pages:
            occupied = set()
            for path, label, row, col, span in page:
                seen.append((path, label))
                cells = {(row+r, col+c) for r in range(span) for c in range(span)}
                self.assertFalse(cells & occupied)
                self.assertTrue(all(r < collector.ROWS and c < collector.COLS for r,c in cells))
                occupied |= cells
        self.assertCountEqual(seen, plots)

    def test_exterior_trim_preserves_faint_nonwhite_annotation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"source.png"
            image = Image.new("RGB", (200, 200), "white")
            image.putpixel((100,100), (0,0,0))
            image.putpixel((20,20), (254,254,254))
            image.save(path)
            result = collector.plot_image(path)
            pixels = np.asarray(result)
            self.assertEqual(np.count_nonzero(np.any(pixels != 255, axis=2)), 2)
            with Image.open(path) as original:
                self.assertEqual(original.size, (200,200))

    def test_failed_fit_keeps_production_flags_and_fills_fallback_from_events(self):
        frame = pd.DataFrame({"mpi0_all": [.135,.135,.135], "mmiss_all": [.93,.93,.93],
                              "pi0_weight": [1.,2.,3.], "scale": [2.,2.,2.],
                              "is_exclusive": [1,0,1]})
        original = frame.copy(deep=True)
        with patch.object(combine, "_fit_combined_2d_mass_cut", return_value=None):
            debug = combine.add_combined_2d_mass_cut(frame)
        pd.testing.assert_frame_equal(original, frame)
        self.assertEqual(debug["params"]["diagnostic_fallback_decorrelation"], 1)
        self.assertAlmostEqual(debug["params"]["legacy_is_exclusive_total_fraction"], 8/12)
        tag = combine.COMBINED_MASS_CUT_TAG
        self.assertEqual(debug["histograms"][f"{tag}_h_legacy_is_exclusive_selected"][0].sum(), 8)
        self.assertEqual(debug["histograms"][f"{tag}_h_mmiss_vs_mpi0_weighted"][0].sum(), 12)
        with tempfile.TemporaryDirectory() as tmp:
            combine.write_combined_mass_cut_canvas(debug, Path(tmp)/"combined.root")
            self.assertTrue(list(Path(tmp).glob("*.pdf")))

    def test_four_bin_fit_is_not_qualified_even_if_solver_reports_valid(self):
        params = {"valid": 1, "ellipse_valid": 1, "fit_subset_bins": 4,
                  "fit_subset_total_fraction": .000303725, "cov_det": 1e-12}
        self.assertIn("too few", combine.ellipse_diagnostic_failure(params))

    def test_combined_physics_overlays_keep_every_existing_selector(self):
        frame = pd.DataFrame({"Q2": [4.,5.,6.], "scale": [1.,1.,1.], "pi0_weight": [1.,2.,3.]})
        for branch in ["is_exclusive", "is_exclusive_ellipse", "is_exclusive_mcd",
                       "is_exclusive_ellipse_combined", "is_exclusive_mcd_combined"]:
            frame[branch] = [1,0,1]
        original = frame.copy(deep=True)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"physics.pdf"
            combine.create_analysis_plots(frame, path)
            self.assertGreater(path.stat().st_size, 0)
        pd.testing.assert_frame_equal(frame, original)


if __name__ == "__main__":
    unittest.main()
