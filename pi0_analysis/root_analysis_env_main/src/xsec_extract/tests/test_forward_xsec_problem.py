"""Deterministic forward-response conservation, support and information tests."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forward_xsec_problem import build_problem, diagnostics, mc_prediction_variance


DATA_FIELDS = ["event_id", "run_number", "q2", "xb", "tprime", "phi", "weight"]
MC_FIELDS = ["event_id", "reco_q2", "reco_xb", "reco_tprime", "reco_phi",
             "truth_q2", "truth_xb", "truth_tprime", "truth_phi", "epsilon", "base_weight"]


def config():
    return {"reco": {"q2_edges": [1, 3], "xb_edges_by_q2": [[.1, .9]],
                      "tprime_edges": [-1, -.5, 0], "phi_edges": np.linspace(0, 2*np.pi, 5).tolist()},
            "truth": {"q2_edges": [1, 3], "xb_edges_by_q2": [[.1, .9]],
                       "tprime_edges": [-1, -.5, 0]},
            "publication_truth_blocks": [0, 1]}


def data(event, t=-.75, phi=45, weight=1, run=1):
    return [event, run, 2, .5, t, np.deg2rad(phi), weight]


def mc(event, reco_t=-.75, truth_t=-.75, phi=45, truth_q=2, truth_x=.5, weight=1e9,
       reco_phi=None, epsilon=.6):
    return [event, 2, .5, reco_t, np.deg2rad(phi if reco_phi is None else reco_phi),
            truth_q, truth_x, truth_t, np.deg2rad(phi), epsilon, weight]


class ForwardProblemTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.cache = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write(self, data_rows, mc_rows):
        for name, fields, rows in (("data_events.csv", DATA_FIELDS, data_rows),
                                   ("mc_events.csv", MC_FIELDS, mc_rows)):
            with (self.cache / name).open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(fields)
                writer.writerows(rows)

    def test_basis_uses_vertex_phi_and_absolute_normalization(self):
        self.write([data(1, phi=90, weight=2), data(2, t=-.2, phi=0, weight=3)],
                   [mc(1, phi=0, reco_phi=90), mc(2, reco_t=-.2, truth_t=-.2, phi=0)])
        problem = build_problem(self.cache, config())
        expected = np.asarray([1, np.sqrt(2*.6*1.6), .6]) / (2*np.pi)
        np.testing.assert_allclose(problem["design"][1, :3], expected)
        np.testing.assert_allclose(problem["mc_basis"][0], expected)
        self.assertEqual(problem["y"].sum(), 5)
        self.assertEqual(problem["sumw2"].sum(), 13)
        self.assertEqual(problem["design"].shape, (8, 6))
        self.assertEqual(np.count_nonzero(problem["y"]), 2)

    def test_reco_refinement_preserves_response_and_truth_coefficients(self):
        self.write([data(1, phi=30), data(2, t=-.2, phi=120)],
                   [mc(1, phi=30), mc(2, reco_t=-.2, truth_t=-.2, phi=120),
                    mc(3, truth_t=-1.2, phi=200)])
        coarse = config()
        coarse["reco"]["phi_edges"] = [0, np.pi, 2*np.pi]
        fine_problem = build_problem(self.cache, config())
        coarse_problem = build_problem(self.cache, coarse)
        np.testing.assert_allclose(fine_problem["design"].reshape(4, 2, -1).sum(axis=1),
                                   coarse_problem["design"])
        np.testing.assert_array_equal(fine_problem["active_global_blocks"], coarse_problem["active_global_blocks"])
        np.testing.assert_allclose(fine_problem["y"].reshape(4, 2).sum(axis=1), coarse_problem["y"])

    def test_unpublished_buffer_remains_independent_nuisance(self):
        self.write([data(1)], [mc(1), mc(2, truth_t=-.2), mc(3, truth_t=-1.2)])
        cfg = config()
        cfg["publication_truth_blocks"] = [1]
        problem = build_problem(self.cache, cfg)
        self.assertEqual(problem["design"].shape[1], 9)
        self.assertEqual(problem["published_blocks"].tolist(), [1])
        self.assertEqual([block["published"] for block in problem["truth_blocks"]], [False, True, False])
        self.assertEqual(problem["truth_blocks"][2]["name"], "tprime_below")

    def test_six_exterior_faces_keep_legacy_corner_priority(self):
        extra = [mc(10, truth_t=-1.1, truth_q=.5, truth_x=.05),
                 mc(11, truth_t=.1, truth_q=4, truth_x=.95),
                 mc(12, truth_q=.5, truth_x=.05), mc(13, truth_q=4, truth_x=.95),
                 mc(14, truth_x=.05), mc(15, truth_x=.95)]
        self.write([data(1)], [mc(1), mc(2, truth_t=-.2)] + extra)
        problem = build_problem(self.cache, config())
        self.assertEqual(problem["mc_blocks"].tolist(), list(range(8)))
        self.assertEqual(problem["design"].shape, (8, 24))

    def test_edges_final_upper_and_phi_wrap(self):
        self.write([data(1, t=0, phi=360)],
                   [mc(1, truth_t=-1, phi=0), mc(2, reco_t=0, truth_t=0, phi=360)])
        problem = build_problem(self.cache, config())
        self.assertEqual(problem["data_rows"].tolist(), [4])
        self.assertEqual(problem["mc_blocks"].tolist(), [0, 1])
        self.assertEqual(problem["mc_rows"].tolist(), [0, 4])

    def test_shifted_phi_origin_wraps_both_endpoints(self):
        self.write([data(1, phi=-45), data(2, phi=315), data(3, phi=0), data(4, phi=45)],
                   [mc(1, phi=-45), mc(2, truth_t=-.2, phi=315),
                    mc(3, phi=0), mc(4, phi=45)])
        cfg = config()
        cfg["reco"]["phi_edges"] = np.linspace(-np.pi/4, 7*np.pi/4, 5).tolist()
        problem = build_problem(self.cache, cfg)
        self.assertEqual(problem["data_rows"].tolist(), [0, 0, 0, 1])
        self.assertEqual(problem["mc_rows"].tolist(), [0, 0, 0, 1])
        cfg["reco"]["phi_edges"] = [-np.pi/4, 0, np.pi]
        with self.assertRaisesRegex(ValueError, "full 2\\*pi period"):
            build_problem(self.cache, cfg)

    def test_nonuniform_xb_bin_counts_have_stable_global_ids(self):
        self.write([], [mc(1, truth_q=1.5, truth_x=.2), mc(2, truth_q=2.5, truth_x=.2),
                        mc(3, truth_q=2.5, truth_x=.7), mc(4, truth_q=2.5, truth_x=.7, truth_t=-.2)])
        cfg = config()
        cfg["truth"]["q2_edges"] = [1, 2, 3]
        cfg["truth"]["xb_edges_by_q2"] = [[.1, .9], [.1, .5, .9]]
        cfg["publication_truth_blocks"] = [0, 1, 2, 5]
        problem = build_problem(self.cache, cfg)
        self.assertEqual(problem["active_global_blocks"].tolist(), [0, 1, 2, 5])

    def test_resampling_universe_is_independent_of_reco_selection(self):
        self.write([data(1), data(2, t=-.2)],
                   [mc(1), mc(2, reco_t=-.2, truth_t=-.2), mc(3, truth_t=-.2)])
        full = build_problem(self.cache, config())
        cfg = config()
        cfg["reco"]["tprime_edges"] = [-1, -.5]
        narrow = build_problem(self.cache, cfg)
        self.assertLess(len(narrow["data_ids"]), len(full["data_ids"]))
        self.assertLess(len(narrow["mc_ids"]), len(full["mc_ids"]))
        np.testing.assert_array_equal(narrow["bootstrap_data_ids"], full["bootstrap_data_ids"])
        np.testing.assert_array_equal(narrow["bootstrap_mc_ids"], full["bootstrap_mc_ids"])

    def test_support_and_cache_completion_are_hard_errors(self):
        self.write([data(1, phi=150)], [mc(1), mc(2, truth_t=-.2)])
        with self.assertRaisesRegex(ValueError, "outside MC response support"):
            build_problem(self.cache, config())
        self.write([], [mc(1)])
        with self.assertRaisesRegex(ValueError, "Published truth blocks"):
            build_problem(self.cache, config())
        (self.cache / "forward_cache_manifest.json").write_text('{"complete":false}')
        with self.assertRaisesRegex(ValueError, "incomplete"):
            build_problem(self.cache, config())

    def test_nonfinite_negative_and_bad_edges_are_rejected(self):
        for bad_weight in (-1, float("nan")):
            self.write([data(1, weight=bad_weight)], [mc(1), mc(2, truth_t=-.2)])
            with self.assertRaises(ValueError):
                build_problem(self.cache, config())
        self.write([data(1)], [mc(1, epsilon=1.1), mc(2, truth_t=-.2)])
        with self.assertRaisesRegex(ValueError, "epsilon"):
            build_problem(self.cache, config())
        cfg = config()
        cfg["truth"]["tprime_edges"] = [-1, -.5, -.5, 0]
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            build_problem(self.cache, cfg)

    def test_independent_event_ids_and_effective_counts(self):
        self.write([data(7, weight=2), data(7, weight=3), data(7, weight=4, run=2)],
                   [mc(8), mc(8), mc(9, truth_t=-.2)])
        problem = build_problem(self.cache, config())
        self.assertEqual(problem["y"][0], 9)
        self.assertEqual(problem["sumw2"][0], 25 + 16)
        self.assertEqual(problem["data_ids"][0], problem["data_ids"][1])
        self.assertNotEqual(problem["data_ids"][0], problem["data_ids"][2])
        report = diagnostics(problem)
        self.assertEqual(report["blocks"][0]["u_response_effective_events"], 1)
        self.assertEqual(report["blocks"][0]["mc_independent_events"], 1)

    def test_mc_variance_preserves_event_and_harmonic_covariance(self):
        self.write([], [mc(8, phi=0), mc(8, phi=0), mc(8, truth_t=-.2, phi=0),
                        mc(9, truth_t=-.2, phi=0)])
        problem = build_problem(self.cache, config())
        coefficients = np.asarray([4, .5, .2, 2, .1, .05])
        predictions = np.sum(problem["mc_basis"] * coefficients.reshape(-1, 3)[problem["mc_blocks"]], axis=1)
        variance = mc_prediction_variance(problem, coefficients)
        self.assertAlmostEqual(variance[0], predictions[:3].sum()**2 + predictions[3]**2)
        self.assertEqual(np.count_nonzero(variance), 1)

    def test_profile_information_matches_orthogonal_projection(self):
        # Publication column 0 has half its information shared with nuisance;
        # other publication columns are orthogonal. No matrix inverse needed.
        design = np.zeros((6, 6))
        design[:3, :3] = np.eye(3)
        design[0, 3] = 1
        design[3, 3] = 1
        design[4, 4] = 1
        design[5, 5] = 1
        problem = {"design": design, "published_blocks": np.asarray([0]),
                   "truth_blocks": [], "mc_blocks": np.asarray([], dtype=int),
                   "mc_basis": np.empty((0, 3)), "mc_ids": np.asarray([])}
        report = diagnostics(problem, row_variance=np.ones(6))
        np.testing.assert_allclose(sorted(report["retained_publication_information_fractions"]), [.5, 1, 1])
        self.assertEqual(report["publication_rank_profiled"], 3)
        json.dumps(report, allow_nan=False)
        problem["design"] = np.eye(3)
        report = diagnostics(problem)
        np.testing.assert_allclose(report["retained_publication_information_fractions"], [1, 1, 1])
        self.assertEqual(report["nuisance_rank"], 0)
        problem["design"] = np.column_stack((np.eye(3), np.eye(3)))
        report = diagnostics(problem)
        self.assertEqual(report["publication_rank_profiled"], 0)


if __name__ == "__main__":
    unittest.main()
