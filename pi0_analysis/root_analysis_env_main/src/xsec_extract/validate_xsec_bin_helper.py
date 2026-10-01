"""Execute real ROOT inputs in the required kernel and check the bin contract.

Run from the repository: /group/nps/singhav/software/python/bin/python
src/xsec_extract/validate_xsec_bin_helper.py --case both (or single).
Artifacts go to output/lt_bin_validation; source notebook is never executed in place.
"""

from pathlib import Path
import argparse
import json
import os
import sys

ROOT = Path(__file__).resolve().parents[2]
PYTHON = Path("/group/nps/singhav/software/python/bin/python")


def kernel_checks(ns, output):
    import numpy as np
    import pandas as pd
    from src.xsec_extract.xsec_bin_dashboard import (
        BinDashboard,
        in_bin,
        _quantile_edges,
    )

    d = ns["dashboard"]
    output = Path(output)
    notebook = json.loads((ROOT / "src/xsec_extract/xsec_bin_helper.ipynb").read_text())
    checks = []

    def audit(label):
        d.require_ready()
        polygon = d._polygon()
        counts = []
        for name, selected in zip(d.settings, d.selected_frames):
            raw = ns["all_data"][name]
            finite = np.logical_and.reduce(
                [
                    np.isfinite(raw[k])
                    for k in [
                        "Q2",
                        "xB",
                        "W",
                        "tprime",
                        "phi_wrapped",
                        "mmiss_all",
                        "data_weight",
                    ]
                ]
            )
            mask = finite & raw.mmiss_all.between(*ns["MMISS_LIMITS"])
            mask &= (
                raw.tprime.between(*d.t_limits.value)
                & raw.Q2.between(*d.q_limits.value)
                & raw.xB.between(*d.x_limits.value)
            )
            if polygon is not None:
                # Independent convex half-plane test, not the dashboard Path mask.
                for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
                    cross = (b[0] - a[0]) * (raw.Q2 - a[1]) - (b[1] - a[1]) * (
                        raw.xB - a[0]
                    )
                    mask &= cross >= -1e-12
            assert set(raw.index[mask]) == set(selected.index), (
                label,
                name,
                "mask mismatch",
            )
            assert f"{len(selected):,}" in d.status.value
            occ = d.occupancy_frame.loc[d.occupancy_frame.setting == name]
            assert occ.events.sum() == len(selected), (
                label,
                name,
                "lost/doubled events",
            )
            # Verify every proposed Q2/xB/tprime/phi bin, including empty bins.
            for row in occ.to_dict("records"):
                iq, ix, it, ip = [
                    row[k] - 1 for k in ["Q2 bin", "xB bin", "t' bin", "phi bin"]
                ]
                expected = (
                    in_bin(selected.Q2, d.edges["Q2"], iq)
                    & in_bin(selected.xB, d.edges["xB_by_Q2"][iq], ix)
                    & in_bin(selected.tprime, d.edges["tprime"], it)
                    & in_bin(selected.phi_wrapped, d.edges["phi"], ip)
                ).sum()
                assert row["events"] == expected
                assert row["sparse"] == (expected < d.sparse_min.value)
            counts.append(len(selected))
        for stats, x, y in zip(d.map_stats, ["Q2", "xB"], ["W", "Q2"]):
            for f, h in zip(d.current_coverage, stats["counts"]):
                np.testing.assert_array_equal(
                    h,
                    np.histogram2d(
                        f[x], f[y], bins=(stats["x_edges"], stats["y_edges"])
                    )[0],
                )
            occupied = stats["counts"] >= d.cell_min.value
            if len(d.settings) > 1:
                np.testing.assert_array_equal(
                    stats["category"] == len(d.settings) + 2, occupied.all(axis=0)
                )
        assert ns["data_selected"].shape[0] == sum(counts)
        np.testing.assert_array_equal(ns["TPRIME_EDGES"], d.edges["tprime"])
        checks.append({"check": label, "counts": counts, "revision": d.revision})

    audit("default")
    baseline = [len(f) for f in d.selected_frames]
    vertices = [(0.285, 3.25), (0.34, 3.3), (0.435, 4.7), (0.375, 4.6)]
    for v in vertices:
        d.place_corner(*v)
    audit("four-corner diamond")
    for imap in range(2):
        for i in range(1, 5):
            assert f'id="corner-{imap}-{i}"' in d.canvas.svg
    d.export_figure(output / "selected_diamond")
    d.figure.savefig(output / "selected_diamond.png", dpi=150)
    original = d.vertices.copy()
    d.corner.value = 3
    d.place_corner(0.37, 4.5)
    audit("move corner 4")
    d._undo()
    assert d.vertices == original
    audit("undo move")
    d._reset()
    audit("clear")
    assert [len(f) for f in d.selected_frames] == baseline
    d._undo()
    assert d.vertices == original

    # Downstream cells must use the same cut and exact edge arrays.
    for i in (13, 15, 17):
        exec("".join(notebook["cells"][i]["source"]), ns)
    np.testing.assert_array_equal(ns["q2_edges_final"], d.edges["Q2"])
    assert ns["config_edges"]["diamond_xb_q2_vertices"] == ns["DIAMOND_XB_Q2_VERTICES"]
    assert (
        ns["bin_contract"]["selection"]["diamond_xb_q2_vertices"]
        == ns["DIAMOND_XB_Q2_VERTICES"]
    )
    expected_reference = len(
        d.selected_frames[d.settings.index(ns["REFERENCE_KINEMATIC"])]
    )
    assert (
        ns["acceptance_weighted_table"].data_events_phi_integrated.sum()
        == expected_reference
    )
    ns["WRITE_EXPORT"] = True
    try:
        exec("".join(notebook["cells"][17]["source"]), ns)
    except RuntimeError as exc:
        assert "Export requires one setting" in str(exc)
    else:
        raise AssertionError("Multi-setting/missing-SIMC export guard bypassed")
    finally:
        ns["WRITE_EXPORT"] = False
    if len(d.settings) > 1:
        saved_settings = ns["KINEMATICS"]
        ns["KINEMATICS"] = [ns["REFERENCE_KINEMATIC"]]
        ns["WRITE_EXPORT"] = True
        try:
            exec("".join(notebook["cells"][17]["source"]), ns)
        except RuntimeError as exc:
            assert "Refusing export" in str(exc)
        else:
            raise AssertionError("Automatic Q2/xB edges were exportable")
        finally:
            ns["KINEMATICS"] = saved_settings
            ns["WRITE_EXPORT"] = False
    rev = d.revision
    d.n_phi.value = 12
    try:
        d.require_ready(rev)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Stale representative coordinates accepted")
    audit("phi bins changed")
    d.slice.value = 0
    for f, h in zip(d.selected_frames, d.occupancy_matrices):
        mask = in_bin(f.Q2, d.edges["Q2"], 0) & in_bin(f.xB, d.edges["xB_by_Q2"][0], 0)
        assert h.sum() == mask.sum()
    d.slice.value = -1

    for mode in ["uniform", "yield_quantile", "count_quantile"]:
        d._busy = True
        d.t_mode.value = d.q_mode.value = d.x_mode.value = mode
        d._busy = False
        d.refresh()
        audit(mode)
        for key in ("tprime", "Q2"):
            e = d.edges[key]
            assert e[0] == d.limits[key][0] and e[-1] == d.limits[key][1]
        if mode == "count_quantile":
            pooled = pd.concat(d.selected_frames)
            np.testing.assert_allclose(
                d.edges["tprime"][1:-1],
                np.quantile(pooled.tprime, np.arange(1, d.n_t.value) / d.n_t.value),
            )
    # Manual arrays control bin cardinality, including conditional uneven xB bins.
    d.manual_fields["tprime"].value = "[-0.7,-0.25,0]"
    d.manual_fields["Q2"].value = "[3,4,5]"
    d.manual_fields["xB"].value = "[[0.25,0.34,0.455],[0.25,0.455]]"
    d._apply_manual()
    audit("manual conditional bins")
    assert d.n_t.value == 2 and [len(e) for e in d.edges["xB_by_Q2"]] == [3, 2]
    d.manual_fields["Q2"].value = "[3,3,5]"
    d._apply_manual()
    assert not d.valid and not d.edges and ns["TPRIME_EDGES"].size == 0
    try:
        d.export_figure(output / "must_not_exist")
    except RuntimeError:
        pass
    else:
        raise AssertionError("Invalid selection exported")
    d.manual_fields["Q2"].value = "[3,4,5]"
    d._apply_manual()
    d._busy = True
    d.t_mode.value = d.q_mode.value = d.x_mode.value = "uniform"
    d._busy = False
    d._reset()
    # Verify expansion uses retained input rather than an irreversibly pre-cut frame.
    applied = d._apply_limits(limits={
        "tprime": ["-1.0", "0.0"],
        "Q2": list(d.q_limits.value),
        "xB": list(d.x_limits.value),
    })
    assert applied and d.limit_fields["tprime"][0].value == -1.0
    audit("expanded tprime")
    d.limit_fields["xB"][0].value = 0.8
    d.limit_fields["xB"][1].value = 0.9
    d._apply_limits()
    assert not d.valid and not d.edges and all(f.empty for f in d.selected_frames)
    d.limit_fields["xB"][0].value = 0.25
    d.limit_fields["xB"][1].value = 0.455
    d._apply_limits()
    d.place_corner(0.3, 3.4)
    try:
        d.require_ready()
    except RuntimeError:
        pass
    else:
        raise AssertionError("Incomplete diamond accepted downstream")
    d._reset()
    # Collapse and weighted-quantile conventions on a small analytic sample.
    d.manual_fields["tprime"].value = json.dumps(np.linspace(-1, 0, 26).tolist())
    d.manual_fields["Q2"].value = "[3,5]"
    d.manual_fields["xB"].value = "[[0.25,0.455]]"
    d._apply_manual()
    audit("25 manual tprime bins")
    assert d.n_t.value == ns["N_TPRIME"] == 25 and len(d.edges["tprime"]) == 26
    np.testing.assert_allclose(
        _quantile_edges([0, 1, 2], [1, 2, 1], (0, 2), 2, "yield_quantile"), [0, 1, 2]
    )
    try:
        _quantile_edges([1, 1, 1], [1, 1, 1], (0, 2), 4, "count_quantile")
    except ValueError:
        pass
    else:
        raise AssertionError("Collapsed quantiles accepted")
    # Multi-setting categorical overlap and W-branch guard.
    f = ns["coverage_frames"][0]
    multi = BinDashboard(
        ["A", "B", "C"],
        [f, f.copy(), f.copy()],
        ns["PROTON_MASS_GEV"],
        {"tprime": (-0.7, 0), "Q2": (3, 5), "xB": (0.25, 0.455)},
        dict.fromkeys(["tprime", "Q2", "xB"], "uniform"),
        {"tprime": 2, "Q2": 1, "xB": 1, "phi": 4},
    )
    assert multi.valid and (multi.map_stats[0]["category"] == 5).any()
    bad = f.copy()
    bad["W"] += 0.01
    try:
        BinDashboard(
            ["bad"],
            [bad],
            ns["PROTON_MASS_GEV"],
            multi.limits,
            dict.fromkeys(["tprime", "Q2", "xB"], "uniform"),
            {"tprime": 2, "Q2": 1, "xB": 1, "phi": 4},
        )
    except ValueError as exc:
        assert "W branch" in str(exc)
    else:
        raise AssertionError("Inconsistent W branch accepted")
    result = {
        "python": sys.executable,
        "settings": d.settings,
        "projection": d.projection_checks,
        "checks": checks,
        "stale_export_blocked": True,
        "invalid_export_blocked": True,
        "incomplete_cut_blocked": True,
        "all_full_bin_counts_exact": True,
        "multi_setting_checked": True,
    }
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print("VALIDATION PASSED", json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=["both", "single"], default="both")
    parser.add_argument(
        "--output", type=Path, default=ROOT / "output/lt_bin_validation"
    )
    args = parser.parse_args()
    if Path(sys.executable).resolve() != PYTHON.resolve():
        raise SystemExit(f"Use {PYTHON}")
    import nbformat
    from nbclient import NotebookClient

    out = args.output.resolve() / args.case
    out.mkdir(parents=True, exist_ok=True)
    n = nbformat.read(ROOT / "src/xsec_extract/xsec_bin_helper.ipynb", as_version=4)
    n.cells.append(
        nbformat.v4.new_code_cell(
            "from src.xsec_extract.validate_xsec_bin_helper import kernel_checks\n"
            f"kernel_checks(globals(), {str(out)!r})"
        )
    )
    os.environ["PATH"] = str(PYTHON.parent) + os.pathsep + os.environ["PATH"]
    os.environ["NPS_BIN_SETTINGS"] = (
        "KinC_x36_4" if args.case == "single" else "KinC_x36_4,KinC_x36_5_407"
    )
    client = NotebookClient(
        n,
        timeout=360,
        kernel_name="python3",
        resources={"metadata": {"path": str(ROOT)}},
    )
    try:
        client.execute()
    finally:
        nbformat.write(n, out / "executed.ipynb")
    print((out / "results.json").read_text())


if __name__ == "__main__":
    main()
