"""Diagnostic regression check of recovery CSV and ROOT contracts."""
import csv
import math
import sys
from pathlib import Path

import ROOT

base = Path(sys.argv[1])
for case, expected in {"healthy": 4, "partial": 3, "nuisance": 3, "rank": 2, "all_failed": 0}.items():
    directory = base / case

    def rows(name):
        with (directory / name).open() as stream:
            return list(csv.DictReader(stream))

    status = rows("fit_status.csv")
    assert len(status) == 4
    assert sum(int(row["fit_ok"]) for row in status) == expected
    assert all(row["fit_scope"] == ("global_migration" if case == "healthy" else "global_migration_subset") for row in status)
    slices = rows("slices.csv")
    assert sum(int(row["fit_xsec_ok"]) for row in slices) == expected * 2
    for row in slices:
        if row["fit_xsec_ok"] == "0":
            assert math.isnan(float(row["fit_xsec_sigmaU"]))
            assert row["fit_failure_reason"]
    for row in rows("summary.csv"):
        if row["fit_xsec_ok"] == "0":
            assert all(math.isnan(float(row[key])) for key in ("sim", "sim_err", "ratio", "xsec", "xsec_err"))
    reco = rows("migration_reco_rows.csv")
    assert sum(row["exclusion"] == "excluded_q2_xb" for row in reco) == (4 - expected) * 24
    for row in reco:
        if row["exclusion"] == "excluded_q2_xb":
            assert row["fit_index"] == "-1" and math.isnan(float(row["prediction"]))
    with ROOT.TFile.Open(str(directory / "results.root")) as result:
        assert result.Get("fit_status").GetEntries() == 4
        assert result.Get("fit_attempts").GetEntries() == len(rows("fit_attempts.csv"))
        assert result.Get("migration_response_cells").GetEntries() == 96 * 14
        assert result.Get("fit_retained_q2_xb_groups").GetVal() == expected
    if case == "nuisance":
        nuisance = [row for row in rows("migration_parameters.csv") if row["truth_block"] == "0"]
        assert len(nuisance) == 3 and all(row["is_nuisance"] == "1" for row in nuisance)
print("PASS recovery CSV/ROOT contracts")
