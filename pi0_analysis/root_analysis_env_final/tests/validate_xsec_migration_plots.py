#!/usr/bin/env python3
"""Validate migration ROOT products against independent CSV records.

Run in the sourced NPS/ROOT environment. Inputs are read-only. An optional
baseline checks that adding diagnostics preserved the numerical extraction.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import ROOT

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('output', type=Path)
parser.add_argument('--baseline', type=Path)
parser.add_argument('--json', type=Path)
args = parser.parse_args()


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def equal(a, b):
    return (math.isnan(a) and math.isnan(b)) or math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-25)


root = ROOT.TFile.Open(str(args.output / 'excl_xsec_pi0_analysis_no_simc_model_output.root'), 'READ')
assert root and not root.IsZombie(), 'Missing output ROOT file'
hist = {name: root.Get('migration/' + name) for name in (
    'response_U', 'response_LT', 'response_TT', 'selected_counts',
    'p_reco_given_truth_selected', 'p_truth_given_reco_selected',
    'vertex_q2_xb_selected', 'reco_q2_xb_selected')}
assert all(hist.values()), 'Missing migration histogram'
cells = rows(args.output / 'migration_response_cells.csv')
reco = rows(args.output / 'migration_reco_rows.csv')
truth = rows(args.output / 'migration_truth_blocks.csv')
row_to_slice = {}
slice_keys = sorted({(int(r['it']), int(r['iq']), int(r['ix'])) for r in reco})
for row in reco:
    row_to_slice[int(row['reco_row'])] = slice_keys.index((int(row['it']), int(row['iq']), int(row['ix'])))
nb, ns, nr = len(truth), len(slice_keys), len(reco)
counts = [[0 for _ in range(nb)] for _ in range(ns)]
for term in ('U', 'LT', 'TT'):
    assert hist['response_' + term].GetNbinsX() == nb
    assert hist['response_' + term].GetNbinsY() == nr
for cell in cells:
    r, b = int(cell['reco_row']), int(cell['truth_block'])
    counts[row_to_slice[r]][b] += int(cell['events'])
    for term in ('U', 'LT', 'TT'):
        assert equal(hist['response_' + term].GetBinContent(b + 1, r + 1), float(cell['basis_' + term])), (term, r, b)
truth_sums = [sum(row[b] for row in counts) for b in range(nb)]
reco_sums = [sum(row) for row in counts]
for s in range(ns):
    for b in range(nb):
        n = counts[s][b]
        assert hist['selected_counts'].GetBinContent(b + 1, s + 1) == n
        assert equal(hist['p_reco_given_truth_selected'].GetBinContent(b + 1, s + 1), n / truth_sums[b] if truth_sums[b] else 0)
        assert equal(hist['p_truth_given_reco_selected'].GetBinContent(b + 1, s + 1), n / reco_sums[s] if reco_sums[s] else 0)
for b, total in enumerate(truth_sums):
    assert equal(sum(hist['p_reco_given_truth_selected'].GetBinContent(b + 1, s + 1) for s in range(ns)), 1.0 if total else 0.0)
for s, total in enumerate(reco_sums):
    assert equal(sum(hist['p_truth_given_reco_selected'].GetBinContent(b + 1, s + 1) for b in range(nb)), 1.0 if total else 0.0)
total = sum(reco_sums)
coverage = {}
for name in ('vertex_q2_xb_selected', 'reco_q2_xb_selected'):
    h = hist[name]
    assert h.GetEntries() == total, name
    assert h.Integral() == total, name
    assert h.Integral(0, h.GetNbinsX() + 1, 0, h.GetNbinsY() + 1) == total, name
    coverage[name] = {'events': total, 'xb_range': [h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()],
                      'q2_range': [h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()]}

compared = []
if args.baseline:
    for name in ('migration_design.csv', 'migration_parameters.csv', 'migration_covariance.csv',
                 'migration_reco_rows.csv', 'migration_response_cells.csv', 'migration_truth_blocks.csv'):
        before, after = rows(args.baseline / name), rows(args.output / name)
        assert len(before) == len(after), name
        for a, b in zip(before, after):
            assert a.keys() == b.keys(), name
            for key in a:
                if a[key] == b[key]:
                    continue
                assert equal(float(a[key]), float(b[key])), (name, key, a[key], b[key])
        compared.append(name)
report = {'result': 'PASS', 'response_shape': [nr, nb], 'reco_slices': ns,
          'selected_mc_entries': total, 'truth_support': truth_sums,
          'coverage': coverage, 'unchanged_numerical_exports': compared}
if args.json:
    args.json.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
root.Close()
