#!/usr/bin/env python3
"""Verify new event-cache normalization against a saved same-binning extraction."""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src/xsec_extract'))
from forward_xsec_problem import build_problem


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    args=parser.parse_args()
    problem=build_problem(args.cache, args.config)
    pars=read(args.reference/'migration_parameters.csv')
    rows=read(args.reference/'migration_reco_rows.csv')
    assert np.array_equal(problem['active_global_blocks'], [int(p['truth_block']) for p in pars[::3]])
    assert problem['design'].shape==(len(rows),len(pars))
    matrix=np.zeros_like(problem['design'])
    for row in read(args.reference/'migration_design.csv'):
        matrix[int(row['reco_row']),int(row['parameter_index'])]=float(row['response'])*1e-9
    report={}
    for name,actual,expected in (
            ('response',problem['design'],matrix),
            ('data_yield',problem['y'],np.array([float(r['data']) for r in rows])),
            ('data_sumw2',problem['sumw2'],np.array([float(r['data_variance']) for r in rows]))):
        absolute=float(np.max(np.abs(actual-expected)))
        scale=float(np.max(np.abs(expected)))
        relative=absolute/scale if scale else absolute
        report[name]={'max_absolute_difference':absolute,'relative_to_global_max':relative}
        if not relative < 1e-11:
            raise AssertionError(f'{name} disagrees with the saved reference: {relative}')
    report['status']='PASS same-response normalization and observation regression; not physics closure'
    rendered=json.dumps(report,indent=2)+'\n'
    if args.output:
        with args.output.open('x') as stream: stream.write(rendered)
    print(rendered,end='')


if __name__=='__main__': main()
