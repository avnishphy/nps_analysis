#!/usr/bin/env python3
"""Read-only C++/ROOT versus NumPy regression on KinC_x60_4b runs 4253/4260/4261.

Run after sourcing the Hall C ROOT environment. Compile artifacts go to /tmp.
Usage: python3 tests/test_dummy_sub.py [--data-base /path/to/workflow]
No production or combined ROOT files are written; all stored physics entries
are used, so this validates subtraction arithmetic, not a final exclusive yield.
"""
import argparse
import csv
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np
import uproot

REPO = Path(__file__).resolve().parents[1]
RUNS = (4253, 4260, 4261)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-base', type=Path, default=REPO)
    args = parser.parse_args()
    with (args.data_base / 'output/efficiency_stuff/efficiency_KinC_x60_4b.csv').open() as f:
        rows = {int(r['run_number']): r for r in csv.DictReader(f)}
    charges = []
    for run in RUNS:
        row = rows[run]
        assert row['run_processing_status'].strip() == 'processed'
        # All three fixture runs have ps4=0 or ps6=0, hence factor 1.
        assert float(row['ps_factor']) == 1.0
        q = (float(row['HEL_charge_after_cut_uC']) / 1000
             * float(row['NewGen_EDTM_livetime'])
             * float(row['HMS_tracking_eff'])
             * float(row['HMS_hodo_3of4_eff']))
        assert np.isfinite(q) and q > 0
        charges.append(q)

    with tempfile.TemporaryDirectory(prefix='nps_dummy_sub_') as tmp:
        exe = str(Path(tmp) / 'test_nps_dummy_sub')
        flags = shlex.split(subprocess.check_output(
            ['root-config', '--cflags', '--libs'], text=True))
        subprocess.run(['c++', '-Wall', '-Wextra', '-pedantic',
                        str(REPO / 'tests/test_nps_dummy_sub.cpp'),
                        '-o', exe] + flags, check=True)
        output = subprocess.check_output(
            [exe, str(args.data_base), *map(str, charges), str(RUNS[0])], text=True)
        print(output, end='')
    observed = np.array([float(x) for x in output.split('RESULT ')[1].split()])
    lh2 = background = variance = 0.0
    for run in RUNS:
        path = args.data_base / f'output/KinC_x60_4b/root/diagnostics_run{run}.root'
        # Chunked reads bound memory even for large production diagnostics.
        for a in uproot.iterate(f'{path}:physics', ['ytar', 'pi0_weight'],
                               library='np', step_size='20 MB'):
            y, w = a['ytar'], a['pi0_weight']
            assert np.isfinite(y).all() and np.isfinite(w).all()
            if run == RUNS[0]:
                weights = w / charges[0]
                lh2 += weights.sum()
            else:
                weights = w * np.where(y < 0, 1 / 8.467,
                                       np.where(y > 0, 1 / 4.256, 0)) / sum(charges[1:])
                background += weights.sum()
            variance += np.square(weights).sum()
    expected = [lh2, background, lh2 - background, variance]
    np.testing.assert_allclose(observed, expected, rtol=1e-10, atol=1e-10)
    print(f'PASS ROOT/NumPy: runs={RUNS}, Qeff_mC={charges}, '
          f'stat_error={np.sqrt(variance):.9g}')


if __name__ == '__main__':
    main()
