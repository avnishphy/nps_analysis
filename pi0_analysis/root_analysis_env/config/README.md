# Configuration Files

The `config/` directory stores run lists, cuts, branch lists, and run metadata
used by the ROOT and Python analysis workflows.

## Run Lists

- `runlist_x60_4a.txt`: x60 4a run list.
- `runlist_x60_4b.txt`: x60 4b run list.
- `runlist_production.txt`: production-data run list.

Run lists should contain one run number per line. Blank lines and comments are
ignored by the wrapper scripts.

## Cuts

- `cuts.conf`: simple key/value cut file.

Current cut categories include:

- HMS PID thresholds
- pi0 invariant-mass window
- vertex placeholders
- skim thresholds

Keep the `KEY = value` format simple so ROOT and shell-side parsers remain
easy to maintain.

## Branch Lists

- `branches_to_read.txt`: branch allow-list used by skim and analysis stages.

This file includes both original skim branch names such as `NPS.cal.*` and HMS
scaler/trigger branches. The waveform-calibrated production workflow maps
several NPS and scaler branches to production equivalents; see
`../README_wfpi0_analysis.md`.

## Run Metadata

- `nps_dvcs_all_kins_main.csv`: main run metadata/configuration table.
- `DataBase_production_runs_newBCMOffset.txt`: production run database with BCM
  offset information.

Python scripts such as `scripts/combine_analysis_branches.py` read the metadata
CSV to determine target, kinematic setting, and prescale information.

## ROOT Cut Artifacts

- `pid_cuts.root`
- `pid_cuts_x60_4a.root`

These are binary cut/config artifacts. Do not edit them manually; regenerate
them with the relevant analysis tooling if their contents need to change.
