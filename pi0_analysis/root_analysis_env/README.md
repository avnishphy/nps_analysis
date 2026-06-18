# NPS Pi0 ROOT Analysis Environment

This workspace contains ROOT/C++, Python, and shell tools for Hall C NPS pi0
analysis. It covers per-run diagnostics, waveform-calibrated production data,
branch combining, SIMC smearing, data/simulation comparisons, and cross-section
extraction.

## Layout

- `src/`: core ROOT analysis macros and helper headers.
- `scripts/`: smearing, SIMC, cross-section, plotting, diagnostics, and
  efficiency tools.
- `config/`: run lists, branch lists, cuts, and run metadata.
- `output/`: generated per-run diagnostics and combined analysis artifacts.
- `output_pi0_xsec/`: generated cross-section plots and summaries.
- `output_pi0_xsec_no_simc_model/`: generated cross-section outputs without
  the SIMC model.
- `tmp_runlists/`: temporary run-list chunks used by parallel jobs.

## Main Workflows

### Original Skim-Based Analysis

```bash
./run_nps_analysis.sh \
  config/runlist_x60_4b.txt \
  /lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/root_analysis_env_skim/x60_4b \
  output/plots/x60_4b \
  10.538
```

This wrapper runs `src/nps_analysis.C` one run at a time, captures the per-run
CSV row emitted by the macro, and rebuilds `summary_all_runs.csv`.

### Waveform-Calibrated Production Analysis

```bash
./run_nps_analysis_wfpi0.sh \
  config/runlist_x60_4b.txt \
  /lustre24/expphy/volatile/hallc/nps/hhuang/farmFile/Production/DVCS \
  output/plots/x60_4b/production_wfpi0 \
  10.538 \
  1800
```

This wrapper runs `src/nps_analysis_wfpi0.C`, which expects production files
with the `t_prod` tree and `NPS.prod.*` branch names.

Key differences from `src/nps_analysis.C`:

- Tree name: `T` -> `t_prod`
- NPS cluster branches: `NPS.cal.*` -> `NPS.prod.*`
- Production files are split by run and loaded with `TChain`
- Production file pattern:

```text
nps_production_<run>_<split>_wf_calib.root
```

Important branch mappings:

```text
NPS.cal.nclust   -> NPS.prod.nclust
NPS.cal.clusE    -> NPS.prod.clusE
NPS.cal.clusX    -> NPS.prod.clusX
NPS.cal.clusY    -> NPS.prod.clusY
NPS.cal.clusT    -> NPS.prod.clusT
H.BCM2.*         -> H.BCM4A.*
H.1MHz.scalerTime -> H.1MHz.scaler
```

Environment overrides supported by both analysis wrappers:

```bash
export NPS_SKIM_DIR="/path/to/input/files"
export NPS_OUTPUT_DIR="output/custom_plots"
export NPS_RUNLIST="config/custom_runlist.txt"
export NPS_EBEAM="10.538"
```

### Parallel Production Analysis

```bash
./run_parallel_nps_analysis_wfpi0.sh \
  config/runlist_x60_4b.txt \
  /lustre24/expphy/volatile/hallc/nps/hhuang/farmFile/Production/DVCS \
  output/plots/x60_4b/production_wfpi0 \
  10.538 \
  1800 \
  8
```

Use a conservative job count when writing to shared filesystems.

### Smearing Pipeline

```bash
./run_smearing_pipeline.sh
```

The pipeline compiles `scripts/nps_sim_smearing_new_try.C`, fits section-wise
smearing parameters, and regenerates smeared SIMC output via
`scripts/simc_pi0_analysis.C`. The smearing fit uses
`is_exclusive_ellipse_combined` for combined data and `is_exclusive_ellipse`
for SIMC by default; regenerate combined data and unsmeared SIMC after changing
the 2D mass-cut logic.

Key diagnostics are written beside `section_map.csv`, including
`smearing_sweep_history.csv`, `smearing_objective_breakdown.csv`,
`smearing_closure_summary.csv`, `smearing_runtime_summary.csv`,
`best_parameters.json`, `optimization_config.json`, and the optimizer
seed/profile summaries. The coupled-sweep optimizer uses an explicit
candidate/current-accepted/best-ever state model, so rejected sweeps are not
allowed to seed later sweeps.

Current validation status is reduced-validation only. Peak-window agreement
improves in the latest reduced tests, but full-width/tail behavior still needs
production-scale physics sign-off; see `smearing.md` for details.

## Configuration

- `config/runlist_x60_4a.txt`, `config/runlist_x60_4b.txt`, and
  `config/runlist_production.txt` list runs, one run number per line.
- `config/cuts.conf` stores simple key/value cuts such as PID and pi0 mass
  windows.
- `config/branches_to_read.txt` lists branches used during skim/analysis.
- `config/nps_dvcs_all_kins_main.csv` stores run metadata used by combining and
  plotting scripts.

## Output Physics Tree

Per-run diagnostics ROOT files contain a `physics` tree with event-level
variables. Important branches include:

```text
mpi0_all
mmiss_all
mmiss_all_corr
pi0_weight
is_exclusive
is_exclusive_ellipse
is_exclusive_mcd
is_weighted
Q2, W, t, tmin, pt, theta, phi, s, xB, z
cluster_x_1, cluster_y_1, cluster_e_1
cluster_x_2, cluster_y_2, cluster_e_2
```

`is_exclusive` is the legacy corrected-missing-mass window flag. The two new
2D mass-cut flags do not remove events from the tree:

```text
is_exclusive_ellipse = 1 if event passes covariance-ellipse cut
is_exclusive_mcd     = 1 if event passes robust MCD ellipse cut
```

Use these branches downstream to compare exclusivity definitions.

## 2D Mass-Cut Exclusivity

The shared implementation is:

```text
src/nps_2d_mass_cut.h
```

The test macro is:

```text
scripts/test_2D_mass_cut.C
```

The cut is built from the weighted `mmiss_all:mpi0_all` distribution:

```text
x = mpi0_all
y = mmiss_all
w = pi0_weight
```

The algorithm first builds a connected high-density core around the maximum
2D bin. In auto mode, `peak_fraction <= 0`, it scans thresholds and chooses the
last compact core before a sideband joins the peak component. The scan outputs
`*_peak_scan.csv` for auditing.

Two smooth branch-level selectors are then evaluated:

1. Standard covariance/Mahalanobis ellipse
2. Robust minimum-covariance-distance (MCD) ellipse

Both use the same acceptance form:

```text
d2 = [x - mu]^T C^-1 [x - mu]
accept if d2 <= d2_cut
```

The MCD method recomputes a robust covariance from the closest weighted subset,
making it less sensitive to sideband/tail bins.

Default debug outputs are written beside analysis plots:

```text
mass_cut_run<run>.root
mass_cut_run<run>.pdf
mass_cut_run<run>.png
mass_cut_run<run>_params.csv
mass_cut_run<run>_peak_scan.csv
```

For systematic-uncertainty details and suggested variations, see:

```text
systematic_uncertainty.md
```

## SIMC Output

`scripts/simc_pi0_analysis.C` writes the same branch names to unsmeared and
smeared simulation trees:

```text
is_exclusive
is_exclusive_ellipse
is_exclusive_mcd
```

`is_exclusive` is the input-sample label (`exclusive` vs `sidis`). The ellipse
and MCD parameters are trained on the `is_exclusive == 1` simulation rows, then
applied to the full output tree using `mpi0`, `mmiss`, and `full_weight`. This
keeps the SIDIS sideband from biasing the learned exclusive region.

The default SIMC inputs can be overridden with:

```bash
export NPS_SIMC_EXCLUSIVE_INPUT="/path/to/exclusive.root"
export NPS_SIMC_SIDIS_INPUT="/path/to/sidis.root"
```

## Generated Outputs

Typical analysis outputs include:

- `diagnostics_run<run>.root`
- per-run diagnostic plots
- `summary_run<run>.txt`
- `summary_all_runs.csv`
- combined branch ROOT files
- cross-section summary CSV files and publication-style PDFs/PNGs

Most generated artifacts are intentionally not good targets for broad edits.
Inspect them when needed, but prefer changing the producing macro/script.

## Requirements

- ROOT and `root-config`
- Bash utilities: `grep`, `sed`, `awk`, `timeout`
- Python packages for selected scripts: `numpy`, `pandas`, `matplotlib`,
  `uproot`
- A Hall C filesystem environment with access to the referenced input ROOT
  files

Set up ROOT in the Hall C/NPS environment with:

```csh
source /group/nps/singhav/setup.csh
```

For non-interactive `tcsh`/`csh` sessions, initialize modules first:

```csh
source /usr/share/Modules/init/csh
source /group/nps/singhav/setup.csh
```

## Notes For Future Work

- Keep the `[CSV_WRITTEN]` marker in ROOT macro output unless wrapper scripts
  are updated.
- Many paths are absolute and specific to the current analysis area.
- Before running full production jobs, test with a one-run run list when
  possible.
- See `AGENTS.md` for agent-specific workflow and safety guidance.
