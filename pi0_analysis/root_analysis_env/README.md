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
with the `t_prod` tree and `NPS.prod.*` branch names. See
`README_wfpi0_analysis.md` for the detailed branch mapping.

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

The pipeline compiles `scripts/nps_sim_smearing_new.C`, fits section-wise
smearing parameters, and regenerates smeared SIMC output via
`scripts/simc_pi0_analysis.C`.

## Configuration

- `config/runlist_x60_4a.txt`, `config/runlist_x60_4b.txt`, and
  `config/runlist_production.txt` list runs, one run number per line.
- `config/cuts.conf` stores simple key/value cuts such as PID and pi0 mass
  windows.
- `config/branches_to_read.txt` lists branches used during skim/analysis.
- `config/nps_dvcs_all_kins_main.csv` stores run metadata used by combining and
  plotting scripts.

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

## Notes For Future Work

- Keep the `[CSV_WRITTEN]` marker in ROOT macro output unless wrapper scripts
  are updated.
- Many paths are absolute and specific to the current analysis area.
- Before running full production jobs, test with a one-run run list when
  possible.
- See `AGENTS.md` for agent-specific workflow and safety guidance.
