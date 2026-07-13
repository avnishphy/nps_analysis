# Repository Guidelines

## Project Structure & Module Organization

This Hall C NPS analysis workspace uses ROOT/C++, Python, shell scripts, and notebooks. Top-level macros and branch lists support quick checks. Domain work is grouped by directory: `pi0_analysis/` for pi0 skims and the nested `root_analysis_env/` production workflow; `elastic_analysis/` for elastic skims, optics, efficiency, and SIMC comparisons; `luminosity_analysis/` for scaler calculations; `det_calibrations/` for detector calibration; `efficiencies/` for report extraction; `bash_scripts/` for utilities; and `swif2/` for batch job generation. Generated plots, ROOT files, summaries, and temporary chunks usually live under `output/`, `*_plots/`, or `tmp_*` paths.

## Build, Test, and Development Commands

There is no top-level build system. Run ROOT macros after loading the Hall C/NPS ROOT environment:

```bash
csh -c 'source /group/nps/singhav/setup.csh; root -l -b -q "Scaler_analysis.C"'
csh -c 'source /group/nps/singhav/setup.csh; root -l -b -q "elastic_analysis/skim_elastic_data.cpp"'
```

Common workflow scripts:

```bash
cd pi0_analysis/root_analysis_env && ./run_nps_analysis.sh
cd pi0_analysis/root_analysis_env && ./run_smearing_pipeline.sh
cd swif2 && ./make_script_runlist.sh 6810 6810 run-lists/runlist_1.dat
```

Use conservative parallel job counts. Do not run production workflows for documentation-only changes.

## Coding Style & Naming Conventions

Match the local style in the file you edit. ROOT/C++ files use `.C`, `.cxx`, or `.cpp`; keep variable names descriptive and physics-facing, such as `mmiss_all`, `pi0_weight`, or `current_mean_uA`. Shell scripts should be executable, start with a shebang, and quote paths and variables. Python scripts should use snake_case. Keep run lists and configs as plain text with one run or key/value setting per line.

## Testing Guidelines

Validation is workflow-specific. For ROOT changes, run the smallest relevant macro or a single-run wrapper before scaling up. For pi0 production work, use `pi0_analysis/root_analysis_env/AGENTS.md` and `README.md`. For Python utilities, run a small known input and confirm expected CSV, plot, or summary output. Avoid committing bulk regenerated outputs unless the change explicitly concerns them.

## Commit & Pull Request Guidelines

Recent history uses conventional-style messages, such as `fix(smearing): complete global-prefit workflow`, and plain descriptive summaries. Prefer short imperative subjects with an optional scope, for example `fix(luminosity): guard missing scaler rows`. Pull requests should describe the area touched, commands run, sample runs used, and generated artifacts intentionally updated. Include representative plots or CSV summaries for physics-output changes.

## Security & Configuration Tips

Do not hard-code private credentials or machine-specific scratch paths in reusable scripts. Prefer environment variables such as `NPS_SKIM_DIR`, `NPS_OUTPUT_DIR`, `NPS_RUNLIST`, and `NPS_EBEAM` when available. Treat large ROOT files, production outputs, and batch logs as reproducible artifacts unless explicitly requested.
