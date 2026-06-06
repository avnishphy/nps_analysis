# Analysis Scripts

The `scripts/` directory contains secondary analysis tools: smearing,
simulation processing, cross-section extraction, plotting, diagnostics, and
efficiency studies.

## Smearing And SIMC

- `nps_sim_smearing_new.C`: C++ source for section-wise smearing fits.
- `nps_sim_smearing_new`: compiled executable generated from the source.
- `simc_pi0_analysis.C`: SIMC pi0 analysis macro and smearing consumer.
- `section_map.csv`: detector section map used by section smearing.
- `TWO_STAGE_SMEARING_EXPLANATION.md`: detailed two-stage smearing notes.
- `NPS_Smearing_Methodology.txt`: additional methodology notes.

The top-level `run_smearing_pipeline.sh` compiles and runs the main smearing
workflow.

## Data/SIMC Comparison And Combining

- `combine_analysis_branches.py`: combines per-run ROOT diagnostic trees into a
  single weighted branch dataset.
- `data_simulation_comparison.py`: data/SIMC comparison plotting.
- `plot_publication_quality.py`: publication-oriented plotting.
- `plot_sigma_coeff_publication.py`: sigma coefficient plotting.
- `create_diagnostic_pdfs.py`: diagnostic PDF generation helper.

## Cross-Section Extraction

- `excl_xsec_pi0_analysis.C`: primary exclusive pi0 cross-section extraction.
- `excl_xsec_pi0_analysis_no_simc_model.C`: variant without the SIMC model.
- `excl_xsec_pi0_analysis_bak_working.C`: retained backup/reference version.
- `excl_xsec_simc_geant.C`: SIMC/GEANT-related cross-section support.

## Efficiency And Livetime Tools

Efficiency-related files live under `scripts/efficiencies/`, including scaler,
livetime, and HMS tracking-efficiency utilities.

## Diagnostics And Notes

The text files in this directory capture previous debugging context:

- `MISSING_MASS_FITTING_ISSUE.txt`
- `SYSTEMATIC_PEAK_SHIFTS_DIAGNOSIS.txt`
- `CONFIGURATION_UPDATES.txt`

Read these before changing related fit or calibration behavior; they contain
analysis-specific breadcrumbs that are easy to lose in code alone.

## Development Notes

- Many scripts use absolute paths into the current Hall C analysis area. Check
  those constants before running from a different environment.
- Keep generated executables, ROOT files, plots, and notebook outputs out of
  code reviews unless the task explicitly concerns those artifacts.
- For Python scripts, prefer making path/config changes near the top of the
  file so batch runs are easy to audit.
