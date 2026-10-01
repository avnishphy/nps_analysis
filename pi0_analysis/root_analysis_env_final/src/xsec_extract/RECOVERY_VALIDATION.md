# Global Q2/xB subset recovery validation (2026-09-17)

Changes were prepared under `/tmp/xsec_bin_recovery/src/xsec_extract`, using
the existing `root_analysis_env_main` files as the baseline. Production inputs
were opened read-only; production outputs were not overwritten.

## Checks completed

- ROOT 6.30.04 executable build; standalone recovery tests and existing linear
  solver tests passed. Recovery tests cover a healthy full fit, missing support
  in a later t' slice, retained cross-bin covariance, excluded truth feed-in as
  free nuisance, rank failure, missing truth support, unchanged target scaling
  across retries, convergence exhaustion, and all-failed CSV/ROOT exports.
- PARTONS-enabled translation unit passed `g++ -fsyntax-only` with the native
  include paths. Native PARTONS model evaluation was not rerun.
- Partial/all-failed plot smoke passed; `pdftotext` confirms the excluded-bin
  message. The all-failed combined PDF has 8 readable pages. ROOT still emits
  the existing combined-PDF object-state warnings described in VALIDATION.md.
- Shell runner syntax and patch whitespace checks passed.

Reproduce the standalone checks from `root_analysis_env_main`:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash tests/run_xsec_bin_recovery_tests.sh /tmp/nps_xsec_recovery_check'
```

The test prints `PASS` for response recovery and CSV/ROOT contracts. It uses
manufactured response matrices, not production event inputs.

## User's KinC_x60_4b input

The ROOT-only executable was run with the user's data, volatile smeared SIMC,
raw vertex SIMC, and missing-mass window. PARTONS and plots were disabled for
this event-level check. Exact executable arguments used:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; /tmp/xsec_bin_recovery/extractor --data-file /work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/KinC_x60_4b/root/combined_branches_LH2.root --sim-file /volatile/hallc/nps/singhav/nps_smearing/smear_x60_4b_test/smearing_output/KinC_x60_4b/root/simc_pi0_analysis_output_smeared.root --vertex_simc_file /work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x60_4b.root --out-dir /tmp/xsec_bin_recovery/real_input --mmiss-lower 0.6 --mmiss-upper 1.1 --no-pdf --no-png --no-diagnostics'
```

Observed result:

- First fit fails at row 106: `it=2, iq=0, ix=0, ip=10`.
- Excludes `(iq=0, ix=0)` and fits the three remaining Q2/xB groups jointly.
- 175 fitted rows, 72 parameters, rank 72, scaled condition 794.455127.
- Chi2/ndf = 94.0839/103; finite-MC variance converges in 13 iterations.
- Seven retained generated t'/Q2/xB bins have negative angular regions and
  remain flagged. This recovery is not a physical-validation claim.

Independent NumPy reproduction was run with the existing validator:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; python3 /work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/tests/validate_migration.py /tmp/xsec_bin_recovery/real_input --json /tmp/xsec_bin_recovery/real_input_validation.json'
```

It reproduced parameters to relative error `2.85e-14` and covariance to
`5.06e-14`. Target covariance/scaling, full off-diagonal covariance, manufactured
algebraic closure and rank rejection checks passed. No new independent-event
or selection-bias physics closure is claimed.
