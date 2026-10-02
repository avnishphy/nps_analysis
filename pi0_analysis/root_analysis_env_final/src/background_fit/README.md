# Background-fit shadow tools

This directory contains opt-in statistical pilots isolated from production.
ALG-002A fits combined timing statistics while retaining run,
acquisition-mode, and multiplicity identity. ALG-002B adds an approved
simultaneous timing/mass shadow fitter restricted to configured production-LH2
runs in `KinC_x36_4`.

It does **not** write `pi0_weight`, read or calculate efficiencies, modify the
legacy timing subtraction, or form a cross section. Output below canonical
`output/` requires an explicit flag and is restricted to the isolated
`output/KinC_x36_4/alg002b/` shadow subtree.

Install the small mass-gradient dependency into the managed analysis Python:

```bash
/group/nps/singhav/software/python/bin/python -m pip install \
  -r src/background_fit/requirements.txt
```

The default `autograd` mass backend differentiates the same profiled likelihood
used for final evaluation. `--mass-gradient-backend finite` retains the slower
finite-difference audit path. Timing gradients continue to use `--nproc` CPU
workers.

Entry points:

- `run_joint_timing_fit.py`: fit ALG-001 `raw_observation` trees;
- `validate_joint_timing_output.py`: write run coverage and explicit
  pass/fail/pending scientific gates;
- `joint_timing_model.py`: model, accepted geometry, fit, and output code.
- `run_joint_timing_mass_fit.py`: run the ALG-002B six-component shadow fit;
- `validate_joint_timing_mass_output.py`: write fail-closed ALG-002B gates;
- `joint_timing_mass_model.py`: joint timing/mass model and isolated outputs.
- `compare_legacy_timing_background.py`: report the original production
  timing-box estimate beside ALG-002B without using it as a fit input.
- `compare_raw_observation_inputs.py`: prove two diagnostic bundles have
  identical `raw_observation` entries and branch buffers before reusing a fit
  initialization.
- `run_joint_timing_mass_campaign.py`: schedule independently dispersed,
  resumable starts and aggregate the approved 20-start reproducibility gate.

Pass the completed campaign directory to
`validate_joint_timing_mass_output.py --start-campaign <directory>` so the
otherwise fail-closed reproducibility gate uses its aggregate result.

See `docs/ALG002A_combined_timing_background_20261001.md` for the full model,
reasoning, commands, output schema, validation history, and current blockers.
See `docs/ALG002B_simultaneous_timing_mass_20261001.md` for the ALG-002B
implementation, strict LH2 manifest, smoke result, and open validation gates.
