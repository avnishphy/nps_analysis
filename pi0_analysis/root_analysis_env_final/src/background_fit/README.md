# Background-fit shadow tools

This directory contains opt-in statistical pilots isolated from production.
ALG-002A fits combined timing statistics while retaining run,
acquisition-mode, and multiplicity identity. ALG-002B adds an approved
simultaneous timing/mass shadow fitter restricted to configured production-LH2
runs in `KinC_x36_4`.

It does **not** write `pi0_weight`, read or calculate efficiencies, modify the
legacy timing subtraction, form a cross section, or write below `output/`.

Entry points:

- `run_joint_timing_fit.py`: fit ALG-001 `raw_observation` trees;
- `validate_joint_timing_output.py`: write run coverage and explicit
  pass/fail/pending scientific gates;
- `joint_timing_model.py`: model, accepted geometry, fit, and output code.
- `run_joint_timing_mass_fit.py`: run the ALG-002B six-component shadow fit;
- `validate_joint_timing_mass_output.py`: write fail-closed ALG-002B gates;
- `joint_timing_mass_model.py`: joint timing/mass model and isolated outputs.

See `docs/ALG002A_combined_timing_background_20261001.md` for the full model,
reasoning, commands, output schema, validation history, and current blockers.
See `docs/ALG002B_simultaneous_timing_mass_20261001.md` for the ALG-002B
implementation, strict LH2 manifest, smoke result, and open validation gates.
