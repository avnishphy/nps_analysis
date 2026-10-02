# Background-fit shadow tools

This directory contains opt-in statistical pilots isolated from production.
ALG-002A fits combined timing statistics while retaining run,
acquisition-mode, and multiplicity identity.

It does **not** write `pi0_weight`, read or calculate efficiencies, modify the
legacy timing subtraction, form a cross section, or write below `output/`.

Entry points:

- `run_joint_timing_fit.py`: fit ALG-001 `raw_observation` trees;
- `validate_joint_timing_output.py`: write run coverage and explicit
  pass/fail/pending scientific gates;
- `joint_timing_model.py`: model, accepted geometry, fit, and output code.

See `docs/ALG002A_combined_timing_background_20261001.md` for the full model,
reasoning, commands, output schema, validation history, and current blockers.

