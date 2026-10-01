# Final workspace operating procedure

## Boundaries

- Editable workspace: `/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final`.
- Frozen source: `/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main` at `e1e7ecc4b18d359e02140c1785df64b69b3b53aa`.
- Read-only audit: `/work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/docs/pi0_uncertainty_audit_2026-10-01.md`.

Never build, cache, log, validate, or write output in the frozen source. Do not
run a migrated launcher until the worklog shows `WF-003` resolved.

## Setup

```bash
cd /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final
./scripts/setup_final_workspace.sh
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; root-config --version'
```

The setup script creates only ignored runtime directories under this workspace.
The ROOT command should report the environment version without compiling or
running analysis code.

## Stage order

1. Verify baseline/copy/config manifests and the `MAIN` freeze.
2. Classify every active absolute path as immutable input, configuration,
   build/cache, log, or output. Resolve all writes into `FINAL`.
3. Map the active call graph, estimators, units, settings, dependencies, and
   output schemas without changing scientific behavior.
4. Run syntax/help checks and bounded behavior-preservation comparisons on
   identified inputs, writing only to `validation/runtime/` or `scratch/`.
5. Prepare dependency-ordered algorithm proposals with prospective validation
   criteria. Wait for explicit user approval before implementation.
6. Implement one approved increment at a time with recovery copy, focused
   validation, diff review, worklog update, and scoped local commit.
7. Run agreed scientific gates, then a resource-bounded production campaign.
8. Generate reviewed publication artifacts under `publication/generated/` and
   promote only explicitly reviewed release products/manifests.

## Stop conditions

Stop execution and record the evidence if a path resolves into `MAIN`, an
input/config hash differs unexpectedly, a process would overwrite an existing
campaign, a required approval is absent, a covariance/fit/completeness check
fails, or an input/resource budget is not established. Do not repair `MAIN` or
silently relax a scientific gate.
