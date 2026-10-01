# Smearing audit fixes, 2026-09-11

Scope: correct the reproduced failures in the post-update audit without changing response formulas, event selection, weighting formulas, fitting strategy, random draws, or interpolation policy.

## Producer input failures

All input branch bindings in the producer now check ROOT's return status and verify the branch entry count against the tree. Incompatible types, missing/incomplete required branches, and failed event reads return nonzero status with the input sample and offending branch/entry. Optional missing `sigcm` for non-exclusive samples retains its existing policy.

Regression: valid synthetic inputs produce 120 validated events with expected `full_weight` values. The former `Float_t`/`Double_t Weight` mismatch now returns status 7 instead of publishing default weights. A 39-entry Weight branch in a 40-entry tree also returns status 7. No new physics cuts or type-coercion policy were introduced.

## Fast evaluator

The bin lookup now follows ROOT's floating-point operation order. All 549 configured-axis boundary and adjacent-value probes match ROOT (the audited version had 18 mismatches). Negative cache budgets, including whitespace-prefixed negative values, are rejected; valid zero still selects the uncached path. Thirteen parser cases pass.

Global and final-section diagnostics only perform the expensive legacy reference evaluation when `VALIDATE_FAST_OBJECTIVE` is enabled. Its default remains enabled, so default validation tolerances/fallback are preserved. This removes redundant work when validation is explicitly disabled; it does not claim a default production speedup.

The existing broader scalar/batch/legacy equivalence harness also passes against the staged changes, including histogram contents/sumw2, cache fallback and map checks. Optional inactive physics modes were not newly toggled or redesigned.

## Input provenance and publication

Raw inputs, normalization files and configuration now use content SHA-256 identities rather than size/integer-second timestamps. Small configuration/normalization inputs are frozen for both producer passes and retained with the archive. Build provenance records dependency, executable and toolchain information; dependency hashes are compared before and after compilation.

Content verification adds full sequential reads of raw ROOT inputs at consistency checkpoints. This deliberately prioritizes verified input identity over the previous cheap stat check; no production I/O speedup is claimed.

The complete validated archive is committed before replacing legacy output files. Publication is serialized, uses unique temporary files and keeps rollback backups/journals. Ordinary publication failure restores previous legacy outputs while preserving the completed archive. Interrupted publication can be recovered automatically before the next publication or explicitly with:

```bash
python3 src/simulation_smearing/publish_smearing_run.py \
  --recover-only --manifest /absolute/output/directory/smearing_latest.json
```

`smearing_latest.json` is updated only after successful publication and identifies the complete archived set. Legacy filenames are still replaced individually: concurrent readers needing a coherent multi-file view must read this manifest once and use its immutable archive paths. No claim of simultaneous atomic replacement of multiple legacy paths is made.

Publication also guards input/source path collisions and conflicting ownership of output destinations. Recovery refuses unsafe cross-manifest reuse rather than risking restoration over another run's results. Recovery state should not be manually deleted while publication remains unresolved.

## Reproducible bounded checks

From the repository root, load Hall C ROOT before compilation/execution. These tests create only temporary outputs:

Validation passed: eight publication/provenance tests, shell syntax, producer input checks, and evaluator boundary/cache-budget tests. The producer and fitter preprocessor/compiler dependency sets also matched under ROOT 6.30.04.

```bash
bash -n src/simulation_smearing/run_smearing_pipeline.sh
python3 tests/test_smearing_publication.py
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ tests/smearing_fast_regression.C `root-config --cflags --libs` -lMathMore -O2 -std=c++17 -fopenmp -o /tmp/smearing_fast_regression && /tmp/smearing_fast_regression'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ src/simulation_smearing/simc_pi0_analysis.C `root-config --cflags --libs` -O2 -std=c++17 -o /tmp/simc_input_check && g++ tests/simc_input_fixture.C `root-config --cflags --libs` -O2 -std=c++17 -o /tmp/simc_input_fixture && python3 tests/test_simc_input_validation.py /tmp/simc_input_check /tmp/simc_input_fixture .'
```

No full production campaign or new cross-section extraction was run. The changes repair failure handling and numerical boundary semantics; they do not recalibrate existing maps or establish new physics results.
