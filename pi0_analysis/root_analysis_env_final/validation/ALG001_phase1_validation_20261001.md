# ALG-001 phase-one validation — 2026-10-01

Scope: passive raw-observation export and run/segment ledgers only. No real data,
SIMC production, background replacement, extraction fit, or efficiency
calculation was run or changed.

## Approval boundary

The user approved examining and implementing ALG-001 except for efficiency
calculations, which are deferred. `git status` and `git diff` over
`src/efficiencies/` and `config/acceptance_cuts.conf` were empty.

Pre-edit versions and matching SHA-256 values are retained under ignored
`recovery/pre_alg001_phase1_20261001/`.

## Passing checks actually run

Working directory for all commands:

```text
/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final
```

Commands:

```bash
bash -n src/analysis/run_parallel_nps_analysis_main.sh
python3 -m py_compile src/analysis/combine_analysis_branches.py tests/test_raw_observation_bundle.py
python3 tests/test_plot_diagnostics.py
g++ -std=c++17 -Wall -Wextra -pedantic tests/test_nps_raw_observation.cpp -o /tmp/test_nps_raw_observation_20261001
/tmp/test_nps_raw_observation_20261001
python3 tests/test_raw_observation_bundle.py
bash src/analysis/run_parallel_nps_analysis_main.sh --help | rg -n 'raw-observation-export|Usage|no-combine'
```

Results:

- Shell and Python syntax passed.
- Six existing plot/diagnostic tests passed. Two pre-existing Matplotlib
  no-label legend warnings remain.
- Pure C++ timing classification passed all prompt, diagonal, horizontal,
  vertical, full-box, outside, strict-boundary, exclusivity, and mask checks.
- Synthetic bundle test passed: run 100 `ready`, run 101 `zero_candidate`, run
  102 `missing_diagnostics`; combined raw tree had two observations and 25
  branches.
- The same synthetic test compared default and opt-in combination. Their 14
  shared columns were exactly equal; `event_id` was absent by default and was
  the sole extra source branch in opt-in mode.
- Launcher help showed the new flag; no analysis started.
- The launcher requires opt-in runs to name an explicit output base different
  from canonical `FINAL/output`. An optional `--efficiency-csv` points only to
  an existing frozen correction input for a one-setting shadow run.

The guard was executed with `ROOT_CMD=true` (a harmless command used only to
pass the launcher's executable-presence check), run 4398 selection, `--run-only`,
and canonical `FINAL/output`. It exited 1 before job construction or analysis:

```text
--raw-observation-export requires an explicit noncanonical --output-base to prevent overwriting legacy outputs.
```

ROOT compilation command:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; root -l -b -q -e '\''gSystem->SetBuildDir("build/alg001_phase1", kTRUE); gROOT->ProcessLine(".L src/analysis/nps_analysis_main.C+");'\'''
```

Result: ROOT 6.30.04 ACLiC created
`build/alg001_phase1/nps_analysis_main_C.so`; exit status 0. The environment
reported its inherited libstdc++ ABI warning and the macro's existing three
signed/unsigned comparison warnings. No new compile error or warning was
reported.

The exact ROOT `std::string` schema used by `raw_observation_segments` was also
written to `/tmp/nps_alg001_string_branch.root` and read through
`read_tree_dataframe`; the recovered row was:

```text
run_number=1, source_tree_number=0, source_tree_name=T,
source_path=/tmp/input.root
```

## Failed check classified as pre-existing

Command:

```bash
python3 tests/test_2d_mass_cut.py
```

One of two tests failed:

```text
0.16964444444444443 not less than 0.15
```

The ALG-001 patch does not edit the mass-cut implementation. To test the actual
pre-edit behavior, the same test module was run after replacing its `combine`
global with the checksum-verified recovery copy
`recovery/pre_alg001_phase1_20261001/src/analysis/combine_analysis_branches.py`.
It produced the identical value and failure. Therefore this is a documented
baseline test failure, not an ALG-001 regression. The threshold or mass-cut
algorithm was not changed.

One final combined syntax command was mistakenly launched from the enclosing
repository while using FINAL-relative paths. `bash` and `py_compile` reported
`No such file or directory`; no file changed. Repeating both commands from the
documented FINAL working directory passed, and the staged `git diff --check`
was separately rerun from the enclosing repository with the correct scoped
path and passed.

## Remaining required validation

- Default-off versus opt-in comparison on one user-selected real run: legacy
  trees/histograms, selected event IDs/order, summaries, weights, and logs.
- Reconcile raw timing-category counts to source histograms by run and
  multiplicity.
- Validate a complete run ledger for one setting, including intended zero-event
  and missing/partial statuses.
- Record input/output checksums and performance.

No production-readiness or publication claim is made until those checks and the
later background/smearing/extraction proposals pass.
