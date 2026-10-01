# ALG-001 phase-one validation — 2026-10-01

Scope: passive raw-observation export and run/segment ledgers, including a
matched default-off/opt-in check on real-data run 5237. No SIMC production,
background replacement, extraction fit, or efficiency calculation was run or
changed.

## Approval boundary

The user approved examining and implementing ALG-001 except for efficiency
calculations, which are deferred. `git status` and `git diff` over
`src/efficiencies/` and `config/acceptance_cuts.conf` were empty.

Pre-edit versions and matching SHA-256 values are retained under ignored
`recovery/pre_alg001_phase1_20261001/`. The pre-real-run-fix macro and edited
documentation are retained under ignored
`recovery/pre_alg001_realrun_fix_20261001/`.

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

## Real-run default-off/opt-in comparison

Run 5237 (`KinC_x50_0a`, LH2, HCANA) was selected because it is marked good
production data, has an existing good-event selection row, and has a single
updated replay segment. Both runs used one worker, a 1800 s timeout, the same
configuration and the same read-only selection report. Outputs were isolated
under ignored `validation/runtime/`; no canonical output was overwritten.

The first default-off attempt exited wrapper status 92 before event processing
because FINAL intentionally does not contain copied efficiency products and
the optional good-event cut looked for a FINAL-local selection report. The
retry set the already-supported `NPS_SELECTION_REPORT_CSV` to the existing
frozen MAIN report. That report was read only; no efficiency was calculated.

Post-fix commands actually run were:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; setenv NPS_SELECTION_REPORT_CSV /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/efficiency_stuff/selection_report_KinC_x50_0a.csv; bash src/analysis/run_parallel_nps_analysis_main.sh --source updated --mode hcana --types production --kin KinC_x50_0a --target LH2 --run 5237 --jobs 1 --run-only --gevnum-cut yes --timeout 1800 --output-base validation/runtime/alg001_run5237_postfix_default'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; setenv NPS_SELECTION_REPORT_CSV /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/efficiency_stuff/selection_report_KinC_x50_0a.csv; bash src/analysis/run_parallel_nps_analysis_main.sh --source updated --mode hcana --types production --kin KinC_x50_0a --target LH2 --run 5237 --jobs 1 --run-only --gevnum-cut yes --timeout 1800 --raw-observation-export --output-base validation/runtime/alg001_run5237_postfix_optin'
python3 tests/compare_alg001_real_run.py --default-dir validation/runtime/alg001_run5237_postfix_default/KinC_x50_0a --optin-dir validation/runtime/alg001_run5237_postfix_optin/KinC_x50_0a --run 5237 --expected-source /lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_5237_0_1_-1.root
mkdir -p validation/runtime/alg001_run5237_pdf_compare/default validation/runtime/alg001_run5237_pdf_compare/optin
pdftoppm -png -r 100 validation/runtime/alg001_run5237_postfix_default/KinC_x50_0a/plots/cut_debug_run5237.pdf validation/runtime/alg001_run5237_pdf_compare/default/cut_debug
pdftoppm -png -r 100 validation/runtime/alg001_run5237_postfix_optin/KinC_x50_0a/plots/cut_debug_run5237.pdf validation/runtime/alg001_run5237_pdf_compare/optin/cut_debug
pdftoppm -png -r 100 validation/runtime/alg001_run5237_postfix_default/KinC_x50_0a/plots/mass_cut_run5237.pdf validation/runtime/alg001_run5237_pdf_compare/default/mass_cut
pdftoppm -png -r 100 validation/runtime/alg001_run5237_postfix_optin/KinC_x50_0a/plots/mass_cut_run5237.pdf validation/runtime/alg001_run5237_pdf_compare/optin/mass_cut
diff -rq validation/runtime/alg001_run5237_pdf_compare/default validation/runtime/alg001_run5237_pdf_compare/optin
```

### Defect found and fixed

The first opt-in output exposed a false segment row for nonexistent
`skim_run5237.root`. ROOT `TChain::Add` accepts a nonexistent exact filename as
a deferred element, so the real replay appeared as segment 1 even though its
event entries were not duplicated. `open_chain_for_run` now skips nonexistent
non-wildcard paths before calling `Add`; wildcard discovery is unchanged. A
direct resolver assertion then returned:

```text
OPEN_OK=1 NFILES=1 TREE=T ENTRIES=653966 ERROR=
```

ROOT ACLiC compilation, the pure C++ timing test, and the synthetic bundle test
passed after the fix. The same three pre-existing signed/unsigned warnings and
the environment's libstdc++ ABI warning remain.

### Matched-output result

`tests/compare_alg001_real_run.py` returned `ALG001_REALRUN_COMPARE=PASS` and
established:

- all 39 `physics` branches and all 121 rows were exactly equal and ordered;
- 183 shared histograms, eight `TParameter<double>` values, two `TNamed`
  objects, and 16 non-canvas auxiliary ROOT objects were exactly equal;
- 23 ordinary outputs, including PNGs, summaries, and diagnostic CSVs, were
  byte-identical;
- eight physics-relevant log facts (source count/tree/mode, timing cuts, entry
  count, helicity source/charge, dead blocks, weighted-event total, and 2D
  mass-cut result) were identical; expected path, opt-in notice, process ID,
  and runtime differences were not treated as physics content;
- separately rendered pages from both differing binary PDFs were
  byte-identical at 100 dpi;
- the only new diagnostic keys were `raw_observation` and
  `raw_observation_segments`;
- all 121 raw rows matched the `physics.event_id` order and shared observable
  values; event IDs spanned 13735 through 651678;
- the segment ledger had exactly one row, segment 0, naming the existing replay
  file, and every raw row referenced it with an in-range source entry and an
  event number exactly equal to source `T.g.evnum`;
- timing categories were outside 28, prompt 54, diagonal 12, horizontal 12,
  vertical 14, full-box-2 1, full-box-1 0, ambiguous 0. Each timing-mask bit
  count exactly matched its source histogram fill count.

Both jobs processed 653966 input entries and reported 44 exclusive events,
ellipse 5, MCD 36, and peak fraction 0.050. Macro runtime was 29.682544 s for
default-off and 27.746202 s for opt-in; this single warm-cache ordering is not
a performance claim. Each isolated output directory used 1.9 MB.

Key SHA-256 values:

| Input/output | SHA-256 |
|---|---|
| run-5237 replay input | `fcc1544c1e19c79d7eccad393380b68ccdbcb012df89ac13db05c3972c51979d` |
| frozen selection report | `c301526074bf5c45bf9101efb6fa971a3b6a6581faf798c0bfb62d43b3da55b1` |
| main run configuration | `f3ee5f17956e0291a5c6f10f0b96b9cfaee222a8155012d9fa7be6814bee61c8` |
| acceptance cuts | `1f2ebf451c8ffd22f3d3f06140a839cf83cb729c9d88d8da0f9ce317136bca05` |
| dead-block configuration | `3c3e35e661428262a44f585b4b4e9319a3ee455535e1dd9d7f52f0af3af086e3` |
| default diagnostic ROOT | `ae5162f4eb453ec7ecc203ec8e5671d5867aff79b1f06dc2d7eac6efe9d410d7` |
| opt-in diagnostic ROOT | `0fea41e0eb31d229ce71abf8f0d34c15933eec000fa630e20f9dc4040d45d81e` |
| identical summary CSV | `648c73fb2edb59328fdcc6d5b6b5babba1931977ce8d2f628ba12422166faf71` |

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

- Validate a complete run ledger for one setting, including intended zero-event
  and missing/partial statuses.
- Exercise combined default-off/opt-in output on that complete setting.
- Run the later pseudoexperiment closure/coverage suite after its likelihood is
  separately approved and implemented.

No production-readiness or publication claim is made until those checks and the
later background/smearing/extraction proposals pass.
