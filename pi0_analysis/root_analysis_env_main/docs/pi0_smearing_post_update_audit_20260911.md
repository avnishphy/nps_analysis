# Post-update smearing audit

Date: 2026-09-11. Target: `root_analysis_env_main/src/simulation_smearing/`. Read-only audit of the four requested updates; no analysis code changed. Findings describe the current source, including some pre-existing behavior exposed by the expanded validation.

## Outcome by requirement

| Requirement | Observed implementation | Assessment |
|---|---|---|
| Always generate fresh nominal simulation | Unconditional `nominal-only` generation into unique staging; fitter receives its path and SHA-256 | Implemented |
| Preserve nominal input through final production | `smeared-only` raw-input replay; nominal SHA checked afterward; no nominal rewrite | Implemented, with input-identity limitations below |
| Propagate failures | Compiled producer returns integer status; staged outputs/provenance validated before publication | Major previous gap fixed; branch-binding/read errors and publication failure still need attention |
| Add interpolated-response summaries | Additional data/section/map overlays and metrics; original section behavior retained; same events, draws and normalization | Implemented; bounded map tests pass |
| Extend fast evaluation and isolate caches | Batched final grids and broader fast objective use; evaluator-owned caches, per-candidate accumulators, rebuilt contexts across sweeps | Implemented substantially; numerical boundary discrepancy and unnecessary duplicate evaluations remain |

## 1. High priority: incompatible input branch types still appear successful

**Location:** [producer:1846](../src/simulation_smearing/simc_pi0_analysis.C#L1846), also bindings at 1786-1788; active-loop entry read at [2161](../src/simulation_smearing/simc_pi0_analysis.C#L2161). Output validation checks branch presence, entry count and provenance at [2652](../src/simulation_smearing/simc_pi0_analysis.C#L2652).

The producer checks that required branch names exist, but ignores `SetBranchAddress` return codes. It also ignores `GetEntry` results. ROOT can report a binding error without the producer returning failure. A failed binding can leave a default value; a failed entry read can leave values from an earlier entry.

**Independently reproduced:** using a freshly compiled current producer, a synthetic fixture was changed only from `Float_t Weight` to `Double_t Weight`, retaining valid selection branches. ROOT printed a type mismatch for each input channel, yet nominal production returned **0**, produced **120 entries**, and `--validate-output` returned **0**. This is a demonstrated failure-propagation defect, not a claim that existing production inputs have the wrong type. Artifacts and log: `/tmp/pi0_weight_type_audit/`.

The output was numerically wrong: input weights `2.0-2.6` with normalization `1/40` should produce `full_weight=0.050-0.065`, but every output weight was **0.025**. Failed binding left the producer's default `Weight=1` in place (declaration at 1354; unchecked binding at 1850).

**Remedy:** centralize checked required-branch binding; fail on incompatible binding status with sample and branch name. Check every expected entry read and fail with sample/index if it fails. Add malformed-type and read-failure tests. This preserves valid-event physics; it must not be implemented by adding new physics cuts or silently coercing incompatible inputs.

## 2. Fast histogram boundary assignment is not ROOT-equivalent

**Location:** [fitter:1400](../src/simulation_smearing/nps_sim_smearing_new.C#L1400).

`FastHistogram1D::binIndex` multiplies by a precomputed reciprocal bin width. ROOT's axis calculation uses a different floating-point operation order. Mathematically equivalent expressions can round to opposite sides of an integer at a bin edge.

**Independently reproduced with actual configured axes:** tests at each boundary and its immediately adjacent representable values found **18 bin-assignment differences**: 5 for Mgg, 7 for missing mass, and 6 for Mpgg2. Missing mass is currently inactive in the objective; the other two are active. Examples, with zero-based bin indices:

```text
Mgg:   x = 0.1172                 ROOT bin 8; fast bin 9
Mpgg2: x = 5.7999999999999998     ROOT bin 8; fast bin 9
```

The existing equivalence harness passes because its boundary-specific histogram uses simple integer-width bins; the tested event sample does not expose this discrepancy. The three-point runtime objective check also cannot guarantee agreement at all subsequent points. This behavior predates the latest extension of the fast path; wider use makes resolving it important.

**Remedy:** reproduce the ROOT axis arithmetic/boundary policy exactly, or use its axis lookup where appropriate; test every actual histogram axis at edges and adjacent representable values. Compare bin contents and sumw2, not only total chi-square. Do not loosen tolerances to conceal a different bin assignment. No production bias or altered fit minimum was measured here.

Reproducer: `/tmp/pi0_latest_fast_regression.C`.

## 3. Publication can leave a mixed result set on failure

**Location:** [pipeline:978](../src/simulation_smearing/run_smearing_pipeline.sh#L978), publication sequence at 992-1017, cleanup at [721](../src/simulation_smearing/run_smearing_pipeline.sh#L721).

Each individual file replacement is atomic, but the collection is not. If a later copy/rename fails after nominal or smeared output has been replaced, some latest files belong to the new run and others to the previous run. The final archive is committed only after these replacements. The EXIT cleanup then deletes both staging directories, removing the completed run's recovery copies. This differs from a failure during generation/fitting, which now preserves previous latest outputs.

**Remedy:** commit a validated immutable run directory first, then atomically update one manifest/pointer selecting that run. If legacy paths must also be maintained, retain the committed run when compatibility-copy publication fails and provide a recoverable transaction/rollback policy. Use uniquely created publication temporary files; `${dst}.tmp.${RUN_TAG}` can collide between concurrent runs sharing a tag.

This is a publication/recovery issue, not a detector-response change. The implementation documentation correctly describes its present guarantee as atomic *per file*.

## 4. Raw/configuration identity does not prove unchanged contents

**Location:** [pipeline:618](../src/simulation_smearing/run_smearing_pipeline.sh#L618).

The input identity combines resolved path, device, inode, size and modification time rounded to integer seconds. Different same-size contents written within one second can retain exactly this identity. The reviewer reproduced this using small temporary configuration files.

The fresh nominal ROOT is protected by SHA-256, but final production rereads raw ROOT/configuration files. Consequently the nominal hash and equal selected counts alone do not prove identical final input contents/events. No input mutation during an actual production run was demonstrated.

**Remedy:** copy/hash small configuration and normalization files into staging and use those snapshots in both passes. Use a content identity or genuinely immutable production manifest for large ROOT inputs. Record a complete build identity as well: the current source hash includes the producer and shared response header but omits other included reconstruction helpers. Fresh compilation avoids stale executable reuse; this remaining issue is provenance completeness.

## 5. Cache lifetime looks sound; invalid budget input bypasses the limit

No cross-event cache contamination or growing shared cache was found in the inspected call paths. Event vectors outlive their evaluators; per-section contexts finish before coupled-sweep updates; shared evaluators used for independent starts are read-only with local accumulators. Construction/destruction and forced-uncached equivalence tests pass. This is not a full leak-sanitizer or production memory-profile result.

One concrete input-validation problem remains at [fitter:485](../src/simulation_smearing/nps_sim_smearing_new.C#L485): `strtoull` accepts `-1` as the largest unsigned value. Setting `NPS_FAST_PULL_CACHE_BUDGET_BYTES=-1` produced **18446744073709551615**, effectively removing the intended limit on this platform.

**Remedy:** reject negative/signed-invalid values before unsigned parsing and test malformed, overflow, zero and valid budgets. Keep the documented default/failure policy explicit. The default 1 GiB setting is unaffected, and it limits pull arrays, not all process memory.

## 6. Remaining performance redundancy

The final 856-point visualization grids now use fast batched evaluation with a legacy validation/fallback path. This addresses the main earlier performance recommendation.

However, global-objective diagnostics ([fitter:5477](../src/simulation_smearing/nps_sim_smearing_new.C#L5477)) and final section breakdowns ([6903](../src/simulation_smearing/nps_sim_smearing_new.C#L6903)) unconditionally calculate a complete legacy result, construct a fast evaluator/cache, then calculate the fast result. The duplicate work occurs even with `VALIDATE_FAST_OBJECTIVE=false`. In these one-shot contexts the cache construction is not amortized over a large scan.

**Remedy:** reuse a compatible existing evaluator/result where its lifetime permits; gate reference evaluation on validation policy and perform it at selected context checks. Retain explicit fallback. For a truly one-shot result, benchmark whether constructing a pull cache is worthwhile. Existing ROOT histogram rendering can remain separate from objective evaluation; replacing it merely to call everything 'fast' is unnecessary.

## Validation performed for this audit

- Rebuilt the existing equivalence harness against the current source: `equivalence_failures=0`.
- Added the actual-axis boundary and negative-budget tests described above; they expose cases absent from that passing harness.
- Independently tested producer-style map lookup against materialized ROOT maps at 792 internal-edge/adjacent-value probes over the default geometry. No coefficient mismatch above `1e-12`; maximum difference `2.22e-16`. Constant/varying-map tests in the rebuilt equivalence harness also pass. Section-versus-interpolated fitting/production behavior was not changed.
- Fresh producer compilation plus the malformed-Weight fixture demonstrated false success. The producer also passes a separate syntax check; the current fitter compiles as part of the rebuilt harness.
- `bash -n` and targeted `git diff --check` pass.
- No full production run, full pipeline failure-recovery run, new production timing measurement, or sanitizer campaign was performed. No existing physics result was shown to be biased.

Reproduce bounded checks after sourcing the Hall C environment:

```bash
# Working directory: root_analysis_env_main
bash -n src/simulation_smearing/run_smearing_pipeline.sh
git diff --check -- src/simulation_smearing
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ /tmp/pi0_smearing_equivalence_20260911.C `root-config --cflags --libs` -lMathMore -O2 -std=c++17 -fopenmp -o /tmp/pi0_latest_audit_equivalence && /tmp/pi0_latest_audit_equivalence'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ /tmp/pi0_latest_fast_regression.C `root-config --cflags --libs` -lMathMore -O2 -std=c++17 -fopenmp -o /tmp/pi0_latest_fast_regression && /tmp/pi0_latest_fast_regression'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ /tmp/pi0_latest_map_edge_audit.C `root-config --cflags --libs` -lMathMore -O2 -std=c++17 -fopenmp -o /tmp/pi0_latest_map_edge_audit && /tmp/pi0_latest_map_edge_audit'
```

Supporting notes: `/tmp/pi0_latest_pipeline_audit.md`, `/tmp/pi0_latest_fast_audit.md`, `/tmp/pi0_latest_update_audit_working.md`. The producer failure reproduction commands/logs are recorded with the pipeline notes and `/tmp/pi0_weight_type_audit/`.
