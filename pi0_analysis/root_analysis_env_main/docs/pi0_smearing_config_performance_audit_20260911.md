# Updated smearing configuration, pipeline and performance audit

Date: 2026-09-11. Reviewed repository: `root_analysis_env_main`, specifically `src/simulation_smearing/`. The sibling `root_analysis_env` still contains an older workflow and was not the audit target. No analysis source or physics logic changed.

## Assessment

The shared configuration refactor is connected correctly for the active response: both macros import `smearing_response_config.h`, including Gaussian energy smearing, the energy-mean convention/floor, position smearing and optional response controls. The producer now uses `ENABLE_POSITION_SMEARING=true`, matching the fitter. Removed fitter model classes/reference cuts/scan constants were unused by the surviving response path. Representative Gaussian regression checks pass.

This confirms the narrow configuration update, not complete equivalence of fitter diagnostics and final production. Two pipeline issues and an existing coefficient-lookup difference remain. Computational savings are available without changing the fitted physics model; the strongest candidate is accelerating the final diagnostic scans using the machinery already used for optimization.

The git diff includes older uncommitted producer changes as well as the latest refactor. Findings below describe current behavior; they are not all attributed to the latest update.

## 1. Findings requiring attention

### A. Existing unsmeared input is reused without checking its provenance

**Evidence:** [pipeline:662](../src/simulation_smearing/run_smearing_pipeline.sh#L662) generates nominal simulation only if its file is missing. After the fit, [pipeline:723](../src/simulation_smearing/run_smearing_pipeline.sh#L723) invokes the producer again. The producer opens both nominal and smeared outputs with `RECREATE` ([producer:1319](../src/simulation_smearing/simc_pi0_analysis.C#L1319)).

**Consequence:** if cuts, geometry, normalization, input production or producer behavior changed, the fitter can use old nominal events while the final nominal file is regenerated with new settings. The final archive then need not contain the nominal input actually fitted. Sharing response constants does not solve this input-consistency problem. A changed smearing-only flag need not invalidate nominal events, but relevant upstream changes must.

**Solution:** fingerprint nominal-producing source/dependencies, input identities, normalization and resolved reconstruction settings. Reuse nominal output only on a match; otherwise regenerate before fitting. Preserve the exact fit-input identity in the run manifest. Avoid rewriting that input after fitting. This is orchestration/provenance work, not a change of selection formulas.

### B. Producer failure can be reported as pipeline success

**Evidence:** configuration exceptions return from the producer's void macro ([producer:1061](../src/simulation_smearing/simc_pi0_analysis.C#L1061)); missing required smearing files also return ([producer:1085](../src/simulation_smearing/simc_pi0_analysis.C#L1085)). The shell invokes ROOT without an explicit producer-status contract ([pipeline:659](../src/simulation_smearing/run_smearing_pipeline.sh#L659)). After the final call, it does not require fresh valid output; archival merely warns if a file is absent ([pipeline:728](../src/simulation_smearing/run_smearing_pipeline.sh#L728)).

**Consequence:** a handled macro failure can leave ROOT returning successfully, allowing an older output to be archived and completion to be printed. Shell `set -e` cannot detect an error that was converted to a normal return.

**Solution:** explicit nonzero failure signaling/completion status; write temporary outputs, validate successful close/tree/schema/current-run identity, then publish. Preflight resolved default production paths before the expensive fit, as well as explicit overrides. File existence alone is insufficient.

### C. Matching settings still leaves different coefficient application

**Evidence:** the fitter's final all-section summary selects discrete fitted section coefficients and interpolates only as fallback ([fitter:7151](../src/simulation_smearing/nps_sim_smearing_new.C#L7151)). The producer preferentially samples interpolated maps ([producer:1470](../src/simulation_smearing/simc_pi0_analysis.C#L1470)).

**Consequence:** for a spatially varying map, the final summary and producer can use different response coefficients at the same photon position. This is pre-existing, not a regression caused by the shared header. Matching summary plots therefore does not prove agreement with the exact production response.

**Solution within the requested boundary:** add a diagnostic comparison using the producer's actual map lookup and quantify the difference. Do not silently change section assignment or interpolation as a cleanup. Any later decision to unify those algorithms should be explicit.

### D. Shared defaults do not validate older fitted artifacts

**Evidence:** interpolated-map loading reads response histograms but does not verify the complete generating response configuration ([producer:1210](../src/simulation_smearing/simc_pi0_analysis.C#L1210)). The shared header is included in the fitter's cache source fingerprint ([fitter:4797](../src/simulation_smearing/nps_sim_smearing_new.C#L4797)), which is useful but does not establish compatibility of arbitrary producer input maps.

**Solution:** record effective response settings with fit artifacts and check them when loading. Distinguish an explicit legacy compatibility policy from a verified match. This matters particularly when applying old maps after changing a response switch. The current header alone does not make those old files compatible.

## 2. Faster computation without changing the physics model

### Highest priority: final diagnostic scans

[fitter:2862](../src/simulation_smearing/nps_sim_smearing_new.C#L2862) generates the final visualization scans through `eval_chi2_selected`, the legacy ROOT-histogram/RNG evaluation path. Current settings request

```text
28*28 + 24 + 24 + 24 = 856 objective evaluations per successful section
```

The last 24 evaluations are the enabled position scan. Each evaluation traverses the section's events and smearing replicas. These scans are regenerated from final fit values, separately from optimization. At 80 replicas the cost scales as `856 * N_section * 80` pair replicas per section. This is a static operation count, not a production timing measurement.

**Recommendation:** reuse one fast evaluator/cache per section and batch the same parameter points, preserving scan grid, event order, weights, random pulls, histogram ranges and normalization. The fitter already has a fast objective and batched Sobol evaluation; this avoids introducing a new approximation. Validate old/new per-observable histograms and chi-square values, especially near bin edges: the existing fast path uses some different arithmetic expressions, so mathematical equivalence alone is not proof of bitwise identity. Keep fitted coefficients unchanged.

**Measured supporting check:** a small harness invoking the actual current evaluator APIs used 200 synthetic events, 80 replicas and 16 parameter points, with seven timing repeats and one thread (ROOT 6.30.04, `-O2`). Both currently active observables and nonzero position response were exercised.

| Evaluator | Median time for 16 points |
|---|---:|
| Legacy | 49.428 ms |
| Fast scalar | 5.901 ms |
| Fast batch | 5.074 ms |

Fast-evaluator construction took another 2.723 ms. Maximum absolute fast/legacy objective difference was `8.53e-14`; batch/scalar difference was zero; all 16 points met existing validation tolerances. This is a synthetic evaluator benchmark, not a production speed estimate or full histogram/fit-equivalence test. The pull cache fit in memory; large-section fallback and inactive modes were not tested.

Optional diagnostic thinning or deferral can also save time, but changes the diagnostic product. Do not confuse it with an exact replacement of the current scans, or reduce the optimizer's replicas/search as a purported logic-preserving optimization.

### Existing acceleration and the cache limit

The fitter already caches event invariants and Gaussian pulls as doubles, uses fixed histogram arrays, batches Sobol candidates and parallelizes suitable fit work. Recommending those features as though absent would be redundant.

The pull cache has a process-wide 1 GiB budget divided among active fits ([fitter:473](../src/simulation_smearing/nps_sim_smearing_new.C#L473), [fitter:1873](../src/simulation_smearing/nps_sim_smearing_new.C#L1873)). Six doubles per event/replica at 80 replicas require **3,840 bytes/event**: approximately 279,620 events for the whole pull budget before its per-fit division. Oversized fits fall back to RNG generation; the batched fast route also depends on cached pulls.

**Recommendation:** first record cache-use/fallback counts, event counts, time per evaluation and peak memory on a representative run. Consider reusing compatible immutable pull arrays or processing fixed event blocks across candidate batches. Preserve each evaluator's existing seed/index ordering and histogram accumulation order. Do not use float pulls, new random sequences or uncontrolled reductions merely to reduce memory. More threads can reduce each fit's cache allowance, so more threads are not automatically faster.

Another bounded optimization is caching prepared response terms for photons outside the section being fitted ([fitter:1530](../src/simulation_smearing/nps_sim_smearing_new.C#L1530)). Their external coefficients are fixed within that section fit. Rebuild this cache when a coupled sweep changes those coefficients; retain any later pair-dependent corrections. Actual benefit needs profiling. Existing sweep `runtime_sec` is repeated in section-summary rows; summing those rows would overcount wall time.

### Pipeline and producer savings

| Opportunity | Current evidence | Safe implementation condition |
|---|---|---|
| Reuse compiled fitter | Pipeline:670-678 builds into a temporary directory and deletes the executable each invocation | Cache against compiler/ROOT versions, flags, source and complete header dependencies; include the shared config header |
| Avoid rewriting nominal output during final smearing application | Producer:1319 and 2438 recreate nominal tree and repeat its 2D mass-cut processing | Add a smeared-only output mode when the exact nominal input is already available and verified; retain current event/RNG ordering |
| Reuse deterministic photon preparation | Producer:2324-2331 computes energy mean for the energy draw and again for the position reference | Prepare mean/width once per photon, then preserve the exact sequence of energy/x/y draws and clamps |

Preparing nominal input and applying the fitted response are distinct necessary stages. The redundant part is repeated nominal output/diagnostics, not the existence of those two stages. Any larger intermediate-event cache must preserve the information needed to reproduce the current reconstruction; changing selection or pairing is outside this performance refactor.

## 3. Smaller cleanup findings

- The new shared configuration header is currently untracked. Include it in any eventual commit/package containing the two changed macros; it is now a required dependency.
- The shell's interpolated-output name derivation treats only `.root` as an extension, whereas the fitter inserts the suffix before any extension ([pipeline:91](../src/simulation_smearing/run_smearing_pipeline.sh#L91), fitter:735). Defaults work; custom `out.dat` produces mismatched expectations. Require `.root` or use a common naming rule.
- Dated archives and latest-output copies serve provenance. They are not duplicate physics computation. Preserve them; filesystem reflinks can reduce copying where supported. Do not hardlink files later opened with `RECREATE`.
- The input-preview approval step is intentional. Binary/input caching can reduce repeated startup work after cancellation without removing it.

## 4. Validation performed and limits

- `bash -n` and targeted `git diff --check`: pass.
- Both macros pass C++17 syntax checks with ROOT 6.30.04 loaded through the Hall C setup.
- Rebuilt and ran the existing small fitter and producer response-regression harnesses against the current source. Both pass: representative Gaussian mean/width and draw sequence agree with the pre-refactor formulas; the producer applies the intentionally enabled position response.
- These tests cover representative active-Gaussian calculations, not all optional response modes, map compatibility, full fit convergence or production closure. Enabling position draws also changes later producer RNG consumption; whole-file equality against the old position-disabled run is not an appropriate acceptance criterion.
- No full production fit or pipeline execution was performed. No production wall-time speedup is claimed.
- The synthetic evaluator benchmark above was executed against the current source; its harness is `/tmp/pi0_smear_perf_bench_20260911.C`.

Commands run from `root_analysis_env_main` (regression files are temporary artifacts):

```bash
bash -n src/simulation_smearing/run_smearing_pipeline.sh
git diff --check -- src/simulation_smearing
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ -fsyntax-only -std=c++17 -fopenmp `root-config --cflags` -Isrc src/simulation_smearing/nps_sim_smearing_new.C; g++ -fsyntax-only -std=c++17 `root-config --cflags` -Isrc src/simulation_smearing/simc_pi0_analysis.C'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ /tmp/nps_smearing_refactor_validation/fitter_response_regression.cpp `root-config --cflags --libs` -lMathMore -O2 -std=c++17 -fopenmp -o /tmp/pi0_reaudit_fitter_regression && /tmp/pi0_reaudit_fitter_regression; g++ /tmp/nps_smearing_refactor_validation/producer_response_regression.cpp `root-config --cflags --libs` -O2 -std=c++17 -o /tmp/pi0_reaudit_producer_regression && /tmp/pi0_reaudit_producer_regression'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ /tmp/pi0_smear_perf_bench_20260911.C `root-config --cflags --libs` -lMathMore -O2 -std=c++17 -fopenmp -o /tmp/pi0_smear_perf_bench_20260911'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; /tmp/pi0_smear_perf_bench_20260911'
```

One sandbox invocation failed because the Hall C module setup could not resolve user ID 11606. Repeating the regression build outside the sandbox successfully loaded ROOT and passed both tests.

Working evidence: `/tmp/pi0_smearing_reaudit_20260911.md`, `/tmp/pi0_pipeline_reaudit_20260911.md`, `/tmp/pi0_perf_reaudit_20260911.md`.
