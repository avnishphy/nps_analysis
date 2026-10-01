# Pi0 smearing pipeline implementation and validation (2026-09-11)

## Implemented data flow

Each pipeline invocation now creates a unique staging directory and follows this flow:

1. Resolve and validate raw exclusive, SIDIS and delta Geant inputs, normalization inputs, acceptance/kinematics/dead-block configuration, and explicit overrides. Record their path/stat identities and resolved settings in a hashed input snapshot.
2. Compile the producer and fitter, then run the producer in `nominal-only` mode into staging. This step is unconditional: neither an existing nominal file nor an old fingerprint can bypass it.
3. Validate the staged nominal ROOT file (`simulation` tree, required branches, positive entries, current run/input provenance), calculate its SHA-256, and give that exact path to the fitter.
4. After fitting and validating the current-run fit ROOT, section CSV, interpolated ROOT maps, comparison CSV and provenance, rerun the producer in `smeared-only` mode from the same raw files and resolved reconstruction settings. Validate the smeared output, equal selected entry count, map/input/run provenance, unchanged raw/config snapshot, and unchanged nominal SHA-256.
5. Copy only current-run staged artifacts into the run archive and publish each final file through a same-directory temporary file plus atomic rename. A failed producer, fitter or validation stage cannot publish or archive a stale prior file.

The nominal output tree is not used as the final producer's event source. It contains only the selected pair and stores cluster quantities as `Float_t`; the producer reconstructs/selects/pairs raw `vector<double>` clusters and may inspect other candidates. Literal tree reuse would therefore change selection information and precision. The smeared-only pass instead replays the immutable raw inputs/configuration and verifies the selected entry count, while the preserved nominal SHA-256 ties fitting and final provenance to the exact nominal artifact.

The producer now returns explicit nonzero status codes, supports `both` (standalone-compatible default), `nominal-only`, and `smeared-only`, and validates its output through `--validate-output`. ROOT write failures, missing trees/branches, zero selected events, missing maps in smeared-only mode, and provenance mismatches are fatal. The wrapper retains fitter status 20 as deliberate input-preview cancellation, atomically retains that preview at its requested path, and treats other nonzero statuses as failures.

## Response comparison

Existing data/unsmeared/section-response summaries are unchanged. A new labelled PDF page, ROOT canvas/histograms, and `smearing_response_lookup_comparison.csv` add Data, section-coefficient simulation, and producer-style interpolated-map simulation for Mgg, Mmiss and Mpgg2. The two simulation curves use copied versions of the same ordered event buffer, identical weights/cuts/Nsmear/normalization, and separately reset `TRandom3(42)` streams, so every event/replica/photon receives the same draws.

The map curve reproduces producer rules: section values are the fallback, inclusive physical map boundaries are accepted, ROOT-style fixed 100x100 bin selection is clamped to an in-range bin, and the persisted bin-center interpolation supplies each of the five response fields. The map is built from the final accepted fitted coefficients. Metrics include integrals, L1/max bin differences, maximum bin pull difference, and data chi-square for both response lookups. No fit or production policy changes automatically from these diagnostics.

## Fast evaluator ownership and invalidation

`FastObjectiveEvaluator` is an RAII object owned by one immutable event/configuration evaluation context. It references an event vector that must outlive it and owns prepared invariants, copied data histogram bins/sumw2 and optional Gaussian pulls. Cached values depend on event contents/order/weights, geometry, histogram definitions, Nsmear, response shape and deterministic seed. Photon external coefficients are part of the immutable event records. Coupled sweeps construct a new evaluator after coefficient updates; section workers own isolated evaluators.

Only parameter-invariant values and random pulls are cached. Parameter-dependent means, resolutions, clamps, positions and kinematics are recomputed per candidate. Candidate accumulators are newly constructed/reset per scalar evaluation or batch chunk. Accumulation remains event, replica, photon order per candidate and uses `double`. The total pull-cache budget defaults to 1 GiB and can be lowered with `NPS_FAST_PULL_CACHE_BUDGET_BYTES`; allocation/budget failure selects the correct uncached scalar path. No global mutable evaluator/cache was added.

Fast evaluation now covers section optimization (pre-existing), independent multistart batches, final visualization grids, optional global electron-scale coarse/refinement grids, optimizer scalar requests, coupled-sweep global diagnostics and final objective breakdowns. Independent grids are batched; sequential optimizer scheduling is unchanged. The legacy evaluator remains the runtime reference and fallback. Non-objective final ROOT histogram rendering still uses the existing histogram-fill path.

## Validation performed

All commands ran from the repository root after:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; <command>'
```

- ROOT 6.30.04: both C++ sources passed `g++ -fsyntax-only`; both full executables compiled. `bash -n` and `git diff --check` passed.
- `/tmp/pi0_smearing_equivalence_20260911.C`: zero failures comparing legacy, fast scalar and fast batch objectives plus every in-range histogram bin content and sumw2. Cases covered histogram min/max/underflow/overflow rules, zero/nonzero position response, external photon coefficients, different event buffers and geometry, repeated evaluator construction/destruction, batch resets, and forced zero-byte cache-budget fallback.
- The same harness showed exact section/map histogram agreement for constant response maps. A spatially varying 2x2 calibration exposed a nonzero section-versus-map histogram difference, while the producer-style lookup matched a materialized ROOT `TH2D::FindBin` map at interior and inclusive boundary points.
- Synthetic producer test: nominal-only created and validated 120 entries; smeared-only replay created and validated 120 entries. The nominal SHA-256 remained `8c2a21e2f10fd66194b5e0c874c10b5cbeefb37960ef9c2f5a55d77deb9799aa` after final production.
- Failure injection: missing raw input exited 1 before compilation; invalid ROOT producer input exited 6; a valid fresh producer followed by an invalid data tree exited 3 at the fitter. In all cases the pre-existing stale nominal checksum (`c398c01827d768ba836ce68ad839bdf4db09f182e6d8fb660f61fe75c091f9c0`) remained unchanged and no current-run archive/output was published. Missing and malformed output validation both exited 10. Input-preview rejection returned the distinct fitter status 20 and wrote the preview without a fit output.
- The wrapper's ROOT validation expression was exercised against readable fitter/interpolated outputs with all five required response maps and provenance, returning 0.
- Bounded deterministic 1x1, 150-data/150-simulation, Nsmear=2 fits were run with the pre-edit validation binary and current binary. Both returned the same minimum and section row: chi2 `26.3645`, `(mu_a,mu_b,mu_c)=(0.1482,0.972038,-0.197444)`, `sigma=0.0142685`, `sigma_pos=0.0593717`, two sweeps ending by repeated rejected candidate, and identical `hesse_ok=0` (no covariance was available in either run). The pre-edit binary is a `/tmp` artifact created earlier on 2026-09-11; its source-to-binary identity was not embedded, so this before/after comparison is supporting evidence rather than a reproducible proof of the precise pre-edit tree.
- Evaluator-only benchmark, one thread, 200 events, 80 replicas, 16 points, seven repeats: legacy median `48.605 ms`, fast scalar `6.023 ms` (8.07x), fast batch `5.226 ms` (9.30x). Maximum legacy/fast objective difference was `8.53e-14`; batch/scalar difference was zero. This does not measure producer, I/O, plotting or full-pipeline speed.

## Remaining validation limits

No full production campaign was run. The active Gaussian/simple-resolution configuration was tested. Compile-time-disabled Landau response, three-term energy resolution, energy-dependent position width, linear Mgg correction and the disabled global electron-scale stage were not toggled; their unchanged legacy fallback paths remain available. The synthetic producer's 2D exclusivity training was intentionally not used as physics validation. These checks demonstrate the stated bounded equivalences and failure behavior, but do not prove the absence of all defects.
