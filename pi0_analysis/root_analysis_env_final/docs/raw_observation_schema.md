# ALG-001 raw-observation schema

Status: opt-in engineering foundation; not a publication estimator.

## Enablement and default

The legacy default is off. Enable the passive exporter through the unified
launcher:

```bash
cd /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash src/analysis/run_parallel_nps_analysis_main.sh --raw-observation-export --output-base /scratch/singhav/alg001_shadow_SETTING --efficiency-csv /absolute/read_only/existing_efficiency_SETTING.csv [normal selection options]'
```

The launcher sets `NPS_RAW_OBSERVATION_EXPORT=yes` for each per-run ROOT job and
passes `--raw-observation-export` to the combiner. Without the flag, the macro
does not create either raw tree, the combiner continues excluding `event_id`,
and no ALG-001 ledger is written.

The command above is a template, not an executed production command. Replace
`SETTING`, the read-only existing correction CSV, and the bracketed selection
only after choosing a validation run. The launcher rejects opt-in execution at
the canonical `FINAL/output` path; a separate explicit output base is required
so no legacy artifact is overwritten. `--efficiency-csv` only selects an
existing input for the unchanged combiner calculations.

## Per-run `raw_observation` tree

One row corresponds to one entry in the existing selected-event `physics` tree,
before timing/combinatorial subtraction or `pi0_weight` assignment.

| Field | Meaning |
|---|---|
| `run_number` | Run key. |
| `event_id` | Zero-based entry in the full input `TChain`; 64-bit. |
| `source_tree_number` | Zero-based chain segment index. |
| `source_entry` | Entry index in that segment's underlying tree. |
| `event_number` | Existing `g.evnum` value when available. |
| `t1_ns`, `t2_ns`, `pair_dt_ns` | Selected pair times and `t1-t2` after the existing mode offset. |
| `timing_region_mask` | Bit mask of all current timing-region predicates. |
| `timing_category` | Exclusive category derived from the mask. |
| `acquisition_mode` | `1=hcana`, `2=waveform`. |
| `shifted_sidebands` | Existing resolved sideband mode, `0/1`. |
| `pair_time_diff_max_ns` | Existing resolved pair-time cut; recorded, not reapplied. |
| `nclust_selected` | Existing selected-cluster multiplicity. |
| `passes_mmiss_exclusive_cut` | Result of the existing missing-mass cut, evaluated in the first pass without requiring positive purity weight. |
| `helicity` | Existing resolved helicity. |
| `mpi0_all`, `mmiss_all`, `mmiss_all_corr`, `Q2`, `W`, `t`, `tmin`, `phi`, `xB` | Existing reconstructed observables copied before purity weighting. |

No `pi0_weight`, subtraction result, fitted background parameter, scale, or
efficiency-derived event weight is present in this tree.

## Timing codes

| Category | Code | Mask bit |
|---|---:|---:|
| Outside the named regions | 0 | none |
| Prompt | 1 | `1` |
| Diagonal | 2 | `2` |
| Horizontal | 3 | `4` |
| Vertical | 4 | `8` |
| Full box 1 | 5 | `16` |
| Full box 2 | 6 | `32` |
| Ambiguous/overlapping predicates | 7 | multiple bits |

Open interval comparisons exactly match the existing histogram-fill code.
`outside` makes the category exhaustive over all selected pairs. An ambiguous
category is retained in the per-run file for diagnosis but blocks combined
ALG-001 export.

## Segment and run ledgers

Each per-run file contains `raw_observation_segments`, mapping
`(run_number, source_tree_number)` to the original tree name and input path.

The combiner writes beside `combined_branches_<target>.root`:

- `combined_branches_<target>_run_ledger.csv`;
- `combined_branches_<target>_segment_ledger.csv`.

The run ledger has a row for every intended config run after setting, target,
type, and explicit run filters, including zero-candidate and missing runs. It
records existing charge/correction central fields without changing or
recalculating efficiency definitions. Statuses outside `ready`,
`zero_candidate`, and `excluded_known_bad` make the opt-in combine fail with
exit status 3 after writing the diagnostic ledgers.

The combined ROOT file retains the legacy `physics` tree, adding `event_id` only
in opt-in mode, and adds a separate unweighted `raw_observation` tree.

## Reproduce the implemented checks

Actually executed from the FINAL workspace:

```bash
bash -n src/analysis/run_parallel_nps_analysis_main.sh
python3 -m py_compile src/analysis/combine_analysis_branches.py tests/test_raw_observation_bundle.py
g++ -std=c++17 -Wall -Wextra -pedantic tests/test_nps_raw_observation.cpp -o /tmp/test_nps_raw_observation_20261001
/tmp/test_nps_raw_observation_20261001
python3 tests/test_raw_observation_bundle.py
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; root -l -b -q -e '\''gSystem->SetBuildDir("build/alg001_phase1", kTRUE); gROOT->ProcessLine(".L src/analysis/nps_analysis_main.C+");'\'''
```

The first command checks shell syntax. `py_compile` checks Python parsing.
The C++ executable exercises every timing category and strict boundaries. The
Python test creates temporary synthetic ROOT files, checks ready/zero/missing
ledger states, proves `event_id` remains excluded by default but is preserved
opt-in, and round-trips the combined raw tree. ACLiC compiles the actual macro
against ROOT 6.30.04. Full results and the known baseline failure are in
`validation/ALG001_phase1_validation_20261001.md`.
