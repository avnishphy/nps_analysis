# FINAL path-isolation validation — 2026-10-01

## Claim and held-fixed behavior

The 18 targets in `manifests/path_relocation_map.tsv` change only the literal
workspace component `root_analysis_env_main` to `root_analysis_env_final`.
Relative input/configuration/output suffixes and all analysis behavior are held
fixed. This is an engineering relocation, not an algorithm or physics change.

Recoverable pre-edit copies are under ignored
`recovery/path_isolation_pre_edit/`. The first 17 targets were verified by
substituting the old prefix in each recovery copy and comparing the result
byte-for-byte with the working file. `tests/run_xsec_tests.sh`, discovered by a
broader follow-up scan, was backed up and verified with the same comparison.

## Checks executed

All commands were run from FINAL. ROOT commands used:

```bash
csh -f -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; root -l -b -q -e "gSystem->SetBuildDir(\"build/path_isolation\", kTRUE); gSystem->CompileMacro(\"<target>\", \"k\");"'
```

The four `<target>` values were:

- `src/analysis/nps_analysis.C`
- `src/analysis/nps_analysis_wfpi0.C`
- `src/analysis/nps_analysis_main.C`
- `src/simulation_smearing/simc_pi0_analysis.C`

Additional checks:

```bash
bash -n <each affected shell file>
PYTHONPYCACHEPREFIX=scratch/pycache_path_isolation python3 -m py_compile <each affected Python file>
python3 -m json.tool <each affected JSON file>
g++ -std=c++17 -fsyntax-only $(root-config --cflags) <affected efficiency source>
cmp src_entries_before_root_compile.tsv src_entries_after_root_compile.tsv
sha256sum --check root_targets_before.SHA256SUMS
git diff --check -- <path-isolation targets>
```

## Results

- Exact-substitution comparisons: 18/18 passed.
- Shell syntax: passed.
- Python byte-compilation: passed; cache confined to `scratch/`.
- JSON parsing: passed.
- Efficiency C++ syntax checks: passed.
- ROOT/ACLiC compilation: 4/4 shared libraries produced under
  `build/path_isolation/`; no compile errors.
- Source immutability: before/after entry manifests compare equal and all four
  source checksums verify.
- Compiler diagnostics retained in runtime logs: ROOT reported a possible C++
  standard-library ABI mismatch for all four macros; `nps_analysis_main.C`
  also reported inherited signed/unsigned comparison warnings.
- No production launcher or analysis event loop ran.

Detailed transient logs are in ignored
`validation/runtime/path_isolation/`. Build products are in ignored
`build/path_isolation/`.
