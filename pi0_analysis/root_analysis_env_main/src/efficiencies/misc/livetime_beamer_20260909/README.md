# KinC_x60_4b LH2 livetime technical presentation

- [Compiled Beamer PDF](KinC_x60_4b_LH2_livetime_Beamer.pdf)
- [Editable LaTeX](livetime.tex)
- [New source investigation and exact definitions](RESEARCH.md)
- [Frozen calculation report](baseline/REPORT.md)
- [All 51 run statuses and numerical results](baseline/run_summary.csv)

56 slides: 27 main slides, followed by 29 backup slides. The main discussion
covers loss definitions, weighting, the four earlier estimators, current
NewGen, common-interval methodology, comparisons, prescales, counter defects,
PRE ambiguity and possible artificial normalization structure. Backup retains
the detailed variants, models, source evidence and all 43 production rows.

The 16:9 layout uses editable text/math, TikZ signal/interval diagrams and
vector analysis plots. The user's original applied multipanel PNG is preserved
unchanged. The old matplotlib-generated presentation remains in its original
archive; this PDF does not embed its slides as images.

## Build and verify

Requirements: TeX Live with pdflatex, Beamer, Latin Modern, TikZ and booktabs;
Poppler's pdftotext. No ROOT environment, network, cache access or Python is
needed to rebuild from the bundled evidence.

```bash
cd /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc/livetime_beamer_20260909
bash build.sh
rg 'Overfull|Error|^!' build/livetime.log
sha256sum -c SHA256SUMS
```

Expected: two successful LaTeX passes, 56 PDF pages, no overfull boxes or TeX
errors. `rg` returns exit 1 for no matches. Rebuilding can change PDF metadata
and checksums; verify the distributed archive before rebuilding. `validation.json`
records layout/numerical checks and `build/previews/` contains rendered pages
and contact sheets used for visual review.

The included `prepare.py` documents the initial evidence-copying and cache
availability scan; **it is not part of the build**. Running it requires the
original Hall C paths and refreshes availability. `fetch_sources.py` similarly
documents the public-source acquisition and is not needed to compile.

## Numerical scope

The baseline is the 2026-09-09 18:56 UTC inventory: 46 cached runs, including
43 production runs, 168 segments and 70,155,126 recorded events. It is not the
full newly staged calculation. The figures and ratios remain unchanged from
that audited baseline. The supplementary cache availability JSON records only
filenames/sizes later in the staging process; file appearance does not establish
complete or stable payloads.

The new +/-2 ns estimate is explicitly diagnostic. Literature supplies a
concrete caution that genuine EDTM can have corrected-time tails in another
configuration; NPS timing acceptance remains unmeasured here. PRE attribution,
run-start TI mode, prescale sampling, anomalous-counter cause and full
good-helicity/yield/charge consistency remain open. No PRE correction, clipping,
inferred scaler repair or production replacement is adopted.

## Reproduce the ROOT analysis later

Use the scripts and full commands in the sibling
`livetime_KinC_x60_4b_lh2_20260909/REPORT.md`. Its frozen manifest specifies the
original selected files. The compact archive excludes large temporary exported
column files, which must be regenerated if absent. Source ROOT before running:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; root-config --version'
```

For a refreshed cache calculation, create a new manifest/output directory and
retain the prior baseline. Never combine alternate updated/production replays
within a run. No jcache or other tape staging is performed by this package.
