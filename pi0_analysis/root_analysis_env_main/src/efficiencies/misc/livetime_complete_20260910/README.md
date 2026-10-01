# Complete KinC_x60_4b LH2 livetime documentation

Edition: 2026-09-10. This edition replaces the earlier presentation for current
discussion. It uses the completed 43-production-run, 183-segment calculation
(77,838,384 events). Including cached controls/junk: 46 runs, 187 segments,
79,564,296 events. No ROOT or raw reread, staging, production correction change,
yield change or efficiency CSV replacement was performed for this edition.

The existing broad-window pE/D estimator remains the core method. Matching
event/scaler intervals and current assignment refines its bookkeeping. Its
interpretation depends on the injected path and the selected trigger. The
EDTM copy bypasses NPS cluster formation; TI4 and TI6 require different physical
interpretations. Counter defects in 4303/4305 remain unresolved. CLT_physics is
only a proposed definition, never a baseline or calibration target.

## Read and edit

| Artifact | Purpose |
|---|---|
| [Technical report PDF](KinC_x60_4b_LH2_livetime_report.pdf) | Complete narrative, all current comparison figures, counts, ratios and references |
| [Technical presentation PDF](KinC_x60_4b_LH2_livetime_presentation.pdf) | Main discussion followed by technical backup, all-run tables and source guides |
| [presentation.tex](presentation.tex) | Directly editable Beamer source; figures remain external |
| [report.tex](report.tex) and [report body](sections/report_body.tex) | Directly editable article source |
| [REPORT.md](REPORT.md) | Full parallel Markdown narrative |
| [RESEARCH_LOG.md](RESEARCH_LOG.md) | Detailed dated research history and documentation-build audit |
| [REFERENCES.md](REFERENCES.md) | Source-by-source role, URL and applicability limitation |
| [PLOT_CATALOG.md](PLOT_CATALOG.md) | Current and historical plots, population and provenance |
| [43-run publication table](data/production_summary.csv) | Exact current count and ratio inputs; raw-subtracted proposed ratio |
| [Trigger interpretation table](data/run_trigger_interpretation.csv) | All 43 metadata choices joined to dated ROC1 TI routing |
| [validation.json](validation.json) | Document build, numeric, link and preservation checks |
| [SHA256SUMS](SHA256SUMS) | Issued package checksums; verify before editing |

The report and presentation deliberately retain unresolved assumptions, failed
alternatives and values above one. The current figures do not substitute a
preferred-looking ratio for missing hardware evidence. Archived names such as
physics_CLT are preserved in data; the text states their conditional meaning.

## Build after editing LaTeX

From this directory, with TeX Live/pdflatex and Poppler installed:

```bash
bash build.sh
```

This compiles both documents three times and writes the two named PDFs. Logs
and extracted text are in `build/`. It does not regenerate prose, reread events,
call ROOT, require network access, or modify evidence. A failed build leaves
the precise pdflatex error in `build/report_pass*.stdout` or
`build/presentation_pass*.stdout`.

Edit `presentation.tex` directly for slides. For the report, edit `report.tex`,
`sections/report_body.tex`, `sections/references.tex` and `sections/tables.tex`.
The Markdown copy is not synchronized automatically with direct LaTeX edits.

Only if deliberately choosing Markdown as the source of truth:

```bash
python3 generate_report.py
bash build.sh
```

`generate_report.py` is a small project-specific converter, not a general
Markdown engine. It overwrites the report body, references and generated
tables from the Markdown/CSV files. Preserve direct LaTeX edits before using it.

## Regenerate the new plots without ROOT

Requires Python 3, numpy and matplotlib:

```bash
python3 make_comparisons.py
bash build.sh
```

This regenerates eleven new figure pairs from the frozen JSON/interval CSVs,
plus `data/method_decomposition.json`. The seven original full-cache diagnostic
figure pairs and unchanged applied NewGen PNG remain supplied copies; their
generator is `evidence/full_cache/make_refresh_figures.py`. See the plot catalog.
No fitted rate correction is produced.

## Verify the issued evidence

```bash
sha256sum --check --quiet SHA256SUMS
```

An unchanged package should exit successfully without output. Editing a listed
file creates an expected checksum difference. Keep one issued copy if you want
to compare later revisions; do not treat a modified copy as the original hash
baseline. `previous_docs/` preserves the five documentation entry points before
this update. `top_docs/` contains their replacement versions as published in
the parent misc directory.

## Evidence and historical scope

| Directory | Contents |
|---|---|
| `evidence/full_cache/` | Latest per-run/per-interval results, phase and clock checks, raw4305 summaries/arrays, user logs, scripts, source snapshots and manifests |
| `evidence/run4398/` | Initial segment0 and full-six-segment investigation, original code/replay snapshots and detailed boundary/current studies |
| `evidence/coda/` | Official CODA/TI/TS source snapshots, log register checks and version limits |
| `evidence/apps_search/` | Bounded local apps search, commands, source comparisons and negative results |
| `evidence/npslib/` | NPSlib source/history, repeated-key model and firmware-interpretation cautions |
| `evidence/routing/` | User-supplied dated logbook text, applicability confirmation and all-run routing join |
| `evidence/refinements/` | Exact original-to-matched 4398 count decomposition |
| `evidence/historical_168_run_audit/` | Earlier all-run numerical archive; historical sample only |
| `evidence/historical_beamer/` | Earlier 56-slide deck and its research/source package; historical sample only |
| `references/initial_sources/`, `references/later_sources/` | Original PDF/HTML/text collections and retrieval metadata |

Historical reports retain the conclusions known at their dates, including
superseded statements about cache coverage and unknown routing. Read the current
report/checkpoint first. Frozen archives are evidence records, not current
instructions to stage data, change corrections or hold presentation updates.
Some archived scripts contain original absolute paths and cannot run from a
relocated copy without explicit path adaptation. Publication scripts inside
old evidence packages are historical provenance and are not build steps.

Full ROOT reproduction requires the same cached replay files and the archived
manifest; large temporary column exports and the 13.3-GB raw file are not copied
into this documentation package. Follow the exact commands in the full-cache
and raw reports, in a fresh writable working copy. Load the Hall C environment
before ROOT/hcana. No additional raw staging is authorized; the single allowance
was already used for `nps_coin_4305.dat.0`.

## Reading the statistical comparisons

The raw-subtracted proposed ratio in `production_summary.csv` uses Eraw. The
legacy `run_summary.csv` CLT_physics column uses Etight. Do not mix them by name.
Paired resampling preserves the shared-count comparison; it is not a final
hardware-path systematic. The original independent-Poisson production error
function has not been replaced. Timing/phase association is not trigger causality.
Normalization-sensitivity plots hold yield and charge fixed; they are not
measurements of a cross-section bias.
