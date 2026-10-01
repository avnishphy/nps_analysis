"""Preserve historical packages and prepare the living-document handoff."""
from pathlib import Path
import json, shutil
P=Path(__file__).resolve().parent
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
origins=json.loads((P/'evidence_origins.json').read_text())
for key,name in [('historical_168_run_audit','livetime_KinC_x60_4b_lh2_20260909'),('historical_beamer','livetime_beamer_20260909')]:
 dst=P/'evidence'/key
 if not dst.exists():shutil.copytree(M/name,dst)
 origins[key]=str(M/name)
(P/'evidence_origins.json').write_text(json.dumps(origins,indent=2)+'\n')
entry='''
## 2026-09-10: complete documentation and LaTeX edition

The user explicitly requested all findings, references, relevant plots and
comparisons, EDTM capabilities/limitations, Hall C treatment, trigger-setting
dependence and updated TeX. This supersedes the earlier hold on presentation
updates. It does not establish a total-livetime prescription or authorize a
production correction change.

New edition: [livetime_complete_20260910/README.md](livetime_complete_20260910/README.md).
Its editable report and Beamer sources compile offline. The report includes
all43 production counts and ratios, nineteen current figures, a complete
reference catalog, failed alternatives, source/version limits and the detailed
reason for separating TI4 singles from TI6 coincidence interpretation. The
slides have a main discussion and technical backup with all-run tables.

No ROOT/raw data were reread and no staging was performed. The new eleven
comparison figure pairs read frozen per-run JSON and interval CSV inputs.
The seven full-cache figure pairs and the applied NewGen PNG are unchanged
copies. Historical168-segment figures and full earlier archives remain separately
labeled. The five pre-update top-level Markdown documents are backed up.

The core method was already the user's broad raw-window pE/D. The nominal
changes are common event/scaler exposure, the same ending-interval current
decision, and compatible segment joins. For4398, same six files:
E112586 ->112263 ->112658 ->112658, D112746 throughout; final raw0.9992194845.
The net+72 tags demonstrate why this is not merely endpoint removal. Narrow
timing remains diagnostic, with pulse-phase-supported rejected tags. CLT_physics
remains only a proposed definition and never supplies a truth baseline.

Working configuration: dated entries4170042 and4171556, user-confirmed applicable.
Five TI6 runs /38 TI4 runs, with matching sole-enabled metadata factors. EDTM
joins the NPS branch at CH downstream of cluster formation. ForTI6 it cannot
measure bypassed NPS formation losses; forTI4 an NPS hardware cluster is not
the selected-input requirement. HMS fanout, sent tap, causal acceptance and
population weighting still need evidence. PRE attribution is still missing.

The documentation explicitly distinguishes conventional electronic/computer/
total LT endpoints, block level from buffer threshold, periodic-pulser from
eligible-event weighting, and single-trigger prescaling from possible
multi-trigger OR/overlap behavior. Actual metadata has one enabled input per
run; mask overlap is not proof of several independent enabled accept paths.

4303/4305 counter defects remain visible. Raw4305 proves the checked deficit
already exists in raw; it does not localize the hardware/ROC mechanism or prove
raw4303 behavior. No repair, missing-pulse addition, yield/charge exclusion,
PRE correction, rate fit, forced unity or total-loss multiplication was adopted.
Same-file old-to-matched-raw fixed-yield/charge normalization sensitivity spans
-0.7666206% (4493) to+0.6910719% (4308), not a measured cross-section bias.

Reproduce the documentation from the edition directory:

```bash
bash build.sh
```

This compiles editable TeX without replacing prose. Optional
`python3 make_comparisons.py` regenerates the eleven new plots from frozen
outputs. `python3 generate_report.py` deliberately replaces generated article
sections from Markdown/CSV and therefore must not be run over unpreserved
direct LaTeX edits. No ROOT environment is needed for document builds.

Build/numeric/link/preservation results are in
[validation.json](livetime_complete_20260910/validation.json). Issued checksums
and file inventory cover the editable source, PDFs and evidence. PDF layouts
were rendered for visual inspection; current table rows and refinements were
checked against the frozen CSV/JSON. Old archive copies are byte-checked against
their originals. Later edits will intentionally change issued checksums.

Remaining scientific work: run-period ROC5 readout/build/configuration, detailed
injection/sent taps, justified pulse/physics weighting, complete yield/charge/
good-helicity matching, and a final uncertainty/total-loss prescription. No
job needs restarting. No additional raw staging is authorized.
'''
old=(P/'previous_docs/LIVETIME_RESEARCH_LOG.md').read_text()
old=old.replace('The investigation has four existing evidence packages:', 'Initial evidence packages (later follow-ups and current edition appear below):')
old=old.replace('Earlier source research and presentation; frozen pending a justified method','Earlier source research and presentation; frozen historical edition')
old=old.replace('- No LaTeX revision until the physical method is established. No numerical\n  repair, event exclusion or production code/CSV change has been adopted.', '- The user has now authorized all documentation and LaTeX updates. The earlier\n  hold is superseded. No numerical repair, event exclusion or production\n  code/CSV change has been adopted.')
marker='## Reading and maintaining this record'
old=old.replace(marker, '## Current documentation edition\n\n[Complete report, presentation, references and plots](livetime_complete_20260910/README.md)\nare now updated at the user\'s explicit request. The final entry below records\nthis documentation build; earlier dated entries retain their historical scope.\n\n'+marker)
(P/'top_docs/LIVETIME_RESEARCH_LOG.md').write_text(old+entry)
checkpoint=(P/'previous_docs/LIVETIME_REFRESH_CHECKPOINT.md').read_text()
checkpoint=checkpoint.replace('- No LaTeX/presentation update until the method is established. Old Beamer\n  archive\'s148 hashes and numerical archive\'s402 hashes verified unchanged.', '- The user now explicitly authorized all documentation and LaTeX updates.\n  The earlier hold is superseded; it is not evidence that the method is final.\n  Earlier immutable archives are preserved. New edition and checks below.')
checkpoint=checkpoint.replace('## Artifacts / execution', '''## Current documentation edition

- [Index/edit/build instructions](livetime_complete_20260910/README.md)
- [Report](livetime_complete_20260910/REPORT.md), [references](livetime_complete_20260910/REFERENCES.md), [plot catalog](livetime_complete_20260910/PLOT_CATALOG.md)
- [Report PDF](livetime_complete_20260910/KinC_x60_4b_LH2_livetime_report.pdf), [slides PDF](livetime_complete_20260910/KinC_x60_4b_LH2_livetime_presentation.pdf)
- Editable report.tex, sections/ and presentation.tex; run `bash build.sh` in
  the edition directory. No ROOT is required. Do not regenerate Markdown-to-TeX
  sections over direct user LaTeX edits without preserving them.
- All43 current count/ratio rows, nineteen current figures, complete frozen
  evidence/source packages, old editions and original top-doc backups included.
- Eleven new plots read existing summaries/intervals; seven current diagnostic
  pairs and the supplied NewGen PNG copied unchanged. No event reread or staging.
- Core pE/D method was already the user's. Same-file refinements concern
  exposure/current/joins. For4398 E112586 ->112263 ->112658 ->112658, D112746.
- No exact total-LT prescription, production replacement or new raw allowance.
- Build, numerical and preservation checks: `livetime_complete_20260910/validation.json`.

## Artifacts / execution''')
(P/'top_docs/LIVETIME_REFRESH_CHECKPOINT.md').write_text(checkpoint)
mapping={Path(src).name:'evidence/'+key for key,src in origins.items()}
def relocate(text):
 text=text.replace('livetime_complete_20260910/','')
 for src,dst in mapping.items():text=text.replace(src+'/',dst+'/')
 text=text.replace('(LIVETIME_REFRESH_CHECKPOINT.md)','(CHECKPOINT.md)')
 text=text.replace('(LIVETIME_RESEARCH_LOG.md)','(RESEARCH_LOG.md)')
 return text
(P/'RESEARCH_LOG.md').write_text(relocate(old+entry))
(P/'CHECKPOINT.md').write_text(relocate(checkpoint))
print('Prepared all five top-level docs, portable companions and historical copies.')
