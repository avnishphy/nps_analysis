# KinC_x60_4b LH2 livetime investigation

Current documentation edition: 2026-09-10. Full calculation: **43 production
runs, 183 updated segments, 77,838,384 events**. The user has now explicitly
authorized all documentation and LaTeX updates. The previous presentation hold
is superseded; the unresolved scientific conditions remain explicit.

## Current complete documentation

- [Package index and editing instructions](livetime_complete_20260910/README.md)
- [Technical report PDF](livetime_complete_20260910/KinC_x60_4b_LH2_livetime_report.pdf)
- [Technical presentation PDF](livetime_complete_20260910/KinC_x60_4b_LH2_livetime_presentation.pdf)
- [Editable presentation LaTeX](livetime_complete_20260910/presentation.tex)
- [Editable report LaTeX](livetime_complete_20260910/report.tex)
- [Complete Markdown report](livetime_complete_20260910/REPORT.md)
- [References and applicability](livetime_complete_20260910/REFERENCES.md)
- [All relevant plots and their scope](livetime_complete_20260910/PLOT_CATALOG.md)
- [43-run counts and ratios](livetime_complete_20260910/data/production_summary.csv)
- [Living research log](LIVETIME_RESEARCH_LOG.md) and [resume checkpoint](LIVETIME_REFRESH_CHECKPOINT.md)

## Current conclusion

The user already implemented the core broad-window pE/D method. The refinements
match event/scaler exposure and ending-interval current, and join only compatible
segments. Tight timing is a diagnostic: its rejected tags have pulse-phase
support. Phase association does not prove that EDTM caused the recorded trigger.

The dated logbook routing is adopted as the working configuration, following
the user's applicability confirmation. The direct EDTM copy bypasses NPS
cluster formation. Five TI6 coincidence runs and 38 TI4 HMS EL-REAL singles runs
therefore need different physical interpretations. Detailed HMS fanout and the
sent-scaler tap remain incompletely specified. PRE40/100/150/200 remain unmapped.

Runs 4303/4305 contain scaler exposure defects; 4305 is verified in the single
authorized raw file. No counter repair or yield/charge exclusion is adopted.
The exact total-livetime prescription remains unresolved. **CLT_physics is only
a proposed definition, never truth, baseline or a calibration target.** No
production correction, physics yield or efficiency CSV has been changed.

## Frozen evidence and previous editions

The complete edition includes copies of the current ROOT/raw/source/routing
evidence and historical numerical/presentation archives. Their original sibling
directories remain unchanged. The earlier 168-segment results and 56-slide deck
are historical; their claims about pending cache, unknown routing and TeX holds
must not be read as the current state. Earlier top-level documents are backed
up in `livetime_complete_20260910/previous_docs/`.

For exact data reproduction use the frozen full-cache manifest, not a fresh
live-cache inventory. No additional raw staging is authorized; the single
allowance has already been used for `nps_coin_4305.dat.0`.
