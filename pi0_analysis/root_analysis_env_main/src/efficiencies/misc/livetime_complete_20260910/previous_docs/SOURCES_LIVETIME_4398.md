# Proposed run-4398 livetime evidence set

## New source extension and LaTeX presentation

The [2026-09-09 research note](livetime_beamer_20260909/RESEARCH.md) catalogs
new primary sources N01--N09, their limits, local analyzer evidence, retrieval
timestamps and hashes. The [Beamer PDF](livetime_beamer_20260909/KinC_x60_4b_LH2_livetime_Beamer.pdf)
and [LaTeX source](livetime_beamer_20260909/livetime.tex) incorporate these
alongside the original approved sources below. The key new constraints are
timing-acceptance validation and recovery of run-specific configuration records.
No source resolves PRE attribution or establishes a final total-livetime model.
The earlier source and numeric archives remain unchanged.

Cache-only extension completed for all KinC_x60_4b LH2 metadata entries:
[`livetime_KinC_x60_4b_lh2_20260909/REPORT.md`](livetime_KinC_x60_4b_lh2_20260909/REPORT.md).
It reuses this source set with its original scope limits, adds measured all-run
count/coverage evidence, and preserves the user-supplied current NewGen plot.
Do not extend the 4398 replay-job provenance to all runs or to both replay variants.

Collected 2026-09-09. User approved the source set and resumption in the new chat
("yeah, resume"). See `livetime_4398_check_20260909/REPORT.md` for the subsequent
six-segment investigation, evidence limits, and unresolved hardware attribution.
Archive: `livetime_4398_sources_20260909/`. Original PDFs, searchable text, metadata snapshots and SHA-256 manifests retained. Equations/diagrams require checking original PDFs.

## Core sources

| ID | Document | Intended role / limit |
|---|---|---|
| D00 | Avnish's [NPS luminosity presentation, May 2025](https://indico.jlab.org/event/946/contributions/16510/attachments/12602/20074/Luminosity_nps_collaboration_meeting_2025.pdf); `singh_nps_luminosity_202505.pdf` | User's prior investigation: starting point for definitions and unresolved discrepancies. |
| D01 | Eric Pooser, [Live Time Calculations, 1022-v1](https://hallcweb.jlab.org/DocDB/0010/001022/001/pooser_live_time.pdf); `doc1022_pooser_live_time.pdf` | Main formalism/checklist: slides 2-10, 14 onward. EDTM versus all-trigger/physics-trigger CLT; decoding, timing, current cuts and event-type pitfalls. F2/commissioning conditions require NPS validation. Metadata/body date May 3, 2019; title slide says 2018. |
| D02 | Dave Mack, [EDTM non-Poissonian bias correction, 1001-v2](https://hallcweb.jlab.org/DocDB/0010/001001/002/EDTMnonPoissonBiasCorrectionv2.pdf), updated September 20, 2019; `doc1001_EDTMnonPoissonBiasCorrectionv2.pdf` | Periodic-pulser and beam-ramp bias. Explicitly unbuffered; correction and assumed readout time cannot be transplanted without validating mode/conditions. |
| D03 | [Hall C Trigger Electronics, 1028-v5](https://hallcweb.jlab.org/doc-public/ShowDocument?docid=1028): [diagrams](https://hallcweb.jlab.org/DocDB/0010/001028/005/trigger_v3.pdf) and Brad Sawatzky's [EL_real update, December 3, 2021](https://hallcweb.jlab.org/DocDB/0010/001028/005/EL_real_trigger_update-03Dec2021.pdf); `doc1028_trigger_v3.pdf`, `doc1028_EL_real_trigger_update-03Dec2021.pdf` | Trace trigger stages, injection and scaler taps. Historical HMS/SHMS layout requires run-period NPS wiring confirmation. |
| D04 | Benjamin Raydo, [NPS Trigger/DAQ, February 2, 2023](https://wiki.jlab.org/cuawiki/images/9/90/NPS_Trigger_Status_Feb_2_2023.pdf); `nps_trigger_status_20230202.pdf` | NPS VTP clustering and waveform-readout architecture. Commissioning description does not establish run 4398 configuration. |
| D05 | Dave Mack, [Non-obvious deadtime sources, 1110-v2](https://hallcweb.jlab.org/DocDB/0011/001110/002/DetailedDeadtimeSourcesv2.pdf) and [Hodo 3of4 livetime combinatorics, 1063-v2](https://hallcweb.jlab.org/DocDB/0010/001063/002/Combinatorics%20in%20Estimating%20the%20Hodo%203of4%20Trigger%20Livetimev2.pdf); `doc1110_DetailedDeadtimeSourcesv2.pdf`, `doc1063_Combinatorics_in_Estimating_the_Hodo_3of4_Trigger_Livetimev2.pdf` | Electronics-loss mechanisms and assumptions in rate-based fallback estimates. Models require independent validation. |
| D06 | [SHMS instrument paper, June 5, 2024 circulated draft](https://mailman.jlab.org/pipermail/hallc/attachments/20240605/2611cf32/attachment-0001.pdf), sections 4.1-4.2 and 6.2; `shms_instrument_draft_20240605.pdf` | EDTM operating principle and buffered/unbuffered DAQ. Archived draft, not verified final journal version; SHMS coverage does not establish NPS coverage. |
| D07 | Joshua P. Crafts, [Semi-Inclusive pi0 Production Analysis, May 6, 2025](https://indico.jlab.org/event/946/contributions/16518/attachments/12601/20084/Edited_SemiInclusive_Pi0_Status_May2025_styled_clean.pdf), slide 13; `nps_pi0_status_20250506.pdf` | NPS precedent showing PRE100/PRE150 formula and TRIG6 CLT expression, EDTM still to be studied. Evidence of a convention, not proof it works. |
| D08 | [nps_replay](https://github.com/JeffersonLab/nps_replay), [hcana](https://github.com/JeffersonLab/hcana), [NPSlib](https://github.com/JeffersonLab/NPSlib) | After approval, establish pass2 revision and trace DEF/CUT/MAP/PARAM/report/decoder semantics. Only nps_replay README archived now; current heads cannot prove historical behavior. |

## Supporting sources (archived)

- [Hall C DAQ](https://hallcweb.jlab.org/wiki/index.php?title=Hall_C_DAQ): `hallc_daq.html`; prescale setting versus factor and run-start coda-flags provenance. Mutable, not run-specific.
- [Trigger Layout](https://hallcweb.jlab.org/wiki/index.php?title=Hall_C_Trigger_Layout), [Scalers](https://hallcweb.jlab.org/wiki/index.php?title=Scalers): `hallc_trigger_layout.html`, `hallc_scalers.html`; provenance and hardware-mapping context.
- [Archived NPS run plan, September 11, 2023](https://hallcweb.jlab.org/wiki/images/archive/1/10/20230911123553%21NPS_DVCS_RunPlan.pdf): `nps_runplan_archive_20230911.pdf`; planned trigger numbering, not proof of enabled triggers in 4398.
- Carlos Yero, [EDTM update, September 6, 2017](https://raw.githubusercontent.com/JeffersonLab/Hall-C-Trigger-Setup/master/edtm_studies/EDTM_report.pdf) and [Livetime studies, November 21, 2017](https://raw.githubusercontent.com/JeffersonLab/Hall-C-Trigger-Setup/master/edtm_studies/EDTM_LT_Studies.pdf): `yero_EDTM_report_20170906.pdf`, `yero_EDTM_LT_Studies_20171121.pdf`; historical commissioning studies linked by archived `hallc_livetime_studies.html`.
- [F2 update, February 17, 2022](https://indico.jlab.org/event/517/contributions/9310/attachments/7536/10488/F2_HallCwinterCollaboration_2022.pdf), slide 7: `f2_livetime_20220217.pdf`; limited EDTM statistics and modeled ELT precedent. Numerical correction not transferable to NPS.

## Unresolved evidence / prior code

- NPS RG1a Computer and Electronic Deadtime link is a redlink: page does not exist. Confirmed in archived `nps_analysis_index.html`; not an analysis prescription.
- Web-tool access to some DocDB pages failed; public pages/attachments subsequently downloaded without authentication using urllib. Download manifests record provenance.
- User's current implementation: `/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies`.
- Earlier implementation: `/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/src/nps_livetimes.h`.
- Both code locations and the earlier scaler/luminosity implementation were
  subsequently inspected; see `livetime_4398_check_20260909/REPORT.md`.
  Archived replay source and embedded DAQ metadata were recovered, but exact
  pass2 executable provenance and actual monitor wiring remain unresolved.

## Reproduce collection / verify

From requested misc directory, with Python 3 and network access:

```bash
python3 livetime_4398_sources_20260909/nps_livetime_collect.py
python3 livetime_4398_sources_20260909/nps_livetime_collect_attachments.py
python3 livetime_4398_sources_20260909/nps_livetime_collect_supplement.py
cd livetime_4398_sources_20260909
sha256sum -c SHA256SUMS
```

Collection commands were run from /tmp during preparation; scripts are preserved here. They write to `/tmp/nps_livetime_4398_sources_20260909`, retain URL/content type/bytes/SHA-256, check PDF magic and refuse differing file overwrites. `downloaded` establishes retrieval, not applicability. Installation verifies copied bytes. SHA-256 command above is provided for independent verification. No ROOT environment needed for this phase.
