# Reference catalog and applicability

References are evidence with stated scope, not authority to import a correction.
The source collections retain original files, text, retrieval metadata and hashes.
The August 2023 logbooks were supplied as text/URLs by the user; the web tool's
expired-authentication response prevented independent page retrieval. All other
retrieval status is recorded in the original manifests, including failed attempts
followed by successful public downloads. No source is silently replaced by a
newer web revision during this documentation build.

## Hall C and the user's earlier analysis

| ID | Reference | Use and limitation |
|---|---|---|
| D00 | [Singh, NPS luminosity, May 2025](https://indico.jlab.org/event/946/contributions/16510/attachments/12602/20074/Luminosity_nps_collaboration_meeting_2025.pdf) | User's earlier four definitions and normalization questions; historical sample differs |
| D01 | [Pooser, Live Time Calculations, DocDB1022-v1](https://hallcweb.jlab.org/DocDB/0010/001022/001/pooser_live_time.pdf) | HallC EDTM concept, ROC-lock and buffered distinction; title2018/body2019 date ambiguity preserved |
| D02 | [Mack, EDTM non-Poisson bias,1001-v2,2019-09-20](https://hallcweb.jlab.org/DocDB/0010/001001/002/EDTMnonPoissonBiasCorrectionv2.pdf) | Periodic pulser and beam-ramp effects; explicitly unbuffered, not an NPS correction |
| D03a | [HallC trigger layout,1028-v5](https://hallcweb.jlab.org/DocDB/0010/001028/005/trigger_v3.pdf) | Historical wiring/stages; not run-specific NPS PRE mapping |
| D03b | [Sawatzky, EL-REAL update,2021-12-03](https://hallcweb.jlab.org/DocDB/0010/001028/005/EL_real_trigger_update-03Dec2021.pdf) | EL-REAL electronics context; later run configuration must be established |
| D04 | [Raydo, NPS trigger status,2023-02-02](https://wiki.jlab.org/cuawiki/images/9/90/NPS_Trigger_Status_Feb_2_2023.pdf) | Commissioning trigger/readout architecture, not Feb 2024 readback |
| D05a | [Mack, non-obvious deadtime sources,1110-v2](https://hallcweb.jlab.org/DocDB/0011/001110/002/DetailedDeadtimeSourcesv2.pdf) | Leading edges, updating logic, recovery and shaping; examples not our measured hardware |
| D05b | [Mack, hodoscope3-of-4 combinatorics,1063-v2](https://hallcweb.jlab.org/DocDB/0010/001063/002/Combinatorics%20in%20Estimating%20the%20Hodo%203of4%20Trigger%20Livetimev2.pdf) | Correlation/model assumptions; no plane-rate correction fitted |
| D06 | [SHMS instrument circulated draft,2024-06-05](https://mailman.jlab.org/pipermail/hallc/attachments/20240605/2611cf32/attachment-0001.pdf) | Sections 4.1-4.2,6.2; draft/SHMS scope, not verified final journal version or NPS coverage |
| D07 | [Crafts, semi-inclusive pi0 status,2025-05-06](https://indico.jlab.org/event/946/contributions/16518/attachments/12601/20084/Edited_SemiInclusive_Pi0_Status_May2025_styled_clean.pdf) | Slide 13 PRE100/150 and TRIG6 precedent; convention is not validation |
| D08a | [JeffersonLab nps_replay](https://github.com/JeffersonLab/nps_replay) | Local archived replay tar/job provenance used; current repository head is not historical executable proof |
| D08b | [JeffersonLab hcana](https://github.com/JeffersonLab/hcana) | Inspected local scaler/event handling; source path/hash in evidence |
| D08c | [JeffersonLab NPSlib](https://github.com/JeffersonLab/NPSlib) | VTP/config/timing source; known installed HEAD does not prove pass 2 binary |

## NPS-period and pulser studies

| ID | Reference | Use and limitation |
|---|---|---|
| N01 | [Murphy, EDTM study,2022-06-28](https://redmine.jlab.org/attachments/download/1499/EDTM_Study_Report.pdf) | PionLT timing tails and competing triggers; do not transplant window/delay values |
| N02 | [Zhang, Deadtime and Efficiency,2024-07-18](https://indico.jlab.org/event/866/contributions/14918/subcontributions/283/attachments/11529/17843/Deadtime%20and%20Efficiency.pdf) | NPS definitions, discrepancies and historical 19-run slope; original units/causality unresolved |
| N03 | [Raydo, NPS Trigger/DAQ,2024-07-17](https://indico.jlab.org/event/866/contributions/14914/subcontributions/234/attachments/11505/17805/NPS_July17_2024_Raydo.pdf) | Trigger versus readout thresholds, EVIO config and VTP decisions; accepted-sample tests cannot recover rejected population |
| N04 | [Zhang, dead-time analysis,2024-02-19](https://redmine.jlab.org/attachments/download/2343/Dead%20time%20analysis_2024-02-19.pdf) | Historical raw-TDC/scaler investigation; examples not all current runs |
| N05 | [Zhang, updates,2024-03-03](https://redmine.jlab.org/attachments/download/2368/Updates_2024-03-03.pdf) | Automatic peak,2-uA cut, inferred scaler repair and terminal handling; not adopted defaults |
| N06 | [NPS issue837](https://redmine.jlab.org/issues/837) | Chronology and attachments; a resolved issue is not validation of this replay |
| N07a | [NPS electronics-deadtime issue866](https://redmine.jlab.org/issues/866) | Task context, not recovered hardware state |
| N07b | [DAQ event-blocking issue836](https://redmine.jlab.org/issues/836) | Task context, not effective whole-DAQ-mode proof |
| N08 | [HallC Trigger History](https://hallcweb.jlab.org/wiki/index.php/Trigger_History) | Links older wiring; historical entry 3532998 content not read |
| N09 | [HallC EDTM Pulser](https://hallcweb.jlab.org/wiki/index.php?title=Hall_C_EDTM_Pulser) | Pulser infrastructure and entry 3893007 link; not this run's measured clock/rate |

## Dated logbooks and CODA

| ID | Reference | Use and limitation |
|---|---|---|
| L01 | [Entry 4170042,2023-08-23](https://logbooks.jlab.org/entry/4170042) | User text: firmware capabilities, local NPS outputs and VLD path; capabilities not active settings |
| L02 | [Entry 4171556,2023-08-29](https://logbooks.jlab.org/entry/4171556) | User text and applicability confirmation: CH delays, EDTM injection and ROC1 TI mapping |
| L03 | User log4303,2024-02-10 | Original supplied attachment preserved under `evidence/full_cache/logbook`; no external entry URL supplied |
| L04 | User log4305,2024-02-10 | Original supplied attachment preserved; TI/ROC counters and interleaved status output |
| C00 | [CODA documentation portal](https://coda.jlab.org/drupal/) | Requested starting point; generic documentation must be tied to run configuration |
| C01 | [CODA single-crate example](https://coda.jlab.org/drupal/content/example-single-crate-configuration) | COOL_HOME and rol1/rol2 fields; example paths are not NPS paths |
| C02 | [CODA3.10 walkthrough](https://coda.jlab.org/drupal/system/files/coda/3.10/walkthrough/index.html) | Compiled ROL1 configuration and reproduction context |
| C03 | [VME module drivers](https://coda.jlab.org/drupal/content/vme-module-drivers) | TI3v11.3 documentation matches logged firmware label, not exact online library build |
| C04a | [TI header3v11.3](https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/tiLib_8h_source.html) | Buffer/selector/busy register masks |
| C04b | [TI configuration API](https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/group__Config.html) | External busy versus buffer-threshold control |
| C04c | [TI status API](https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/group__Status.html) | Timer units and latch/reset readout options; no proof of actual tiStatus execution |
| C04d | [TI master configuration API](https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/group__MasterConfig.html) | Latch requirement before timer getters |
| C05 | [TI hardware manual](https://coda.jlab.org/drupal/system/files/pdfs/HardwareManual/TI/TI.pdf) | Saved update2026-08-27; later features not assumed in 2024 |
| C06 | [TS hardware manual](https://coda.jlab.org/drupal/system/files/pdfs/HardwareManual/TS/TS.pdf) | Saved2017 Version 4 architecture; actual NPS master not identified |
| C07 | [CODA downloads](https://coda.jlab.org/drupal/content/downloads) | Public files also at /group/da/distribution/coda; matching local copies checked |
| C08 | [ROC landing page](https://coda.jlab.org/drupal/content/readout-controller-roc) | Inspected, no detailed implementation there; exact ROC5 source remains missing |

## Supporting historical sources retained

- [HallC DAQ](https://hallcweb.jlab.org/wiki/index.php?title=Hall_C_DAQ): prescale setting/factor and coda-flags context, mutable and not run-specific.
- [HallC Trigger Layout](https://hallcweb.jlab.org/wiki/index.php?title=Hall_C_Trigger_Layout) and [Scalers](https://hallcweb.jlab.org/wiki/index.php?title=Scalers): map context, not channel attribution by rate alone.
- [Archived NPS run plan,2023-09-11](https://hallcweb.jlab.org/wiki/images/archive/1/10/20230911123553%21NPS_DVCS_RunPlan.pdf): planned TI4/TI6 labels, now supplemented by L02.
- [Yero, EDTM report,2017-09-06](https://raw.githubusercontent.com/JeffersonLab/Hall-C-Trigger-Setup/master/edtm_studies/EDTM_report.pdf) and [livetime studies,2017-11-21](https://raw.githubusercontent.com/JeffersonLab/Hall-C-Trigger-Setup/master/edtm_studies/EDTM_LT_Studies.pdf): commissioning precedents, not current calibration.
- [F2 update,2022-02-17](https://indico.jlab.org/event/517/contributions/9310/attachments/7536/10488/F2_HallCwinterCollaboration_2022.pdf): limited EDTM statistics/model precedent, not transferable correction.
- [Historical wiring entry 3532998](https://logbooks.jlab.org/entry/3532998) and [pulser entry 3893007](https://logbooks.jlab.org/entry/3893007): linked leads, original content not inspected. These are not substitutes for L01/L02.
- The NPS RG1a Computer and Electronic Deadtime wiki entry was a redlink in the archived index; it supplied no prescription.

## Local implementation and data references

The full local source manifests retain paths, bytes and SHA256; selected
original Git identifiers are documented in the research log. These references
are local files, not web-derived claims:

- Original current NewGen: `evidence/full_cache/snapshot/newgen_edtm_livetime.h`, `compute_efficiencies_stuff.cxx`, `prescale_beamtime_helper.h`, `good_event_selection_helper.h`, `root_file_discovery.h` and the original plot/CSV.
- Original older luminosity and `nps_livetimes.h` source: `evidence/run4398/code_provenance.json` and its archived source snapshots; later `DT_Analyzer_new.C` under `references/later_sources/local`.
- Full-cache methodology: `evidence/full_cache/analyze.py`, `audited_results.json`, results, `acceptance_tests.py/CSV`, timestamp/phase/clock checks and all validation manifests.
- Raw 4305: `evidence/full_cache/raw_extract.C`, `raw_validate.py`, `raw_config_probe.py`, raw 4305 arrays/configs, logbook originals and raw checksum manifest.
- Actual-map evidence: `evidence/run4398/replay_snapshot/MAPS_db_HScalevt.dat` and other archived replay files; job/tar attribution remains limited to its documented scope.
- NPSlib versions 0269b41,bf8ec57,fddb8e9: evidence/npslib source snapshots, Git history and repeated-key model.
- Apps/hcana search: evidence/apps_search commands, negative results, file inventory and source hashes; offline sources are not online ROC5 implementations.
- Final routing interpretation and 43-run metadata join: evidence/routing; exact original-to-matched4398 comparison: evidence/refinements.

No source gives permission to identify PRE channels by resemblance to a pulse
width or fit a correction to CLT_physics. Applicability is assessed separately
from successful download, and missing hardware evidence remains visible.
