"""Build the investigation presentation from measured full-run results."""
from pathlib import Path
import csv,json,textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
P=Path(__file__).resolve().parent
r=json.loads((P/'full_run_results.json').read_text()); prod=r['estimators'][2]
pdf=PdfPages(P/'Run4398_livetime_investigation.pdf');page=0
def canvas(title,ref='Run 4398 | updated replay | 2026-09-09'):
    global page
    page+=1;fig=plt.figure(figsize=(12.8,7.2),facecolor='white')
    heading=fig.text(.055,.91,title,fontsize=25,weight='bold',color='#153b59')
    fig.canvas.draw()
    while heading.get_window_extent(fig.canvas.get_renderer()).width > .89*fig.bbox.width:
        heading.set_fontsize(heading.get_fontsize()-1)
        fig.canvas.draw()
    fig.text(.055,.035,ref,fontsize=9,color='#555555');fig.text(.945,.035,str(page),ha='right',fontsize=10)
    return fig
def save(fig):
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    for obj in fig.texts:
        box=obj.get_window_extent(renderer)
        if box.x0<0 or box.x1>fig.bbox.width or box.y0<0 or box.y1>fig.bbox.height:
            raise RuntimeError('Text outside slide: '+obj.get_text()[:100])
    pdf.savefig(fig);plt.close(fig)
def textslide(title,paragraphs,ref='',size=18):
    fig=canvas(title,ref or 'Run 4398 | analysis findings and remaining evidence')
    y=.80
    for p in paragraphs:
        lines=textwrap.wrap(p,width=86 if size<18 else 82,break_long_words=True,break_on_hyphens=False)
        fig.text(.065,y,'\n'.join(lines),fontsize=size,va='top',linespacing=1.35)
        y-=.044*len(lines)+.035
    if y<.055:raise RuntimeError('Slide text overflows: '+title)
    save(fig)
def image_slide(title,path,caption,ref):
    fig=canvas(title,ref);ax=fig.add_axes([.07,.17,.86,.64]);ax.imshow(plt.imread(P/path));ax.axis('off')
    fig.text(.06,.095,caption,fontsize=12);save(fig)
def table_slide(title,columns,rows,foot):
    fig=canvas(title);ax=fig.add_axes([.055,.22,.89,.57]);ax.axis('off')
    tab=ax.table(cellText=rows,colLabels=columns,cellLoc='center',loc='center')
    tab.auto_set_font_size(False);tab.set_fontsize(13);tab.scale(1,2.2)
    for (i,j),cell in tab.get_celld().items():
        if i==0:cell.set_facecolor('#dceaf1');cell.set_text_props(weight='bold')
    fig.text(.065,.10,foot,fontsize=12);save(fig)

textslide('Run 4398: what the evidence establishes',[
    'All six updated replay segments were inspected: 3,222,394 contiguous recorded events.',
    'The segment-0 EDTM ratio above one has an event/scaler endpoint mismatch. Joining real snapshots across files recovers intermediate segment tails.',
    f'Fixed production-current intervals: EDTM survival = {prod["EDTM_tight"]:.8f}; physics computer livetime = {prod["CLTP_tight"]:.8f} under the stated timing/boundary convention.',
    'Total physics livetime remains conditional on actual EDTM coverage and DAQ configuration. PRE monitor trigger attribution is unresolved.'
])
textslide('Question and scope',[
    'Determine a defensible livetime from the updated run-4398 replay, despite concerns about EDTM and PRE monitor interpretation.',
    'The user clarified that the PRE40/100/150/200 issue is which trigger they monitor. Presence and nonzero values do not establish applicability.',
    'This presentation reports measured component ratios, implementation findings, and exact missing evidence. Production efficiency code and CSVs were not changed.'
])
textslide('Data and replay provenance',[
    'Requested input: /cache/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_4398_*_1_-1.root.',
    'Initially only segment 0 was cached; tape request 100230264 retrieved segments 1--5. All six exports completed.',
    'Archived SWIF job identifies nps_replay_pass2_github.tar.gz. Archived HEAD: 2d462a85b860c54785d8287e644619a1cb1acc19. Relevant file bytes and hashes are retained.',
    'Original job stdout/stderr were not found. The exact pass2 executable build is therefore not fully recovered.'
],size=17)
table_slide('Segment coverage',['Segment','Events','First event','Last event','Real snapshots'],[
    [str(q['segment']),f"{q['events']:,}",f"{q['first_event']:,}",f"{q['last_event']:,}",str(q['real_scaler_rows'])] for q in r['segments']
], 'No event-number gaps or duplicates across the six exported event trees.')
textslide('Definitions: keep stage boundaries explicit',[
    'N = recorded trigger events; E = accepted EDTM subset; S = trigger-input scaler counts; D = delivered EDTM scaler pulses.',
    'For one enabled trigger with factor 1: pulser survival E/D; all-trigger computer LT N/S; physics computer LT (N-E)/(S-D).',
    'Intentional prescaling, trigger/detector efficiency, electronics losses, DAQ busy, and offline cuts are different effects.',
    'An EDTM estimate that already includes DAQ loss must not be multiplied by computer livetime again.'
], 'D01 Pooser, slides 7--10, 14--15; run-specific assumptions stated in REPORT.md')
textslide('Counter provenance and signal paths',[
    'Archived map: EDTM = ROC 5 / slot 12 / channel 14; TRIG4 = ROC 5 / slot 11 / channel 13; L1 accept = ROC 5 / slot 10 / channel 16.',
    'EDTM_CP is a copy monitor, not an accepted-pulser channel. The two clock channels are also copies.',
    'The trigger-input/L1 ratio probes losses between those taps. EDTM probes only stages traversed after its actual injection points.',
    'The historical NPS plan identifies trigger 4 as HMS EL-REAL. That is context, not proof of run-specific wiring or NPS VTP coverage.'
], 'Archived scaler map; D03 printed page 14; archived NPS run plan printed page 5',size=17)
textslide('Embedded metadata: useful, but incomplete',[
    'Run_Data stores trigger-4 factor 1. Its prescale set/read flags are zero, so this is not independent proof of a decoded TI prescale event.',
    'Every recorded event has g.evtyp=1. Trigger masks vary and all include bit 3; event type alone cannot classify EDTM.',
    'Ten embedded FADC/VTP configuration strings were saved. No TI buffer level, block level or ROC-lock setting was located.',
    'VTP prescale settings are not TI prescale factors. Do not transplant an unbuffered 250-us correction into unknown DAQ conditions.'
], 'run_metadata.log; daq_config.txt; D02 Mack 1001-v2, slides 2, 4--7',size=17)
textslide('Why segment 0 gave EDTM > 1',[
    'Existing result: 19,341 accepted candidates / 19,310 delivered pulses = 1.00160539.',
    'The final TSH row advances evNumber to 542585, while its clock and all relevant scaler counts repeat the preceding snapshot at event 541321.',
    'Interval [541321,542585) contains 1,264 recorded events and 49 EDTM candidates, but no new scaler counts.',
    'Removing that uncovered tail gives 19,292 / 19,310 = 0.99906784 for the original raw-time window.'
], 'Segment-0 counts; inspected THcScalerEvtHandler::End(); endpoint_mismatch.png')
image_slide('The endpoint error is visible in cumulative counts','endpoint_mismatch.png',
    'The terminal jump adds accepted events without any new denominator counts.', 'Segment 0; all currents; raw-window EDTM')
textslide('Match intervals before forming ratios',[
    'Use common event boundaries: TSH.evNumber, TSHelH.evNumber, and T.g.evnum. The two scaler-tree evcount variables are separate record counters.',
    'Difference real cumulative snapshots. Join snapshots across adjacent segments and reject duplicated terminal records.',
    'Apply the ending-snapshot current to the same covered event interval. Event-stored scaler current is from the preceding snapshot.',
    'Keep charge integration and accepted physics yield on the same intervals. A livetime-only tail cut does not authorize unmatched yield normalization.'
], 'full_results.py; full_run_intervals.csv; boundary convention documented in REPORT.md',size=17)
image_slide('Joining six segments removes artificial endpoint jumps','full_run_alignment.png',
    'Real snapshots span 21.963170--4160.053707 s; primary cumulative counters have no resets.',
    'Full updated run; one terminal non-EDTM event excluded by the half-open convention')
table_slide('Measured full-run ratios', ['Selection','N','E, +/-2 ns','D','EDTM','Physics CLT'],[
    [q['selection'].replace('fixed production current','Production I').replace('all covered','All I').replace('stable production current','Stable I'),
     f"{q['N']:,}",f"{q['E_tight']:,}",f"{q['D']:,}",f"{q['EDTM_tight']:.8f}",f"{q['CLTP_tight']:.8f}"] for q in r['estimators']
], 'Production I is the fixed diagnostic window 33.3625--45.1375 uA; not a final total-LT prescription.')
textslide('Timing: channels and nanoseconds are different',[
    'The current implementation cuts +/-500 raw TDC channels, although its comment says ns. The archived conversion is 0.09766 ns/channel: +/-48.83 ns.',
    'Reference-subtracted EDTM timing peaks near 245.11508 ns. Segment 0 has two raw-window timing outliers beyond +/-2 ns.',
    'Full production-current numerator: 112658 in the raw window, 112650 at +/-2 ns, and 112652 at +/-5 and +/-10 ns.',
    'Use this observed timing sensitivity in the uncertainty study. Calling the outliers accidentals requires additional evidence.'
], 'THcTrigDet.cxx; archived thms_nps23.param; full_run_timing_sensitivity.csv',size=17)
image_slide('Raw timing is broad; corrected timing resolves outliers','edtm_timing.png',
    'Illustrative segment-0 timing plot; full-run cut sensitivity is in the accompanying CSV/PDF.',
    'Covered segment-0 production-current intervals')
textslide('Statistical precision and time variation',[
    f'Conditional independent-binomial scales: EDTM {prod["conditional_sigma_EDTM"]:.2g}; physics computer LT {prod["conditional_sigma_CLTP"]:.2g}. These are not final uncertainties.',
    'Time-block resampling (10, 20, 60 s) gives EDTM spread about 0.00013 and physics-CLT spread about 0.00011, revealing substantial time variation.',
    'Near clock 820.1256 s, one ~2-s interval loses 283 of 2082 trigger inputs and 11 of 80 EDTM pulses. Retain real losses in the run average.',
    'Changing the half-open boundary convention shifts EDTM by 0.00001774 and CLT by 0.00000072. None of these checks diagnoses unmonitored electronics.'
], 'full_run_block_bootstrap.json; full_run_intervals.csv; bootstrap is a sensitivity diagnostic',size=17)
image_slide('Losses are not uniform in time','full_run_loss_bursts.png',
    'Coincident trigger and EDTM losses support busy sensitivity, not proof of all upstream coverage.',
    'Fixed production-current intervals; no loss bursts removed')
textslide('PRE scalers: trigger identity remains unresolved',[
    'All PRE40/100/150/200 channels are populated. That answers availability, not which physical trigger each channel samples.',
    'Segment-0 HMS PRE150/PRE100 = 0.99067720; pPRE150/pPRE100 = 0.95547866. Different monitor choices imply very different corrections.',
    'Required before use: physical input trigger, cable/module routing, actual widths, overlap/retrigger behavior, and applicability to the desired physics population.',
    'A close rate correlation can suggest wiring candidates; it cannot establish trigger identity. No PRE-based correction is adopted.'
], 'User clarification; archived maps; D07 slide 13 is precedent, not run-4398 validation',size=17)
textslide('Implementation findings to address in a controlled change',[
    'Event/scaler coverage: reject unsupported tails or recover them by joining real snapshots across segments.',
    'Current alignment: event-stored current belongs to a prior interval; match boundaries explicitly.',
    'Counter mismatch: TSHelH evcount ranges are applied to TSH evcount. In segment 0 this gives 163.950216 s instead of 483.315153 s of covered current intervals.',
    'Completeness: any cached updated file is currently accepted as the run. Also correct timing-unit comments and justify numerator/denominator covariance.'
], 'code_evidence/; source files retained unmodified',size=17)
textslide('What to use, and what still needs evidence',[
    f'A defensible measured computer-livetime candidate for the stated fixed-current sample is {prod["CLTP_tight"]:.8f}, under the trigger/prescale and boundary assumptions.',
    f'EDTM pulse survival is {prod["EDTM_tight"]:.8f} with the stated cut. Agreement with CLT supports consistency of observed DAQ losses.',
    'A final total physics correction needs run-specific EDTM/PRE routing, TI buffering/readout information, and exact matching to the production charge/yield selection.',
    'Do not assign an unverified PRE correction, silently label computer LT as total LT, or clip a ratio above one.'
], 'Conclusions are conditional; unresolved hardware evidence is listed in REPORT.md',size=17)
textslide('Reproducibility and validation',[
    'Read REPORT.md first. Exact counts and interval records are in full_run_estimators.csv and full_run_intervals.csv.',
    'export_full_run.C reads the six updated ROOT files into diagnostic column caches. full_results.py recreates the combined tables, plots and sensitivity checks.',
    'Independent ROOT/uproot comparison agreed exactly for segment 0: 542584 x 92 event values and 355 x 696 scaler values.',
    'Source snapshots, original run-specific CSV extracts, DAQ metadata, archive checksums and presentation source are retained. No production source or CSV was changed.'
], 'ROOT environment: source /usr/share/Modules/init/csh, then /group/nps/singhav/setup.csh',size=17)
textslide('Reference and evidence index',[
    'D00: Avnish, NPS luminosity presentation, May 2025. D01: Pooser, Live Time Calculations, 1022-v1, especially slides 7--10 and 14--15.',
    'D02: Mack, non-Poissonian EDTM correction, 1001-v2, September 2019. D03: Hall C Trigger Electronics, 1028-v5. Applicability is conditional.',
    'Archived September 2023 NPS run plan: planned trigger numbering. D07: Crafts, May 2025, slide 13: legacy PRE-monitor convention.',
    'Local evidence: replay_provenance.json, code_provenance.json, inventory.json, run_metadata.log, daq_config.txt and all count/sensitivity tables.'
], 'Original PDFs and URL/SHA-256 manifests remain in the adjacent source archive.',size=17)
pdf.close();print(f'Wrote {page} slides to {P / "Run4398_livetime_investigation.pdf"}')
