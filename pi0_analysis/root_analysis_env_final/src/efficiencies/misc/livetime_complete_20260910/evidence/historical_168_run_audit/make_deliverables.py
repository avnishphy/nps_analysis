"""Generate the report and presentation from the audited run results."""
from pathlib import Path
import csv,json,textwrap,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
P=Path(__file__).resolve().parent
R=json.loads((P/'audited_results.json').read_text());V=json.loads((P/'validation.json').read_text());valid=[r for r in R if r['status']=='CALCULATED'];prod=[r for r in valid if r['run_type']=='production'];byrun={r['run']:r for r in valid}
summary=list(csv.DictReader((P/'run_summary.csv').open()));cat=json.loads((P/'catalog_coverage.json').read_text())['runs'];catby={r['run']:r for r in cat}
suspects=list(csv.DictReader((P/'counter_suspect_intervals.csv').open())) if (P/'counter_suspect_intervals.csv').exists() else []
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])
def fmt(x):return f'{x:.8f}' if isinstance(x,(int,float)) else str(x)
ps=[r for r in prod if r['prescale_factor']>1]
replacements={
'SUMMARY':f"Completed **{V['calculated']} cached runs**, including all 43 production runs, from **{V['exported_segments']} segments and {V['events']:,} recorded events**. Five junk runs have no selected cache files: {', '.join(map(str,V['no_cache']))}. The same-coverage control reproduces all {V['saved_count_checks']} saved production count pairs and peaks exactly. {V['counter_suspect_intervals']} selected intervals in runs {', '.join(map(str,V['counter_suspect_runs']))} have recorded-event increments exceeding L1 increments by more than four. These defects remain visibly flagged; no ratio was forced to one.",
'COVERAGE_TABLE':table(['Run','Source','Cached','Missing catalogued segments'],[[r['run'],r['source'],','.join(map(str,r['cached_segments'])) or 'none',','.join(map(str,r['missing_catalog_segments']))] for r in cat if r['missing_catalog_segments']]),
'PRESCALES':table(['Run','Trigger / p','pE/D','pN/S','p(N-E)/(S-D)'],[[r['run'],f"{r['trigger']} / {r['prescale_factor']}",fmt(r['nominal']['EDTM_tight']),fmt(r['nominal']['CLT_all']),fmt(r['nominal']['CLT_physics'])] for r in ps]),
'COUNTER_AUDIT':table(['Run / component / interval','dt (s)','N','A','S','D','E'],[[f"{int(float(r['run']))} / {int(float(r['component']))} / {int(float(r['index']))}",f"{float(r['dt']):.6f}",int(float(r['N'])),int(float(r['A'])),int(float(r['S'])),int(float(r['D'])),int(float(r['E_tight']))] for r in suspects])+'\n\n'+table(['Run','Matched EDTM','If suspect intervals excluded','Matched physics CLT','If excluded'],[[r['run'],fmt(r['nominal']['EDTM_tight']),fmt(r['closure_exclusion_diagnostic']['EDTM']),fmt(r['nominal']['CLT_physics']),fmt(r['closure_exclusion_diagnostic']['CLT_physics'])] for r in valid if r['closure_exclusion_diagnostic']['removed_intervals']]),
'SENSITIVITIES':f"Across all calculated runs, the maximum absolute EDTM boundary-convention change is {V['max_boundary_delta']:.8g}; the maximum absolute +/-2 to +/-10-ns change is {V['max_width10_delta']:.8g}. These are observed variations, not adopted systematic errors. Every per-run variation is in `run_summary.csv` and the JSON/timing tables.",
'ALL_RUN_TABLE':'## All-run ratio table\n\n'+table(['Run','Type','Source','Segments','p','Saved NewGen','Same-cache NewGen','Matched EDTM','Physics CLT','Flags'],[[r['run'],r['run_type'],r['source'],r['cached_segments'] or 'none',r['prescale_factor'],r['saved_NewGen'],r['old_same_cache'],r['matched_tight'],r['CLT_physics'],r['flags'] or r['status']] for r in summary])
}
report=(P/'REPORT_TEMPLATE.md').read_text()
for k,v in replacements.items():report=report.replace('@@'+k+'@@',v)
if '@@' in report:raise RuntimeError('Unexpanded report placeholder')
(P/'REPORT.md').write_text(report)

pdf=PdfPages(P/'KinC_x60_4b_LH2_livetime_investigation.pdf');page=0;checks=[]
def canvas(title,source='Cache-only diagnostic | KinC_x60_4b / LH2 | 2026-09-09'):
    global page
    page+=1;fig=plt.figure(figsize=(12.8,7.2),facecolor='white')
    fig.add_artist(plt.Line2D([.055,.945],[.873,.873],transform=fig.transFigure,lw=2,color='#25859a'))
    h=fig.text(.055,.92,title,fontsize=26,weight='bold',color='#173b50',va='center')
    fig.canvas.draw()
    while h.get_window_extent(fig.canvas.get_renderer()).width>.89*fig.bbox.width:
        h.set_fontsize(h.get_fontsize()-1);fig.canvas.draw()
    fig.text(.055,.033,source,fontsize=9,color='#59636b');fig.text(.945,.033,str(page),ha='right',fontsize=10)
    return fig
def save(fig):
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    for obj in fig.texts:
        b=obj.get_window_extent(renderer)
        if b.x0<0 or b.x1>fig.bbox.width or b.y0<0 or b.y1>fig.bbox.height:raise RuntimeError('Text exceeds slide: '+obj.get_text())
    checks.append(dict(page=page,title=fig.texts[0].get_text(),text_bounds_pass=True))
    pdf.savefig(fig);plt.close(fig)
def textslide(title,paragraphs,source='',size=19):
    fig=canvas(title,source or 'Measured facts, conditional interpretation and unresolved coverage are distinguished')
    y=.80
    for p in paragraphs:
        lines=textwrap.wrap(p,width=78 if size<=19 else 72,break_long_words=False,break_on_hyphens=False)
        fig.text(.066,y,'\n'.join(lines),fontsize=size,va='top',linespacing=1.3)
        y-=.046*len(lines)+.035
    if y<.06:raise RuntimeError('Too many lines: '+title)
    save(fig)
def image_slide(title,name,caption,source='Figures generated from audited interval counts; ratios are not clipped'):
    fig=canvas(title,source);ax=fig.add_axes([.055,.19,.89,.65]);ax.imshow(plt.imread(P/'figures'/(name+'.png')));ax.axis('off')
    fig.text(.065,.105,caption,fontsize=12,va='center');save(fig)
def table_slide(title,headers,rows,caption,source='run_summary.csv and per-run interval tables'):
    fig=canvas(title,source);ax=fig.add_axes([.055,.24,.89,.58]);ax.axis('off')
    wrapped=[[textwrap.fill(str(c),width=max(10,102//len(headers)),break_long_words=False) for c in row] for row in rows]
    tab=ax.table(cellText=wrapped,colLabels=headers,cellLoc='center',bbox=[0,0,1,1]);tab.auto_set_font_size(False);tab.set_fontsize(13)
    for (i,j),c in tab.get_celld().items():
        if i==0:c.set_facecolor('#dfedf2');c.set_text_props(weight='bold',color='#173b50')
        elif i%2==0:c.set_facecolor('#f4f7f8')
    fig.text(.065,.12,textwrap.fill(caption,125),fontsize=12,va='center');save(fig)

textslide('Livetime can create artificial physics structure',[
    f"All 43 production runs checked; {V['calculated']} cached LH2 runs in total. {V['exported_segments']} selected segments, {V['events']:,} recorded events. No tape staging.",
    'Endpoint/current alignment changes the run dependence. Counter defects and prescale sampling can survive that alignment.',
    'The goal is a defensible normalization with explicit coverage, not a curve forced toward one.'
])
textslide('Scope: all LH2 entries, with their run types visible',[
    'Master list: 43 production, six junk and two efficiency runs. Five junk runs have no cache files; they receive a missing-data status, not a livetime.',
    '168 files selected from 253 cached replay variants. Match the existing updated-first, production-fallback policy; do not mix replay versions within a run.',
    'Main plots contain production runs. All 51 statuses and available-data results are in the appendix and CSV.'
],source='inventory.json; catalog_coverage.json; snapshotted master and efficiency CSVs')
table_slide('Partial cache coverage is explicit',['Run','Source','Cached segments','Missing'],[[r['run'],r['source'],','.join(map(str,r['cached_segments'])) or 'none',','.join(map(str,r['missing_catalog_segments']))] for r in cat if r['missing_catalog_segments']],
    'Read-only catalog listing; no retrieval. Unlisted uncached runs have no selected-source catalog entry.')
textslide('The original study tested four definitions',[
    'May 2025: whole-file EDTM survival; scaler L1/input computer LT; trigger-TDC computer LT for all events; and trigger-TDC computer LT with EDTM removed.',
    'The TDC computer estimates used scaler beam-on fractions to weight whole-file event counts. That does not establish direct interval matching.',
    'That talk concerns earlier luminosity scans. Its reported CODA problems are historical evidence, not a diagnosis of these runs.'
],source='Singh, NPS collaboration meeting May 2025, slides 3--11, 20, 23; prior REPORT.md')
table_slide('Previous estimators: the population matters',['Method','Numerator','Denominator'],[
    ['EDTM, old talk','Whole-file accepted EDTM','Accumulated EDTM'],
    ['CLT all, scalers','Current-cut L1 accepts','Current-cut trigger inputs'],
    ['CLT all, TDC','Whole-file trigger TDC x beam fraction','Beam-on trigger inputs'],
    ['CLT physics, TDC','Non-EDTM trigger TDC x beam fraction','Trigger inputs - sent EDTM']],
    'Prescale factors and exact cuts must also match. Proxy beam fractions do not fix endpoints.',
    source='Original talk slide 23; compute_luminosity_scaler.cxx inspected in run-4398 archive')
textslide('The older database CLT is a component comparison',[
    'The active older nps_livetimes.h reads CPU_LT from a database; its previous calculator is commented out.',
    'For 4398, database CPU_LT=0.9992 agrees at stored precision with the earlier >2-uA matched physics CLT=0.9992068726.',
    'That agreement does not certify total electronic/NPS coverage. The master CSV\'s separate Computer_live_time=0.992 has different unresolved provenance.'
],source='Prior run-4398 REPORT.md and archived nps_livetimes.h; do not conflate database columns')
textslide('Current NewGen: reproduced before changing it',[
    'NewGen = p E_old / D_old. E_old uses event-stored current and a run-level raw-TDC peak window. D_old uses ending-scaler current and per-file differences.',
    'All 42 production runs with unchanged saved coverage reproduce numerator, denominator and peak exactly. Run 4398 now has more cached segments.',
    'The division is not the endpoint bug. Numerator and denominator must first describe the same covered intervals.'
],source='snapshot/compute_efficiencies_stuff.cxx; saved_reproduction_checks.csv')
image_slide('Saved versus new: the full production comparison','comparison_all',
    'Shading: prescaled runs. Red rings: selected counter-closure defects. No forced upper bound.')
image_slide('Factor-one runs reveal the smaller structures','comparison_unprescaled',
    'This panel excludes prescaled runs, so their large sampling residuals do not compress the scale.')
textslide('A common interval is the unit of calculation',[
    'Use TSH.evNumber and T.g.evnum. Difference consecutive compatible snapshots; count events in [b_i, b_next).',
    'Join adjacent segments only with continuous events and compatible cumulative counters. Missing files or counter/clock restarts create separate components.',
    'Remove unchanged terminal snapshots. Keep zero-clock records with changing counters visible as anomalies; they are not ordinary duplicates.'
],source='analyze.py; per-run interval CSVs; component-break and anomaly records')
textslide('Initial counts and final tails are different problems',[
    'A first cumulative snapshot can already contain pre-Go pulses and trigger inputs. Treat it as the baseline for differences.',
    'A terminal scaler row may advance its event number without adding any counters or time. Events beyond the last real snapshot lack a matching denominator.',
    'A compatible next file can recover that tail. Do not drop every file tail blindly, and never bridge a missing segment.'
],source='Run 4398: prior endpoint investigation and THcScalerEvtHandler::End() evidence')
table_slide('What fixed the original run-4398 segment-0 excess',['Count population','Accepted E','Sent D','E/D'],[
    ['Original NewGen','19341','19310','1.00160539'],
    ['Covered, original raw window','19292','19310','0.99906784'],
    ['Covered, corrected +/-2 ns','19290','19310','0.99896427']],
    'The unmatched tail contained 1,264 events, including 49 raw-window EDTM candidates.',
    source='Earlier segment-0 endpoint check: real boundary 541321; duplicate terminal row 542585')
r=byrun[4398]
table_slide('Run 4398: separate a coverage change from a method change',['Calculation','Segments','EDTM ratio'],[
    ['Saved NewGen','0 only',fmt(float(r['saved']['NewGen_EDTM_livetime']))],
    ['NewGen recomputed','0--5',fmt(r['old_ratio'])],
    ['Matched raw timing','0--5',fmt(r['nominal']['EDTM_raw'])],
    ['Matched +/-2 ns','0--5',fmt(r['nominal']['EDTM_tight'])]],
    'The saved-to-new difference includes extra data. The same-cache comparison isolates the method.')
textslide('Current alignment and helicity are separate issues',[
    'The current stored on T normally belongs to the preceding scaler snapshot. Use the ending snapshot current for the interval being counted.',
    'Each segment retains the current helper\'s +/-15% peak window. At a cross-file interval, the ending file defines that window.',
    'TSHelH.evcount and TSH.evcount are different counters. The existing beam-time/gated diagnostic mixes them; main NewGen and this study remain current-only.'
],source='good_event_selection_helper.h; prescale_beamtime_helper.h; selection_checks.csv')
image_slide('Raw channels are not corrected nanoseconds','timing_peaks',
    'NewGen +/-500 raw channels is about +/-48.83 ns. The new corrected-time peak is found per run.')
table_slide('State the prescale and the counted population',['Quantity','Expression','Necessary interpretation'],[
    ['Pulser survival','p E / D','Representative delivered/accepted pulse sample'],
    ['All-trigger CLT','p N / S','One configured trigger; compatible coverage'],
    ['Physics CLT','p (N-E) / (S-D)','One EDTM input per D; same prescale population'],
    ['L1/input ratio','p A / S','Only the stages between those scaler taps']],
    'N = recorded; E = accepted EDTM subset; S = trigger input; D = sent EDTM; A = L1 accepts.')
image_slide('Prescaled runs require their own interpretation','prescales',
    'Factors come from the metadata convention; stored parameters do not independently prove TI settings.')
textslide('Stored prescales are useful evidence, with a clear limit',[
    '40 embedded Run_Data records contain populated prescales agreeing with the master metadata. The other 128 contain all-zero placeholders.',
    'All prescale set/read flags are zero. A later segment\'s all-zero record does not mean that its trigger was disabled.',
    'Every recorded event contains the configured trigger bit. This supports the selection, but is not independent proof of run-start TI settings or EDTM routing.'
],source='embedded_prescale_checks.csv; validation.json; per-run trigger-mask checks')
r=byrun[4259];n=r['nominal']
table_slide('Run 4259: endpoint matching is not the whole explanation',['Count / ratio','Value'],[
    ['TRIG / prescale factor','6 / 5'],['E / D',f"{n['E_tight']} / {int(n['D'])}"],['p E / D',fmt(n['EDTM_tight'])],['p N / S',fmt(n['CLT_all'])],['p (N-E) / (S-D)',fmt(n['CLT_physics'])]],
    'E/D is below one. The excess appears after multiplying by p; do not call it E > D.')
textslide('These are correlated ratios, not three independent tests',[
    'L_all = (1 - D/S) L_physics + (D/S) L_EDTM. This identity follows from shared counts.',
    'If prescaled EDTM is overrepresented, the subtracted physics estimate moves downward at fixed all-trigger ratio.',
    'Periodic pulser/prescale correlations and sampling fluctuations remain possible. The observed excess alone does not identify their cause.'
],source='Exact count identity; run4259.json; historical periodic-pulser literature is conditional')
table_slide('Counter defects: a plausible run average can hide them',['Run','dt (s)','N','A','S','D','E'],[[int(float(r['run'])),f"{float(r['dt']):.6f}",int(float(r['N'])),int(float(r['A'])),int(float(r['S'])),int(float(r['D'])),int(float(r['E_tight']))] for r in suspects],
    'Selected intervals with N-A > 4. Tiny scaler time accompanies many recorded events.')
image_slide('Run 4305: this is a counter mismatch, not ordinary busy','intervals_4305',
    'The discrepant interval adds 1,047 recorded events and 80 pulses with almost no scaler exposure.')
textslide('A diagnostic exclusion is not an authorized physics repair',[
    'The report shows the effect of excluding intervals with N-A > 4. This tolerance flags a mismatch; it is not a calibrated hardware cut.',
    'Those exclusions are not adopted. Ratios above one remain visible, and negative N-A is retained because it can represent loss after L1.',
    'If any interval is ultimately removed, its physics yield and charge must be handled consistently. Removing it only from livetime creates another bias.'
])
textslide('Run 4350: monotonicity alone is insufficient',[
    'A zero-clock record adds 35,540 EDTM counts. The next record jumps by 4,294,931,836, producing a cumulative value near 2^32.',
    'Both occur below the nominal current window. Keep their audit evidence; do not interpret them as physical pulser bursts.',
    'EDTM cumulative state restarts across segments 3->4; the clock restarts across 4->5. These boundaries are not stitched.'
],source='run4350.json: scaler_anomalies and component_breaks; recorded values, cause unresolved')
image_slide('Real loss bursts must remain in the average','intervals_4398',
    'Example: S=2082, N=A=1799, D=80, E=69 near 820.13 s. This is different from N exceeding A.')
image_slide('Early production runs: inspect the run dependence','comparison_early',
    'Matched ratios are conditional component measurements; defects and prescales remain labeled.')
image_slide('Later production runs: inspect the run dependence','comparison_late',
    'A change toward one can be explained by coverage, but proximity to one is not validation.')
image_slide('Separate the contributions to the method change','decomposition',
    'Differences are sequential and expressed in percentage points; their sum is the old-to-new change.')
image_slide('A livetime error becomes a normalization error','normalization_effect',
    'Hold yield, charge and other factors fixed. This illustrates sensitivity, not a measured cross-section bias.')
textslide('How artificial structures enter physics results',[
    'Endpoint exposure depends on run length and file boundaries. Current lag depends on beam trips and ramps. Prescale sampling depends on configuration.',
    'These features can vary systematically across run groups, rates or kinematics. Dividing by an inconsistent L can create false slopes and steps.',
    'Do not choose a correction because it flattens a yield. First match populations, validate counters and establish which physical losses are covered.'
])
image_slide('Timing selection must survive a sensitivity check','timing_sensitivity',
    'The width scan measures sensitivity. Outlying candidates are not proven accidentals by this plot.')
image_slide('Boundary, timing and current changes are distinct','sensitivity_all',
    'Changing current stability changes the sample. Single-event latch ambiguity remains explicit.')
image_slide('Do not mistake a convenient error bar for an uncertainty','uncertainties',
    'Bernoulli and block spreads are conditional diagnostics. Neither includes unknown hardware coverage.')
textslide('PRE attribution remains unresolved',[
    'PRE40/100/150/200 channels are populated. The issue is which input trigger each h/p monitor represents, and its actual width and routing.',
    'Nonzero counts or rate correlation do not establish applicability. No PRE-based correction is used.',
    'EDTM_CP is a copy, not an accepted-pulser monitor. A historical unbuffered 250-us model also requires the actual TI configuration before use.'
],source='User clarification; archived map and earlier 4398 PRE diagnostics; Mack 1001-v2 is conditional')
textslide('What is needed before adopting total livetime',[
    'Establish run-period EDTM/PRE signal paths and TI prescale/buffering/readout settings; verify the exact replay behavior at anomalous intervals.',
    'Match the final physics yield and charge to the chosen good beam/helicity intervals. Current-only pulser survival is not that final integration.',
    'Determine upstream and NPS-path coverage. If EDTM already includes computer loss, multiplying it by CLT counts that loss twice.'
])
textslide('What this calculation delivers',[
    f"{V['calculated']} cached LH2 runs calculated; all 51 metadata entries accounted for. {V['saved_count_checks']} exact original-count reproductions on unchanged coverage.",
    'Common-interval ratios, timing/current/boundary studies, counter anomalies and coverage differences are reproducible from the archived tables and scripts.',
    'Production corrections remain unchanged. Hardware attribution and consistent yield/charge integration remain explicit conditions for adoption.'
],source='REPORT.md; run_summary.csv; validation.json; source_provenance.json')
textslide('Reproduce without staging',[
    'Copy this archive to writable scratch. Source the Hall C/NPS environment, then run export_cached.C for its three manifest batches.',
    'Run analyze.py, audit.py, make_figures.py and make_deliverables.py. REPORT.md contains the exact commands and prerequisites.',
    'Use the frozen manifest. Missing cache files must be reported, not retrieved automatically. SHA256SUMS verifies the compact installed artifacts.'
])
fig=canvas('Appendix: supplied applied-efficiency reference','Original user-supplied PNG, copied unchanged; numeric comparison uses the saved CSV')
ax=fig.add_axes([.06,.09,.88,.75]);ax.imshow(plt.imread(P/'snapshot/efficiency_multipanel_vs_run_KinC_x60_4b_lh2.png'));ax.axis('off');save(fig)
for start in range(0,len(summary),9):
    group=summary[start:start+9];rows=[]
    for r in group:
        def small(k):return f"{float(r[k]):.6f}" if r[k]!='' else '--'
        notes=[]
        if r['status']=='NO_CACHED_FILES':notes.append('missing')
        if int(r['closure_suspect_intervals'] or 0)>0:notes.append('check')
        if r['missing_catalog_segments'] and r['status']!='NO_CACHED_FILES':notes.append('partial')
        if 'cross_segment_counter_or_clock_restart' in r['flags']:notes.append('restart')
        if float(r['prescale_factor'] or 1)>1:notes.append('p sample')
        note=', '.join(notes)
        rows.append([r['run'],r['run_type'][:4],r['prescale_factor'] or '--',small('saved_NewGen'),small('old_same_cache'),small('matched_tight'),small('CLT_physics'),note])
    table_slide(f'Appendix: all LH2 entries ({start+1}--{start+len(group)})',['Run','Type','p','Saved','Old/cache','EDTM','Phys. CLT','Status'],rows,
        'EDTM=pE/D, +/-2 ns. "check" = closure defect. Full counts, precision, coverage and flags in run_summary.csv.')
pdf.close();(P/'slide_checks.json').write_text(json.dumps(checks,indent=2));print('Wrote report and',page,'slides')
