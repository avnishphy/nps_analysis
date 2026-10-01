"""Validate completed outputs and retain unsupported counter intervals visibly.

The supplementary exclusion is NOT the adopted result: N-A>4 indicates a
counter/event coverage defect or unexplained latch behavior. Real S-N losses
are not excluded. Negative N-A can be actual loss after L1 and is retained.
"""
from pathlib import Path
import csv,json,re,hashlib
import numpy as np
from analyze import P,O,inventory,writecsv,div,read
summary=[];suspects=[];checks=[];selections=[];anomalies=[]
catalog={r['run']:r for r in json.loads((P/'catalog_coverage.json').read_text())['runs']}
savedsel={ (int(r['run_number']),int(r['segment_number'])):r for r in csv.DictReader((P/'snapshot/selection_report_KinC_x60_4b.csv').open())}
all_results=[]
for q in inventory:
    r=json.loads((O/f"run{q['run']}.json").read_text());all_results.append(r)
    if r['status'] not in ['CALCULATED','NO_CACHED_FILES','NO_VALID_SELECTION']:raise RuntimeError('Unfinished/failed run '+str(r))
    common=dict(run=q['run'],run_type=q['run_type'],status=r['status'],source=q['source'],cached_segments=';'.join(str(f['segment']) for f in q['files']),catalog_segments=';'.join(map(str,catalog[q['run']]['catalog_segments'])),missing_catalog_segments=';'.join(map(str,catalog[q['run']]['missing_catalog_segments'])),same_saved_files=q['same_saved_files'])
    z=dict(**common,trigger='',prescale_factor='',events='',D='',E='',S='',N='',A='',saved_NewGen='',old_same_cache='',matched_raw='',matched_tight='',CLT_all='',CLT_physics='',L1_over_input='',closure_exclusion_EDTM='',closure_exclusion_CLT='',closure_suspect_intervals='',EDTM_boundary_delta='',EDTM_width10_delta='',EDTM_stable_delta='',EDTM_block20_sigma='',current_low_min='',current_high_max='',raw_peak_channels='',corrected_peak_ns='',flags='')
    if r['status']!='CALCULATED':summary.append(z);continue
    anomalies.extend(dict(run=r['run'],**x) for x in r.get('scaler_anomalies',[]))
    ii=[{k:float(v) for k,v in a.items()} for a in csv.DictReader((O/f"run{r['run']}_intervals.csv").open())]
    nom=[x for x in ii if x['nominal']];bad=[x for x in nom if x['N']-x['A']>4];clean=[x for x in nom if x['N']-x['A']<=4]
    p=r['prescale_factor'];n=r['nominal'];badsum={k:sum(x[k] for x in bad) for k in ['N','A','D','S','E_tight','dt','charge_uC']}
    c={k:sum(x[k] for x in clean) for k in ['N','A','D','S','E_tight','dt','charge_uC']}
    c['EDTM']=div(p*c['E_tight'],c['D']);c['CLT_physics']=div(p*(c['N']-c['E_tight']),c['S']-c['D'])
    c['removed_intervals']=len(bad);c['removed_counts']=badsum
    r['closure_exclusion_diagnostic']=c
    r['max_interval_N_minus_A']=max(x['N']-x['A'] for x in nom)
    r['min_interval_N_minus_A']=min(x['N']-x['A'] for x in nom)
    if bad:r['flags'].append('recorded_exceeds_L1_interval')
    if catalog[q['run']]['missing_catalog_segments']:r['flags'].append('partial_catalog_coverage')
    for x in bad:suspects.append(x)
    # Check replay metadata evidence without treating stored parameters as TI proof.
    for d in r['segment_details']:
        old=savedsel.get((r['run'],d['segment']))
        selection_match=bool(old) and abs(d['low_uA']-float(old['current_min_uA']))<1e-7 and abs(d['high_uA']-float(old['current_max_uA']))<1e-7
        selections.append(dict(run=r['run'],segment=d['segment'],saved_selection_present=bool(old),current_window_match=selection_match,low=d['low_uA'],high=d['high_uA'],event_clock_previous_fraction=d['event_clock_matches_previous']/d['events'],event_current_previous_fraction=d['event_current_matches_previous']/d['events']))
        if old and not selection_match:raise RuntimeError('Saved selection mismatch '+str(d))
    old=r['saved']
    if old and r['same_saved_files']:
        ok=r['old_E']==int(old['NewGen_EDTM_num']) and abs(r['old_D']-float(old['NewGen_EDTM_den']))<1e-7 and r['raw_peak_channels']==float(old['NewGen_EDTM_peak'])
        checks.append(dict(run=r['run'],saved_counts_reproduced=ok,E=r['old_E'],D=r['old_D']))
        if not ok:raise RuntimeError('Saved original count mismatch '+str(r['run']))
    z.update(trigger=r['trigger'],prescale_factor=p,events=r['events'],D=n['D'],E=n['E_tight'],S=n['S'],N=n['N'],A=n['A'],saved_NewGen=float(old['NewGen_EDTM_livetime']) if old else '',old_same_cache=r['old_ratio'],matched_raw=n['EDTM_raw'],matched_tight=n['EDTM_tight'],CLT_all=n['CLT_all'],CLT_physics=n['CLT_physics'],L1_over_input=n['L1_over_input'],closure_exclusion_EDTM=c['EDTM'],closure_exclusion_CLT=c['CLT_physics'],closure_suspect_intervals=len(bad),EDTM_boundary_delta=n['EDTM_alternative']-n['EDTM_tight'],EDTM_width10_delta=next(t['ratio'] for t in r['timing'] if t['width_ns']==10)-n['EDTM_tight'],EDTM_stable_delta=r['stable']['EDTM_tight']-n['EDTM_tight'],EDTM_block20_sigma=next((t['EDTM_sigma'] for t in r['block_bootstrap'] if t['block_seconds']==20),float('nan')),current_low_min=min(s['low_uA'] for s in r['segment_details']),current_high_max=max(s['high_uA'] for s in r['segment_details']),raw_peak_channels=r['raw_peak_channels'],corrected_peak_ns=r['corrected_peak_ns'],flags=';'.join(r['flags']))
    summary.append(z)
writecsv(P/'run_summary.csv',summary);writecsv(P/'counter_suspect_intervals.csv',suspects);writecsv(P/'saved_reproduction_checks.csv',checks);writecsv(P/'selection_checks.csv',selections);writecsv(P/'scaler_anomalies.csv',anomalies)
(P/'audited_results.json').write_text(json.dumps(all_results,indent=2))
# Verify the frozen files still exist with the inventoried size and mtime.
filechecks=[]
for q in inventory:
    for f in q['files']:
        st=Path(f['path']).stat();same=(st.st_size==f['size'] and st.st_mtime==f['mtime'])
        filechecks.append(dict(run=q['run'],segment=f['segment'],path=f['path'],same_size_mtime=same))
        if not same:raise RuntimeError('Cache file changed '+f['path'])
writecsv(P/'file_stability_checks.csv',filechecks)
metadata={}
for log in [P/'pilot4259.log']+sorted(P.glob('export_batch*.log')):
    for block in re.split(r'(?=^START )',log.read_text(),flags=re.M):
        start=re.match(r'START (\d+) (\d+) (\S+)',block)
        ps=re.search(r'Prescale factors: \(1-12\)\s*([\d/\.\-]+)',block)
        fl=re.search(r'Prescales set/rd/req:\s*(\d+)\s+(\d+)\s+(\d+)',block)
        if start and ps and fl:
            run,seg=int(start[1]),int(start[2]);stored=list(map(float,ps[1].split('/')))
            result=next(r for r in all_results if r['run']==run)
            agrees=stored[result['trigger']-1]==result['prescale_factor'] and sum(x>0 for x in stored[:6])==1
            state='empty_placeholder' if all(x==0 for x in stored[:6]) else 'populated_agrees' if agrees else 'populated_conflict'
            metadata[(run,seg)]=dict(run=run,segment=seg,stored_prescales=ps[1],set_flag=int(fl[1]),read_flag=int(fl[2]),required_flag=int(fl[3]),agrees_with_metadata=agrees,evidence_state=state)
if len(metadata)!=len(filechecks):raise RuntimeError('Missing Run_Data metadata prints')
writecsv(P/'embedded_prescale_checks.csv',list(metadata.values()))
if any(x['evidence_state']=='populated_conflict' for x in metadata.values()):raise RuntimeError('Populated stored prescale conflicts with metadata')
v=dict(runs=len(inventory),calculated=sum(r['status']=='CALCULATED' for r in all_results),no_cache=[r['run'] for r in all_results if r['status']=='NO_CACHED_FILES'],exported_segments=len(filechecks),events=sum(r.get('events',0) for r in all_results),saved_count_checks=len(checks),counter_suspect_intervals=len(suspects),counter_suspect_runs=sorted({int(x['run']) for x in suspects}),max_boundary_delta=max(abs(float(x['EDTM_boundary_delta'])) for x in summary if x['EDTM_boundary_delta']!=''),max_width10_delta=max(abs(float(x['EDTM_width10_delta'])) for x in summary if x['EDTM_width10_delta']!=''))
v.update(embedded_prescale_checks=len(metadata),all_prescale_read_flags_zero=all(x['read_flag']==0 for x in metadata.values()),all_prescale_set_flags_zero=all(x['set_flag']==0 for x in metadata.values()),all_events_have_configured_trigger_bit=all(r.get('triggerbit_fraction',1)==1 for r in all_results),scaler_anomalies=len(anomalies),selected_scaler_anomalies=sum(x['in_current_window'] for x in anomalies))
v.update(populated_prescale_agreements=sum(x['evidence_state']=='populated_agrees' for x in metadata.values()),empty_prescale_placeholders=sum(x['evidence_state']=='empty_placeholder' for x in metadata.values()))
(P/'validation.json').write_text(json.dumps(v,indent=2));print(json.dumps(v,indent=2))
