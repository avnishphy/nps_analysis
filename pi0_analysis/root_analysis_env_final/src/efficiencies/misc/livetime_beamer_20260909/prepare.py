"""Copy frozen evidence and create editable numerical slide inputs; no ROOT reads."""
from pathlib import Path
import shutil,json,hashlib,datetime
P=Path(__file__).resolve().parent
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
B=M/'livetime_KinC_x60_4b_lh2_20260909'
manifest=[]
def copy(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    manifest.append({'source':str(src),'file':str(dst.relative_to(P)),'sha256':hashlib.sha256(dst.read_bytes()).hexdigest()})
for name in ['audited_results.json','run_summary.csv','inventory.json','catalog_coverage.json','REPORT.md','validation.json','independent_reader_checks.json']:
    copy(B/name,P/'baseline'/name)
for f in (B/'figures').glob('*.pdf'):copy(f,P/'figures'/f.name)
copy(B/'snapshot/efficiency_multipanel_vs_run_KinC_x60_4b_lh2.png',P/'figures/applied_NewGen.png')
for name in ['singh_nps_luminosity_202505','doc1022_pooser_live_time','doc1001_EDTMnonPoissonBiasCorrectionv2','doc1110_DetailedDeadtimeSourcesv2','doc1063_Combinatorics_in_Estimating_the_Hodo_3of4_Trigger_Livetimev2','doc1028_trigger_v3']:
    for ext in ['pdf','txt']:
        f=M/'livetime_4398_sources_20260909'/f'{name}.{ext}'
        if f.exists():copy(f,P/'sources'/f.name)
for name in ['DT_Analyzer_new.C','DT_Analyzer.C']:
    copy(Path('/group/nps/yaopeng/Analysis/Deadtime/mass-production-pass1')/name,P/'sources/local'/name)
copy(M/'livetime_4398_check_20260909/code_evidence/compute_luminosity_scaler.cxx',P/'sources/local/compute_luminosity_scaler.cxx')
(P/'sources/local_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
rows=json.loads((B/'audited_results.json').read_text());prod=[r for r in rows if r['run_type']=='production' and r['status']=='CALCULATED']
def table(rs):
    s=r'\begin{tabular}{rrrrrr}\toprule Run & Trig./$p$ & NewGen & matched raw & matched $\pm2$ ns & physics CLT\\\midrule'+'\n'
    for r in rs:
        flag=r'\textsuperscript{!}' if r['closure_exclusion_diagnostic']['removed_intervals'] else ''
        n=r['nominal']
        s+=f"{r['run']}{flag} & {r['trigger']}/{r['prescale_factor']} & {r['old_ratio']:.6f} & {n['EDTM_raw']:.6f} & {n['EDTM_tight']:.6f} & {n['CLT_physics']:.6f}"+r'\\'+'\n'
    return s+r'\bottomrule\end{tabular}'
(P/'tables').mkdir(exist_ok=True)
for i in range(0,len(prod),11):(P/f'tables/production_{i//11+1}.tex').write_text(table(prod[i:i+11]))
copy(B/'snapshot/nps_dvcs_all_kins_main.csv',P/'baseline/master.csv') if (B/'snapshot/nps_dvcs_all_kins_main.csv').exists() else None
# Availability check only; preserve the analysis snapshot and selection policy.
inv=json.loads((B/'inventory.json').read_text());availability=[]
for r in inv['runs']:
    variants={}
    for v in ['updated','production']:
        fs=sorted(Path('/cache/hallc/c-nps/analysis/pass2/replays',v).glob(f"nps_hms_coin_{r['run']}_*_1_-1.root"))
        variants[v]=[{'path':str(f),'bytes':f.stat().st_size,'segment':int(f.name.split('_')[-3])} for f in fs]
    availability.append({'run':r['run'],'type':r['run_type'],'baseline_source':r['source'],'baseline_segments':[f['segment'] for f in r['files']],'current_cache':variants})
(P/'cache_availability_only.json').write_text(json.dumps({'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'note':'Read-only filenames/stat. No reanalysis or staging; prior ratios retain frozen coverage.','runs':availability},indent=2)+'\n')
print('Frozen production rows:',len(prod),'copied evidence files:',len(manifest))
print('Current files by replay variant:',{v:sum(len(r['current_cache'][v]) for r in availability) for v in ['updated','production']})
