"""Assemble a self-contained documentation edition from frozen evidence."""
from pathlib import Path
import hashlib,json,shutil
P=Path(__file__).resolve().parent
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
for d in ['figures','data','references','evidence','previous_docs','sections','build']:(P/d).mkdir(exist_ok=True)
def copytree(src,dst):
 if not dst.exists():shutil.copytree(src,dst)
archives={
 'full_cache':'livetime_refresh_20260909',
 'run4398':'livetime_4398_check_20260909',
 'coda':'coda_livetime_research_20260910',
 'apps_search':'nps_roc_apps_research_20260910',
 'npslib':'npslib_roc_research_20260910',
 'routing':'nps_trigger_routing_20260910',
 'refinements':'nps_livetime_refinements_20260910'}
for dest,src in archives.items():copytree(M/src,P/'evidence'/dest)
copytree(M/'livetime_4398_sources_20260909',P/'references/initial_sources')
copytree(M/'livetime_beamer_20260909/sources',P/'references/later_sources')
copytree(M/'livetime_beamer_20260909/figures',P/'figures/historical_168_segments')
for f in (M/'livetime_refresh_20260909/figures').iterdir():shutil.copyfile(f,P/'figures'/f.name)
shutil.copyfile(M/'livetime_refresh_20260909/snapshot/efficiency_multipanel_vs_run_KinC_x60_4b_lh2.png',P/'figures/applied_NewGen.png')
for n in ['production_summary.csv','acceptance_tests.csv','saved_reproduction_checks.csv','selection_checks.csv','source_provenance.json','combined_residual_diagnostic.json']:
 shutil.copyfile(M/'livetime_refresh_20260909'/n,P/'data'/n)
shutil.copyfile(M/'nps_trigger_routing_20260910/run_trigger_interpretation.csv',P/'data/run_trigger_interpretation.csv')
docs=['LIVETIME_4398_INVESTIGATION.md','LIVETIME_KinC_x60_4b_LH2.md','LIVETIME_REFRESH_CHECKPOINT.md','LIVETIME_RESEARCH_LOG.md','SOURCES_LIVETIME_4398.md']
old={}
for n in docs:
 f=M/n;b=f.read_bytes();(P/'previous_docs'/n).write_bytes(b);old[n]=hashlib.sha256(b).hexdigest()
(P/'original_doc_hashes.json').write_text(json.dumps(old,indent=2)+'\n')
for src,dest in [('SOURCES_LIVETIME_4398.md','original_source_catalog.md'),('livetime_beamer_20260909/RESEARCH.md','later_source_catalog.md')]:shutil.copyfile(M/src,P/'references'/dest)
origin={dest:str(M/src) for dest,src in archives.items()}
(P/'evidence_origins.json').write_text(json.dumps(origin,indent=2)+'\n')
print('Prepared',len(archives),'frozen evidence packages; source collections, plots, tables and previous-doc backups.')
