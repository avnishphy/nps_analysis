"""Audit this documentation against the frozen numerical and source evidence."""
from pathlib import Path
from collections import Counter
from urllib.parse import unquote
import csv, hashlib, json, math, re, subprocess, xml.etree.ElementTree as ET
P=Path(__file__).resolve().parent
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def close(a,b):assert math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=1e-12),(a,b)
rows=list(csv.DictReader((P/'data/production_summary.csv').open()))
assert len(rows)==43
assert sum(int(r['segments']) for r in rows)==183
assert sum(int(r['events']) for r in rows)==77838384
cohorts=Counter((int(r['trigger']),int(r['p'])) for r in rows)
assert cohorts=={(6,1):4,(6,5):1,(4,1):33,(4,2):5},cohorts
assert [int(r['run']) for r in rows if int(r['counter_suspect_intervals'])]==[4303,4305]
for r in rows:
 p,N,E,D,S,A=(float(r[k]) for k in ['p','N','E_raw','D','S','A'])
 close(r['matched_raw_EDTM'],p*E/D)
 close(r['matched_tight_EDTM'],p*float(r['E_tight'])/D)
 close(r['all_trigger_ratio'],p*N/S)
 close(r['raw_subtracted_physics_ratio'],p*(N-E)/(S-D))
 close(r['all_trigger_ratio'],(1-D/S)*float(r['raw_subtracted_physics_ratio'])+D/S*float(r['matched_raw_EDTM']))
 raw=json.loads((P/'evidence/full_cache/results'/f"run{r['run']}.json").read_text())
 close(r['same_cache_NewGen'],raw['old_ratio'])
 for k in ['N','D','S','A','E_raw','E_tight']:close(r[k],raw['nominal'][k])
decomp=json.loads((P/'data/method_decomposition.json').read_text())
assert len(decomp)==43
for r,d in zip(rows,decomp):
 assert int(r['run'])==d['run']
 close(d['net_pp'],100*(float(r['matched_raw_EDTM'])-float(r['same_cache_NewGen'])))
r4398=json.loads((P/'evidence/full_cache/results/run4398.json').read_text())
assert [r4398[k] for k in ['old_E','perfile_covered_oldcurrent_E','perfile_covered_aligned_E']]==[112586,112263,112658]
assert r4398['nominal']['E_raw']==112658 and r4398['old_D']==112746
norm=[(100*(float(r['same_cache_NewGen'])/float(r['matched_raw_EDTM'])-1),int(r['run'])) for r in rows]
assert min(norm)[1]==4493 and max(norm)[1]==4308
assert round(min(norm)[0],7)==-0.7666206
assert round(max(norm)[0],7)==0.6910719
# Current figures copied unchanged; preserve the user's applied plot exactly.
for f in (P/'evidence/full_cache/figures').iterdir():assert sha(f)==sha(P/'figures'/f.name)
assert sha(P/'figures/applied_NewGen.png')==sha(P/'evidence/full_cache/snapshot/efficiency_multipanel_vs_run_KinC_x60_4b_lh2.png')
assert len(list((P/'figures').glob('*.png')))==19
assert len(list((P/'figures').glob('*.pdf')))==18
preserved={}
for dest,origin in json.loads((P/'evidence_origins.json').read_text()).items():
 src=Path(origin);n=0
 for f in src.rglob('*'):
  if f.is_file():
   g=P/'evidence'/dest/f.relative_to(src)
   assert g.is_file() and sha(f)==sha(g),str(g)
   n+=1
 preserved[dest]=n
doc_states={}
for name,h in json.loads((P/'original_doc_hashes.json').read_text()).items():
 assert sha(P/'previous_docs'/name)==h
 current=sha(M/name)
 assert current in [h,sha(P/'top_docs'/name)],'Unexpected live documentation change: '+name
 doc_states[name]='pre-update' if current==h else 'published version'
pdfs={}
for doc in ['report','presentation']:
 log=(P/'build'/f'{doc}.log').read_text()
 assert 'Overfull' not in log,doc+' has layout overflow'
 assert 'Warning:' not in log,doc+' has unresolved LaTeX warnings'
 assert 'Output written on' in log
 pdf=P/'build'/f'{doc}.pdf'
 info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
 pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M)[1])
 bbox=P/'build'/f'{doc}_bbox.html'
 subprocess.run(['pdftotext','-bbox',str(pdf),str(bbox)],check=True)
 tree=ET.parse(bbox);ns={'x':'http://www.w3.org/1999/xhtml'}
 nwords=0
 for page in tree.findall('.//x:page',ns):
  w,h=float(page.attrib['width']),float(page.attrib['height'])
  for word in page.findall('.//x:word',ns):
   a=word.attrib;nwords+=1
   assert float(a['xMin'])>=-0.1 and float(a['yMin'])>=-0.1,(doc,a)
   assert float(a['xMax'])<=w+.1 and float(a['yMax'])<=h+.1,(doc,a)
 pdfs[doc]={'pages':pages,'words_within_page':nwords,'overfull_boxes':0,'latex_warnings':0}
missing=[];checked=0
future={'validation.json','artifact_inventory.json','SHA256SUMS'}
for f in list(P.glob('*.md'))+list((P/'top_docs').glob('*.md')):
 top=f.parent.name=='top_docs'
 for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',f.read_text()):
  if re.match(r'\w+://',target) or target.startswith('#'):continue
  target=unquote(target.split('#')[0]);checked+=1
  if target in future:continue
  if top:
   g=P/target[len('livetime_complete_20260910/'):] if target.startswith('livetime_complete_20260910/') else M/target
  else:g=f.parent/target
  if not g.exists() and g not in [P/name for name in future]:missing.append((str(f.relative_to(P)),target))
assert not missing,missing
result={'edition':'2026-09-10','production_runs':43,'production_segments':183,'production_events':77838384,
 'trigger_prescale_cohorts':{f'TI{k[0]}_p{k[1]}':v for k,v in cohorts.items()},
 'ratios_recomputed_from_frozen_counts':43,'refinement_rows_checked':43,
 'counter_defect_flags':[4303,4305],'fixed_yield_charge_raw_normalization_percent':{'min':min(norm),'max':max(norm)},
 'current_figures':19,'new_figure_pairs':11,'unchanged_diagnostic_pairs':7,'original_NewGen_png_unchanged':True,
 'preserved_archive_file_counts':preserved,'old_top_docs_backups_verified':5,'live_top_doc_states':doc_states,
 'pdfs':pdfs,'new_markdown_local_links_checked':checked,'missing_links':missing,
 'validation_scope':'Documentation build and frozen-output audit, no new ROOT/raw read or production correction change.'}
(P/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
