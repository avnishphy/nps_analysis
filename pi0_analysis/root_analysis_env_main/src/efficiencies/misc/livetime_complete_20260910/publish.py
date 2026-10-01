"""One-time publication of an audited, immutable edition and five doc entry points."""
from pathlib import Path
import hashlib,json,os,shutil
P=Path(__file__).resolve().parent
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
D=M/'livetime_complete_20260910'
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
assert not D.exists(),'Destination already exists; refusing to overwrite an edition.'
v=json.loads((P/'validation.json').read_text())
assert v['production_runs']==43 and not v['missing_links']
assert all(d['overfull_boxes']==0 and d['latex_warnings']==0 for d in v['pdfs'].values())
for name,h in json.loads((P/'original_doc_hashes.json').read_text()).items():
 assert sha(M/name)==h,'Live documentation changed since snapshot: '+name
 assert sha(P/'previous_docs'/name)==h
manifest=[]
for line in (P/'SHA256SUMS').read_text().splitlines():
 h,rel=line.split('  ',1);f=P/rel
 assert not Path(rel).is_absolute() and '..' not in Path(rel).parts
 assert sha(f)==h,'Staging file changed after validation: '+rel
 manifest.append((h,rel))
incoming=M/'.livetime_complete_20260910.incoming'
assert not incoming.exists(),'Unfinished publication exists; inspect before resuming.'
incoming.mkdir()
for h,rel in manifest:
 dst=incoming/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(P/rel,dst)
 assert sha(dst)==h,'Copy verification failed: '+rel
shutil.copy2(P/'SHA256SUMS',incoming/'SHA256SUMS')
os.replace(incoming,D)
for name in json.loads((P/'original_doc_hashes.json').read_text()):
 dst=M/name;tmp=M/('.'+name+'.livetime-doc-update')
 assert not tmp.exists()
 with tmp.open('xb') as out:out.write((P/'top_docs'/name).read_bytes())
 shutil.copymode(dst,tmp);os.replace(tmp,dst)
 assert sha(dst)==sha(D/'top_docs'/name)
print(json.dumps({'published':str(D),'verified_package_files':len(manifest),'updated_docs':list(json.loads((P/'original_doc_hashes.json').read_text())),'pdfs':v['pdfs']},indent=2))
