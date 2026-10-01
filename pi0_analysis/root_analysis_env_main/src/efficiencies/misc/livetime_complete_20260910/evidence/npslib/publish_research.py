"""Install source-search evidence and guarded living Markdown updates only."""
from pathlib import Path
import hashlib,json,shutil
P=Path(__file__).resolve().parent
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
D=M/P.name
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
old=json.loads((P/'original_doc_hashes.json').read_text())
assert not D.exists(),D
for n,h in old.items():
 assert (not (M/n).exists()) if h is None else sha(M/n)==h,('changed',n)
files=sorted(f for f in P.rglob('*') if f.is_file() and f.name not in ['SHA256SUMS','publish_receipt.json'])
(P/'SHA256SUMS').write_text(''.join(f'{sha(f)}  {f.relative_to(P)}\n' for f in files))
shutil.copytree(P,D)
for f in files:assert sha(D/f.relative_to(P))==sha(f)
for n,h in old.items():
 assert (not (M/n).exists()) if h is None else sha(M/n)==h,('changed',n)
 shutil.copyfile(P/n,M/n)
 assert sha(M/n)==sha(P/n)
receipt=dict(archive=str(D),verified_artifacts=len(files),markdown_updated=list(old))
(P/'publish_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
