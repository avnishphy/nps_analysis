"""Install the prepared new archive and two authorized Markdown updates.

Existing main note is hash-guarded; previous numerical and Beamer archives
are never modified. Run only after inspecting the concrete prepared package.
"""
from pathlib import Path
import json,hashlib,shutil
P=Path(__file__).resolve().parent;A=P/'archive'
D=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
dest=D/'livetime_refresh_20260909';note=D/'LIVETIME_KinC_x60_4b_LH2.md'
expected=json.loads((P/'publish_expected.json').read_text())
assert not dest.exists(),'New archive path already exists'
assert hashlib.sha256(note.read_bytes()).hexdigest()==expected['main_doc_sha256'],'Main note changed; merge required'
manifest=(A/'SHA256SUMS').read_text().splitlines()
for line in manifest:
    digest,rel=line.split('  ',1);assert hashlib.sha256((A/rel).read_bytes()).hexdigest()==digest
shutil.copytree(A,dest)
for line in manifest:
    digest,rel=line.split('  ',1);assert hashlib.sha256((dest/rel).read_bytes()).hexdigest()==digest
shutil.copy2(P/'MAIN_DOC.md',note)
shutil.copy2(P/'CHECKPOINT_FINAL.md',D/'LIVETIME_REFRESH_CHECKPOINT.md')
result=dict(archive=str(dest),files_verified=len(manifest),main_note=str(note),checkpoint=str(D/'LIVETIME_REFRESH_CHECKPOINT.md'))
(P/'publish_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
