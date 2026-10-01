"""Freeze distributable checksums after build/validation and visual inspection."""
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent
v=json.loads((P/'validation.json').read_text())
assert v['pages']==56 and not v['out_of_page_text'] and not v['overfull_or_tex_errors']
v['visual_review']={'all_pages_in_contact_sheets':True,'additional_high_resolution_pages':[4,15,55],'checks':'Readable math, table rows, plot captions, diagram labels and page margins; no content overlap observed.'}
B=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc/livetime_KinC_x60_4b_lh2_20260909')
count=0
for line in (B/'SHA256SUMS').read_text().splitlines():
    sha,name=line.split('  ',1)
    assert hashlib.sha256((B/name).read_bytes()).hexdigest()==sha,name
    count+=1
v['original_analysis_archive_checksums_verified']=count
(P/'validation.json').write_text(json.dumps(v,indent=2)+'\n')
files=sorted(f for f in P.rglob('*') if f.is_file() and f.name!='SHA256SUMS' and '__pycache__' not in f.parts)
(P/'SHA256SUMS').write_text(''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+str(f.relative_to(P))+'\n' for f in files))
print('New package artifacts:',len(files),'unchanged original archive artifacts:',count)
