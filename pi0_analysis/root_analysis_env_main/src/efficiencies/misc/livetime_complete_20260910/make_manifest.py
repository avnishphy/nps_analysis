"""Seal the issued documentation; excludes local rendering/cache scratch."""
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
def included(f):
 r=f.relative_to(P)
 return f.is_file() and '__pycache__' not in r.parts and not any(str(r).startswith(x) for x in ['build/render/','build/mpl/'])
def digest(f):return hashlib.sha256(f.read_bytes()).hexdigest()
fs=[f for f in sorted(P.rglob('*')) if included(f) and f not in [P/'artifact_inventory.json',P/'SHA256SUMS']]
items=[{'path':str(f.relative_to(P)),'bytes':f.stat().st_size,'sha256':digest(f)} for f in fs]
(P/'artifact_inventory.json').write_text(json.dumps({'edition':'2026-09-10','scope':'All issued files except this inventory and its checksum manifest; local visual-review rasters and matplotlib cache excluded.','files':items},indent=2)+'\n')
fs.append(P/'artifact_inventory.json');fs.sort()
(P/'SHA256SUMS').write_text(''.join(digest(f)+'  '+str(f.relative_to(P))+'\n' for f in fs))
print('Sealed',len(fs),'files;',sum(f.stat().st_size for f in fs),'bytes.')
