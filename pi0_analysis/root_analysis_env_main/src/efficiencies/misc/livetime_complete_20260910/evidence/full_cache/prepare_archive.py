"""Prepare a compact immutable result package; never include large ROOT columns."""
from pathlib import Path
import shutil,hashlib,json
P=Path(__file__).resolve().parent;A=P/'archive'
if A.exists():raise RuntimeError('Archive staging already exists; inspect before replacing')
A.mkdir()
skip={'MAIN_DOC.md','CHECKPOINT_FINAL.md','CHECKPOINT.md','publish_expected.json','make_figures.py'}
for f in P.iterdir():
    if f.is_file() and f.name not in skip and f.suffix in ['.py','.C','.md','.csv','.json','.jsonl','.tsv','.log']:
        shutil.copy2(f,A/f.name)
for name in ['snapshot','results','timestamp','figures','logbook','raw4305']:
    shutil.copytree(P/name,A/name)
shutil.copy2(P/'CHECKPOINT_FINAL.md',A/'CHECKPOINT.md')
shutil.copy2(P/'CHECKPOINT.md',A/'CHECKPOINT_PROGRESS.md')
(A/'selection_records').mkdir()
for f in (P/'columns').iterdir():
    if f.name.endswith('_selection.tsv') or f.name.endswith('_columns.txt'):shutil.copy2(f,A/'selection_records'/f.name)
files=sorted(f for f in A.rglob('*') if f.is_file())
manifest=''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+str(f.relative_to(A))+'\n' for f in files)
(A/'SHA256SUMS').write_text(manifest)
print(json.dumps(dict(files=len(files),bytes=sum(f.stat().st_size for f in files),staging=str(A))))
