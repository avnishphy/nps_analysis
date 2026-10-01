"""Stop only this directory's export controller/readers, retain completed files."""
from pathlib import Path
import os,signal,json,time
P=Path(__file__).resolve().parent;targets=[]
for d in Path('/proc').glob('[0-9]*'):
    try:
        if Path(os.readlink(d/'cwd')).resolve()!=P:continue
        args=(d/'cmdline').read_bytes().decode().split('\0');cmd=' '.join(args)
        controller=any(Path(a).name=='run_exports.py' for a in args)
        reader=any('export_cached.C(' in a for a in args) and Path(args[0]).name in ['hcana','hcana.exe','csh','tcsh']
        if controller or reader:targets.append(dict(pid=int(d.name),command=cmd))
    except (FileNotFoundError,PermissionError,ProcessLookupError):pass
(P/'restart_targets.json').write_text(json.dumps(targets,indent=2)+'\n')
for t in targets:
    try:os.kill(t['pid'],signal.SIGTERM)
    except ProcessLookupError:pass
print('Stopped only diagnostic export processes:',[t['pid'] for t in targets])
