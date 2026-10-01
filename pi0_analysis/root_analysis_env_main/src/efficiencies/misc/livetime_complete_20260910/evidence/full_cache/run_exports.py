"""Cache-only exports: all other runs first, refreshed 4397/4399/4402 last."""
from pathlib import Path
import json,subprocess,datetime,time
P=Path(__file__).resolve().parent;TAIL={4397,4399,4402}
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def note(state,**kw):
    row=dict(time_utc=now(),state=state,**kw)
    with (P/'dispatch.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    print(json.dumps(row),flush=True)
def manifest(rows):
    (P/'manifest.tsv').write_text(''.join(f"{r['run']}\t{f['segment']}\t{f['path']}\n" for r in rows for f in r['files']))
def spawn(arg,name):
    log=(P/name).open('a')
    cmd=f'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; hcana -l -b -q "export_cached.C({arg})"'
    p=subprocess.Popen(['csh','-c',cmd],cwd=P,stdout=log,stderr=subprocess.STDOUT)
    note('worker_started',pid=p.pid,argument=arg,log=name);return p,log
def complete(workers):
    for p,log in workers:
        status=p.wait();log.close()
        if status:raise RuntimeError(f'Export worker {p.pid} exited {status}')
def check_done(rows):
    missing=[(r['run'],f['segment']) for r in rows for f in r['files'] if not (P/f"columns/run{r['run']}_seg{f['segment']}.done").exists()]
    if missing:raise RuntimeError('Missing export completions: '+str(missing))
def scan(run):
    found={}
    for source in ['updated','production']:
        fs=Path('/cache/hallc/c-nps/analysis/pass2/replays',source).glob(f'nps_hms_coin_{run}_*_1_-1.root')
        found[source]=sorted([dict(path=str(f),segment=int(f.name.split('_')[4]),size=f.stat().st_size,mtime=f.stat().st_mtime) for f in fs if f.is_file()],key=lambda f:f['segment'])
    return found
inv=json.loads((P/'inventory.json').read_text());(P/'inventory_initial.json').write_text(json.dumps(inv,indent=2)+'\n')
early=[r for r in inv['runs'] if r['run'] not in TAIL]
manifest(early);(P/'manifest_early.tsv').write_text((P/'manifest.tsv').read_text())
note('early_started',runs=[r['run'] for r in early],segments=sum(len(r['files']) for r in early))
complete([spawn(f'0,{b},3',f'export_batch{b}.log') for b in range(3)])
check_done(early);note('early_finished')
# Freeze deferred runs only after the early exports have all finished.
for r in inv['runs']:
    if r['run'] not in TAIL:continue
    found=scan(r['run']);time.sleep(10)
    stable=scan(r['run'])
    if stable!=found:raise RuntimeError(f"Deferred run {r['run']} files still changing; retry tail after staging stabilizes")
    source='updated' if found['updated'] else 'production' if found['production'] else 'none'
    r.update(source=source,files=found.get(source,[]),other_cached_sources=found,dispatch_snapshot_utc=now())
    identity=lambda s:(Path(s).parent.name,Path(s).name)
    r['same_saved_files']={identity(s) for s in r['saved_files']}=={identity(f['path']) for f in r['files']}
    note('tail_inventory',run=r['run'],source=source,segments=[f['segment'] for f in r['files']])
inv['deferred_refresh_utc']=now();inv['deferred_runs']=sorted(TAIL)
(P/'inventory.json').write_text(json.dumps(inv,indent=2)+'\n')
manifest(inv['runs']);tail=[r for r in inv['runs'] if r['run'] in TAIL]
note('tail_started')
complete([spawn(str(r),f'export_batch_tail{r}.log') for r in sorted(TAIL)])
check_done(inv['runs']);note('exports_complete',segments=sum(len(r['files']) for r in inv['runs']))
