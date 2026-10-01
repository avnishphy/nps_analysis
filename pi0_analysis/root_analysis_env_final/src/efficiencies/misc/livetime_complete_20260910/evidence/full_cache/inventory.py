"""Read cache directory entries only; never request tape staging."""
from pathlib import Path
import csv,json,hashlib,shutil,re,datetime
P=Path(__file__).resolve().parent
MAIN=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main')
CACHE=Path('/cache/hallc/c-nps/analysis/pass2/replays')
def rows(p):
    return [{k.strip():(v or '').strip() for k,v in r.items() if k is not None} for r in csv.DictReader(p.open())]
meta=[r for r in rows(MAIN/'config/nps_dvcs_all_kins_main.csv') if r['Kin_old']=='KinC_x60_4b' and r['target'].lower()=='lh2']
meta.sort(key=lambda r:int(r['run_number']))
(P/'snapshot').mkdir(exist_ok=True)
sources=['config/nps_dvcs_all_kins_main.csv','output/efficiency_stuff/efficiency_KinC_x60_4b.csv','output/efficiency_stuff/selection_report_KinC_x60_4b.csv','output/efficiency_stuff/plots/efficiency_multipanel_vs_run_KinC_x60_4b_lh2.png']
sources += [str(p.relative_to(MAIN)) for p in (MAIN/'src/efficiencies').glob('*.h')]
sources += ['src/efficiencies/compute_efficiencies_stuff.cxx']
prov=[]
for rel in sources:
    src=MAIN/rel;dst=P/'snapshot'/src.name
    shutil.copy2(src,dst)
    prov.append(dict(source=str(src),snapshot=str(dst.relative_to(P)),sha256=hashlib.sha256(dst.read_bytes()).hexdigest(),mtime=src.stat().st_mtime))
(P/'source_provenance.json').write_text(json.dumps(prov,indent=2))
base={int(r['run_number']):r for r in rows(P/'snapshot/efficiency_KinC_x60_4b.csv')}
savedsel=rows(P/'snapshot/selection_report_KinC_x60_4b.csv')
inventory=[];manifest=[]
for r in meta:
    run=int(r['run_number']);found={}
    for source in ['updated','production']:
        fs=list((CACHE/source).glob(f'nps_hms_coin_{run}_*_1_-1.root'))
        found[source]=sorted([dict(path=str(p),segment=int(p.name.split('_')[4]),size=p.stat().st_size,mtime=p.stat().st_mtime) for p in fs if p.is_file()],key=lambda x:x['segment'])
    # Same source preference as current analysis; never combine replay versions.
    source='updated' if found['updated'] else 'production' if found['production'] else 'none'
    files=found.get(source,[])
    oldfiles=[x['segment_file'] for x in savedsel if int(x['run_number'])==run]
    def identity(s):
        p=Path(s);return (p.parent.name,p.name)
    q=dict(run=run,run_type=r['Type'],prescale_token=r['prescale'],metadata=r,source=source,files=files,other_cached_sources=found,saved=base.get(run),saved_files=oldfiles,same_saved_files={identity(s) for s in oldfiles}=={identity(f['path']) for f in files})
    inventory.append(q)
    for f in files:manifest.append(f"{run}\t{f['segment']}\t{f['path']}")
(P/'inventory.json').write_text(json.dumps(dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),policy='cache only; updated preferred per run; no staging commands; do not bridge missing segments',runs=inventory),indent=2))
(P/'manifest.tsv').write_text('\n'.join(manifest)+'\n')
print('Runs',len(meta),'files',len(manifest),'cached runs',sum(bool(q['files']) for q in inventory),'same saved coverage',sum(q['same_saved_files'] for q in inventory))
for q in inventory:print(q['run'],q['run_type'],q['prescale_token'],q['source'],[f['segment'] for f in q['files']],'saved',len(q['saved_files']),'same',q['same_saved_files'])
