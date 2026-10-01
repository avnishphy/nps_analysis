"""Analyze each run once its frozen scalar exports are complete."""
from pathlib import Path
import subprocess,json,time,datetime
P=Path(__file__).resolve().parent;finished=set()
while True:
    inv=json.loads((P/'inventory.json').read_text())['runs']
    for r in inv:
        run=r['run']
        if run in finished:continue
        if all((P/f'columns/run{run}_seg{f["segment"]}.done').exists() for f in r['files']):
            subprocess.run(['python3','analyze.py',str(run)],cwd=P,check=True)
            result=json.loads((P/f'results/run{run}.json').read_text())
            if result['status'] not in ['CALCULATED','NO_CACHED_FILES','NO_VALID_SELECTION']:raise RuntimeError(str(result))
            finished.add(run)
    (P/'analysis_progress.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'finished':sorted(finished),'pending':[r['run'] for r in inv if r['run'] not in finished]},indent=2)+'\n')
    if len(finished)==len(inv):break
    time.sleep(5)
