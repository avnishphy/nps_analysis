"""Run the phase diagnostic as each run's ordinary analysis becomes available."""
from pathlib import Path
import subprocess,json,time
P=Path(__file__).resolve().parent
while True:
    inv=json.loads((P/'inventory.json').read_text())['runs'];pending=[]
    for r in inv:
        if not r['files']:continue
        run=r['run'];dst=P/f'timestamp/run{run}.json';src=P/f'results/run{run}.json'
        if dst.exists():continue
        if not src.exists():pending.append(run);continue
        result=json.loads(src.read_text())
        if result['status']=='CALCULATED':subprocess.run(['python3','timestamp_probe.py',str(run)],cwd=P,check=True)
        else:pending.append(run)
    if not pending:break
    time.sleep(5)
