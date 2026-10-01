"""Finish refreshed diagnostics after every cached run has completed."""
from pathlib import Path
import json,subprocess,time,concurrent.futures
P=Path(__file__).resolve().parent
while True:
    inv=json.loads((P/'inventory.json').read_text())['runs']
    progress=json.loads((P/'analysis_progress.json').read_text())
    if not progress['pending'] and all((P/f'timestamp/run{r["run"]}.json').exists() for r in inv if r['files']):break
    time.sleep(5)
def execute(script):
    print('START',script,flush=True)
    with (P/(Path(script).stem+'.log')).open('w') as f:subprocess.run(['python3',script],cwd=P,stdout=f,stderr=subprocess.STDOUT,check=True)
    print('DONE',script,flush=True)
execute('catalog_coverage.py');execute('audit.py')
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    for _ in pool.map(execute,['acceptance_tests.py','clock_scaler_checks.py','event_clock_checks.py','phase_timing_bridge.py']):pass
execute('make_refresh_figures.py')
(P/'finalize.done').write_text('All final diagnostics and focused figures completed.\n')
