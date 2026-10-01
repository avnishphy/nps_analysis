"""Test, rather than assume, cancellation of EDTM timing tails against trigger time."""
from pathlib import Path
import json,numpy as np
from timestamp_probe import read
P=Path(__file__).resolve().parent;out=[]
for run in [4253,4255,4259,4301,4305,4350]:
 r=json.loads((P/f'results/run{run}.json').read_text());values=[];tails=[]
 for d in r['segment_details']:
  e=read(run,d['segment'],'T');raw=e['T.hms.hEDTM_tdcTimeRaw'];tc=e['T.hms.hEDTM_tdcTime'];tr=e[f'T.hms.hTRIG{r["trigger"]}_tdcTime']
  good=(raw>1)&(abs(raw-r['raw_peak_channels'])<=500)&(e[f'T.hms.hTRIG{r["trigger"]}_tdcTimeRaw']>1)
  core=good&(abs(tc-r['corrected_peak_ns'])<=2);tail=good&~core
  values.extend((tc-tr)[core].tolist());tails.extend((tc-tr)[tail].tolist())
 center=float(np.median(values));out.append(dict(run=run,selection='whole available files, raw-window candidates with selected-trigger raw TDC >1',relative_center_ns=center,core_quantiles_ns=np.quantile(values,[.001,.5,.999]).tolist(),absolute_time_tail_candidates=len(tails),tails_recovered_with_relative_2ns=int((abs(np.array(tails)-center)<=2).sum()),tail_relative_times_ns=tails))
(P/'relative_timing_check.json').write_text(json.dumps(out,indent=2)+'\n')
print([(r['run'],r['absolute_time_tail_candidates'],r['tails_recovered_with_relative_2ns']) for r in out])
