"""Read archived logs; verify TI settings without assigning interleaved ROC text."""
import csv, hashlib, json, re
from pathlib import Path
P=Path(__file__).resolve().parent
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
rows=[]; summary=[]
for run,total in [(4303,3310432),(4305,560647)]:
 f=M/'livetime_refresh_20260909/logbook'/f'run{run}_user_log.txt'
 s=f.read_text(); boards=list(re.finditer(r'boardID\s+\(0x0000\) = (0x[0-9a-fA-F]+)',s))
 assert len(boards)==10
 for i,m in enumerate(boards):
  part=s[m.end():boards[i+1].start() if i+1<len(boards) else len(s)]
  vals={k:re.findall(k+r'\s+\(0x[0-9a-f]+\) = (0x[0-9a-fA-F]+)',part) for k in ['dataFormat','vmeControl','livetime','busytime']}
  # Local register pair follows its board ID. Timers can be interleaved:
  # refuse ambiguous chunks instead of inventing a board assignment.
  assert len(vals['dataFormat'])==len(vals['vmeControl'])==1
  df=int(vals['dataFormat'][0],16);vc=int(vals['vmeControl'][0],16)
  assert (df>>24)==5 and not vc&(1<<22) and vc&(1<<23)
  timer_ok=len(vals['livetime'])==len(vals['busytime'])==1
  rows.append(dict(run=run,phase='initial' if i<5 else 'end',boardID=m[1],line=s[:m.start()].count('\n')+1,
   dataFormat=hex(df),vmeControl=hex(vc),broadcast_buffer=df>>24,use_local=bool(vc&(1<<22)),busy_on_buffer=bool(vc&(1<<23)),
   timers_unambiguous=timer_ok,live_hex=vals['livetime'][0] if timer_ok else '',busy_hex=vals['busytime'][0] if timer_ok else ''))
 counts={k:[int(x) for x in re.findall(k+r':\s*(\d+)',s)] for k in ['Readout Count','Ack Count','L1A Count']}
 for k,v in counts.items():assert len(v)==10 and v[-5:]==[total]*5,(run,k,v)
 summary.append(dict(run=run,log_path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),end_counts=counts,expected_recorded=total))
with (P/'ti_log_registers.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(P/'log_checks.json').write_text(json.dumps(summary,indent=2)+'\n')
print('PASS: 20 TI register pairs: broadcast5, local disabled, busy enabled; all 10 end TI count triplets match recorded totals.')
print('Two interleaved timer chunks left unassigned; no livetime correction calculated.')
