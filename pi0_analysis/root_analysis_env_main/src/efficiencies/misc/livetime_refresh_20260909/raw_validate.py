"""Validate the one authorized raw file against ROOT and audit scaler banks.

Unassigned channels are not adopted as pulser denominators. Matching counts
do not establish cabling. No counter repair is performed.
"""
from pathlib import Path
import json,csv,re,numpy as np
from timestamp_probe import read
P=Path(__file__).resolve().parent;D=P/'raw4305'
e=np.fromfile(D/'physics.bin',dtype=np.uint64).reshape(-1,4);r=read(4305,0,'T');s=read(4305,0,'TSH')
physics={k:bool(np.array_equal(e[:,j],r[k])) for j,k in [(1,'g.evnum'),(2,'g.evtime'),(3,'g.trigbits')]};assert all(physics.values()) and len(e)==560647
headers=[0x620,0x200720,0x400820,0x600920,0x800a20,0xa00b20,0xc00c20]
rows=[];data=[]
for row in map(json.loads,(D/'scaler_banks.jsonl').read_text().splitlines()):
    v=np.fromfile(D/row['file'],dtype=np.uint32);a=[]
    for slot,h in zip(range(6,13),headers):
        pos=2+33*(slot-6);assert v[pos]==h
        a.append(v[pos+1:pos+33])
    data.append(a);rows.append(row)
a=np.array(data,dtype=np.int64);b=np.array([x['last_physics'] for x in rows]);assert np.array_equal(b,s['evNumber'][:-1]) and s['evNumber'][-1]==b[-1]+1
mapping={'H.1MHz.scalerTime':(12,31,1e-6),'H.EDTM.scaler':(12,14,1),'H.hL1ACCP.scaler':(10,16,1)}
for trig in range(1,7):mapping[f'H.hTRIG{trig}.scaler']=(11,9+trig,1)
checks=[]
for k,(slot,ch,factor) in mapping.items():
    x=a[:,slot-6,ch]*factor;ok=bool(np.allclose(x,s[k][:-1],rtol=0,atol=1e-9));assert ok,k
    assert abs(s[k][-1]-x[-1])<1e-9
    checks.append(dict(branch=k,slot=slot,channel=ch,all_raw_snapshots_match_ROOT=ok,last_ROOT_row_duplicates_final_raw_counter_with_event_boundary_incremented=True))
diff=a[1:]-a[:-1];i=int(np.flatnonzero((b[:-1]==105080)&(b[1:]==106127))[0]);lead=a[:,1,31]-a[:,6,14]
np.savez_compressed(D/'scaler_modules.npz',event=b,values=a)
with (D/'selected_raw_counters.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['raw_bank','event','clock_slot12_ch31','EDTM_slot12_ch14','L1_slot10_ch16','TRIG4_slot11_ch13','unassigned_slot7_ch31','unassigned_slot9_ch31'])
    for j,x in enumerate(a):w.writerow([j+1,b[j],x[6,31],x[6,14],x[4,16],x[5,13],x[1,31],x[3,31]])
changes=np.flatnonzero(np.diff(lead)!=0)
result=dict(raw_physics_events=len(e),physics_branches_identical=physics,raw_scaler_banks=len(a),ROOT_scaler_rows=len(s['evNumber']),scaler_checks=checks,anomaly=dict(raw_banks=[i+1,i+2],raw_records=[rows[i]['record'],rows[i+1]['record']],event_boundaries=[int(b[i]),int(b[i+1])],per_slot_channel_increments={str(slot):diff[i,slot-6].tolist() for slot in range(6,13)}),unassigned_channels=dict(slot7_ch31_equals_slot9_ch31=bool(np.array_equal(a[:,1,31],a[:,3,31])),difference_from_EDTM_values_counts={str(k):int(v) for k,v in zip(*np.unique(lead,return_counts=True))},difference_steps=[dict(before_event=int(b[k]),after_event=int(b[k+1]),before_difference=int(lead[k]),after_difference=int(lead[k+1])) for k in changes],interpretation='Pulse-like count agreement is diagnostic only; map labels these channels Empty, not EDTM. No replacement denominator or repair adopted.'))
(D/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
