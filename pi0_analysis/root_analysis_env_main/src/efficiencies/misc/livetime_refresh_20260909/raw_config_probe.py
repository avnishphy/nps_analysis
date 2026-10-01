"""Preserve raw config text and search named settings; no inferred routing."""
from pathlib import Path
import json,re,collections,numpy as np,shutil,hashlib
P=Path(__file__).resolve().parent;D=P/'raw4305';O=D/'config_text';O.mkdir(exist_ok=True)
matches=collections.defaultdict(set);configs=[]
for row in map(json.loads,(D/'index.jsonl').read_text().splitlines()):
    if row['tag'] not in [137,182,183,184,185]:continue
    v=np.fromfile(D/row['file'],dtype=np.uint32)
    assert ((v[1]>>8)&63)==16 and ((v[3]>>8)&63)==3
    text=v[4:].tobytes().decode('latin1').split('\x00')[0]
    if row['tag']==137:
        dst=O/(Path(row['file']).stem+'.txt');dst.write_text(text)
        configs.append(dict(record=row['record'],ROC_bank_tag=int(v[3]>>16),text_file=str(dst.relative_to(D))))
    for line in text.splitlines():
        if re.search(r'EDTM|PRE40|PRE100|PRE150|PRE200|prescale|\bps[1-6]\b|\bTI_|\bTS_|TRG|TRIG|BLOCK',line,re.I):matches[(row['tag'],line)].add(row['last_physics'])
(D/'config_search.json').write_text(json.dumps(dict(configs=configs,matches=[dict(tag=k[0],text=k[1],distinct_event_snapshots=len(v)) for k,v in matches.items()]),indent=2)+'\n')
S=D/'source_evidence';S.mkdir(exist_ok=True);sources=[]
for src in [Path('/u/group/halla/apps/analyzer/1.7.12/src/hana_decode')/n for n in ['CodaDecoder.cxx','THaCodaFile.cxx','GenScaler.cxx','Scaler3801.cxx']]+[Path('/group/nps/singhav/nps_replay/MAPS/db_cratemap.dat'),Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc/livetime_4398_check_20260909/replay_snapshot/MAPS_db_HScalevt.dat')]:
    dst=S/src.name;shutil.copy2(src,dst);sources.append(dict(source=str(src),file=str(dst.relative_to(D)),sha256=hashlib.sha256(dst.read_bytes()).hexdigest()))
(D/'source_manifest.json').write_text(json.dumps(sources,indent=2)+'\n')
print('Saved',len(configs),'FADC/VTP configuration text banks. No master prescale readback inferred from VTP settings.')
