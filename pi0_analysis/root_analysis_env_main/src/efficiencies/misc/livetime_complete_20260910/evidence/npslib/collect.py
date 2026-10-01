"""Freeze NPSlib source/history and audit previously extracted raw config text."""
from pathlib import Path
import datetime,hashlib,json,os,re,shlex,subprocess
P=Path(__file__).resolve().parent;E=P/'evidence';E.mkdir(exist_ok=True)
B=Path('/group/nps/apps/NPSlib')
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
commands=[];sources=[]
def run(name,args):
 r=subprocess.run(args,capture_output=True,timeout=40)
 for suffix,data in [('stdout',r.stdout),('stderr',r.stderr)]: (E/(name+'.'+suffix+'.txt')).write_bytes(data)
 commands.append(dict(name=name,command=shlex.join(args),argv=args,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),exit=r.returncode))
 return r
for v in ['0269b41','bf8ec57','fddb8e9']:
 repo=B/v/'src';git=['git','-c','safe.directory='+str(repo),'-C',str(repo)]
 run(v+'_head',git+['log','-1','--format=%H %cI %s'])
 run(v+'_status',git+['status','--short','--','README.md','src'])
 run(v+'_search',['rg','-n','-i',r'edtm|pre40|pre100|pre150|pre200|livetime|deadtime|sis3801|\brol[12]\b|coin_sparse',str(repo/'src')])
 for name in ['README.md','src/THcNPSCoinTime.h','src/THcNPSCoinTime.cxx','src/THcNPSConfigEvtHandler.h','src/THcNPSConfigEvtHandler.cxx','src/VTPModule.cxx','src/VTPModule.h','src/THcNPSCalorimeter.cxx']:
  f=repo/name
  if not f.exists():continue
  b=f.read_bytes();out=E/(v+'_'+f.name);out.write_bytes(b)
  sources.append(dict(path=str(f),snapshot=str(out.relative_to(P)),sha256=hashlib.sha256(b).hexdigest(),bytes=len(b)))
repo=B/'fddb8e9/src';git=['git','-c','safe.directory='+str(repo),'-C',str(repo)]
run('edtm_blame',git+['blame','-L','78,78','--','src/THcNPSCoinTime.h'])
run('config_history',git+['log','--all','--format=%H %aI %s','--','src/THcNPSConfigEvtHandler.cxx'])
run('history_search',git+['log','--all','--oneline','--extended-regexp','--regexp-ignore-case','--grep=EDTM|ROC5|livetime|deadtime|prescale|trigger|config'])
rows=[];first={}
cfgdir=M/'livetime_refresh_20260909/raw4305/config_text'
for f in sorted(cfgdir.glob('record*_tag137.txt'),key=lambda f:int(re.search(r'record(\d+)',f.name)[1])):
 data=f.read_bytes();keys={}
 for lineno,line in enumerate(data.decode().splitlines(),1):
  a=line.split()
  if not a:continue
  val=' '.join(a[1:])
  if val=='end':continue
  keys.setdefault(a[0],[]).append(dict(line=lineno,value=val))
  # Model unique-key map emplace only, not a runtime replay of this handler.
  first.setdefault(a[0],val)
 rows.append(dict(file=str(f),sha256=hashlib.sha256(data).hexdigest(),keys={k:keys[k] for k in ['VTP_NPS_TRIG_LATENCY','VTP_NPS_TRIG_PRESCALE'] if k in keys}))
audit=dict(same_directory=os.path.samefile(B,'/u/group/nps/apps/NPSlib'),config_records=rows,first_key_model={k:first[k] for k in ['VTP_NPS_TRIG_LATENCY','VTP_NPS_TRIG_PRESCALE']},runtime_handler_execution=False)
(E/'config_key_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
(E/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
(E/'sources.json').write_text(json.dumps(sources,indent=2)+'\n')
print('same directory:',audit['same_directory'],'source snapshots:',len(sources),'config records:',len(rows))
print('VTP prescale rows per record:',[len(r['keys'].get('VTP_NPS_TRIG_PRESCALE',[])) for r in rows if r['keys']])
print('First-key model:',audit['first_key_model'])
