"""Bounded, read-only source discovery; output is an auditable local record."""
import datetime,hashlib,json,shlex,subprocess
from pathlib import Path
P=Path(__file__).resolve().parent
E=P/'evidence';E.mkdir(exist_ok=True)
A='/u/group/nps/apps'
M=Path('/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/efficiencies/misc')
records=[]
def run(name,args):
 result=subprocess.run(args,capture_output=True,timeout=50)
 for key,data in [('stdout',result.stdout),('stderr',result.stderr)]:
  (E/(name+'.'+key+'.txt')).write_bytes(data)
 records.append(dict(name=name,argv=args,command=shlex.join(args),utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),exit_code=result.returncode,stdout_bytes=len(result.stdout),stderr_bytes=len(result.stderr)))
 print(name,result.returncode,len(result.stdout),len(result.stderr))
 return result
run('apps_listing',['ls','-la',A])
inventory=run('apps_files',['rg','--files','--hidden','--no-ignore','-g','!**/.git/**',A])
files=inventory.stdout.decode().splitlines()
run('apps_roc_text',['rg','-n','-i','--hidden','--no-ignore','-g','!**/.git/**','-g','!**/BUILD/**','-g','!*.so*','-g','!*.o','-g','!*.a','-g','!*.pcm','-g','!*.root','coin_sparse|COOL_HOME|cdaqfs|\brol[12]\b|sis3801|tiLib|/site/coda',A])
run('apps_roc_text_precise',['rg','-n','-i','--hidden','--no-ignore','-g','!**/.git/**','-g','!**/BUILD*/**','-g','!**/build/**','-g','!*.so*','-g','!*.o','-g','!*.a','-g','!*.pcm','-g','!*.root',r'coin_sparse|COOL_HOME|cdaqfs|\brol[12]\b|sis3801|\btiLib\b|/site/coda',A])
run('apps_decoder_positive_control',['rg','-n','Scaler3801|AnalyzeBuffer|scal_overflows',A+'/hcana/1.0.1/src/src/THcScalerEvtHandler.cxx'])
run('modules_paths',['rg','-n','cdaqfs|module load hcana|set topdir|REPLAYDIR','/u/group/nps/modulefiles/nps_replay'])
run('online_locations',['ls','-ld','/cdaqfs1','/u/cdaqfs1','/home/coda','/u/home/coda','/adaqfs'])
for version in ['1.0.0','1.0.1','1.0.2']:
 repo=A+'/hcana/'+version+'/src'
 git=['git','-c','safe.directory='+repo,'-C',repo]
 run('hcana_'+version+'_head',git+['log','-1','--format=%H %cI %s'])
 run('hcana_'+version+'_scaler_status',git+['status','--short','--','src/THcScalerEvtHandler.cxx'])
for runnum in [4303,4305]:
 run('run'+str(runnum)+'_roc_identity',['rg','-n','channel ROC5|channel npsvme5|ROC5 INFO: ended after|ROC1 INFO: ended after|ROC3 INFO: ended after',str(M/'livetime_refresh_20260909/logbook'/f'run{runnum}_user_log.txt')])
sources=[Path(A)/'setup.sh',Path('/u/group/nps/modulefiles/nps_replay/11.11.23'),Path('/u/group/nps/modulefiles/nps_replay/2.20.24')]
sources += [Path(A)/'hcana'/v/'src/src/THcScalerEvtHandler.cxx' for v in ['1.0.0','1.0.1','1.0.2']]
source_rows=[]
for i,f in enumerate(sources):
 b=f.read_bytes();rel='evidence/source_'+str(i)+'_'+f.name;(P/rel).write_bytes(b)
 source_rows.append(dict(path=str(f),resolved_path=str(f.resolve()),snapshot=rel,bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
(E/'source_manifest.json').write_text(json.dumps(source_rows,indent=2)+'\n')
(E/'commands.json').write_text(json.dumps(records,indent=2)+'\n')
(E/'inventory_summary.json').write_text(json.dumps(dict(files_listed=len(files),follow_directory_symlinks=False,hidden_included=True,git_internal_excluded=True),indent=2)+'\n')
print('inventory',len(files),'files; source hashes:',[(r['path'].split('/')[-4],r['sha256'][:12]) for r in source_rows[-3:]])
