import math
import csv
import subprocess
import sys
from array import array
import ROOT
from pathlib import Path

out=Path(sys.argv[1])
out.mkdir(exist_ok=True)
builder=Path(__file__).resolve().parent/'build_xsec_test_binary.py'
def build_extractor(nphi):
    return Path(subprocess.check_output([
        sys.executable, str(builder), '--out-dir', str(out/f'build_{nphi}'),
        '--phi-bins', str(nphi), '--t-edges=-0.6,-0.375,-0.15',
        '--q2-edges=3,5', '--xb-edges=.25,.455'
    ], text=True).strip())
extractors={6:build_extractor(6), 12:build_extractor(12)}
plus=[1.60077,-0.01523,37.08142,-4.11060,23.26192,.00983,.87073,-5.77115,-271.08678,.13766,-.00855,.27885,-1.13212,-1.50415,-6.34766,.55769,-.01709]
minus=[1.75169,.11144,47.35877,-4.69434,1.60552,.008,.44194,-2.29188,-41.67194,.69475,.02527,-.50178,-1.22825,-1.16878,5.75825,-1.00355,.05055]
truth=plus.copy()
truth[4]*=.88; truth[6]*=1.20; truth[8]*=.82
print('truth plus.p5,p7,p9',truth[4],truth[6],truth[8])
mp=.9382720813; mpi=.1349768; q2=4.; w=2.8
xb=q2/(w*w-mp*mp+q2)
nu=(w*w-mp*mp+q2)/(2*mp); ebeam=10.538
y=nu/ebeam; z=q2/(4*ebeam*ebeam)
eps=(1-y-z)/(1-y+.5*y*y+z)
def sigma(t,phi,p):
    ws=w*w; eg=(ws-mp*mp-q2)/(2*w); pg=math.sqrt(eg*eg+q2)
    epi=(ws+mpi*mpi-mp*mp)/(2*w); pp=math.sqrt(epi*epi-mpi*mpi)
    ct=(t+q2-mpi*mpi+2*eg*epi)/(2*pg*pp)
    st=math.sqrt(1-ct*ct); norm=8.539/(ws-.938*.938)**2/2/3.1415928/1e6
    v=0
    for a in (p,minus):
        T=a[4]/q2*math.exp(a[5]*q2*q2)/(ws**a[11]+w**a[15])*math.exp(a[13]*abs(t))
        LT=a[6]/(1+a[9]*q2)*math.exp(a[7]*abs(t))*st/ws**a[12]
        TT=a[8]/(1+q2)*math.exp(-7*abs(t))*st*st
        v+=.5*norm*(T+math.sqrt(2*eps*(1+eps))*LT*math.cos(phi)+eps*TT*math.cos(2*phi))
    return v
def tree(path,name,fields):
    f=ROOT.TFile(str(path),'RECREATE'); t=ROOT.TTree(name,name)
    buf={}
    for key,typ in fields.items():
        code={'F':'f','D':'d','I':'i'}[typ]
        buf[key]=array(code,[0])
        t.Branch(key,buf[key],f'{key}/{typ}')
    return f,t,buf
df,dt,d=tree(out/'data.root','physics',{
 'Q2':'D','t':'D','tmin':'D','xB':'D','phi':'D','mmiss_all':'D',
 'pi0_weight':'D','scale':'F','charge_uC':'F','run_number':'I','W':'D'})
sf,st,m=tree(out/'sim.root','simulation',{
 'Q2':'F','t':'F','tmin':'F','xB':'F','phi':'F','mmiss':'F',
 'full_weight':'F','sigcm':'F','is_exclusive':'I','W':'F',
 'Q2i':'F','Wi':'F','ti':'F','phipqi':'F','epsilon_i':'F'})
def fill(buf,t,values):
    for k,v in values.items(): buf[k][0]=v
    t.Fill()
for it,tc in enumerate([-.5,-.25]):
    for ip in range(6):
        expected=0.
        for j in range(50):
            # Vary event kinematics inside each reconstructed bin.
            t=tc+.012*math.sin(j*.8)
            phi=(ip+.5)*math.pi/3+.08*math.sin(j*.6)
            base=1.e7*(1+.08*math.cos(j*.4))
            original=sigma(t,phi,plus)
            expected+=base*sigma(t,phi,truth)
            fill(m,st,dict(Q2=q2,t=t,tmin=-.1,xB=xb,phi=phi,mmiss=.9,
                full_weight=base*original,sigcm=original,is_exclusive=1,W=w,
                Q2i=q2,Wi=w,ti=-t,phipqi=phi,
                epsilon_i=(-999. if it==0 and ip==0 and j==0 else eps)))
        for j in range(200):
            fill(d,dt,dict(Q2=q2,t=tc, tmin=-.1,xB=xb,
                phi=(ip+.5)*math.pi/3,mmiss_all=.9,
                pi0_weight=expected/200,scale=1,charge_uC=1,run_number=1,W=w))
# Invalid denominator and invalid full weight: selected geometry, rejected before division.
for sig,weight in [(0.,1.),(float('nan'),1.),(1.e-8,float('nan'))]:
    fill(m,st,dict(Q2=q2,t=-.5,tmin=-.1,xB=xb,phi=.5,mmiss=.9,
        full_weight=weight,sigcm=sig,is_exclusive=1,W=w,
        Q2i=q2,Wi=w,ti=.5,phipqi=.5,epsilon_i=eps))
df.cd(); dt.Write(); df.Close()
sf.cd(); st.Write(); sf.Close()


def run(args, expect_success=True):
    p=subprocess.run([str(extractors[int(args[1])]), '--data-file',str(out/'data.root'),
                      '--sim-file',str(out/'sim.root'), '--out-dir',str(out/args[0]),
                      '--target-contam','1','--target-contam-err','0',
                      '--no-diagnostics','--no-png','--no-pdf']+args[2:],
                     text=True,capture_output=True)
    if (p.returncode==0)!=expect_success:
        raise RuntimeError(f'Unexpected extractor status {p.returncode}: {p.stderr} {p.stdout}')
    return p
p=run(['fit','6'])
assert '[MODEL_FIT] status=converged' in p.stdout, p.stdout
assert 'bad_full_weight=1 bad_sigcm=2' in p.stdout
assert 'epsilon_fallback=1' in p.stdout
rows=list(csv.DictReader(open(out/'fit/excl_xsec_pi0_analysis_simc_model_summary.csv')))
assert len(rows)==12
rf=ROOT.TFile.Open(str(out/'fit/excl_xsec_pi0_analysis_simc_model_output.root'))
best=plus.copy()
for item in rf.Get('model_parameters'):
    name=str(item.name)
    if name.startswith('plus.p'):
        best[int(name.split('p')[-1])-1]=item.value
for j in (4,6,8):
    assert abs(best[j]/truth[j]-1)<.02,(j,best[j],truth[j])
assert rf.Get('model_fit_covariance').GetNrows()==3
assert rf.Get('n_sim_bad_sigcm').GetVal()==2
assert rf.Get('n_sim_bad_full_weight').GetVal()==1
assert rf.Get('n_sim_cached').GetVal()==600
assert rf.Get('n_sim_epsilon_fallback').GetVal()==1
sf=ROOT.TFile.Open(str(out/'sim.root')); sim=sf.Get('simulation')
max_var_rel=max_point_rel=max_q2_rel=max_pull=max_xsec_rel=0.
for row in rows:
    it=int(row['it']); lo=float(row['t_lo']); hi=float(row['t_hi'])
    plo=float(row['phi_lo']); phi_hi=float(row['phi_hi'])
    weights=[]
    for ev in sim:
        if not (math.isfinite(ev.sigcm) and ev.sigcm>0 and math.isfinite(ev.full_weight)): continue
        if not (lo<=ev.t<(hi if it==0 else hi+1e-9) and plo<=ev.phi<phi_hi): continue
        weights.append(ev.full_weight/ev.sigcm*sigma(-ev.ti,ev.phipqi,best))
    v=sum(x*x for x in weights)
    max_var_rel=max(max_var_rel,abs(float(row['sim_err'])**2/v-1))
    assert abs(float(row['sim'])/sum(weights)-1)<2e-6
    center=sigma(float(row['t_center']),float(row['phi_center']),best)
    max_point_rel=max(max_point_rel,abs(float(row['model_xsec_phi_center'])/center-1))
    truth_center=sigma(float(row['t_center']),float(row['phi_center']),truth)
    max_xsec_rel=max(max_xsec_rel,abs(float(row['xsec'])/truth_center-1))
    wref=float(row['W_ref']); xbref=float(row['xb_ref'])
    q2ref=xbref*(wref*wref-mp*mp)/(1-xbref)
    max_q2_rel=max(max_q2_rel,abs(float(row['q2_ref'])/q2ref-1))
    assert abs(float(row['t_center'])-
               .5*(float(row['t_lo'])+float(row['t_hi'])))<1e-8
    assert abs(float(row['phi_center'])-
               .5*(float(row['phi_lo'])+float(row['phi_hi'])))<1e-8
    max_pull=max(max_pull,abs(float(row['ratio'])-1)/float(row['ratio_err']))
    assert row['bin_status']=='measured'
assert max_var_rel<1e-5 and max_point_rel<1e-5 and max_q2_rel<1e-7
assert max_pull<.5 and max_xsec_rel<.03
bad=run(['bad_norm','6','--normalize_mmiss'],False)
assert 'conflicts with iterative absolute' in bad.stderr
failed=run(['failed','6','--model-max-evaluations','1'],False)
assert 'fit failed' in failed.stderr
inactive=run(['inactive','6','--model-free','plus.p1'],False)
assert 'inactive when pi0 fpifact=0' in inactive.stderr
fixed=run(['empty','12','--fixed-default-model'])
empty=list(csv.DictReader(open(out/'empty/excl_xsec_pi0_analysis_simc_model_summary.csv')))
assert any(r['bin_status']!='measured' and math.isnan(float(r['xsec']))
           for r in empty)
print(f'PASS closure: fitted p5={best[4]:.6g} p7={best[6]:.6g} p9={best[8]:.6g}')
print(f'PASS max ratio pull={max_pull:.3g}, xsec relative error={max_xsec_rel:.3g}')
print(f'PASS MC sumw2={max_var_rel:.3g}, model center={max_point_rel:.3g}, Q2_ref={max_q2_rel:.3g}')
print('PASS invalid inputs, normalization conflict, convergence failure, empty bins')
