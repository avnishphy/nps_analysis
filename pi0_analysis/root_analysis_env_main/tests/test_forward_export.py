#!/usr/bin/env python3
"""Synthetic ROOT integration test; source the Hall C ROOT environment first."""
from array import array
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

import ROOT

ROOT.gROOT.SetBatch(True)
REPO = Path(__file__).resolve().parents[1]
SRC = REPO / 'src/xsec_extract'


def make_tree(path, name, schema, entries):
    output = ROOT.TFile(str(path), 'RECREATE')
    tree = ROOT.TTree(name, name)
    buffers = {}
    types = {'D': 'd', 'F': 'f', 'I': 'i', 'l': 'Q'}
    for key, kind in schema.items():
        buffers[key] = array(types[kind], [0])
        tree.Branch(key, buffers[key], f'{key}/{kind}')
    for event in entries:
        for key in schema:
            buffers[key][0] = event[key]
        tree.Fill()
    tree.Write()
    output.Close()


def main():
    out = Path(tempfile.mkdtemp(prefix='pi0_forward_export_test_'))
    config = json.loads((SRC/'xsec_config/xsec_config_x36_4.json').read_text())
    config['diamond_xb_q2_vertices'] = None
    (out/'config.json').write_text(json.dumps(config))
    subprocess.run([sys.executable, str(SRC/'generate_xsec_config.py'),
                    str(out/'config.json'), str(out/'xsec_config.h')], check=True)
    flags = subprocess.check_output(['root-config','--cflags','--libs'], text=True).split()
    subprocess.run(['g++','-std=c++17','-O2',str(SRC/'excl_xsec_pi0_analysis_no_simc_model.C'),
                    '-I'+str(out), '-I'+str(SRC), *flags, '-lMinuit2','-o',str(out/'extractor')], check=True)
    mp, mpi, q2, xb = .9382720813, .1349768, 4., .36
    W = math.sqrt(mp*mp+q2*(1/xb-1))
    q0=(W*W-mp*mp-q2)/(2*W); epi=(W*W+mpi*mpi-mp*mp)/(2*W)
    tmin=mpi*mpi-q2-2*(q0*epi-math.sqrt(q0*q0+q2)*math.sqrt(epi*epi-mpi*mpi))
    data=[]; mc=[]; truth=[]
    for i, tp in enumerate((-.8,-.3,-.1,.2,-.2)):
        mass = .95 if i < 4 else 1.5
        data.append(dict(Q2=q2,t=tp-.1,tmin=-.1,xB=xb,phi=.5*i,pi0_weight=.5,
                         scale=2.,charge_uC=100.,run_number=6415,mmiss_all=mass,W=W))
        mc.append(dict(Q2=q2,t=tp-.1,tmin=-.1,xB=xb,phi=.5*i,full_weight=1e-7,
                       is_exclusive=1,mmiss=mass,W=W,sigcm=1e-8,event_id=i))
        truth.append(dict(Q2i=q2,Wi=W,ti=-(min(tp,-.05)+tmin),phipqi=.5*i,
                          sigcm=1e-8,hsxptari=0.,hsyptari=0.))
    ds={k:'D' for k in data[0]};ds.update(scale='F',charge_uC='F',run_number='I')
    ms={k:'F' for k in mc[0]};ms.update(is_exclusive='I',event_id='l')
    make_tree(out/'data.root','physics',ds,data)
    make_tree(out/'mc.root','simulation',ms,mc)
    make_tree(out/'truth.root','h10',{k:'F' for k in truth[0]},truth)
    command=[str(out/'extractor'),'--data-file',str(out/'data.root'),'--sim-file',str(out/'mc.root'),
             '--vertex_simc_file',str(out/'truth.root'),'--out-dir',str(out/'cache'),'--prepare-forward-inputs']
    result=subprocess.run(command,text=True,capture_output=True)
    if result.returncode: raise AssertionError(result.stdout+result.stderr)
    with (out/'cache/data_events.csv').open() as f: de=list(csv.DictReader(f))
    with (out/'cache/mc_events.csv').open() as f: me=list(csv.DictReader(f))
    manifest=json.loads((out/'cache/forward_cache_manifest.json').read_text())
    assert len(de)==len(me)==4 and manifest['complete']
    assert float(de[0]['tprime']) < -.75 and float(de[3]['tprime']) > 0
    for row in de: assert math.isclose(float(row['weight']),1/.584,rel_tol=1e-7)
    for i,row in enumerate(me):
        assert int(row['event_id'])==i
        assert math.isclose(float(row['base_weight']),10.,rel_tol=1e-6)
        assert math.isclose(float(row['truth_tprime']),min((-.8,-.3,-.1,.2)[i],-.05),abs_tol=1e-6)
    assert manifest['rectangular_kinematic_cuts_applied'] is False
    assert manifest['target_divisor']==.584
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (out/'cache').iterdir()}
    result=subprocess.run(command,text=True,capture_output=True)
    assert result.returncode and 'refuses existing' in result.stderr
    assert hashes=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (out/'cache').iterdir()}
    print('PASS forward export: rectangular cuts deferred; mass cut, units, truth matching, target once, overwrite refusal')
    print(out)


if __name__=='__main__': main()
