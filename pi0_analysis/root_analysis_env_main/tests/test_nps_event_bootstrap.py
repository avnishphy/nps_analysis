"""Statistical identities and actual fitted-bootstrap checks; ROOT environment required."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from bootstrap_pi0_data import Fitter


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--library',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    fit=Fitter(args.library)
    rng=np.random.default_rng(729)
    mult=rng.poisson(1,(20000,300))
    count=mult[:,:200].sum(axis=1)
    sub=count-.25*mult[:,200:].sum(axis=1)
    counting_ratio=count.var(ddof=1)/200
    sideband_ratio=sub.var(ddof=1)/(200+.25**2*100)
    if abs(counting_ratio-1)>.05 or abs(sideband_ratio-1)>.05:
        raise RuntimeError('Poisson or fixed-sideband variance identity failed')
    x=(np.arange(200)+.5)*.002
    zero_classes=[];zero_y=[];failed=0;positive_bg=[];positive_y=[]
    for replica in range(2000):
        local=np.random.default_rng(np.random.SeedSequence([55,replica]))
        # Physical signal events plus disjoint coin/sideband accidentals.
        c=np.zeros(200);b=np.zeros(200)
        c[67]=local.poisson(200);c[25]=local.poisson(10);b[25]=local.poisson(60)
        y=c-b/6;variance=c+b/36
        code,final,info=fit.fit(y,variance)
        if code:failed+=1
        else:zero_classes.append(bool(info[0]));zero_y.append(final.sum())
        # Well-populated positive continuum; same disjoint Poisson event
        # counts feed both numerator and its observed variance.
        if replica>=250:continue
        background=12/(1+np.exp((x-.145)/.025))
        signal=30*np.exp(-.5*((x-.135)/.006)**2)
        c=local.poisson(background+signal+2).astype(float)
        b=local.poisson(np.full(200,12.)).astype(float)
        code,final,info=fit.fit(c-b/6,c+b/36)
        if not code:positive_bg.append(info[6]);positive_y.append(final.sum())
    if not any(zero_classes) or all(zero_classes):
        raise RuntimeError('Boundary bootstrap did not visit both zero and positive branches')
    if len(positive_bg)<240 or np.std(positive_bg,ddof=1)<=0:
        raise RuntimeError('Positive-background fits did not propagate fitted-background fluctuations')
    if failed or abs(np.var(zero_y,ddof=1)/(210+60/36)-1)>.13:
        raise RuntimeError('Zero-background toy fails timing-statistics closure')
    result=dict(poisson_variance_ratio=counting_ratio,sideband_variance_ratio=sideband_ratio,
                zero_background_replicas=2000,zero_branches=sum(zero_classes),positive_branches=len(zero_classes)-sum(zero_classes),
                zero_case_failures=failed,zero_case_variance=float(np.var(zero_y,ddof=1)),
                zero_case_timing_variance=210+60/36,
                positive_case_successful=len(positive_bg),positive_background_sd=float(np.std(positive_bg,ddof=1)),
                positive_signal_sd=float(np.std(positive_y,ddof=1)))
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
