"""Direct raw-bank exposure comparison; unassigned counters stay unassigned."""
from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;D=P/'raw4305';F=P/'figures'
a=np.load(D/'scaler_modules.npz');b=a['event'];v=a['values'];e=np.fromfile(D/'physics.bin',dtype=np.uint64).reshape(-1,4)
freq=next(r['relative_ticks_per_scaler_second'] for r in json.loads((P/'event_clock_checks.json').read_text())['components'] if r['run']==4305)
t=e[b-1,2].astype(float);event=(t-t[0])/freq;scaler=(v[:,6,31]-v[0,6,31])*1e-6
fig,axs=plt.subplots(3,1,figsize=(11.8,8),sharex=True,layout='constrained')
axs[0].plot(b,event-scaler,color='#176b9b',lw=1.7);axs[0].set_ylabel('Event elapsed - scaler elapsed (s)')
axs[1].step(b,v[:,1,31]-v[:,6,14],where='post',color='#cf6a18');axs[1].set_ylabel('Slot7/ch31 - EDTM slot12/ch14\n(raw counts)')
axs[2].step(b,b-v[:,4,16],where='post',color='#23764a');axs[2].set_ylabel('Recorded event number - raw L1');axs[2].set_xlabel('Recorded event number at raw scaler snapshot')
for ax in axs:
    ax.axvspan(105080,106127,color='#b3294d',alpha=.3);ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.15)
fig.suptitle('4305: the exposure discontinuity is already in the raw banks\nUnassigned slot7/ch31 equals slot9/ch31 throughout; their signal attribution remains unverified.',fontsize=13)
fig.savefig(F/'raw4305_exposure_step.png',dpi=180);fig.savefig(F/'raw4305_exposure_step.pdf');plt.close(fig)
