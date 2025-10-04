import csv, sys
from pathlib import Path
import matplotlib.pyplot as plt
from math import pi, gamma

outdir = Path(sys.argv[1] if len(sys.argv)>1 else "result")

def load(d):
    xs,v,se=[],[],[]
    with (outdir/f"results_d{d}.csv").open() as f:
        r=csv.DictReader(f)
        for row in r:
            xs.append(float(row["sqrtN"]))
            v.append(float(row["volume"]))
            se.append(float(row["se"]))
    return xs,v,se

def Vd(d):
    return (pi**(d/2))/gamma(d/2+1)

fig,axs=plt.subplots(2,1,figsize=(7,8),sharex=True)

for d,mark in [(3,"o"),(5,"s"),(10,"^")]:
    xs,v,se=load(d)
    axs[0].errorbar(xs,v,yerr=se,fmt=mark,label=f"d={d}",capsize=3,linewidth=1)
    axs[0].plot(xs,[Vd(d)]*len(xs),"--",linewidth=1)
axs[0].set_ylabel("Estimated volume")
axs[0].legend()

for d,mark in [(3,"o"),(5,"s"),(10,"^")]:
    xs,v,se=load(d)
    axs[1].plot(xs,se,mark,label=f"d={d}")
axs[1].set_xlabel(r"$\sqrt{N}$")
axs[1].set_ylabel("Statistical uncertainty (SE)")
axs[1].legend()

plt.tight_layout()
plt.savefig(outdir/"convergence.png",dpi=150)
