import sys, csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

outdir = Path(sys.argv[1] if len(sys.argv)>1 else "sweep_out")
r1 = sys.argv[2] if len(sys.argv)>2 else None
r2 = sys.argv[3] if len(sys.argv)>3 else None
a  = sys.argv[4] if len(sys.argv)>4 else None

def load(d):
    xs,v,se=[],[],[]
    with (outdir/f"results_d{d}.csv").open() as f:
        r=csv.DictReader(f)
        for row in r:
            xs.append(float(row["sqrtN"]))
            v.append(float(row["volume"]))
            se.append(float(row["se"]))
    xs=np.array(xs); v=np.array(v); se=np.array(se)
    idx=np.argsort(xs)
    return xs[idx], v[idx], se[idx]

fig,axs=plt.subplots(2,1,figsize=(9,10),sharex=True)

series=[(3,"o","tab:blue",  -0.7),
        (5,"s","tab:orange",  0.0),
        (10,"^","tab:green",  0.7)]

for d,mark,col,xshift in series:
    xs,v,se=load(d)
    x = xs + xshift
    axs[0].errorbar(x, v, yerr=se, fmt=mark, ms=3.5, mfc="white",
                    ecolor=col, mec=col, mew=1, color=col,
                    capsize=2, elinewidth=0.9, alpha=0.95, label=f"d={d}")
    axs[0].plot(x, v, color=col, linewidth=1.0, alpha=0.6)
axs[0].set_title("Estimated volume vs $\\sqrt{N}$")
axs[0].set_ylabel("Estimated volume")
axs[0].grid(True, alpha=0.25)
axs[0].legend(ncol=3, frameon=False, loc="upper right")

for d,mark,col,xshift in series:
    xs,v,se=load(d)
    x = xs + xshift
    axs[1].plot(x, se, mark, ms=3.5, mfc="white", mec=col, mew=1,
                color=col, alpha=0.95, label=f"d={d}")
    axs[1].plot(x, se, color=col, linewidth=1.0, alpha=0.6)
axs[1].set_title("Statistical uncertainty (SE) vs $\\sqrt{N}$")
axs[1].set_xlabel("$\\sqrt{N}$")
axs[1].set_ylabel("SE")
axs[1].grid(True, alpha=0.25)
axs[1].legend(ncol=3, frameon=False, loc="upper right")

if r1 and r2 and a:
    fig.suptitle(f"Monte Carlo (stone-throwing) for intersection volume — r1={r1}, r2={r2}, a={a}", y=0.99, fontsize=12)

plt.tight_layout(rect=(0,0,1,0.98))
(fig if hasattr(fig, "savefig") else plt).savefig(outdir/"convergence.png", dpi=180)
plt.savefig(outdir/"convergence.pdf", dpi=180)

