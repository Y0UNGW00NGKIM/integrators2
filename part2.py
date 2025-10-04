import sys, math, csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import qmc
from scipy.special import gamma

def V5():
    return (math.pi**(2.5))/gamma(3.5)

def parse_list_int(s, default):
    if s is None or s=="default":
        return default
    return [int(t) for t in s.split(",") if t.strip()]

def grid_estimate(m):
    h=2.0/m
    xs=(-1.0+(np.arange(m)+0.5)*h)**2
    ps=(xs[:,None]+xs[None,:]).ravel()
    ps.sort()
    x2=xs
    c=0
    for i in range(m):
        si=x2[i]
        for j in range(m):
            sj=si+x2[j]
            for k in range(m):
                s3=sj+x2[k]
                rem=1.0-s3
                if rem<0.0:
                    continue
                c+=int(np.searchsorted(ps, rem, side="right"))
    vol=c*(h**5)
    Neff=m**5
    return vol, Neff

def pseudo_estimate(N, seed=12345):
    rng=np.random.default_rng(seed)
    b=1_000_000
    hits=0
    done=0
    while done<N:
        n=min(b, N-done)
        x=rng.uniform(-1.0,1.0,size=(n,5))
        s=np.sum(x*x,axis=1)
        hits+=int(np.count_nonzero(s<=1.0))
        done+=n
    vol=(hits/float(N))*(2.0**5)
    return vol, N

def sobol_estimate(N):
    if qmc is None:
        raise RuntimeError("scipy.stats.qmc.Sobol not available")
    eng=qmc.Sobol(d=5, scramble=False)
    b=1_000_000
    hits=0
    done=0
    while done<N:
        n=min(b, N-done)
        u=eng.random(n)
        x=2.0*u-1.0
        s=np.sum(x*x,axis=1)
        hits+=int(np.count_nonzero(s<=1.0))
        done+=n
    vol=(hits/float(N))*(2.0**5)
    return vol, N

def relerr(est, true): 
    return abs(est-true)/true

def run_all(outdir, Ns_mc, Ms_grid, seed=12345):
    outdir=Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    VT=V5()
    with (outdir/"grid.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["sqrtN","estimate","relerr"])
        for m in Ms_grid:
            est,Neff=grid_estimate(m)
            w.writerow([math.sqrt(Neff), est, relerr(est,VT)])
    with (outdir/"pseudo.csv").open("w",newline="") as f:
        w=csv.writer(f); w.writerow(["sqrtN","estimate","relerr"])
        for N in Ns_mc:
            est,Neff=pseudo_estimate(int(N), seed)
            w.writerow([math.sqrt(Neff), est, relerr(est,VT)])
    sobol_ok=False
    if qmc is not None:
        with (outdir/"sobol.csv").open("w",newline="") as f:
            w=csv.writer(f); w.writerow(["sqrtN","estimate","relerr"])
            for N in Ns_mc:
                est,Neff=sobol_estimate(int(N))
                w.writerow([math.sqrt(Neff), est, relerr(est,VT)])
        sobol_ok=True
    return sobol_ok

def loadcsv(p):
    a=np.genfromtxt(p, delimiter=",", names=True)
    return a["sqrtN"], a["estimate"], a["relerr"]

def plot_two_panel(outdir, title_note=""):
    outdir=Path(outdir)
    fig,axs=plt.subplots(2,1,figsize=(9,10),sharex=True)
    series=[]
    p=outdir/"grid.csv"
    if p.exists(): series.append(("Grid","tab:purple","D",p,0.0))
    p=outdir/"pseudo.csv"
    if p.exists(): series.append(("Pseudo","tab:blue","o",p,-0.6))
    p=outdir/"sobol.csv"
    if p.exists(): series.append(("Sobol","tab:orange","s",p,0.6))
    for name,col,mark,p,xshift in series:
        xs,est,err=loadcsv(p)
        idx=np.argsort(xs); xs=xs[idx]; est=est[idx]; err=err[idx]
        x=xs+xshift
        axs[0].plot(x, est, mark, ms=3.5, mfc="white", mec=col, mew=1, color=col, alpha=0.95, label=name)
        axs[0].plot(x, est, color=col, linewidth=1.0, alpha=0.6)
        axs[1].plot(x, err, mark, ms=3.5, mfc="white", mec=col, mew=1, color=col, alpha=0.95, label=name)
        axs[1].plot(x, err, color=col, linewidth=1.0, alpha=0.6)
    axs[0].axhline(V5(), ls="--", lw=1.0, color="gray")
    axs[0].set_title("5D unit ball volume vs $\\sqrt{N}$"+(" — "+title_note if title_note else ""))
    axs[0].set_ylabel("Estimated volume")
    axs[0].grid(True, alpha=0.25)
    axs[0].legend(ncol=len(series), frameon=False, loc="best")
    axs[1].set_title("Relative error vs $\\sqrt{N}$")
    axs[1].set_xlabel("$\\sqrt{N}$")
    axs[1].set_ylabel("Relative error")
    axs[1].grid(True, alpha=0.25)
    axs[1].legend(ncol=len(series), frameon=False, loc="best")
    plt.tight_layout()
    plt.savefig(outdir/"methods.png", dpi=180)

def main():
    outdir = sys.argv[1] if len(sys.argv)>1 else "part2_out"
    Ns_mc = parse_list_int(sys.argv[2] if len(sys.argv)>2 else "default",
                           [64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304,8388608,16777216])
    Ms_grid = parse_list_int(sys.argv[3] if len(sys.argv)>3 else "default",
                             [8,10,12,16,20,24,28,32,36,40])
    seed = int(sys.argv[4]) if len(sys.argv)>4 else 12345
    sobol_ok = run_all(outdir, Ns_mc, Ms_grid, seed)
    note = "Grid / Pseudo" + (" / Sobol" if sobol_ok else " (Sobol unavailable)")
    plot_two_panel(outdir, note)
if __name__=="__main__":
    main()

