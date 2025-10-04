import sys, math, random

def sq(x): return x*x

if len(sys.argv) != 6:
    sys.stderr.write(f"Usage: {sys.argv[0]} <d> <N> <r1> <r2> <a>\n")
    sys.exit(1)

d = int(sys.argv[1])
N = int(sys.argv[2])
r1 = float(sys.argv[3])
r2 = float(sys.argv[4])
a = float(sys.argv[5])

if d < 1 or N <= 0 or r1 <= 0.0 or r2 <= 0.0:
    sys.stderr.write("Bad inputs: d>=1, N>0, r1>0, r2>0\n")
    sys.exit(1)

maxr = r1 if r1 >= r2 else r2
xmin = (-r1) if (-r1) < (a - r2) else (a - r2)
xmax = ( r1) if ( r1) > (a + r2) else (a + r2)
Lx = xmax - xmin
Lperp = 2.0 * maxr

Vbox = Lx
for _ in range(1, d): Vbox *= Lperp

rng = random.Random(0xC0FFEE)
u0 = lambda: rng.uniform(xmin, xmax)
u1 = lambda: rng.uniform(-maxr, maxr)

r1sq = r1*r1
r2sq = r2*r2

hits = 0
for _ in range(N):
    x0 = u0()
    d1 = x0*x0
    t = x0 - a
    d2 = t*t
    for _ in range(1, d):
        xj = u1()
        y = xj*xj
        d1 += y
        d2 += y
    if d1 <= r1sq and d2 <= r2sq:
        hits += 1

phat = hits / float(N)
volume = Vbox * phat
stdev = Vbox * math.sqrt(phat * (1.0 - phat) / float(N)) if 0.0 < phat < 1.0 else 0.0

print(f"(r1,r2): {r1} {r2}")
print(f"(d,N,a): {d} {N} {a}")
print(f"volume: {volume}")
print(f"stat uncertainty: {stdev}")
