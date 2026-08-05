"""Exclusion map: defect families at several heights; does fine-delta self-consistency
always repair to floor? Families: pair-merge, triple-merge, block-shift (5 zeros +0.4sp)."""
import numpy as np
from scipy.special import loggamma as slg
from sympy import primerange
g = np.load("zeros10k.npy")[:5000]
N = len(g); k = np.arange(1, N+1)
th_of = lambda tt: np.imag(slg(0.25+0.5j*np.asarray(tt))) - np.asarray(tt)/2*np.log(np.pi)
sp = np.gradient(g)
osc = []
for p in primerange(2, 500000):
    lp = np.log(p); m = 1
    while p**m < 500000:
        osc.append((lp*m, 1.0/(np.pi*np.sqrt(p**m)*m))); m += 1
Ls = np.array([o[0] for o in osc]); As = np.array([o[1] for o in osc])
W12 = As*np.exp(-0.5*(Ls*0.12)**2); keep = W12 > 1e-8
Lk, Ak = Ls[keep], W12[keep]
def S_Sp(t):
    S = np.zeros_like(t); Sp = np.zeros_like(t)
    for L, a in zip(Lk, Ak):
        S += -a*np.sin(t*L); Sp += -a*L*np.cos(t*L)
    return S, Sp
def solve(t, iters=6):
    for _ in range(iters):
        S, Sp = S_Sp(t)
        thp = 0.5*np.log(t/(2*np.pi))/np.pi
        Jc = np.maximum(thp + Sp, 0.3*thp)
        t = t - np.clip(((th_of(t)/np.pi + 1 + S) - (k-0.5))/Jc, -1.5, 1.5)
    return t
floor = 0.0059
for idx in (500, 1500, 2500, 3500, 4500):
    for fam in ("pair", "triple", "block"):
        t = g.copy()
        if fam == "pair":
            m = 0.5*(g[idx]+g[idx+1]); t[idx] = t[idx+1] = m
        elif fam == "triple":
            m = (g[idx]+g[idx+1]+g[idx+2])/3; t[idx] = t[idx+1] = t[idx+2] = m
        else:
            t[idx:idx+5] += 0.4*sp[idx:idx+5]
        t = solve(t)
        e = (t-g)/sp
        loc = np.abs(e[max(0,idx-2):idx+7]).max()
        print(f"h~{g[idx]:.0f} {fam:>6}: max local err after solve {loc:.4f} ({'REPAIRED' if loc < 3*floor else 'STUCK'})", flush=True)
print("DONE")
