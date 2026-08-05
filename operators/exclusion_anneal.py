"""Exclusion map v2: same 15 defects, but solved with the FULL annealing schedule
(delta 0.5 -> 0.12) that the global recovery uses. Do the stuck points dissolve?"""
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
banks = {}
for delta in (0.5, 0.35, 0.25, 0.18, 0.12):
    w = As*np.exp(-0.5*(Ls*delta)**2); keep = w > 1e-8
    banks[delta] = (Ls[keep], w[keep])
def S_Sp(t, delta):
    Lk, Ak = banks[delta]
    S = np.zeros_like(t); Sp = np.zeros_like(t)
    for L, a in zip(Lk, Ak):
        S += -a*np.sin(t*L); Sp += -a*L*np.cos(t*L)
    return S, Sp
def solve_anneal(t):
    for delta in (0.5, 0.35, 0.25, 0.18, 0.12):
        for _ in range(3):
            S, Sp = S_Sp(t, delta)
            thp = 0.5*np.log(t/(2*np.pi))/np.pi
            Jc = np.maximum(thp + Sp, 0.3*thp)
            t = t - np.clip(((th_of(t)/np.pi + 1 + S) - (k-0.5))/Jc, -1.5, 1.5)
    return t
floor = 0.0059
res = []
for idx in (500, 1500, 2500, 3500, 4500):
    for fam in ("pair", "triple", "block"):
        t = g.copy()
        if fam == "pair":
            m = 0.5*(g[idx]+g[idx+1]); t[idx] = t[idx+1] = m
        elif fam == "triple":
            m = (g[idx]+g[idx+1]+g[idx+2])/3; t[idx] = t[idx+1] = t[idx+2] = m
        else:
            t[idx:idx+5] += 0.4*sp[idx:idx+5]
        t = solve_anneal(t)
        e = (t-g)/sp
        loc = np.abs(e[max(0,idx-2):idx+7]).max()
        ok = loc < 3*floor
        res.append(ok)
        print(f"h~{g[idx]:.0f} {fam:>6}: max local err {loc:.4f} ({'REPAIRED' if ok else 'STUCK'})", flush=True)
print("TOTAL REPAIRED: %d/15" % sum(res))
print("DONE")
