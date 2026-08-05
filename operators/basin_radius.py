"""Basin radius: perturb ALL 5000 zeros by iid noise of amplitude a spacings; annealed solve;
recovery vs a. Where does the single basin end?"""
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
for delta in (0.8, 0.5, 0.35, 0.25, 0.18, 0.12):
    w = As*np.exp(-0.5*(Ls*delta)**2); keep = w > 1e-8
    banks[delta] = (Ls[keep], w[keep])
def S_Sp(t, delta):
    Lk, Ak = banks[delta]
    S = np.zeros_like(t); Sp = np.zeros_like(t)
    for L, a in zip(Lk, Ak):
        S += -a*np.sin(t*L); Sp += -a*L*np.cos(t*L)
    return S, Sp
def solve_anneal(t):
    for delta in (0.8, 0.5, 0.35, 0.25, 0.18, 0.12):
        for _ in range(3):
            S, Sp = S_Sp(t, delta)
            thp = 0.5*np.log(t/(2*np.pi))/np.pi
            Jc = np.maximum(thp + Sp, 0.3*thp)
            t = t - np.clip(((th_of(t)/np.pi + 1 + S) - (k-0.5))/Jc, -2.5, 2.5)
    return t
rng = np.random.default_rng(3)
for a in (0.5, 1.0, 2.0, 4.0, 8.0):
    t0 = g + a*sp*rng.standard_normal(N)
    t0 = np.sort(t0)
    t = solve_anneal(t0.copy())
    e = (t-g)/sp
    rms = np.sqrt(np.mean(e**2)); worst = np.abs(e).max()
    frac = np.mean(np.abs(e) < 0.05)
    print(f"a={a:.1f} sp: recovered rms {rms:.4f} | worst {worst:.3f} | frac<0.05sp {frac:.3f}", flush=True)
print("DONE")
