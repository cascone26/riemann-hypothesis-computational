import numpy as np, os
from scipy.signal import argrelmax

# Is the twisted-Landau comb the MOBIUS/squarefree comb?
# G_t(u) = (1/N) sum_k e^{-2i arg zeta'(rho_k)} e^{i u gamma_k}, 500 zeros.
# Tests: (1) top-peak hit rates vs three mark sets (prime powers / squarefree
# composites / non-squarefree>3), with circular-shift p-values.
# (2) sign structure: Re(e^{i c} G_t(log n)) vs mu(n) for squarefree n.

DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros500.npy"))
ph = np.load(os.path.join(DIR, "argzp500.npy"))
us = np.linspace(0.5, 4.6, 8201)
E = np.exp(1j * np.outer(us, Z))
Gt = (E * np.exp(-2j * ph)[None, :]).sum(axis=1) / len(Z)
A = np.abs(Gt)

def factor(n):
    f = {}
    d = 2; m = n
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f

def mu(n):
    f = factor(n)
    if any(e > 1 for e in f.values()): return 0
    return -1 if len(f) % 2 else 1

NMAX = 99
marks_pp = [n for n in range(2, NMAX) if len(factor(n)) == 1]                       # p^m
marks_sf = [n for n in range(2, NMAX) if mu(n) != 0 and len(factor(n)) >= 2]        # squarefree composite
marks_nsf = [n for n in range(4, NMAX) if mu(n) == 0]                               # non-squarefree
rng = np.random.default_rng(13)

ix = argrelmax(A, order=40)[0]
top = ix[np.argsort(A[ix])[::-1][:20]]
peak_us = np.sort(us[top])

def hitstats(name, nums):
    mk = np.array([np.log(n) for n in nums if 0.55 < np.log(n) < 4.55])
    hits = sum(1 for u in peak_us if np.min(np.abs(mk - u)) < 0.012)
    # circular-shift control on the mark set
    ctrl = []
    for _ in range(2000):
        sh = (mk - 0.55 + rng.uniform(0.05, 3.9)) % 4.0 + 0.55
        ctrl.append(sum(1 for u in peak_us if np.min(np.abs(sh - u)) < 0.012))
    ctrl = np.array(ctrl)
    p = float(np.mean(ctrl >= hits))
    print(f"  {name:>24}: {len(mk):2d} marks | hits {hits}/20 | ctrl {ctrl.mean():.1f}±{ctrl.std():.1f} | p={p:.4f}")

print("top-20 twisted peaks at u =", np.round(peak_us, 3))
hitstats("prime powers", marks_pp)
hitstats("squarefree composites", marks_sf)
hitstats("non-squarefree (4,8,9..)", marks_nsf)
hitstats("ALL squarefree (incl p)", sorted(set(marks_pp) & set(range(2, NMAX)) - set()) )

sf_all = [n for n in range(2, NMAX) if mu(n) != 0]
hitstats("all squarefree n", sf_all)

print("\nsign test: best global rotation c, corr( Re(e^{ic}G(log n)), mu(n)/sqrt(n)-style ) over squarefree n")
lgn = np.log(np.array(sf_all, dtype=float))
Gn = (np.exp(1j * np.outer(lgn, Z)) * np.exp(-2j * ph)[None, :]).sum(axis=1) / len(Z)
mus = np.array([mu(n) for n in sf_all], dtype=float)
best = (0, -2)
for c in np.linspace(-np.pi, np.pi, 361):
    r = np.corrcoef(np.real(np.exp(1j * c) * Gn), mus)[0, 1]
    if r > best[1]: best = (float(c), float(r))
print(f"  vs mu(n):        best c={best[0]:+.2f} r={best[1]:+.3f}")
best2 = (0, -2)
for c in np.linspace(-np.pi, np.pi, 361):
    r = np.corrcoef(np.real(np.exp(1j * c) * Gn), mus / np.sqrt(sf_all))[0, 1]
    if r > best2[1]: best2 = (float(c), float(r))
print(f"  vs mu(n)/sqrt(n): best c={best2[0]:+.2f} r={best2[1]:+.3f}")
shuf = mus.copy()
ctrl = []
for _ in range(500):
    rng.shuffle(shuf)
    m = max(np.corrcoef(np.real(np.exp(1j * c) * Gn), shuf)[0, 1] for c in np.linspace(-np.pi, np.pi, 37))
    ctrl.append(m)
print(f"  shuffled-mu control: mean {np.mean(ctrl):.3f} max {np.max(ctrl):.3f}")
