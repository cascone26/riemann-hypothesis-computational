import numpy as np

# Phase-only Mobius claim at 10^4 zeros (gamma up to ~9878):
#  - weight ladder r vs mu(n)/sqrt(n) at n<99 / n<300 / n<1000 / n<3000
#  - squarefree comb peak test in u in [4.5, 8.0] (n ~ 90..3000), tight tolerance
# Convergence vs N in {2000, 5000, 10000}.

Z = np.load("zeros10k.npy"); ph = np.load("argzp10k.npy"); mag = np.load("abszp10k.npy")

def factor(n):
    f = {}; d = 2; m = n
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f

def mu(n):
    f = factor(n)
    return 0 if any(e > 1 for e in f.values()) else (-1 if len(f) % 2 else 1)

print(f"zeros loaded: {len(Z)} (gamma_max={Z[-1]:.1f})")
for NMAX in [99, 300, 1000, 3000]:
    sf = [n for n in range(2, NMAX) if mu(n) != 0]
    target = np.array([mu(n) / np.sqrt(n) for n in sf])
    lgn = np.log(np.array(sf, float))
    print(f"--- n<{NMAX} ({len(sf)} sf targets) ---")
    print(f"{'N':>7} {'gonek':>8} {'phase1':>8} {'gap':>8}")
    for N in [2000, 5000, 10000]:
        E = np.exp(1j * np.outer(lgn, Z[:N]))
        rs = []
        for w in [np.exp(-1j * ph[:N]) / mag[:N], np.exp(-1j * ph[:N])]:
            Gn = (E * w[None, :]).sum(axis=1) / np.sum(np.abs(w))
            rs.append(max(float(np.corrcoef(np.real(np.exp(1j * c) * Gn), target)[0, 1])
                          for c in np.linspace(-np.pi, np.pi, 241)))
        print(f"{N:>7} {rs[0]:>8.4f} {rs[1]:>8.4f} {rs[0] - rs[1]:>8.4f}", flush=True)

from scipy.signal import argrelmax
rng = np.random.default_rng(53)
us = np.linspace(4.5, 8.0, 70001)
E = np.exp(1j * np.outer(us, Z))
Gt = np.abs((E * np.exp(-2j * ph)[None, :]).sum(axis=1)) / len(Z)
ix = argrelmax(Gt, order=80)[0]
top = np.sort(ix[np.argsort(Gt[ix])[::-1][:40]])
sfm = np.log(np.array([n for n in range(90, 3000) if mu(n) != 0], float))
nsfm = np.log(np.array([n for n in range(90, 3000) if mu(n) == 0], float))
tol = 0.00025
h_sf = sum(1 for i in top if np.min(np.abs(sfm - us[i])) < tol)
h_nsf = sum(1 for i in top if np.min(np.abs(nsfm - us[i])) < tol)
ctrl = []
for _ in range(1000):
    sh = 4.5 + (sfm - 4.5 + rng.uniform(0.05, 3.4)) % 3.5
    ctrl.append(sum(1 for i in top if np.min(np.abs(sh - us[i])) < tol))
ctrl = np.array(ctrl)
print(f"\ncomb u in [4.5,8.0] (n~90-3000), top-40 peaks, tol {tol}:")
print(f"  squarefree {h_sf}/40 (ctrl {ctrl.mean():.1f}±{ctrl.std():.1f}, p={float(np.mean(ctrl >= h_sf)):.4f}) | non-squarefree {h_nsf}/40")
print("DONE", flush=True)
