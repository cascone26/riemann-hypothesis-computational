#!/usr/bin/env python3
"""
lattice_highN_verify.py — Higher-N local verification of key lattice points.
N = min(3*lsq, 1200) — better than previous min(4*lsq,600) cap.
"""
import time, json
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 50000
THRESH = 0.02
OUTFILE = '/Users/scones/projects/riemann/results/lattice_highN_mac.out'

def vm_table(kmax):
    tab = []
    for k in range(2, int(kmax) + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab

def measure(lsq, N):
    L  = float(mp.log(lsq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lsq)
    idx = np.arange(-N, N + 1, dtype=float)
    ns  = idx
    R   = TWO_PI * N / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]; w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI
    d = (TWO_PI * ns / L)[:, None] - grid[None, :]
    tiny = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L),
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L)
                 * 2 * np.sin(0.5 * d * L) / (d + 1e-300 * tiny))
    arch = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
    B = (32.0 * L * sh * (L ** 2 - 16.0 * np.pi ** 2 * np.outer(ns, ns)) /
         np.outer(L ** 2 + 16.0 * np.pi ** 2 * ns ** 2,
                  L ** 2 + 16.0 * np.pi ** 2 * ns ** 2))
    if VM:
        lkpp = np.array([np.log(k) for k, _ in VM])
        wts  = np.array([lk * k ** (-0.5) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))
    AB = arch + B
    tau = (AB - P); tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    lam0 = float(ev[0])
    null_mask = np.abs(ev) < THRESH
    dim = int(null_mask.sum())
    if dim > 0:
        Vn = V[:, null_mask]
        ABp = Vn.T @ AB @ Vn
        Pp  = Vn.T @ P  @ Vn
        sc = float(np.linalg.norm(ABp - Pp, 'fro') / (np.linalg.norm(ABp, 'fro') + 1e-100))
    else:
        sc = float('nan')
    return sc, lam0, dim

tests = [
    (80,  'control         '),
    (160, 'γ₁,n=-1 [SPIKE] '),
    (180, 'γ₂,n=-6 [sub]   '),
    (200, 'γ₁,n=0  [MIN]   '),
    (209, 'γ₂,n=-5 [MIN?]  '),
    (243, 'γ₂,n=-4 [SPIKE] '),
    (250, 'γ₁,n=1  [MIN]   '),
    (282, 'γ₂,n=-3 [sub]   '),
    (312, 'γ₁,n=2  [MIN]   '),
]

with open(OUTFILE, 'w') as f:
    f.write(f"Lattice highN verification (N=min(3*lsq,1200)) {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"{'lsq':>5} {'label':20} {'N':>5} {'sc':>12} {'pred':>12} {'ratio':>7} {'type':>6} {'dt':>5}\n")
    f.write("-" * 85 + "\n")

print(f"Lattice highN verification — N=min(3*lsq,1200)")
print(f"{'lsq':>5} {'label':20} {'N':>5} {'sc':>12} {'pred':>12} {'ratio':>7} {'type':>6} {'dt':>5}")
print("-" * 85)

results = []
t0 = time.time()
for lsq, label in tests:
    N = min(3 * lsq, 1200)
    L = float(mp.log(lsq))
    t1 = time.time()
    sc, lam0, dim = measure(lsq, N)
    dt = time.time() - t1
    sc_pred = 3.073e-3 / L ** 1.9151
    ratio = sc / sc_pred if sc_pred > 0 and not np.isnan(sc) else float('nan')
    typ = "SPIKE" if ratio > 1.5 else ("MIN" if ratio < 0.2 else "sub")
    line = f"{lsq:5d} {label:20} {N:5d}: sc={sc:.4e}  pred={sc_pred:.4e}  ratio={ratio:7.3f}  {typ:5}  ({dt:.0f}s)"
    print(line, flush=True)
    with open(OUTFILE, 'a') as f:
        f.write(line + "\n")
    results.append({'lsq': lsq, 'label': label.strip(), 'N': N, 'sc': sc, 'ratio': ratio, 'dt': round(dt)})

print(f"\nTotal: {time.time()-t0:.1f}s")
with open(OUTFILE.replace('.out', '.json'), 'w') as f:
    json.dump(results, f, indent=2)
print(f"Done → {OUTFILE}")
