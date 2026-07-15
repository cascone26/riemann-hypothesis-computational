#!/usr/bin/env python3
"""
gamma3_dense.py — Dense scan around γ₃ resonance at lsq=625.5.

γ₃ = 25.0109, γ₃² = 625.5  → main minimum at lsq ≈ 625

At N=1000 (Mac), we have N ≈ 1.6*lsq at lsq=625 → undersaturated.
BUT: γ₃ scan showed lam0 rapidly dropping as we approach 625:
  lsq=552: lam0=-1.07e-6
  lsq=575: lam0=-1.35e-7
These lam0 values are 1000x smaller than lsq=486 (-7.7e-4).

Dense scan from lsq=580 to lsq=680 at step=5, N=1000.
Watch for:
  - lam0 minimum at lsq≈625 (deep resonance)
  - sc ratio minimum (even at undersaturation, sc may be low at exact resonance)
"""
import time, sys
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 50000
THRESH = 0.02
N      = 1000
OUTFILE = '/Users/scones/projects/riemann/results/gamma3_dense.out'

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

def measure(lsq):
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

def main():
    lsq_vals = list(range(580, 685, 5))  # 580, 585, ..., 680

    sc_pred_ref = 3.073e-3 / float(mp.log(625)) ** 1.9151

    print(f"γ₃ dense scan: lsq=580-680, N={N}")
    print(f"γ₃²=625.5 → main resonance at lsq≈625")
    print(f"pred sc at lsq=625: {sc_pred_ref:.4e}")
    print(f"{'lsq':>5} {'sc':>12} {'ratio':>7} {'lam0':>12} {'dim':>5} {'dt':>5}")
    print("-" * 55)
    with open(OUTFILE, 'w') as f:
        f.write(f"gamma3_dense N={N} {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"γ₃²=625.5\n")

    for lsq in lsq_vals:
        t0 = time.time()
        sc, lam0, dim = measure(lsq)
        dt = time.time() - t0
        sc_pred = 3.073e-3 / float(mp.log(lsq)) ** 1.9151
        ratio = sc / sc_pred if sc_pred > 0 and not np.isnan(sc) else float('nan')
        label = " ← γ₃ n=0!" if abs(lsq - 625) <= 3 else ""
        line = f"{lsq:5d} {sc:12.4e} {ratio:7.3f} {lam0:12.4e} {dim:5d}  ({dt:.0f}s){label}"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")

if __name__ == "__main__":
    main()
