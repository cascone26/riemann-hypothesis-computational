#!/usr/bin/env python3
"""
sc_oscillation_map.py — Dense lsq scan to map the sub_cancel oscillation.

We found sc(lsq=200) = 4.16e-6 (WAY below L^{-1.92} power law).
Is this a deep local minimum? What's the oscillation pattern?

Scan: lsq = 50, 70, 80, 90, 100, 120, 140, 150, 170, 190, 200, 220, 250, 300
N = 5*lsq for each (except lsq=300: N=1200 to save time).
Also measure lam0 and dim to track saturation.
"""
import os, json, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
THRESH = 0.02
OUTFILE = os.path.expanduser("~/riemann/results/sc_oscillation_map.json")


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
    exp_dim = 2 * (lsq - 2) + 1
    if dim > 0:
        Vn = V[:, null_mask]
        ABp = Vn.T @ AB @ Vn
        Pp  = Vn.T @ P  @ Vn
        sc = float(np.linalg.norm(ABp - Pp, 'fro') / (np.linalg.norm(ABp, 'fro') + 1e-100))
    else:
        sc = float('nan')
    return sc, lam0, dim, exp_dim


def main():
    # Dense scan — include known interesting points
    lsq_list = [50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
                150, 160, 170, 180, 190, 200, 210, 220, 230, 250, 300]
    # N = 5*lsq except for large lsq
    def N_for(lsq):
        return min(5 * lsq, 1500)  # cap at N=1500 (matrix 3001x3001)

    results = []
    t0 = time.time()
    for lsq in lsq_list:
        N = N_for(lsq)
        L = float(mp.log(lsq))
        t1 = time.time()
        sc, lam0, dim, exp_dim = measure(lsq, N)
        dt = time.time() - t1
        # Power law prediction
        sc_pred = 3.073e-3 / L ** 1.9151
        ratio = sc / sc_pred if sc_pred > 0 and not np.isnan(sc) else float('nan')
        print(f"lsq={lsq:4d} L={L:.4f} N={N:5d}: sc={sc:.4e}  pred={sc_pred:.4e}  ratio={ratio:.3f}  lam0={lam0:.2e}  dim={dim}/{exp_dim}  ({dt:.1f}s)", flush=True)
        results.append({'lsq': lsq, 'L': round(L,5), 'N': N, 'sc': sc,
                        'sc_pred': round(sc_pred, 8), 'ratio': round(ratio, 5) if not np.isnan(ratio) else None,
                        'lam0': lam0, 'dim': dim, 'exp_dim': exp_dim, 'dt': round(dt,1)})

    out = {'lsq_list': lsq_list, 'results': results, 'total_time': round(time.time()-t0,1)}
    with open(OUTFILE,'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main()
