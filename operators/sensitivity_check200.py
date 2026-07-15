#!/usr/bin/env python3
"""Quick: lsq=200 N=1000 (5x) sensitivity check."""
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
GAMMA1 = 14.134725
OUTFILE = os.path.expanduser("~/riemann/results/sensitivity_check200.json")


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


def run(lsq, N):
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
    Vn = V[:, null_mask]
    ABp = Vn.T @ AB @ Vn
    Pp  = Vn.T @ P  @ Vn
    sc = float(np.linalg.norm(ABp - Pp, 'fro') / (np.linalg.norm(ABp, 'fro') + 1e-100))
    # Sensitivity
    s0 = 0.25 + 0.5j * grid
    def contrib(sigma):
        rho = sigma + 1j * GAMMA1
        rho_bar = sigma - 1j * GAMMA1
        return (np.real(1.0 / (s0 - rho) + 1.0 / rho) +
                np.real(1.0 / (s0 - rho_bar) + 1.0 / rho_bar))
    for ds in [0.001, 0.002, 0.005]:
        dOm = contrib(0.5 + ds) - contrib(0.5)
        dA = ((E * (w * dOm)) @ E.conj().T).real / (2.0 * np.pi)
        tau2 = (AB + dA - P); tau2 = (tau2 + tau2.T) / 2.0
        lam_g = float(scipy.linalg.eigvalsh(tau2, subset_by_index=[0,0])[0])
        print(f"  Δσ={ds}: lam_ghost={lam_g:.4e}  slope={( lam_g - lam0)/ds:.4e}", flush=True)
    print(f"lam0={lam0:.4e}  sc={sc:.4e}  dim={dim}/{2*(lsq-2)+1}", flush=True)
    return {'lsq': lsq, 'N': N, 'L': round(L,4), 'lam0': lam0, 'sc': sc, 'dim': dim}


t0 = time.time()
lsq = 200; N = 5 * lsq
print(f"lsq={lsq} N={N} (5x)")
r = run(lsq, N)
print(f"Done in {time.time()-t0:.1f}s")
with open(OUTFILE, 'w') as f:
    json.dump(r, f, indent=2)
