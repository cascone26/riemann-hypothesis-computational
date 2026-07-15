#!/usr/bin/env python3
"""
sigma_sharpening.py — Does the sub_cancel minimum at σ=0.5 get SHARPER with lsq?

From sigma_fine_scan.py at lsq=50:
  sc(0.499)/sc(0.500) = 11.6
  sc(0.501)/sc(0.500) = 11.9

Question: does this ratio GROW with lsq (sharpening phase transition)?
If ratio → ∞ as lsq → ∞, it means σ=0.5 is THE unique minimizer.
If ratio → constant, it's still a minimum but the contrast stays bounded.

Also tests σ = 0.5 ± 0.001, 0.5 ± 0.005, 0.5 ± 0.01 at lsq = 20, 30, 50, 80, 150.
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

OUTFILE = os.path.expanduser("~/riemann/results/sigma_sharpening.json")


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


def measure_sc(lambda_sq, N_use, sigma):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N_use, N_use + 1, dtype=float)
    ns  = idx
    R   = TWO_PI * N_use / L + 80.0
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
    AB = arch + B
    if VM:
        lkpp = np.array([np.log(k) for k, _ in VM])
        wts  = np.array([lk * (k ** (-sigma)) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))
    tau = AB - P
    tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    null_mask = np.abs(ev) < THRESH
    dim = null_mask.sum()
    if dim == 0:
        return float('nan'), int(dim), float(ev[0])
    V_null = V[:, null_mask]
    AB_proj = V_null.T @ AB @ V_null
    P_proj  = V_null.T @ P  @ V_null
    sc = float(np.linalg.norm(AB_proj - P_proj, 'fro') / (np.linalg.norm(AB_proj, 'fro') + 1e-100))
    return sc, int(dim), float(ev[0])


def main():
    # lsq values to test
    lsq_list = [20, 30, 50, 80, 150]
    # sigma values: center ± offsets
    deltas = [0.0, 0.001, 0.005, 0.01, 0.05]

    results = []
    t0 = time.time()

    print(f"sigma_sharpening: lsq={lsq_list}")
    print(f"Testing sigma offsets: {deltas}")

    for lsq in lsq_list:
        N_use = int(np.ceil(3.0 * lsq))
        L = float(mp.log(lsq))
        print(f"\nlsq={lsq}, N={N_use}, L={L:.4f}", flush=True)

        sc_by_delta = {}
        for delta in deltas:
            for sign in ([0] if delta == 0.0 else [+1, -1]):
                sigma = 0.5 + sign * delta
                t1 = time.time()
                sc, dim, lam_min = measure_sc(lsq, N_use, sigma)
                dt = time.time() - t1
                key = f'{delta:+.3f}' if sign >= 0 else f'{-delta:+.3f}'
                print(f"  σ=0.5{sign*delta:+.3f}: sc={sc:.4e}  dim={dim}  lam={lam_min:.3e}  ({dt:.1f}s)", flush=True)
                sc_by_delta[(delta, sign)] = {'sigma': sigma, 'sc': sc, 'dim': dim, 'lam_min': lam_min}

        # Compute ratios
        sc_0 = sc_by_delta[(0.0, 0)]['sc']
        print(f"  Ratios sc(σ=0.5±Δ) / sc(σ=0.5):")
        for delta in deltas:
            if delta == 0.0:
                continue
            sc_plus  = sc_by_delta[(delta, +1)]['sc']
            sc_minus = sc_by_delta[(delta, -1)]['sc']
            r_plus  = sc_plus  / sc_0 if sc_0 > 0 else float('inf')
            r_minus = sc_minus / sc_0 if sc_0 > 0 else float('inf')
            print(f"    Δ={delta}: ratio(+)={r_plus:.2f}  ratio(-)={r_minus:.2f}", flush=True)

        results.append({'lsq': lsq, 'N': N_use, 'L': round(L, 4),
                       'sc_by_sigma': {f'{s}:{sg}': v for (s,sg),v in sc_by_delta.items()}})

    out = {
        'lsq_list': lsq_list,
        'deltas': deltas,
        'results': results,
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
