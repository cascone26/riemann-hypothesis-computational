#!/usr/bin/env python3
"""
pnt_scaling.py — Does the PNT vs RH sensitivity gap grow with L?

From pnt_vs_rh.py:
  lsq=50: true ratio=0.599, random ratio=6.694  → gap = 11.2
  lsq=100: true ratio=0.263, random ratio=33.0   → gap = 125.5

If gap ~ L^α, fit α.

Also test: what is the OPTIMAL shift φ that maximizes sc (worst case shift)?
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
OUTFILE = '/Users/scones/projects/riemann/results/pnt_scaling.out'

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

def measure_sc(lsq, N, phase_offset=0.0):
    """Build tau with prime phases shifted by phase_offset × L / (2π)."""
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
        lkpp = np.array([float(mp.log(k)) for k, _ in VM])
        lkpp_mod = lkpp + phase_offset * L / (2 * np.pi)  # shift phases
        wts  = np.array([lk * k ** (-0.5) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp_mod) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp_mod / L)))
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
    return sc, lam0, dim, L

def main():
    # Measure true vs shifted at multiple lsq
    configs = [
        (20, 100),
        (25, 125),
        (35, 175),
        (50, 250),
        (75, 375),
        (100, 500),
    ]

    print(f"PNT sensitivity scaling: gap = ratio(π-shift) / ratio(true) vs L")
    print(f"{'lsq':>5} {'L':>7} {'sc_true':>12} {'sc_π':>12} {'ratio_true':>11} {'ratio_π':>11} {'gap':>8} {'dt':>5}")
    print("-" * 80)
    with open(OUTFILE, 'w') as f:
        f.write(f"pnt_scaling {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Sensitivity gap: ratio(π-shift)/ratio(true)\n\n")

    results = []
    for lsq, N in configs:
        t0 = time.time()
        sc_true, lam0_true, dim_true, L = measure_sc(lsq, N, 0.0)
        sc_pi, lam0_pi, dim_pi, _ = measure_sc(lsq, N, np.pi)
        dt = time.time() - t0
        sc_pred = 3.073e-3 / L ** 1.9151
        ratio_true = sc_true / sc_pred if not np.isnan(sc_true) else float('nan')
        ratio_pi   = sc_pi / sc_pred if not np.isnan(sc_pi) else float('nan')
        gap = ratio_pi / ratio_true if not (np.isnan(ratio_true) or ratio_true < 1e-10) else float('nan')
        line = f"{lsq:5d} {L:7.4f} {sc_true:12.4e} {sc_pi:12.4e} {ratio_true:11.3f} {ratio_pi:11.3f} {gap:8.2f}  ({dt:.0f}s)"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")
        results.append({'lsq': lsq, 'L': L, 'sc_true': sc_true, 'sc_pi': sc_pi,
                        'ratio_true': ratio_true, 'ratio_pi': ratio_pi, 'gap': gap})

    # Fit gap ~ L^α
    Ls = np.array([r['L'] for r in results if not np.isnan(r['gap'])])
    gaps = np.array([r['gap'] for r in results if not np.isnan(r['gap'])])
    if len(Ls) >= 3:
        logL = np.log(Ls); logG = np.log(gaps)
        A_fit = np.column_stack([logL, np.ones_like(logL)])
        coef, _, _, _ = np.linalg.lstsq(A_fit, logG, rcond=None)
        alpha, logC = coef
        C = np.exp(logC)
        print(f"\nGap power law: gap ~ {C:.2f} × L^{alpha:.2f}")
        with open(OUTFILE, 'a') as f:
            f.write(f"\ngap ~ {C:.2f} × L^{alpha:.2f}\n")

    json.dump(results, open(OUTFILE.replace('.out', '.json'), 'w'), indent=2, default=str)

if __name__ == "__main__":
    main()
