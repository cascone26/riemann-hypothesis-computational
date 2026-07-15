#!/usr/bin/env python3
"""
witness_highL.py — Measure sensitivity at larger L values.

From witness_scaling (lsq=20..75), sensitivity ~ L (linear).
Extend to lsq=150,200,300 to verify the scaling at larger L.
Also check: does sensitivity at N=5*lsq saturate or keep growing?

If sensitivity ~ C × L^alpha with alpha ≈ 1.0, then:
  |dlam_min/dΔσ| = C × log(lsq)

Combined with high-precision lam_min at c=13 (L=13) of +2.077e-59:
  Bound on first zero: Δσ₁ < lam_min / sensitivity(L=13) = 2.077e-59 / (C×13)

This is a NUMERICAL BOUND on zero location, not a proof.
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
GAMMA1 = 14.134725
OUTFILE = '/Users/scones/projects/riemann/results/witness_highL.out'

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

def sensitivity_at(lsq, N, dsigma=0.001):
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

    A = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
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

    tau_base = (A + B - P); tau_base = (tau_base + tau_base.T) / 2.0

    xi2 = grid / 2.0
    a1 = 0.25 + dsigma; b1 = xi2 - GAMMA1
    a2 = 0.25 - dsigma; b2 = xi2 + GAMMA1
    a0 = 0.25
    dOm = ((-a1)/(a1**2+b1**2) + (-a2)/(a2**2+b2**2)
           - (-a0)/(a0**2+b1**2) - (-a0)/(a0**2+b2**2))
    dA = ((E * (w * dOm)) @ E.conj().T).real / (2.0 * np.pi)
    tau_tilde = (tau_base + dA + (tau_base + dA).T) / 2.0

    lam0 = float(scipy.linalg.eigvalsh(tau_base, subset_by_index=[0,0])[0])
    lam1 = float(scipy.linalg.eigvalsh(tau_tilde, subset_by_index=[0,0])[0])
    return lam0, lam1, abs(lam1 - lam0) / dsigma, L

def main():
    # lsq=100 at both N=5*lsq and N=3*lsq to check N-dependence
    # lsq=150,200 at N=min(5*lsq, 800)
    configs = [
        (100, 300, "N=3*lsq"),
        (100, 500, "N=5*lsq"),
        (150, 600, "N=4*lsq (near-sat)"),
        (200, 600, "N=3*lsq"),
        (200, 800, "N=4*lsq"),
    ]

    print(f"Witness sensitivity at higher L")
    print(f"{'lsq':>5} {'N':>5} {'L':>7} {'lam0':>12} {'lam1':>12} {'sens':>12} {'dt':>5}")
    print("-" * 70)

    with open(OUTFILE, 'w') as f:
        f.write(f"witness_highL {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    results = []
    for lsq, N, label in configs:
        t0 = time.time()
        lam0, lam1, sens, L = sensitivity_at(lsq, N)
        dt = time.time() - t0
        ratio = sens / L
        line = f"{lsq:5d} {N:5d} {L:7.4f} {lam0:12.4e} {lam1:12.4e} {sens:12.4e}  ({dt:.0f}s)  sens/L={ratio:.4f}"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")
        results.append({'lsq': lsq, 'N': N, 'L': L, 'sens': sens, 'ratio': ratio})

    json.dump(results, open(OUTFILE.replace('.out', '.json'), 'w'), indent=2)

    # Summary: sens/L convergence
    print(f"\nKey: sens/L = {results[-1]['ratio']:.4f}  (converging to C=?)")
    print(f"High-precision bound at L=13: Δσ < 2.077e-59 / (C×13)")
    with open(OUTFILE, 'a') as f:
        ratios_str = str([round(r['ratio'], 4) for r in results])
        f.write(f"\nsens/L values: {ratios_str}\n")

if __name__ == "__main__":
    main()
