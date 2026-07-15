#!/usr/bin/env python3
"""
witness_scaling.py — Measure |dlam_min/dΔσ| vs L to get sensitivity scaling.

From witness_exact:
  lsq=50 (L=3.91): |dlam/dΔσ| ≈ 3.91
  lsq=100 (L=4.61): |dlam/dΔσ| ≈ 4.68

Ratio = 1.20, L ratio = 1.18 → roughly L^1.

Extend to lsq=25,50,100,200 at N=5*lsq (safe saturation).
Use Δσ=0.001 (small enough for linear regime).

FIT: |dlam/dΔσ| ~ C × L^α
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
GAMMA1 = 14.134725  # First Riemann zero
OUTFILE = '/Users/scones/projects/riemann/results/witness_scaling.out'

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

def build_tau_and_delta(lsq, N, dsigma):
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

    # Exact delta_Omega for moving gamma1 off-line by dsigma
    xi2 = grid / 2.0
    a1 = 0.25 + dsigma; b1 = xi2 - GAMMA1
    a2 = 0.25 - dsigma; b2 = xi2 + GAMMA1
    a0 = 0.25
    dOm = ((-a1)/(a1**2+b1**2) + (-a2)/(a2**2+b2**2)
           - (-a0)/(a0**2+b1**2) - (-a0)/(a0**2+b2**2))

    dA = ((E * (w * dOm)) @ E.conj().T).real / (2.0 * np.pi)
    tau_tilde = tau_base + dA; tau_tilde = (tau_tilde + tau_tilde.T) / 2.0

    lam_base = float(scipy.linalg.eigvalsh(tau_base, subset_by_index=[0,0])[0])
    lam_tilde = float(scipy.linalg.eigvalsh(tau_tilde, subset_by_index=[0,0])[0])
    return lam_base, lam_tilde, L

def main():
    # Test cases: small N (2*lsq) for speed check, then N=min(5*lsq,800)
    configs = [
        (20, 100, "N=5*lsq"),
        (25, 125, "N=5*lsq"),
        (35, 175, "N=5*lsq"),
        (50, 250, "N=5*lsq"),
        (75, 375, "N=5*lsq"),
        (100, 500, "N=5*lsq"),
        (150, 750, "N=5*lsq"),
        (200, 800, "N=4*lsq (Mac limit)"),
    ]
    dsigma = 0.001  # Small: linear regime

    print(f"Witness sensitivity scaling: |dlam/dΔσ| vs L")
    print(f"Δσ = {dsigma}")
    print(f"{'lsq':>5} {'L':>8} {'lam_base':>12} {'lam_tilde':>12} {'|Δlam/Δσ|':>14} {'dt':>5}")
    print("-" * 65)
    with open(OUTFILE, 'w') as f:
        f.write(f"witness_scaling {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Δσ={dsigma}\n")
        f.write(f"{'lsq':>5} {'L':>8} {'lam_base':>12} {'lam_tilde':>12} {'|Δlam/Δσ|':>14}\n")
        f.write("-" * 55 + "\n")

    results = []
    for lsq, N, label in configs:
        t0 = time.time()
        lam_base, lam_tilde, L = build_tau_and_delta(lsq, N, dsigma)
        dt = time.time() - t0
        sensitivity = abs(lam_tilde - lam_base) / dsigma
        line = f"{lsq:5d} {L:8.4f} {lam_base:12.4e} {lam_tilde:12.4e} {sensitivity:14.4e}  ({dt:.0f}s)"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(f"{lsq:5d} {L:8.4f} {lam_base:12.4e} {lam_tilde:12.4e} {sensitivity:14.4e}\n")
        results.append({'lsq': lsq, 'L': L, 'sens': sensitivity})

    # Fit power law |sens| ~ C × L^α
    Ls = np.array([r['L'] for r in results])
    Ss = np.array([r['sens'] for r in results])
    logL = np.log(Ls)
    logS = np.log(Ss)
    A_fit = np.column_stack([logL, np.ones_like(logL)])
    coef, _, _, _ = np.linalg.lstsq(A_fit, logS, rcond=None)
    alpha, logC = coef
    C = np.exp(logC)
    print(f"\nPower law fit: |dlam/dΔσ| ~ {C:.4f} × L^{alpha:.4f}")
    print(f"α = {alpha:.4f}  (expect ~1.0 if sensitivity ~ L)")
    with open(OUTFILE, 'a') as f:
        f.write(f"\nPower law: |sens| ~ {C:.4f} × L^{alpha:.4f}\n")
        f.write(f"α = {alpha:.4f}\n")

    json.dump(results, open(OUTFILE.replace('.out', '.json'), 'w'), indent=2)
    print(f"Output: {OUTFILE}")

if __name__ == "__main__":
    main()
