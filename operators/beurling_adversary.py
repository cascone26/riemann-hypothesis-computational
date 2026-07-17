#!/usr/bin/env python3
"""
beurling_adversary.py — PRISM item #1: Beurling adversary diff test.

Tests whether the CvS tau_N operator DIFFERENTIATES between:
  (A) Real prime comb (Riemann zeta → zeros on critical line)
  (B) Beurling prime combs where RH is KNOWN to fail

If the operator is genuinely selective (not circular), it should:
  - Stay PSD for real primes
  - Go non-PSD (negative eigenvalue) or have significantly larger min_ev for Beurling combs

Three Beurling systems tested:
  B1: {2} only — ζ_B(s)=1/(1-2^{-s}), zeros at s=2πin/log2 (on σ=0, not σ=1/2)
  B2: {2, 3} only — sparse, missing all primes ≥ 5
  B3: Real primes + fake prime at 4.1 — known not to be a prime power
  B4: Evenly-spaced "primes" at 2,3,5,7,11... shifted by 0.5 — breaks PNT
"""

import numpy as np
import scipy.linalg
from scipy.special import digamma as _digamma
import mpmath as mp
import time

mp.mp.dps = 30
LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
G = 30000  # arch quadrature points

def omega_vec(r):
    return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI

def e_hat_vec(r, n):
    s = TWO_PI * n / L
    d = s - r
    out = np.empty_like(d, dtype=complex)
    small = np.abs(d) < 1e-12
    out[small] = np.sqrt(L) + 0j
    dd = d[~small]
    out[~small] = (1/np.sqrt(L)) * np.exp(0.5j*dd*L) * 2*np.sin(0.5*dd*L) / dd
    return out

def build_arch(idx, G=G):
    D = len(idx)
    R = TWO_PI * max(abs(n) for n in idx) / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5*h
    Om = omega_vec(grid)
    E = np.array([e_hat_vec(grid, n) for n in idx])
    A = E * (w * Om)
    return (A @ E.conj().T).real / (2*np.pi)

def build_B(idx):
    D = len(idx)
    W = np.zeros((D,D))
    sh = np.sinh(L/4)**2
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            num = L**2 - 16*np.pi**2 * m * n
            den = (L**2 + 16*np.pi**2*m**2) * (L**2 + 16*np.pi**2*n**2)
            W[ii,jj] = 32*L*sh*num/den
    return W

def q_kernel(n, m, y):
    if n == m:
        return 2*(1 - abs(y)/L) * np.cos(TWO_PI*n*y/L)
    return (np.sin(TWO_PI*m*y/L) - np.sin(TWO_PI*n*y/L)) / (np.pi*(n-m))

def build_P(idx, comb):
    """comb: list of (k_real_or_float, Lambda_value) pairs."""
    D = len(idx)
    P = np.zeros((D,D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            s = sum(Lam * k**(-0.5) * q_kernel(n, m, np.log(k)) for k, Lam in comb)
            P[ii,jj] = P[jj,ii] = s
    return P

def build_tau(N, comb):
    idx = list(range(-N, N+1))
    arch = build_arch(idx)
    B = build_B(idx)
    P = build_P(idx, comb)
    tau = arch + B - P
    return (tau + tau.T)/2

def min_ev(tau):
    evs = scipy.linalg.eigvalsh(tau, subset_by_index=[0,1])
    return min(abs(evs[0]), abs(evs[1])), evs[0]

# ── Real prime comb ───────────────────────────────────────────────────────────
def real_comb():
    tab = []
    for k in range(2, LAMBDA_SQ+1):
        n, p, f, ok = k, 0, 2, True
        while f*f <= n:
            if n%f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n%f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((float(k), float(mp.log(p))))
    return tab

# ── Beurling combs ────────────────────────────────────────────────────────────
def beurling_powers_of_2():
    """B1: only prime is 2. ζ_B(s) = 1/(1-2^{-s}), zeros at σ=0."""
    return [(float(2**j), float(np.log(2))) for j in range(1, int(np.log2(LAMBDA_SQ))+1)]

def beurling_2_3():
    """B2: primes = {2,3} only."""
    tab = []
    for p in [2, 3]:
        j = 1
        while p**j <= LAMBDA_SQ:
            tab.append((float(p**j), float(np.log(p))))
            j += 1
    return tab

def beurling_plus_fake():
    """B3: Real primes + fake prime power at 4.1 (not a prime power)."""
    tab = real_comb()
    k_fake = 4.1
    if k_fake <= LAMBDA_SQ:
        tab.append((k_fake, float(np.log(k_fake))))
    return sorted(tab)

def beurling_shifted():
    """B4: Real prime powers shifted by +0.5 (none are integers → breaks PNT)."""
    return [(float(k) + 0.5, lam) for k, lam in real_comb() if float(k)+0.5 <= LAMBDA_SQ]


if __name__ == "__main__":
    N = 8  # small but enough to see sign differences
    print(f"CvS tau_N Beurling adversary test  N={N}  LAMBDA_SQ={LAMBDA_SQ}  L={L:.4f}")
    print(f"{'Comb':<30} {'min_ev':>12} {'ev[0]':>14} {'PSD?':>6}  {'vs real':>10}  time")
    print("-"*90)

    combs = [
        ("Real primes (baseline)",   real_comb()),
        ("B1: {2} only (σ=0 zeros)", beurling_powers_of_2()),
        ("B2: {2,3} only",           beurling_2_3()),
        ("B3: real + fake@4.1",      beurling_plus_fake()),
        ("B4: shifted primes +0.5",  beurling_shifted()),
    ]

    real_min = None
    results = []
    for name, comb in combs:
        t0 = time.time()
        tau = build_tau(N, comb)
        me, ev0 = min_ev(tau)
        dt = time.time() - t0
        psd = "YES" if ev0 >= -1e-12 else f"NO({ev0:.2e})"
        if real_min is None:
            real_min = me
            ratio = "  baseline"
        else:
            ratio = f"{me/real_min:>8.2f}×"
        print(f"{name:<30} {me:>12.4e} {ev0:>14.4e} {psd:>6}  {ratio}  ({dt:.1f}s)")
        results.append({"name": name, "min_ev": me, "ev0": ev0, "psd": ev0 >= -1e-12})

    print("\nKey question: do any Beurling combs produce ev0 < 0 (non-PSD)?")
    non_psd = [r for r in results if not r["psd"]]
    if non_psd:
        print(f"  YES — {len(non_psd)} non-PSD cases: {[r['name'] for r in non_psd]}")
    else:
        print("  All PSD at N=8. Try larger N or check if discrimination is via min_ev RATIO.")
    print(f"\nReal prime min_ev: {real_min:.4e}")
    for r in results[1:]:
        print(f"  {r['name']}: ratio = {r['min_ev']/real_min:.3f}× ({'LARGER' if r['min_ev']>real_min else 'SMALLER'})")
