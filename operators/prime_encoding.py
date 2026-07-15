#!/usr/bin/env python3
"""
prime_encoding.py — What in the primes encodes RH?

From pnt_vs_rh.py: shifting all prime phases by π gives sc 36× larger.
This means the specific prime phases (log(p_k) mod 2π/L) are critical.

QUESTION: Is the critical structure:
  (a) The GAPS between consecutive primes (prime gap distribution)?
  (b) The exact prime VALUES (log spacing)?
  (c) The PRIME ORBIT structure (sum over prime powers)?

TEST APPROACH:
Build 5 variants of the prime sum P:
1. true: actual primes
2. permuted: same primes, random permutation of weights
3. scaled: all logs scaled by constant α (shifts all zeros' contributions)
4. smoothed: replace primes with their linear-interpolated positions
5. poisson: prime-like Poisson process with same density

If 'permuted' ≈ 'true': order doesn't matter, only individual values
If 'scaled' ≈ 'true': multiplicative structure preserved
If 'smoothed' >> 'true': prime clumping is essential
If 'poisson' >> 'true': true randomness destroys cancellation
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
OUTFILE = '/Users/scones/projects/riemann/results/prime_encoding.out'

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

def build_tau_with_prime_table(lsq, N, log_k_vals, weights):
    """Build tau with custom prime table (log_k_vals, weights)."""
    L  = float(mp.log(lsq))
    sh = float(np.sinh(L / 4) ** 2)
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

    if len(log_k_vals) > 0:
        lkpp = np.array(log_k_vals)
        wts  = np.array(weights)
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
    lsq, N = 50, 250
    L = float(mp.log(lsq))
    VM = vm_table(lsq)
    sc_pred = 3.073e-3 / L ** 1.9151

    # True prime table: phases use log(prime_power)=log(k), weights use log(prime_base)=log_k
    lk_true = np.array([float(mp.log(k)) for k, _ in VM])   # log(prime_power) for phases
    lk_base = np.array([log_k for _, log_k in VM])           # log(prime_base) for weights
    wts_true = np.array([log_k * k ** (-0.5) for k, log_k in VM])

    print(f"Prime encoding test: lsq={lsq}, N={N}, {len(VM)} prime powers")
    print(f"pred sc = {sc_pred:.4e}")
    print()

    with open(OUTFILE, 'w') as f:
        f.write(f"prime_encoding {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    variants = []

    # 1. True primes
    variants.append(('true', lk_true, wts_true))

    # 2. Permuted weights (same log values, randomized weights)
    rng = np.random.RandomState(42)
    wts_perm = wts_true[rng.permutation(len(wts_true))]
    variants.append(('wt_perm', lk_true, wts_perm))

    # 3. Permuted log positions (same weights, randomized log values within range)
    lk_perm = lk_true[rng.permutation(len(lk_true))]
    variants.append(('lk_perm', lk_perm, wts_true))

    # 4. Scaled logs by √2 (stretches all positions by constant factor)
    lk_scaled = lk_true * np.sqrt(2)
    variants.append(('lk_√2', lk_scaled, wts_true))

    # 5. Smooth approximation: replace log(k) with smooth function log(k) → continuous spacing
    # Space primes uniformly in log-space
    n_primes = len(VM)
    lk_smooth = np.linspace(lk_true[0], lk_true[-1], n_primes)
    wts_smooth = np.array([wts_true[i] for i in range(n_primes)])
    variants.append(('smooth', lk_smooth, wts_smooth))

    # 6. Poisson: prime positions replaced by Poisson process with same density
    mean_gap = (lk_true[-1] - lk_true[0]) / (n_primes - 1)
    lk_poisson_raw = np.sort(rng.exponential(mean_gap, size=n_primes * 3).cumsum() + lk_true[0])
    lk_poisson = lk_poisson_raw[:n_primes]  # take exactly n_primes points
    variants.append(('poisson', lk_poisson, wts_smooth))

    # 7. GUE-spaced: prime positions follow GUE level spacing statistics
    u = rng.uniform(0, 1, size=n_primes)
    gue_gaps = np.sqrt(-4 / np.pi * np.log(1 - u + 1e-10)) * mean_gap
    lk_gue = np.cumsum(gue_gaps) + lk_true[0]
    lk_gue = lk_gue[:n_primes]
    variants.append(('GUE', lk_gue, wts_smooth))

    print(f"{'variant':12s} {'sc':>12} {'ratio':>8} {'lam0':>12} {'dt':>5}")
    print("-" * 55)
    with open(OUTFILE, 'a') as f:
        f.write(f"{'variant':12s} {'sc':>12} {'ratio':>8} {'lam0':>12}\n")
        f.write("-" * 48 + "\n")

    for name, lk, wts in variants:
        t0 = time.time()
        sc, lam0, dim = build_tau_with_prime_table(lsq, N, lk[:n_primes], wts[:n_primes])
        dt = time.time() - t0
        ratio = sc / sc_pred if not np.isnan(sc) else float('nan')
        line = f"{name:12s} {sc:12.4e} {ratio:8.3f} {lam0:12.4e}  ({dt:.0f}s)"
        print(line)
        with open(OUTFILE, 'a') as f:
            f.write(f"{name:12s} {sc:12.4e} {ratio:8.3f} {lam0:12.4e}\n")

    print()
    print("INTERPRETATION:")
    print("  true: baseline (ratio~0.6)")
    print("  wt_perm: permuting weights only")
    print("  lk_perm: permuting log positions only")
    print("  lk_√2: stretching all positions by sqrt(2)")
    print("  smooth: smooth log-spacing (no prime clumping)")
    print("  poisson: random prime-like positions")
    print("  GUE: positions with GUE level statistics")
    print("  Small ratio = cancellation preserved  |  Large ratio = cancellation destroyed")

if __name__ == "__main__":
    main()
