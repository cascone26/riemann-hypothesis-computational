#!/usr/bin/env python3
"""
zero_crossing_scan.py — Fine scan of lambda^2 across the first Riemann zero.

First zero of zeta: gamma_1 = 14.134725141734693...
Key question: does tau's spectrum "see" the zero? Do dim, lam_min, gap change at gamma_1?

Scans both fine (near gamma_1) and coarse (wide range).
"""
import sys, json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
N      = 64       # large enough to be reliable, fast enough to scan many values
THRESH_TIGHT = 1e-3
THRESH_WIDE  = 1e-2
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/zero_crossing_scan.json")

# First 20 nontrivial zeros of zeta (imaginary parts)
ZEROS = [
    14.134725141734693,
    21.022039638771555,
    25.010857580145688,
    30.424876125859513,
    32.935061587739189,
    37.586178158825671,
    40.918719012147495,
    43.327073280914999,
    48.005150881167159,
    49.773832477672302,
]


def vm_table(kmax):
    """Prime powers up to floor(kmax). Returns (prime_power, log(prime_base))."""
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


def build_tau(lambda_sq, N):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N, N + 1, dtype=float)
    D   = len(idx)
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
        log_k_pp = np.array([np.log(k) for k, _ in VM])
        wts      = np.array([lk * k ** (-0.5) for k, lk in VM])
        phases   = TWO_PI * np.outer(ns, log_k_pp) / L
        S        = np.sin(phases) @ wts
        dn       = ns[:, None] - ns[None, :]   # dn[i,j] = ns[i]-ns[j]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(phases) @ (2 * wts * (1 - log_k_pp / L)))
    else:
        P = np.zeros((D, D))
    tau = arch + B - P
    return (tau + tau.T) / 2.0


def analyze(lambda_sq, N):
    t0 = time.time()
    tau = build_tau(lambda_sq, N)
    ev  = scipy.linalg.eigvalsh(tau)
    dt  = time.time() - t0

    results = {}
    for label, thr in [('tight', THRESH_TIGHT), ('wide', THRESH_WIDE)]:
        mask = np.abs(ev) < thr
        dim  = int(mask.sum())
        lam_min = float(ev[0])
        above = ev[~mask & (ev > 0)]
        if dim > 0 and len(above) > 0:
            gap = float(above.min() - ev[mask].max())
        else:
            gap = float('nan')
        results[label] = {'dim': dim, 'lam_min': lam_min, 'gap': gap}
    results['t'] = round(dt, 2)

    # Distance to nearest zero
    nearest_zero = min(ZEROS, key=lambda z: abs(z - lambda_sq))
    dist_to_zero = lambda_sq - nearest_zero

    results['lambda_sq']     = lambda_sq
    results['L']             = float(mp.log(lambda_sq))
    results['n_pp']          = len(vm_table(lambda_sq))
    results['nearest_zero']  = nearest_zero
    results['dist_to_zero']  = dist_to_zero
    results['N']             = N
    return results


def main():
    # Scan plan:
    # 1. Fine scan around gamma_1 = 14.135
    # 2. Coarse scan out to gamma_5 = 32.9 to check first 5 zeros
    # 3. Wide scan to 100

    gamma1 = 14.134725141734693

    lambda_sqs = sorted(set([
        # Below gamma_1
        10.0, 12.0, 13.0, 13.5, 13.8, 14.0, 14.1,
        14.12, 14.13, 14.134, 14.1347, 14.1348,
        # Right at gamma_1
        gamma1 - 0.001, gamma1, gamma1 + 0.001,
        # Just above
        14.135, 14.14, 14.15, 14.2, 14.5,
        # Integer values
        15.0, 16.0, 18.0, 20.0,
        # Second zero region: gamma_2 = 21.02
        21.0, 21.022039638771555, 21.5, 22.0,
        # Third: gamma_3 = 25.01
        25.0, 25.01, 25.5,
        # Wide
        30.0, 36.0, 49.0, 64.0, 100.0,
    ]))

    results = []
    print(f"Zero-crossing scan | N={N} | G={G}")
    print(f"First zero: gamma_1 = {gamma1:.6f}")
    print()
    print(f"{'lam_sq':>10}  {'dist_zero':>10}  {'n_pp':>4}  {'dim_tight':>9}  "
          f"{'dim_wide':>9}  {'lam_min':>14}  {'gap_tight':>12}  {'t':>6}")

    for lsq in lambda_sqs:
        try:
            r = analyze(lsq, N)
            dt  = r['tight']['dim']
            dw  = r['wide']['dim']
            lm  = r['tight']['lam_min']
            gap = r['tight']['gap']
            dist = r['dist_to_zero']
            print(f"{lsq:>10.6f}  {dist:>+10.4f}  {r['n_pp']:>4}  "
                  f"{dt:>9}  {dw:>9}  {lm:>14.4e}  {gap:>12.4e}  {r['t']:>6.2f}s")
            sys.stdout.flush()
            results.append(r)
        except Exception as e:
            print(f"lam_sq={lsq:.4f} ERROR: {e}")
            sys.stdout.flush()

    with open(OUTFILE, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
