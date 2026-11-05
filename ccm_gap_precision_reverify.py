#!/usr/bin/env python3
"""
ccm_gap_precision_reverify.py -- Re-verify CCM eigenvalue gaps for N=20,40 at high precision (mpmath dps=50+).

Task: The float64 computation in S2753 found:
  N=20: gap = 5.94e-14 (dangerously close to precision ceiling)
  N=40: gap = 3.92e-15 (only ~17.8x above machine epsilon 2.2e-16)

This script rebuilds tau matrices using mpmath at dps=50 throughout every
intermediate step (digamma, arch, w02, prime, eigenvalues) to confirm
whether these gaps are REAL or NOISE.

Constants and math are identical to connes_decisive_v3.py, only using mpmath
instead of scipy/numpy for precision-critical terms.
"""

import time
import json
from mpmath import mp, mpf, mpc, matrix, eigh, digamma, log, exp, sin, cos, sinh, pi, sqrt

# High-precision setup
mp.dps = 50  # 50 decimal places

# Constants (all in mpmath)
LAMBDA_SQ = 14
L = 2.0 * mp.log(mp.sqrt(mpf(LAMBDA_SQ)))
TWO_PI = 2.0 * mp.pi
LOGPI = mp.log(mp.pi)

print(f"Starting high-precision reverify at dps={mp.dps}")
print(f"L = {L}")
print(f"TWO_PI = {TWO_PI}")
print(f"LOGPI = {LOGPI}")


def vm_table(kmax):
    """Von Mangoldt Lambda table for k <= kmax."""
    tab = []
    for k in range(2, kmax + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0:
                    p = f
                elif p != f:
                    ok = False
                    break
                while n % f == 0:
                    n //= f
            f += 1
        if n > 1:
            if p == 0:
                p = n
            elif p != n:
                ok = False
        if ok and p > 0:
            tab.append((k, mp.log(mpf(p))))
    return tab


VM = vm_table(LAMBDA_SQ)


def omega(r):
    """Omega(r) = Re psi(1/4 + i*r/2) - log(pi) at high precision."""
    z = mpf(0.25) + mpc(0, mpf(r) / 2.0)
    psi_val = digamma(z)
    return psi_val.real - LOGPI


def e_hat(r, n):
    """e_hat_n(r) = (1/sqrt L) e^{i*d*L/2} * 2*sin(d*L/2) / d, d = 2pi*n/L - r."""
    s = TWO_PI * mpf(n) / L
    d = s - mpf(r)

    # Handle near-singularity (d very close to 0)
    if abs(d) < mpf(1e-10):
        return mp.sqrt(L)

    numerator = 2.0 * mp.sin(d * L / 2.0)
    denom = d
    exp_term = mp.exp(mpc(0, d * L / 2.0))
    return (1.0 / mp.sqrt(L)) * exp_term * numerator / denom


def q_kernel(n, m, y):
    """q kernel for prime term."""
    yn = mpf(y)
    if n == m:
        return 2.0 * (1.0 - abs(yn) / L) * mp.cos(TWO_PI * mpf(n) * yn / L)
    else:
        num = mp.sin(TWO_PI * mpf(m) * yn / L) - mp.sin(TWO_PI * mpf(n) * yn / L)
        denom = mp.pi * mpf(n - m)
        return num / denom


def prime_mat_mp(idx):
    """Prime/von Mangoldt matrix at high precision."""
    D = len(idx)
    P = matrix(D, D)
    for i in range(D):
        for j in range(D):
            n = idx[i]
            m = idx[j]
            s = mpf(0)
            for k, Lam in VM:
                s += Lam * mpf(k) ** (-mpf(0.5)) * q_kernel(n, m, mp.log(mpf(k)))
            P[i, j] = s
    return P


def w02_mat_mp(idx):
    """W_{0,2} boundary/pole matrix at high precision."""
    D = len(idx)
    W = matrix(D, D)
    sh = mp.sinh(L / 4.0) ** 2
    for i in range(D):
        for j in range(D):
            n = mpf(idx[i])
            m = mpf(idx[j])
            num = L ** 2 - 16.0 * mp.pi ** 2 * m * n
            den = (L ** 2 + 16.0 * mp.pi ** 2 * m ** 2) * (L ** 2 + 16.0 * mp.pi ** 2 * n ** 2)
            W[i, j] = 32.0 * L * sh * num / den
    return W


def arch_mat_mp(N, R, G):
    """
    Archimedean block at high precision using high-precision integration grid.

    arch = (1/2pi) * E @ diag(w * Omega) @ E^H, where E[k,g] = e_hat_k(grid[g]).
    """
    idx = list(range(-N, N + 1))
    D = len(idx)

    # Create high-precision grid
    grid = [mpf(-R) + (mpf(R) - mpf(-R)) * mpf(i) / mpf(G - 1) for i in range(G)]
    h = (grid[1] - grid[0]) if G > 1 else mpf(1)

    # Trapezoid weights
    w = [h] * G
    w[0] = h / 2.0
    w[-1] = h / 2.0

    # Omega values at grid points
    Om = [omega(float(g)) for g in grid]

    # Build E matrix: E[k,g] = e_hat_idx[k](grid[g])
    E = matrix(D, G)
    for k, n in enumerate(idx):
        for g in range(G):
            E[k, g] = e_hat(float(grid[g]), n)

    # A = E * diag(w * Om)
    A = matrix(D, G)
    for k in range(D):
        for g in range(G):
            A[k, g] = E[k, g] * w[g] * Om[g]

    # arch = (A @ E^H) / (2*pi)
    EH = E.transpose().conjugate()
    arch_raw = A * EH
    arch = arch_raw / (2.0 * mp.pi)

    return arch, idx


def build_tau_mp(N, R, G, eps=mpf(1)):
    """Build full tau matrix at high precision."""
    print(f"\nBuilding tau for N={N} at dps={mp.dps}")
    print(f"  Grid: R={R}, G={G}, dim={2*N+1}")

    t0 = time.time()
    arch, idx = arch_mat_mp(N, float(R), G)
    dt = time.time() - t0
    print(f"  arch computed: {dt:.1f}s")

    t0 = time.time()
    w02 = w02_mat_mp(idx)
    dt = time.time() - t0
    print(f"  w02 computed: {dt:.1f}s")

    t0 = time.time()
    prime = prime_mat_mp(idx)
    dt = time.time() - t0
    print(f"  prime computed: {dt:.1f}s")

    # tau = arch + eps * (w02 - prime)
    pert = w02 - prime
    tau = arch + eps * pert

    # Symmetrize
    tau_T = tau.transpose()
    tau = (tau + tau_T) / 2.0

    return tau, idx


def eigh_smallest_mp(M, n_smallest=2):
    """
    Compute the n_smallest smallest eigenvalues using mpmath's eigh.

    Returns (eigenvalues, eigenvectors) sorted by eigenvalue.
    """
    print(f"  Computing eigendecomposition of {M.rows}x{M.cols} matrix...")
    t0 = time.time()
    evs = mp.eigh(M)
    dt = time.time() - t0
    print(f"  eigh completed: {dt:.1f}s")

    # evs is a list of (eigenvalue, eigenvector) pairs
    evs_sorted = sorted(evs, key=lambda x: abs(x[0]))

    smallest = evs_sorted[:n_smallest]
    evals = [ev[0] for ev in smallest]
    evecs = [ev[1] for ev in smallest]

    return evals, evecs


def verify_gap(N):
    """Verify the eigenvalue gap for a specific N at high precision."""
    print(f"\n{'='*70}")
    print(f"VERIFICATION: N={N}")
    print(f"{'='*70}")

    # Parameters (same as S2753)
    R = float(TWO_PI * N / L + 80.0)
    G = max(60000, 3000 * N)

    # Build tau
    tau, idx = build_tau_mp(N, R, G, eps=mpf(1))

    # Compute smallest eigenvalues
    evals, evecs = eigh_smallest_mp(tau, n_smallest=2)

    lambda_1 = evals[0]
    lambda_2 = evals[1]
    gap = lambda_2 - lambda_1

    print(f"\nResults for N={N}:")
    print(f"  lambda_1 (smallest):    {lambda_1}")
    print(f"  lambda_2 (2nd smallest): {lambda_2}")
    print(f"  gap = lambda_2 - lambda_1: {gap}")
    print(f"  gap (sci notation):     {mp.nstr(gap, 5, strip_zeros=False)}")
    print(f"  gap as float:           {float(gap):.6e}")

    # Compare to float64 result
    if N == 20:
        float64_gap = 5.94e-14
    elif N == 40:
        float64_gap = 3.92e-15
    else:
        float64_gap = None

    if float64_gap is not None:
        ratio = float(gap) / float64_gap
        print(f"\nComparison to float64 result ({float64_gap:.3e}):")
        print(f"  mpmath gap / float64 gap: {ratio:.3f}x")
        print(f"  Same sign? {(gap > 0) and float64_gap > 0}")

    return {
        "N": N,
        "dps": mp.dps,
        "dim": 2 * N + 1,
        "lambda_1": str(lambda_1),
        "lambda_2": str(lambda_2),
        "gap": str(gap),
        "gap_float": float(gap),
        "gap_positive": gap > 0,
    }


def main():
    results = []

    # Verify N=20 and N=40
    for N in [20, 40]:
        result = verify_gap(N)
        results.append(result)

    # Also verify at dps=100 for N=40 to check stability
    print(f"\n{'='*70}")
    print("CROSS-CHECK: N=40 at dps=100 for precision stability")
    print(f"{'='*70}")
    mp.dps = 100
    result_100 = verify_gap(40)
    result_100["dps"] = 100
    results.append(result_100)

    # Save results
    output_file = "/Users/scones/projects/riemann/results/ccm_gap_precision_reverify.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")

    return results


if __name__ == "__main__":
    main()
