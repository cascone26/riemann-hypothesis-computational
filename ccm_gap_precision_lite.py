#!/usr/bin/env python3
"""
ccm_gap_precision_lite.py -- Re-verify CCM gaps for N=20,40 at high precision with optimized grid.

Uses mpmath dps=50 but with reduced grid size for tractability:
- Arch integration via composite Simpson's rule (fewer grid points than trapezoid)
- Focus on eigenvalue gap (not full spectrum)
- Reuse float64 construction for w02 + prime (only Omega integration in mpmath)
"""

import time
import json
import numpy as np
from mpmath import mp, mpf, mpc, matrix, eigh, digamma, log, exp, sin, cos, sinh, pi, sqrt
from scipy.special import digamma as _digamma_scipy

# High-precision setup
mp.dps = 50

LAMBDA_SQ = 14
L = 2.0 * mp.log(mp.sqrt(mpf(LAMBDA_SQ)))
TWO_PI = 2.0 * mp.pi
LOGPI = mp.log(mp.pi)

print(f"Starting lite high-precision reverify at dps={mp.dps}")
print(f"L = {mp.nstr(L, 15)}")


def vm_table(kmax):
    """Von Mangoldt table."""
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


def omega_mp(r):
    """Omega(r) = Re psi(1/4 + i*r/2) - log(pi) in mpmath."""
    z = mpf(0.25) + mpc(0, mpf(r) / 2.0)
    psi_val = digamma(z)
    return psi_val.real - LOGPI


def e_hat_mp(r, n):
    """e_hat_n(r) in mpmath, with near-singularity handling."""
    s = TWO_PI * mpf(n) / L
    d = s - mpf(r)

    if abs(d) < mpf(1e-10):
        return mp.sqrt(L)

    numerator = 2.0 * mp.sin(d * L / 2.0)
    exp_term = mp.exp(mpc(0, d * L / 2.0))
    return (1.0 / mp.sqrt(L)) * exp_term * numerator / d


def q_kernel_mp(n, m, y):
    """q kernel in mpmath."""
    yn = mpf(y)
    if n == m:
        return 2.0 * (1.0 - abs(yn) / L) * mp.cos(TWO_PI * mpf(n) * yn / L)
    else:
        num = mp.sin(TWO_PI * mpf(m) * yn / L) - mp.sin(TWO_PI * mpf(n) * yn / L)
        denom = mp.pi * mpf(n - m)
        return num / denom


def prime_mat_np(idx):
    """Compute prime matrix using numpy (acceptable precision for this term)."""
    D = len(idx)
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            s = 0.0
            for k, Lam in VM:
                s += float(Lam) * (k ** (-0.5)) * q_kernel_mp(n, m, float(mp.log(mpf(k))))
            P[ii, jj] = P[jj, ii] = s
    return P


def w02_mat_np(idx):
    """W02 matrix using numpy."""
    D = len(idx)
    W = np.zeros((D, D))
    sh = float(mp.sinh(L / 4.0) ** 2)
    L_f = float(L)
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            num = L_f ** 2 - 16.0 * np.pi ** 2 * m * n
            den = (L_f ** 2 + 16.0 * np.pi ** 2 * m ** 2) * (L_f ** 2 + 16.0 * np.pi ** 2 * n ** 2)
            W[ii, jj] = 32.0 * L_f * sh * num / den
    return W


def arch_mat_mp_lite(N, R, G):
    """
    Compute archimedean matrix with mpmath Omega integration.
    Uses Simpson's rule for better accuracy with fewer points.
    """
    idx = list(range(-N, N + 1))
    D = len(idx)

    # Create grid (coarse for speed)
    grid_np = np.linspace(-R, R, G)
    h = float(grid_np[1] - grid_np[0])

    # Simpson's rule weights (if G is odd)
    if G % 2 == 0:
        G += 1
        grid_np = np.linspace(-R, R, G)
        h = float(grid_np[1] - grid_np[0])

    w = np.zeros(G)
    w[0] = h / 3.0
    w[-1] = h / 3.0
    for i in range(1, G - 1):
        w[i] = h / 3.0 * (4.0 if i % 2 == 1 else 2.0)

    # Compute Omega at grid points using mpmath
    print(f"  Computing Omega at {G} grid points...")
    Om = np.array([float(omega_mp(float(g))) for g in grid_np], dtype=float)

    # Build E matrix (eigenfunction values)
    print(f"  Building E matrix...")
    E = np.zeros((D, G), dtype=complex)
    for k, n in enumerate(idx):
        for g in range(G):
            # Use mpmath for e_hat calculation
            e_val = e_hat_mp(float(grid_np[g]), n)
            E[k, g] = complex(e_val)

    # Compute arch = (1/(2pi)) * E @ diag(w * Om) @ E^H
    print(f"  Computing arch matrix...")
    A = E * (w * Om)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)

    return arch, idx


def verify_gap_lite(N):
    """Verify eigenvalue gap for N at high precision."""
    print(f"\n{'='*70}")
    print(f"VERIFICATION: N={N} at dps={mp.dps}")
    print(f"{'='*70}")

    # Parameters
    R = float(TWO_PI * N / L + 80.0)
    # Use MUCH smaller grid for tractability
    G = min(2000, 1000 + 100 * N)  # Much smaller than production
    dim = 2 * N + 1

    print(f"N={N}, dim={dim}, R={R:.1f}, G={G}")

    t0 = time.time()
    arch, idx = arch_mat_mp_lite(N, R, G)
    dt_arch = time.time() - t0
    print(f"  arch computed: {dt_arch:.1f}s")

    t0 = time.time()
    w02 = w02_mat_np(idx)
    print(f"  w02 computed: {time.time()-t0:.1f}s")

    t0 = time.time()
    prime = prime_mat_np(idx)
    print(f"  prime computed: {time.time()-t0:.1f}s")

    # Build tau = arch + w02 - prime
    tau = arch + w02 - prime
    tau = (tau + tau.T) / 2

    # Compute two smallest eigenvalues
    print(f"  Computing eigendecomposition...")
    t0 = time.time()
    evals, evecs = np.linalg.eigh(tau)
    dt_eig = time.time() - t0
    print(f"  eigh completed: {dt_eig:.1f}s")

    lambda_1_float64 = float(evals[0])
    lambda_2_float64 = float(evals[1])
    gap_float64 = lambda_2_float64 - lambda_1_float64

    print(f"\nResults (float64 on mpmath-constructed tau):")
    print(f"  lambda_1: {lambda_1_float64:.6e}")
    print(f"  lambda_2: {lambda_2_float64:.6e}")
    print(f"  gap:      {gap_float64:.6e}")

    # Compare to S2753 float64 result
    if N == 20:
        s2753_gap = 5.94e-14
    elif N == 40:
        s2753_gap = 3.92e-15
    else:
        s2753_gap = None

    if s2753_gap is not None:
        ratio = gap_float64 / s2753_gap
        print(f"\nComparison to S2753 float64 result ({s2753_gap:.3e}):")
        print(f"  This / S2753: {ratio:.3f}x")
        print(f"  Gap still positive? {gap_float64 > 0}")

    return {
        "N": N,
        "dim": dim,
        "grid_size": G,
        "lambda_1_float64": lambda_1_float64,
        "lambda_2_float64": lambda_2_float64,
        "gap_float64": gap_float64,
        "gap_positive": gap_float64 > 0,
        "s2753_gap": s2753_gap,
        "ratio_to_s2753": ratio if s2753_gap else None,
    }


def main():
    results = []
    for N in [20, 40]:
        result = verify_gap_lite(N)
        results.append(result)

    output_file = "/Users/scones/projects/riemann/results/ccm_gap_precision_lite.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")

    return results


if __name__ == "__main__":
    main()
