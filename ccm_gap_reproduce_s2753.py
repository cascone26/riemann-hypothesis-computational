#!/usr/bin/env python3
"""
Reproduce S2753 results for N=20 and N=40 exactly.

Verify the reported gaps and compute residuals/condition numbers to assess reliability.
"""

import numpy as np
from scipy.special import digamma as _digamma
import time

# Constants from S2753
LAMBDA_SQ = 14
L = float(2 * np.log(np.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))

print(f"Reproducing S2753 results")
print(f"L = {L}")
print(f"TWO_PI = {TWO_PI}")

def vm_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
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
            tab.append((k, float(np.log(p))))
    return tab

VM = vm_table(LAMBDA_SQ)

def omega_vec(r):
    return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI

def e_hat_vec(r, n):
    s = TWO_PI * n / L
    d = s - r
    out = np.empty_like(d, dtype=complex)
    small = np.abs(d) < 1e-12
    out[small] = np.sqrt(L) + 0j
    dd = d[~small]
    out[~small] = (1.0 / np.sqrt(L)) * np.exp(0.5j * dd * L) * 2.0 * np.sin(0.5 * dd * L) / dd
    return out

def q_kernel(n, m, y):
    if n == m:
        return 2.0 * (1.0 - abs(y) / L) * np.cos(TWO_PI * n * y / L)
    return (np.sin(TWO_PI * m * y / L) - np.sin(TWO_PI * n * y / L)) / (np.pi * (n - m))

def prime_mat(idx):
    D = len(idx)
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            s = 0.0
            for k, Lam in VM:
                s += Lam * k ** (-0.5) * q_kernel(n, m, np.log(k))
            P[ii, jj] = P[jj, ii] = s
    return P

def w02_mat(idx):
    D = len(idx)
    W = np.zeros((D, D))
    sh = np.sinh(L / 4) ** 2
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            num = (L ** 2 - 16 * np.pi ** 2 * m * n)
            den = (L ** 2 + 16 * np.pi ** 2 * m ** 2) * (L ** 2 + 16 * np.pi ** 2 * n ** 2)
            W[ii, jj] = 32 * L * sh * num / den
    return W

def build_tau(N):
    """Build tau matrix exactly as in S2753."""
    idx = list(range(-N, N + 1))
    D = len(idx)

    # Use production grid sizes from S2753
    R = TWO_PI * N / L + 80.0
    G = max(60000, 3000 * N)

    print(f"\nN={N}: dim={D}, R={R:.1f}, G={G}")

    t0 = time.time()
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h)
    w[0] = w[-1] = 0.5 * h

    Om = omega_vec(grid)
    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        E[k] = e_hat_vec(grid, n)

    A = E * (w * Om)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)
    dt = time.time() - t0
    print(f"  arch built: {dt:.1f}s")

    t0 = time.time()
    w02 = w02_mat(idx)
    print(f"  w02 built: {time.time()-t0:.2f}s")

    t0 = time.time()
    prime = prime_mat(idx)
    print(f"  prime built: {time.time()-t0:.2f}s")

    tau = arch + w02 - prime
    tau = (tau + tau.T) / 2

    return tau, idx

def analyze_gap(N):
    """Compute and analyze the eigenvalue gap."""
    tau, idx = build_tau(N)

    # Compute eigendecomposition
    print(f"\n  Computing eigendecomposition...")
    t0 = time.time()
    evals, evecs = np.linalg.eigh(tau)
    dt = time.time() - t0
    print(f"  eigh completed: {dt:.2f}s")

    # Analyze smallest eigenvalues
    lambda_1 = evals[0]
    lambda_2 = evals[1]
    gap = lambda_2 - lambda_1

    print(f"\n  Results:")
    print(f"    lambda_1: {lambda_1:.6e}")
    print(f"    lambda_2: {lambda_2:.6e}")
    print(f"    gap (λ₂ - λ₁): {gap:.6e}")

    # Compute residuals to assess reliability
    # res = ||τ v - λ v||
    v1 = evecs[:, 0]
    v2 = evecs[:, 1]

    res1 = np.linalg.norm(tau @ v1 - lambda_1 * v1)
    res2 = np.linalg.norm(tau @ v2 - lambda_2 * v2)

    # Condition number of tau
    cond = np.linalg.cond(tau)

    # Relative gaps
    rel_gap_1 = abs(gap) / abs(lambda_1) if lambda_1 != 0 else gap
    rel_gap_2 = abs(gap) / abs(lambda_2) if lambda_2 != 0 else gap

    print(f"\n  Reliability assessment:")
    print(f"    Residual (λ₁): {res1:.3e}")
    print(f"    Residual (λ₂): {res2:.3e}")
    print(f"    Matrix condition #: {cond:.3e}")
    print(f"    Relative gap / λ₁: {rel_gap_1:.3e}")
    print(f"    Relative gap / λ₂: {rel_gap_2:.3e}")
    print(f"    Float64 epsilon: {np.finfo(float).eps:.3e}")
    print(f"    Gap > epsilon? {gap > np.finfo(float).eps}")
    print(f"    Gap / epsilon: {gap / np.finfo(float).eps:.1f}x")

    # Compare to S2753
    if N == 20:
        s2753_gap = 5.94e-14
    elif N == 40:
        s2753_gap = 3.92e-15
    else:
        s2753_gap = None

    if s2753_gap:
        print(f"\n  Comparison to S2753:")
        print(f"    S2753 reported gap: {s2753_gap:.3e}")
        print(f"    This run gap: {gap:.3e}")
        print(f"    Match? {abs(gap - s2753_gap) / s2753_gap < 0.1}")
        print(f"    Ratio: {gap / s2753_gap:.3f}x")

    return {
        "N": N,
        "lambda_1": lambda_1,
        "lambda_2": lambda_2,
        "gap": gap,
        "residual_1": res1,
        "residual_2": res2,
        "condition": cond,
        "s2753_gap": s2753_gap if s2753_gap else None,
    }

if __name__ == "__main__":
    results = {}
    for N in [20, 40]:
        results[N] = analyze_gap(N)

    import json
    with open("/Users/scones/projects/riemann/results/ccm_gap_reproduce_s2753.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\n" + "="*70)
    print("Saved to /Users/scones/projects/riemann/results/ccm_gap_reproduce_s2753.json")
