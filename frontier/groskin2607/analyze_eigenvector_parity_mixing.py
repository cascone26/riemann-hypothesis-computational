#!/usr/bin/env python3
"""
Confirm the coupling mechanism: the global minimum eigenvector
mixes BOTH even and odd parity components significantly.
"""
import sys
import numpy as np
from flint import arb, ctx
import mpmath as mp

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from arb_ldlt_certify import build_arb_tau

def analyze_min_eigenvector_parity(c, N, prec):
    """Compute global minimum eigenvector and analyze its parity composition."""
    ctx.prec = prec
    tau_arb, dim = build_arb_tau(c, N, prec)

    # Convert to mpmath
    mp.mp.dps = int(prec * 0.29)
    tau_mp = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            tau_mp[i, j] = mp.mpf(tau_arb[i, j].str(int(prec*0.28), radius=False))

    # Convert to numpy for eigsh
    tau_f64 = np.array([[float(tau_mp[i, j]) for j in range(dim)] for i in range(dim)])

    # Compute all eigenvalues/eigenvectors
    evals, evecs = np.linalg.eigh(tau_f64)

    # Minimum eigenvalue and corresponding eigenvector
    lambda_min = evals[0]
    psi_min = evecs[:, 0]

    print("="*100)
    print(f"EIGENVECTOR PARITY ANALYSIS: c={c}, N={N}")
    print("="*100)
    print(f"\nGlobal minimum eigenvalue: λ_min = {lambda_min:.6e}")
    print(f"Eigenvector norm: ||ψ|| = {np.linalg.norm(psi_min):.6f}\n")

    # Decompose into even/odd components
    even_idx = np.array([i for i in range(dim) if (i - N) % 2 == 0])
    odd_idx = np.array([i for i in range(dim) if (i - N) % 2 != 0])

    psi_even = psi_min[even_idx]
    psi_odd = psi_min[odd_idx]

    norm_even = np.linalg.norm(psi_even)
    norm_odd = np.linalg.norm(psi_odd)
    norm_total = np.linalg.norm(psi_min)

    frac_even = (norm_even / norm_total) ** 2  # Energy fraction
    frac_odd = (norm_odd / norm_total) ** 2

    print(f"Parity decomposition (energy fractions):")
    print(f"  Even modes: {frac_even:.4f} ({norm_even:.6f})")
    print(f"  Odd modes:  {frac_odd:.4f} ({norm_odd:.6f})")
    print(f"  Total:      {frac_even + frac_odd:.4f}")

    if frac_even > 0.3 and frac_odd > 0.3:
        print("\n✗ STRONG MIXING: Both parities are significantly present in ψ_min")
        print("  → The minimum eigenstate is a superposition of even and odd modes")
        print("  → This proves the sectors are COUPLED and cannot be treated independently")
    elif frac_even > 0.7 or frac_odd > 0.7:
        print("\n~ PARTIAL COUPLING: Dominated by one parity but other is not negligible")
    else:
        print("\n✓ NEARLY DECOUPLED: Eigenvector is mostly in one sector")

    print("\nEigenvector components (first 10 modes):")
    print(f"{'i':>3} {'n':>3} {'ψ[i]':>15} {'parity':>8} {'|ψ[i]|²':>15}")
    print("-"*50)
    for i in range(min(10, dim)):
        n = i - N
        parity = "even" if n % 2 == 0 else "odd"
        print(f"{i:>3} {n:>3} {psi_min[i]:>15.6e} {parity:>8} {psi_min[i]**2:>15.6e}")

if __name__ == '__main__':
    analyze_min_eigenvector_parity(30, 4, 300)
    print("\n" + "="*100 + "\n")
    analyze_min_eigenvector_parity(30, 8, 300)
