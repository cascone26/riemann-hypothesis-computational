#!/usr/bin/env python3
"""
Check: Is τ actually block-diagonal, or are there off-diagonal blocks
between even and odd frequency components?

This is the critical check to understand why min(λ_even, λ_odd) < min(λ_full).
"""
import sys
import numpy as np
from flint import arb, ctx
import mpmath as mp

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from arb_ldlt_certify import build_arb_tau

def check_block_structure(c, N, prec):
    """Build τ and check for coupling between even/odd frequency modes."""
    ctx.prec = prec
    tau_arb, dim = build_arb_tau(c, N, prec)

    # Convert to float for inspection
    mp.mp.dps = int(prec * 0.29)
    tau_mp = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            tau_mp[i, j] = mp.mpf(tau_arb[i, j].str(int(prec*0.28), radius=False))

    # Convert to numpy float64 for display
    tau_f64 = np.array([[float(tau_mp[i, j]) for j in range(dim)] for i in range(dim)])

    # Check: are there significant entries between even and odd rows/cols?
    even_idx = [i for i in range(dim) if (i - N) % 2 == 0]
    odd_idx = [i for i in range(dim) if (i - N) % 2 != 0]

    print("="*100)
    print(f"BLOCK DIAGONALITY CHECK: c={c}, N={N}, dim={dim}")
    print("="*100)

    # Extract off-diagonal blocks
    tau_even_odd = tau_f64[np.ix_(even_idx, odd_idx)]
    tau_odd_even = tau_f64[np.ix_(odd_idx, even_idx)]

    max_even_odd = np.max(np.abs(tau_even_odd))
    max_odd_even = np.max(np.abs(tau_odd_even))
    max_full = np.max(np.abs(tau_f64))

    print(f"\nEven-Odd coupling (even rows × odd cols):")
    print(f"  Max absolute value: {max_even_odd:.2e}")
    print(f"  As fraction of max(τ): {max_even_odd / max_full if max_full > 0 else 0:.2e}")
    print(f"  Sample entries:")
    for i_idx in range(min(3, len(even_idx))):
        for j_idx in range(min(3, len(odd_idx))):
            i, j = even_idx[i_idx], odd_idx[j_idx]
            print(f"    τ[{i},{j}] = {tau_f64[i, j]:+.6e}")

    print(f"\nOdd-Even coupling (odd rows × even cols):")
    print(f"  Max absolute value: {max_odd_even:.2e}")
    print(f"  As fraction of max(τ): {max_odd_even / max_full if max_full > 0 else 0:.2e}")
    print(f"  Sample entries:")
    for i_idx in range(min(3, len(odd_idx))):
        for j_idx in range(min(3, len(even_idx))):
            i, j = odd_idx[i_idx], even_idx[j_idx]
            print(f"    τ[{i},{j}] = {tau_f64[i, j]:+.6e}")

    print("\n" + "="*100)
    if max_even_odd < 1e-14 * max_full and max_odd_even < 1e-14 * max_full:
        print("✓ BLOCK-DIAGONAL: Off-diagonal blocks are machine-noise only")
        return True
    elif max_even_odd < 1e-6 * max_full and max_odd_even < 1e-6 * max_full:
        print("~ NEARLY BLOCK-DIAGONAL: Small coupling exists but might be numerical artifact")
        return True
    else:
        print("✗ NOT BLOCK-DIAGONAL: Significant coupling between even and odd sectors!")
        return False

if __name__ == '__main__':
    check_block_structure(30, 4, 300)
    check_block_structure(30, 8, 300)
