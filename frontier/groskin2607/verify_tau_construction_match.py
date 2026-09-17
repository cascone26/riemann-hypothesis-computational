#!/usr/bin/env python3
"""
VERIFICATION: Do build_arb_tau (arb_ldlt_certify.py) and build_tau_and_terms
(cascade_generalizes_finiteN.py) construct the SAME τ matrix?

Test at small (c=13, N=4) using both methods.
"""
import sys
import numpy as np
from flint import arb, acb, arb_mat, ctx
import mpmath as mp

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from arb_ldlt_certify import build_arb_tau
from cascade_generalizes_finiteN import build_tau_and_terms

def test_match():
    c, N = 13, 4
    prec_arb = 300  # bits

    print("="*100)
    print(f"CONSTRUCTION MATCH TEST: c={c}, N={N}")
    print("="*100)

    # Method 1: Arb (arb_ldlt_certify.py)
    print("\n1. build_arb_tau (arb_ldlt_certify.py) at 300-bit Arb precision:")
    ctx.prec = prec_arb
    tau_arb, dim_arb = build_arb_tau(c, N, prec_arb)

    # Convert to mpmath for comparison
    mp.mp.dps = int(prec_arb * 0.29)  # Same as margin_law.py does
    tau_arb_mp = mp.matrix(dim_arb, dim_arb)
    for i in range(dim_arb):
        for j in range(dim_arb):
            tau_arb_mp[i,j] = mp.mpf(tau_arb[i,j].str(int(prec_arb*0.28), radius=False))

    print(f"   Dimension: {dim_arb}")
    print(f"   τ[0,0] = {tau_arb_mp[0,0]}")
    print(f"   τ[4,0] = {tau_arb_mp[4,0]}")
    print(f"   τ[4,4] = {tau_arb_mp[4,4]}")

    # Method 2: cascade_generalizes_finiteN.py (float64 + mpmath)
    print("\n2. build_tau_and_terms (cascade_generalizes_finiteN.py) at float64 + mpmath:")
    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau_cascade = W02_WR - Wp_total

    print(f"   Dimension: {2*N+1}")
    print(f"   τ[0,0] = {tau_cascade[0,0]}")
    print(f"   τ[4,0] = {tau_cascade[4,0]}")
    print(f"   τ[4,4] = {tau_cascade[4,4]}")

    # Compare entry-by-entry
    print("\n3. ENTRY-BY-ENTRY COMPARISON:")
    print(f"{'(i,j)':>6} {'arb_ldlt':>20} {'cascade':>20} {'rel_diff':>15} {'status':>8}")
    print("-"*72)

    max_rel_diff = 0
    test_indices = [(0,0), (0,4), (1,0), (1,4), (4,0), (4,4), (2,3), (3,2)]

    for (i,j) in test_indices:
        arb_val = float(tau_arb_mp[i,j])
        cascade_val = tau_cascade[i,j]

        # Relative difference
        scale = max(abs(arb_val), abs(cascade_val), 1e-100)
        rel_diff = abs(arb_val - cascade_val) / scale if scale > 0 else 0
        max_rel_diff = max(max_rel_diff, rel_diff)

        status = "✓" if rel_diff < 1e-8 else "✗" if rel_diff > 1e-4 else "~"
        print(f"({i},{j}) {arb_val:>20.10e} {cascade_val:>20.10e} {rel_diff:>15.2e} {status:>8}")

    print("\n4. SUMMARY:")
    print(f"   Max relative difference: {max_rel_diff:.2e}")
    if max_rel_diff < 1e-10:
        print("   ✓ CONSTRUCTIONS MATCH (numeric precision difference only)")
        return True
    elif max_rel_diff < 1e-4:
        print("   ~ CONSTRUCTIONS ROUGHLY MATCH (within float64 tolerance)")
        return True
    else:
        print("   ✗ CONSTRUCTIONS DIFFER SIGNIFICANTLY (possible bug or formula mismatch)")
        return False

if __name__ == '__main__':
    match = test_match()
    sys.exit(0 if match else 1)
