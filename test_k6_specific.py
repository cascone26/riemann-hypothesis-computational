#!/usr/bin/env python3
"""
Test K=6 cubic variants with SPECIFIC parameter structures (not abstract free variables).

The key insight from SESSION_K6_CUBIC_FIX files: the hsum constraint structure is what
determines whether the polynomial identity can hold. Test:

1. Current cubic (degree 3) with K=6 parameters → should NOT vanish (verified)
2. Degree-4 cubic for K=6 block (while keeping K-1 blocks as quartics)
3. K-dependent degree (e.g., degree = K-2 = 4 for K=6)
"""

import sympy as sp
from sympy import symbols, expand, collect, simplify, Rational

def test_k6_cubic_variant(variant_name, cubic_deg, other_degs):
    """
    Test a K=6 variant where the designated block has degree cubic_deg,
    and other blocks have degrees other_degs (list of 5 degrees, or a constant).
    
    Return whether the sum-of-products can plausibly vanish (rough heuristic test).
    """
    K = 6
    delta = symbols('delta', real=True)
    
    # For simplicity, test with SPECIFIC numerical coefficient patterns
    # (matching the actual CrossK5Key structure in spirit)
    
    # Quartic form: delta^4 + a*delta^2 + b
    def make_poly(deg, block_idx):
        # Use specific coefficients that would appear in the actual problem
        # (These are illustrative; the real values come from the hsum constraint)
        coeffs = [Rational(1, 1) if k == deg else symbols(f'c{block_idx}_{k}') for k in range(deg+1)]
        return sum(coeffs[k] * delta**k for k in range(deg+1))
    
    # Build sum-of-products
    total = 0
    designated_deg = cubic_deg if isinstance(cubic_deg, int) else cubic_deg
    
    for i in range(1, K+1):
        # Block i uses the test degree
        if isinstance(cubic_deg, int):
            poly_i = make_poly(cubic_deg, i)
        else:
            poly_i = cubic_deg  # already a polynomial
        
        # Other blocks use their specified degrees
        product_others = 1
        for j in range(1, K+1):
            if j != i:
                if isinstance(other_degs, list):
                    other_deg = other_degs[j-1] if j-1 < len(other_degs) else 4
                else:
                    other_deg = other_degs
                product_others *= make_poly(other_deg, j)
        
        total += poly_i * product_others
    
    # Try to simplify
    total_expanded = expand(total)
    
    # Check highest degree term
    poly_total = sp.Poly(total_expanded, delta)
    highest_deg = poly_total.degree()
    
    print(f"\n{variant_name}:")
    print(f"  Highest degree in sum-of-products: {highest_deg}")
    print(f"  (Expected max degree: cubic_deg={cubic_deg} + 4*(K-1)={4*(K-1)} = {cubic_deg + 4*(K-1)})")
    
    # Quick heuristic: can it vanish?
    # In reality, we'd need to solve the full hsum constraint system
    # For now, just report the structure
    return highest_deg


print("=" * 80)
print("K=6 CUBIC VARIANT TESTS")
print("=" * 80)

# Variant 1: Current (degree 3 cubic)
test_k6_cubic_variant("Variant 1: Degree-3 cubic (current)", cubic_deg=3, other_degs=4)

# Variant 2: Degree-4 cubic (all-quartic style)
test_k6_cubic_variant("Variant 2: Degree-4 cubic (symmetric)", cubic_deg=4, other_degs=4)

# Variant 3: Degree-5 cubic
test_k6_cubic_variant("Variant 3: Degree-5 cubic", cubic_deg=5, other_degs=4)

# Variant 4: K-scaled (degree = K-2 = 4 for K=6)
test_k6_cubic_variant("Variant 4: Degree K-2 = 4 (K-scaled)", cubic_deg=4, other_degs=4)

print("\n" + "=" * 80)
print("INTERPRETATION")
print("=" * 80)
print("""
The symbolic structure alone doesn't tell us if a variant works.
The key is: what does the K-dependent hsum constraint allow?

For K=3,4,5: The cubic-based structure works with the specific hsum constraint.
For K=6: The hsum constraint structure changes, breaking the cubic identity.

To actually test if a variant works, we would need to:
1. Build the hsum constraint for K=6 with the variant polynomial
2. Check if the constraint can be satisfied (null space exists)
3. Verify the sum-of-products vanishes when the constraint is applied

This requires either:
- Lean compilation (definitive, Handle-tier)
- Symbolic constraint solving in sympy (research-scale, may be ill-conditioned)
""")

