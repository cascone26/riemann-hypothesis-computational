#!/usr/bin/env python3
"""
Test degree-4+ polynomial generalization for K=6+ induction.

Two main hypotheses:
A. All-quartic: Replace cubic_i with a degree-4 "quartic_cubic_i" for the i-th block,
   making ALL K blocks symmetric (K-1 of them already quartic, now K-th one is also degree-4).
B. K-scaled degree: Find a K-dependent degree formula that works for K=3,4,5 and fix K=6.
   E.g.: degree = K-2 (gives 3 for K=5, 4 for K=6) or degree = f(K) for some pattern.

Test structure:
1. For each K in [3,4,5,6], for each hypothesis:
   - Generate the symbolic polynomial identity
   - Check if it vanishes (all coefficients in delta = 0)
   - Report findings

2. If a candidate works for K=3,4,5, verify it ALSO works for K=6 (fix condition)
"""

import sympy as sp
from sympy import symbols, expand, collect, Poly, simplify
import sys

def build_sum_of_products_cubic_based(K):
    """
    Current K-induction setup: K-1 quartics + 1 cubic.

    Symbolic structure (ignoring A_i, z_i details, focusing on delta = 2π(n-N)/L):
    - Each block j (j ≠ i) contributes a quartic: q_j(delta) = delta^4 + ...
    - Block i contributes a cubic: c_i(delta) = delta^3 + ...
    - Sum-of-products: ∑_i c_i(delta) * ∏_{j≠i} q_j(delta)

    For symbolic testing, use simplified forms:
    - quartic_j = delta^4 + a_j*delta^2 + b_j  (generic biquadratic form)
    - cubic_i = delta^3 + c_i*delta  (generic cubic form, symmetric in structure)
    """
    delta = symbols('delta', real=True)

    # Generic quartic (degree 4): represents the K-1 "other" blocks
    def quartic_generic(j):
        a, b = symbols(f'a{j} b{j}', real=True)
        return delta**4 + a*delta**2 + b

    # Current cubic structure (degree 3)
    def cubic_generic(i):
        c = symbols(f'c{i}', real=True)
        return delta**3 + c*delta

    # Sum-of-products: each term i = cubic_i * ∏_{j≠i} quartic_j
    total = 0
    for i in range(1, K+1):
        cubic_i = cubic_generic(i)
        product_others = 1
        for j in range(1, K+1):
            if j != i:
                product_others *= quartic_generic(j)
        total += cubic_i * product_others

    return expand(total), delta


def build_sum_of_products_all_quartic(K):
    """
    Option A: All-quartic (symmetric degree-4 for all K blocks).

    Replace cubic_i with a degree-4 polynomial q_i, making all K blocks symmetric.
    Each block gets the SAME structure, a quartic.

    q_i(delta) = delta^4 + alpha_i*delta^2 + beta_i  (biquadratic, degree 4)
    Sum-of-products: ∑_i q_i(delta) * ∏_{j≠i} q_j(delta)

    This is now FULLY SYMMETRIC in all K blocks.
    """
    delta = symbols('delta', real=True)

    def quartic_generic(i):
        alpha, beta = symbols(f'alpha{i} beta{i}', real=True)
        return delta**4 + alpha*delta**2 + beta

    total = 0
    for i in range(1, K+1):
        q_i = quartic_generic(i)
        product_others = 1
        for j in range(1, K+1):
            if j != i:
                product_others *= quartic_generic(j)
        total += q_i * product_others

    return expand(total), delta


def build_sum_of_products_k_scaled_degree(K, degree):
    """
    Option B: K-scaled degree.

    Parameterize the designated block i's polynomial by a K-dependent degree.
    For K=3,4,5, degree should equal 3 (current form).
    For K=6, try degree=4 (K-2 formula) or other patterns.

    poly_i(delta) = sum_{k=0}^{degree} p_{i,k} * delta^k
    For simplicity, test with "generic degree-N" where we set all coefficients as free symbols.
    """
    delta = symbols('delta', real=True)

    # For K-1 blocks, use standard quartic
    def quartic_generic(j):
        a, b = symbols(f'a{j} b{j}', real=True)
        return delta**4 + a*delta**2 + b

    # For block i, use degree-N polynomial
    # To keep it manageable, use a specific "canonical" form with specific coefficients
    poly_i_coeffs = [symbols(f'p{k}', real=True) for k in range(degree+1)]
    poly_i = sum(poly_i_coeffs[k] * delta**k for k in range(degree+1))

    total = 0
    for i in range(1, K+1):
        # Use same generic form for all i (they differ only in the free parameters)
        product_others = 1
        for j in range(1, K+1):
            if j != i:
                product_others *= quartic_generic(j)
        total += poly_i * product_others

    # Note: this is a REPRESENTATIVE test (all designated blocks have the same degree,
    # parameterized by p_k). In a full test, we'd need to set p_k values, but for
    # vanishing-coefficient check, the structure is what matters.
    return expand(total), delta


def test_vanishing_coefficients(expr, delta, K, name):
    """
    Check if the expression vanishes by collecting coefficients in delta.
    Return a dict of {power: coefficient} for all nonzero terms.
    """
    poly = Poly(expr, delta)
    coeffs = poly.all_coeffs()  # ordered from highest to lowest degree

    nonzero = {}
    for deg, coeff in enumerate(reversed(coeffs)):
        coeff_simplified = simplify(coeff)
        if coeff_simplified != 0:
            nonzero[deg] = coeff_simplified

    return nonzero


def main():
    print("=" * 80)
    print("DEGREE-4+ GENERALIZATION TEST FOR K-INDUCTION")
    print("=" * 80)

    results = {}

    # Test Option A: All-quartic (symmetric)
    print("\n### OPTION A: ALL-QUARTIC (Symmetric Degree-4 for All Blocks) ###\n")
    for K in [3, 4, 5, 6]:
        expr, delta = build_sum_of_products_all_quartic(K)
        nonzero = test_vanishing_coefficients(expr, delta, K, f"Option A (K={K})")
        results[f"OptA_K{K}"] = nonzero
        if nonzero:
            print(f"K={K}: NONZERO coefficients in delta (not vanishing)")
            for deg in sorted(nonzero.keys(), reverse=True):
                print(f"  delta^{deg}: {nonzero[deg]}")
        else:
            print(f"K={K}: ✓ VANISHES (all coefficients zero)")

    # Test Option B: K-scaled degree
    print("\n### OPTION B: K-SCALED DEGREE ###\n")

    # Hypothesis B1: degree = K-2 (gives 3 for K=5, 4 for K=6)
    print("Hypothesis B1: degree = K-2")
    for K in [3, 4, 5, 6]:
        deg = K - 2
        print(f"K={K}: degree = {deg}")
        # This is complex to test symbolically, so we'll note the structure
        # and skip full computation (too many free parameters)
        print(f"  [Complex to test fully; structure would have degree-{deg} polynomial in block i]")

    # Hypothesis B2: degree = 3 (constant, current form)
    print("\nHypothesis B2: degree = 3 (constant, current form)")
    for K in [3, 4, 5, 6]:
        expr, delta = build_sum_of_products_cubic_based(K)
        nonzero = test_vanishing_coefficients(expr, delta, K, f"Current cubic (K={K})")
        results[f"OptB_current_K{K}"] = nonzero
        status = "✓ VANISHES" if not nonzero else f"NONZERO at {len(nonzero)} powers"
        print(f"K={K}: {status}")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    all_a_vanish = all(not results[f"OptA_K{K}"] for K in [3,4,5,6])
    all_b_vanish = all(not results[f"OptB_current_K{K}"] for K in [3,4,5])
    b6_vanish = not results[f"OptB_current_K6"]

    print(f"Option A (all-quartic): {'✓ All K vanish' if all_a_vanish else '✗ Some K nonzero'}")
    if all_a_vanish:
        print("  → ALL-QUARTIC WORKS for K=3,4,5,6!")
        print("  → Recommendation: Test all-quartic in Lean for K=5,6 to verify fix")

    print(f"Option B (current cubic): K=3,4,5 {'vanish ✓' if all_b_vanish else 'nonzero ✗'}, K=6 {'vanishes ✓' if b6_vanish else 'nonzero ✗'}")
    if not b6_vanish:
        print("  → Current cubic fails at K=6 (confirmed)")

    if not all_a_vanish:
        print("\n⚠ Option A (all-quartic) does not universally vanish. See details above.")

    return results


if __name__ == "__main__":
    try:
        results = main()
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
