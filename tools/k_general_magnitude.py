"""
General-K magnitude polynomial: closed-form solution to the "quadruple/quintuple-product
expand" combinatorial wall that blocked K=6+ symbolic verification (S2583/S2585, this
session's S2681 first failed attempt, then S2682's correct derivation).

VERIFIED (S2682): each "other" block j contributes a biquadratic factor
    quartic_j(X) = X^4 - B_j*X^2 + A_j = (X^2 - z_j^2)(X^2 - zc_j^2)
where A_j = z_j^2*zc_j^2, B_j = z_j^2 + zc_j^2 (matching the real Lean source's Q-polynomial
factor structure, confirmed directly in CrossK5SixthCoeff.lean).

The (K-1)-fold product of these quartics, substituting Y=X^2, is
    M_K(Y) = prod_{j in others} (Y^2 - B_j*Y + A_j)
a polynomial of degree 2*(K-1) in Y (degree 4*(K-1) in X). The degree-b (X-degree, b even)
magnitude term used in Cross{K}Coeff's coefficient bricks is exactly the coefficient of
Y^(b/2) in M_K(Y).

VERIFIED EXACTLY against K=4's real, independently machine-verified magnitude_sympy formulas
(tools/build_k4reduced_generic.py) for ALL SIX of its b-values (0,2,4,6,8,10) -- exact term-
for-term match including sign. Also verified against K=5's real compiled Lean source (Q1's
constant term = product of all 4 A's, matching this formula's Y^0 coefficient).

Fast: sympy Poly arithmetic in a SINGLE variable Y, not a raw multi-variable expansion --
confirmed sub-10ms even at K=8 (7 other factors), where the old expansion-based approach
(S2583/S2585) hit a hard performance wall at K=6.
"""
import sympy as sp

def magnitude_polynomial(K, other_labels=None):
    """
    Returns the full (K-1)-fold magnitude polynomial M_K(Y) = prod_j (Y^2 - B_j*Y + A_j),
    as a sympy Poly in Y. other_labels: optional list of K-1 label strings for the symbols
    (default: numeric 0..K-2).
    """
    n_others = K - 1
    if other_labels is None:
        other_labels = [str(i) for i in range(n_others)]
    assert len(other_labels) == n_others
    syms = [sp.symbols(f'A{lab} B{lab}') for lab in other_labels]
    Y = sp.Symbol('Y')
    poly = sp.Poly(1, Y)
    for A, B in syms:
        poly = poly * sp.Poly(Y**2 - B*Y + A, Y)
    return poly, syms

def magnitude_coeff(K, b, other_labels=None, factor=False):
    """
    Returns the degree-b (X-degree, b even, 0<=b<=4*(K-1)) magnitude coefficient for a
    K-block construction, i.e. the coefficient of Y^(b/2) in the (K-1)-fold product.
    """
    assert b % 2 == 0 and 0 <= b <= 4 * (K - 1)
    poly, _ = magnitude_polynomial(K, other_labels)
    c = poly.coeff_monomial(sp.Symbol('Y') ** (b // 2))
    return sp.factor_terms(c) if factor else c

if __name__ == "__main__":
    # Self-verification against K=4's real, independently-proven formulas.
    Aj,Ak,Al,Bj,Bk,Bl = sp.symbols('Aj Ak Al Bj Bk Bl')
    known_k4 = {
        10: (-1, Bj+Bk+Bl),
        8:  (1, Bj*Bk+Bj*Bl+Bk*Bl+Aj+Ak+Al),
        6:  (-1, Aj*Bk+Aj*Bl+Ak*Bj+Ak*Bl+Al*Bj+Al*Bk+Bj*Bk*Bl),
        4:  (1, Aj*Ak+Aj*Al+Ak*Al+Aj*Bk*Bl+Ak*Bj*Bl+Al*Bj*Bk),
        2:  (-1, Aj*Ak*Bl+Aj*Al*Bk+Ak*Al*Bj),
        0:  (1, Aj*Ak*Al),
    }
    all_match = True
    for b, (osign, mag) in known_k4.items():
        computed = magnitude_coeff(4, b, other_labels=['j','k','l'])
        expected = osign * mag
        match = sp.simplify(computed - expected) == 0
        all_match = all_match and match
        print(f"K=4 b={b}: match={match}")
    print(f"\nALL MATCH: {all_match}")

    # Generate K=6 formulas (previously blocked)
    print("\n--- K=6 formulas (b=0..20), newly derivable ---")
    for b in range(0, 4*5+1, 2):
        c = magnitude_coeff(6, b, other_labels=['j','k','l','m','n'], factor=True)
        print(f"b={b}: {c}")
