#!/usr/bin/env python3
"""
Hybrid approach: numerical nullspace extraction + high-precision numerical verification
+ symbolic cubic evaluation on exact coefficients extracted from numerical data.
"""

import numpy as np
from sympy import (
    Matrix, Rational, I, pi, sqrt, sin, conjugate, simplify, expand,
    N as sym_N, nsimplify
)

N_freq = 10
L_val = 2 * pi

configs = {
    7: {'name': 'Symmetric (-2..2) + 4i',
        'z_vals': [-2+4*I, -1+4*I, 0+4*I, 1+4*I, 2+4*I]},
    9: {'name': 'SmallReal (±0.01) + 7i',
        'z_vals': [-Rational(1,100)+7*I, -Rational(5,100)+7*I, 0+7*I,
                   Rational(5,100)+7*I, Rational(1,100)+7*I]},
    14: {'name': 'NegReal (-5..3) + 6i',
         'z_vals': [-5+6*I, -3+6*I, -1+6*I, 1+6*I, 3+6*I]}
}


def v_plus_numeric(n, z):
    L = 2 * np.pi
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z + freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom


def v_minus_numeric(n, z):
    L = 2 * np.pi
    freq_shift = 2 * np.pi * (n - N_freq) / L
    denom = z - freq_shift
    num = np.sqrt(2) * np.sin(L * z / 2)
    return num / denom


def build_matrix_numeric(z_sym_vals):
    """Build constraint matrix numerically."""
    z_vals = [complex(float(sym_N(z.as_real_imag()[0])),
                      float(sym_N(z.as_real_imag()[1]))) for z in z_sym_vals]
    L = 2 * np.pi
    M = []
    for n in range(2 * N_freq + 1):
        row = []
        for z in z_vals:
            z_conj = np.conj(z)
            row.extend([v_plus_numeric(n, z), v_minus_numeric(n, z),
                       v_plus_numeric(n, z_conj), v_minus_numeric(n, z_conj)])
        M.append(row)
    return np.array(M, dtype=complex)


def extract_nullspace_numeric(M):
    """Extract nullspace using SVD."""
    from scipy.linalg import svd
    _, s, Vh = svd(M)
    rank = np.sum(s > 1e-10)
    dim = M.shape[1] - rank
    null_vecs = Vh[rank:, :].T
    return null_vecs, dim


def v_plus_sym(n, L, z):
    freq_shift = 2 * pi * (n - N_freq) / L
    return sqrt(2) * sin(L * z / 2) / (z + freq_shift)


def v_minus_sym(n, L, z):
    freq_shift = 2 * pi * (n - N_freq) / L
    return sqrt(2) * sin(L * z / 2) / (z - freq_shift)


def verify_hsum_numeric(coeffs_exact, z_vals):
    """Verify hsum constraint numerically with high precision."""
    z_numeric = [complex(float(sym_N(z.as_real_imag()[0])),
                         float(sym_N(z.as_real_imag()[1]))) for z in z_vals]
    M = build_matrix_numeric(z_vals)

    coeffs_num = [complex(float(sym_N(c.as_real_imag()[0])) if hasattr(c, 'as_real_imag') else float(sym_N(c)),
                          float(sym_N(c.as_real_imag()[1])) if hasattr(c, 'as_real_imag') else 0)
                  for c in coeffs_exact]

    residual = M @ np.array(coeffs_num)
    residual_norm = np.linalg.norm(residual)
    return residual_norm < 1e-9, residual_norm


def cubic_sym(z_vals, coeffs_exact, n_val):
    """Evaluate cubic symbolically using exact coefficients."""
    t = 2 * pi * (n_val - N_freq) / L_val

    result = 0
    for i, z_i in enumerate(z_vals):
        z_conj = conjugate(z_i)
        A_i = sqrt(2) * sin(L_val * z_i / 2)
        A_conj = conjugate(A_i)

        p_i, q_i, r_i, s_i = coeffs_exact[4*i:4*i+4]

        cubic_i = (
            (p_i + q_i) * A_i * (z_i * z_conj**2)
            + (r_i + s_i) * A_conj * (z_i**2 * z_conj)
            + (q_i - p_i) * A_i * (z_conj**2) * t
            + (s_i - r_i) * A_conj * (z_i**2) * t
            - (p_i + q_i) * A_i * z_i * (t**2)
            - (r_i + s_i) * A_conj * z_conj * (t**2)
            + (p_i - q_i) * A_i * (t**3)
            + (r_i - s_i) * A_conj * (t**3)
        )

        quartic_prod = 1
        for j, z_j in enumerate(z_vals):
            if i != j:
                z_conj_j = conjugate(z_j)
                quartic_prod *= (z_j + t) * (z_j - t) * (z_conj_j + t) * (z_conj_j - t)

        result += cubic_i * quartic_prod

    return result


def test_config_hybrid(idx):
    config = configs[idx]
    z_vals = config['z_vals']

    print(f"\nCONFIG {idx}: {config['name']}")
    print("-" * 60)

    # Numeric extraction
    print("  Extracting nullspace numerically...")
    M = build_matrix_numeric(z_vals)
    null_vecs, dim = extract_nullspace_numeric(M)
    print(f"    Nullspace dim: {dim}")

    # Use first vector
    coeffs_num = null_vecs[:, 0]

    # Convert to exact rationals carefully
    print("  Converting to exact rationals...")
    coeffs_exact = []
    for c in coeffs_num:
        if abs(c.imag) < 1e-14:
            coeffs_exact.append(nsimplify(c.real, rational=True, tolerance=1e-12))
        else:
            re_rational = nsimplify(c.real, rational=True, tolerance=1e-12)
            im_rational = nsimplify(c.imag, rational=True, tolerance=1e-12)
            coeffs_exact.append(re_rational + im_rational*I)

    # Verify constraint
    print("  Verifying hsum constraint...")
    is_satisfied, residual_norm = verify_hsum_numeric(coeffs_exact, z_vals)
    print(f"    Residual norm: {residual_norm:.2e}")

    if residual_norm > 1e-6:
        print("    WARNING: Large residual, solution may not be accurate")
        return {'status': 'UNCERTAIN', 'residual': residual_norm}

    # Test cubic at n=10
    print("  Testing cubic at n=10...")
    try:
        poly_val = cubic_sym(z_vals, coeffs_exact, N_freq)
        poly_simp = simplify(expand(poly_val))

        if poly_simp == 0:
            print("    ✓ CUBIC = 0")
            return {'status': 'PASS'}
        else:
            print(f"    ✗ CUBIC ≠ 0: {str(poly_simp)[:80]}...")
            return {'status': 'FAIL', 'value': poly_simp}
    except Exception as e:
        print(f"    ERROR: {str(e)[:60]}")
        return {'status': 'ERROR', 'error': str(e)}


if __name__ == '__main__':
    print("=" * 60)
    print("K=5 CUBIC TEST (HYBRID NUMERIC-SYMBOLIC)")
    print("=" * 60)

    results = {}
    for idx in [7, 9, 14]:
        results[idx] = test_config_hybrid(idx)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for idx in [7, 9, 14]:
        print(f"Config {idx}: {results[idx]['status']}")
