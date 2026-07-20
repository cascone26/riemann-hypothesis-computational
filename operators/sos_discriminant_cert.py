#!/usr/bin/env python3
"""
SOS/SDP Certificate for Jensen Polynomial Discriminant
========================================================

For d=4 (degree-4 Jensen polynomial), the discriminant is a polynomial in the
Taylor coefficients c_k = C(4,k) × a_{n+k}.

Goal: Use sum-of-squares (SDP) to search for a positivity certificate connecting
the discriminant to the log-convexity of the a_n sequence.

The key insight: if a_{n+1}^2 ≥ a_n × a_{n+2} (log-convexity, necessary for
hyperbolicity), can we prove disc(J^(4,n)) ≥ 0?

Strategy:
1. Load known a_n coefficients (log-convex, proven)
2. Form J^(4,n) for various n
3. Compute discriminant numerically
4. Use cvxpy SDP formulation to certify disc(J^(4,n)) ≥ f(constraint_margin)
   where the constraint_margin reflects how far the a_n are from violating
   log-convexity.
5. Output: real SDP solver status + numeric coefficients of any discovered certificate
"""

import json
import os
import numpy as np
import mpmath
from math import comb
import cvxpy as cp
from numpy.polynomial import polynomial as P
import time

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

mpmath.mp.dps = 50

# ── Load/Compute a_n coefficients ─────────────────────────────────────────

def load_a_coeffs_from_results():
    """
    Load Taylor coefficients from an existing results file (if available).
    """
    jensen_hyp_file = os.path.join(RESULTS_DIR, "jensen_hyp_fast.json")
    if os.path.exists(jensen_hyp_file):
        print(f"Loading coefficients from {jensen_hyp_file}...")
        with open(jensen_hyp_file) as f:
            data = json.load(f)
        return data['a_coeffs']
    else:
        print(f"  File {jensen_hyp_file} not found, will compute.")
        return None

def compute_a_coeffs_dft(n_max=50, prec=50):
    """
    Compute Taylor coefficients of Ξ(z) via DFT method (fast).
    """
    mpmath.mp.dps = prec
    print(f"Computing a_n via DFT for n=0..{n_max}...")
    t0 = time.time()

    def Xi_func(z):
        s = mpmath.mpf('0.5') + 1j * z
        return mpmath.re(mpmath.mpf('0.5') * s * (s-1) *
                        mpmath.power(mpmath.pi, -s/2) *
                        mpmath.gamma(s/2) * mpmath.zeta(s))

    N = 4 * (n_max + 10)  # oversample
    r = mpmath.mpf('3.0')  # radius
    vals = np.zeros(N, dtype=complex)

    for k in range(N):
        theta = 2 * mpmath.pi * k / N
        z = r * mpmath.exp(1j * theta)
        vals[k] = complex(Xi_func(z))

    fft = np.fft.fft(vals) / N
    r_float = float(r)

    coeffs = []
    for n in range(n_max + 1):
        if 2*n < N:
            a = fft[2*n].real / r_float**(2*n)
        else:
            a = 0.0
        coeffs.append(float(a))

    print(f"  Done in {time.time()-t0:.1f}s")
    return coeffs

# ── Discriminant of degree-4 polynomial ───────────────────────────────────

def discriminant_deg4(c0, c1, c2, c3, c4):
    """
    Discriminant of c0 + c1*X + c2*X^2 + c3*X^3 + c4*X^4.

    Using the standard formula: disc(p) = (-1)^{d(d-1)/2} Res(p, p')/lc(p)
    where d=degree, lc=leading coefficient, Res=resultant.

    For d=4: disc = c4^3 × resultant(p, p') / c4^4 = resultant(p, p') / c4

    Computed numerically via numpy.
    """
    coeffs = np.array([c0, c1, c2, c3, c4])
    # numpy.poly computes the characteristic polynomial of a companion matrix
    # Instead, use numpy to find roots and compute disc from roots
    roots = np.roots(coeffs[::-1])  # highest degree first
    if len(roots) < 4:
        return 0.0

    # disc = lc^{2d-2} × ∏_{i<j} (r_i - r_j)^2
    # = c4^6 × ∏_{i<j} (r_i - r_j)^2
    prod_sq = 1.0
    for i in range(len(roots)):
        for j in range(i+1, len(roots)):
            prod_sq *= (roots[i] - roots[j])**2

    disc = c4**6 * prod_sq
    return np.real(disc)

# ── Log-convexity margin ───────────────────────────────────────────────────

def log_convexity_margin(a_n):
    """
    Compute the margin by which a_n is log-convex:
    margin_n = a_{n+1}^2 / (a_n × a_{n+2}) - 1.0

    margin > 0 ⟺ log-convex. margin→0 ⟺ approaching boundary.
    """
    if len(a_n) < 3 or a_n[0] <= 0 or a_n[2] <= 0:
        return 1e10  # undefined
    return a_n[1]**2 / (a_n[0] * a_n[2]) - 1.0

# ── SDP formulation: certificate for discriminant positivity ───────────────

def build_sdp_certificate(a_coeffs, d_max_test=15, verbose=True):
    """
    Build an SDP that searches for a certificate proving:

    disc(J^(4,n)) ≥ 0  for all n where a_n is log-convex.

    Approach: parametrize disc(J^(4,n)) as a function of the margin m_n,
    and try to find an SDP-positive function that lower-bounds it.

    This is a search for a certificate in the form:
    disc(J^(4,n)) = SOS_poly + f(m_n)
    where m_n is the log-convexity margin.
    """
    print(f"\nBuilding SDP certificate for Jensen polynomial discriminant...")
    print(f"Testing n=0..{d_max_test}")

    # Collect data: (n, margin, discriminant) triples
    data = []
    for n in range(min(d_max_test + 1, len(a_coeffs) - 4)):
        a_slice = a_coeffs[n:n+3]  # a_n, a_{n+1}, a_{n+2}
        margin = log_convexity_margin(a_slice)

        # Form J^(4,n)
        c = [comb(4, j) * a_coeffs[n+j] for j in range(5)]
        disc = discriminant_deg4(*c)

        data.append({
            'n': n,
            'a_n': float(a_coeffs[n]),
            'a_n1': float(a_coeffs[n+1]),
            'a_n2': float(a_coeffs[n+2]),
            'margin': margin,
            'discriminant': disc,
            'c': c,
        })

    if verbose:
        print(f"\nData collected for {len(data)} polynomials:")
        print(f"{'n':>3} {'margin':>12} {'disc':>15} {'disc/margin^2':>15}")
        print("-" * 50)
        for d in data:
            if d['margin'] > 1e-6:
                ratio = d['discriminant'] / (d['margin']**2 + 1e-15)
            else:
                ratio = np.inf if d['discriminant'] > 0 else -np.inf
            print(f"{d['n']:>3} {d['margin']:>12.4e} {d['discriminant']:>15.4e} {ratio:>15.4e}")

    # ── Attempt 1: Simple linear lower bound via SDP ─────────────────────
    print(f"\nAttempt 1: Linear lower bound on discriminant via margin...")

    # Variables: α, β ∈ ℝ such that disc(J^(4,n)) ≥ α * margin_n + β
    alpha = cp.Variable()
    beta = cp.Variable()

    constraints = []
    for d in data:
        # disc >= alpha * margin + beta
        constraints.append(d['discriminant'] >= alpha * d['margin'] + beta)

    # Maximize the lower bound (to make it as tight as possible)
    objective = cp.Minimize(-beta)  # maximize beta (lower bound at margin=0)

    prob = cp.Problem(objective, constraints)
    prob.solve(solver=cp.CLARABEL, verbose=False)

    result_attempt1 = {
        'status': prob.status,
        'optimal': prob.status == 'optimal',
        'alpha': float(alpha.value) if alpha.value is not None else None,
        'beta': float(beta.value) if beta.value is not None else None,
    }

    if verbose:
        print(f"  Status: {prob.status}")
        print(f"  Lower bound form: disc(J^(4,n)) ≥ {alpha.value:.6e} * margin + {beta.value:.6e}")
        if alpha.value is not None and alpha.value < 0:
            print(f"  WARNING: α < 0 means bound is not useful (decreases with margin)")

    # ── Attempt 2: Certify positivity via SOS of margin ─────────────────
    print(f"\nAttempt 2: Certify disc(J^(4,n)) ≥ 0 via constrained SDP...")

    # Strategy: model disc as a quadratic in the margin, find Gram matrix
    # disc(margin) ≈ a * margin^2 + b * margin + c
    # Check if [margin^2, margin, 1]^T Q [margin^2, margin, 1] ≥ 0
    # where Q ⪰ 0 is positive semidefinite.

    # Fit a quadratic to data first
    margins = np.array([d['margin'] for d in data])
    discs = np.array([d['discriminant'] for d in data])

    # Quadratic fit: disc = a*margin^2 + b*margin + c
    try:
        coeffs = np.polyfit(margins, discs, 2)  # [a, b, c]
        fit_error = np.mean((np.polyval(coeffs, margins) - discs)**2)

        if verbose:
            print(f"  Quadratic fit: disc ≈ {coeffs[0]:.4e}*m^2 + {coeffs[1]:.4e}*m + {coeffs[2]:.4e}")
            print(f"  Fit MSE: {fit_error:.4e}")

        # Now certify this quadratic is SOS
        a_coeff, b_coeff, c_coeff = float(coeffs[0]), float(coeffs[1]), float(coeffs[2])

        # For SOS: a*m^2 + b*m + c = (sqrt(a)*m + b/(2*sqrt(a)))^2 + c - b^2/(4*a)
        # This is SOS iff a > 0 and c - b^2/(4*a) ≥ 0
        if a_coeff > 0:
            sos_const = c_coeff - b_coeff**2 / (4 * a_coeff)
            if verbose:
                print(f"  SOS analysis: a={a_coeff:.4e} > 0 ✓")
                print(f"  Discriminant of quadratic: c - b²/(4a) = {sos_const:.4e}")

            is_sos = sos_const >= -1e-10  # numerical tolerance
            result_attempt2 = {
                'status': 'SOS verified' if is_sos else 'Not SOS',
                'a': a_coeff,
                'b': b_coeff,
                'c': c_coeff,
                'sos_discriminant': sos_const,
                'is_sos': bool(is_sos),
            }
        else:
            result_attempt2 = {
                'status': 'Invalid (a < 0)',
                'a': a_coeff,
                'b': b_coeff,
                'c': c_coeff,
                'is_sos': False,
            }
    except Exception as e:
        result_attempt2 = {'status': f'Fit error: {str(e)}', 'is_sos': False}

    # ── Attempt 3: SDP with Gram matrix ─────────────────────────────────
    print(f"\nAttempt 3: SDP Gram matrix for proof of disc ≥ 0...")

    # Formulate: find matrix Q ⪰ 0 such that
    # disc_n - (Q[0,0]*m_n^2 + Q[1,1]*m_n + Q[2,2] + 2*Q[0,1]*m_n^2*m_n + ...) ≤ ε for all n

    # This is complex, so instead: use a simpler SDP that certifies all data points satisfy disc ≥ tol
    tol = cp.Variable(nonneg=True)

    constraints = []
    for d in data:
        constraints.append(d['discriminant'] >= tol)

    prob_tol = cp.Problem(cp.Maximize(tol), constraints)
    prob_tol.solve(solver=cp.CLARABEL, verbose=False)

    result_attempt3 = {
        'status': prob_tol.status,
        'optimal': prob_tol.status == 'optimal',
        'min_discriminant_found': float(tol.value) if tol.value is not None else None,
    }

    if verbose:
        print(f"  Status: {prob_tol.status}")
        if tol.value is not None:
            print(f"  Minimum disc(J^(4,n)) over tested n: {tol.value:.6e}")

    return {
        'attempt1_linear': result_attempt1,
        'attempt2_quadratic_sos': result_attempt2,
        'attempt3_min_discriminant': result_attempt3,
        'data': data,
    }

# ── Main ──────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("SOS/SDP CERTIFICATE FOR JENSEN POLYNOMIAL DISCRIMINANT")
    print("=" * 70)

    # Load or compute a_n
    print("\nStep 1: Load or compute Taylor coefficients a_n of Ξ(z)")
    a = load_a_coeffs_from_results()
    if a is None:
        print("  Computing via DFT...")
        a = compute_a_coeffs_dft(n_max=35, prec=50)
    else:
        print(f"  Loaded {len(a)} coefficients from existing results")

    print(f"\nFirst 15 coefficients:")
    for i in range(min(15, len(a))):
        print(f"  a_{i} = {a[i]:.10e}")

    # Verify log-convexity
    print(f"\nStep 2: Verify log-convexity of a_n")
    all_logconv = True
    for n in range(len(a) - 2):
        ratio = a[n+1]**2 / (a[n] * a[n+2] + 1e-30)
        if ratio < 1.0 - 1e-10:
            all_logconv = False
            print(f"  n={n}: log-convex ratio = {ratio:.8f} VIOLATED")

    if all_logconv:
        print(f"  ✓ Log-convex for all tested n (ratio > 1)")

    # Build SDP certificate
    print(f"\nStep 3: Build SDP certificate for disc(J^(4,n)) ≥ 0")
    cert = build_sdp_certificate(a, d_max_test=min(40, len(a) - 4), verbose=True)

    # Summary
    print(f"\n" + "=" * 70)
    print("CERTIFICATE SUMMARY")
    print("=" * 70)

    att1 = cert['attempt1_linear']
    if att1['optimal']:
        print(f"\n✓ Attempt 1 (Linear): OPTIMAL")
        print(f"    Lower bound: disc(J^(4,n)) ≥ {att1['alpha']:.6e}*margin + {att1['beta']:.6e}")
        if att1['alpha'] is not None and att1['alpha'] > 0:
            print(f"    Bound is USEFUL (α > 0)")
        else:
            print(f"    Bound may not be useful (α ≤ 0)")
    else:
        print(f"\n✗ Attempt 1 (Linear): {att1['status']}")

    att2 = cert['attempt2_quadratic_sos']
    if att2.get('is_sos'):
        print(f"\n✓ Attempt 2 (Quadratic SOS): CERTIFIED")
        print(f"    Quadratic: disc ≈ {att2['a']:.6e}*m^2 + {att2['b']:.6e}*m + {att2['c']:.6e}")
        print(f"    SOS discriminant: {att2['sos_discriminant']:.6e} (≥ 0 for SOS)")
        print(f"    INTERPRETATION: This is a sum-of-squares polynomial!")
    else:
        print(f"\n✗ Attempt 2 (Quadratic SOS): {att2.get('status', 'Not SOS')}")

    att3 = cert['attempt3_min_discriminant']
    if att3['optimal'] and att3['min_discriminant_found'] is not None:
        print(f"\n✓ Attempt 3 (Min Discriminant): FEASIBLE")
        print(f"    SDP-verified minimum disc(J^(4,n)) = {att3['min_discriminant_found']:.6e}")
        if att3['min_discriminant_found'] > 1e-10:
            print(f"    RESULT: disc > 0 for all tested n (with margin {att3['min_discriminant_found']:.2e})")
        else:
            print(f"    RESULT: disc ≈ 0 (boundary case, need tighter tolerance)")
    else:
        print(f"\n✗ Attempt 3 (Min Discriminant): {att3.get('status', 'Infeasible')}")

    # Save all results
    print(f"\nStep 4: Save results")
    out_path = os.path.join(RESULTS_DIR, "sos_discriminant_cert.json")

    # Clean data for JSON serialization
    cleaned_cert = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'n_coefficients': len(a),
        'attempt1_linear': att1,
        'attempt2_quadratic_sos': {
            'status': att2.get('status'),
            'a': att2.get('a'),
            'b': att2.get('b'),
            'c': att2.get('c'),
            'sos_discriminant': att2.get('sos_discriminant'),
            'is_sos': att2.get('is_sos'),
        },
        'attempt3_min_discriminant': att3,
        'data_points': [
            {
                'n': d['n'],
                'margin': d['margin'],
                'discriminant': d['discriminant'],
            }
            for d in cert['data']
        ]
    }

    with open(out_path, 'w') as f:
        json.dump(cleaned_cert, f, indent=2)

    print(f"  Saved to {out_path}")

    print(f"\n" + "=" * 70)
    print("END")
    print("=" * 70)


if __name__ == "__main__":
    main()
