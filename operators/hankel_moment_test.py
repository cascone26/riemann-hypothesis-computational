"""
Hankel Matrix / Stieltjes Moment Test for Li coefficients

If the sequence (lambda_n) is a Stieltjes moment sequence — i.e., there exists
a positive Borel measure mu on [0,inf) such that lambda_n = integral t^n dmu(t) —
then ALL lambda_n > 0 automatically for all n, unconditionally.

Test: the Hamburger moment condition requires all Hankel matrices
  H_N = [lambda_{i+j+1}]_{i,j=0..N-1}
to be positive definite (all eigenvalues > 0).

The Stieltjes condition additionally requires positivity of the shifted Hankel matrices
  H_N^{(1)} = [lambda_{i+j+2}]_{i,j=0..N-1}

If both families are PD for all N, then (lambda_n) is a Stieltjes moment sequence.
"""

import json
import numpy as np
import os

DATA_PATH = "/Users/scones/projects/riemann/results/li_unconditional_final.json"

def load_lambda():
    with open(DATA_PATH) as f:
        data = json.load(f)
    # data is a list of {n: ..., lambda: ...} or similar
    # Let's check structure
    if isinstance(data, list):
        if isinstance(data[0], dict):
            # Extract lambda values sorted by n
            pairs = [(d['n'], d['lambda']) for d in data]
            pairs.sort()
            ns = [p[0] for p in pairs]
            lams = [p[1] for p in pairs]
        else:
            lams = data
            ns = list(range(1, len(lams)+1))
    elif isinstance(data, dict):
        # Might be {"1": val, "2": val, ...} or {"lambdas": [...]}
        if "lambdas" in data:
            lams = data["lambdas"]
            ns = list(range(1, len(lams)+1))
        else:
            pairs = [(int(k), v) for k, v in data.items()]
            pairs.sort()
            ns = [p[0] for p in pairs]
            lams = [p[1] for p in pairs]
    return np.array(ns), np.array(lams, dtype=np.float64)

def build_hankel(lams, shift, N):
    """
    Build Hankel matrix H[i,j] = lams[i+j+shift] for i,j = 0..N-1
    shift=0: H[i,j] = lambda_{i+j+1} (1-indexed, lambda_1 at H[0,0])
    shift=1: H[i,j] = lambda_{i+j+2}
    """
    H = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            idx = i + j + shift
            if idx < len(lams):
                H[i, j] = lams[idx]
            else:
                H[i, j] = np.nan
    return H

def check_pd(H, label):
    """Check positive definiteness via eigenvalues."""
    # Remove any NaN rows/cols
    valid = ~np.any(np.isnan(H), axis=1)
    H_clean = H[valid][:, valid]
    if H_clean.shape[0] == 0:
        return None, "empty"
    eigvals = np.linalg.eigvalsh(H_clean)
    min_eig = eigvals.min()
    max_eig = eigvals.max()
    is_pd = bool(min_eig > 0)
    return is_pd, min_eig, max_eig

def main():
    print("Loading Li coefficients...")
    ns, lams = load_lambda()
    print(f"Loaded {len(lams)} values, n={ns[0]}..{ns[-1]}")
    print(f"lambda_1 = {lams[0]:.10f}")
    print(f"lambda_10 = {lams[9]:.10f}")
    print(f"lambda_100 = {lams[99]:.10f}")
    print(f"lambda_1000 = {lams[999]:.10f}" if len(lams) > 999 else "")
    print()

    # lambda array is 0-indexed: lams[0] = lambda_1, lams[k] = lambda_{k+1}

    print("="*60)
    print("HAMBURGER MOMENT TEST: Hankel matrices H_N[i,j] = lambda_{i+j+1}")
    print("(PD for all N ↔ (lambda_n) is a Hamburger moment sequence)")
    print("="*60)

    all_pd_hamburger = True
    for N in [2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]:
        if 2*N - 1 >= len(lams):
            print(f"N={N}: not enough data (need {2*N-1} values, have {len(lams)})")
            break
        H = build_hankel(lams, shift=0, N=N)
        result = check_pd(H, f"H_{N}")
        is_pd, min_eig, max_eig = result
        status = "PD ✓" if is_pd else "NOT PD ✗"
        print(f"  N={N:5d}: {status}  min_eig={min_eig:.6e}  max_eig={max_eig:.6e}")
        if not is_pd:
            all_pd_hamburger = False
            print(f"  *** HAMBURGER CONDITION FAILS AT N={N} ***")

    print()
    if all_pd_hamburger:
        print("RESULT: All tested Hankel matrices are POSITIVE DEFINITE")
        print("→ Strong evidence (lambda_n) is a Hamburger moment sequence")
        print("→ IF true for all N: exists positive measure μ with lambda_n = ∫t^n dμ(t)")
        print("→ This would imply ALL lambda_n > 0 unconditionally")
    else:
        print("RESULT: Hamburger condition FAILS — moment sequence approach ruled out")

    print()
    print("="*60)
    print("STIELTJES MOMENT TEST: Also check shifted H^(1)_N[i,j] = lambda_{i+j+2}")
    print("(Both families PD ↔ (lambda_n) is a Stieltjes moment sequence on [0,∞))")
    print("="*60)

    all_pd_stieltjes = True
    for N in [2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]:
        if 2*N >= len(lams):
            print(f"N={N}: not enough data")
            break
        H1 = build_hankel(lams, shift=1, N=N)
        result = check_pd(H1, f"H1_{N}")
        is_pd, min_eig, max_eig = result
        status = "PD ✓" if is_pd else "NOT PD ✗"
        print(f"  N={N:5d}: {status}  min_eig={min_eig:.6e}  max_eig={max_eig:.6e}")
        if not is_pd:
            all_pd_stieltjes = False
            print(f"  *** STIELTJES CONDITION FAILS AT N={N} ***")

    print()
    if all_pd_stieltjes and all_pd_hamburger:
        print("RESULT: BOTH families PD — (lambda_n) CONSISTENT WITH Stieltjes moment sequence")
        print("→ There may exist a positive measure μ on [0,∞) with lambda_n = ∫t^n dμ(t)")
        print("→ Finding μ analytically = proof that ALL lambda_n > 0 = PROOF OF RH")
    elif all_pd_hamburger and not all_pd_stieltjes:
        print("RESULT: Hamburger OK but Stieltjes FAILS — moment on R but not on [0,∞)")
    else:
        print("RESULT: Moment sequence approach ruled out")

    print()
    print("="*60)
    print("HAUSDORFF TEST: Finite differences (-1)^k Δ^k lambda_n ≥ 0?")
    print("(Would imply moment on [0,1] — probably fails given lambda_n ~ n log n)")
    print("="*60)

    # Just check first few finite differences
    for k in range(1, 6):
        # Δ^k lambda_0 (0-indexed: lambda_1)
        coeffs = np.array([(-1)**j * int(np.math.comb(k, j)) for j in range(k+1)])
        val = sum(coeffs[j] * lams[j] for j in range(k+1))
        sign = "≥0 ✓" if val >= 0 else "<0 ✗"
        print(f"  (-1)^{k} Δ^{k} λ_1 = {val:.6e}  {sign}")

    # Check at n=100
    print("  ... at n=100:")
    for k in range(1, 6):
        if 100 + k > len(lams): break
        coeffs = np.array([(-1)**j * int(np.math.comb(k, j)) for j in range(k+1)])
        val = sum(coeffs[j] * lams[99+j] for j in range(k+1))
        sign = "≥0 ✓" if val >= 0 else "<0 ✗"
        print(f"  (-1)^{k} Δ^{k} λ_100 = {val:.6e}  {sign}")

    print()
    print("="*60)
    print("MONOTONICITY OF DIFFERENCES (structure of growth):")
    print("="*60)
    delta1 = np.diff(lams[:1000])  # Δλ_n = λ_{n+1} - λ_n
    delta2 = np.diff(delta1)       # Δ²λ_n
    delta3 = np.diff(delta2)       # Δ³λ_n
    print(f"  Δλ_n  (n=1..999): min={delta1.min():.6e}, max={delta1.max():.6e}, all>0: {bool(delta1.min()>0)}")
    print(f"  Δ²λ_n (n=1..998): min={delta2.min():.6e}, max={delta2.max():.6e}, all>0: {bool(delta2.min()>0)}")
    print(f"  Δ³λ_n (n=1..997): min={delta3.min():.6e}, max={delta3.max():.6e}, all>0: {bool(delta3.min()>0)}")

    # Check sign pattern of higher differences
    print()
    print("  Sign pattern of Δ^k λ_1:")
    for k in range(1, 10):
        if k >= len(lams): break
        from math import comb
        val = sum((-1)**j * comb(k, j) * lams[j] for j in range(k+1))
        print(f"    k={k}: Δ^k λ_1 = {val:.4e}  sign={'+'if val>0 else '-'}")

if __name__ == "__main__":
    main()
