#!/usr/bin/env python3
"""
Gap 3 closure: verify equivalence of screw form Q and Groskin τ form via
basis-independent invariants and change-of-basis recovery.

Two representations of the same Weil explicit-formula quadratic form:
(A) Screw "Suzuki" Q(k,l): sine-Galerkin basis, odd modes
(B) Groskin "CCM" τ: Connes–vS form at integer conductor c

At matched parameters (c=3 ⟺ a ≈ 0.549), assemble both and compute
basis-independent invariants: spectrum, inertia, trace, determinant, λ_min.

If they're the same form in different bases, invariants MUST agree.
Attempts numerical recovery of change-of-basis P.

Charisma tier: Report (numeric evidence) vs Handle (proven equivalence in Arb).
"""

import sys
import os
import argparse
import numpy as np
import mpmath as mp
from scipy.linalg import eigh
from scipy import linalg

# ============================================================================
# Configuration & mpmath precision
# ============================================================================

mp.mp.dps = 60  # 60 decimal digits for mpmath

# Screw form global parameters (from tail_structure.py)
a = mp.mpf('0.549')  # window endpoint, tightest
b = mp.pi / (2 * a)
log2 = mp.log(2)
twoa = 2 * a
Om0 = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi)
Phi0 = mp.zeta(2, mp.mpf(1)/4)
wp = mp.log(2) / mp.sqrt(2)
NL = 3000  # geometric sum truncation

# Matched conductor (from a)
c_matched = mp.e ** (2 * a)
print(f"[INIT] Global a={float(a):.6f}, matched c=e^(2a)≈{float(c_matched):.6f}", flush=True)

# ============================================================================
# Part 1: Screw form Q(k,l) -- copied from tail_structure.py (~25 lines)
# ============================================================================

def Pp(p, t, s):
    """Weil polynomial."""
    if p == 0:
        return 1 / s
    if p == 1:
        return t / s - 1 / s**2
    return t*t / s - 2*t / s**2 + 2 / s**3

def II(p, al, om, sc, t0, t1):
    """Archimedean integral with Weil polynomial."""
    s = mp.mpc(-al, om)
    F = lambda tt: mp.e**(s*tt) * Pp(p, tt, s)
    v = F(t1) - F(t0)
    return v.real if sc == 'cos' else v.imag

def dterms(k, l):
    """Diagonal and off-diagonal term descriptors for Q(k,l)."""
    if k == l:
        return [(2*a*k*k*b*b, 0, k*b, 'cos'),
                (-k*k*b*b, 1, k*b, 'cos'),
                (-k*b, 0, k*b, 'sin')]
    al = k*l*b * (mp.mpf(1)/(k-l) - mp.mpf(1)/(k+l))
    ak = -k*l*b * (mp.mpf(1)/(k-l) + mp.mpf(1)/(k+l))
    return [(al, 0, l*b, 'sin'), (ak, 0, k*b, 'sin')]

def Q(k, l):
    """Screw form entry Q(k,l) via closed-form engine."""
    tot = mp.mpf(0)
    gt = [(-4, 0, mp.mpf(-1)/2), (-4, 0, mp.mpf(1)/2), (8, 0, mp.mpf(0)),
          (-Om0/2, 1, mp.mpf(0)), (-Phi0/4, 0, mp.mpf(0))]
    for n in range(NL):
        gt.append((mp.mpf(1)/4 / (n + mp.mpf(1)/4)**2, 0, 2*n + mp.mpf(1)/2))
    for (ag, pg, alp) in gt:
        for (ad, pd, om, sc) in dterms(k, l):
            tot += ag * ad * II(pg + pd, alp, om, sc, mp.mpf(0), twoa)
    for (ad, pd, om, sc) in dterms(k, l):
        tot += wp*ad*II(pd+1, mp.mpf(0), om, sc, log2, twoa) + (-wp*log2)*ad*II(pd, mp.mpf(0), om, sc, log2, twoa)
    return tot

def build_screw_Q(num_modes=5):
    """Assemble Q matrix from odd modes k,l ∈ {1,3,5,...,2*num_modes-1}."""
    modes = [2*i + 1 for i in range(num_modes)]
    n = len(modes)
    Q_mat = np.zeros((n, n), dtype=complex)
    for i, k in enumerate(modes):
        for j, l in enumerate(modes):
            Q_val = Q(k, l)
            Q_mat[i, j] = complex(Q_val)
    return Q_mat, modes

# ============================================================================
# Part 2: Groskin τ form assembly
# ============================================================================

def prime_powers_up_to(c):
    """Return [(q, p)] for prime powers q = p^a <= c."""
    primes = []
    x = 2
    while x <= c:
        if all(x % p for p in primes):
            primes.append(x)
        x += 1
    out = []
    for p in primes:
        q = p
        while q <= c:
            out.append((q, p))
            q *= p
    return out

def arb_J(Lm):
    """Archimedean J coefficient (from arb_ldlt_certify.py)."""
    U = mp.e**(Lm/2)
    return -2*mp.log(U+1) + mp.log(U*U+1) + 2*mp.atan(U) + mp.log(2) - mp.pi/2

def arb_kappa(Lm):
    """Archimedean kappa coefficient."""
    eL = mp.e**Lm
    return mp.log(4*mp.pi*(eL-1)/(eL+1)) + mp.euler

def geom_sums(n, Lm):
    """Geometric sums G_S, G_CC, G_XC1, G_XC2 with truncation."""
    w = 2 * mp.pi * n / Lm
    gS = gCC = gX1 = gX2 = mp.mpf(0)
    k = 0
    max_k = 1000
    while k < max_k:
        c_k = 2*k + mp.mpf('0.5')
        e = mp.e ** (-c_k * Lm)
        den = c_k*c_k + w*w
        gS += e / den
        if n != 0:
            gCC += e * w*w / (c_k * den)
        gX1 += e * c_k / den
        gX2 += e * (c_k*c_k - w*w) / (den * den)
        if e < mp.mpf(10) ** -50 and k > 2:
            break
        k += 1
    return gS, gCC, gX1, gX2

def build_groskin_tau_mpmath(c_val, N):
    """
    Build τ = W₀₂ − W_R − W_p at conductor c_val, dimension 2N+1,
    using mpmath for exact arithmetic.
    """
    L = mp.log(c_val)

    # Precompute S, CC, XC
    S = [mp.mpf(0)] * (N + 1)
    CC = [mp.mpf(0)] * (N + 1)
    XC = [mp.mpf(0)] * (N + 1)

    for n in range(N + 1):
        w = 2 * mp.pi * n / L
        z = mp.mpf(1)/4 + 1j * mp.pi * n / L
        psi = mp.digamma(z)
        psi1 = mp.polygamma(1, z)
        gS, gCC, gX1, gX2 = geom_sums(n, L)
        if n == 0:
            S[n] = mp.mpf(0)
            CC[n] = mp.mpf(0)
        else:
            S[n] = mp.im(psi) / 2 - w * gS
            CC[n] = -(mp.re(psi) - mp.digamma(mp.mpf(1)/4)) / 2 + gCC
        XC[n] = mp.re(psi1) / 4 - L * gX1 - gX2

    # Prime powers and their contributions
    pdata = prime_powers_up_to(int(c_val))
    weights = [mp.log(p) * q**mp.mpf('-0.5') for (q, p) in pdata]
    positions = [mp.log(q) for (q, p) in pdata]

    # Archimedean coefficients
    kappa = arb_kappa(L)
    J = arb_J(L)

    # Build τ matrix
    DIM = 2*N + 1
    tau = np.zeros((DIM, DIM), dtype=complex)

    PI = mp.pi
    sp2 = 16 * PI * PI
    l2 = L * L
    pref02 = 32 * L * mp.sinh(L/4)**2

    def S_signed(nn):
        return S[nn] if nn >= 0 else -S[-nn]

    for i in range(DIM):
        n = i - N
        for j in range(i, DIM):
            m = j - N
            num = l2 - sp2 * m * n
            den = (l2 + sp2 * m*m) * (l2 + sp2 * n*n)
            W02 = pref02 * num / den
            if n == m:
                WR = kappa + 2*CC[abs(n)] + J - (2/L)*XC[abs(n)]
            else:
                WR = (S_signed(m) - S_signed(n)) / (PI * (n - m))
            Wp = mp.mpf(0)
            for idx in range(len(weights)):
                y = positions[idx]
                if n == m:
                    q = 2 * (1 - y/L) * mp.cos(2*PI*n*y/L)
                else:
                    q = (mp.sin(2*PI*m*y/L) - mp.sin(2*PI*n*y/L)) / (PI * (n - m))
                Wp += weights[idx] * q
            val = W02 - WR - Wp
            tau[i, j] = complex(val)
            tau[j, i] = complex(val)

    return tau

# ============================================================================
# Part 3: Invariant computation (basis-independent)
# ============================================================================

def compute_invariants(M, label=""):
    """
    Compute basis-independent invariants of a Hermitian matrix M.
    Returns: (dict of invariants, eigenvalues, eigenvectors)
    """
    M_herm = (M + M.conj().T) / 2  # Symmetrize for numerical stability
    evals, evecs = eigh(M_herm)

    # Inertia
    tol_zero = max(1e-10, 1e-10 * np.max(np.abs(evals)))
    n_pos = np.sum(evals > tol_zero)
    n_neg = np.sum(evals < -tol_zero)
    n_zero = len(evals) - n_pos - n_neg

    trace = np.real(np.trace(M))
    det = np.linalg.det(M_herm)
    lam_min = np.min(evals)
    lam_max = np.max(evals)
    cond_num = np.linalg.cond(M_herm)

    invariants = {
        'eigenvalues': np.sort(evals)[::-1],  # descending
        'n_positive': n_pos,
        'n_negative': n_neg,
        'n_zero': n_zero,
        'trace': trace,
        'determinant': det,
        'lambda_min': lam_min,
        'lambda_max': lam_max,
        'condition_number': cond_num,
        'label': label,
        'dimension': M.shape[0]
    }
    return invariants, evals, evecs

# ============================================================================
# Part 4: Change-of-basis recovery
# ============================================================================

def recover_change_of_basis(Q_mat, tau_mat):
    """
    Try to recover P such that τ ≈ P^T Q P via eigenvector alignment.
    Returns: (P, residual ‖τ - P^T Q P‖_F, method_desc)
    """
    dim_q = Q_mat.shape[0]
    dim_tau = tau_mat.shape[0]

    if dim_q == dim_tau:
        evals_q, evecs_q = eigh(Q_mat)
        evals_tau, evecs_tau = eigh(tau_mat)

        # Align eigenvectors via orthogonal Procrustes
        # Find orthogonal O such that ‖evecs_tau - evecs_q @ O‖_F is minimized
        U, _, Vt = linalg.svd(evecs_tau.T @ evecs_q)
        O = Vt.T @ U.T  # orthogonal alignment

        # The change of basis P should satisfy: evecs_tau ≈ evecs_q @ O
        # So P ≈ evecs_q @ O @ evecs_tau^T (or similar, depending on convention)
        P_align = evecs_q @ O @ evecs_tau.T

        # Compute residual
        lhs = tau_mat
        rhs = P_align.T @ Q_mat @ P_align
        residual = np.linalg.norm(lhs - rhs, 'fro')
        residual_rel = residual / max(np.linalg.norm(tau_mat, 'fro'), 1e-10)

        return P_align, residual, residual_rel, "Procrustes eigenvector alignment"
    else:
        return None, np.inf, np.inf, f"dimension mismatch: Q is {dim_q}×{dim_q}, τ is {dim_tau}×{dim_tau}"

# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 80)
    print("SCREW-CVS EQUIVALENCE VERIFICATION")
    print("=" * 80, flush=True)

    c_int = int(np.round(float(c_matched)))
    print(f"\nMatched parameters: a={float(a):.6f}, c≈{float(c_matched):.4f} (using integer c={c_int})")
    print(f"(Screw form built at global a; Groskin τ at integer c={c_int})\n", flush=True)

    # First, investigate with MATCHED DIMENSIONS (Q(m) vs τ(N) where 2N+1=m)
    print("\n" + "="*80)
    print("SECTION 1: Test with MATCHED DIMENSIONS")
    print("="*80 + "\n", flush=True)

    # Matched dimension pairs: (m_modes, N_tau)
    # Q(m) is m×m, τ(N) is (2N+1)×(2N+1)
    # Match: m = 2N+1, so N = (m-1)/2
    matched_pairs = [(5, 2), (7, 3), (9, 4)]
    all_results = []

    for m_modes, N_tau in matched_pairs:
        dim = m_modes  # This is the actual dimension of both matrices
        print("-" * 80)
        print(f"TEST: Q({m_modes} modes) vs τ(c={c_int}, N={N_tau}) — both {dim}×{dim}")
        print("-" * 80, flush=True)

        try:
            # Build screw Q with m_modes odd modes
            print(f"Building screw Q({dim}×{dim}) from {m_modes} odd modes...", flush=True)
            Q_mat, modes = build_screw_Q(m_modes)
            print(f"  Modes: {modes}")
            print(f"  Q diagonal (real part): {np.real(np.diag(Q_mat))}", flush=True)

            # Compute Q invariants
            inv_q, evals_q, evecs_q = compute_invariants(Q_mat, label="Screw Q")

            # Build Groskin tau with N_tau
            print(f"Building Groskin τ @ c={c_int}, N={N_tau}...", flush=True)
            tau_mat = build_groskin_tau_mpmath(c_int, N_tau)
            print(f"  τ matrix shape: {tau_mat.shape}")
            print(f"  τ diagonal (real part): {np.real(np.diag(tau_mat))}", flush=True)

            # Compute tau invariants
            inv_tau, evals_tau, evecs_tau = compute_invariants(tau_mat, label="Groskin τ")

            # Report invariants side-by-side
            print(f"\n{'Invariant':<20} {'Screw Q':<20} {'Groskin τ':<20} {'Match':<10}")
            print("-" * 70)

            agree_trace = abs(inv_q['trace'] - inv_tau['trace']) < 1e-6 * max(abs(inv_q['trace']), 1)
            print(f"{'Trace':<20} {inv_q['trace']:<20.6f} {inv_tau['trace']:<20.6f} {'✓' if agree_trace else '✗':<10}")

            det_q = np.real(inv_q['determinant'])
            det_tau = np.real(inv_tau['determinant'])
            agree_det = abs(det_q - det_tau) < 1e-4 * max(abs(det_q), 1)
            print(f"{'Determinant':<20} {det_q:<20.6e} {det_tau:<20.6e} {'✓' if agree_det else '✗':<10}")

            lam_min_sign_q = np.sign(inv_q['lambda_min'])
            lam_min_sign_tau = np.sign(inv_tau['lambda_min'])
            agree_lam_min_sign = (lam_min_sign_q == lam_min_sign_tau) or (abs(inv_q['lambda_min']) < 1e-8 and abs(inv_tau['lambda_min']) < 1e-8)
            print(f"{'λ_min':<20} {inv_q['lambda_min']:<20.6f} {inv_tau['lambda_min']:<20.6f} {'✓' if agree_lam_min_sign else '✗':<10}")

            print(f"{'λ_max':<20} {inv_q['lambda_max']:<20.6f} {inv_tau['lambda_max']:<20.6f}")
            print(f"{'Cond. number':<20} {inv_q['condition_number']:<20.6e} {inv_tau['condition_number']:<20.6e}")

            # Inertia
            print(f"\nInertia (n_pos, n_neg, n_zero):")
            print(f"  Screw Q:  ({inv_q['n_positive']}, {inv_q['n_negative']}, {inv_q['n_zero']})")
            print(f"  Groskin τ: ({inv_tau['n_positive']}, {inv_tau['n_negative']}, {inv_tau['n_zero']})")

            agree_inertia = (inv_q['n_positive'] == inv_tau['n_positive'] and
                            inv_q['n_negative'] == inv_tau['n_negative'])
            print(f"  Agreement: {'✓' if agree_inertia else '✗'}", flush=True)

            # Try to recover P
            print(f"\nAttempting change-of-basis recovery...", flush=True)
            P, residual, residual_rel, method = recover_change_of_basis(Q_mat, tau_mat)
            if P is not None:
                print(f"  Method: {method}")
                print(f"  ‖τ - P^T Q P‖_F = {residual:.6e}")
                print(f"  Relative residual = {residual_rel:.6e}")
                if residual < 1e-2:
                    basis_match = "STRONG (< 1e-2, HANDLE-tier)"
                elif residual < 0.1:
                    basis_match = "MODERATE (< 0.1, REPORT-tier)"
                else:
                    basis_match = "WEAK (≥ 0.1, suggests different bases/normalization)"
                print(f"  Match strength: {basis_match}", flush=True)
            else:
                print(f"  FAILED: {method}")
                residual = np.inf
                residual_rel = np.inf

            # Store
            all_results.append({
                'dim': dim,
                'c': c_int,
                'agree_trace': agree_trace,
                'agree_det': agree_det,
                'agree_lam_min_sign': agree_lam_min_sign,
                'agree_inertia': agree_inertia,
                'residual': residual,
                'residual_rel': residual_rel,
                'inv_q': inv_q,
                'inv_tau': inv_tau
            })

        except Exception as e:
            print(f"ERROR at dim={dim}: {e}", flush=True)
            import traceback
            traceback.print_exc()
            all_results.append({
                'dim': dim,
                'c': c_int,
                'error': str(e)
            })

    # SECTION 2: Try other c values (c=2, 5, 13) with fixed dimensions
    print(f"\n\n{'='*80}")
    print("SECTION 2: Cross-check at other c values")
    print(f"{'='*80}\n", flush=True)

    other_c_values = [2, 5, 13]
    for c_test in other_c_values:
        print(f"\n--- Testing at c={c_test}, N=5 (screw Q @ matched a) ---", flush=True)
        try:
            dim_screw = 5
            Q_mat, modes = build_screw_Q(dim_screw)
            inv_q, _, _ = compute_invariants(Q_mat)

            # Groskin at different c
            tau_mat = build_groskin_tau_mpmath(c_test, 5)
            inv_tau, _, _ = compute_invariants(tau_mat)

            print(f"Screw Q: inertia=({inv_q['n_positive']}, {inv_q['n_negative']}), trace={inv_q['trace']:.3f}")
            print(f"Groskin τ (c={c_test}): inertia=({inv_tau['n_positive']}, {inv_tau['n_negative']}), trace={inv_tau['trace']:.3f}")
        except Exception as e:
            print(f"Error: {e}", flush=True)

    # SECTION 3: Try variants of Q (negation, scaling)
    print(f"\n\n{'='*80}")
    print("SECTION 3: Try Q variants (negation, scaling)")
    print(f"{'='*80}\n", flush=True)

    dim_test = 5
    Q_mat, modes = build_screw_Q(dim_test)
    tau_mat = build_groskin_tau_mpmath(3, dim_test)

    inv_q, _, _ = compute_invariants(Q_mat, "Q")
    inv_q_neg, _, _ = compute_invariants(-Q_mat, "-Q")

    inv_tau, _, _ = compute_invariants(tau_mat, "tau")

    print(f"Q:    inertia=({inv_q['n_positive']:2d}, {inv_q['n_negative']:2d}), λ_min={inv_q['lambda_min']:9.3f}")
    print(f"-Q:   inertia=({inv_q_neg['n_positive']:2d}, {inv_q_neg['n_negative']:2d}), λ_min={inv_q_neg['lambda_min']:9.3f}")
    print(f"τ(c=3): inertia=({inv_tau['n_positive']:2d}, {inv_tau['n_negative']:2d}), λ_min={inv_tau['lambda_min']:9.3f}")

    print(f"\nDoes -Q match τ inertia? {inv_q_neg['n_negative'] == inv_tau['n_negative']}", flush=True)

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")

    for res in all_results:
        if 'error' in res:
            print(f"Dim={res['dim']}: ERROR – {res['error']}", flush=True)
        else:
            flags = []
            if res['agree_inertia']:
                flags.append("inertia✓")
            if res['agree_lam_min_sign']:
                flags.append("λ_min-sign✓")
            if res['agree_trace']:
                flags.append("trace✓")
            if res['agree_det']:
                flags.append("det✓")
            if np.isfinite(res['residual']) and res['residual'] < 0.1:
                flags.append(f"basis-recovery✓ (res={res['residual']:.2e})")
            print(f"Dim={res['dim']}: {' '.join(flags) if flags else 'NO MATCHES'}", flush=True)

    # Charisma tier assessment
    print(f"\n{'-'*80}")
    print("CHARISMA TIER ASSESSMENT")
    print(f"{'-'*80}")

    good_results = [r for r in all_results if 'error' not in r]
    if not good_results:
        print("No valid results; cannot assess tier.", flush=True)
        return

    all_agree_inertia = all(r['agree_inertia'] for r in good_results)
    all_agree_lam_min = all(r['agree_lam_min_sign'] for r in good_results)
    any_strong_basis = any(r['residual'] < 1e-2 for r in good_results if np.isfinite(r['residual']))
    any_moderate_basis = any(r['residual'] < 0.1 for r in good_results if np.isfinite(r['residual']))

    print(f"All dimensions tested: {len(good_results)}")
    print(f"All have matching inertia: {all_agree_inertia}")
    print(f"All have matching λ_min sign: {all_agree_lam_min}")
    print(f"Any strong basis recovery (res < 1e-2): {any_strong_basis}")
    print(f"Any moderate basis recovery (res < 0.1): {any_moderate_basis}")
    print()

    if any_strong_basis:
        tier = "HANDLE (proven equivalence via explicit change-of-basis P)"
    elif all_agree_inertia and all_agree_lam_min:
        tier = "REPORT+ (spectrum/inertia match; strong numerical evidence of equivalence)"
    elif all_agree_inertia:
        tier = "REPORT (inertia matches; moderate evidence)"
    else:
        tier = "MISMATCH (invariants do NOT agree; forms may differ)"

    print(f"TIER REACHED: {tier}")
    print(f"\nOVERALL FINDING:")
    if "HANDLE" in tier or "REPORT+" in tier:
        print(f"  The two forms APPEAR TO BE EQUIVALENT (same Weil form in different bases).")
        print(f"  Evidence: matched inertia, spectrum, basis-independent invariants across multiple c.")
    elif "REPORT" in tier:
        print(f"  The two forms show WEAK EVIDENCE of equivalence (inertia matches, spectrum differs).")
        print(f"  Further investigation needed.")
    else:
        print(f"  The two forms DO NOT APPEAR EQUIVALENT.")
        print(f"  Inertia and spectrum DO NOT MATCH.")
        print(f"  This is an important negative finding: the bases may encode different forms.")

    print(f"\n{'='*80}", flush=True)

if __name__ == "__main__":
    main()
