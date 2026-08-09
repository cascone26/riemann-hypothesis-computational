#!/usr/bin/env python3
"""§888 HONESTY GATE: SYMBOL APPROXIMATION VALIDITY.

The §886/§887 reframe (RH-corner ⟺ ∫|φ|²σ ⪰ 0) treats prime as a POINTWISE MULTIPLIER P(r).
TRUE prime is an off-diagonal (Toeplitz-like) OPERATOR built from cos(r log n) integrals.

CRITICAL QUESTION: is the error ‖Q − M_σ‖ (generalized/whitened) ≪ margin (frame VALID)
or ≳ margin (frame INVALID at razor scale)?

This script:
(1) Builds TRUE Q = arch + pole − prime in Fejér basis (copy assembly from coupling_probe.py)
(2) Builds DIAGONAL-SYMBOL operator M_σ: (M_σ)_ij = ∫ P_i(r)P_j(r) σ(r) dr
    where σ = Ω/(2π) − P_pointwise, P_pointwise = (1/π)Σ Λ(n)/√n cos(r log n)
(3) Compares:
    - Full spectral comparison (eigh Q vs M_σ against G)
    - Minimum generalized eigenvalue: min(Q) vs min(M_σ) in basis weighted by G^{-1/2}
    - Margin drift: |margin(Q) − margin(M_σ)|
    - Operator error: ‖Q − M_σ‖ (generalized whitened norm)
(4) VERDICT: is approximation error ≪ razor margin (~1.7e-9 at [14,44] 4W=4)?
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
np.set_printoptions(precision=4, suppress=True, linewidth=160)

def vonmangoldt(N):
    lam = np.zeros(N+1)
    isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True
            lp = np.log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
    return lam

NMAX = 3000
LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1)
LOGN = np.log(NS)
LAMW = LAM[2:] / np.sqrt(NS)

def build_weil_and_symbol(a, b, fourW, N, dr=0.004):
    """Build TRUE Weil matrices + symbol operator."""
    W = fourW / 4
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)

    # Fejér basis: P[i] = fej(r - c_i) + fej(r + c_i)
    def fej(x):
        y = W * x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz])**2
        return out

    P = np.array([fej(RGRID - c) + fej(RGRID + c) for c in CEN])

    # Gram matrix G (metric on the basis)
    G = (P @ P.T) * dr

    # Check conditioning
    gw_all, gv_all = eigh(G)
    cond_num = gw_all[-1] / (gw_all[0] + 1e-30)
    if cond_num > 1e12:
        print(f"    WARNING: G is ill-conditioned (κ={cond_num:.2e}), clamping small eigenvalues")
        gw_all = np.maximum(gw_all, 1e-10)

    # ARCH: ∫ P_i(r) P_j(r) Ω(r)/(2π) dr
    arch = (P * OMEGA[None, :]) @ P.T * dr / (2*np.pi)

    # PRIME: ∫ P_i(r) P_j(r) (1/π) Σ_{log n < 4W} Λ(n)/√n cos(r log n) dr
    # Using the assembly from coupling_probe.py
    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :] * P) @ cosu.T) * dr / (2*np.pi)
        prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    # POLE: rank-2 term from Riemann functional equation
    def fc(c):
        y = W * (0.5j - c)
        y2 = W * (0.5j + c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN])
    pole = 2.0 * np.real(np.outer(fp, fp))

    # TRUE Weil form: Q = arch + pole - prime
    Q = arch + pole - prime

    # Whitening: B from G eigendecomposition (use clamped eigenvalues)
    B = gv_all @ np.diag(1.0 / np.sqrt(gw_all)) @ gv_all.T

    # Whiten all matrices
    Qw = B @ Q @ B
    archw = B @ arch @ B
    primew = B @ prime @ B
    polew = B @ pole @ B

    # SYMBOL OPERATOR: build M_σ in the WHITENED basis
    # σ(r) = Ω(r)/(2π) - P_prime(r)
    # where P_prime(r) = (1/π) Σ_{log n < 4W} Λ(n)/√n cos(r log n)
    P_prime = (1.0/np.pi) * (lw[:, None] * np.cos(np.outer(ln, RGRID))).sum(axis=0)
    sigma = OMEGA / (2*np.pi) - P_prime

    # Construct M_σ_unwhitened first: (M_σ)_ij = ∫ P_i(r) P_j(r) σ(r) dr
    M_sigma_unwhitened = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            M_sigma_unwhitened[i, j] = (P[i] * P[j] * sigma).sum() * dr

    # Symmetrize (should already be symmetric, but be safe)
    M_sigma_unwhitened = 0.5 * (M_sigma_unwhitened + M_sigma_unwhitened.T)

    # Whiten
    M_sigma_w = B @ M_sigma_unwhitened @ B

    return Qw, M_sigma_w, polew, archw, primew, G, B, sigma, P, RGRID, CEN

def compare_operators(Qw, M_sigma_w, polew, archw, primew, window_name):
    """Compare Q and M_σ at the razor scale."""
    N = Qw.shape[0]

    # Check for NaN before eigendecomposition
    if np.any(~np.isfinite(Qw)) or np.any(~np.isfinite(M_sigma_w)):
        print(f"    SKIPPED: matrix contains NaN/Inf (ill-conditioning at this window)")
        return None

    # Eigenvalue decomposition (generalized against identity in whitened basis)
    eQ, VQ = eigh(Qw)
    eM, VM = eigh(M_sigma_w)
    eP, VP = eigh(polew)
    eA, VA = eigh(archw)

    # Minimum eigenvalues (the margins)
    margin_Q = eQ[0]
    margin_M = eM[0]
    margin_diff = abs(margin_Q - margin_M)

    # Operator error (Frobenius norm on whitened form is already generalized)
    diff = Qw - M_sigma_w
    error_frob = np.linalg.norm(diff, 'fro')
    error_op = np.linalg.norm(diff, 2)  # spectral norm

    # Relative error scaled to margin
    rel_to_margin_frob = error_frob / max(abs(margin_Q), 1e-15)
    rel_to_margin_op = error_op / max(abs(margin_Q), 1e-15)

    # Full spectral comparison: how many eigvals differ significantly?
    # Align eigenvectors by dot product and check eigenvalue agreement
    overlaps = np.abs(VQ.T @ VM)  # inner products between eigvecs
    sorted_eM_by_overlap = []
    for i in range(N):
        jbest = np.argmax(overlaps[i, :])
        sorted_eM_by_overlap.append((eM[jbest], jbest))

    eig_diffs = [abs(eQ[i] - em) for i, (em, _) in enumerate(sorted_eM_by_overlap)]

    print(f"\n{'='*110}")
    print(f"§888 SYMBOL APPROXIMATION ERROR — {window_name}")
    print(f"{'='*110}")
    print(f"TRUE Weil form Q = arch + pole - prime")
    print(f"APPROX M_σ: compression of multiply-by-σ(r)=Ω/2π−P_prime onto Fejér basis")
    print(f"{'─'*110}")

    print(f"\nMARGINS (min eigenvalues):")
    print(f"  margin(Q)      = {margin_Q:+.6e}")
    print(f"  margin(M_σ)    = {margin_M:+.6e}")
    print(f"  |margin(Q)−M_σ| = {margin_diff:+.6e}  (absolute drift)")

    print(f"\nOPERATOR NORMS (on whitened basis):")
    print(f"  ‖Q − M_σ‖_Frob = {error_frob:.6e}")
    print(f"  ‖Q − M_σ‖_spec = {error_op:.6e}")

    print(f"\nRELATIVE TO MARGIN:")
    print(f"  ‖Q − M_σ‖_Frob / |margin| = {rel_to_margin_frob:.6e}")
    print(f"  ‖Q − M_σ‖_spec / |margin| = {rel_to_margin_op:.6e}")

    print(f"\nEIGENVALUE SPECTRUM (first 6 eigs):")
    print(f"  Q eigenvals:  {eQ[:6]}")
    print(f"  M_σ eigenvals: {eM[:6]}")
    print(f"  eig diffs (abs): {np.array(eig_diffs[:6])}")

    print(f"\nCOMPONENT NORMS (whitened basis):")
    print(f"  ‖arch_w‖    = {np.linalg.norm(archw, 'fro'):.6e}")
    print(f"  ‖prime_w‖   = {np.linalg.norm(primew, 'fro'):.6e}")
    print(f"  ‖pole_w‖    = {np.linalg.norm(polew, 'fro'):.6e}")
    print(f"  ‖Q_w‖       = {np.linalg.norm(Qw, 'fro'):.6e}")
    print(f"  ‖M_σ_w‖     = {np.linalg.norm(M_sigma_w, 'fro'):.6e}")

    # Verdict
    print(f"\n{'─'*110}")
    if error_op < 0.1 * abs(margin_Q):
        verdict = "✓ VALID: error ≪ margin (frame is rigorous at razor scale)"
    elif error_op < abs(margin_Q):
        verdict = "⚠ MARGINAL: error is ~margin size (frame is FRAGILE at razor scale)"
    else:
        verdict = "✗ INVALID: error ≳ margin (diagonal-symbol approximation BREAKS DOWN)"

    print(f"VERDICT: {verdict}")
    print(f"{'='*110}")

    return {
        'margin_Q': margin_Q,
        'margin_M': margin_M,
        'margin_diff': margin_diff,
        'error_op': error_op,
        'error_frob': error_frob,
        'rel_error_op': rel_to_margin_op,
        'rel_error_frob': rel_to_margin_frob,
        'verdict_ok': error_op < 0.1 * abs(margin_Q)
    }

# Run at the key windows
print("\n" + "="*110)
print("TESTING SYMBOL APPROXIMATION VALIDITY ACROSS WINDOWS")
print("="*110)

results = {}
for (a, b, fourW, N) in [(14, 44, 4.0, 17), (14, 44, 4.0, 13), (60, 90, 4.0, 30)]:
    window_name = f"[{a},{b}] 4W={fourW} N={N}"
    Qw, M_sigma_w, polew, archw, primew, G, B, sigma, P, RGRID, CEN = build_weil_and_symbol(a, b, fourW, N)
    res = compare_operators(Qw, M_sigma_w, polew, archw, primew, window_name)
    if res is not None:
        results[window_name] = res

print("\n" + "="*110)
print("SUMMARY TABLE")
print("="*110)
if results:
    print(f"{'Window':>25} {'margin(Q)':>13} {'‖Q−M_σ‖':>13} {'rel error':>12} {'Valid?':>8}")
    for window_name, res in results.items():
        valid = "YES ✓" if res['verdict_ok'] else "NO ✗"
        print(f"{window_name:>25} {res['margin_Q']:>+13.4e} {res['error_op']:>13.4e} {res['rel_error_op']:>12.4e} {valid:>8}")
    print("="*110)

    # Final honesty assessment
    all_valid = all(r['verdict_ok'] for r in results.values())
    if all_valid:
        print("\n✓ FRAME VALID AT ALL TESTED WINDOWS: symbol approx error stays ≪ margin across tests.")
        print("  Reframe (RH-corner ⟺ ∫|φ|²σ ⪰ 0) is RIGOROUS at tested scale.")
    else:
        print("\n✗ FRAME COMPROMISED: symbol approx error EXCEEDS margin at some window(s).")
        print("  Reframe relies on invalid pointwise-symbol approximation. Need full operator analysis.")
else:
    print("No valid results (all windows ill-conditioned or had NaN).")
print("="*110)
