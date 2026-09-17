#!/usr/bin/env python3
"""Candidate (B): Augment pole⊥-PSD form toward full m-function.

Core task: Given M_σ|pole⊥ ⪰ 0 (verified HANDLE §976), and the pole direction v
with M_σ v ≈ μv (μ < 0, single negative eigenvalue), derive the PRECISE linear-algebra
condition for "restoring v produces exactly κ=1 negative eigenvalue" — not 0, not 2.

This is a classical bordered matrix problem:
  M_σ = [[A, b], [b^T, c]]  where A ≈ M_σ|v⊥ ⪰ 0 (PSD)
                            b ≈ M_σ(rest) coupling to v
                            c ≈ <v, M_σ v> = μ (negative)

The exact condition for exactly one negative eigenvalue: Schur complement tells us
  κ(M_σ) = κ(A) + κ(c - b^T A^+ b) = 0 + κ(c - b^T A^+ b)

where A^+ is the pseudoinverse (or Moore-Penrose). We need:
  κ(c - b^T A^+ b) = 1  ⟺  c - b^T A^+ b < 0  (exactly one negative)

The challenge: this is still circular (depends on knowing A's structure). But we can
phrase it WITHOUT the pseudoinverse via the generalized eigenvalue reformulation (§980).

Implementation:
1. Compute M_σ from the exact Suzuki form
2. Extract v (principal eigenvector of pole term)
3. Compute M_σ|v⊥ explicitly
4. Derive the Schur-complement based condition
5. Test it: does it correctly predict κ=1 at true primes and fail at perturbed primes?
6. Derive whether this can be restated in m-function language (Nevanlinna pole)
"""
import numpy as np
from scipy.linalg import eigh
import sys, os
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
os.chdir('/Users/scones/projects/riemann/operators')
from suzuki_exact_Q import build_Q, vonmangoldt_upto

def perturb_primes(K, delta, seed=None):
    """Perturb log-frequencies to create RH-violating form."""
    if seed is not None:
        np.random.seed(seed)
    Lam = vonmangoldt_upto(int(np.floor(np.exp(2*2.0)))+1)
    primes = np.arange(2, len(Lam))[Lam > 0]
    log_primes_true = np.log(primes)
    log_primes_pert = log_primes_true * (1 + delta * np.random.uniform(-1, 1, len(primes)))
    return log_primes_true, log_primes_pert

def build_M_sigma(a, K, NG=320):
    """Extract M_σ = arch - prime from the exact Suzuki form.

    In the exact Suzuki code:
      Q = Qpole + Qcore + Qtail + Qprime
    where Qprime = Σ(-Λ(n)/√n * Sn), which is NEGATIVE.

    The Weil form is: Q = pole + arch - prime
    where arch = Qcore + Qtail, and prime = -Qprime (the positive magnitude).

    Therefore: M_σ = arch - prime = (core + tail) - (-Qprime) = core + tail + Qprime
    """
    Q, pieces = build_Q(a, K, NG, return_pieces=True)
    arch = pieces['core'] + pieces['tail']
    # pieces['prime'] is already negative (the code accumulates -Λ(n)... terms)
    # So M_σ = arch - prime_magnitude = arch - (-pieces['prime']) = arch + pieces['prime']
    M_sigma = arch + pieces['prime']
    return M_sigma, arch, pieces['prime'], pieces

def extract_pole_principal_direction(Qpole):
    """Extract the principal eigenvector of the pole term (rank-2 symmetric).

    The pole term has rank 2, so we extract its principal (largest) direction.
    """
    eigvals, eigvecs = eigh(Qpole)
    # Last eigenvector is the principal direction (largest eigenvalue)
    v_pole = eigvecs[:, -1]
    return v_pole / np.linalg.norm(v_pole), eigvals[-1]

def extract_negative_direction(M_sigma):
    """Extract the principal (most negative) eigenvector of M_σ."""
    eigvals, eigvecs = eigh(M_sigma)
    # First eigenvector corresponds to the most negative eigenvalue
    v_neg = eigvecs[:, 0]
    return v_neg / np.linalg.norm(v_neg), eigvals[0]

def restrict_to_orthogonal_complement(M, v):
    """Restrict matrix M to the orthogonal complement of v."""
    v_norm = v / np.linalg.norm(v)
    P = np.eye(len(v)) - np.outer(v_norm, v_norm)  # Projection onto v⊥
    M_restricted = P @ M @ P
    return M_restricted, P

def test_bordering_condition(a=1.6, K=120, NG=320, delta_perturb=None):
    """
    Core test: compute M_σ, extract pole direction v, restrict to v⊥,
    and test the Schur-complement based condition for κ=1.
    """
    print(f"\n=== Bordering Analysis: a={a}, K={K} ===")

    # Build M_σ
    M_sigma, arch, prime, pieces = build_M_sigma(a, K, NG)

    # Full spectral analysis
    eigs_full, eigvecs_full = eigh(M_sigma)
    kappa_full = int((eigs_full < -1e-10).sum())
    mu_min = eigs_full[0]
    print(f"M_σ full: κ={kappa_full}, λ_min={mu_min:.4e}, λ_2={eigs_full[1]:.4e}")

    # Extract pole principal direction
    v_pole, eig_pole = extract_pole_principal_direction(pieces['pole'])

    # Extract M_σ's negative direction
    v_neg, eig_neg = extract_negative_direction(M_sigma)

    # Alignment
    alignment = abs(np.dot(v_pole, v_neg))
    print(f"Pole principal-negative alignment: cos = {alignment:.6f}")

    # Restrict M_σ to v⊥
    M_perp, P = restrict_to_orthogonal_complement(M_sigma, v_pole)
    eigs_perp = np.linalg.eigvalsh(M_perp)
    kappa_perp = int((eigs_perp < -1e-10).sum())
    print(f"M_σ|v⊥: κ={kappa_perp}, λ_min={eigs_perp[0]:.4e}")

    # ===== BORDERED MATRIX ANALYSIS =====
    # Think of M_σ as [[A, b], [b^T, c]] where:
    # A = M_σ|v⊥ (the (K-1)×(K-1) restriction)
    # b = projection of M_σ v onto v⊥ (should be ~0 if aligned)
    # c = <v, M_σ v> = μ_min

    Mv = M_sigma @ v_pole
    c = np.dot(v_pole, Mv)  # Should be ≈ μ_min
    b = P @ Mv - c * v_pole  # Coupling: M_σ v restricted to v⊥
    b_norm = np.linalg.norm(b)

    print(f"\nBordered structure:")
    print(f"  c = <v, M_σ v> = {c:.6e} (should match λ_min={mu_min:.6e})")
    print(f"  b (coupling) norm = {b_norm:.6e} (should be ~0 if v is aligned)")

    # The Schur complement (simplest form, if A is invertible):
    # schur_c = c - b^T A^{-1} b
    # For PSD A with possible nullspace, use pseudoinverse:
    A = M_perp
    try:
        A_inv = np.linalg.inv(A)
        schur_c_from_inv = c - b @ A_inv @ b
        print(f"  Schur complement (via A^-1): {schur_c_from_inv:.6e}")
    except np.linalg.LinAlgError:
        print(f"  A is singular, using pseudoinverse...")
        A_pinv = np.linalg.pinv(A)
        schur_c_from_pinv = c - b @ A_pinv @ b
        print(f"  Schur complement (via A^+): {schur_c_from_pinv:.6e}")

    # ===== REFORMULATION AS GENERALIZED EIGENVALUE =====
    # RH ⟺ exactly one λ_top of the generalized problem (prime, arch) > 1.
    # At the true primes, κ = 1 always. At perturbation, κ explodes.

    if delta_perturb is not None:
        print(f"\n=== Perturbation Test: δ={delta_perturb} ===")
        # Perturb the primes (this is a REPORT-level test, not rigorous proof)
        # For now, just check the robustness of the pole-alignment
        print("(Perturbation test deferred to multi-seed verification pass)")

    return {
        'kappa_full': kappa_full,
        'kappa_perp': kappa_perp,
        'mu_min': mu_min,
        'alignment': alignment,
        'b_norm': b_norm,
        'c': c,
    }

def analyze_generalized_eigenvalue_lens(a=1.6, K=120, NG=320):
    """
    Alternative lens: the generalized eigenvalue problem (prime, arch).
    κ = #{λ > 1} in the generalized-eigenvalue sense.

    This is a different formulation of the same condition, based on §980.
    """
    print(f"\n=== Generalized Eigenvalue Lens: a={a}, K={K} ===")

    M_sigma, arch, prime, pieces = build_M_sigma(a, K, NG)

    # Generalized eigenvalue problem: prime x = λ arch x
    try:
        from scipy.linalg import eigh as eigh_gen
        eigenvals = eigh_gen(prime, arch)[0]
    except Exception as e:
        print(f"  Warning: generalized eigh failed ({e}), skipping")
        return None, None

    kappa_gen = int((eigenvals > 1.0 + 1e-10).sum())
    print(f"Generalized eigenvalues (prime, arch):")
    print(f"  λ_top = {eigenvals[-1]:.6e} (>1? {eigenvals[-1] > 1.0})")
    print(f"  λ_2nd = {eigenvals[-2]:.6e} (>1? {eigenvals[-2] > 1.0})")
    print(f"  κ = #{λ > 1} = {kappa_gen}")
    print(f"  Margin λ_top - 1 = {eigenvals[-1] - 1.0:.6e}")

    return eigenvals, kappa_gen

if __name__ == "__main__":
    # Test at the truth point (real primes)
    result_truth = test_bordering_condition(a=1.6, K=120, NG=320)
    eigenvals_gen, kappa_gen = analyze_generalized_eigenvalue_lens(a=1.6, K=120, NG=320)

    print("\n" + "="*80)
    print("FINDINGS (this run):")
    print(f"  κ(M_σ) = {result_truth['kappa_full']} (should be 1 at truth)")
    print(f"  κ(M_σ|v⊥) = {result_truth['kappa_perp']} (should be 0 — PSD)")
    print(f"  Pole-negative alignment = {result_truth['alignment']:.6f} (should be ~1)")
    if kappa_gen is not None:
        print(f"  κ (generalized eig) = {kappa_gen} (should be 1)")
    print("\nNext steps:")
    print("  1. Test faithfulness: does the condition break under prime perturbation?")
    print("  2. Derive the Schur-complement condition as an explicit inequality about arch+pole vs prime")
    print("  3. Can this be restated without needing the full M_σ spectrum?")
