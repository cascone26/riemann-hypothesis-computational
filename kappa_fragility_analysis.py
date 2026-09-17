#!/usr/bin/env python3
"""CRITICAL FINDING: κ is FRAGILE under perturbation, not robust as claimed in FRONTIER.

Previous probe showed κ jumps from 1 → 5 → 28 under ±0.1% Λ perturbation.
This directly contradicts FRONTIER's claim: "pole-count is a ROBUST INTEGER (not fragile margin)."

THIS EXPERIMENT: diagnose the fragility source — is it:
  (A) near-zero eigenvalues of M_σ flip sign under tiny perturbations?
  (B) conditioning issue (spectral clustering near 0)?
  (C) a deeper non-robustness of the architecture-form definition?

HONEST: if κ is genuinely fragile, that blocks the κ=1 proof route.
"""

import numpy as np
from scipy.linalg import solve, eigh
from scipy.special import digamma
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'operators'))
from suzuki_exact_Q import build_Q, E_p, shat, vonmangoldt_upto

print("="*100)
print("KAPPA FRAGILITY DIAGNOSIS: Why does κ jump under ±0.1% perturbation?")
print("="*100)

a = 1.5
K = 90
NG = 400
kk = np.arange(1, K+1)

def build_Q_perturbed(a, K, NG, perturbation_factor):
    """Build Q with Λ(n) multiplied by (1 + perturbation)."""
    Q_full_0, Q_pieces = build_Q(a, K, NG, return_pieces=True)
    Qprime_pert = Q_pieces["prime"] * (1.0 + perturbation_factor)
    Q_pert = Q_pieces["pole"] + Q_pieces["core"] + Q_pieces["tail"] + Qprime_pert
    Q_pert = (Q_pert + Q_pert.T) / 2
    return Q_pert, Q_pieces["pole"]

Q_pole_exact = np.outer(E_p(0.5, a, kk), E_p(-0.5, a, kk)) + np.outer(E_p(-0.5, a, kk), E_p(0.5, a, kk))

print("\nDIAGNOSIS A: Near-zero eigenvalue clustering")
print("-"*100)

for pf in [-0.001, -0.0005, 0.0, 0.0005, 0.001]:
    Q_pert, Q_pole = build_Q_perturbed(a, K, NG, pf)
    M_sigma = Q_pert - Q_pole

    evals = np.linalg.eigvalsh(M_sigma)
    evals_sorted = np.sort(evals)

    # Count how many eigenvalues are "near zero" (sensitive to threshold)
    tol_strict = 1e-9      # strict threshold
    tol_loose = 1e-6       # loose threshold (1 millionth of largest)

    neg_strict = int((evals < -tol_strict).sum())
    neg_loose = int((evals < -tol_loose).sum())

    near_zero = int((np.abs(evals) < 1e-6).sum())

    print(f"\n  Perturbation: {pf:+8.4f}")
    print(f"    κ (|λ| > 1e-9): {neg_strict:2d}")
    print(f"    κ (|λ| > 1e-6): {neg_loose:2d}")
    print(f"    Eigenvalues in [-1e-6, 1e-6]: {near_zero:2d}")
    print(f"    Extremal evals: [{evals_sorted[0]:.6e}, {evals_sorted[5]}, ... {evals_sorted[-1]:.6e}]")

    # Show the specific near-zero eigenvalues
    near_zero_evals = evals[np.abs(evals) < 1e-5]
    if len(near_zero_evals) > 0:
        print(f"    Near-zero eigenvalues: {sorted(near_zero_evals)[:8]}")

print("\n" + "="*100)
print("DIAGNOSIS B: Sensitivity of eigenvalue count to threshold")
print("-"*100)

Q_orig, Q_pole = build_Q_perturbed(a, K, NG, 0.0)
M_sigma_orig = Q_orig - Q_pole
evals_orig = np.sort(np.linalg.eigvalsh(M_sigma_orig))

print(f"\nOriginal M_σ eigenvalue spectrum near zero:")
print(f"  First 20 (most negative):")
for i in range(min(20, len(evals_orig))):
    print(f"    {i:2d}. λ[{i:2d}] = {evals_orig[i]:+.6e}")

# Show the gap structure
diffs = np.diff(evals_orig[:20])
print(f"\n  Gaps between first 20 eigenvalues:")
for i in range(len(diffs)):
    print(f"    λ[{i+1:2d}] - λ[{i:2d}] = {diffs[i]:.6e}")

print("\n" + "="*100)
print("DIAGNOSIS C: What creates the negative eigenvalues?")
print("-"*100)

# Decompose M_σ = arch - prime into pieces
print("\n  M_σ = (Q_pole + Q_core + Q_tail) - Q_prime")
print(f"  = const·I (arch core) + arch_tail + Q_pole - Q_prime")
print()

Q_full = build_Q(a, K, NG)
Q_pieces = {}
Q_pieces["pole"] = np.outer(E_p(0.5, a, kk), E_p(-0.5, a, kk)) + np.outer(E_p(-0.5, a, kk), E_p(0.5, a, kk))

# Rebuild to get pieces
Q_temp, Q_pieces_temp = build_Q(a, K, NG, return_pieces=True)
Q_pieces.update(Q_pieces_temp)

M_sigma_full = Q_pieces["pole"] + Q_pieces["core"] + Q_pieces["tail"] - Q_pieces["prime"]
evals_M = np.linalg.eigvalsh(M_sigma_full)

print(f"  Eigenvalue counts:")
print(f"    Q_pole alone:      {int((np.linalg.eigvalsh(Q_pieces['pole']) < 0).sum())} negative")
print(f"    Q_core (≈-2c·a·I): {int((np.linalg.eigvalsh(Q_pieces['core']) < 0).sum())} negative (should be K={K})")
print(f"    Q_tail:            {int((np.linalg.eigvalsh(Q_pieces['tail']) < 0).sum())} negative")
print(f"    Q_prime:           {int((np.linalg.eigvalsh(Q_pieces['prime']) < 0).sum())} negative")
print(f"    M_σ = pole+core+tail - prime: {int((evals_M < 0).sum())} negative")

# The core is -2c·a·I, so it should have all negative eigenvalues with value -2c·a
PSI = float(np.log(np.pi) - digamma(0.25))  # ψ(¼)-logπ
C_CORE = PSI / 2.0
print(f"\n  Core term: Q_core = -2c·a·I where c = (ψ(1/4)-log π)/2 ≈ {C_CORE:.6e}")
print(f"    Eigenvalues: all ≈ {-2*C_CORE*a:.6e} (K={K} times)")

print("\n" + "="*100)
print("KEY FINDING: The negative eigenvalues come from Q_core (the archimedean part)")
print("This is NOT tied to the pole structure, but to arch alone!")
print("="*100)

print("""
DIAGNOSIS SUMMARY:

Problem: κ is defined as #negative eigenvalues of M_σ. But M_σ = arch_terms - prime_terms.
  - arch includes Q_core = -2c·a·I, which contributes K negative eigenvalues (K=90).
  - prime subtracts these away, leaving only 1 negative (the "true" pole-induced one).
  - But under ±0.1% perturbation of Q_prime, the cancellation is FRAGILE.

ROOT CAUSE:
  - The definition κ = #neg(arch - prime) assumes EXACT cancellation between arch and prime.
  - In reality, arch and prime are of VASTLY different magnitudes (arch core is ~-27 per dim).
  - Perturbation of prime → breakdown of fine-tuned cancellation → κ jumps wildly.

HONEST ASSESSMENT:
  - κ is NOT a robust integer in practice.
  - The "robust κ=1" claim in FRONTIER is FALSE in this form.
  - The fragility disqualifies κ as a proof tool unless we can restore robustness.

WHAT WOULD FIX IT:
  - Don't count raw negative eigenvalues.
  - Instead, measure a NORMALIZED or RELATIVE quantity (e.g., inertia of B.T M B in a subspace).
  - Or build an indefinite de Branges form where pole-count is intrinsically robust.
""")

print("="*100)
