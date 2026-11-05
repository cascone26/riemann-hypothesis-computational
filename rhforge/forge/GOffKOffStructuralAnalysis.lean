import GOff0Symm
import DLNKOff0
import Mathlib

open scoped MatrixOrder Matrix

/-!
# PART 2: G_off0 and K_off0 Structural Relationship Analysis

INVESTIGATION: Can the characteristic polynomial identity
  `det(w·1 - G_off0⁻¹·K_off0) = Q_off_K(w)`
be proven via a different approach that leverages the ALREADY-PROVEN matrix-null-vector machinery?

KEY OBSERVATION (STRUCTURAL DUALITY):
=====================================

Both G_off0 and K_off0 are congruences using the SAME C_matrix:

  G_off0 = C^T · S_T_off_real · C                    (4K × 4K matrix)
  K_off0 = C^T · (S_T_off_real · D_LN) · C           (4K × 4K matrix)

where:
  - C_matrix is the isometric contrast embedding (injective, full column rank)
  - S_T_off_real is the indefinite symmetric matrix (rank 4K, with one null vector e_off)
  - D_LN is the frequency diagonal matrix with ν_{n,L} on diagonal

CRITICAL STRUCTURAL INSIGHT:
============================

The product structure K_off0 = C^T · (S_T_off_real · D_LN) · C shows that K_off0 is
NOT independent of G_off0 — it is derived from the SAME base matrix S_T_off_real,
but weighted by the DIAGONAL frequency matrix D_LN.

This suggests that the eigenvalues of G_off0^{-1}·K_off0 might be directly related
to the FREQUENCY NODES (the diagonal entries of D_LN) rather than purely to z_k values.

HYPOTHESIS (TO BE INVESTIGATED):
================================

If we write S_T_off_real in the basis where D_LN acts naturally (i.e., in the
frequency/Fourier basis), the relationship between G_off0 and K_off0 might
decompose into a DIRECT SCALING by the frequency nodes — giving:

  (G_off0)^{-1}·(K_off0) = (C^T·S·C)^{-1}·(C^T·S·D_LN·C)
                         = (C^T·S·C)^{-1}·(C^T·S·C)·(something involving D_LN)

provided certain compatibility conditions hold (which may follow from the
already-proven null-vector relationships).

CURRENT KNOWLEDGE (from LEDGER §2782):
======================================

1. G_off0 invertibility reduces to ONE condition: `C^T *_v e_off = 0`
2. This condition is "almost certainly FALSE" based on numerical tests (dimension argument)
3. Therefore, G_off0 invertibility cannot be deduced from the orthogonality claim
4. BUT: this does NOT rule out G_off0 being invertible via OTHER mechanisms

POSSIBLE RESOLUTION PATHS:
==========================

**Path A: G_off0 may be invertible despite the orthogonality claim being false**
  - The matrix S_T_off_real might have a different structural property that ensures
    C^T·S_T_off_real·C is nonsingular even when e_off is NOT orthogonal to Im(C)
  - The fact that S_T_off_real has RANK EXACTLY 4K (one null vector in a (4K+1)×(4K+1) matrix)
    and C has FULL COLUMN RANK might combine in a way that makes the congruence invertible
  - This would be a "generic" invertibility property of congruences of rank-deficient matrices

**Path B: The characteristic polynomial identity follows from the congruence structure directly**
  - Rather than proving G_off0 is invertible first, we might directly relate
    det(w·1 - G_off0^{-1}·K_off0) to Q_off_K by analyzing the congruence relationships
  - The fact that both come from the same S_T_off_real suggests there's a UNIFIED
    algebraic identity that doesn't require explicit inversion
  - For example: maybe there's a Schur complement or matrix pencil identity that relates
    the two determinants without requiring G_off0^{-1} to be explicitly defined

**Path C: Use explicit hypothesis for G_off0 invertibility**
  - Rather than prove it from first principles, state eq_7_17 with the assumption
    `(hG : IsUnit (G_off0 ...))`
  - This matches how Shi's original paper treats it (implicitly assumed)
  - This unblocks further progress while leaving the invertibility question for future work
  - The problem becomes well-scoped: assuming invertibility, prove the characteristic polynomial

MATHEMATICAL BRIDGE (HYPOTHETICAL):
===================================

If S_T_off_real has a spectral decomposition S = Q·Λ·Q^T where Λ is diagonal
with the generalized eigenvalues as diagonal entries, and if C*Q relates to the
frequency basis in some way, then the congruence relationships might simplify.

This is speculative without explicit computation, but it suggests the path
forward is through:
1. Spectral analysis of S_T_off_real (via already-proven null-vector properties)
2. Understanding how C_matrix aligns with or relates to that spectral structure
3. Using the diagonal frequency structure of D_LN to connect eigenvalues to ν_n

NEXT IMMEDIATE STEPS (TO UNBLOCK eq_7_17):
===========================================

**Option 1 (Recommended for immediate progress):**
  Add explicit hypothesis `(hG : IsUnit (G_off0 N L K z).det)` to eq_7_17_statement
  This is mathematically sound and unblocks the rest of the identity proof.

**Option 2 (Research-grade, higher effort):**
  Prove G_off0 invertibility via spectral/geometric argument using S_T_off_real properties
  (potentially using the already-proven null-vector machinery from §2368 onwards)

**Option 3 (Alternative characterization):**
  Prove the characteristic polynomial identity WITHOUT requiring explicit G_off0 inversion,
  using a direct matrix pencil or Schur complement approach on the base matrices

FORMALIZATION STATUS:
====================

This file is analytical/exploratory — it documents the investigation without
claiming to have resolved the question. The structural insights are real and
useful for future work, even if the hypothesis about spectral decomposition
remains to be verified.

-/

/-- Structural fact: G_off0 and K_off0 share the same C_matrix congruence basis. -/
lemma G_off0_K_off0_congruence_structure (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) :
    G_off0 N L K z = (C_matrix (2 * N + 1) (by omega))ᵀ * (S_T_off_real N L K z) * (C_matrix (2 * N + 1) (by omega)) ∧
    K_off0 N L K z = (C_matrix (2 * N + 1) (by omega))ᵀ * ((S_T_off_real N L K z) * (D_LN N L)) * (C_matrix (2 * N + 1) (by omega)) := by
  exact ⟨rfl, rfl⟩
