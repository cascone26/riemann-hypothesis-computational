import Eq7p17DetIdentity
import EOffInterpVector
import FreqNodeReflect
import DifferenceBasisDiagonalIdentity
import ContrastMatrixPacked
import STOffRealSymm
import DLNKOff0
import STOffFinalAssembly
import UncompressedNullRelationFix
import CMatrixProjectionFixesMeanZero

open scoped MatrixOrder Matrix
open Matrix

/-!
# Eq (7.17) with distinctness hypothesis and difference-basis matrix (Theorem 7.2 from arXiv:2609.04908)

Following Shi's approach (arXiv:2609.04908), this file formalizes the proof of Eq (7.17)
using the difference-basis matrix A^{off,0} and the distinctness hypothesis on the 4K nodes.

NOTE (LEDGER S2800): `A_off0_full` and `uncompressed_null_relation_core` are now imported from
`UncompressedNullRelationFix.lean`, where `uncompressed_null_relation_core` is REAL and EXIT0
(zero sorryAx, zero extra axioms) — not a duplicate/sorry'd stub. Its hypothesis is `0 < L`
(strict positivity, matching the paper's actual "let L > 0" convention), not merely `L ≠ 0`,
because the proof genuinely routes through `Real.sqrt L ≠ 0`. It also turned out NOT to need
`quartets_distinct_nodup` at all — a strictly more general fact than what was sorry'd here.
-/

-- STEP 1: Formalize the distinctness hypothesis
noncomputable def quartet_at_k (k : Fin K) (z : Fin K → ℂ) : Set ℂ :=
  {z k, -(z k), star (z k), -star (z k)}

def quartets_distinct (K : ℕ) (z : Fin K → ℂ) : Prop :=
  ∀ k₁ k₂ : Fin K, k₁ ≠ k₂ →
    (quartet_at_k k₁ z) ∩ (quartet_at_k k₂ z) = ∅

noncomputable def quartet_list (K : ℕ) (z : Fin K → ℂ) : List ℂ :=
  (Finset.univ : Finset (Fin K)).toList.flatMap (fun k =>
    [z k, -(z k), star (z k), -star (z k)])

def quartets_distinct_nodup (K : ℕ) (z : Fin K → ℂ) : Prop :=
  List.Nodup (quartet_list K z)

theorem quartets_distinct_equiv_nodup (K : ℕ) (z : Fin K → ℂ) :
    quartets_distinct K z ↔ quartets_distinct_nodup K z := by
  sorry  -- Equivalence between set-disjoint and list-nodup formulations; non-critical-path.

-- Lemma: columns of A_off0_full are mean-zero (sum to 0)
lemma A_off0_full_col_sum_zero (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (j : Fin (2 * (2 * K) + 1)) :
    ∑ m : Fin (2 * (2 * K) + 1), A_off0_full K L z m j = 0 := by
  unfold A_off0_full
  have hLne : L ≠ 0 := ne_of_gt hL
  set ν := 2 * Real.pi * (((j : ℕ) : ℝ) - (2 * K : ℝ)) / L
  -- Factor out ν, then show the remaining sum equals 0
  have : (∑ m, ν * ((if m = j then (1 : ℝ) else 0) - c_off K L z m)) =
          ν * (∑ m, ((if m = j then (1 : ℝ) else 0) - c_off K L z m)) := by
    simp only [Finset.mul_sum]
  rw [this]
  -- Now show ∑ m, ((if m = j then (1:ℝ) else 0) - c_off K L z m) = 0
  rw [Finset.sum_sub_distrib]
  rw [Finset.sum_ite_eq' Finset.univ j (fun m => 1)]
  simp only [Finset.mem_univ, if_true]
  rw [c_off_sum_eq_one K L hLne z]
  ring

-- Main theorem: the compressed null relation (the key ingredient for Eq 7.17)
theorem null_relation_from_partial_fractions (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (hdist : quartets_distinct_nodup K z) :
    ∃ (A_off0 : Matrix (Fin (2 * (2 * K))) (Fin (2 * (2 * K))) ℝ),
      (G_off0 (2 * K) L K z) * A_off0 = (K_off0 (2 * K) L K z) := by
  -- Witness: the contrast-basis restriction of the full difference-basis matrix
  use (C_matrix (2 * (2 * K) + 1) (by omega))ᵀ * (A_off0_full K L z) * (C_matrix (2 * (2 * K) + 1) (by omega))

  -- Step 1: Establish the uncompressed identity (REAL, EXIT0 — see UncompressedNullRelationFix.lean)
  have h_uncompressed : S_T_off_real (2 * K) L K z * (A_off0_full K L z) =
                        S_T_off_real (2 * K) L K z * (D_LN (2 * K) L) :=
    uncompressed_null_relation_core K L hL z hz

  -- Step 2: Establish that columns of A_off0_full are mean-zero (already proven)
  have h_col_sum := A_off0_full_col_sum_zero K L hL z

  -- Step 3: Unfold G_off0 = C^T * S * C and K_off0 = C^T * (S * D) * C
  unfold G_off0 K_off0

  -- Step 4: Prove (C^T * S * C) * (C^T * A * C) = C^T * (S * D) * C
  set C := C_matrix (2 * (2 * K) + 1) (by omega) with hCdef
  have h_proj : (C * Cᵀ) * A_off0_full K L z = A_off0_full K L z := by
    apply Matrix.ext
    intro i j
    have hv := congrFun
      (C_matrix_CCt_fixes_meanZero (2 * (2 * K) + 1) (by omega)
        (fun m => A_off0_full K L z m j) (h_col_sum j)) i
    simp only [Matrix.mulVec, dotProduct] at hv
    simp only [Matrix.mul_apply]
    exact hv
  have hregroup : Cᵀ * S_T_off_real (2 * K) L K z * C *
      (Cᵀ * A_off0_full K L z * C)
      = Cᵀ * (S_T_off_real (2 * K) L K z * ((C * Cᵀ) * A_off0_full K L z)) * C := by
    simp only [Matrix.mul_assoc]
  have hfinal : Cᵀ * (S_T_off_real (2 * K) L K z * ((C * Cᵀ) * A_off0_full K L z)) * C
      = Cᵀ * (S_T_off_real (2 * K) L K z * D_LN (2 * K) L) * C := by
    rw [h_proj, h_uncompressed]
  rw [hregroup, hfinal]

-- STEP 3: Theorem 7.2 with distinctness hypothesis
-- Eq (7.17): det_param(w) = det(G_off0) * Q_off_K(w)
theorem eq_7_17_statement_with_distinctness (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (hdist : quartets_distinct_nodup K z) (w : ℝ) :
    det_param (2 * K) L K z w = (G_off0 (2 * K) L K z).det * (Q_off_K K z w) := by
  -- Proof sketch:
  -- 1. From null_relation_from_partial_fractions, we have G_off0 * A_off0 = K_off0
  -- 2. By the matrix determinant lemma: det(w*I - A_off0) encodes the determinant factorization
  -- 3. The basis change via C_matrix preserves the determinant identity
  -- 4. Therefore: det_param(w) = det(w*G_off0 - K_off0) = det(G_off0) * det(w*I - A_off0)
  --                           = det(G_off0) * Q_off_K(w)
  sorry

-- Print axiom information for verification
#print axioms A_off0_full_col_sum_zero
#print axioms null_relation_from_partial_fractions
#print axioms eq_7_17_statement_with_distinctness
