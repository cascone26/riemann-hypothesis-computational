import NormPartialSumCharLe
import Mathlib.Algebra.BigOperators.Module

open Finset

/-!
# §1515 — Abel bound for the twisted tail block: summation-by-parts + bounded partial sums

Applying Mathlib's `sum_Ioc_by_parts` to `f i = i^{-s}`, `g i = χ(↑i)` (partial sums `G n = ∑_{range n} χ(↑i)`
bounded by `M` via §1514), the twisted block sum obeys

  `‖∑_{i∈Ioc m n} i^{-s}·χ(↑i)‖ ≤ M·(‖n^{-s}‖ + ‖(m+1)^{-s}‖ + ∑_{i∈Ioc m (n-1)} ‖(i+1)^{-s}-i^{-s}‖)`.

The remaining `∑‖(i+1)^{-s}-i^{-s}‖` is the total variation of `i^{-s}`, bounded next by `|s|·∑ i^{-σ-1}`.  This is
the core Abel step toward the truncated `L(χ)` representation and its explicit growth bound.  EXIT0 / std-axioms.
-/

theorem norm_sum_Ioc_twist_by_parts_le {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M}
    (hχ : χ ≠ 1) {s : ℂ} {m n : ℕ} (hmn : m < n) :
    ‖∑ i ∈ Finset.Ioc m n, (i : ℂ) ^ (-s) • χ (i : ZMod M)‖
      ≤ ‖(n : ℂ) ^ (-s)‖ * M + ‖((m + 1 : ℕ) : ℂ) ^ (-s)‖ * M
        + ∑ i ∈ Finset.Ioc m (n - 1), ‖((i + 1 : ℕ) : ℂ) ^ (-s) - (i : ℂ) ^ (-s)‖ * M := by
  set f : ℕ → ℂ := fun i => (i : ℂ) ^ (-s) with hf
  set g : ℕ → ℂ := fun i => χ (i : ZMod M) with hg
  have hG : ∀ k, ‖∑ i ∈ range k, g i‖ ≤ (M : ℝ) := fun k => norm_partial_sum_char_le hχ k
  rw [Finset.sum_Ioc_by_parts f g hmn]
  refine le_trans (norm_sub_le _ _) ?_
  refine add_le_add (le_trans (norm_sub_le _ _) ?_) ?_
  · refine add_le_add ?_ ?_
    · rw [norm_smul]; exact mul_le_mul_of_nonneg_left (hG _) (norm_nonneg _)
    · rw [norm_smul]; exact mul_le_mul_of_nonneg_left (hG _) (norm_nonneg _)
  · refine le_trans (norm_sum_le _ _) ?_
    refine Finset.sum_le_sum (fun i _ => ?_)
    rw [norm_smul]
    exact mul_le_mul_of_nonneg_left (hG _) (norm_nonneg _)

#print axioms norm_sum_Ioc_twist_by_parts_le
