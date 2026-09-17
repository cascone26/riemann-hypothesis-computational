import Mathlib.Tactic

/-!
# §1810 — Weighted form PSD when all real roots lie in [−2,2] (the ⟸ of the interval iff)

Multiset companion of §1791: for a real root multiset `R` with every `r² ≤ 4` (roots in `[−2,2]`) and real
coefficients, the `(4−r²)`-weighted Hankel form is nonnegative — each term is a nonneg weight times a square.
With §1809 (the ⟹) this gives the interval iff.  EXIT0 / std-axioms.
-/

open Finset

variable {d : ℕ}

/-- **All roots in [−2,2] ⟹ weighted form PSD.** -/
theorem interval_weighted_nonneg (R : Multiset ℝ) (c : Fin (d + 1) → ℝ)
    (hR : ∀ r ∈ R, r ^ 2 ≤ 4) :
    0 ≤ (R.map (fun r => (4 - r ^ 2) * (∑ i : Fin (d + 1), c i * r ^ (i : ℕ)) ^ 2)).sum := by
  apply Multiset.sum_nonneg
  intro y hy
  simp only [Multiset.mem_map] at hy
  obtain ⟨r, hr, rfl⟩ := hy
  exact mul_nonneg (by linarith [hR r hr]) (sq_nonneg _)
