import Mathlib

open Matrix

/-- **Linear-algebra core of the Weil criterion's backward direction (off-line zeros ⟹ not PSD).**

An off-line zero of an L-function contributes to the Weil quadratic form a rank-2 term
`vecMulVec a a - vecMulVec b b`, with `a = Re w`, `b = Im w` for the (complex) test-function
evaluation vector `w` at the zero (on the line `w` is real, `b = 0`, and the term is the positive
Gram piece `vecMulVec a a`). Whenever some direction `x` annihilates `a` but not `b` — `Im w` not
parallel to `Re w`, the generic off-line case — this term is NOT positive semidefinite:
`xᵀ(aaᵀ - bbᵀ)x = (a⬝x)² - (b⬝x)² = -(b⬝x)² < 0`.

Machine-checked heart of the numerical degree-2 ghost-zero detector (cycle 253): moving a zero off
the critical line injects exactly this indefinite rank-2 term into the Weil form. -/
theorem rank2_diff_not_posSemidef {n : Type*} [Fintype n] [DecidableEq n]
    (a b x : n → ℝ) (h1 : a ⬝ᵥ x = 0) (h2 : b ⬝ᵥ x ≠ 0) :
    ¬ (vecMulVec a a - vecMulVec b b).PosSemidef := by
  intro hpsd
  have key : ∀ v : n → ℝ, star x ⬝ᵥ (vecMulVec v v *ᵥ x) = (v ⬝ᵥ x) * (v ⬝ᵥ x) := by
    intro v
    simp only [star_trivial, dotProduct, mulVec, vecMulVec_apply]
    rw [Finset.sum_mul_sum]
    simp only [Finset.mul_sum]
    apply Finset.sum_congr rfl; intro i _
    apply Finset.sum_congr rfl; intro j _
    ring
  have hx := hpsd.dotProduct_mulVec_nonneg x
  rw [sub_mulVec, dotProduct_sub, key a, key b, h1] at hx
  nlinarith [mul_self_pos.mpr h2]
