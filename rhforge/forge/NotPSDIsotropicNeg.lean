import Mathlib

open Matrix

/-- **Composition primitive for the backward Weil criterion.**

If a direction `x` is isotropic for `P` (`xᵀP x = 0` — e.g. `x` in the null space of the PSD Gram
sum over the on-line zeros) and strictly negative for `D` (`xᵀD x < 0` — e.g. `D = aaᵀ - bbᵀ` the
indefinite rank-2 term contributed by one off-line zero, via `rank2_diff_not_posSemidef`), then the
combined form `P + D` is NOT positive semidefinite: `xᵀ(P+D)x = 0 + (neg) < 0`.

Composed with `WeilFormPSD` (the on-line Gram sum is PSD) and `rank2_diff_not_posSemidef`, this is
the backward direction of the Weil criterion in Lean: a single off-line zero destroys Weil
positivity. -/
theorem not_posSemidef_of_isotropic_neg {n : Type*} [Fintype n] [DecidableEq n]
    (P D : Matrix n n ℝ) (x : n → ℝ)
    (hP : x ⬝ᵥ (P *ᵥ x) = 0) (hD : x ⬝ᵥ (D *ᵥ x) < 0) :
    ¬ (P + D).PosSemidef := by
  intro hpsd
  have hx := hpsd.dotProduct_mulVec_nonneg x
  rw [add_mulVec, dotProduct_add] at hx
  simp only [star_trivial] at hx
  rw [hP] at hx
  linarith [hD]
