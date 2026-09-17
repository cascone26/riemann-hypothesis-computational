import Mathlib

/-!
# Nonneg-weighted combination of PSD is PSD (RH session 2026-08-29)

The algebraic core of the Euler-product discriminator's jump-form positivity (cycle 199, the ζ /
trivial-character case): the prime side is `∑ₙ wₙ Dₙ` with each single-jump form `Dₙ ⪰ 0`; when the
weights `wₙ = Λ(n)/√n ≥ 0` (which for ζ is the Euler product), the whole form is `⪰ 0`.
-/

open Matrix

/-- A nonneg-weighted sum of positive-semidefinite matrices is positive semidefinite. -/
theorem nonneg_combo_posSemidef {ι n : Type*} [Fintype ι] [Fintype n] [DecidableEq n]
    (w : ι → ℝ) (D : ι → Matrix n n ℝ)
    (hw : ∀ i, 0 ≤ w i) (hD : ∀ i, (D i).PosSemidef) :
    (∑ i, w i • D i).PosSemidef := by
  apply Matrix.posSemidef_sum Finset.univ
  intro i _
  exact (hD i).smul (hw i)
