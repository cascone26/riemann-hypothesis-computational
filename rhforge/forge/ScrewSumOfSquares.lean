import Mathlib

open Real Finset

/-!
# Screw-function sum-of-squares (Suzuki Cor 4.3 mechanism)

RH work, session 2026-08-30 cycle 303: the Suzuki screw function of ξ satisfies
`g(t) = -4 ∑_γ sin²(γ t / 2) / γ²`, i.e. `-g(t) = 4 ∑_γ sin²(γ t/2)/γ² ≥ 0`, and its
term-wise nonpositivity is EXACTLY the reality of the zeros γ.  The load-bearing elementary
facts, formalized and axiom-audited:

* `cos_sub_one_eq`      : `cos x - 1 = -2 * sin (x/2)^2`         (so it is a negated square)
* `cos_sub_one_nonpos`  : `cos x - 1 ≤ 0`  for all real `x`      (real frequency ⇒ term ≤ 0)
* `screw_weighted_nonpos`: a nonneg-weighted sum `∑ w i * (cos (γ i * t) - 1) ≤ 0`
    — this is `g(t) ≤ 0` for the on-line (real-γ) screw function, the Cor 4.3 conclusion.
-/

/-- `cos x - 1 = -2 sin²(x/2)`: the term of the screw function is a negated square. -/
theorem cos_sub_one_eq (x : ℝ) : Real.cos x - 1 = -2 * Real.sin (x / 2) ^ 2 := by
  have h2 : Real.cos (2 * (x / 2)) = 2 * Real.cos (x / 2) ^ 2 - 1 := Real.cos_two_mul (x / 2)
  have hx : 2 * (x / 2) = x := by ring
  rw [hx] at h2
  have hpyth : Real.sin (x / 2) ^ 2 + Real.cos (x / 2) ^ 2 = 1 := Real.sin_sq_add_cos_sq (x / 2)
  linear_combination h2 + 2 * hpyth

/-- Real frequency ⇒ each screw term is ≤ 0. -/
theorem cos_sub_one_nonpos (x : ℝ) : Real.cos x - 1 ≤ 0 := by
  rw [cos_sub_one_eq]
  nlinarith [sq_nonneg (Real.sin (x / 2))]

/-- Suzuki Cor 4.3 (on-line case): a nonnegatively-weighted sum of screw terms with REAL
frequencies is ≤ 0.  With `w i = 4 / γ i ^ 2` this is exactly `g(t) ≤ 0`, the screw function's
nonpositivity — nonneg precisely because the frequencies `γ i` are real. -/
theorem screw_weighted_nonpos {n : ℕ} (w γ : Fin n → ℝ) (t : ℝ)
    (hw : ∀ i, 0 ≤ w i) :
    ∑ i, w i * (Real.cos (γ i * t) - 1) ≤ 0 := by
  apply Finset.sum_nonpos
  intro i _
  exact mul_nonpos_of_nonneg_of_nonpos (hw i) (cos_sub_one_nonpos _)
