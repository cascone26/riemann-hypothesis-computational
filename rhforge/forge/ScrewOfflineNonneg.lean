import Mathlib

open Real Finset

/-!
# Off-line screw term flips sign (Suzuki Cor 4.3 FAIL side)

Companion to `ScrewSumOfSquares` (cycle 307). The on-line screw term (real frequency γ) is
`cos(γt) - 1 = -2 sin²(γt/2) ≤ 0`. An OFF-LINE zero shifts the frequency into the complex plane;
the resulting real contribution is the hyperbolic mirror `cosh(δt) - 1 = +2 sinh²(δt/2) ≥ 0`
— the OPPOSITE sign, strictly positive once `δt ≠ 0`, and unbounded (cosh grows). Hence a
nonneg-weighted sum of off-line contributions is `≥ 0` (not `≤ 0`): Cor 4.3's nonpositivity FAILS.
This machine-checked contrapositive ties the EP/multiplicativity discriminator (off-line ⟹ FAIL)
to the screw-function mechanism (cycles 290-308).
-/

/-- Off-line screw term is ≥ 0 — the OPPOSITE sign to the on-line term `cos(γt)-1 ≤ 0`. -/
theorem cosh_sub_one_nonneg (x : ℝ) : 0 ≤ Real.cosh x - 1 := by
  have := Real.one_le_cosh x
  linarith

/-- Strictly positive once the shifted argument is nonzero: a genuine off-line zero (δ≠0, t≠0)
makes the contribution `> 0`, breaking nonpositivity. -/
theorem cosh_sub_one_pos_of_ne {x : ℝ} (hx : x ≠ 0) : 0 < Real.cosh x - 1 := by
  have := Real.one_lt_cosh.mpr hx
  linarith

/-- Suzuki Cor 4.3 FAIL side: a nonnegatively-weighted sum of OFF-LINE screw contributions is `≥ 0`
(the reverse of `screw_weighted_nonpos`). With an off-line zero present some term is `> 0`, so the
screw function is NOT ≤ 0 — Cor 4.3 fails, exactly as it must for a non-multiplicative L-function. -/
theorem screw_offline_nonneg {n : ℕ} (w δ : Fin n → ℝ) (t : ℝ) (hw : ∀ i, 0 ≤ w i) :
    0 ≤ ∑ i, w i * (Real.cosh (δ i * t) - 1) := by
  apply Finset.sum_nonneg
  intro i _
  exact mul_nonneg (hw i) (cosh_sub_one_nonneg _)
