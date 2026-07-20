import Mathlib

/-!
# Route-6 tail-safety lemma — elementary kernel (formalized, zero `sorry`)

The RH route-6 tail lemma (LEDGER §21) reduces, via the exact closed form
`F_n(σ,γ) = 4 - 4 cos(nθ) cosh(nu)`, to the scalar inequality `cos(nθ) cosh(nu) ≤ 1`
under `θ ≥ u ≥ 0`, `n(θ+u) ≤ 2π`, `0 < nθ < 2π`.  Writing `x = nθ`, `w = nu`:
the two nontrivial regions are `x ∈ (0, π/2]` (uses `w ≤ x`) and `x ∈ [3π/2, 2π)`
(uses `x + w ≤ 2π` via `τ = 2π - x ≥ w`); the middle `x ∈ [π/2, 3π/2]` has `cos x ≤ 0`.

The two analytic inputs are TEXTBOOK facts, taken here as explicit hypotheses so the
**novel combinatorial kernel** (the case-split producing the `2πγ` threshold) is
machine-checked with NO `sorry`:
  (B1)  `cos t ≤ exp(-(t²)/2)`  for `t ∈ [0, π/2]`   (⇔ `log cos t ≤ -t²/2`, since `tan t ≥ t`)
  (B2)  `cosh z ≤ exp(z²/2)`    for all `z`          (term-by-term: `(2n)! ≥ 2ⁿ·n!`)
-/

open Real

/-- **Tail-lemma kernel** (fully proven from the two textbook analytic inputs B1, B2).
Given `x = nθ ∈ (0, 2π)`, `w = nu ≥ 0`, the structural hypotheses `w ≤ x` (from `θ ≥ u`)
and `x + w ≤ 2π` (from `n(θ+u) ≤ 2π`), the product `cos x · cosh w ≤ 1`; hence
`F = 4 - 4·cos x·cosh w ≥ 0`. -/
theorem tail_kernel {x w : ℝ} (hx0 : 0 < x) (hx2 : x < 2 * π) (hw0 : 0 ≤ w)
    (hlow : w ≤ x) (hhigh : x + w ≤ 2 * π)
    (B1 : ∀ t : ℝ, 0 ≤ t → t ≤ π / 2 → cos t ≤ exp (-(t ^ 2) / 2))
    (B2 : ∀ z : ℝ, cosh z ≤ exp (z ^ 2 / 2)) :
    cos x * cosh w ≤ 1 := by
  rcases le_or_gt (cos x) 0 with hcos | hcos
  · -- cos x ≤ 0 : product ≤ 0 ≤ 1
    have : cos x * cosh w ≤ 0 := mul_nonpos_of_nonpos_of_nonneg hcos (by positivity)
    linarith
  · rcases le_or_gt x (π / 2) with hxle | hxgt
    · -- low branch:  cos x ≤ e^{-x²/2},  cosh w ≤ e^{w²/2},  w ≤ x  ⇒  product ≤ e^{-(x²-w²)/2} ≤ 1
      have hb1 : cos x ≤ exp (-(x ^ 2) / 2) := B1 x (le_of_lt hx0) hxle
      have hb2 : cosh w ≤ exp (w ^ 2 / 2) := B2 w
      have hprod : cos x * cosh w ≤ exp (-(x ^ 2) / 2) * exp (w ^ 2 / 2) :=
        mul_le_mul hb1 hb2 (by positivity) (by positivity)
      rw [← Real.exp_add] at hprod
      have hexp : -(x ^ 2) / 2 + w ^ 2 / 2 ≤ 0 := by nlinarith [hlow, hw0, hx0]
      calc cos x * cosh w ≤ exp (-(x ^ 2) / 2 + w ^ 2 / 2) := hprod
        _ ≤ exp 0 := Real.exp_le_exp.mpr hexp
        _ = 1 := Real.exp_zero
    · -- high branch:  cos x > 0 with x ∈ (π/2, 2π) forces x > 3π/2.  τ = 2π - x ∈ (0, π/2], cos x = cos τ.
      have hx3 : 3 * π / 2 < x := by
        by_contra h
        push_neg at h
        have : cos x ≤ 0 :=
          Real.cos_nonpos_of_pi_div_two_le_of_le (le_of_lt hxgt) (by linarith [Real.pi_pos])
        linarith
      have hτpos : 0 < 2 * π - x := by linarith
      have hτle : 2 * π - x ≤ π / 2 := by linarith
      have hcosτ : cos x = cos (2 * π - x) := by rw [Real.cos_two_pi_sub]
      rw [hcosτ]
      have hb1 : cos (2 * π - x) ≤ exp (-((2 * π - x) ^ 2) / 2) := B1 _ (le_of_lt hτpos) hτle
      have hb2 : cosh w ≤ exp (w ^ 2 / 2) := B2 w
      have hprod : cos (2 * π - x) * cosh w ≤ exp (-((2 * π - x) ^ 2) / 2) * exp (w ^ 2 / 2) :=
        mul_le_mul hb1 hb2 (by positivity) (by positivity)
      rw [← Real.exp_add] at hprod
      have hexp : -((2 * π - x) ^ 2) / 2 + w ^ 2 / 2 ≤ 0 := by nlinarith [hhigh, hw0, hτpos]
      calc cos (2 * π - x) * cosh w ≤ exp (-((2 * π - x) ^ 2) / 2 + w ^ 2 / 2) := hprod
        _ ≤ exp 0 := Real.exp_le_exp.mpr hexp
        _ = 1 := Real.exp_zero

/-- Restated as the nonnegativity of `F = 4 - 4·cos x·cosh w`. -/
theorem tail_F_nonneg {x w : ℝ} (hx0 : 0 < x) (hx2 : x < 2 * π) (hw0 : 0 ≤ w)
    (hlow : w ≤ x) (hhigh : x + w ≤ 2 * π)
    (B1 : ∀ t : ℝ, 0 ≤ t → t ≤ π / 2 → cos t ≤ exp (-(t ^ 2) / 2))
    (B2 : ∀ z : ℝ, cosh z ≤ exp (z ^ 2 / 2)) :
    0 ≤ 4 - 4 * (cos x * cosh w) := by
  have := tail_kernel hx0 hx2 hw0 hlow hhigh B1 B2
  linarith
