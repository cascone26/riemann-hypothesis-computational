import RiemannLean.TailAnalytic

/-!
# Route-6 tail-safety lemma — elementary kernel (fully self-contained, zero `sorry`, zero hypotheses)

The RH route-6 tail lemma (LEDGER §21) reduces, via the exact closed form
`F_n(σ,γ) = 4 - 4 cos(nθ) cosh(nu)`, to `cos(nθ) cosh(nu) ≤ 1` under `θ ≥ u ≥ 0`, `n(θ+u) ≤ 2π`,
`0 < nθ < 2π`.  Writing `x = nθ`, `w = nu`, the three regions are `cos x ≤ 0` (trivial),
`x ∈ (0,π/2]` (uses `w ≤ x`), `x ∈ [3π/2,2π)` (uses `x+w ≤ 2π` via `τ = 2π−x ≥ w`).

The two analytic inputs are now PROVEN in `TailAnalytic.lean`:
`cos_le_exp_neg_half_sq` and `cosh_le_exp_half_sq`.  Hence this whole file is self-contained.
-/

open Real

/-- **Tail-lemma kernel** (fully self-contained).  For `x = nθ ∈ (0,2π)`, `w = nu ≥ 0`,
`w ≤ x` (θ≥u) and `x+w ≤ 2π` (n(θ+u)≤2π):  `cos x · cosh w ≤ 1`. -/
theorem tail_kernel {x w : ℝ} (hx0 : 0 < x) (hx2 : x < 2 * π) (hw0 : 0 ≤ w)
    (hlow : w ≤ x) (hhigh : x + w ≤ 2 * π) :
    cos x * cosh w ≤ 1 := by
  rcases le_or_gt (cos x) 0 with hcos | hcos
  · have : cos x * cosh w ≤ 0 := mul_nonpos_of_nonpos_of_nonneg hcos (by positivity)
    linarith
  · rcases le_or_gt x (π / 2) with hxle | hxgt
    · have hb1 : cos x ≤ exp (-(x ^ 2) / 2) := _root_.cos_le_exp_neg_half_sq (le_of_lt hx0) hxle
      have hb2 : cosh w ≤ exp (w ^ 2 / 2) := _root_.cosh_le_exp_half_sq w
      have hprod : cos x * cosh w ≤ exp (-(x ^ 2) / 2) * exp (w ^ 2 / 2) :=
        mul_le_mul hb1 hb2 (by positivity) (by positivity)
      rw [← Real.exp_add] at hprod
      have hexp : -(x ^ 2) / 2 + w ^ 2 / 2 ≤ 0 := by nlinarith [hlow, hw0, hx0]
      calc cos x * cosh w ≤ exp (-(x ^ 2) / 2 + w ^ 2 / 2) := hprod
        _ ≤ exp 0 := Real.exp_le_exp.mpr hexp
        _ = 1 := Real.exp_zero
    · have hx3 : 3 * π / 2 < x := by
        by_contra h
        push_neg at h
        have : cos x ≤ 0 :=
          Real.cos_nonpos_of_pi_div_two_le_of_le (le_of_lt hxgt) (by linarith [Real.pi_pos])
        linarith
      have hτpos : 0 < 2 * π - x := by linarith
      have hτle : 2 * π - x ≤ π / 2 := by linarith
      have hcosτ : cos x = cos (2 * π - x) := by rw [Real.cos_two_pi_sub]
      rw [hcosτ]
      have hb1 : cos (2 * π - x) ≤ exp (-((2 * π - x) ^ 2) / 2) :=
        cos_le_exp_neg_half_sq (le_of_lt hτpos) hτle
      have hb2 : cosh w ≤ exp (w ^ 2 / 2) := _root_.cosh_le_exp_half_sq w
      have hprod : cos (2 * π - x) * cosh w ≤ exp (-((2 * π - x) ^ 2) / 2) * exp (w ^ 2 / 2) :=
        mul_le_mul hb1 hb2 (by positivity) (by positivity)
      rw [← Real.exp_add] at hprod
      have hexp : -((2 * π - x) ^ 2) / 2 + w ^ 2 / 2 ≤ 0 := by nlinarith [hhigh, hw0, hτpos]
      calc cos (2 * π - x) * cosh w ≤ exp (-((2 * π - x) ^ 2) / 2 + w ^ 2 / 2) := hprod
        _ ≤ exp 0 := Real.exp_le_exp.mpr hexp
        _ = 1 := Real.exp_zero

/-- Nonnegativity of `F = 4 - 4·cos x·cosh w` (= the tail-safety `F_n ≥ 0`), fully self-contained. -/
theorem tail_F_nonneg {x w : ℝ} (hx0 : 0 < x) (hx2 : x < 2 * π) (hw0 : 0 ≤ w)
    (hlow : w ≤ x) (hhigh : x + w ≤ 2 * π) :
    0 ≤ 4 - 4 * (cos x * cosh w) := by
  have := tail_kernel hx0 hx2 hw0 hlow hhigh
  linarith
