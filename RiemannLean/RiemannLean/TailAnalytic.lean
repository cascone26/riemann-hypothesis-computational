import Mathlib

/-!
# The two analytic inputs (B1), (B2) for the route-6 tail lemma — self-contained proofs.

(B1)  `cos t ≤ exp(-(t²)/2)` on `[0, π/2]`, via `g x = log (cos x) + x²/2` strictly antitone
      (`g' = x - tan x < 0`, using `Real.lt_tan`).
(B2)  `cosh z ≤ exp(z²/2)` for all `z`, via `h x = log (cosh x) - x²/2` (even; `h' = tanh x - x < 0`).
-/

open Real Set

lemma cos_le_exp_neg_half_sq {t : ℝ} (h0 : 0 ≤ t) (h1 : t ≤ π / 2) :
    cos t ≤ exp (-(t ^ 2) / 2) := by
  rcases eq_or_lt_of_le h0 with h | hpos
  · simp [← h]
  rcases eq_or_lt_of_le h1 with h | hlt
  · subst h
    rw [Real.cos_pi_div_two]; positivity
  · have hcospos : 0 < cos t := Real.cos_pos_of_mem_Ioo ⟨by linarith [Real.pi_pos], hlt⟩
    have hlog : Real.log (cos t) ≤ -(t ^ 2) / 2 := by
      set g : ℝ → ℝ := fun x => Real.log (cos x) + x ^ 2 / 2 with hg
      have hg0 : g 0 = 0 := by simp [hg]
      have hcont : ContinuousOn g (Icc 0 t) := by
        apply ContinuousOn.add _ (by fun_prop)
        apply ContinuousOn.log (Real.continuous_cos.continuousOn)
        intro x hx
        exact ne_of_gt (Real.cos_pos_of_mem_Ioo ⟨by linarith [hx.1, Real.pi_pos], by linarith [hx.2]⟩)
      have hanti : StrictAntiOn g (Icc 0 t) := by
        apply strictAntiOn_of_hasDerivWithinAt_neg (convex_Icc 0 t) hcont
        · intro x hx
          rw [interior_Icc] at hx
          have hxc : 0 < cos x :=
            Real.cos_pos_of_mem_Ioo ⟨by linarith [hx.1, Real.pi_pos], by linarith [hx.2]⟩
          have hc : HasDerivAt (fun x => Real.log (cos x)) (-tan x) x := by
            have := (Real.hasDerivAt_log (ne_of_gt hxc)).comp x (Real.hasDerivAt_cos x)
            simp only [Function.comp] at this
            convert this using 1
            rw [Real.tan_eq_sin_div_cos]; ring
          have h2 : HasDerivAt (fun x => x ^ 2 / 2) x x := by
            simpa using ((hasDerivAt_pow 2 x).div_const 2)
          have : HasDerivAt g (-tan x + x) x := hc.add h2
          exact this.hasDerivWithinAt
        · intro x hx
          rw [interior_Icc] at hx
          have := Real.lt_tan hx.1 (by linarith [hx.2] : x < π / 2)
          linarith
      have hlt' : g t < g 0 := hanti (left_mem_Icc.mpr h0) (right_mem_Icc.mpr h0) hpos
      rw [hg0] at hlt'
      simp only [hg] at hlt'
      linarith
    calc cos t = exp (Real.log (cos t)) := (Real.exp_log hcospos).symm
      _ ≤ exp (-(t ^ 2) / 2) := Real.exp_le_exp.mpr hlog

lemma cosh_le_exp_half_sq (z : ℝ) : cosh z ≤ exp (z ^ 2 / 2) := by
  wlog hz : 0 ≤ z generalizing z with H
  · have := H (-z) (by linarith [not_le.mp hz])
    simpa [Real.cosh_neg, even_two.neg_pow] using this
  have hlog : Real.log (cosh z) ≤ z ^ 2 / 2 := by
    set h : ℝ → ℝ := fun x => Real.log (cosh x) - x ^ 2 / 2 with hh
    have hh0 : h 0 = 0 := by simp [hh]
    have hcont : ContinuousOn h (Icc 0 z) := by
      apply ContinuousOn.sub _ (by fun_prop)
      exact ContinuousOn.log Real.continuous_cosh.continuousOn
        (fun x _ => ne_of_gt (Real.cosh_pos x))
    have hanti : StrictAntiOn h (Icc 0 z) := by
      apply strictAntiOn_of_hasDerivWithinAt_neg (convex_Icc 0 z) hcont
      · intro x hx
        rw [interior_Icc] at hx
        have hc : HasDerivAt (fun x => Real.log (cosh x)) (tanh x) x := by
          have := (Real.hasDerivAt_log (ne_of_gt (Real.cosh_pos x))).comp x (Real.hasDerivAt_cosh x)
          simp only [Function.comp] at this
          convert this using 1
          rw [Real.tanh_eq_sinh_div_cosh]; ring
        have h2 : HasDerivAt (fun x => x ^ 2 / 2) x x := by
          simpa using ((hasDerivAt_pow 2 x).div_const 2)
        exact (hc.sub h2).hasDerivWithinAt
      · intro x hx
        rw [interior_Icc] at hx
        have hx0 : 0 < x := hx.1
        have hkey : sinh x < x * cosh x := by
          have hcontφ : ContinuousOn (fun y => y * cosh y - sinh y) (Icc 0 x) :=
            ((continuous_id.mul Real.continuous_cosh).sub Real.continuous_sinh).continuousOn
          have hmono : StrictMonoOn (fun y => y * cosh y - sinh y) (Icc 0 x) := by
            apply strictMonoOn_of_hasDerivWithinAt_pos (convex_Icc 0 x) hcontφ
            · intro y hy
              rw [interior_Icc] at hy
              have := (((hasDerivAt_id y).mul (Real.hasDerivAt_cosh y)).sub
                (Real.hasDerivAt_sinh y))
              exact this.hasDerivWithinAt
            · intro y hy
              rw [interior_Icc] at hy
              have : 0 < y * sinh y := mul_pos hy.1 (Real.sinh_pos_iff.mpr hy.1)
              simp only [id_eq]
              nlinarith [this]
          have := hmono (left_mem_Icc.mpr hx0.le) (right_mem_Icc.mpr hx0.le) hx0
          simpa using this
        have htanh : tanh x < x := by
          rw [Real.tanh_eq_sinh_div_cosh, div_lt_iff₀ (Real.cosh_pos x)]; exact hkey
        linarith
    rcases eq_or_lt_of_le hz with h | hzpos
    · simp [← h]
    · have hlt : h z < h 0 := hanti (left_mem_Icc.mpr hz) (right_mem_Icc.mpr hz) hzpos
      rw [hh0] at hlt; simp only [hh] at hlt; linarith
  calc cosh z = exp (Real.log (cosh z)) := (Real.exp_log (Real.cosh_pos z)).symm
    _ ≤ exp (z ^ 2 / 2) := Real.exp_le_exp.mpr hlog
