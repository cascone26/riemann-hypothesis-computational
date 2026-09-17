import LandauZeroDistance
import Mathlib.Analysis.SpecialFunctions.Complex.LogDeriv

set_option maxHeartbeats 2000000

open Complex Metric

/-!
# §1696a — integrability along the left line, LOCAL-T zero-freeness signature (matches §1695)

Two facts enabling the left-line error integral: (a) every point of the shifted line
`Re = 1 - c₁/log(T+2)` (heights `≤ T+1`) is ζ-free — the slick argument: apply the §1686
zero-distance bound with `ρ := z` ITSELF; a zero at `z` would be at distance `0 ≥ c₁/log(T+2) > 0`
from itself. (b) Hence `t ↦ F(σ₁ + it)` (the Perron integrand) is continuous, so interval-integrable
on every `[a, b]` with heights `≤ T+1`.
-/

theorem left_line_F_intervalIntegrable' {c₁ : ℝ} (hc₁ : 0 < c₁)
    {T x a b : ℝ} (hT : 2 ≤ T) (hx : 1 ≤ x) (hc₁L : c₁ / Real.log (T + 2) < 1)
    (hzf : ∀ z : ℂ,
      z.re = 1 - c₁ / Real.log (T + 2) → |z.im| ≤ T + 1 → riemannZeta z ≠ 0)
    (hab : ∀ t ∈ Set.uIcc a b, |t| ≤ T + 1) :
    IntervalIntegrable (fun t : ℝ =>
      -logDeriv riemannZeta (((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I) *
          (x : ℂ) ^ ((((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I) + 1) /
        ((((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I) *
          ((((1 - c₁ / Real.log (T + 2) : ℝ) : ℂ) + (t : ℝ) * Complex.I) + 1)))
      MeasureTheory.volume a b := by
  set σ₁ : ℝ := 1 - c₁ / Real.log (T + 2) with hσ₁def
  have hL : 0 < Real.log (T + 2) := Real.log_pos (by linarith)
  have hσ₁pos : 0 < σ₁ := by rw [hσ₁def]; linarith
  have hσ₁lt : σ₁ < 1 := by
    rw [hσ₁def]
    have : 0 < c₁ / Real.log (T + 2) := by positivity
    linarith
  apply ContinuousOn.intervalIntegrable
  intro t ht
  apply ContinuousAt.continuousWithinAt
  set s : ℂ := ((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I with hsdef
  have hsre : s.re = σ₁ := by rw [hsdef]; simp
  have hsim : s.im = t := by rw [hsdef]; simp
  have hs1 : s ≠ 1 := by
    intro h
    have := congrArg Complex.re h
    rw [hsre] at this
    simp at this
    linarith
  have hs0 : s ≠ 0 := by
    intro h
    have := congrArg Complex.re h
    rw [hsre] at this
    simp at this
    linarith
  have hsm1 : s + 1 ≠ 0 := by
    intro h
    have := congrArg Complex.re h
    rw [Complex.add_re, hsre] at this
    simp at this
    linarith
  have hζs : riemannZeta s ≠ 0 := hzf s hsre (by rw [hsim]; exact hab t ht)
  -- continuity of the line map
  have hline : ContinuousAt (fun t : ℝ => ((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) t := by
    fun_prop
  -- continuity of logDeriv ζ at s
  have hld : ContinuousAt (logDeriv riemannZeta) s := by
    rw [logDeriv]
    apply ContinuousAt.div
    · have han : AnalyticAt ℂ riemannZeta s := by
        apply analyticOn_riemannZeta s
        simp [hs1]
      exact han.deriv.continuousAt
    · exact (differentiableAt_riemannZeta hs1).continuousAt
    · exact hζs
  -- assemble
  have h1 : ContinuousAt (fun t : ℝ =>
      logDeriv riemannZeta (((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I)) t :=
    ContinuousAt.comp (by rw [← hsdef]; exact hld) hline
  have h2 : ContinuousAt (fun t : ℝ =>
      (x : ℂ) ^ ((((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) + 1)) t := by
    apply ContinuousAt.const_cpow
    · exact hline.add continuousAt_const
    · left
      intro h
      rw [Complex.ofReal_eq_zero] at h
      linarith
  have h3 : ContinuousAt (fun t : ℝ =>
      (((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) *
        ((((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) + 1)) t :=
    hline.mul (hline.add continuousAt_const)
  have h4 : (fun t : ℝ => (((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) *
      ((((σ₁ : ℝ) : ℂ) + (t : ℝ) * Complex.I) + 1)) t ≠ 0 := by
    simp only []
    exact mul_ne_zero (by rw [← hsdef]; exact hs0) (by rw [← hsdef]; exact hsm1)
  exact ((h1.neg.mul h2).div h3 h4)

#print axioms left_line_F_intervalIntegrable'
