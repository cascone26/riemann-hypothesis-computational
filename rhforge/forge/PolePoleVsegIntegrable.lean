import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic
import Mathlib.Analysis.SpecialFunctions.Pow.Complex

set_option maxHeartbeats 1000000

/-!
# §1623 — vertical-segment integrability of the pole term `(x²/2)/(s-1)`

Along a vertical segment `t ↦ σ + it` with `σ ≠ 1`, the pole term `(x²/2)/((σ+it)-1)` is continuous (its denominator
never vanishes, having real part `σ - 1 ≠ 0`), hence interval-integrable.  This is the pole-side integrability input
for the Perron line-integral split `∫ F = ∫ pole + ∫ R`.  EXIT0 / std-axioms.
-/

open Complex MeasureTheory

theorem intervalIntegrable_pole_vseg {x : ℝ} {σ a b : ℝ} (hσ : σ ≠ 1) :
    IntervalIntegrable
      (fun t : ℝ => (x : ℂ) ^ (2 : ℂ) / 2 / (((σ : ℂ) + (t : ℝ) * I) - 1)) volume a b := by
  apply ContinuousOn.intervalIntegrable
  intro t _ht
  have hden : (σ : ℂ) + (t : ℝ) * I - 1 ≠ 0 := by
    intro h
    have hre : σ - 1 = 0 := by have := congrArg Complex.re h; simpa using this
    exact hσ (by linarith)
  have hg : ContinuousAt (fun t : ℝ => (σ : ℂ) + (t : ℝ) * I - 1) t := by fun_prop
  exact (continuousAt_const.div hg hden).continuousWithinAt

#print axioms intervalIntegrable_pole_vseg
