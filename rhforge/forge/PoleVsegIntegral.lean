import PolePoleVsegIntegrable
import Mathlib.Analysis.SpecialFunctions.Complex.LogDeriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

set_option maxHeartbeats 1000000

/-!
# §1626 — the pole-term integral in closed form: `∫ (x²/2)/(s-1) = (x²/2)(-i)·[log w(b) - log w(a)]`

For `σ > 1` the path `t ↦ (σ+it)-1` stays in the right half-plane (`Re = σ-1 > 0`), inside the slit plane where
`Complex.log` is holomorphic, so `G(t) = (x²/2)(-i)·log((σ+it)-1)` is an antiderivative of `(x²/2)/((σ+it)-1)`
(`HasDerivAt.clog_real`, `(-i)·i = 1`).  The fundamental theorem of calculus then evaluates the pole-term integral
exactly.  As `T→∞` this closed form supplies the `iπ`-jump (change of argument) that becomes the residue in the
contour shift.  EXIT0 / std-axioms.
-/

open Complex MeasureTheory intervalIntegral

theorem pole_vseg_integral {x : ℝ} {σ a b : ℝ} (hσ : 1 < σ) :
    (∫ t : ℝ in a..b, (x : ℂ) ^ (2 : ℂ) / 2 / (((σ : ℂ) + (t : ℝ) * I) - 1))
      = (x : ℂ) ^ (2 : ℂ) / 2 * (-I) * Complex.log ((σ : ℂ) + (b : ℝ) * I - 1)
        - (x : ℂ) ^ (2 : ℂ) / 2 * (-I) * Complex.log ((σ : ℂ) + (a : ℝ) * I - 1) := by
  have hII : (-I) * I = 1 := by rw [neg_mul, Complex.I_mul_I, neg_neg]
  -- antiderivative has the pole term as derivative at every point
  have hG : ∀ t : ℝ, HasDerivAt
      (fun t : ℝ => (x : ℂ) ^ (2 : ℂ) / 2 * (-I) * Complex.log ((σ : ℂ) + (t : ℝ) * I - 1))
      ((x : ℂ) ^ (2 : ℂ) / 2 / (((σ : ℂ) + (t : ℝ) * I) - 1)) t := by
    intro t
    have hw : HasDerivAt (fun t : ℝ => (σ : ℂ) + (t : ℝ) * I - 1) I t := by
      have hoRt : HasDerivAt (fun t : ℝ => ((t : ℝ) : ℂ)) 1 t := by
        simpa using (hasDerivAt_id t).ofReal_comp
      have := ((hoRt.mul_const I).const_add (σ : ℂ)).sub_const 1
      simpa using this
    have hslit : (σ : ℂ) + (t : ℝ) * I - 1 ∈ slitPlane := by
      apply mem_slitPlane_iff.mpr; left
      simp [Complex.sub_re, Complex.add_re, Complex.mul_re]
      linarith
    have hclog := hw.clog_real hslit
    have h := hclog.const_mul ((x : ℂ) ^ (2 : ℂ) / 2 * (-I))
    have heq : (x : ℂ) ^ (2 : ℂ) / 2 * (-I) * (I / ((σ : ℂ) + (t : ℝ) * I - 1))
        = (x : ℂ) ^ (2 : ℂ) / 2 / (((σ : ℂ) + (t : ℝ) * I) - 1) := by
      rw [mul_assoc, ← mul_div_assoc, hII, mul_one_div]
    rw [heq] at h
    exact h
  rw [integral_eq_sub_of_hasDerivAt (fun t _ => hG t) (intervalIntegrable_pole_vseg hσ.ne')]

#print axioms pole_vseg_integral
