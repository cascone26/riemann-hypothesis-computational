import Psi1Def
import PsiRightContinuous
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

set_option maxHeartbeats 1000000

/-!
# §1593 — the complex primitive `Ψ(x) = ∫_1^x ψ(t)dt : ℂ` and its RIGHT-derivative `ψ(b)` at every `b ≥ 0`

Working ℂ-valued from the start (integrand `((ψ t:ℝ):ℂ)`) dodges `ofReal`-composition lemmas.  The right-FTC
derivative (`integral_hasDerivWithinAt_right` + §1592 right-continuity) holds at EVERY `b` — the input for
integration by parts `∫_1^N ψ(t)t^{-s-1} = Ψ(N)N^{-s-1} + (s+1)∫_1^N Ψ(t)t^{-s-2}`.  EXIT0 / std-axioms.
-/

open Chebyshev MeasureTheory intervalIntegral

/-- The ℂ-valued smoothed Chebyshev primitive based at `1`. -/
noncomputable def psi1C (x : ℝ) : ℂ := ∫ t in (1 : ℝ)..x, ((ψ t : ℝ) : ℂ)

theorem psi_intervalIntegrable_C (a b : ℝ) :
    IntervalIntegrable (fun t : ℝ => ((ψ t : ℝ) : ℂ)) MeasureTheory.volume a b :=
  ⟨(psi_integrableOn_Ioc a b).ofReal, (psi_integrableOn_Ioc b a).ofReal⟩

theorem continuous_psi1C : Continuous psi1C :=
  intervalIntegral.continuous_primitive (fun a b => psi_intervalIntegrable_C a b) 1

theorem psi1C_hasDerivWithinAt_right {b : ℝ} (hb : 0 ≤ b) :
    HasDerivWithinAt psi1C ((ψ b : ℝ) : ℂ) (Set.Ici b) b := by
  apply intervalIntegral.integral_hasDerivWithinAt_right
    (psi_intervalIntegrable_C 1 b) (t := Set.Ioi b)
  · exact ⟨Set.univ, Filter.univ_mem,
      (Complex.measurable_ofReal.comp psi_measurable).aestronglyMeasurable⟩
  · exact Complex.continuous_ofReal.continuousAt.comp_continuousWithinAt
      (psi_continuousWithinAt_Ioi hb)

#print axioms psi1C_hasDerivWithinAt_right
