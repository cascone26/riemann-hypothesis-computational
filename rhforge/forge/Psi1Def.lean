import PsiCpowIntervalIntegrable
import Mathlib.MeasureTheory.Integral.DominatedConvergence

set_option maxHeartbeats 1000000

/-!
# §1590 — the smoothed Chebyshev function `ψ₁(x) = ∫_0^x ψ(t) dt`: definition, integrability, CONTINUITY

Foundation of the SMOOTHED PERRON arc.  `ψ` is integrable on every finite interval (measurable §1582 +
monotone-bounded + finite measure), so its primitive `ψ₁` is well-defined and CONTINUOUS everywhere
(`continuous_primitive`) — exactly the `ContinuousAt` hypothesis Mathlib's Mellin inversion needs at every `x`
(where the unsmoothed step function `ψ` fails at prime powers).  EXIT0 / std-axioms.
-/

open Chebyshev MeasureTheory

theorem psi_integrableOn_Ioc (a b : ℝ) :
    IntegrableOn (fun t : ℝ => ψ t) (Set.Ioc a b) := by
  by_cases hab : b ≤ a
  · rw [Set.Ioc_eq_empty (not_lt.mpr hab)]
    exact integrableOn_empty
  · have hfin : IsFiniteMeasure (MeasureTheory.volume.restrict (Set.Ioc a b)) :=
      ⟨by rw [MeasureTheory.Measure.restrict_apply_univ]; exact measure_Ioc_lt_top⟩
    apply MeasureTheory.Integrable.mono' (MeasureTheory.integrable_const (ψ b))
    · exact psi_measurable.aestronglyMeasurable
    · rw [MeasureTheory.ae_restrict_iff' measurableSet_Ioc]
      refine Filter.Eventually.of_forall (fun t ht => ?_)
      rw [Set.mem_Ioc] at ht
      rw [Real.norm_eq_abs, abs_of_nonneg (psi_nonneg t)]
      exact psi_mono ht.2

theorem psi_intervalIntegrable (a b : ℝ) :
    IntervalIntegrable (fun t : ℝ => ψ t) MeasureTheory.volume a b :=
  ⟨psi_integrableOn_Ioc a b, psi_integrableOn_Ioc b a⟩

/-- The smoothed Chebyshev function `ψ₁(x) = ∫_0^x ψ(t) dt`. -/
noncomputable def psi1 (x : ℝ) : ℝ := ∫ t in (0 : ℝ)..x, ψ t

theorem continuous_psi1 : Continuous psi1 :=
  intervalIntegral.continuous_primitive (fun a b => psi_intervalIntegrable a b) 0

#print axioms continuous_psi1
