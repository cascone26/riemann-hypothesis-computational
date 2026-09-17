import PsiCpowIntervalIntegrable
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals

set_option maxHeartbeats 1000000

/-!
# §1587 — global integrability: `t ↦ ψ(t)·t^{-s-1}` is integrable on `(1, ∞)` for `Re s > 1`

`‖ψ(t)t^{-s-1}‖ = ψ(t)·t^{-σ-1} ≤ (log 4+4)·t·t^{-σ-1} = (log 4+4)·t^{-σ}` (Chebyshev `psi_le_const_mul_self`),
and `t^{-σ}` is integrable on `(1,∞)` for `σ > 1` (`integrableOn_Ioi_rpow_of_lt`).  Measurability as in §1582.
The input `intervalIntegral_tendsto_integral_Ioi` needs for `∫_1^N → ∫_1^∞`.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory Set

theorem psi_cpow_integrableOn_Ioi {s : ℂ} (hs : 1 < s.re) :
    IntegrableOn (fun t : ℝ => ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1)) (Set.Ioi (1 : ℝ)) := by
  have hfm : AEStronglyMeasurable (fun t : ℝ => ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1))
      (MeasureTheory.volume.restrict (Set.Ioi (1 : ℝ))) := by
    apply AEStronglyMeasurable.mul
    · exact (Complex.measurable_ofReal.comp psi_measurable).aestronglyMeasurable
    · refine (ContinuousOn.aestronglyMeasurable ?_ measurableSet_Ioi)
      intro t ht
      rw [Set.mem_Ioi] at ht
      have ht0 : t ≠ 0 := by linarith
      exact (continuousAt_ofReal_cpow_const t (-s - 1) (Or.inr ht0)).continuousWithinAt
  have hg : IntegrableOn (fun t : ℝ => (Real.log 4 + 4) * t ^ (-s.re)) (Set.Ioi (1 : ℝ)) :=
    (integrableOn_Ioi_rpow_of_lt (by linarith : -s.re < -1) one_pos).const_mul _
  apply MeasureTheory.Integrable.mono' hg hfm
  rw [MeasureTheory.ae_restrict_iff' measurableSet_Ioi]
  refine Filter.Eventually.of_forall (fun t ht => ?_)
  rw [Set.mem_Ioi] at ht
  have ht0 : (0 : ℝ) < t := by linarith
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (psi_nonneg t),
    Complex.norm_cpow_eq_rpow_re_of_pos ht0]
  have hψ : ψ t ≤ (Real.log 4 + 4) * t := psi_le_const_mul_self ht0.le
  have hexp : ((-s - 1 : ℂ)).re = -s.re - 1 := by
    simp [Complex.sub_re, Complex.neg_re]
  calc ψ t * t ^ ((-s - 1 : ℂ)).re
      ≤ ((Real.log 4 + 4) * t) * t ^ ((-s - 1 : ℂ)).re :=
        mul_le_mul_of_nonneg_right hψ (Real.rpow_nonneg ht0.le _)
    _ = (Real.log 4 + 4) * t ^ (-s.re) := by
        rw [hexp, show -s.re - 1 = -s.re + (-1 : ℝ) by ring, Real.rpow_add ht0,
          Real.rpow_neg_one]
        field_simp
#print axioms psi_cpow_integrableOn_Ioi
