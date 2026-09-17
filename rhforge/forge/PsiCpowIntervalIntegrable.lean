import Mathlib.NumberTheory.Chebyshev
import Mathlib.MeasureTheory.Function.Floor
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Complex
import Mathlib.Analysis.SpecialFunctions.Pow.Continuity

set_option maxHeartbeats 1000000

/-!
# §1582 — strip integrability: `t ↦ ψ(t)·t^{-s-1}` is interval-integrable on `[k, k+1]` (`k ≥ 1`)

Fifth Mellin brick.  `ψ` is measurable (`ψ = (ψ ∘ ↑) ∘ ⌊·⌋₊`, `Nat.measurable_floor`, discrete source); the cpow
factor is continuous on the strip (positive real base, `slitPlane`); the product is a.e. strongly measurable on
the strip and bounded by `ψ(k+1)` (`ψ` monotone, `t^{-σ-1} ≤ 1` for `t ≥ 1`, `σ ≥ 0`), so it is integrable on the
finite-measure strip.  This is the integrability input for gluing the strips into `∫_1^N`.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory

theorem psi_measurable : Measurable (fun t : ℝ => ψ t) := by
  have heq : (fun t : ℝ => ψ t) = (fun n : ℕ => ψ (n : ℝ)) ∘ Nat.floor := by
    funext t
    simp only [Function.comp_apply]
    rw [psi_eq_psi_coe_floor]
  rw [heq]
  exact measurable_from_top.comp Nat.measurable_floor

theorem psi_cpow_intervalIntegrable {s : ℂ} (hσ : 0 ≤ s.re) {k : ℕ} (hk : 0 < k) :
    IntervalIntegrable (fun t : ℝ => ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1)) MeasureTheory.volume
      (k : ℝ) ((k : ℝ) + 1) := by
  have hkR : (1 : ℝ) ≤ (k : ℝ) := by exact_mod_cast hk
  rw [intervalIntegrable_iff_integrableOn_Ioc_of_le (by linarith : (k : ℝ) ≤ (k : ℝ) + 1)]
  have hcont : ContinuousOn (fun t : ℝ => (t : ℂ) ^ (-s - 1)) (Set.Ioc (k : ℝ) ((k : ℝ) + 1)) := by
    intro t ht
    rw [Set.mem_Ioc] at ht
    have ht0 : t ≠ 0 := by intro h; rw [h] at ht; linarith [ht.1]
    exact (continuousAt_ofReal_cpow_const t (-s - 1) (Or.inr ht0)).continuousWithinAt
  have hfm : AEStronglyMeasurable (fun t : ℝ => ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1))
      (MeasureTheory.volume.restrict (Set.Ioc (k : ℝ) ((k : ℝ) + 1))) := by
    apply AEStronglyMeasurable.mul
    · exact (Complex.measurable_ofReal.comp psi_measurable).aestronglyMeasurable
    · exact hcont.aestronglyMeasurable measurableSet_Ioc
  have : IsFiniteMeasure (MeasureTheory.volume.restrict (Set.Ioc (k : ℝ) ((k : ℝ) + 1))) :=
    ⟨by rw [MeasureTheory.Measure.restrict_apply_univ]; exact measure_Ioc_lt_top⟩
  apply MeasureTheory.Integrable.mono'
    (MeasureTheory.integrable_const (ψ ((k : ℝ) + 1)))
  · exact hfm
  · rw [MeasureTheory.ae_restrict_iff' measurableSet_Ioc]
    refine Filter.Eventually.of_forall (fun t ht => ?_)
    rw [Set.mem_Ioc] at ht
    have ht0 : 0 < t := by linarith [ht.1]
    have ht1 : 1 ≤ t := by linarith [ht.1]
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (psi_nonneg t),
      Complex.norm_cpow_eq_rpow_re_of_pos ht0]
    have h1 : ψ t ≤ ψ ((k : ℝ) + 1) := psi_mono ht.2
    have h2 : t ^ ((-s - 1).re) ≤ 1 := by
      apply Real.rpow_le_one_of_one_le_of_nonpos ht1
      simp only [Complex.sub_re, Complex.neg_re, Complex.one_re]
      linarith
    calc ψ t * t ^ ((-s - 1).re)
        ≤ ψ ((k : ℝ) + 1) * 1 :=
          mul_le_mul h1 h2 (Real.rpow_nonneg ht0.le _) (psi_nonneg _)
      _ = ψ ((k : ℝ) + 1) := mul_one _

#print axioms psi_cpow_intervalIntegrable
