import Psi1CDeriv
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

set_option maxHeartbeats 1000000

/-!
# §1594 ★★ — BY-PARTS: `∫_1^N ψ(t)t^{-s-1}dt = Ψ(N)N^{-s-1} + (s+1)∫_1^N Ψ(t)t^{-s-2}dt`

The truncated SMOOTHED identity (`Ψ = psi1C`).  Right-FTC (`integral_eq_sub_of_hasDeriv_right`) applied to
`F(u) = Ψ(u)·u^{-s-1}`: `F` is continuous on `[1,N]`, has RIGHT-derivative `ψ(u)u^{-s-1} + Ψ(u)(-s-1)u^{-s-2}`
at every interior point (§1593 + `hasDerivAt_ofReal_cpow_const`, product rule — no exceptional set!), and the
derivative is interval-integrable.  `F(1) = 0` since `Ψ(1) = 0`.  Rearranged, this trades one power of `1/|s|`
for the `1/(s+1)` the Mellin-inversion `VerticalIntegrable` hypothesis needs.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory

theorem psi_cpow_intervalIntegrable_gen {r : ℂ} (hr : r.re ≤ 0) {a b : ℝ}
    (h1 : 1 ≤ a) (hab : a ≤ b) :
    IntervalIntegrable (fun t : ℝ => ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ r) MeasureTheory.volume a b := by
  constructor
  · have hfin : IsFiniteMeasure (MeasureTheory.volume.restrict (Set.Ioc a b)) :=
      ⟨by rw [MeasureTheory.Measure.restrict_apply_univ]; exact measure_Ioc_lt_top⟩
    have hfm : AEStronglyMeasurable (fun t : ℝ => ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ r)
        (MeasureTheory.volume.restrict (Set.Ioc a b)) := by
      apply AEStronglyMeasurable.mul
      · exact (Complex.measurable_ofReal.comp psi_measurable).aestronglyMeasurable
      · refine (ContinuousOn.aestronglyMeasurable ?_ measurableSet_Ioc)
        intro t ht
        rw [Set.mem_Ioc] at ht
        have ht0 : t ≠ 0 := ne_of_gt (by linarith [ht.1])
        exact (continuousAt_ofReal_cpow_const t r (Or.inr ht0)).continuousWithinAt
    apply MeasureTheory.Integrable.mono' (MeasureTheory.integrable_const (ψ b)) hfm
    rw [MeasureTheory.ae_restrict_iff' measurableSet_Ioc]
    refine Filter.Eventually.of_forall (fun t ht => ?_)
    rw [Set.mem_Ioc] at ht
    have ht0 : (0 : ℝ) < t := by linarith [ht.1]
    have ht1 : (1 : ℝ) ≤ t := by linarith [ht.1]
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (psi_nonneg t),
      Complex.norm_cpow_eq_rpow_re_of_pos ht0]
    calc ψ t * t ^ r.re
        ≤ ψ b * 1 := mul_le_mul (psi_mono ht.2)
          (Real.rpow_le_one_of_one_le_of_nonpos ht1 hr) (Real.rpow_nonneg ht0.le _) (psi_nonneg _)
      _ = ψ b := mul_one _
  · rw [Set.Ioc_eq_empty (not_lt.mpr hab)]
    exact integrableOn_empty

theorem truncated_smoothed {s : ℂ} (hs : 1 < s.re) {N : ℕ} (hN : 1 ≤ N) :
    ∫ t in (1:ℝ)..(N:ℝ), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1)
      = psi1C (N : ℝ) * (N : ℂ) ^ (-s - 1)
        + (s + 1) * ∫ t in (1:ℝ)..(N:ℝ), psi1C t * (t : ℂ) ^ (-s - 2) := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hs1ne : (-s - 1 : ℂ) ≠ 0 := by
    intro h
    have h' := congrArg Complex.re h
    simp only [Complex.sub_re, Complex.neg_re, Complex.one_re, Complex.zero_re] at h'
    linarith
  have hcpow_cont : ∀ r : ℂ, ContinuousOn (fun u : ℝ => (u : ℂ) ^ r) (Set.uIcc 1 (N : ℝ)) := by
    intro r t ht
    rw [Set.uIcc_of_le hNR, Set.mem_Icc] at ht
    have ht0 : t ≠ 0 := ne_of_gt (by linarith [ht.1])
    exact (continuousAt_ofReal_cpow_const t r (Or.inr ht0)).continuousWithinAt
  have hint1 : IntervalIntegrable (fun t : ℝ => ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1))
      MeasureTheory.volume 1 (N : ℝ) :=
    psi_cpow_intervalIntegrable_gen
      (by simp only [Complex.sub_re, Complex.neg_re, Complex.one_re]; linarith) le_rfl hNR
  have hint2 : IntervalIntegrable (fun u : ℝ => psi1C u * ((-s - 1) * (u : ℂ) ^ (-s - 2)))
      MeasureTheory.volume 1 (N : ℝ) :=
    (continuous_psi1C.continuousOn.mul (continuousOn_const.mul (hcpow_cont (-s - 2)))).intervalIntegrable
  have hcont : ContinuousOn (fun u : ℝ => psi1C u * (u : ℂ) ^ (-s - 1)) (Set.uIcc 1 (N : ℝ)) :=
    continuous_psi1C.continuousOn.mul (hcpow_cont (-s - 1))
  have hderiv : ∀ x ∈ Set.Ioo (min 1 (N : ℝ)) (max 1 (N : ℝ)),
      HasDerivWithinAt (fun u : ℝ => psi1C u * (u : ℂ) ^ (-s - 1))
        (((ψ x : ℝ) : ℂ) * (x : ℂ) ^ (-s - 1) + psi1C x * ((-s - 1) * (x : ℂ) ^ (-s - 2)))
        (Set.Ioi x) x := by
    intro x hx
    rw [min_eq_left hNR, max_eq_right hNR] at hx
    obtain ⟨hx1, _⟩ := hx
    have hd1 : HasDerivWithinAt psi1C ((ψ x : ℝ) : ℂ) (Set.Ioi x) x :=
      (psi1C_hasDerivWithinAt_right (by linarith)).mono Set.Ioi_subset_Ici_self
    have hd2 : HasDerivAt (fun u : ℝ => (u : ℂ) ^ (-s - 1))
        ((-s - 1) * (x : ℂ) ^ (-s - 1 - 1)) x :=
      hasDerivAt_ofReal_cpow_const (ne_of_gt (by linarith)) hs1ne
    have hexp : (-s - 1 - 1 : ℂ) = -s - 2 := by ring
    rw [hexp] at hd2
    exact hd1.mul hd2.hasDerivWithinAt
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDeriv_right hcont hderiv (hint1.add hint2)
  rw [intervalIntegral.integral_add hint1 hint2] at hFTC
  have hpsi1C1 : psi1C (1 : ℝ) = 0 := intervalIntegral.integral_same
  rw [hpsi1C1, zero_mul, sub_zero] at hFTC
  have hp2 : ∫ u in (1:ℝ)..(N:ℝ), psi1C u * ((-s - 1) * (u : ℂ) ^ (-s - 2))
      = (-s - 1) * ∫ u in (1:ℝ)..(N:ℝ), psi1C u * (u : ℂ) ^ (-s - 2) := by
    rw [← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr
    intro u _
    ring
  rw [hp2] at hFTC
  linear_combination hFTC

#print axioms truncated_smoothed
