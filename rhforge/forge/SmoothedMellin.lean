import TruncatedSmoothed
import Psi1CBoundary
import MellinPillar
import PsiCpowIntegrableOnIoi

set_option maxHeartbeats 1000000

/-!
# §1596 ★★★ — THE SMOOTHED MELLIN IDENTITY: `-ζ'/ζ(s) = s(s+1)·∫₁^∞ Ψ(t)t^{-s-2} dt` (`Re s > 1`)

`N → ∞` in the by-parts identity §1594: boundary dies (§1595), both integrals converge (domination by
`C·t^{-σ}`), limits are unique, and §1588 identifies the left side.  The `1/(s(s+1))` decay is exactly what
Mellin inversion's `VerticalIntegrable` needs — with `Ψ` continuous, the smoothed Perron formula is now within
reach of `mellinInv_mellin_eq`.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory Filter

theorem psi1C_norm_le_real {t : ℝ} (ht : 1 ≤ t) :
    ‖psi1C t‖ ≤ (Real.log 4 + 4) * t ^ (2 : ℝ) := by
  have h : ‖∫ u in (1:ℝ)..t, ((ψ u : ℝ) : ℂ)‖ ≤ (Real.log 4 + 4) * t * |t - 1| := by
    apply intervalIntegral.norm_integral_le_of_norm_le_const
    intro u hu
    rw [Set.uIoc_of_le ht, Set.mem_Ioc] at hu
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (psi_nonneg u)]
    calc ψ u ≤ ψ t := psi_mono hu.2
      _ ≤ (Real.log 4 + 4) * t := psi_le_const_mul_self (by linarith)
  calc ‖psi1C t‖
      ≤ (Real.log 4 + 4) * t * |t - 1| := h
    _ ≤ (Real.log 4 + 4) * t * t := by
        have habs : |t - 1| ≤ t := by rw [abs_of_nonneg (by linarith)]; linarith
        have hc : (0 : ℝ) ≤ (Real.log 4 + 4) * t := by
          have := Real.log_nonneg (by norm_num : (1:ℝ) ≤ 4)
          positivity
        exact mul_le_mul_of_nonneg_left habs hc
    _ = (Real.log 4 + 4) * t ^ (2 : ℝ) := by
        rw [show ((2:ℝ)) = ((2:ℕ):ℝ) by norm_num, Real.rpow_natCast]
        ring

theorem psi1C_cpow_integrableOn_Ioi {s : ℂ} (hs : 1 < s.re) :
    IntegrableOn (fun t : ℝ => psi1C t * (t : ℂ) ^ (-s - 2)) (Set.Ioi (1 : ℝ)) := by
  have hfm : AEStronglyMeasurable (fun t : ℝ => psi1C t * (t : ℂ) ^ (-s - 2))
      (MeasureTheory.volume.restrict (Set.Ioi (1 : ℝ))) := by
    apply AEStronglyMeasurable.mul
    · exact continuous_psi1C.aestronglyMeasurable
    · refine (ContinuousOn.aestronglyMeasurable ?_ measurableSet_Ioi)
      intro t ht
      rw [Set.mem_Ioi] at ht
      have ht0 : t ≠ 0 := by linarith
      exact (continuousAt_ofReal_cpow_const t (-s - 2) (Or.inr ht0)).continuousWithinAt
  have hg : IntegrableOn (fun t : ℝ => (Real.log 4 + 4) * t ^ (-s.re)) (Set.Ioi (1 : ℝ)) :=
    (integrableOn_Ioi_rpow_of_lt (by linarith : -s.re < -1) one_pos).const_mul _
  apply MeasureTheory.Integrable.mono' hg hfm
  rw [MeasureTheory.ae_restrict_iff' measurableSet_Ioi]
  refine Filter.Eventually.of_forall (fun t ht => ?_)
  rw [Set.mem_Ioi] at ht
  have ht0 : (0 : ℝ) < t := by linarith
  rw [norm_mul, Complex.norm_cpow_eq_rpow_re_of_pos ht0]
  have hexp : ((-s - 2 : ℂ)).re = -s.re - 2 := by
    simp [Complex.sub_re, Complex.neg_re]
  rw [hexp]
  calc ‖psi1C t‖ * t ^ (-s.re - 2)
      ≤ ((Real.log 4 + 4) * t ^ (2 : ℝ)) * t ^ (-s.re - 2) :=
        mul_le_mul_of_nonneg_right (psi1C_norm_le_real (le_of_lt ht)) (Real.rpow_nonneg ht0.le _)
    _ = (Real.log 4 + 4) * t ^ (-s.re) := by
        rw [mul_assoc, ← Real.rpow_add ht0]
        congr 2
        ring

theorem neg_logDeriv_zeta_eq_smoothed {s : ℂ} (hs : 1 < s.re) :
    -logDeriv riemannZeta s
      = s * (s + 1) * ∫ t in Set.Ioi (1 : ℝ), psi1C t * (t : ℂ) ^ (-s - 2) := by
  have hI2 := MeasureTheory.intervalIntegral_tendsto_integral_Ioi (1 : ℝ)
    (psi1C_cpow_integrableOn_Ioi hs) tendsto_natCast_atTop_atTop
  have hLHS := MeasureTheory.intervalIntegral_tendsto_integral_Ioi (1 : ℝ)
    (psi_cpow_integrableOn_Ioi hs) tendsto_natCast_atTop_atTop
  have hRHS : Tendsto (fun N : ℕ => psi1C (N : ℝ) * (N : ℂ) ^ (-s - 1)
      + (s + 1) * ∫ t in (1:ℝ)..(N:ℝ), psi1C t * (t : ℂ) ^ (-s - 2)) atTop
      (nhds (0 + (s + 1) * ∫ t in Set.Ioi (1 : ℝ), psi1C t * (t : ℂ) ^ (-s - 2))) :=
    (psi1C_boundary_tendsto_zero hs).add (hI2.const_mul (s + 1))
  have hcongr : Tendsto (fun N : ℕ => ∫ t in (1:ℝ)..(N:ℝ), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1))
      atTop (nhds (0 + (s + 1) * ∫ t in Set.Ioi (1 : ℝ), psi1C t * (t : ℂ) ^ (-s - 2))) := by
    apply hRHS.congr'
    filter_upwards [eventually_ge_atTop 1] with N hN
    exact (truncated_smoothed hs hN).symm
  have huniq := tendsto_nhds_unique hLHS hcongr
  rw [neg_logDeriv_zeta_eq_mellin hs, huniq]
  ring

#print axioms neg_logDeriv_zeta_eq_smoothed
