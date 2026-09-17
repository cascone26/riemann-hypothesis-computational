import MellinPillar
import Mathlib.Analysis.MellinTransform

set_option maxHeartbeats 1000000

/-!
# §1589 — the pillar in Mathlib `mellin` form: `-ζ'/ζ(s) = s·mellin ψ (-s)` for `Re s > 1`

Rewrites §1588 through Mathlib's `mellin` (integral over `Ioi 0`): since `ψ(t) = 0` for `t < 2`, the integrand
vanishes on `(0,1]` (indicator congruence), so `∫_{Ioi 0} = ∫_{Ioi 1}`.  This plugs the pillar into Mathlib's
Mellin ecosystem (`mellinInv_mellin_eq` inversion) — the road to the SMOOTHED explicit formula: direct `ψ`
inversion fails `VerticalIntegrable` (`~1/|y|`), but `ψ₁(x) = ∫_0^x ψ` gives `~1/|y|²` decay and continuity
everywhere, exactly Mathlib's hypotheses.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory Set

theorem neg_logDeriv_zeta_eq_smul_mellin {s : ℂ} (hs : 1 < s.re) :
    -logDeriv riemannZeta s = s * mellin (fun t : ℝ => ((ψ t : ℝ) : ℂ)) (-s) := by
  rw [neg_logDeriv_zeta_eq_mellin hs]
  congr 1
  rw [mellin]
  have hcongr : ∀ t ∈ Set.Ioi (0 : ℝ),
      (t : ℂ) ^ (-s - 1) • ((ψ t : ℝ) : ℂ)
      = (Set.Ioi (1 : ℝ)).indicator (fun t : ℝ => (t : ℂ) ^ (-s - 1) • ((ψ t : ℝ) : ℂ)) t := by
    intro t ht
    rw [Set.mem_Ioi] at ht
    by_cases h1 : 1 < t
    · rw [Set.indicator_of_mem (Set.mem_Ioi.mpr h1)]
    · rw [Set.indicator_of_notMem (by rw [Set.mem_Ioi]; exact h1)]
      rw [not_lt] at h1
      rw [psi_eq_zero_of_lt_two (by linarith)]
      simp
  rw [MeasureTheory.setIntegral_congr_fun measurableSet_Ioi hcongr,
    MeasureTheory.setIntegral_indicator measurableSet_Ioi]
  have hinter : Set.Ioi (0 : ℝ) ∩ Set.Ioi (1 : ℝ) = Set.Ioi (1 : ℝ) := by
    rw [Set.Ioi_inter_Ioi]
    congr 1
    exact max_eq_right (by norm_num)
  rw [hinter]
  apply setIntegral_congr_fun measurableSet_Ioi
  intro t _
  simp only [smul_eq_mul]
  ring

#print axioms neg_logDeriv_zeta_eq_smul_mellin
