import SmoothedMellin
import Mathlib.Analysis.MellinTransform

set_option maxHeartbeats 1000000

/-!
# §1597-§1598 — `Ψ` vanishes on `(0,1]`, and the smoothed identity in Mathlib `mellin` form

`Ψ(t) = ∫_1^t ψ = 0` for `t ≤ 1` (`ψ ≡ 0` on `[t,1] ⊆ [0,2)`), so `mellin Ψ (-s-1) = ∫_{Ioi 0} t^{-s-2}Ψ(t)`
has integrand vanishing on `(0,1]` — it collapses to the `Ioi 1` integral of §1596.  Hence

  `s(s+1)·mellin Ψ (-s-1) = -ζ'/ζ(s)`   (`Re s > 1`),

plugging the pillar into Mathlib's Mellin-inversion ecosystem.  EXIT0 / std-axioms.
-/

open Complex Chebyshev MeasureTheory Set

theorem psi1C_eq_zero_of_le_one {t : ℝ} (ht : t ≤ 1) : psi1C t = 0 := by
  rw [psi1C, ← intervalIntegral.integral_zero (a := (1:ℝ)) (b := t)]
  apply intervalIntegral.integral_congr
  intro u hu
  rw [Set.uIcc_comm, Set.uIcc_of_le ht, Set.mem_Icc] at hu
  simp only [Complex.ofReal_eq_zero]
  rw [psi_eq_zero_of_lt_two (by linarith [hu.2])]

theorem mellin_psi1C_eq {s : ℂ} (hs : 1 < s.re) :
    s * (s + 1) * mellin psi1C (-s - 1) = -logDeriv riemannZeta s := by
  rw [neg_logDeriv_zeta_eq_smoothed hs]
  congr 1
  rw [mellin]
  have hcongr : ∀ t ∈ Set.Ioi (0 : ℝ),
      (t : ℂ) ^ ((-s - 1) - 1) • psi1C t
      = (Set.Ioi (1 : ℝ)).indicator (fun t : ℝ => psi1C t * (t : ℂ) ^ (-s - 2)) t := by
    intro t ht
    rw [Set.mem_Ioi] at ht
    by_cases h1 : 1 < t
    · rw [Set.indicator_of_mem (Set.mem_Ioi.mpr h1), smul_eq_mul, mul_comm]
      congr 2
      ring
    · rw [Set.indicator_of_notMem (by rw [Set.mem_Ioi]; exact h1)]
      rw [not_lt] at h1
      rw [psi1C_eq_zero_of_le_one h1, smul_zero]
  rw [MeasureTheory.setIntegral_congr_fun measurableSet_Ioi hcongr,
    MeasureTheory.setIntegral_indicator measurableSet_Ioi]
  have hinter : Set.Ioi (0 : ℝ) ∩ Set.Ioi (1 : ℝ) = Set.Ioi (1 : ℝ) := by
    rw [Set.Ioi_inter_Ioi]; congr 1; exact max_eq_right (by norm_num)
  rw [hinter]

#print axioms mellin_psi1C_eq
