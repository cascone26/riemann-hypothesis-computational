import CCMWeightGram

/-!
# S2943 — shape of the CCM/Shi weight: Fejér form, periodicity, quadratic small-γ bound (lane 2)

w_γ = 4 sin²(γL/2)/L.  Here:
* `ccmW_eq_cos` : w_γ = 2(1 − cos γL)/L  (the Fejér / 1−cos kernel form);
* `ccmW_periodic` : w is (2π/L)-periodic in γ (resonance lattice of S2939);
* `ccmW_le_sq` : w_γ ≤ L·γ²  (quadratic vanishing at γ = 0, from |sin x| ≤ |x|);
* `ccmW_le_min` : w_γ ≤ min(4/L, Lγ²) — the two-regime envelope used to split low-lying vs high zeros
  in the zero-side Gram (S2931/S2937 brackets use only the 4/L half).
-/

namespace CCMWeight

theorem ccmW_eq_cos (L γ : ℝ) : ccmW L γ = 2 * (1 - Real.cos (γ * L)) / L := by
  unfold ccmW
  have h := Real.cos_two_mul (γ * L / 2)
  rw [show 2 * (γ * L / 2) = γ * L by ring] at h
  have hs := Real.sin_sq_add_cos_sq (γ * L / 2)
  congr 1
  nlinarith [h, hs]

theorem ccmW_periodic {L : ℝ} (hL : L ≠ 0) (γ : ℝ) :
    ccmW L (γ + 2 * Real.pi / L) = ccmW L γ := by
  rw [ccmW_eq_cos, ccmW_eq_cos]
  congr 3
  rw [show (γ + 2 * Real.pi / L) * L = γ * L + 2 * Real.pi by field_simp]
  exact Real.cos_add_two_pi _

theorem ccmW_le_sq {L : ℝ} (hL : 0 < L) (γ : ℝ) : ccmW L γ ≤ L * γ ^ 2 := by
  unfold ccmW
  rw [div_le_iff₀ hL]
  have h := Real.abs_sin_le_abs (x := γ * L / 2)
  have h2 : Real.sin (γ * L / 2) ^ 2 ≤ (γ * L / 2) ^ 2 := by
    rw [← sq_abs, ← sq_abs (γ * L / 2)]
    exact pow_le_pow_left₀ (abs_nonneg _) h 2
  nlinarith [h2]

theorem ccmW_le_min {L : ℝ} (hL : 0 < L) (γ : ℝ) : ccmW L γ ≤ min (4 / L) (L * γ ^ 2) :=
  le_min (ccmW_le hL γ) (ccmW_le_sq hL γ)

end CCMWeight

#print axioms CCMWeight.ccmW_eq_cos
#print axioms CCMWeight.ccmW_periodic
#print axioms CCMWeight.ccmW_le_sq
#print axioms CCMWeight.ccmW_le_min
