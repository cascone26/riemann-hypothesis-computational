import Mathlib

/-!
# S3022 — Blaschke factor toolkit (lane 1, κ 2 → 1 program, brick 1)

Plan (Report): replace the two-radii Titchmarsh construction (zeros in the half disk, BC on R/2, sphere bound on R) by the
Jensen/Blaschke one-radius construction: divide ζ by b_ρ(z) = R(w − a)/(R² − conj(a)·w) (w = z − c, a = ρ − c) for every
zero in the disk, so g = ζ/∏b is zero-free on the WHOLE disk with |g| = |ζ| on the sphere; then sharp BC (S3018) on
radius R gives growth (2/R)·log(M/|ζ(c)|) instead of (2/(R/2))·…, i.e. κ = 2 → 1.  The zero terms become
Re 1/(c−ρ) − (σ−β)/R² ≥ 0 (since |c − ρ| ≤ R), so all but the chosen zero may still be dropped.
This brick — the algebra of one factor:
* `blaschke_identity`  — |R² − ā w|² − R²|w − a|² = (R² − |a|²)(R² − |w|²).
* `blaschke_den_ne`    — |a| < R, |w| ≤ R ⇒ R² − ā w ≠ 0.
* `blaschke_norm_le_one` — |a| < R, |w| ≤ R ⇒ ‖b w‖ ≤ 1.
* `blaschke_norm_sphere` — |a| < R, |w| = R ⇒ ‖b w‖ = 1.
* `blaschke_logDeriv_zero` — logDeriv b 0 = −1/a + ā/R²  (a ≠ 0).
* `blaschke_zero_term_nonneg` — |a| ≤ R, a ≠ 0, Re a ≤ 0 ⇒ 0 ≤ Re(1/(−a)) + Re a/R².
-/

open Complex

namespace Blaschke

noncomputable def b (R : ℝ) (a w : ℂ) : ℂ := (R : ℂ) * (w - a) / ((R : ℂ) ^ 2 - (starRingEnd ℂ) a * w)

theorem blaschke_identity (R : ℝ) (a w : ℂ) :
    ‖(R : ℂ) ^ 2 - (starRingEnd ℂ) a * w‖ ^ 2 - R ^ 2 * ‖w - a‖ ^ 2
      = (R ^ 2 - ‖a‖ ^ 2) * (R ^ 2 - ‖w‖ ^ 2) := by
  simp only [Complex.sq_norm, Complex.normSq_apply]
  simp [pow_two]
  ring

theorem blaschke_den_ne {R : ℝ} {a w : ℂ} (ha : ‖a‖ < R) (hw : ‖w‖ ≤ R) :
    (R : ℂ) ^ 2 - (starRingEnd ℂ) a * w ≠ 0 := by
  intro h
  have h1 : ‖(starRingEnd ℂ) a * w‖ = ‖(R : ℂ) ^ 2‖ := by rw [sub_eq_zero.mp h]
  rw [norm_mul, Complex.norm_conj, norm_pow, Complex.norm_real, Real.norm_eq_abs] at h1
  have hR : 0 < R := lt_of_le_of_lt (norm_nonneg _) ha
  rw [abs_of_pos hR] at h1
  have : ‖a‖ * ‖w‖ < R ^ 2 := by
    calc ‖a‖ * ‖w‖ ≤ ‖a‖ * R := mul_le_mul_of_nonneg_left hw (norm_nonneg _)
      _ < R * R := mul_lt_mul_of_pos_right ha hR
      _ = R ^ 2 := by ring
  linarith

theorem blaschke_norm_le_one {R : ℝ} {a w : ℂ} (ha : ‖a‖ < R) (hw : ‖w‖ ≤ R) : ‖b R a w‖ ≤ 1 := by
  have hR : 0 < R := lt_of_le_of_lt (norm_nonneg _) ha
  have hden := blaschke_den_ne ha hw
  unfold b
  rw [norm_div, div_le_one (norm_pos_iff.mpr hden), norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hR]
  have hid := blaschke_identity R a w
  have hprod : 0 ≤ (R ^ 2 - ‖a‖ ^ 2) * (R ^ 2 - ‖w‖ ^ 2) := by
    apply mul_nonneg <;> nlinarith [norm_nonneg a, norm_nonneg w]
  have hsq : (R * ‖w - a‖) ^ 2 ≤ ‖(R : ℂ) ^ 2 - (starRingEnd ℂ) a * w‖ ^ 2 := by nlinarith
  exact (sq_le_sq₀ (by positivity) (norm_nonneg _)).mp hsq

theorem blaschke_norm_sphere {R : ℝ} {a w : ℂ} (ha : ‖a‖ < R) (hw : ‖w‖ = R) : ‖b R a w‖ = 1 := by
  have hR : 0 < R := lt_of_le_of_lt (norm_nonneg _) ha
  have hden := blaschke_den_ne ha hw.le
  unfold b
  rw [norm_div, div_eq_one_iff_eq (norm_ne_zero_iff.mpr hden), norm_mul, Complex.norm_real, Real.norm_eq_abs,
    abs_of_pos hR]
  have hid := blaschke_identity R a w
  rw [hw, sub_self, mul_zero] at hid
  have hsq : (R * ‖w - a‖) ^ 2 = ‖(R : ℂ) ^ 2 - (starRingEnd ℂ) a * w‖ ^ 2 := by nlinarith
  exact (sq_eq_sq₀ (by positivity) (norm_nonneg _)).mp hsq

theorem blaschke_logDeriv_zero {R : ℝ} {a : ℂ} (ha : ‖a‖ < R) (ha0 : a ≠ 0) :
    logDeriv (b R a) 0 = -1 / a + (starRingEnd ℂ) a / (R : ℂ) ^ 2 := by
  have hR : 0 < R := lt_of_le_of_lt (norm_nonneg _) ha
  have hR2 : (R : ℂ) ^ 2 ≠ 0 := by exact_mod_cast (by positivity : R ^ 2 ≠ 0)
  have hRc : (R : ℂ) ≠ 0 := by exact_mod_cast hR.ne'
  have hnum : HasDerivAt (fun w : ℂ => (R : ℂ) * (w - a)) (R : ℂ) 0 := by
    simpa using ((hasDerivAt_id (0 : ℂ)).sub_const a).const_mul (R : ℂ)
  have hden : HasDerivAt (fun w : ℂ => (R : ℂ) ^ 2 - (starRingEnd ℂ) a * w) (-(starRingEnd ℂ) a) 0 := by
    simpa using ((hasDerivAt_id (0 : ℂ)).const_mul ((starRingEnd ℂ) a)).const_sub ((R : ℂ) ^ 2)
  have hd0 : (R : ℂ) ^ 2 - (starRingEnd ℂ) a * 0 ≠ 0 := by simpa using hR2
  have hb := hnum.div hden hd0
  have hb' : HasDerivAt (b R a) (((R : ℂ) * ((R : ℂ) ^ 2 - (starRingEnd ℂ) a * 0)
      - (R : ℂ) * (0 - a) * -(starRingEnd ℂ) a) / ((R : ℂ) ^ 2 - (starRingEnd ℂ) a * 0) ^ 2) 0 := hb
  rw [logDeriv_apply, hb'.deriv]
  unfold b
  simp only [mul_zero, sub_zero, zero_sub]
  field_simp
  ring

/-- Per-zero contribution sign: with a = ρ − c and Re a = β − σ ≤ 0, |a| ≤ R, the Blaschke-corrected term
Re(1/(c − ρ)) + Re(ρ − c)/R² = −Re a·(1/|a|² − 1/R²) is ≥ 0, so zeros may be dropped from −Re ζ'/ζ(c). -/
theorem blaschke_zero_term_nonneg {R : ℝ} {a : ℂ} (hR : 0 < R) (ha : ‖a‖ ≤ R) (ha0 : a ≠ 0) (hre : a.re ≤ 0) :
    0 ≤ (1 / (-a)).re + a.re / R ^ 2 := by
  have hn : 0 < Complex.normSq a := Complex.normSq_pos.mpr ha0
  have h1 : (1 / (-a)).re = -a.re / Complex.normSq a := by
    simp [Complex.div_re]
    ring
  rw [h1]
  have hsq : Complex.normSq a ≤ R ^ 2 := by
    rw [← Complex.sq_norm]; exact pow_le_pow_left₀ (norm_nonneg _) ha 2
  have h2 : -a.re / R ^ 2 ≤ -a.re / Complex.normSq a := div_le_div_of_nonneg_left (by linarith) hn hsq
  have h3 : a.re / R ^ 2 = -(-a.re / R ^ 2) := by ring
  linarith

end Blaschke

#print axioms Blaschke.blaschke_identity
#print axioms Blaschke.blaschke_den_ne
#print axioms Blaschke.blaschke_norm_le_one
#print axioms Blaschke.blaschke_norm_sphere
#print axioms Blaschke.blaschke_logDeriv_zero
#print axioms Blaschke.blaschke_zero_term_nonneg
