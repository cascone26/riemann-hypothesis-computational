import ZetaConvexityZeta
import ZetaReGeOne
import NormZetaGeExplicit

/-!
# S3044 — pointwise log|ζ| bounds on the Poisson circle (lane 1, Poisson–Landau P5b-1)

c = σ + iγ, 1 < σ ≤ 3/2, |γ| ≥ 3, ‖w‖ ≤ 1/2, L := log(2+|γ|):
* `log_zeta_upper` — ∃ A₁ ≥ 0: log‖ζ(c+w)‖ ≤ A₁ + log(L+6) + (L/2)·max(−Re w, 0)
  (left of Re = 1: convexity S3002, exponent (1−σ')/2 ≤ max(−Re w,0)/2; right: ‖ζ‖ ≤ log|t|+6, S3006);
* `log_zeta_lower` — Re w ≥ 0 ⇒ −log(3 + 1/(σ−1)) ≤ log‖ζ(c+w)‖  (NormZetaGeExplicit).
These are exactly `hup`/`hlo` of PoissonAvg.neg_avg_le (S3041) with A = A₁ + log(L+6), D = log(3+1/(σ−1)).
-/

open Complex

namespace ZetaCircleLog

theorem log_zeta_upper : ∃ A₁ : ℝ, 0 ≤ A₁ ∧ ∀ (σ γ : ℝ) (w : ℂ), 1 < σ → σ ≤ 3 / 2 → 3 ≤ |γ| → ‖w‖ ≤ 1 / 2 →
    Real.log ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ) + w)‖
      ≤ A₁ + Real.log (Real.log (2 + |γ|) + 6) + Real.log (2 + |γ|) / 2 * max (-w.re) 0 := by
  obtain ⟨C, hC0, hC⟩ := ZetaConvexityZeta.zeta_convexity
  refine ⟨Real.log (C + 1), Real.log_nonneg (by linarith), ?_⟩
  intro σ γ w hσ hσ2 hγ hw
  set L := Real.log (2 + |γ|) with hL
  have hγ2 : (1 : ℝ) ≤ 2 + |γ| := by linarith [abs_nonneg γ]
  have hL0 : 0 ≤ L := Real.log_nonneg hγ2
  have hLL : 0 ≤ Real.log (L + 6) := Real.log_nonneg (by linarith)
  have hA0 : 0 ≤ Real.log (C + 1) := Real.log_nonneg (by linarith)
  have hwre : |w.re| ≤ 1 / 2 := (Complex.abs_re_le_norm w).trans hw
  have hwim : |w.im| ≤ 1 / 2 := (Complex.abs_im_le_norm w).trans hw
  set σ' := σ + w.re with hσ'
  set T := γ + w.im with hT
  have hz : (σ : ℂ) + I * (γ : ℝ) + w = (σ' : ℂ) + (T : ℂ) * I := by
    apply Complex.ext <;> simp [hσ', hT] <;> ring
  rw [hz]
  have hTabs1 : |T| ≤ |γ| + 1 / 2 := by
    rw [hT]; calc |γ + w.im| ≤ |γ| + |w.im| := abs_add_le _ _
      _ ≤ |γ| + 1 / 2 := by linarith
  have hTabs2 : |γ| - 1 / 2 ≤ |T| := by
    have : |γ| ≤ |T| + |w.im| := by
      calc |γ| = |T - w.im| := by rw [hT]; ring_nf
        _ ≤ |T| + |w.im| := abs_sub _ _
    linarith
  have hmax0 : 0 ≤ max (-w.re) 0 := le_max_right _ _
  set Z := riemannZeta ((σ' : ℂ) + (T : ℂ) * I) with hZ
  rcases (norm_nonneg Z).eq_or_lt with h0 | hpos
  · rw [← h0, Real.log_zero]; positivity
  have hRHS0 : 0 ≤ L / 2 * max (-w.re) 0 := by positivity
  by_cases hle : σ' ≤ 1
  · have hb := hC T σ' (by linarith) (by rw [hσ']; linarith [abs_le.mp hwre]) hle
    have hX1 : (1 : ℝ) ≤ 1 + |T| := by linarith [abs_nonneg T]
    have hXle : 1 + |T| ≤ 2 + |γ| := by linarith
    have hexp : (1 - σ') / 2 ≤ max (-w.re) 0 / 2 := by
      have : 1 - σ' ≤ -w.re := by rw [hσ']; linarith
      have := this.trans (le_max_left (-w.re) 0)
      linarith
    have hexp0 : 0 ≤ (1 - σ') / 2 := by linarith
    have hpow : (1 + |T|) ^ ((1 - σ') / 2) ≤ (2 + |γ|) ^ (max (-w.re) 0 / 2) :=
      (Real.rpow_le_rpow (by linarith) hXle hexp0).trans (Real.rpow_le_rpow_of_exponent_le hγ2 hexp)
    have hlog6 : Real.log (1 + |T|) + 6 ≤ L + 6 := by
      have := Real.log_le_log (by linarith) hXle; rw [hL]; linarith
    have hl6pos : 0 < Real.log (1 + |T|) + 6 := by
      have := Real.log_nonneg hX1; linarith
    have hbound : ‖Z‖ ≤ (C + 1) * (2 + |γ|) ^ (max (-w.re) 0 / 2) * (L + 6) := by
      calc ‖Z‖ ≤ C * (1 + |T|) ^ ((1 - σ') / 2) * (Real.log (1 + |T|) + 6) := hb
        _ ≤ (C + 1) * (2 + |γ|) ^ (max (-w.re) 0 / 2) * (L + 6) := by
          apply mul_le_mul (mul_le_mul (by linarith) hpow (by positivity) (by linarith)) hlog6 hl6pos.le
          positivity
    have hlog := Real.log_le_log hpos hbound
    rw [Real.log_mul (by positivity) (by linarith), Real.log_mul (by positivity) (by positivity),
      Real.log_rpow (by linarith)] at hlog
    rw [← hL] at hlog
    linarith
  · push Not at hle
    have hz2 : ((σ' : ℂ) + (T : ℂ) * I).re = σ' := by simp
    have hzi : ((σ' : ℂ) + (T : ℂ) * I).im = T := by simp
    have hb := ZetaReGeOne.norm_zeta_le_log_right (z := (σ' : ℂ) + (T : ℂ) * I)
      (by rw [hz2]; exact hle) (by rw [hz2, hσ']; linarith [abs_le.mp hwre]) (by rw [hzi]; linarith)
    rw [hzi] at hb
    have hTpos : 0 < |T| := by linarith
    have hlogT : Real.log |T| ≤ L := Real.log_le_log hTpos (by linarith)
    have hlog := Real.log_le_log hpos (hb.trans (by linarith : Real.log |T| + 6 ≤ L + 6))
    linarith

theorem log_zeta_lower {σ γ : ℝ} {w : ℂ} (hσ : 1 < σ) (hw : 0 ≤ w.re) :
    -Real.log (3 + 1 / (σ - 1)) ≤ Real.log ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ) + w)‖ := by
  have hre : ((σ : ℂ) + I * (γ : ℝ) + w).re = σ + w.re := by simp
  have hcre : 1 < ((σ : ℂ) + I * (γ : ℝ) + w).re := by rw [hre]; linarith
  have hlb := norm_riemannZeta_ge_explicit hcre
  rw [hre] at hlb
  have hD : 0 < 3 + 1 / (σ - 1) := by have : 0 < 1 / (σ - 1) := one_div_pos.mpr (by linarith); linarith
  have hD' : 0 < 3 + 1 / (σ + w.re - 1) := by
    have : 0 < 1 / (σ + w.re - 1) := one_div_pos.mpr (by linarith); linarith
  have hmono : 3 + 1 / (σ + w.re - 1) ≤ 3 + 1 / (σ - 1) := by
    have : 1 / (σ + w.re - 1) ≤ 1 / (σ - 1) := one_div_le_one_div_of_le (by linarith) (by linarith)
    linarith
  have h1 : (3 + 1 / (σ - 1))⁻¹ ≤ ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ) + w)‖ :=
    (inv_anti₀ hD' hmono).trans hlb
  have := Real.log_le_log (inv_pos.mpr hD) h1
  rw [Real.log_inv] at this
  exact this

end ZetaCircleLog

#print axioms ZetaCircleLog.log_zeta_upper
#print axioms ZetaCircleLog.log_zeta_lower
