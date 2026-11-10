import InvGammaStrip
import CompletedZeta0StripReflect
import GammaRNorm
import ZetaEntireE

/-!
# S2997 — growth of E(s) = (s−1)ζ(s) on the closed critical strip (lane 1: hB for three-lines)

* `strip_term_le` — for 1/2 ≤ Re u ≤ 1, |Im u| ≥ 1: the forge Λ₀ strip RHS at u is ≤ 14 + 12|Im u|.
* `norm_lam0_le` — 0 ≤ σ ≤ 1, |t| ≥ 1 ⇒ ‖Λ₀(s)‖ ≤ 14 + 12|t|  (direct for σ ≥ 1/2, reflected for σ < 1/2).
* `norm_zeta_strip_le` — ‖ζ(s)‖ ≤ (16 + 12|t|)·2·e^{π|t|/2}   (ζ = Λ/Γℝ, |1/Γ(s/2)| ≤ (2/π)e^{π|t|/2} by S2995).
* `norm_zetaE_strip_le` — ‖E(s)‖ ≤ 2(1 + |t|)(16 + 12|t|)e^{π|t|/2}.
Then |E·G_T| ≤ poly(t)·e^{π|t|/2 + σ² − (t−T)²} is bounded on the strip ⇒ S2984's hB (next brick: the BddAbove packaging).
-/

open Complex

namespace ZetaEStrip

theorem strip_term_le (u : ℂ) (hu0 : 1 / 2 ≤ u.re) (hu1 : u.re ≤ 1) (hut : 1 ≤ |u.im|) :
    (Real.pi ^ (-u.re / 2) * Real.Gamma (u.re / 2)) * (‖u‖ / ‖u - 1‖ + ‖u‖ / u.re)
      + 1 / ‖u‖ + 1 / ‖1 - u‖ ≤ 14 + 12 * |u.im| := by
  have hpi1 : (1 : ℝ) ≤ Real.pi := by linarith [Real.pi_gt_three]
  have hp : Real.pi ^ (-u.re / 2) ≤ 1 := Real.rpow_le_one_of_one_le_of_nonpos hpi1 (by linarith)
  have hp0 : 0 ≤ Real.pi ^ (-u.re / 2) := Real.rpow_nonneg Real.pi_pos.le _
  have hG : Real.Gamma (u.re / 2) ≤ 4 := by
    refine le_trans (InvGammaStrip.Gamma_le_inv (by linarith) (by linarith)) ?_
    rw [div_le_iff₀ (by linarith)]; linarith
  have hG0 : 0 ≤ Real.Gamma (u.re / 2) := (Real.Gamma_pos_of_pos (by linarith)).le
  have hun : ‖u‖ ≤ 1 + |u.im| := by
    have := Complex.norm_le_abs_re_add_abs_im u
    rw [abs_of_pos (by linarith : (0 : ℝ) < u.re)] at this; linarith
  have hu1n : 1 ≤ ‖u‖ := le_trans hut (Complex.abs_im_le_norm u)
  have hum1 : 1 ≤ ‖u - 1‖ := by
    have := Complex.abs_im_le_norm (u - 1); simp at this; linarith
  have h1mu : 1 ≤ ‖1 - u‖ := by rw [norm_sub_rev]; exact hum1
  have hA : ‖u‖ / ‖u - 1‖ ≤ 1 + |u.im| := by
    rw [div_le_iff₀ (by linarith)]; nlinarith [norm_nonneg u]
  have hB : ‖u‖ / u.re ≤ 2 * (1 + |u.im|) := by
    rw [div_le_iff₀ (by linarith)]; nlinarith [norm_nonneg u]
  have hC : 1 / ‖u‖ ≤ 1 := by rw [div_le_one (by linarith)]; exact hu1n
  have hD : 1 / ‖1 - u‖ ≤ 1 := by rw [div_le_one (by linarith)]; exact h1mu
  have hPG : Real.pi ^ (-u.re / 2) * Real.Gamma (u.re / 2) ≤ 4 := by nlinarith
  have hPG0 : 0 ≤ Real.pi ^ (-u.re / 2) * Real.Gamma (u.re / 2) := mul_nonneg hp0 hG0
  have hsum0 : 0 ≤ ‖u‖ / ‖u - 1‖ + ‖u‖ / u.re := by positivity
  have : (Real.pi ^ (-u.re / 2) * Real.Gamma (u.re / 2)) * (‖u‖ / ‖u - 1‖ + ‖u‖ / u.re)
      ≤ 4 * (3 * (1 + |u.im|)) := mul_le_mul hPG (by linarith) hsum0 (by norm_num)
  linarith

theorem norm_lam0_le {s : ℂ} (h0 : 0 ≤ s.re) (h1 : s.re ≤ 1) (ht : 1 ≤ |s.im|) :
    ‖completedRiemannZeta₀ s‖ ≤ 14 + 12 * |s.im| := by
  have hs0 : s ≠ 0 := by intro h; rw [h] at ht; simp at ht; linarith
  have hs1 : s ≠ 1 := by intro h; rw [h] at ht; simp at ht; linarith
  rcases le_or_gt (1 / 2) s.re with hh | hl
  · exact le_trans (norm_completedRiemannZeta0_le_strip (by linarith) hs1)
      (strip_term_le s hh h1 ht)
  · have h := norm_completedRiemannZeta0_le_strip_reflect (by linarith) hs0
    have hre : (1 - s).re = 1 - s.re := by simp
    have him : (1 - s).im = -s.im := by simp
    have := strip_term_le (1 - s) (by rw [hre]; linarith) (by rw [hre]; linarith)
      (by rw [him, abs_neg]; exact ht)
    rw [him, abs_neg] at this
    exact le_trans h this

theorem norm_zeta_strip_le {s : ℂ} (h0 : 0 ≤ s.re) (h1 : s.re ≤ 1) (ht : 1 ≤ |s.im|) :
    ‖riemannZeta s‖ ≤ (16 + 12 * |s.im|) * (2 * Real.exp (Real.pi * |s.im| / 2)) := by
  have hs0 : s ≠ 0 := by intro h; rw [h] at ht; simp at ht; linarith
  have hu1n : 1 ≤ ‖s‖ := le_trans ht (Complex.abs_im_le_norm s)
  have h1s : 1 ≤ ‖1 - s‖ := by
    have := Complex.abs_im_le_norm (1 - s); simp at this; linarith
  have hLam : ‖completedRiemannZeta s‖ ≤ 16 + 12 * |s.im| := by
    rw [completedRiemannZeta_eq]
    have hl := norm_lam0_le h0 h1 ht
    have e1 : ‖(1 : ℂ) / s‖ ≤ 1 := by rw [norm_div, norm_one, div_le_one (by linarith)]; exact hu1n
    have e2 : ‖(1 : ℂ) / (1 - s)‖ ≤ 1 := by
      rw [norm_div, norm_one, div_le_one (by linarith)]; exact h1s
    calc ‖completedRiemannZeta₀ s - 1 / s - 1 / (1 - s)‖
        ≤ ‖completedRiemannZeta₀ s - 1 / s‖ + ‖(1 : ℂ) / (1 - s)‖ := norm_sub_le _ _
      _ ≤ ‖completedRiemannZeta₀ s‖ + ‖(1 : ℂ) / s‖ + ‖(1 : ℂ) / (1 - s)‖ := by
          gcongr; exact norm_sub_le _ _
      _ ≤ 16 + 12 * |s.im| := by linarith
  rw [riemannZeta_def_of_ne_zero hs0, div_eq_mul_inv, norm_mul, norm_inv, norm_Gammaℝ]
  have hz : (s / 2).re = s.re / 2 := by simp
  have hzi : (s / 2).im = s.im / 2 := by simp
  have hinvG := InvGammaStrip.norm_inv_Gamma_le (z := s / 2) (by rw [hz]; linarith)
    (by rw [hz]; linarith)
  rw [hzi, norm_inv] at hinvG
  have hpow : (Real.pi ^ (-s.re / 2) * ‖Complex.Gamma (s / 2)‖)⁻¹
      = Real.pi ^ (s.re / 2) * ‖Complex.Gamma (s / 2)‖⁻¹ := by
    rw [mul_inv, ← Real.rpow_neg Real.pi_pos.le]; ring_nf
  rw [hpow]
  have hpi2 : Real.pi ^ (s.re / 2) ≤ 2 := by
    have : Real.pi ^ (s.re / 2) ≤ Real.pi ^ (1 / 2 : ℝ) :=
      Real.rpow_le_rpow_of_exponent_le (by linarith [Real.pi_gt_three]) (by linarith)
    refine le_trans this ?_
    rw [← Real.sqrt_eq_rpow, Real.sqrt_le_left (by norm_num)]
    linarith [Real.pi_lt_four]
  have hpi0 : 0 ≤ Real.pi ^ (s.re / 2) := Real.rpow_nonneg Real.pi_pos.le _
  have habs : |s.im / 2| = |s.im| / 2 := by rw [abs_div, abs_two]
  rw [habs] at hinvG
  have h2pi : 2 / Real.pi ≤ 1 := by rw [div_le_one Real.pi_pos]; linarith [Real.pi_gt_three]
  have hexp0 : 0 ≤ Real.exp (Real.pi * (|s.im| / 2)) := (Real.exp_pos _).le
  have hinv' : ‖Complex.Gamma (s / 2)‖⁻¹ ≤ Real.exp (Real.pi * |s.im| / 2) := by
    rw [show Real.pi * |s.im| / 2 = Real.pi * (|s.im| / 2) by ring]
    nlinarith
  have hL0 : 0 ≤ ‖completedRiemannZeta s‖ := norm_nonneg _
  calc ‖completedRiemannZeta s‖ * (Real.pi ^ (s.re / 2) * ‖Complex.Gamma (s / 2)‖⁻¹)
      ≤ (16 + 12 * |s.im|) * (2 * Real.exp (Real.pi * |s.im| / 2)) := by
        apply mul_le_mul hLam _ (by positivity) (by positivity)
        exact mul_le_mul hpi2 hinv' (by positivity) (by norm_num)

theorem norm_zetaE_strip_le {s : ℂ} (h0 : 0 ≤ s.re) (h1 : s.re ≤ 1) (ht : 1 ≤ |s.im|) :
    ‖ZetaEntire.zetaE s‖
      ≤ (1 + |s.im|) * ((16 + 12 * |s.im|) * (2 * Real.exp (Real.pi * |s.im| / 2))) := by
  have hs1 : s ≠ 1 := by intro h; rw [h] at ht; simp at ht; linarith
  rw [ZetaEntire.zetaE_eq hs1, norm_mul]
  have hsm1 : ‖s - 1‖ ≤ 1 + |s.im| := by
    have := Complex.norm_le_abs_re_add_abs_im (s - 1)
    simp only [Complex.sub_re, Complex.one_re, Complex.sub_im, Complex.one_im, sub_zero] at this
    have : |s.re - 1| ≤ 1 := by rw [abs_le]; constructor <;> linarith
    linarith
  exact mul_le_mul hsm1 (norm_zeta_strip_le h0 h1 ht) (norm_nonneg _) (by positivity)

end ZetaEStrip

#print axioms ZetaEStrip.strip_term_le
#print axioms ZetaEStrip.norm_lam0_le
#print axioms ZetaEStrip.norm_zeta_strip_le
#print axioms ZetaEStrip.norm_zetaE_strip_le
