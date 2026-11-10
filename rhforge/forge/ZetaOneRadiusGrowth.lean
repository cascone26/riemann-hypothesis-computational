import ZetaBlaschkeLandau
import NormZetaGeExplicit
import AnalyticDiffContOnCl

/-!
# S3031 — one-radius growth & zero bounds for ζ, scalar form (lane 1, κ 2 → 1 program, brick 5a)

* `growth_one_radius` — no tracked zero: Re(−ζ'/ζ(c)) ≤ 2(log Mζ − log‖ζ c‖)/R
  (zero in ball: S3030 + nonneg bracket; zero-free ball: max modulus + sharp BC directly);
* `growth_one_radius_scalar` — Mζ = Ms + 1, ‖ζ c‖ ≥ (3 + 1/(Re c − 1))⁻¹:  ≤ 2·log((Ms+1)(3 + 1/(Re c − 1)))/R;
* `zero_one_radius_scalar` — with a zero ρ = β + iγ, c = σ + iγ:
  Re(−ζ'/ζ(c)) ≤ 2·log((Ms+1)(3+1/(σ−1)))/R + (σ − β)/R² − 1/(σ − β).
-/

open Complex Metric

namespace OneRadius

theorem growth_one_radius {c : ℂ} {R R' Ms Mζ : ℝ} (hR : 0 < R) (hRR' : R < R') (hcre : 1 < c.re)
    (hball : (1 : ℂ) ∉ ball c R') (hζsphere : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Ms) (hsl : Ms < Mζ) :
    (-logDeriv riemannZeta c).re ≤ 2 * (Real.log Mζ - Real.log ‖riemannZeta c‖) / R := by
  have hζc : riemannZeta c ≠ 0 := riemannZeta_ne_zero_of_one_lt_re hcre
  by_cases hz : ∃ ρ ∈ ball c R, riemannZeta ρ = 0
  · obtain ⟨ρ, hρ, hζρ⟩ := hz
    have key := zeta_landau_blaschke hR hRR' hcre hball hζsphere hsl hρ hζρ
    have hρre : ρ.re < 1 := by
      by_contra h; push_neg at h; exact riemannZeta_ne_zero_of_one_le_re h hζρ
    have hX := Blaschke.blaschke_zero_term_nonneg hR (by rw [← dist_eq_norm]; exact (mem_ball.mp hρ).le)
      (fun h => hζc (by rw [sub_eq_zero.mp h] at hζρ; exact hζρ)) (by simp; linarith)
    rw [show -(ρ - c) = c - ρ by ring] at hX
    linarith
  · push_neg at hz
    have hR' : 0 < R' := hR.trans hRR'
    have han : AnalyticOnNhd ℂ riemannZeta (ball c R') := fun z hz' =>
      analyticAt_riemannZeta (fun h => hball (h ▸ hz'))
    have hd : DiffContOnCl ℂ riemannZeta (ball c R) := analyticOnNhd_diffContOnCl hR hRR' han
    have hbd : ∀ z ∈ ball c R, ‖riemannZeta z‖ ≤ Ms := by
      intro z hzb
      refine Complex.norm_le_of_forall_mem_frontier_norm_le isBounded_ball hd ?_ (subset_closure hzb)
      rw [frontier_ball _ hR.ne']; exact hζsphere
    have hbc := SharpBC.re_logDeriv_ge_of_bounded_sharp_at hR hd.differentiableOn hz
      (fun z hzb => (hbd z hzb).trans hsl.le) (lt_of_le_of_lt (hbd c (mem_ball_self hR)) hsl)
    rw [Complex.neg_re]; linarith

theorem log_bound_scalar {c : ℂ} {Ms : ℝ} (hcre : 1 < c.re) (hMs : 0 ≤ Ms) :
    Real.log (Ms + 1) - Real.log ‖riemannZeta c‖ ≤ Real.log ((Ms + 1) * (3 + 1 / (c.re - 1))) := by
  have hlb := norm_riemannZeta_ge_explicit hcre
  have hD : 0 < 3 + 1 / (c.re - 1) := by
    have : 0 < 1 / (c.re - 1) := one_div_pos.mpr (by linarith)
    linarith
  have hζpos : 0 < ‖riemannZeta c‖ := norm_pos_iff.mpr (riemannZeta_ne_zero_of_one_lt_re hcre)
  rw [Real.log_mul (by linarith) hD.ne']
  have : -Real.log ‖riemannZeta c‖ ≤ Real.log (3 + 1 / (c.re - 1)) := by
    have h1 := Real.log_le_log (inv_pos.mpr hD) hlb
    rw [Real.log_inv] at h1; linarith
  linarith

theorem growth_one_radius_scalar {c : ℂ} {R R' Ms : ℝ} (hR : 0 < R) (hRR' : R < R') (hcre : 1 < c.re)
    (hball : (1 : ℂ) ∉ ball c R') (hζsphere : ∀ z ∈ sphere c R, ‖riemannZeta z‖ ≤ Ms) :
    (-logDeriv riemannZeta c).re ≤ 2 * Real.log ((Ms + 1) * (3 + 1 / (c.re - 1))) / R := by
  have hMs : 0 ≤ Ms := by
    have hz0 : c + (R : ℂ) ∈ sphere c R := by
      rw [mem_sphere, dist_eq_norm, add_sub_cancel_left, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hR]
    exact le_trans (norm_nonneg _) (hζsphere _ hz0)
  have h := growth_one_radius hR hRR' hcre hball hζsphere (lt_add_one Ms)
  have hl := log_bound_scalar hcre hMs
  have : 2 * (Real.log (Ms + 1) - Real.log ‖riemannZeta c‖) / R ≤ 2 * Real.log ((Ms + 1) * (3 + 1 / (c.re - 1))) / R :=
    div_le_div_of_nonneg_right (by linarith) hR.le
  linarith

theorem zero_one_radius_scalar {σ β γ R R' Ms : ℝ} (hR : 0 < R) (hRR' : R < R') (hσ : 1 < σ)
    (hball : (1 : ℂ) ∉ ball ((σ : ℂ) + I * (γ : ℝ)) R')
    (hζsphere : ∀ z ∈ sphere ((σ : ℂ) + I * (γ : ℝ)) R, ‖riemannZeta z‖ ≤ Ms)
    (hρ : ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) R)
    (hζρ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0) (hβσ : β < σ) :
    (-logDeriv riemannZeta ((σ : ℂ) + I * (γ : ℝ))).re
      ≤ 2 * Real.log ((Ms + 1) * (3 + 1 / (σ - 1))) / R + (σ - β) / R ^ 2 - 1 / (σ - β) := by
  have hcre : 1 < ((σ : ℂ) + I * (γ : ℝ)).re := by simp; exact hσ
  have hMs : 0 ≤ Ms := by
    have hz0 : ((σ : ℂ) + I * (γ : ℝ)) + (R : ℂ) ∈ sphere ((σ : ℂ) + I * (γ : ℝ)) R := by
      rw [mem_sphere, dist_eq_norm, add_sub_cancel_left, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hR]
    exact le_trans (norm_nonneg _) (hζsphere _ hz0)
  have key := zeta_landau_blaschke hR hRR' hcre hball hζsphere (lt_add_one Ms) hρ hζρ
  have hl := log_bound_scalar hcre hMs
  have hre1 : ((σ : ℂ) + I * (γ : ℝ)).re = σ := by simp
  rw [hre1] at hl
  have e1 : (1 / (((σ : ℂ) + I * (γ : ℝ)) - ((β : ℂ) + I * (γ : ℝ)))).re = 1 / (σ - β) := by
    rw [show ((σ : ℂ) + I * (γ : ℝ)) - ((β : ℂ) + I * (γ : ℝ)) = ((σ - β : ℝ) : ℂ) by push_cast; ring]
    rw [← Complex.ofReal_one, ← Complex.ofReal_div, Complex.ofReal_re]
  have e2 : (((β : ℂ) + I * (γ : ℝ)) - ((σ : ℂ) + I * (γ : ℝ))).re = β - σ := by simp
  rw [e1, e2] at key
  have h2 : 2 * (Real.log (Ms + 1) - Real.log ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ))‖) / R
      ≤ 2 * Real.log ((Ms + 1) * (3 + 1 / (σ - 1))) / R := div_le_div_of_nonneg_right (by linarith) hR.le
  have e3 : (β - σ) / R ^ 2 = -((σ - β) / R ^ 2) := by ring
  linarith

end OneRadius

#print axioms OneRadius.growth_one_radius
#print axioms OneRadius.growth_one_radius_scalar
#print axioms OneRadius.zero_one_radius_scalar
