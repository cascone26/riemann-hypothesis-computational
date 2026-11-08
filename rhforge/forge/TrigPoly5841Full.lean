import TrigPoly5841Log
import HBexplicitOfGe
import ZetaSphereN1
import PoleAvoidGeom
import AbsLeNormCenterSubOne

/-!
# S2945 — the UNCONDITIONAL 5-8-4-1 zero-free region (lane 1)

Discharges the last hypothesis of S2942 (the zero-input sphere data at σ+iγ) exactly as §1448 did for 3-4-1:
R = 1/4, R' = 1/2, Ms := 2 + 1/(‖c−1‖−R) + (‖c‖+R)/(c.re−R) (§1447 N=1 sphere bound), Mζ := Ms+1,
`hBexplicit_of_ge` (§1444) + §1452 + S2940 `cleanRHS_quarter_le` + S2930 `log_sharp_close` (d = 3) give
B = 64·log|γ|.  A zero with β < σ − R/2 is already left of the claimed line (𝒦 ≥ 832·16), and β ≥ σ > 1 is
impossible.  Result:

  ∃ A₀ ≥ 0, ∀ γ β, |γ| ≥ e¹⁶ ∧ 2·log(16(3+6A₀+640)) ≤ log|γ| ∧ ζ(β+iγ) = 0  ⟹
      β ≤ 1 − (9/26)/(5A₀ + 832·log|γ|).

No analytic hypothesis remains: pole (§1572), zero (S2934), growth (S2928) and all constants are discharged.
-/

open Complex Metric

namespace TrigPoly5841

theorem norm_sub_center {β σ γ : ℝ} :
    ‖((β : ℂ) + I * (γ : ℝ)) - ((σ : ℂ) + I * (γ : ℝ))‖ = |β - σ| := by
  rw [show ((β : ℂ) + I * (γ : ℝ)) - ((σ : ℂ) + I * (γ : ℝ)) = ((β - σ : ℝ) : ℂ) by push_cast; ring]
  exact Complex.norm_real _

theorem riemannZeta_zero_free_5841 :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ,
      Real.exp 16 ≤ |γ| →
      2 * Real.log (16 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - (9 / 26) / (5 * A₀ + 832 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841_log
  refine ⟨A₀, hA0, ?_⟩
  intro γ β hxe hthr hζρ₀
  set L := Real.log |γ| with hLdef
  have hL16 : 16 ≤ L := by
    calc (16 : ℝ) = Real.log (Real.exp 16) := (Real.log_exp 16).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 16) hxe
  have hγ1 : 1 < |γ| := lt_of_lt_of_le (by
      have := Real.add_one_le_exp (16 : ℝ); linarith) hxe
  set 𝒦 := 5 * A₀ + 832 * L with h𝒦
  have hKpos : 0 < 𝒦 := by rw [h𝒦]; nlinarith
  set σ := 1 + (3 / 2) / 𝒦 with hσdef
  have hq : 0 < (3 / 2) / 𝒦 := by positivity
  have hq1 : (3 / 2) / 𝒦 ≤ 1 := by rw [div_le_one hKpos]; rw [h𝒦]; nlinarith
  have hσ1 : 1 < σ := by rw [hσdef]; linarith
  have hσ2 : σ ≤ 2 := by rw [hσdef]; linarith
  -- β < σ, else ζ ≠ 0
  have hβσ : β < σ := by
    by_contra hcon
    push Not at hcon
    have hre : 1 < ((β : ℂ) + I * (γ : ℝ)).re := by simp; linarith
    exact riemannZeta_ne_zero_of_one_lt_re hre hζρ₀
  -- far case
  by_cases hfar : β < σ - 1 / 8
  · have h1 : (3 / 2) / 𝒦 + (9 / 26) / 𝒦 ≤ 1 / 8 := by
      rw [← add_div, div_le_iff₀ hKpos]; rw [h𝒦]; nlinarith
    linarith
  push Not at hfar
  -- close case: sphere data at c = σ + iγ with R = 1/4, R' = 1/2
  set c : ℂ := (σ : ℂ) + I * (γ : ℝ) with hc_def
  have hcre_eq : c.re = σ := by rw [hc_def]; simp
  have hcre : 1 < c.re := by rw [hcre_eq]; exact hσ1
  have hR : (0 : ℝ) < 1 / 4 := by norm_num
  have hcR : 0 < c.re - 1 / 4 := by rw [hcre_eq]; linarith
  have hc1 : (1 / 4 : ℝ) < ‖c - 1‖ :=
    lt_of_lt_of_le (show (1 / 4 : ℝ) < |γ| by linarith) (by rw [hc_def]; exact abs_le_norm_center_sub_one)
  set Ms : ℝ := 2 + 1 / (‖c - 1‖ - 1 / 4) + (‖c‖ + 1 / 4) / (c.re - 1 / 4) with hMs_def
  have hden1 : (0 : ℝ) < ‖c - 1‖ - 1 / 4 := by linarith
  have hMspos : 0 < Ms + 1 := by rw [hMs_def]; positivity
  have hdist : |β - σ| ≤ 1 / 4 / 2 := by rw [abs_sub_comm, abs_of_pos (by linarith)]; linarith
  have hp : ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall c (1 / 4 / 2) := by
    rw [mem_closedBall, dist_eq_norm, hc_def, norm_sub_center]; exact hdist
  have hρball : ((β : ℂ) + I * (γ : ℝ)) ∈ ball c (1 / 2) := by
    rw [mem_ball, dist_eq_norm, hc_def, norm_sub_center]; linarith
  have hball : (1 : ℂ) ∉ ball c (1 / 2) := by
    rw [hc_def]; exact one_notMem_ball_test (by linarith)
  -- B = 64 L bound
  have hinv : 1 / (σ - 1) = 𝒦 / (3 / 2) := by
    rw [hσdef, show 1 + 3 / 2 / 𝒦 - 1 = 3 / 2 / 𝒦 by ring, one_div_div]
  have hKε : 𝒦 / (3 / 2) ≤ 6 * A₀ + 640 * L := by
    rw [div_le_iff₀ (by norm_num), h𝒦]; nlinarith
  have hclean : cleanRHS σ γ (1 / 4) ≤ 64 * L := by
    refine le_trans (cleanRHS_quarter_le hσ1 hσ2 hγ1) ?_
    rw [hinv]
    refine le_trans ?_ (log_sharp_close (a := A₀) (x := |γ|) (d := 3) hA0 (by norm_num) (by norm_num) hxe hthr)
    have hdx : (0 : ℝ) < 3 * |γ| + 10 := by positivity
    have hf : 3 + 𝒦 / (3 / 2) ≤ 3 + 6 * A₀ + 640 * L := by linarith
    have hfp : 0 < 3 + 𝒦 / (3 / 2) := by positivity
    have := Real.log_le_log (by positivity) (mul_le_mul_of_nonneg_left hf hdx.le)
    linarith
  have hgr : growthRHS σ γ (1 / 4) ≤ 64 * L :=
    growthRHS_le_clean hσ1 hR (by norm_num) (by linarith) hclean
  have hBexp : 4 * Real.log ((Ms + 1) / ‖riemannZeta c‖) / (1 / 4 / 2) ≤ 64 * L :=
    hBexplicit_of_ge hcre hMspos hR hgr
  exact H hxe hthr hσdef hβσ hR (by norm_num) hball hρball hζρ₀ hp
    (fun z hz => norm_riemannZeta_le_sphere_N1 hR hcR hc1 hz) (by linarith) hMspos hBexp

end TrigPoly5841

#print axioms TrigPoly5841.riemannZeta_zero_free_5841
