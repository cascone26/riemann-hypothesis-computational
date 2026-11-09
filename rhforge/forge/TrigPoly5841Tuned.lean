import TrigPoly5841Uniform
import LogSharpClose40

/-!
# S2950 — TUNED unconditional 5-8-4-1 region: c₀ = 221/332640 ≈ 6.64e-4 (lane 1)

Re-runs the S2942→S2945→S2947 chain with the tighter S2949 closing estimate (40·log x instead of 64·log x) and the
near-optimal ε = 13/10 (maximizer of ε(3−ε)/(5+ε) is −5+√40 ≈ 1.3246; 13/10 also makes 520/ε = 400 ≤ 640):
  B = C₂ = C₃ = 40·log|γ|,  𝒦 = 5A₀ + 520·log|γ|,  ε(3−ε)/(5+ε) = 221/630.
* `zeta_region_5841_log40`, `riemannZeta_zero_free_5841_40` : β ≤ 1 − (221/630)/(5A₀ + 520·log|γ|) for |γ| ≥ e²⁵⁶ past
  the threshold 8·log(16(3+6A₀+640)) ≤ log|γ|;
* `riemannZeta_zero_free_5841_clean40` : ∃T₀, |γ| ≥ T₀ ⟹ β ≤ 1 − (221/332640)/log|γ|  (66% better than S2947's 9/22464,
  5.95× the forge 3-4-1 clean constant 1/8960).
-/

open Complex Metric

namespace TrigPoly5841


theorem zeta_region_5841_log40 :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ {γ β σ R R' Mζ Ms : ℝ},
      Real.exp 256 ≤ |γ| →
      8 * Real.log (16 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| →
      σ = 1 + (13 / 10) / (5 * A₀ + 520 * Real.log |γ|) →
      β < σ →
      0 < R → R < R' →
      (1 : ℂ) ∉ ball ((σ : ℂ) + I * (γ : ℝ)) R' →
      ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) R' →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall ((σ : ℂ) + I * (γ : ℝ)) (R / 2) →
      (∀ z ∈ sphere ((σ : ℂ) + I * (γ : ℝ)) R, ‖riemannZeta z‖ ≤ Ms) →
      Ms < Mζ → 0 < Mζ →
      4 * Real.log (Mζ / ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ 40 * Real.log |γ| →
      β ≤ 1 - (221 / 630) / (5 * A₀ + 520 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841_scalar
  refine ⟨A₀, hA0, ?_⟩
  intro γ β σ R R' Mζ Ms hxe hthr hσdef hβσ hR hRR' hball hρ₀ball hζρ₀ hp hζsphere hsl hMζ hBexp
  set L := Real.log |γ| with hLdef
  have hL16 : 256 ≤ L := by
    calc (256 : ℝ) = Real.log (Real.exp 256) := (Real.log_exp 256).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 256) hxe
  have hγ1 : 1 < |γ| := lt_of_lt_of_le (by
      have := Real.add_one_le_exp (256 : ℝ); linarith) hxe
  have hK : 5 * A₀ + 8 * (40 * L) + 4 * (40 * L) + 40 * L = 5 * A₀ + 520 * L := by ring
  have hKpos : 0 < 5 * A₀ + 520 * L := by nlinarith
  have hKε : (5 * A₀ + 520 * L) / (13 / 10) ≤ 6 * A₀ + 640 * L := by
    rw [div_le_iff₀ (by norm_num)]; nlinarith
  have hfac : 3 + (5 * A₀ + 520 * L) / (13 / 10) ≤ 3 + 6 * A₀ + 640 * L := by linarith
  have hfacpos : 0 < 3 + (5 * A₀ + 520 * L) / (13 / 10) := by positivity
  have hgrow : ∀ d : ℝ, 0 ≤ d →
      32 * Real.log ((d * |γ| + 10) * (3 + 6 * A₀ + 640 * L)) ≤ 40 * L →
      32 * Real.log ((d * |γ| + 10) * (3 + (5 * A₀ + 520 * L) / (13 / 10))) ≤ 40 * L := by
    intro d hd h
    refine le_trans ?_ h
    have hdx : 0 < d * |γ| + 10 := by positivity
    have := Real.log_le_log (by positivity) (mul_le_mul_of_nonneg_left hfac hdx.le)
    linarith
  have e2 : 3 * |2 * γ| = 6 * |γ| := by rw [abs_mul]; norm_num; ring
  have e3 : 3 * |3 * γ| = 9 * |γ| := by rw [abs_mul]; norm_num; ring
  have g6 := hgrow 6 (by norm_num) (log_sharp_close_40 hA0 (by norm_num) (by norm_num) hxe hthr)
  have g9 := hgrow 9 (by norm_num) (log_sharp_close_40 hA0 (by norm_num) (by norm_num) hxe hthr)
  have hgap := H (γ := γ) (β := β) (ε := 13 / 10) (B := 40 * L) (C₂ := 40 * L) (C₃ := 40 * L)
    (σ := σ) (R := R) (R' := R') (Mζ := Mζ) (Ms := Ms)
    (by norm_num) (by positivity) (by positivity) (by positivity)
    (by rw [hK]; nlinarith) (by rw [hK]; exact hσdef) hβσ hR hRR' hball hρ₀ball hζρ₀ hp hζsphere
    hsl hMζ hBexp hγ1
    (by rw [hK, e2]; exact g6) (by rw [hK, e3]; exact g9)
  rw [hK] at hgap
  have e : (13 / 10 : ℝ) * (3 - 13 / 10) / ((5 + 13 / 10) * (5 * A₀ + 520 * L))
      = (221 / 630) / (5 * A₀ + 520 * L) := by
    field_simp; ring
  rw [e] at hgap
  exact hgap

theorem riemannZeta_zero_free_5841_40 :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ,
      Real.exp 256 ≤ |γ| →
      8 * Real.log (16 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - (221 / 630) / (5 * A₀ + 520 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841_log40
  refine ⟨A₀, hA0, ?_⟩
  intro γ β hxe hthr hζρ₀
  set L := Real.log |γ| with hLdef
  have hL16 : 256 ≤ L := by
    calc (256 : ℝ) = Real.log (Real.exp 256) := (Real.log_exp 256).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 256) hxe
  have hγ1 : 1 < |γ| := lt_of_lt_of_le (by
      have := Real.add_one_le_exp (256 : ℝ); linarith) hxe
  set 𝒦 := 5 * A₀ + 520 * L with h𝒦
  have hKpos : 0 < 𝒦 := by rw [h𝒦]; nlinarith
  set σ := 1 + (13 / 10) / 𝒦 with hσdef
  have hq : 0 < (13 / 10) / 𝒦 := by positivity
  have hq1 : (13 / 10) / 𝒦 ≤ 1 := by rw [div_le_one hKpos]; rw [h𝒦]; nlinarith
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
  · have h1 : (13 / 10) / 𝒦 + (221 / 630) / 𝒦 ≤ 1 / 8 := by
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
  -- B = 40 L bound
  have hinv : 1 / (σ - 1) = 𝒦 / (13 / 10) := by
    rw [hσdef, show 1 + 13 / 10 / 𝒦 - 1 = 13 / 10 / 𝒦 by ring, one_div_div]
  have hKε : 𝒦 / (13 / 10) ≤ 6 * A₀ + 640 * L := by
    rw [div_le_iff₀ (by norm_num), h𝒦]; nlinarith
  have hclean : cleanRHS σ γ (1 / 4) ≤ 40 * L := by
    refine le_trans (cleanRHS_quarter_le hσ1 hσ2 hγ1) ?_
    rw [hinv]
    refine le_trans ?_ (log_sharp_close_40 (a := A₀) (x := |γ|) (d := 3) hA0 (by norm_num) (by norm_num) hxe hthr)
    have hdx : (0 : ℝ) < 3 * |γ| + 10 := by positivity
    have hf : 3 + 𝒦 / (13 / 10) ≤ 3 + 6 * A₀ + 640 * L := by linarith
    have hfp : 0 < 3 + 𝒦 / (13 / 10) := by positivity
    have := Real.log_le_log (by positivity) (mul_le_mul_of_nonneg_left hf hdx.le)
    linarith
  have hgr : growthRHS σ γ (1 / 4) ≤ 40 * L :=
    growthRHS_le_clean hσ1 hR (by norm_num) (by linarith) hclean
  have hBexp : 4 * Real.log ((Ms + 1) / ‖riemannZeta c‖) / (1 / 4 / 2) ≤ 40 * L :=
    hBexplicit_of_ge hcre hMspos hR hgr
  exact H hxe hthr hσdef hβσ hR (by norm_num) hball hρball hζρ₀ hp
    (fun z hz => norm_riemannZeta_le_sphere_N1 hR hcR hc1 hz) (by linarith) hMspos hBexp

theorem riemannZeta_zero_free_5841_clean40 :
    ∃ T₀ : ℝ, ∀ γ β : ℝ, T₀ ≤ |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - (221 / 332640) / Real.log |γ| := by
  obtain ⟨A₀, hA₀nn, hmain⟩ := riemannZeta_zero_free_5841_40
  refine ⟨max (max (Real.exp 256) (Real.exp (8 * Real.log (16 * (3 + 6 * A₀ + 640)))))
      (Real.exp (5 * A₀ / 8)), ?_⟩
  intro γ β hγ hζ
  have h1 : Real.exp 256 ≤ |γ| := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hγ
  have hexp2 : Real.exp (8 * Real.log (16 * (3 + 6 * A₀ + 640))) ≤ |γ| :=
    le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hγ
  have hexp3 : Real.exp (5 * A₀ / 8) ≤ |γ| := le_trans (le_max_right _ _) hγ
  have h2 : 8 * Real.log (16 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| := by
    calc 8 * Real.log (16 * (3 + 6 * A₀ + 640))
        = Real.log (Real.exp (8 * Real.log (16 * (3 + 6 * A₀ + 640)))) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2
  have hL16 : 256 ≤ Real.log |γ| := by
    calc (256 : ℝ) = Real.log (Real.exp 256) := (Real.log_exp 256).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 256) h1
  have hLpos : 0 < Real.log |γ| := by linarith
  have hL5A : 5 * A₀ / 8 ≤ Real.log |γ| := by
    calc 5 * A₀ / 8 = Real.log (Real.exp (5 * A₀ / 8)) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp3
  have hmain' := hmain γ β h1 h2 hζ
  have hKpos : 0 < 5 * A₀ + 520 * Real.log |γ| := by nlinarith
  have hcmp : 5 * A₀ + 520 * Real.log |γ| ≤ 528 * Real.log |γ| := by nlinarith
  have hdiv : (221 / 630) / (528 * Real.log |γ|) ≤ (221 / 630) / (5 * A₀ + 520 * Real.log |γ|) :=
    div_le_div_of_nonneg_left (by norm_num) hKpos hcmp
  have he : (221 / 332640 : ℝ) / Real.log |γ| = (221 / 630) / (528 * Real.log |γ|) := by
    field_simp; ring
  rw [he]; linarith

end TrigPoly5841

#print axioms TrigPoly5841.zeta_region_5841_log40
#print axioms TrigPoly5841.riemannZeta_zero_free_5841_40
#print axioms TrigPoly5841.riemannZeta_zero_free_5841_clean40
