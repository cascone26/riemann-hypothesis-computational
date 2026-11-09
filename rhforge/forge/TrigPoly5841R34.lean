import TrigPoly5841Tuned

/-!
# S2952 — 5-8-4-1 region at Jensen radius R = 3/4: c₀ = 221/110880 ≈ 1.99e-3 (lane 1)

The Borel–Carathéodory/Jensen step costs 4/(R/2) = 8/R per logarithm; the chain so far used R = 1/4 (cost 32).
At R = 3/4 the cost is 32/3, while the sphere bound only worsens inside the log:
  `cleanRHS_34_le` : 1<σ≤2, |t| ≥ 5 ⟹ cleanRHS σ t (3/4) ≤ (32/3)·log((5|t|+10)(3+1/(σ−1))).
Heights γ, 2γ, 3γ give d = 5, 10, 15 — all inside S2949 `log_sharp_close_40` (d ≤ 15), so
B = C₂ = C₃ = (40/3)·log|γ|, 𝒦 = 5A₀ + (520/3)·log|γ|, ε = 13/10:
* `riemannZeta_zero_free_5841_R34` : β ≤ 1 − (221/630)/(5A₀ + (520/3)·log|γ|);
* `riemannZeta_zero_free_5841_clean_R34` : β ≤ 1 − (221/110880)/log|γ| for |γ| ≥ T₀  (×3 over S2950).
-/

open Complex Metric

namespace TrigPoly5841

theorem cleanRHS_34_le {σ t : ℝ} (hσ1 : 1 < σ) (hσ2 : σ ≤ 2) (ht : 5 ≤ |t|) :
    cleanRHS σ t (3 / 4) ≤ 32 / 3 * Real.log ((5 * |t| + 10) * (3 + 1 / (σ - 1))) := by
  unfold cleanRHS
  have hden_pos : 0 < σ - 3 / 4 := by linarith
  have htpos : 0 < |t| - 3 / 4 := by linarith
  have ha : 1 / (|t| - 3 / 4) ≤ 1 := by rw [div_le_iff₀ htpos]; linarith
  have hb : (σ + |t| + 3 / 4) / (σ - 3 / 4) ≤ 4 * |t| + 11 := by
    rw [div_le_iff₀ hden_pos]
    nlinarith [hσ2, hden_pos, abs_nonneg t, ht,
      mul_nonneg (show (0 : ℝ) ≤ 4 * |t| + 11 by positivity)
        (show (0 : ℝ) ≤ (σ - 3 / 4) - 1 / 4 by linarith)]
  have hinner : 3 + 1 / (|t| - 3 / 4) + (σ + |t| + 3 / 4) / (σ - 3 / 4) ≤ 5 * |t| + 10 := by
    nlinarith [ha, hb, abs_nonneg t, ht]
  have hInnerPos : 0 < 3 + 1 / (|t| - 3 / 4) + (σ + |t| + 3 / 4) / (σ - 3 / 4) := by positivity
  have hfac : 0 < 3 + 1 / (σ - 1) := by
    have : 0 < σ - 1 := by linarith
    positivity
  have hlog : Real.log ((3 + 1 / (|t| - 3 / 4) + (σ + |t| + 3 / 4) / (σ - 3 / 4)) * (3 + 1 / (σ - 1)))
      ≤ Real.log ((5 * |t| + 10) * (3 + 1 / (σ - 1))) :=
    Real.log_le_log (by positivity) (by nlinarith [hinner, hfac, hInnerPos])
  have h32 : 4 * Real.log ((3 + 1 / (|t| - 3 / 4) + (σ + |t| + 3 / 4) / (σ - 3 / 4))
        * (3 + 1 / (σ - 1))) / (3 / 4 / 2)
      = 32 / 3 * Real.log ((3 + 1 / (|t| - 3 / 4) + (σ + |t| + 3 / 4) / (σ - 3 / 4))
        * (3 + 1 / (σ - 1))) := by ring
  rw [h32]; linarith [hlog]

theorem riemannZeta_zero_free_5841_R34 :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ,
      Real.exp 256 ≤ |γ| →
      8 * Real.log (16 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - (221 / 630) / (5 * A₀ + 520 / 3 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841_clean
  refine ⟨A₀, hA0, ?_⟩
  intro γ β hxe hthr hζρ₀
  set L := Real.log |γ| with hLdef
  have hL : 256 ≤ L := by
    calc (256 : ℝ) = Real.log (Real.exp 256) := (Real.log_exp 256).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 256) hxe
  have hγ5 : 5 ≤ |γ| := le_trans (by
      have := Real.add_one_le_exp (256 : ℝ); linarith) hxe
  set 𝒦 := 5 * A₀ + 520 / 3 * L with h𝒦
  have hK : 5 * A₀ + 8 * (40 / 3 * L) + 4 * (40 / 3 * L) + 40 / 3 * L = 𝒦 := by rw [h𝒦]; ring
  have hKpos : 0 < 𝒦 := by rw [h𝒦]; nlinarith
  set σ := 1 + (13 / 10) / 𝒦 with hσdef
  have hq : 0 < (13 / 10) / 𝒦 := by positivity
  have hq1 : (13 / 10) / 𝒦 ≤ 1 := by rw [div_le_one hKpos, h𝒦]; nlinarith
  have hσ1 : 1 < σ := by rw [hσdef]; linarith
  have hσ2 : σ ≤ 2 := by rw [hσdef]; linarith
  have hβσ : β < σ := by
    by_contra hcon
    push Not at hcon
    have hre : 1 < ((β : ℂ) + I * (γ : ℝ)).re := by simp; linarith
    exact riemannZeta_ne_zero_of_one_lt_re hre hζρ₀
  by_cases hfar : β < σ - 3 / 8
  · have h1 : (13 / 10) / 𝒦 + (221 / 630) / 𝒦 ≤ 3 / 8 := by
      rw [← add_div, div_le_iff₀ hKpos, h𝒦]; nlinarith
    linarith
  push Not at hfar
  have hinv : 1 / (σ - 1) = 𝒦 / (13 / 10) := by
    rw [hσdef, show 1 + 13 / 10 / 𝒦 - 1 = 13 / 10 / 𝒦 by ring, one_div_div]
  have hKε : 𝒦 / (13 / 10) ≤ 6 * A₀ + 640 * L := by
    rw [div_le_iff₀ (by norm_num), h𝒦]; nlinarith
  have habs : ∀ t : ℝ, 5 ≤ |t| → 5 * |t| ≤ 15 * |γ| →
      cleanRHS σ t (3 / 4) ≤ 40 / 3 * L := by
    intro t ht htγ
    refine le_trans (cleanRHS_34_le hσ1 hσ2 ht) ?_
    rw [hinv]
    have hγpos : 0 < |γ| := by linarith
    have hd : 0 ≤ 5 * |t| / |γ| := by positivity
    have hd15 : 5 * |t| / |γ| ≤ 15 := by rw [div_le_iff₀ hγpos]; linarith
    have key := log_sharp_close_40 (a := A₀) (x := |γ|) (d := 5 * |t| / |γ|) hA0 hd hd15 hxe hthr
    rw [show 5 * |t| / |γ| * |γ| = 5 * |t| by field_simp] at key
    have hdx : (0 : ℝ) < 5 * |t| + 10 := by positivity
    have hf : 3 + 𝒦 / (13 / 10) ≤ 3 + 6 * A₀ + 640 * L := by linarith
    have hfp : 0 < 3 + 𝒦 / (13 / 10) := by positivity
    have := Real.log_le_log (by positivity) (mul_le_mul_of_nonneg_left hf hdx.le)
    linarith
  have h2γ : |2 * γ| = 2 * |γ| := by rw [abs_mul]; norm_num
  have h3γ : |3 * γ| = 3 * |γ| := by rw [abs_mul]; norm_num
  have hG2 := habs (2 * γ) (by rw [h2γ]; linarith) (by rw [h2γ]; linarith)
  have hG3 := habs (3 * γ) (by rw [h3γ]; linarith) (by rw [h3γ]; linarith)
  have hGB := habs γ hγ5 (by linarith)
  set c : ℂ := (σ : ℂ) + I * (γ : ℝ) with hc_def
  have hcre_eq : c.re = σ := by rw [hc_def]; simp
  have hcre : 1 < c.re := by rw [hcre_eq]; exact hσ1
  have hR : (0 : ℝ) < 3 / 4 := by norm_num
  have hcR : 0 < c.re - 3 / 4 := by rw [hcre_eq]; linarith
  have hc1 : (3 / 4 : ℝ) < ‖c - 1‖ :=
    lt_of_lt_of_le (show (3 / 4 : ℝ) < |γ| by linarith)
      (by rw [hc_def]; exact abs_le_norm_center_sub_one)
  set Ms : ℝ := 2 + 1 / (‖c - 1‖ - 3 / 4) + (‖c‖ + 3 / 4) / (c.re - 3 / 4) with hMs_def
  have hden1 : (0 : ℝ) < ‖c - 1‖ - 3 / 4 := by linarith
  have hMspos : 0 < Ms + 1 := by rw [hMs_def]; positivity
  have hdist : |β - σ| ≤ 3 / 4 / 2 := by rw [abs_sub_comm, abs_of_pos (by linarith)]; linarith
  have hp : ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall c (3 / 4 / 2) := by
    rw [mem_closedBall, dist_eq_norm, hc_def, norm_sub_center]; exact hdist
  have hρball : ((β : ℂ) + I * (γ : ℝ)) ∈ ball c 1 := by
    rw [mem_ball, dist_eq_norm, hc_def, norm_sub_center]; linarith
  have hball : (1 : ℂ) ∉ ball c 1 := by
    rw [hc_def]; exact one_notMem_ball_test (by linarith)
  have hgr : growthRHS σ γ (3 / 4) ≤ 40 / 3 * L :=
    growthRHS_le_clean hσ1 hR (by norm_num) (by linarith) hGB
  have hBexp : 4 * Real.log ((Ms + 1) / ‖riemannZeta c‖) / (3 / 4 / 2) ≤ 40 / 3 * L :=
    hBexplicit_of_ge hcre hMspos hR hgr
  have hgap := H (γ := γ) (β := β) (ε := 13 / 10) (B := 40 / 3 * L) (C₂ := 40 / 3 * L)
    (C₃ := 40 / 3 * L) (σ := σ) (R := 3 / 4) (R' := 1) (Mζ := Ms + 1) (Ms := Ms)
    (Rg := 3 / 4) (Rg' := 1)
    (by norm_num) (by positivity) (by positivity) (by positivity)
    (by rw [hK, h𝒦]; nlinarith) (by rw [hK]) hβσ hR (by norm_num) hball hρball hζρ₀ hp
    (fun z hz => norm_riemannZeta_le_sphere_N1 hR hcR hc1 hz) (by linarith) hMspos hBexp
    hR (by norm_num) (by norm_num) (by rw [h2γ]; linarith) hG2 hG3
  rw [hK] at hgap
  have e : (13 / 10 : ℝ) * (3 - 13 / 10) / ((5 + 13 / 10) * 𝒦) = (221 / 630) / 𝒦 := by
    field_simp; ring
  rw [e] at hgap
  exact hgap

theorem riemannZeta_zero_free_5841_clean_R34 :
    ∃ T₀ : ℝ, ∀ γ β : ℝ, T₀ ≤ |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - (221 / 110880) / Real.log |γ| := by
  obtain ⟨A₀, hA₀nn, hmain⟩ := riemannZeta_zero_free_5841_R34
  refine ⟨max (max (Real.exp 256) (Real.exp (8 * Real.log (16 * (3 + 6 * A₀ + 640)))))
      (Real.exp (15 * A₀ / 8)), ?_⟩
  intro γ β hγ hζ
  have h1 : Real.exp 256 ≤ |γ| := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hγ
  have hexp2 : Real.exp (8 * Real.log (16 * (3 + 6 * A₀ + 640))) ≤ |γ| :=
    le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hγ
  have hexp3 : Real.exp (15 * A₀ / 8) ≤ |γ| := le_trans (le_max_right _ _) hγ
  have h2 : 8 * Real.log (16 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| := by
    calc 8 * Real.log (16 * (3 + 6 * A₀ + 640))
        = Real.log (Real.exp (8 * Real.log (16 * (3 + 6 * A₀ + 640)))) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2
  have hL : 256 ≤ Real.log |γ| := by
    calc (256 : ℝ) = Real.log (Real.exp 256) := (Real.log_exp 256).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 256) h1
  have hA : 15 * A₀ / 8 ≤ Real.log |γ| := by
    calc 15 * A₀ / 8 = Real.log (Real.exp (15 * A₀ / 8)) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp3
  have hmain' := hmain γ β h1 h2 hζ
  have hKpos : 0 < 5 * A₀ + 520 / 3 * Real.log |γ| := by nlinarith
  have hcmp : 5 * A₀ + 520 / 3 * Real.log |γ| ≤ 176 * Real.log |γ| := by nlinarith
  have hdiv : (221 / 630) / (176 * Real.log |γ|)
      ≤ (221 / 630) / (5 * A₀ + 520 / 3 * Real.log |γ|) :=
    div_le_div_of_nonneg_left (by norm_num) hKpos hcmp
  have he : (221 / 110880 : ℝ) / Real.log |γ| = (221 / 630) / (176 * Real.log |γ|) := by
    field_simp; ring
  rw [he]; linarith

end TrigPoly5841

#print axioms TrigPoly5841.cleanRHS_34_le
#print axioms TrigPoly5841.riemannZeta_zero_free_5841_R34
#print axioms TrigPoly5841.riemannZeta_zero_free_5841_clean_R34
