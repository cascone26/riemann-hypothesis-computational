import TrigPoly5841R34

/-!
# S2955 — 5-8-4-1 region at Jensen radius R = 15/16: c₀ = 221/88200 ≈ 2.51e-3 (lane 1)

S2954 `log_sharp_close_40K`: the closing estimate with coefficient room d ≤ 63 (base 64 instead of 16).
`cleanRHS_1516_le`: 1<σ≤2, |t| ≥ 41 ⟹ cleanRHS σ t (15/16) ≤ (128/15)·log((17|t|+10)(3+1/(σ−1))) (Jensen cost 8/R = 128/15).
Heights γ,2γ,3γ ⟹ d = 17,34,51 ≤ 63 ⟹ B = C₂ = C₃ = (32/3)·log|γ|, 𝒦 = 5A₀ + (416/3)·log|γ|, ε = 13/10.
* `riemannZeta_zero_free_5841_R1516` : β ≤ 1 − (221/630)/(5A₀ + (416/3)·log|γ|);
* `riemannZeta_zero_free_5841_clean_R1516` : β ≤ 1 − (221/88200)/log|γ| for |γ| ≥ T₀.
-/

open Complex Metric

theorem log_sharp_close_40K {a x d : ℝ} (ha : 0 ≤ a) (hd : 0 ≤ d) (hd63 : d ≤ 63)
    (hxe : Real.exp 256 ≤ x)
    (hthr : 8 * Real.log (64 * (3 + 6 * a + 640)) ≤ Real.log x) :
    32 * Real.log ((d * x + 10) * (3 + 6 * a + 640 * Real.log x)) ≤ 40 * Real.log x := by
  have hx10 : (10 : ℝ) ≤ x := le_trans (by nlinarith [Real.add_one_le_exp (256 : ℝ)]) hxe
  have hx1 : (1 : ℝ) ≤ x := by linarith
  have hxpos : 0 < x := by linarith
  have hL : 256 ≤ Real.log x := by
    calc (256 : ℝ) = Real.log (Real.exp 256) := (Real.log_exp 256).symm
      _ ≤ Real.log x := Real.log_le_log (Real.exp_pos 256) hxe
  have hLpos : 0 < Real.log x := by linarith
  have hbase : 0 < 64 * (3 + 6 * a + 640) := by nlinarith [ha]
  have hslack : (0:ℝ) ≤ (64 - d) * x - 10 := by nlinarith [hx10, hd63]
  have h1 : d * x + 10 ≤ 64 * x := by linarith
  have h2 : 3 + 6 * a + 640 * Real.log x ≤ (3 + 6 * a + 640) * Real.log x := by nlinarith
  have hprod : (d * x + 10) * (3 + 6 * a + 640 * Real.log x)
      ≤ 64 * (3 + 6 * a + 640) * (x * Real.log x) := by
    have hA : 0 ≤ d * x + 10 := by nlinarith [mul_nonneg hd hxpos.le]
    have hB : 0 ≤ 3 + 6 * a + 640 * Real.log x := by nlinarith
    calc (d * x + 10) * (3 + 6 * a + 640 * Real.log x)
        ≤ (64 * x) * ((3 + 6 * a + 640) * Real.log x) := mul_le_mul h1 h2 hB (by positivity)
      _ = 64 * (3 + 6 * a + 640) * (x * Real.log x) := by ring
  have hargpos : 0 < (d * x + 10) * (3 + 6 * a + 640 * Real.log x) :=
    mul_pos (by nlinarith [mul_nonneg hd (le_of_lt hxpos)]) (by nlinarith [ha, hLpos])
  have hlog1 : Real.log ((d * x + 10) * (3 + 6 * a + 640 * Real.log x))
      ≤ Real.log (64 * (3 + 6 * a + 640) * (x * Real.log x)) := Real.log_le_log hargpos hprod
  have hlog2 : Real.log (64 * (3 + 6 * a + 640) * (x * Real.log x))
      = Real.log (64 * (3 + 6 * a + 640)) + Real.log x + Real.log (Real.log x) := by
    rw [Real.log_mul (ne_of_gt hbase) (by positivity), Real.log_mul (ne_of_gt hxpos) (ne_of_gt hLpos)]
    ring
  have hloglog : Real.log (Real.log x) ≤ Real.log x / 8 := by
    have h := mul_log_le_self (C := 8) (x := Real.log x) (by norm_num) (by norm_num; linarith)
    linarith
  linarith [hlog1, hlog2, hthr, hloglog, hL]

namespace TrigPoly5841

theorem cleanRHS_1516_le {σ t : ℝ} (hσ1 : 1 < σ) (hσ2 : σ ≤ 2) (ht : 41 ≤ |t|) :
    cleanRHS σ t (15 / 16) ≤ 128 / 15 * Real.log ((17 * |t| + 10) * (3 + 1 / (σ - 1))) := by
  unfold cleanRHS
  have hden_pos : 0 < σ - 15 / 16 := by linarith
  have htpos : 0 < |t| - 15 / 16 := by linarith
  have ha : 1 / (|t| - 15 / 16) ≤ 1 := by rw [div_le_iff₀ htpos]; linarith
  have hb : (σ + |t| + 15 / 16) / (σ - 15 / 16) ≤ 16 * |t| + 47 := by
    rw [div_le_iff₀ hden_pos]
    nlinarith [hσ2, hden_pos, abs_nonneg t, ht,
      mul_nonneg (show (0 : ℝ) ≤ 16 * |t| + 47 by positivity)
        (show (0 : ℝ) ≤ (σ - 15 / 16) - 1 / 16 by linarith)]
  have hinner : 3 + 1 / (|t| - 15 / 16) + (σ + |t| + 15 / 16) / (σ - 15 / 16) ≤ 17 * |t| + 10 := by
    nlinarith [ha, hb, abs_nonneg t, ht]
  have hInnerPos : 0 < 3 + 1 / (|t| - 15 / 16) + (σ + |t| + 15 / 16) / (σ - 15 / 16) := by positivity
  have hfac : 0 < 3 + 1 / (σ - 1) := by
    have : 0 < σ - 1 := by linarith
    positivity
  have hlog : Real.log ((3 + 1 / (|t| - 15 / 16) + (σ + |t| + 15 / 16) / (σ - 15 / 16)) * (3 + 1 / (σ - 1)))
      ≤ Real.log ((17 * |t| + 10) * (3 + 1 / (σ - 1))) :=
    Real.log_le_log (by positivity) (by nlinarith [hinner, hfac, hInnerPos])
  have h32 : 4 * Real.log ((3 + 1 / (|t| - 15 / 16) + (σ + |t| + 15 / 16) / (σ - 15 / 16))
        * (3 + 1 / (σ - 1))) / (15 / 16 / 2)
      = 128 / 15 * Real.log ((3 + 1 / (|t| - 15 / 16) + (σ + |t| + 15 / 16) / (σ - 15 / 16))
        * (3 + 1 / (σ - 1))) := by ring
  rw [h32]; linarith [hlog]

theorem riemannZeta_zero_free_5841_R1516 :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ,
      Real.exp 256 ≤ |γ| →
      8 * Real.log (64 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - (221 / 630) / (5 * A₀ + 416 / 3 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841_clean
  refine ⟨A₀, hA0, ?_⟩
  intro γ β hxe hthr hζρ₀
  set L := Real.log |γ| with hLdef
  have hL : 256 ≤ L := by
    calc (256 : ℝ) = Real.log (Real.exp 256) := (Real.log_exp 256).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 256) hxe
  have hγ5 : 41 ≤ |γ| := le_trans (by
      have := Real.add_one_le_exp (256 : ℝ); linarith) hxe
  set 𝒦 := 5 * A₀ + 416 / 3 * L with h𝒦
  have hK : 5 * A₀ + 8 * (32 / 3 * L) + 4 * (32 / 3 * L) + 32 / 3 * L = 𝒦 := by rw [h𝒦]; ring
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
  by_cases hfar : β < σ - 15 / 32
  · have h1 : (13 / 10) / 𝒦 + (221 / 630) / 𝒦 ≤ 15 / 32 := by
      rw [← add_div, div_le_iff₀ hKpos, h𝒦]; nlinarith
    linarith
  push Not at hfar
  have hinv : 1 / (σ - 1) = 𝒦 / (13 / 10) := by
    rw [hσdef, show 1 + 13 / 10 / 𝒦 - 1 = 13 / 10 / 𝒦 by ring, one_div_div]
  have hKε : 𝒦 / (13 / 10) ≤ 6 * A₀ + 640 * L := by
    rw [div_le_iff₀ (by norm_num), h𝒦]; nlinarith
  have habs : ∀ t : ℝ, 41 ≤ |t| → 17 * |t| ≤ 51 * |γ| →
      cleanRHS σ t (15 / 16) ≤ 32 / 3 * L := by
    intro t ht htγ
    refine le_trans (cleanRHS_1516_le hσ1 hσ2 ht) ?_
    rw [hinv]
    have hγpos : 0 < |γ| := by linarith
    have hd : 0 ≤ 17 * |t| / |γ| := by positivity
    have hd15 : 17 * |t| / |γ| ≤ 63 := by rw [div_le_iff₀ hγpos]; linarith
    have key := log_sharp_close_40K (a := A₀) (x := |γ|) (d := 17 * |t| / |γ|) hA0 hd hd15 hxe hthr
    rw [show 17 * |t| / |γ| * |γ| = 17 * |t| by field_simp] at key
    have hdx : (0 : ℝ) < 17 * |t| + 10 := by positivity
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
  have hR : (0 : ℝ) < 15 / 16 := by norm_num
  have hcR : 0 < c.re - 15 / 16 := by rw [hcre_eq]; linarith
  have hc1 : (15 / 16 : ℝ) < ‖c - 1‖ :=
    lt_of_lt_of_le (show (15 / 16 : ℝ) < |γ| by linarith)
      (by rw [hc_def]; exact abs_le_norm_center_sub_one)
  set Ms : ℝ := 2 + 1 / (‖c - 1‖ - 15 / 16) + (‖c‖ + 15 / 16) / (c.re - 15 / 16) with hMs_def
  have hden1 : (0 : ℝ) < ‖c - 1‖ - 15 / 16 := by linarith
  have hMspos : 0 < Ms + 1 := by rw [hMs_def]; positivity
  have hdist : |β - σ| ≤ 15 / 16 / 2 := by rw [abs_sub_comm, abs_of_pos (by linarith)]; linarith
  have hp : ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall c (15 / 16 / 2) := by
    rw [mem_closedBall, dist_eq_norm, hc_def, norm_sub_center]; exact hdist
  have hρball : ((β : ℂ) + I * (γ : ℝ)) ∈ ball c 1 := by
    rw [mem_ball, dist_eq_norm, hc_def, norm_sub_center]; linarith
  have hball : (1 : ℂ) ∉ ball c 1 := by
    rw [hc_def]; exact one_notMem_ball_test (by linarith)
  have hgr : growthRHS σ γ (15 / 16) ≤ 32 / 3 * L :=
    growthRHS_le_clean hσ1 hR (by norm_num) (by linarith) hGB
  have hBexp : 4 * Real.log ((Ms + 1) / ‖riemannZeta c‖) / (15 / 16 / 2) ≤ 32 / 3 * L :=
    hBexplicit_of_ge hcre hMspos hR hgr
  have hgap := H (γ := γ) (β := β) (ε := 13 / 10) (B := 32 / 3 * L) (C₂ := 32 / 3 * L)
    (C₃ := 32 / 3 * L) (σ := σ) (R := 15 / 16) (R' := 1) (Mζ := Ms + 1) (Ms := Ms)
    (Rg := 15 / 16) (Rg' := 1)
    (by norm_num) (by positivity) (by positivity) (by positivity)
    (by rw [hK, h𝒦]; nlinarith) (by rw [hK]) hβσ hR (by norm_num) hball hρball hζρ₀ hp
    (fun z hz => norm_riemannZeta_le_sphere_N1 hR hcR hc1 hz) (by linarith) hMspos hBexp
    hR (by norm_num) (by norm_num) (by rw [h2γ]; linarith) hG2 hG3
  rw [hK] at hgap
  have e : (13 / 10 : ℝ) * (3 - 13 / 10) / ((5 + 13 / 10) * 𝒦) = (221 / 630) / 𝒦 := by
    field_simp; ring
  rw [e] at hgap
  exact hgap

theorem riemannZeta_zero_free_5841_clean_R1516 :
    ∃ T₀ : ℝ, ∀ γ β : ℝ, T₀ ≤ |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - (221 / 88200) / Real.log |γ| := by
  obtain ⟨A₀, hA₀nn, hmain⟩ := riemannZeta_zero_free_5841_R1516
  refine ⟨max (max (Real.exp 256) (Real.exp (8 * Real.log (64 * (3 + 6 * A₀ + 640)))))
      (Real.exp (15 * A₀ / 4)), ?_⟩
  intro γ β hγ hζ
  have h1 : Real.exp 256 ≤ |γ| := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hγ
  have hexp2 : Real.exp (8 * Real.log (64 * (3 + 6 * A₀ + 640))) ≤ |γ| :=
    le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hγ
  have hexp3 : Real.exp (15 * A₀ / 4) ≤ |γ| := le_trans (le_max_right _ _) hγ
  have h2 : 8 * Real.log (64 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| := by
    calc 8 * Real.log (64 * (3 + 6 * A₀ + 640))
        = Real.log (Real.exp (8 * Real.log (64 * (3 + 6 * A₀ + 640)))) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2
  have hL : 256 ≤ Real.log |γ| := by
    calc (256 : ℝ) = Real.log (Real.exp 256) := (Real.log_exp 256).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 256) h1
  have hA : 15 * A₀ / 4 ≤ Real.log |γ| := by
    calc 15 * A₀ / 4 = Real.log (Real.exp (15 * A₀ / 4)) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp3
  have hmain' := hmain γ β h1 h2 hζ
  have hKpos : 0 < 5 * A₀ + 416 / 3 * Real.log |γ| := by nlinarith
  have hcmp : 5 * A₀ + 416 / 3 * Real.log |γ| ≤ 140 * Real.log |γ| := by nlinarith
  have hdiv : (221 / 630) / (140 * Real.log |γ|)
      ≤ (221 / 630) / (5 * A₀ + 416 / 3 * Real.log |γ|) :=
    div_le_div_of_nonneg_left (by norm_num) hKpos hcmp
  have he : (221 / 88200 : ℝ) / Real.log |γ| = (221 / 630) / (140 * Real.log |γ|) := by
    field_simp; ring
  rw [he]; linarith

end TrigPoly5841

#print axioms log_sharp_close_40K
#print axioms TrigPoly5841.riemannZeta_zero_free_5841_R1516
#print axioms TrigPoly5841.riemannZeta_zero_free_5841_clean_R1516
