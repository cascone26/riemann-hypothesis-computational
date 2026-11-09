import TrigPoly5841AllSharp
import MulLogLeSelf

/-!
# S2969 — explicit c₀ from the all-sharp 5-8-4-1 region: c₀ = 221/68040 ≈ 3.248e-3 (lane 1)

S2967's region with R = 1/2, R' = 1, ε = 13/10 and B = C₂ = C₃ = (33/4)·log|γ|.  For 1 < σ ≤ 5/4:
`1/(1-σ+1/2) ≤ 4`, `4/(σ-1/2) ≤ 8`, so sharpRHS(σ, kγ, 1/2) ≤ 8·log|γ| + 8 log 3 + 16 log 15
+ 16·log(3+4A₀+83) + 16·log log|γ|, and 16·log L ≤ L/8 (mul_log_le_self, C = 128, L ≥ 65536).
* `riemannZeta_zero_free_5841_sharpc0` — β ≤ 1 − (221/630)/(5A₀ + (429/4)·log|γ|)  (|γ| ≥ e^65536 + threshold).
* `riemannZeta_zero_free_5841_clean_sharpc0` — β ≤ 1 − (221/68040)/log|γ| for |γ| ≥ T₀.
Ladder: 221/74340 ≈ 2.973e-3 (S2957) → 221/68040 ≈ 3.248e-3 (+9.3%), past the old N=1 pipeline ceiling ≈ 3.16e-3.
-/

open Complex Metric TrigPoly5841Sharp

namespace TrigPoly5841SharpC0

theorem norm_sub_center' {β σ γ : ℝ} :
    ‖((β : ℂ) + I * (γ : ℝ)) - ((σ : ℂ) + I * (γ : ℝ))‖ = |β - σ| := by
  rw [show ((β : ℂ) + I * (γ : ℝ)) - ((σ : ℂ) + I * (γ : ℝ)) = ((β - σ : ℝ) : ℂ) by push_cast; ring]
  exact Complex.norm_real _

theorem riemannZeta_zero_free_5841_sharpc0 :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ γ β : ℝ,
      Real.exp 65536 ≤ |γ| →
      8 * (8 * Real.log 3 + 16 * Real.log 15 + 16 * Real.log (3 + 4 * A₀ + 83)) ≤ Real.log |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      β ≤ 1 - (221 / 630) / (5 * A₀ + 429 / 4 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := TrigPoly5841AllSharp.zeta_region_5841_allsharp
  refine ⟨A₀, hA0, ?_⟩
  intro γ β hxe hthr hζρ₀
  set L := Real.log |γ| with hLdef
  have hL : 65536 ≤ L := by
    calc (65536 : ℝ) = Real.log (Real.exp 65536) := (Real.log_exp 65536).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 65536) hxe
  have hγ3 : 3 ≤ |γ| := le_trans (by
      have := Real.add_one_le_exp (65536 : ℝ); linarith) hxe
  have hγpos : 0 < |γ| := by linarith
  set 𝒦 := 5 * A₀ + 429 / 4 * L with h𝒦
  have hK : 5 * A₀ + 8 * (33 / 4 * L) + 4 * (33 / 4 * L) + 33 / 4 * L = 𝒦 := by rw [h𝒦]; ring
  have hKpos : 0 < 𝒦 := by rw [h𝒦]; nlinarith
  set σ := 1 + (13 / 10) / 𝒦 with hσdef
  have hq : 0 < (13 / 10) / 𝒦 := by positivity
  have hq4 : (13 / 10) / 𝒦 ≤ 1 / 4 := by rw [div_le_iff₀ hKpos, h𝒦]; nlinarith
  have hσ1 : 1 < σ := by rw [hσdef]; linarith
  have hσ54 : σ ≤ 5 / 4 := by rw [hσdef]; linarith
  have hβσ : β < σ := by
    by_contra hcon
    push Not at hcon
    have hre : 1 < ((β : ℂ) + I * (γ : ℝ)).re := by simp; linarith
    exact riemannZeta_ne_zero_of_one_lt_re hre hζρ₀
  by_cases hfar : β < σ - 1 / 4
  · have h1 : (13 / 10) / 𝒦 + (221 / 630) / 𝒦 ≤ 1 / 4 := by
      rw [← add_div, div_le_iff₀ hKpos, h𝒦]; nlinarith
    linarith
  push Not at hfar
  have hinv : 1 / (σ - 1) = 𝒦 / (13 / 10) := by
    rw [hσdef, show 1 + 13 / 10 / 𝒦 - 1 = 13 / 10 / 𝒦 by ring, one_div_div]
  have hKε : 𝒦 / (13 / 10) ≤ 4 * A₀ + 83 * L := by
    rw [div_le_iff₀ (by norm_num), h𝒦]; nlinarith
  have hll : 16 * Real.log L ≤ L / 8 := by
    have := mul_log_le_self (C := 128) (x := L) (by norm_num) (by norm_num; linarith)
    linarith
  have hLpos : 0 < L := by linarith
  have hCA : 0 < 3 + 4 * A₀ + 83 := by linarith
  have hD : Real.log (3 + 1 / (σ - 1)) ≤ Real.log (3 + 4 * A₀ + 83) + Real.log L := by
    rw [hinv, ← Real.log_mul hCA.ne' hLpos.ne']
    have hpos : 0 < 3 + 𝒦 / (13 / 10) := by positivity
    refine Real.log_le_log hpos ?_
    nlinarith
  have hK15 : Real.log (1 / (1 - (σ - 1 / 2)) + 3 + 4 / (σ - 1 / 2)) ≤ Real.log 15 := by
    have ha : 0 < 1 - (σ - 1 / 2) := by linarith
    have hb : 0 < σ - 1 / 2 := by linarith
    have h1 : 1 / (1 - (σ - 1 / 2)) ≤ 4 := by rw [div_le_iff₀ ha]; linarith
    have h2 : 4 / (σ - 1 / 2) ≤ 8 := by rw [div_le_iff₀ hb]; linarith
    refine Real.log_le_log (by positivity) ?_
    linarith
  have habs : ∀ t : ℝ, 3 ≤ |t| → |t| ≤ 3 * |γ| → sharpRHS σ t (1 / 2) ≤ 33 / 4 * L := by
    intro t ht htγ
    unfold sharpRHS
    have htpos : 0 < |t| := by linarith
    have hlt : Real.log |t| ≤ Real.log 3 + L := by
      rw [hLdef, ← Real.log_mul (by norm_num) hγpos.ne']
      exact Real.log_le_log htpos htγ
    rw [show (8 : ℝ) / (1 / 2) = 16 by norm_num]
    nlinarith [hlt, hK15, hD, hll, hthr]
  have h2γ : |2 * γ| = 2 * |γ| := by rw [abs_mul]; norm_num
  have h3γ : |3 * γ| = 3 * |γ| := by rw [abs_mul]; norm_num
  have hG2 := habs (2 * γ) (by rw [h2γ]; linarith) (by rw [h2γ]; linarith)
  have hG3 := habs (3 * γ) (by rw [h3γ]; linarith) (by rw [h3γ])
  have hGB := habs γ hγ3 (by linarith)
  have hp : ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall ((σ : ℂ) + I * (γ : ℝ)) (1 / 2 / 2) := by
    rw [mem_closedBall, dist_eq_norm, norm_sub_center', abs_sub_comm, abs_of_pos (by linarith)]
    linarith
  have hρball : ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) 1 := by
    rw [mem_ball, dist_eq_norm, norm_sub_center', abs_sub_comm, abs_of_pos (by linarith)]
    linarith
  have hgap := H (γ := γ) (β := β) (ε := 13 / 10) (B := 33 / 4 * L) (C₂ := 33 / 4 * L)
    (C₃ := 33 / 4 * L) (σ := σ) (R := 1 / 2) (R' := 1)
    (by norm_num) (by positivity) (by positivity) (by positivity)
    (by rw [hK, h𝒦]; nlinarith) (by rw [hK]) hβσ (by norm_num) (by norm_num) (by linarith)
    (by norm_num) (by linarith) (by linarith) hγ3 hρball hζρ₀ hp hGB hG2 hG3
  rw [hK] at hgap
  have e : (13 / 10 : ℝ) * (3 - 13 / 10) / ((5 + 13 / 10) * 𝒦) = (221 / 630) / 𝒦 := by
    field_simp; ring
  rw [e] at hgap
  exact hgap

theorem riemannZeta_zero_free_5841_clean_sharpc0 :
    ∃ T₀ : ℝ, ∀ γ β : ℝ, T₀ ≤ |γ| →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 → β ≤ 1 - (221 / 68040) / Real.log |γ| := by
  obtain ⟨A₀, hA₀nn, hmain⟩ := riemannZeta_zero_free_5841_sharpc0
  set Cst := 8 * (8 * Real.log 3 + 16 * Real.log 15 + 16 * Real.log (3 + 4 * A₀ + 83)) with hCst
  refine ⟨max (max (Real.exp 65536) (Real.exp Cst)) (Real.exp (20 * A₀ / 3)), ?_⟩
  intro γ β hγ hζ
  have h1 : Real.exp 65536 ≤ |γ| := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hγ
  have hexp2 : Real.exp Cst ≤ |γ| := le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hγ
  have hexp3 : Real.exp (20 * A₀ / 3) ≤ |γ| := le_trans (le_max_right _ _) hγ
  have h2 : Cst ≤ Real.log |γ| := by
    calc Cst = Real.log (Real.exp Cst) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp2
  have hL : 65536 ≤ Real.log |γ| := by
    calc (65536 : ℝ) = Real.log (Real.exp 65536) := (Real.log_exp 65536).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos 65536) h1
  have hA : 20 * A₀ / 3 ≤ Real.log |γ| := by
    calc 20 * A₀ / 3 = Real.log (Real.exp (20 * A₀ / 3)) := (Real.log_exp _).symm
      _ ≤ Real.log |γ| := Real.log_le_log (Real.exp_pos _) hexp3
  have hmain' := hmain γ β h1 h2 hζ
  have hKpos : 0 < 5 * A₀ + 429 / 4 * Real.log |γ| := by nlinarith
  have hcmp : 5 * A₀ + 429 / 4 * Real.log |γ| ≤ 108 * Real.log |γ| := by nlinarith
  have hdiv : (221 / 630) / (108 * Real.log |γ|)
      ≤ (221 / 630) / (5 * A₀ + 429 / 4 * Real.log |γ|) :=
    div_le_div_of_nonneg_left (by norm_num) hKpos hcmp
  have he : (221 / 68040 : ℝ) / Real.log |γ| = (221 / 630) / (108 * Real.log |γ|) := by
    field_simp; ring
  rw [he]; linarith

end TrigPoly5841SharpC0

#print axioms TrigPoly5841SharpC0.riemannZeta_zero_free_5841_sharpc0
#print axioms TrigPoly5841SharpC0.riemannZeta_zero_free_5841_clean_sharpc0
