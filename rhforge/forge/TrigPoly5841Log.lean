import TrigPoly5841Scalar
import LogSharpClose
import LogSharpCloseWide

/-!
# S2942 — 5-8-4-1 region with the GROWTH constants discharged: C₂ = C₃ = 64·log|γ| (lane 1)

Instantiate S2940 at ε = 3/2, B = C₂ = C₃ = 64·L (L = log|γ|), so 𝒦 = 5A₀ + 832·L and 𝒦/ε ≤ 6A₀ + 640·L.
The two growth hypotheses 32·log((3|kγ|+10)(3+𝒦/ε)) ≤ 64·L (k = 2, 3, i.e. d = 6, 9) are then closed by
S2930 `log_sharp_close` (d = 6) and `log_sharp_close_wide` (d = 9) for |γ| ≥ e¹⁶ past the A₀-threshold.
Result: every ζ zero β+iγ in that range satisfying the zero-input sphere data with B = 64·log|γ| obeys
  β ≤ 1 − (9/26)/(5A₀ + 832·log|γ|)        — the 5-8-4-1 σ > 1 − c/log|γ| shape.
Only the zero-input sphere data (the iγ term) remains as hypothesis.
-/

open Complex Metric

namespace TrigPoly5841

theorem zeta_region_5841_log :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ {γ β σ R R' Mζ Ms : ℝ},
      Real.exp 16 ≤ |γ| →
      2 * Real.log (16 * (3 + 6 * A₀ + 640)) ≤ Real.log |γ| →
      σ = 1 + (3 / 2) / (5 * A₀ + 832 * Real.log |γ|) →
      β < σ →
      0 < R → R < R' →
      (1 : ℂ) ∉ ball ((σ : ℂ) + I * (γ : ℝ)) R' →
      ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) R' →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall ((σ : ℂ) + I * (γ : ℝ)) (R / 2) →
      (∀ z ∈ sphere ((σ : ℂ) + I * (γ : ℝ)) R, ‖riemannZeta z‖ ≤ Ms) →
      Ms < Mζ → 0 < Mζ →
      4 * Real.log (Mζ / ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ 64 * Real.log |γ| →
      β ≤ 1 - (9 / 26) / (5 * A₀ + 832 * Real.log |γ|) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841_scalar
  refine ⟨A₀, hA0, ?_⟩
  intro γ β σ R R' Mζ Ms hxe hthr hσdef hβσ hR hRR' hball hρ₀ball hζρ₀ hp hζsphere hsl hMζ hBexp
  set L := Real.log |γ| with hLdef
  have hL16 : 16 ≤ L := by
    calc (16 : ℝ) = Real.log (Real.exp 16) := (Real.log_exp 16).symm
      _ ≤ L := Real.log_le_log (Real.exp_pos 16) hxe
  have hγ1 : 1 < |γ| := lt_of_lt_of_le (by
      have := Real.add_one_le_exp (16 : ℝ); linarith) hxe
  have hK : 5 * A₀ + 8 * (64 * L) + 4 * (64 * L) + 64 * L = 5 * A₀ + 832 * L := by ring
  have hKpos : 0 < 5 * A₀ + 832 * L := by nlinarith
  have hKε : (5 * A₀ + 832 * L) / (3 / 2) ≤ 6 * A₀ + 640 * L := by
    rw [div_le_iff₀ (by norm_num)]; nlinarith
  have hfac : 3 + (5 * A₀ + 832 * L) / (3 / 2) ≤ 3 + 6 * A₀ + 640 * L := by linarith
  have hfacpos : 0 < 3 + (5 * A₀ + 832 * L) / (3 / 2) := by positivity
  have hgrow : ∀ d : ℝ, 0 ≤ d →
      32 * Real.log ((d * |γ| + 10) * (3 + 6 * A₀ + 640 * L)) ≤ 64 * L →
      32 * Real.log ((d * |γ| + 10) * (3 + (5 * A₀ + 832 * L) / (3 / 2))) ≤ 64 * L := by
    intro d hd h
    refine le_trans ?_ h
    have hdx : 0 < d * |γ| + 10 := by positivity
    have := Real.log_le_log (by positivity) (mul_le_mul_of_nonneg_left hfac hdx.le)
    linarith
  have e2 : 3 * |2 * γ| = 6 * |γ| := by rw [abs_mul]; norm_num; ring
  have e3 : 3 * |3 * γ| = 9 * |γ| := by rw [abs_mul]; norm_num; ring
  have g6 := hgrow 6 (by norm_num) (log_sharp_close hA0 (by norm_num) (by norm_num) hxe hthr)
  have g9 := hgrow 9 (by norm_num) (log_sharp_close_wide hA0 (by norm_num) (by norm_num) hxe hthr)
  have hgap := H (γ := γ) (β := β) (ε := 3 / 2) (B := 64 * L) (C₂ := 64 * L) (C₃ := 64 * L)
    (σ := σ) (R := R) (R' := R') (Mζ := Mζ) (Ms := Ms)
    (by norm_num) (by positivity) (by positivity) (by positivity)
    (by rw [hK]; nlinarith) (by rw [hK]; exact hσdef) hβσ hR hRR' hball hρ₀ball hζρ₀ hp hζsphere
    hsl hMζ hBexp hγ1
    (by rw [hK, e2]; exact g6) (by rw [hK, e3]; exact g9)
  rw [hK] at hgap
  have e : (3 / 2 : ℝ) * (3 - 3 / 2) / ((5 + 3 / 2) * (5 * A₀ + 832 * L))
      = (9 / 26) / (5 * A₀ + 832 * L) := by
    field_simp; ring
  rw [e] at hgap
  exact hgap

end TrigPoly5841

#print axioms TrigPoly5841.zeta_region_5841_log
