import TrigPoly5841Clean

/-!
# S2940 — 5-8-4-1 region with the growth inputs in log-absorption form (lane 1)

σ-generic version of §1460 `hBexp_rhs_le`: for 1 < σ ≤ 2 and |t| > 1, at radius R = 1/4,
  cleanRHS σ t (1/4) ≤ 32·log((3|t| + 10)·(3 + 1/(σ−1))).
At σ = 1 + ε/𝒦 we have 1/(σ−1) = 𝒦/ε, so the 5-8-4-1 growth hypotheses become
  32·log((3|2γ|+10)(3+𝒦/ε)) ≤ C₂,   32·log((3|3γ|+10)(3+𝒦/ε)) ≤ C₃,
i.e. exactly the `(d·x+10)(3+…)` shape of S2930 `log_sharp_close` (d = 6) / `log_sharp_close_wide` (d = 9).
-/

open Complex Metric

namespace TrigPoly5841

theorem cleanRHS_quarter_le {σ t : ℝ} (hσ1 : 1 < σ) (hσ2 : σ ≤ 2) (ht : 1 < |t|) :
    cleanRHS σ t (1 / 4) ≤ 32 * Real.log ((3 * |t| + 10) * (3 + 1 / (σ - 1))) := by
  unfold cleanRHS
  have hden_pos : 0 < σ - 1 / 4 := by linarith
  have htpos : 0 < |t| - 1 / 4 := by linarith
  have ha : 1 / (|t| - 1 / 4) ≤ 4 / 3 := by rw [div_le_iff₀ htpos]; nlinarith [ht]
  have hb : (σ + |t| + 1 / 4) / (σ - 1 / 4) ≤ (4 / 3) * (|t| + 9 / 4) := by
    rw [div_le_iff₀ hden_pos]
    nlinarith [hσ2, hden_pos, abs_nonneg t, ht,
      mul_nonneg (show (0 : ℝ) ≤ (4 / 3) * (|t| + 9 / 4) by positivity)
        (show (0 : ℝ) ≤ (σ - 1 / 4) - 3 / 4 by linarith)]
  have hinner : 3 + 1 / (|t| - 1 / 4) + (σ + |t| + 1 / 4) / (σ - 1 / 4) ≤ 3 * |t| + 10 := by
    nlinarith [ha, hb, abs_nonneg t, ht]
  have hInnerPos : 0 < 3 + 1 / (|t| - 1 / 4) + (σ + |t| + 1 / 4) / (σ - 1 / 4) := by positivity
  have hfac : 0 < 3 + 1 / (σ - 1) := by have : 0 < σ - 1 := by linarith
                                        positivity
  have hlog : Real.log ((3 + 1 / (|t| - 1 / 4) + (σ + |t| + 1 / 4) / (σ - 1 / 4)) * (3 + 1 / (σ - 1)))
      ≤ Real.log ((3 * |t| + 10) * (3 + 1 / (σ - 1))) :=
    Real.log_le_log (by positivity) (by nlinarith [hinner, hfac, hInnerPos])
  have h32 : 4 * Real.log ((3 + 1 / (|t| - 1 / 4) + (σ + |t| + 1 / 4) / (σ - 1 / 4))
        * (3 + 1 / (σ - 1))) / (1 / 4 / 2)
      = 32 * Real.log ((3 + 1 / (|t| - 1 / 4) + (σ + |t| + 1 / 4) / (σ - 1 / 4))
        * (3 + 1 / (σ - 1))) := by ring
  rw [h32]; linarith [hlog]

theorem zeta_region_5841_scalar :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ {γ β ε B C₂ C₃ σ R R' Mζ Ms : ℝ},
      0 < ε → 0 ≤ B → 0 ≤ C₂ → 0 ≤ C₃ →
      ε ≤ 5 * A₀ + 8 * B + 4 * C₂ + C₃ →
      σ = 1 + ε / (5 * A₀ + 8 * B + 4 * C₂ + C₃) →
      β < σ →
      0 < R → R < R' →
      (1 : ℂ) ∉ ball ((σ : ℂ) + I * (γ : ℝ)) R' →
      ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) R' →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall ((σ : ℂ) + I * (γ : ℝ)) (R / 2) →
      (∀ z ∈ sphere ((σ : ℂ) + I * (γ : ℝ)) R, ‖riemannZeta z‖ ≤ Ms) →
      Ms < Mζ → 0 < Mζ →
      4 * Real.log (Mζ / ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B →
      1 < |γ| →
      32 * Real.log ((3 * |2 * γ| + 10) * (3 + (5 * A₀ + 8 * B + 4 * C₂ + C₃) / ε)) ≤ C₂ →
      32 * Real.log ((3 * |3 * γ| + 10) * (3 + (5 * A₀ + 8 * B + 4 * C₂ + C₃) / ε)) ≤ C₃ →
      β ≤ 1 - ε * (3 - ε) / ((5 + ε) * (5 * A₀ + 8 * B + 4 * C₂ + C₃)) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841_clean
  refine ⟨A₀, hA0, ?_⟩
  intro γ β ε B C₂ C₃ σ R R' Mζ Ms hε hB hC2 hC3 hεK hσdef hβσ hR hRR' hball hρ₀ball hζρ₀
    hp hζsphere hsl hMζ hBexp hγ hG2 hG3
  set 𝒦 := 5 * A₀ + 8 * B + 4 * C₂ + C₃ with h𝒦def
  have hKpos : 0 < 𝒦 := lt_of_lt_of_le hε hεK
  have hq : 0 < ε / 𝒦 := div_pos hε hKpos
  have hq1 : ε / 𝒦 ≤ 1 := (div_le_one hKpos).mpr hεK
  have hσ1 : 1 < σ := by rw [hσdef]; linarith
  have hσ2 : σ ≤ 2 := by rw [hσdef]; linarith
  have hinv : 1 / (σ - 1) = 𝒦 / ε := by
    rw [hσdef, show 1 + ε / 𝒦 - 1 = ε / 𝒦 by ring, one_div_div]
  have h2 : 1 < |2 * γ| := by rw [abs_mul]; norm_num; linarith [abs_nonneg γ]
  have h3 : 1 < |3 * γ| := by rw [abs_mul]; norm_num; linarith [abs_nonneg γ]
  have hRg' : (1 / 2 : ℝ) ≤ |2 * γ| := by linarith
  refine H (Rg := 1 / 4) (Rg' := 1 / 2) hε hB hC2 hC3 hεK hσdef hβσ hR hRR' hball hρ₀ball hζρ₀ hp
    hζsphere hsl hMζ hBexp (by norm_num) (by norm_num) (by norm_num) hRg' ?_ ?_
  · refine le_trans (cleanRHS_quarter_le hσ1 hσ2 h2) ?_
    rw [hinv]; exact hG2
  · refine le_trans (cleanRHS_quarter_le hσ1 hσ2 h3) ?_
    rw [hinv]; exact hG3

end TrigPoly5841

#print axioms TrigPoly5841.cleanRHS_quarter_le
#print axioms TrigPoly5841.zeta_region_5841_scalar
