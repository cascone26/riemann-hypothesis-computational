import TrigPoly5841Region
import GrowthSharpClean

/-!
# S2965 — 5-8-4-1 zero-free region driven by the SHARP growth input (lane 1)

* `zeta_region_5841_Tz` — the S2936 region with the growth input stated directly as `Tz σ γ 2 ≤ C₂`,
  `Tz σ γ 3 ≤ C₃` (decoupled from any particular sphere bound).
* `sharpRHS σ t R := 8·log|t| + (8/R)·(log(1/(1-σ+R)+3+4/(σ-R)) + log(3+1/(σ-1)))`  (S2963).
* `zeta_region_5841_sharp` — region with `sharpRHS σ (2γ) R ≤ C₂`, `sharpRHS σ (3γ) R ≤ C₃`,
  for `σ ≤ 2`, `0 < R < 1`, `σ - R < 1`, `3 ≤ |2γ|`.  Leading growth coefficient 8 (vs 8/R before).
-/

open Complex Metric

namespace TrigPoly5841Sharp

noncomputable def sharpRHS (σ t R : ℝ) : ℝ :=
  8 * Real.log |t| + 8 / R * (Real.log (1 / (1 - (σ - R)) + 3 + 4 / (σ - R))
    + Real.log (3 + 1 / (σ - 1)))

theorem zeta_region_5841_Tz :
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
      TrigPolyRegion.Tz σ γ 2 ≤ C₂ → TrigPolyRegion.Tz σ γ 3 ≤ C₃ →
      β ≤ 1 - ε * (3 - ε) / ((5 + ε) * (5 * A₀ + 8 * B + 4 * C₂ + C₃)) := by
  obtain ⟨A₀, hA0, H⟩ := TrigPolyRegion.zeta_region_P_pole
  refine ⟨A₀, hA0, ?_⟩
  intro γ β ε B C₂ C₃ σ R R' Mζ Ms hε hB hC2 hC3 hεK hσdef hβσ hR hRR' hball hρ₀ball hζρ₀
    hp hζsphere hsl hMζ hBexp hT2 hT3
  set 𝒦 := 5 * A₀ + 8 * B + 4 * C₂ + C₃ with h𝒦def
  have hKpos : 0 < 𝒦 := lt_of_lt_of_le hε hεK
  have hσ : 1 < σ := by rw [hσdef]; have : 0 < ε / 𝒦 := div_pos hε hKpos; linarith
  have h𝒦 : (![5, 8, 4, 1] : Fin 4 → ℝ) 0 * A₀ + (![5, 8, 4, 1] : Fin 4 → ℝ) 1 * B
      + ∑ k : Fin 2, (![5, 8, 4, 1] : Fin 4 → ℝ) k.succ.succ * (![C₂, C₃] : Fin 2 → ℝ) k = 𝒦 := by
    simp [Fin.sum_univ_two, h𝒦def]; ring
  have hσ' : σ = 1 + ε / ((![5, 8, 4, 1] : Fin 4 → ℝ) 0 * A₀ + (![5, 8, 4, 1] : Fin 4 → ℝ) 1 * B
      + ∑ k : Fin 2, (![5, 8, 4, 1] : Fin 4 → ℝ) k.succ.succ * (![C₂, C₃] : Fin 2 → ℝ) k) := by
    rw [h𝒦]; exact hσdef
  have hgap := H (γ := γ) (β := β) (ε := ε) (B := B) (K := 2) ![5, 8, 4, 1] ![C₂, C₃]
    (fun θ => by
      have := TrigPolyDLVP.p5841_nonneg θ
      simp [Fin.sum_univ_four]; linarith)
    (by norm_num) (by norm_num) (fun k => by fin_cases k <;> norm_num) hε hB
    (fun k => by fin_cases k <;> simpa)
    (by rw [h𝒦]; exact hKpos) (by rw [h𝒦]; exact hεK) (by rw [← hσ']; exact hβσ)
    (by
      rw [← hσ']
      exact TrigPolyRegion.hzero_Tz hR hRR' hσ hβσ hball hρ₀ball hζρ₀ hp hζsphere hsl hMζ hBexp)
    (by
      rw [← hσ']
      intro k
      fin_cases k
      · exact le_trans (by simpa using hT2) le_rfl
      · exact le_trans (by simpa using hT3) le_rfl)
  rw [h𝒦] at hgap
  have e1 : (![5, 8, 4, 1] : Fin 4 → ℝ) 1 - (![5, 8, 4, 1] : Fin 4 → ℝ) 0 - ε = 3 - ε := by
    simp; norm_num
  have e0 : (![5, 8, 4, 1] : Fin 4 → ℝ) 0 + ε = 5 + ε := by simp
  rw [e1, e0] at hgap
  exact hgap

theorem Tz_le_sharpRHS {σ γ R C : ℝ} {k : ℕ} (hσ : 1 < σ) (hσ2 : σ ≤ 2) (hR : 0 < R)
    (hR1 : R < 1) (haR : σ - R < 1) (ht : 3 ≤ |(k : ℝ) * γ|)
    (hC : sharpRHS σ ((k : ℝ) * γ) R ≤ C) : TrigPolyRegion.Tz σ γ k ≤ C :=
  le_trans (GrowthSharpClean.growth_at_sharp_clean (t := (k : ℝ) * γ) hσ hσ2 hR hR1 haR ht) hC

theorem zeta_region_5841_sharp :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ {γ β ε B C₂ C₃ σ R R' Mζ Ms Rg : ℝ},
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
      σ ≤ 2 → 0 < Rg → Rg < 1 → σ - Rg < 1 → 3 ≤ |2 * γ| →
      sharpRHS σ (2 * γ) Rg ≤ C₂ → sharpRHS σ (3 * γ) Rg ≤ C₃ →
      β ≤ 1 - ε * (3 - ε) / ((5 + ε) * (5 * A₀ + 8 * B + 4 * C₂ + C₃)) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841_Tz
  refine ⟨A₀, hA0, ?_⟩
  intro γ β ε B C₂ C₃ σ R R' Mζ Ms Rg hε hB hC2 hC3 hεK hσdef hβσ hR hRR' hball hρ₀ball hζρ₀
    hp hζsphere hsl hMζ hBexp hσ2 hRg hRg1 haRg h2γ hG2 hG3
  have hKpos : 0 < 5 * A₀ + 8 * B + 4 * C₂ + C₃ := lt_of_lt_of_le hε hεK
  have hσ : 1 < σ := by
    rw [hσdef]; have : 0 < ε / (5 * A₀ + 8 * B + 4 * C₂ + C₃) := div_pos hε hKpos; linarith
  have e2 : ((2 : ℕ) : ℝ) * γ = 2 * γ := by norm_num
  have e3 : ((3 : ℕ) : ℝ) * γ = 3 * γ := by norm_num
  have h3γ : 3 ≤ |3 * γ| := by
    refine le_trans h2γ ?_
    rw [abs_mul, abs_mul]
    exact mul_le_mul_of_nonneg_right (by norm_num) (abs_nonneg γ)
  have hT2 : TrigPolyRegion.Tz σ γ 2 ≤ C₂ :=
    Tz_le_sharpRHS hσ hσ2 hRg hRg1 haRg (by rw [e2]; exact h2γ) (by rw [e2]; exact hG2)
  have hT3 : TrigPolyRegion.Tz σ γ 3 ≤ C₃ :=
    Tz_le_sharpRHS hσ hσ2 hRg hRg1 haRg (by rw [e3]; exact h3γ) (by rw [e3]; exact hG3)
  exact H hε hB hC2 hC3 hεK hσdef hβσ hR hRR' hball hρ₀ball hζρ₀ hp hζsphere hsl hMζ hBexp hT2 hT3

end TrigPoly5841Sharp

#print axioms TrigPoly5841Sharp.zeta_region_5841_Tz
#print axioms TrigPoly5841Sharp.Tz_le_sharpRHS
#print axioms TrigPoly5841Sharp.zeta_region_5841_sharp
