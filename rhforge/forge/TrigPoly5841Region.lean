import HzeroAtPoint
import GrowthAtKGamma
import TrigPolyTermNonneg

/-!
# S2936 — the 5-8-4-1 zero-free region as ONE composite theorem (lane 1)

P = 5 + 8cosθ + 4cos2θ + cos3θ = (1+cosθ)(1+2cosθ)² ≥ 0 (S2921 `p5841_nonneg`) plugged into the
pole-discharged P-generic region (S2932 `zeta_region_P_pole`), with the ZERO input supplied by
S2934 `hzero_Tz` (BC sphere data at σ+iγ) and the two GROWTH inputs (2iγ, 3iγ) supplied by
S2928 `growth_Tz`.  There is one absolute A₀ ≥ 0 such that, with 𝒦 = 5A₀ + 8B + 4C₂ + C₃ and
σ = 1 + ε/𝒦, every ζ zero β+iγ (β < σ) satisfying the sphere data obeys
  β ≤ 1 − ε(3 − ε)/((5 + ε)𝒦).
No abstract `hpole/hzero/hgrowth` hypothesis remains — only the explicit sphere data and the scalar
inequalities growthRHS(σ, kγ, R) ≤ C_k (to be discharged by S2930 log absorption).
-/

open Complex Metric

namespace TrigPoly5841

/-- The explicit right-hand side of S2928 `growth_at` at σ + it. -/
noncomputable def growthRHS (σ t R : ℝ) : ℝ :=
  4 * Real.log ((2 + 1 / (‖(σ : ℂ) + I * (t : ℂ) - 1‖ - R)
          + (‖(σ : ℂ) + I * (t : ℂ)‖ + R) / (((σ : ℂ) + I * (t : ℂ)).re - R) + 1)
          * (3 + 1 / (((σ : ℂ) + I * (t : ℂ)).re - 1))) / (R / 2)

theorem growth_Tz' {σ γ R R' : ℝ} (k : ℕ) (hσ : 1 < σ) (hR : 0 < R) (hR1 : R < 1)
    (hRR' : R < R') (hR't : R' ≤ |(k : ℝ) * γ|) :
    TrigPolyRegion.Tz σ γ k ≤ growthRHS σ ((k : ℝ) * γ) R :=
  GrowthK.growth_Tz k hσ hR hR1 hRR' hR't

theorem zeta_region_5841 :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ {γ β ε B C₂ C₃ σ R R' Mζ Ms Rg Rg' : ℝ},
      0 < ε → 0 ≤ B → 0 ≤ C₂ → 0 ≤ C₃ →
      ε ≤ 5 * A₀ + 8 * B + 4 * C₂ + C₃ →
      σ = 1 + ε / (5 * A₀ + 8 * B + 4 * C₂ + C₃) →
      β < σ →
      -- zero input: BC sphere data at σ + iγ
      0 < R → R < R' →
      (1 : ℂ) ∉ ball ((σ : ℂ) + I * (γ : ℝ)) R' →
      ((β : ℂ) + I * (γ : ℝ)) ∈ ball ((σ : ℂ) + I * (γ : ℝ)) R' →
      riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0 →
      ((β : ℂ) + I * (γ : ℝ)) ∈ closedBall ((σ : ℂ) + I * (γ : ℝ)) (R / 2) →
      (∀ z ∈ sphere ((σ : ℂ) + I * (γ : ℝ)) R, ‖riemannZeta z‖ ≤ Ms) →
      Ms < Mζ → 0 < Mζ →
      4 * Real.log (Mζ / ‖riemannZeta ((σ : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B →
      -- growth input at 2iγ and 3iγ
      0 < Rg → Rg < 1 → Rg < Rg' → Rg' ≤ |2 * γ| →
      growthRHS σ (2 * γ) Rg ≤ C₂ → growthRHS σ (3 * γ) Rg ≤ C₃ →
      β ≤ 1 - ε * (3 - ε) / ((5 + ε) * (5 * A₀ + 8 * B + 4 * C₂ + C₃)) := by
  obtain ⟨A₀, hA0, H⟩ := TrigPolyRegion.zeta_region_P_pole
  refine ⟨A₀, hA0, ?_⟩
  intro γ β ε B C₂ C₃ σ R R' Mζ Ms Rg Rg' hε hB hC2 hC3 hεK hσdef hβσ hR hRR' hball hρ₀ball hζρ₀
    hp hζsphere hsl hMζ hBexp hRg hRg1 hRgRg' hRg't hG2 hG3
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
      · have h := growth_Tz' (γ := γ) 2 hσ hRg hRg1 hRgRg' (by exact_mod_cast hRg't)
        simp only [Matrix.cons_val_zero]
        have e : ((2 : ℕ) : ℝ) * γ = 2 * γ := by norm_num
        rw [e] at h
        exact le_trans (by simpa using h) hG2
      · have h3 : Rg' ≤ |((3 : ℕ) : ℝ) * γ| := by
          refine le_trans hRg't ?_
          rw [abs_mul, abs_mul]; norm_num
          exact mul_le_mul_of_nonneg_right (by norm_num) (abs_nonneg γ)
        have h := growth_Tz' (γ := γ) 3 hσ hRg hRg1 hRgRg' h3
        have e : ((3 : ℕ) : ℝ) * γ = 3 * γ := by norm_num
        rw [e] at h
        simp only [Matrix.cons_val_one, Matrix.head_cons]
        exact le_trans (by simpa using h) hG3)
  rw [h𝒦] at hgap
  have e1 : (![5, 8, 4, 1] : Fin 4 → ℝ) 1 - (![5, 8, 4, 1] : Fin 4 → ℝ) 0 - ε = 3 - ε := by
    simp; norm_num
  have e0 : (![5, 8, 4, 1] : Fin 4 → ℝ) 0 + ε = 5 + ε := by simp
  rw [e1, e0] at hgap
  exact hgap

end TrigPoly5841

#print axioms TrigPoly5841.zeta_region_5841
