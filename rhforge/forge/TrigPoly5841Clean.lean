import TrigPoly5841Region
import HByscalarOfExplicit

/-!
# S2938 — 5-8-4-1 region with NORM-FREE growth hypotheses (lane 1)

S2936 left the growth inputs as `growthRHS σ (kγ) Rg ≤ C_k`, a statement about complex norms.  Here §1452
(`hByscalar_of_explicit`, the same norm→scalar bridge the 3-4-1 chain used) turns each into the clean real
inequality
  4·log((3 + 1/(|kγ|−Rg) + (σ+|kγ|+Rg)/(σ−Rg))·(3 + 1/(σ−1)))/(Rg/2) ≤ C_k,   k = 2, 3,
which is the exact shape that S2930's log absorption (`log_sharp_close`, `log_sharp_close_wide`) discharges.
-/

open Complex Metric

namespace TrigPoly5841

/-- Clean scalar growth bound at height t. -/
noncomputable def cleanRHS (σ t R : ℝ) : ℝ :=
  4 * Real.log ((3 + 1 / (|t| - R) + (σ + |t| + R) / (σ - R)) * (3 + 1 / (σ - 1))) / (R / 2)

theorem growthRHS_le_clean {σ t R C : ℝ} (hσ : 1 < σ) (hR : 0 < R) (hR1 : R < 1) (hRt : R < |t|)
    (h : cleanRHS σ t R ≤ C) : growthRHS σ t R ≤ C :=
  hByscalar_of_explicit hσ hR hRt (by linarith) h

theorem zeta_region_5841_clean :
    ∃ A₀ : ℝ, 0 ≤ A₀ ∧ ∀ {γ β ε B C₂ C₃ σ R R' Mζ Ms Rg Rg' : ℝ},
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
      0 < Rg → Rg < 1 → Rg < Rg' → Rg' ≤ |2 * γ| →
      cleanRHS σ (2 * γ) Rg ≤ C₂ → cleanRHS σ (3 * γ) Rg ≤ C₃ →
      β ≤ 1 - ε * (3 - ε) / ((5 + ε) * (5 * A₀ + 8 * B + 4 * C₂ + C₃)) := by
  obtain ⟨A₀, hA0, H⟩ := zeta_region_5841
  refine ⟨A₀, hA0, ?_⟩
  intro γ β ε B C₂ C₃ σ R R' Mζ Ms Rg Rg' hε hB hC2 hC3 hεK hσdef hβσ hR hRR' hball hρ₀ball hζρ₀
    hp hζsphere hsl hMζ hBexp hRg hRg1 hRgRg' hRg't hG2 hG3
  have hKpos : 0 < 5 * A₀ + 8 * B + 4 * C₂ + C₃ := lt_of_lt_of_le hε hεK
  have hσ : 1 < σ := by
    rw [hσdef]; have : 0 < ε / (5 * A₀ + 8 * B + 4 * C₂ + C₃) := div_pos hε hKpos; linarith
  have h2 : Rg < |2 * γ| := lt_of_lt_of_le hRgRg' hRg't
  have h3 : Rg < |3 * γ| := by
    refine lt_of_lt_of_le h2 ?_
    rw [abs_mul, abs_mul]
    exact mul_le_mul_of_nonneg_right (by norm_num) (abs_nonneg γ)
  exact H hε hB hC2 hC3 hεK hσdef hβσ hR hRR' hball hρ₀ball hζρ₀ hp hζsphere hsl hMζ hBexp
    hRg hRg1 hRgRg' hRg't (growthRHS_le_clean hσ hRg hRg1 h2 hG2)
    (growthRHS_le_clean hσ hRg hRg1 h3 hG3)

end TrigPoly5841

#print axioms TrigPoly5841.zeta_region_5841_clean
