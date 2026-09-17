import RegionExplicit
import HfarHolds

/-!
# §1436 ★★ — corrected region with the zero-geometry case-split (`hp`/`hρ₀ball` removed)

Same case-split as §1433 but over the CORRECTED, satisfiable §1435: CLOSE (`σK-β ≤ R/2`) → §1435 (`hρ₀ball` from
`hp`); FAR (`σK-β > R/2`) → `β < σK - R/2 ≤ 1-1/(14K)` via the scalar `hfar` (discharged by §1434 for `8 ≤ 7KR`).
Every remaining hypothesis is scalar/geometry or an explicit sphere bound — all satisfiable.  EXIT0 / std-axioms.
-/

open Complex Metric DirichletCharacter

theorem region_riemannZeta_explicit_or_far {γ β A B C R R' Mζ Ms R₂ R'₂ Mζ₂ Ms₂ : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hfar : (1 + 1 / (2 * (3 * A + 4 * B + C))) - R / 2 ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)))
    (hball : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hζsphere : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R,
      ‖riemannZeta z‖ ≤ Ms)
    (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 4 * Real.log (Mζ / ‖riemannZeta
      (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    (hR₂ : 0 < R₂) (hRR'₂ : R₂ < R'₂)
    (hball₂ : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)) R'₂)
    (hζsphere₂ : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
        + I * ((2 * γ : ℝ) : ℂ)) R₂, ‖riemannZeta z‖ ≤ Ms₂)
    (hsl₂ : Ms₂ < Mζ₂) (hMζ₂ : 0 < Mζ₂)
    (hBexplicit₂ : 4 * Real.log (Mζ₂ / ‖riemannZeta
      (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))‖) / (R₂ / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  by_cases hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2)
  · have hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
        ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R' :=
      mem_ball.mpr (lt_of_le_of_lt (mem_closedBall.mp hp) (by linarith))
    exact region_riemannZeta_explicit hR hRR' hK hA hB hC hβ hball hρ₀ball hζρ₀ hp
      hζsphere hsl hMζ hBexplicit hA_pole hR₂ hRR'₂ hball₂ hζsphere₂ hsl₂ hMζ₂ hBexplicit₂
  · rw [mem_closedBall, not_le] at hp
    have hdist : dist ((β : ℂ) + I * (γ : ℝ))
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))
        = |β - (1 + 1 / (2 * (3 * A + 4 * B + C)))| := by
      rw [Complex.dist_eq]
      have hsub : ((β : ℂ) + I * (γ : ℝ)) - (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))
          = (((β - (1 + 1 / (2 * (3 * A + 4 * B + C))) : ℝ)) : ℂ) := by push_cast; ring
      rw [hsub, Complex.norm_real, Real.norm_eq_abs]
    rw [hdist] at hp
    have hβlt : β - (1 + 1 / (2 * (3 * A + 4 * B + C))) < 0 := by linarith
    rw [abs_of_neg hβlt] at hp
    linarith

#print axioms region_riemannZeta_explicit_or_far
