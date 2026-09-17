import RegionTwoGrowthPole

/-!
# §1433 ★★ — case-split removing the zero-geometry hypothesis `hp`

Either the zero `β+iγ` is CLOSE to the test centre (`σK-β ≤ R/2`, i.e. `hp` holds → apply the capstone §1432), or
it is FAR (`σK-β > R/2`), in which case `β < σK - R/2 ≤ 1-1/(14K)` directly from the scalar threshold `hfar`.
Either way `β ≤ 1-1/(14K)`.  This removes `hp` from the interface, replacing it with the pure scalar condition
`hfar : σK - R/2 ≤ 1-1/(14K)` (which holds once `K` is large enough, `K ≥ ~5`).  EXIT0 / std-axioms.
-/

open Complex Metric DirichletCharacter

theorem region_riemannZeta_or_far {γ β A B C R R' R'' R₂ R'₂ R''₂ : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hfar : (1 + 1 / (2 * (3 * A + 4 * B + C))) - R / 2 ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)))
    (hball : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hR''0 : 0 < R'') (hR''1 : R'' < 1) (hR''t : R'' < |γ|)
    (hincl : (1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 + R < R'')
    (hBexplicit : ∀ Ms : ℝ, (∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R,
        ‖riemannZeta z‖ ≤ Ms) →
      4 * Real.log ((Ms + 1) / ‖riemannZeta
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    (hR₂ : 0 < R₂) (hRR'₂ : R₂ < R'₂)
    (hball₂ : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)) R'₂)
    (hR''₂0 : 0 < R''₂) (hR''₂1 : R''₂ < 1) (hR''₂t : R''₂ < |2 * γ|)
    (hincl₂ : (1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 + R₂ < R''₂)
    (hBexplicit₂ : ∀ Ms₂ : ℝ, (∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
        + I * ((2 * γ : ℝ) : ℂ)) R₂, ‖riemannZeta z‖ ≤ Ms₂) →
      4 * Real.log ((Ms₂ + 1) / ‖riemannZeta
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))‖) / (R₂ / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  by_cases hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2)
  · have hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
        ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R' :=
      mem_ball.mpr (lt_of_le_of_lt (mem_closedBall.mp hp) (by linarith))
    exact region_riemannZeta_scalar hR hRR' hK hA hB hC hβ hball hρ₀ball hζρ₀ hp
      hR''0 hR''1 hR''t hincl hBexplicit hA_pole hR₂ hRR'₂ hball₂ hR''₂0 hR''₂1 hR''₂t hincl₂ hBexplicit₂
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

#print axioms region_riemannZeta_or_far
