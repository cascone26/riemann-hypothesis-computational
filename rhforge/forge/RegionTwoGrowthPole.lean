import RegionTwoGrowth
import HpoleRelaxed

/-!
# §1432 ★★★ — quantitative ζ region reduced to PURE scalar/geometry data (`hpole` discharged)

Discharges §1431's `hpole` via §1427 (`hpole_relaxed`): any `A` above the pole constant
`‖logDeriv (LFunctionTrivChar₁ 1) σK‖` gives the pole value bound.  The quantitative de la Vallée Poussin
zero-free region `β ≤ 1-1/(14K)` for real ζ now follows from an entirely **scalar + geometry** hypothesis set — no
`g₀`, no `Mζ` existential, no complex-analytic object anywhere in the interface:

* constants `hK/hA/hB/hC/hβ`; `hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1) σK‖ ≤ A`;
* the kept-zero data `hζρ₀`, `hρ₀ball`, `hp`; the two disks' geometry; and the two explicit BC scalars.

The ONLY thing between this and an unconditional region is the scalar constant-choice fixed-point.  EXIT0 / std.
-/

open Complex Metric DirichletCharacter

theorem region_riemannZeta_scalar {γ β A B C R R' R'' R₂ R'₂ R''₂ : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hball : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2))
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
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  exact region_riemannZeta_from_two_growth hR hRR' hK hA hB hC hβ hball hρ₀ball hζρ₀ hp
    hR''0 hR''1 hR''t hincl hBexplicit (hpole_relaxed hσK hA_pole)
    hR₂ hRR'₂ hball₂ hR''₂0 hR''₂1 hR''₂t hincl₂ hBexplicit₂

#print axioms region_riemannZeta_scalar
