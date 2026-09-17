import RegionExplicitOrFar
import PoleAvoidGeom

/-!
# §1437 ★★ — corrected region with pole-avoidance geometry discharged (`hball`/`hball₂` removed)

Discharges §1436's `hball` and `hball₂` (`1 ∉ ball (σK+iγ) R'`, `1 ∉ ball (σK+2iγ) R'₂`) via §1371
(`one_notMem_ball_test`): since `|c-1| = √((σK-1)²+γ²) ≥ |γ|`, any radius `≤ |γ|` avoids the pole.  Replaces the two
`ball` non-membership hypotheses with the scalar radius bounds `R' ≤ |γ|`, `R'₂ ≤ |2γ|`.  The corrected region now
rests on: `ζ(β+iγ)=0` + constants (`hK/hA/hB/hC/hβ/hfar/hA_pole`) + scalar radius bounds + explicit sphere bounds
(`Ms/Ms₂`) + the two BC scalars.  EXIT0 / std-axioms.
-/

open Complex Metric DirichletCharacter

theorem region_riemannZeta_explicit_geom {γ β A B C R R' Mζ Ms R₂ R'₂ Mζ₂ Ms₂ : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hfar : (1 + 1 / (2 * (3 * A + 4 * B + C))) - R / 2 ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)))
    (hR'γ : R' ≤ |γ|)
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hζsphere : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R,
      ‖riemannZeta z‖ ≤ Ms)
    (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 4 * Real.log (Mζ / ‖riemannZeta
      (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    (hR₂ : 0 < R₂) (hRR'₂ : R₂ < R'₂) (hR'₂γ : R'₂ ≤ |2 * γ|)
    (hζsphere₂ : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
        + I * ((2 * γ : ℝ) : ℂ)) R₂, ‖riemannZeta z‖ ≤ Ms₂)
    (hsl₂ : Ms₂ < Mζ₂) (hMζ₂ : 0 < Mζ₂)
    (hBexplicit₂ : 4 * Real.log (Mζ₂ / ‖riemannZeta
      (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))‖) / (R₂ / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hball : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R' :=
    one_notMem_ball_test hR'γ
  have hball₂ : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
      + I * ((2 * γ : ℝ) : ℂ)) R'₂ := by
    have := one_notMem_ball_test (σ := 1 + 1 / (2 * (3 * A + 4 * B + C))) (γ := 2 * γ) (R := R'₂) hR'₂γ
    simpa using this
  exact region_riemannZeta_explicit_or_far hR hRR' hK hA hB hC hβ hfar hball hζρ₀
    hζsphere hsl hMζ hBexplicit hA_pole hR₂ hRR'₂ hball₂ hζsphere₂ hsl₂ hMζ₂ hBexplicit₂

#print axioms region_riemannZeta_explicit_geom
