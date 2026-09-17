import RegionRiemannZetaFromSphere
import ZetaBoundSphereCenter

/-!
# §1426 ★★ — quantitative ζ region with the sphere bound discharged (§1398)

Discharges §1425's sphere-bound trio (`hζsphere`, `hsl`, `hMζ`) using §1398 `exists_norm_riemannZeta_le_sphere`:
obtain some `Ms` bounding `‖ζ‖` on `sphere(σK+iγ) R`, take `Mζ := Ms + 1` (so `Ms < Mζ` and `0 < Mζ` since a
sphere point gives `Ms ≥ 0`).  The region `β ≤ 1 - 1/(14K)` now follows from: the constant inequalities, the
factorization geometry (`R < R'`, `1 ∉ ball c R'`), the §1398 growth geometry (`R'' < 1`, `R'' < |γ|`,
`σK-1+R < R''`), the explicit BC scalar `hBexplicit`, the kept-zero geometry, and the pole/growth value bounds.
The interface is now free of any `Mζ` existential.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem region_riemannZeta_from_growth {γ β A B C R R' R'' : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
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
    (hpole : (-logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ))).re
      ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : (-logDeriv riemannZeta
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))).re ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  obtain ⟨Ms, hMs⟩ := exists_norm_riemannZeta_le_sphere (R := R) hσK hR''0 hR''1 hR''t hincl
  have hs0mem : (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) + (R : ℝ)
      ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R := by
    simp only [mem_sphere_iff_norm, add_sub_cancel_left, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos hR]
  have hMs0 : 0 ≤ Ms := (norm_nonneg _).trans (hMs _ hs0mem)
  exact region_riemannZeta_from_sphere hR hRR' hK hA hB hC hβ hball hρ₀ball hζρ₀ hp
    hMs (by linarith) (by linarith) (hBexplicit Ms hMs) hpole hgrowth

#print axioms region_riemannZeta_from_growth
