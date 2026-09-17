import RegionRiemannZetaFromGrowth
import ReNegLogDerivLeBOfSphere

/-!
# §1430 ★★ — quantitative ζ region with BOTH disks' data (`hgrowth` discharged via §1429)

Builds §1426's `hgrowth` by applying §1429 (`re_neg_logDeriv_le_B_of_sphere`) at the growth centre `σK+2iγ`
(a point with `Re = σK > 1`), so the growth value bound `(-logDeriv ζ(σK+2iγ)).re ≤ C` becomes the explicit BC
scalar `hBexplicit₂` at the `2iγ` disk.  The region `β ≤ 1-1/(14K)` now follows from purely `g₀`-free data on TWO
Borel–Carathéodory disks (`iγ` for `hgB`, `2iγ` for `hgrowth`) plus the pole bound and geometry — the last
non-scalar input (`hgrowth`) is gone.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem region_riemannZeta_from_growth_disks {γ β A B C R R' R'' R₂ R'₂ Mζ₂ Ms₂ : ℝ}
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
    (hpole : (-logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ))).re
      ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    -- growth disk (centre σK+2iγ):
    (hR₂ : 0 < R₂) (hRR'₂ : R₂ < R'₂)
    (hball₂ : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)) R'₂)
    (hζsphere₂ : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
        + I * ((2 * γ : ℝ) : ℂ)) R₂, ‖riemannZeta z‖ ≤ Ms₂)
    (hsl₂ : Ms₂ < Mζ₂) (hMζ₂ : 0 < Mζ₂)
    (hBexplicit₂ : 4 * Real.log (Mζ₂ / ‖riemannZeta
      (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))‖) / (R₂ / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hcre₂ : 1 < (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re := by
    have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
      have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
      linarith
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
    linarith
  have hgrowth := re_neg_logDeriv_le_B_of_sphere hR₂ hRR'₂ hcre₂ hball₂ hζsphere₂ hsl₂ hMζ₂ hBexplicit₂
  exact region_riemannZeta_from_growth hR hRR' hK hA hB hC hβ hball hρ₀ball hζρ₀ hp
    hR''0 hR''1 hR''t hincl hBexplicit hpole hgrowth

#print axioms region_riemannZeta_from_growth_disks
