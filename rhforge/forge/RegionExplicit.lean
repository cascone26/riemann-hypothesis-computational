import RegionRiemannZetaFromSphere
import ReNegLogDerivLeBOfSphere
import HpoleRelaxed

/-!
# §1435 ★★★ — quantitative ζ region from EXPLICIT sphere bounds (corrected, satisfiable interface)

Rebuilds the two-disk region on §1425 (`region_riemannZeta_from_sphere`, explicit `Ms`) — NOT the flawed ∀Ms
discharge of §1426.  `hgrowth` is supplied by §1429 at the growth centre `σK+2iγ` with an explicit `Ms₂`, and
`hpole` by §1427 (`A` free).  Every hypothesis is now genuinely satisfiable: explicit sphere bounds `Ms`, `Ms₂` (to
be supplied by an explicit ζ growth estimate), the two explicit BC scalars, and the pole constant.  Conclusion:
`β ≤ 1-1/(14K)`.  EXIT0 / std-axioms.
-/

open Complex Metric DirichletCharacter

theorem region_riemannZeta_explicit {γ β A B C R R' Mζ Ms R₂ R'₂ Mζ₂ Ms₂ : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hball : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2))
    (hζsphere : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R,
      ‖riemannZeta z‖ ≤ Ms)
    (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 4 * Real.log (Mζ / ‖riemannZeta
      (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    -- growth disk (centre σK+2iγ), explicit Ms₂:
    (hR₂ : 0 < R₂) (hRR'₂ : R₂ < R'₂)
    (hball₂ : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)) R'₂)
    (hζsphere₂ : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
        + I * ((2 * γ : ℝ) : ℂ)) R₂, ‖riemannZeta z‖ ≤ Ms₂)
    (hsl₂ : Ms₂ < Mζ₂) (hMζ₂ : 0 < Mζ₂)
    (hBexplicit₂ : 4 * Real.log (Mζ₂ / ‖riemannZeta
      (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))‖) / (R₂ / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  have hcre₂ : 1 < (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
    linarith
  have hgrowth := re_neg_logDeriv_le_B_of_sphere hR₂ hRR'₂ hcre₂ hball₂ hζsphere₂ hsl₂ hMζ₂ hBexplicit₂
  exact region_riemannZeta_from_sphere hR hRR' hK hA hB hC hβ hball hρ₀ball hζρ₀ hp
    hζsphere hsl hMζ hBexplicit (hpole_relaxed hσK hA_pole) hgrowth

#print axioms region_riemannZeta_explicit
