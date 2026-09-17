import HdecompHgBRelaxedOfSphere
import RegionRiemannZetaConnect
import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1425 ★★★ — quantitative ζ region from `g₀`-free data (only the scalar fixed-point remains)

Specializes the crux §1424 at the DLVP test centre `c = σK+iγ` (`σK = 1+1/(2K)`, `K = 3A+4B+C`) and threads the
result through the connector §1409.  Every hypothesis is now `g₀`-free — sphere bound `Ms`, slack `Ms < Mζ`, the
explicit BC scalar `4·log(Mζ/‖ζ c‖)/(R/2) ≤ B`, the pole and growth value bounds, the kept-zero geometry, and the
constant inequalities.  Conclusion: `β ≤ 1 - 1/(14K)`.  This is the quantitative de la Vallée Poussin zero-free
region for ζ with the factorization fully internalized — the ONLY thing between this and an unconditional statement
is the scalar constant-choice fixed-point (choosing `A,B,C,R,R',Mζ,Ms,γ` self-consistently).  EXIT0 / std-axioms.
-/

open Complex Metric

theorem region_riemannZeta_from_sphere {γ β A B C R R' Mζ Ms : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
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
    (hpole : (-logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ))).re
      ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : (-logDeriv riemannZeta
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))).re ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  have hcre : 1 < (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, zero_sub, neg_zero, add_zero]
    linarith
  obtain ⟨g, hdec, hgB⟩ :=
    hdecomp_and_hgB_relaxed_of_sphere hR hRR' hcre hball hζsphere hsl hMζ hBexplicit
  exact region_riemannZeta_connect hK hA hB hC hβ hball hρ₀ball hζρ₀
    (mem_ball_self (lt_trans hR hRR')) (riemannZeta_ne_zero_of_one_lt_re hcre) hp hgB hdec hpole hgrowth

#print axioms region_riemannZeta_from_sphere
