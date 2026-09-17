import HdecompHgBRelaxedOfSphereLFunction
import RegionLFunctionConnect
import BetaLtSigmaKLFunction
import Mathlib.NumberTheory.LSeries.Nonvanishing

/-!
# §1509 ★★★ — quantitative L(·,χ) region from `g₀`-free sphere data (only the scalar fixed-point remains)

L-analogue of ζ §1425.  Specializes the CRUX §1508 at the DLVP centre `c = σK+iγ` (`σK = 1+1/(2K)`,
`K = 3A+4B+C`) and threads it through the connector §1498.  Every hypothesis is now `g₀`-free — sphere bound `Ms`,
slack `Ms < Mζ`, the explicit BC scalar `4·log(Mζ/‖L c‖)/(R/2) ≤ B`, the pole/growth value bounds, and the
kept-zero geometry.  `hβ`, `hz₀`, `hz₀ne` are auto-discharged.  Conclusion: `β ≤ 1 - 1/(14K)`.  The ONLY thing
between this and an unconditional statement is the scalar constant-choice fixed-point.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem region_LFunction_from_sphere {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {γ β A B C R R' Mζ Ms : ℝ} (hR : 0 < R) (hRR' : R < R')
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C) (hK : 0 < 3 * A + 4 * B + C)
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2))
    (hζsphere : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R,
      ‖LFunction χ z‖ ≤ Ms)
    (hsl : Ms < Mζ) (hMζ : 0 < Mζ)
    (hBexplicit : 4 * Real.log (Mζ / ‖LFunction χ
      (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B)
    (hpole : (-logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ))).re
      ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : (-logDeriv (LFunction (χ ^ 2))
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))).re ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hcre : 1 < (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)).re := by
    have hpos : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
    linarith
  obtain ⟨g, hdec, hgB⟩ :=
    hdecomp_and_hgB_relaxed_of_sphere_LFunction hχ hR hRR' hcre hζsphere hsl hMζ hBexplicit
  have hβ := beta_lt_sigmaK_of_LFunction_zero hχ hK hLρ₀
  have hz₀ne : LFunction χ (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) ≠ 0 :=
    LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) hcre.le
  exact region_LFunction_connect hχ hK hA hB hC hβ hρ₀ball hLρ₀
    (mem_ball_self (lt_trans hR hRR')) hz₀ne hp hgB hdec hpole hgrowth

#print axioms region_LFunction_from_sphere
