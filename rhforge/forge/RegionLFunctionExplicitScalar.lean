import RegionLFunctionDischarge
import HByscalarOfExplicitLFunction

set_option maxHeartbeats 1000000

/-!
# §1543 — quantitative L(·,χ) region from CLEAN explicit scalar bounds (no norms)

L-analogue of ζ §1453.  Discharges §1541's two BC conditions via §1542, so the region `β ≤ 1-1/(14K)` for
`LFunction χ` (χ≠1, χ²≠1) now follows from the norm-free explicit inequalities

  `B ≥ 4·log((2+N+(σK+|γ|+R)·N/(σK-R))·(3+1/(σK-1)))/(R/2)`     (χ, at iγ)
  `C ≥ 4·log((2+N+(σK+|2γ|+R)·N/(σK-R))·(3+1/(σK-1)))/(R/2)`    (χ², at 2iγ)

together with `LFunction χ (β+iγ)=0`, the constant/geometry scalars, and `hA_pole`.  `σK = 1+1/(2K)`,
`K = 3A+4B+C`, `N` = character modulus.  Fully explicit scalar form — last brick before the numeric constant
choice.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem region_LFunction_explicit_scalar {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2 : χ ^ 2 ≠ 1) {γ β A B C R R' : ℝ}
    (hR : 0 < R) (hR1 : R < 1) (hRR' : R < R')
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C) (hK : 0 < 3 * A + 4 * B + C)
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2))
    (hBexp : 4 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |γ| + R) * (N : ℝ)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - R))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (R / 2) ≤ B)
    (hCexp : 4 * Real.log ((2 + (N : ℝ)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |2 * γ| + R) * (N : ℝ)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - R))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (R / 2) ≤ C)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  have hRσK : R < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by linarith
  exact region_LFunction_discharge hχ hχ2 hR hR1 hRR' hA hB hC hK hρ₀ball hLρ₀ hp
    (hByscalar_of_explicit_LFunction hσK hR hRσK hBexp)
    (hByscalar_of_explicit_LFunction (γ := 2 * γ) hσK hR hRσK hCexp)
    hA_pole

#print axioms region_LFunction_explicit_scalar
