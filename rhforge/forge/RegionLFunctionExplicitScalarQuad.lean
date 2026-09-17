import RegionLFunctionDischargeQuad
import HByscalarOfExplicitLFunction
import HByscalarOfExplicit
import ZetaSphereN1
import AbsLeNormCenterSubOne
import PoleAvoidGeom

set_option maxHeartbeats 1000000

/-!
# §1558 — quantitative L(·,χ) region (QUADRATIC `χ²=1`) from CLEAN explicit (σ,γ) scalars

Quadratic analogue of §1543.  Discharges §1557's data: the χ side `hByscalar` via §1542 (character-free), and the
χ² side (ζ-form, `LFunction χ²` has a pole) via the ζ N=1 sphere bound §1447 + the ζ collapse §1452 + `hball₂`
(§1371).  Leaves the norm-free explicit inequalities:

  `B ≥ 4·log((2+N+(σK+|γ|+R)·N/(σK-R))·(3+1/(σK-1)))/(R/2)`                       (χ)
  `C - ∑_{p|N} log p ≥ 4·log((3+1/(|2γ|-R)+(σK+|2γ|+R)/(σK-R))·(3+1/(σK-1)))/(R/2)`  (χ², ζ-form + correction)

EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem region_LFunction_explicit_scalar_quad {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N}
    (hχ : χ ≠ 1) (hχ2eq : χ ^ 2 = 1) {γ β A B C R R' : ℝ}
    (hR : 0 < R) (hR1 : R < 1) (hRR' : R < R') (hRγ2 : R < |2 * γ|) (hR'2γ : R' ≤ |2 * γ|)
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
    (hCexp : 4 * Real.log ((3 + 1 / (|2 * γ| - R)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |2 * γ| + R)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - R))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (R / 2)
      ≤ C - ∑ p ∈ N.primeFactors, Real.log p)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  have hRσK : R < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by linarith
  have hc2re : (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re
      = 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  have hc2R : 0 < (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re - R := by
    rw [hc2re]; linarith
  have hc1₂ : R < ‖(((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)) - 1‖ :=
    lt_of_lt_of_le hRγ2 (abs_le_norm_center_sub_one (σ := 1 + 1 / (2 * (3 * A + 4 * B + C)))
      (γ := 2 * γ))
  refine region_LFunction_discharge_quad
    (Ms₂ := 2 + 1 / (‖(((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)) - 1‖ - R)
      + (‖(((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))‖ + R)
        / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re - R))
    hχ hχ2eq hR hR1 hRR' hA hB hC hK hρ₀ball hLρ₀ hp
    (hByscalar_of_explicit_LFunction hσK hR hRσK hBexp)
    (one_notMem_ball_test hR'2γ)
    (fun z hz => norm_riemannZeta_le_sphere_N1 hR hc2R hc1₂ hz)
    (hByscalar_of_explicit hσK hR hRγ2 hRσK hCexp)
    hA_pole

#print axioms region_LFunction_explicit_scalar_quad
