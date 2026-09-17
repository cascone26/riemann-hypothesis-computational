import RegionLFunctionFromSphere
import LFunctionSphereN1
import HBexplicitOfGeLFunction
import HgrowthLFunctionQuad
import HpoleRiemannZeta

set_option maxHeartbeats 1000000

/-!
# §1557 — quantitative L(·,χ) region for the QUADRATIC case (`χ²=1`) from scalar data

Quadratic variant of §1541.  The χ side (χ≠1, entire) is IDENTICAL: `hζsphere` via §1535, `hBexplicit` via §1539,
`hpole` via §1405.  The only change is the `hgrowth` (T₂ = `LFunction χ²` at `c₂ = σK+2iγ`): since `χ²=1`,
`LFunction χ²` has a ζ-like pole, so it is discharged via §1556 `hgrowth_LFunction_quad` (ζ sphere bound `Ms₂` +
`hball₂` pole-avoidance + ζ-BC scalar `hCyscalar` + the finite-product correction `∑_{p|N} log p`) instead of §1540.
Result: `β ≤ 1 - 1/(14K)` for a zero of `LFunction χ` with `χ²=1`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem region_LFunction_discharge_quad {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    (hχ2eq : χ ^ 2 = 1) {γ β A B C R R' Ms₂ : ℝ} (hR : 0 < R) (hR1 : R < 1) (hRR' : R < R')
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C) (hK : 0 < 3 * A + 4 * B + C)
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2))
    (hByscalar : 4 * Real.log ((1 + ((N : ℝ)
          * (2 : ℝ) ^ (-((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)).re - R))
          + (‖(((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))‖ + R) * (N : ℝ)
            / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)).re - R)) + 1)
        * (3 + 1 / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)).re - 1)))
        / (R / 2) ≤ B)
    (hball₂ : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
      + I * ((2 * γ : ℝ) : ℂ)) R')
    (hζsphere₂ : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
        + I * ((2 * γ : ℝ) : ℂ)) R, ‖riemannZeta z‖ ≤ Ms₂)
    (hCyscalar : 4 * Real.log ((Ms₂ + 1)
        * (3 + 1 / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
            + I * ((2 * γ : ℝ) : ℂ)).re - 1))) / (R / 2)
      ≤ C - ∑ p ∈ N.primeFactors, Real.log p)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  set c : ℂ := ((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ) with hc_def
  set c₂ : ℂ := ((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ) with hc2_def
  have hcre : c.re = 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    rw [hc_def]; simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re,
      Complex.I_im, Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  have hc2re : c₂.re = 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    rw [hc2_def]; simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re,
      Complex.I_im, Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  have hcre1 : 1 < c.re := by rw [hcre]; exact hσK
  have hc2re1 : 1 < c₂.re := by rw [hc2re]; exact hσK
  have hcR : 0 < c.re - R := by rw [hcre]; linarith
  set Ms : ℝ := 1 + ((N : ℝ) * (2 : ℝ) ^ (-(c.re - R)) + (‖c‖ + R) * (N : ℝ) / (c.re - R)) with hMs
  have hMspos : 0 < Ms + 1 := by
    have ha : (0 : ℝ) ≤ ‖c‖ + R := add_nonneg (norm_nonneg c) hR.le
    have hb : (0 : ℝ) ≤ (N : ℝ) := Nat.cast_nonneg N
    have hnn : 0 ≤ (‖c‖ + R) * (N : ℝ) / (c.re - R) := div_nonneg (mul_nonneg ha hb) hcR.le
    have hnn2 : 0 ≤ (N : ℝ) * (2 : ℝ) ^ (-(c.re - R)) := mul_nonneg hb (Real.rpow_nonneg (by norm_num) _)
    rw [hMs]; linarith
  refine region_LFunction_from_sphere (Mζ := Ms + 1) (Ms := Ms) hχ hR hRR' hA hB hC hK
    hρ₀ball hLρ₀ hp ?_ (by linarith) hMspos ?_ ?_ ?_
  · intro z hz
    exact norm_LFunction_le_sphere_N1 hχ hR hcR hz
  · exact hBexplicit_of_ge_LFunction hχ hcre1 hMspos hR hByscalar
  · exact le_trans (hpole_riemannZeta hσK) (by linarith [hA_pole])
  · -- hgrowth via §1556 (quadratic χ²=1)
    have hq := hgrowth_LFunction_quad hχ2eq hR hRR' hc2re1 hball₂ hζsphere₂ hCyscalar
    linarith [hq]

#print axioms region_LFunction_discharge_quad
