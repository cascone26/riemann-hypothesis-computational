import RegionLFunctionFromSphere
import LFunctionSphereN1
import HBexplicitOfGeLFunction
import HgrowthLFunctionScalar
import HpoleRiemannZeta

set_option maxHeartbeats 1000000

/-!
# §1541 — quantitative L(·,χ) region from PURE SCALAR data (all analytic hyps discharged)

L-analogue of ζ §1448.  Threads §1509 (`region_LFunction_from_sphere`) at `c = σK+iγ`, `R`, discharging EVERY
complex-analytic hypothesis by a scalar dischargER:
  * `hζsphere` via §1535 (`Ms := 1 + (N·2^{-(c.re-R)} + (‖c‖+R)·N/(c.re-R))`, `Mζ := Ms+1`),
  * `hBexplicit` via §1539 (BC scalar from the L lower bound §1538) — reduces to `hByscalar`,
  * `hgrowth` via §1540 at `χ²`, `c₂ = σK+2iγ` (needs `χ²≠1`, the non-quadratic case) — reduces to `hCyscalar`,
  * `hpole` via §1405 (ζ pole) — reduces to `hA_pole` (the ζ pole constant, §1360).
Leaves only PURE SCALAR inequalities `hByscalar`/`hCyscalar`/`hA_pole` + the zero-geometry `hp`/`hρ₀ball`.
Result: `β ≤ 1 - 1/(14K)`, `K = 3A+4B+C`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric

theorem region_LFunction_discharge {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    (hχ2 : χ ^ 2 ≠ 1) {γ β A B C R R' : ℝ} (hR : 0 < R) (hR1 : R < 1) (hRR' : R < R')
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
    (hCyscalar : 4 * Real.log ((1 + ((N : ℝ)
          * (2 : ℝ) ^ (-((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
              + I * ((2 * γ : ℝ) : ℂ)).re - R))
          + (‖(((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))‖ + R) * (N : ℝ)
            / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re - R)) + 1)
        * (3 + 1 / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
            + I * ((2 * γ : ℝ) : ℂ)).re - 1)))
        / (R / 2) ≤ C)
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
  have hc2R : 0 < c₂.re - R := by rw [hc2re]; linarith
  -- centre-χ scalar constants
  set Ms : ℝ := 1 + ((N : ℝ) * (2 : ℝ) ^ (-(c.re - R)) + (‖c‖ + R) * (N : ℝ) / (c.re - R)) with hMs
  have hMspos : 0 < Ms + 1 := by
    have ha : (0 : ℝ) ≤ ‖c‖ + R := add_nonneg (norm_nonneg c) hR.le
    have hb : (0 : ℝ) ≤ (N : ℝ) := Nat.cast_nonneg N
    have hnn : 0 ≤ (‖c‖ + R) * (N : ℝ) / (c.re - R) := div_nonneg (mul_nonneg ha hb) hcR.le
    have hnn2 : 0 ≤ (N : ℝ) * (2 : ℝ) ^ (-(c.re - R)) := mul_nonneg hb (Real.rpow_nonneg (by norm_num) _)
    rw [hMs]; linarith
  refine region_LFunction_from_sphere (Mζ := Ms + 1) (Ms := Ms) hχ hR hRR' hA hB hC hK
    hρ₀ball hLρ₀ hp ?_ (by linarith) hMspos ?_ ?_ ?_
  · -- hζsphere via §1535
    intro z hz
    exact norm_LFunction_le_sphere_N1 hχ hR hcR hz
  · -- hBexplicit via §1539
    exact hBexplicit_of_ge_LFunction hχ hcre1 hMspos hR hByscalar
  · -- hpole via §1405 + hA_pole
    exact le_trans (hpole_riemannZeta hσK) (by linarith [hA_pole])
  · -- hgrowth via §1540 at χ², c₂
    exact hgrowth_LFunction_scalar hχ2 hR hc2re1 hc2R hCyscalar

#print axioms region_LFunction_discharge
