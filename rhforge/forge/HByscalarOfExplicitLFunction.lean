import NormCenterLe
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# §1542 — discharge §1541's L `hByscalar`/`hCyscalar` from an explicit `(σ,γ,R)` bound

L-analogue of ζ §1452 — but for the L-shaped BC argument `Ms+1 = 2 + N·2^{-(σ-R)} + (‖c‖+R)·N/(σ-R)` (from the
L N=1 sphere bound §1535: no pole term, but the modulus `N` and the `2^{-(σ-R)}` factor).  Bounds it by the clean
explicit `M' := 2 + N + (σ+|γ|+R)·N/(σ-R)` using `2^{-(σ-R)} ≤ 1` (`σ-R>0`) and `‖c‖ ≤ σ+|γ|` (§1450); then by
log-monotonicity §1541's `hByscalar` (and `hCyscalar`, at `γ→2γ`) holds whenever
`B ≥ 4·log(M'·(3+1/(σ-1)))/(R/2)`.  EXIT0 / std-axioms.
-/

open Complex

theorem hByscalar_of_explicit_LFunction {σ γ R : ℝ} {N : ℕ} {B : ℝ} (hσ : 1 < σ) (hR : 0 < R)
    (hRσ : R < σ)
    (hB : 4 * Real.log ((2 + (N : ℝ) + (σ + |γ| + R) * (N : ℝ) / (σ - R)) * (3 + 1 / (σ - 1)))
        / (R / 2) ≤ B) :
    4 * Real.log ((1 + ((N : ℝ) * (2 : ℝ) ^ (-((((σ : ℝ) : ℂ) + I * (γ : ℝ)).re - R))
          + (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) * (N : ℝ) / ((((σ : ℝ) : ℂ) + I * (γ : ℝ)).re - R)) + 1)
        * (3 + 1 / ((((σ : ℝ) : ℂ) + I * (γ : ℝ)).re - 1))) / (R / 2) ≤ B := by
  have hcre : (((σ : ℝ) : ℂ) + I * (γ : ℝ)).re = σ := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  rw [hcre]
  have hσR : (0 : ℝ) < σ - R := by linarith
  have hσ1 : (0 : ℝ) < σ - 1 := by linarith
  have hNnn : (0 : ℝ) ≤ (N : ℝ) := Nat.cast_nonneg N
  have hcn : ‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ ≤ σ + |γ| := by
    have h := norm_center_le (σ := σ) (γ := γ)
    rwa [abs_of_pos (by linarith : (0 : ℝ) < σ)] at h
  -- inner-factor bound: Ms+1 ≤ M'
  have hpow : (2 : ℝ) ^ (-(σ - R)) ≤ 1 := by
    rw [Real.rpow_neg (by norm_num), inv_le_one₀ (by positivity)]
    exact Real.one_le_rpow (by norm_num) hσR.le
  have hNpow : (N : ℝ) * (2 : ℝ) ^ (-(σ - R)) ≤ (N : ℝ) := by
    calc (N : ℝ) * (2 : ℝ) ^ (-(σ - R)) ≤ (N : ℝ) * 1 :=
          mul_le_mul_of_nonneg_left hpow hNnn
      _ = (N : ℝ) := mul_one _
  have hterm2 : (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) * (N : ℝ) / (σ - R)
      ≤ (σ + |γ| + R) * (N : ℝ) / (σ - R) := by
    gcongr
  have hinner : 1 + ((N : ℝ) * (2 : ℝ) ^ (-(σ - R))
        + (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) * (N : ℝ) / (σ - R)) + 1
      ≤ 2 + (N : ℝ) + (σ + |γ| + R) * (N : ℝ) / (σ - R) := by linarith
  have hfac_pos : (0 : ℝ) < 3 + 1 / (σ - 1) := by positivity
  have hinner_pos : (0 : ℝ) < 1 + ((N : ℝ) * (2 : ℝ) ^ (-(σ - R))
        + (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) * (N : ℝ) / (σ - R)) + 1 := by
    have t1 : (0 : ℝ) ≤ (N : ℝ) * (2 : ℝ) ^ (-(σ - R)) :=
      mul_nonneg hNnn (Real.rpow_nonneg (by norm_num) _)
    have t2 : (0 : ℝ) ≤ (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) * (N : ℝ) / (σ - R) :=
      div_nonneg (mul_nonneg (add_nonneg (norm_nonneg _) hR.le) hNnn) hσR.le
    linarith
  have hprod := mul_le_mul_of_nonneg_right hinner hfac_pos.le
  have hlog := Real.log_le_log (mul_pos hinner_pos hfac_pos) hprod
  have hRhalf : (0 : ℝ) < R / 2 := by linarith
  calc 4 * Real.log ((1 + ((N : ℝ) * (2 : ℝ) ^ (-(σ - R))
            + (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) * (N : ℝ) / (σ - R)) + 1)
          * (3 + 1 / (σ - 1))) / (R / 2)
      ≤ 4 * Real.log ((2 + (N : ℝ) + (σ + |γ| + R) * (N : ℝ) / (σ - R)) * (3 + 1 / (σ - 1)))
          / (R / 2) := by gcongr
    _ ≤ B := hB

#print axioms hByscalar_of_explicit_LFunction
