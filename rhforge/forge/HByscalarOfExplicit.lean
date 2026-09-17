import InnerFactorLe

/-!
# §1452 — discharge §1449's `hByscalar` from an explicit `(σ,γ,R)` bound on `B`

Using §1451 (`Ms+1 ≤ M' := 3+1/(|γ|-R)+(σ+|γ|+R)/(σ-R)`) and `c.re = σ`, the giant BC log-argument
`(Ms+1)·(3+1/(c.re-1))` is bounded by `M'·(3+1/(σ-1))`; by log-monotonicity §1449's `hByscalar` holds whenever
`B ≥ 4·log(M'·(3+1/(σ-1)))/(R/2)` — a clean explicit inequality in `σ = σK` and `γ`.  EXIT0 / std-axioms.
-/

open Complex

theorem hByscalar_of_explicit {σ γ R B : ℝ} (hσ : 1 < σ) (hR : 0 < R) (hRγ : R < |γ|) (hRσ : R < σ)
    (hB : 4 * Real.log ((3 + 1 / (|γ| - R) + (σ + |γ| + R) / (σ - R)) * (3 + 1 / (σ - 1))) / (R / 2) ≤ B) :
    4 * Real.log ((2 + 1 / (‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ - R)
        + (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) / ((((σ : ℝ) : ℂ) + I * (γ : ℝ)).re - R) + 1)
        * (3 + 1 / ((((σ : ℝ) : ℂ) + I * (γ : ℝ)).re - 1))) / (R / 2) ≤ B := by
  have hcre : (((σ : ℝ) : ℂ) + I * (γ : ℝ)).re = σ := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  have hif := inner_factor_le hσ hR hRγ hRσ
  rw [hcre] at hif ⊢
  have hden2 : (0 : ℝ) < σ - R := by linarith
  have hden1 : (0 : ℝ) < |γ| - R := by linarith
  have hMs1_pos : 0 < 2 + 1 / (‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ - R)
      + (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) / (σ - R) + 1 := by
    have h1 : (0 : ℝ) ≤ 1 / (‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ - R) := by
      have : (0 : ℝ) < ‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ - R := by
        have := abs_le_norm_center_sub_one (σ := σ) (γ := γ); linarith
      positivity
    have h2 : (0 : ℝ) ≤ (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) / (σ - R) := by positivity
    linarith
  have hfac_pos : (0 : ℝ) < 3 + 1 / (σ - 1) := by
    have : (0 : ℝ) < σ - 1 := by linarith
    positivity
  have hprod := mul_le_mul_of_nonneg_right hif (le_of_lt hfac_pos)
  have hlog := Real.log_le_log (by positivity) hprod
  calc 4 * Real.log ((2 + 1 / (‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ - R)
          + (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) / (σ - R) + 1) * (3 + 1 / (σ - 1))) / (R / 2)
      ≤ 4 * Real.log ((3 + 1 / (|γ| - R) + (σ + |γ| + R) / (σ - R)) * (3 + 1 / (σ - 1))) / (R / 2) := by
        gcongr
    _ ≤ B := hB

#print axioms hByscalar_of_explicit
