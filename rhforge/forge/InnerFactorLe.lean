import NormCenterLe
import AbsLeNormCenterSubOne

/-!
# §1451 — bound §1449's BC inner factor by an explicit `(σ,γ)` expression

`Ms+1 = 2 + 1/(‖c-1‖-R) + (‖c‖+R)/(c.re-R) + 1` with `c = σ+iγ` is bounded above using `‖c-1‖ ≥ |γ|` (§1446),
`‖c‖ ≤ σ+|γ|` (§1450, `σ>0`), and `c.re = σ`:

  `Ms+1 ≤ 3 + 1/(|γ|-R) + (σ+|γ|+R)/(σ-R)`.

Fully explicit in `σ,γ,R` — the giant norm terms are gone.  Feeds the log step of the constant choice.
EXIT0 / std-axioms.
-/

open Complex

theorem inner_factor_le {σ γ R : ℝ} (hσ : 1 < σ) (hR : 0 < R) (hRγ : R < |γ|) (hRσ : R < σ) :
    2 + 1 / (‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ - R)
        + (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) / ((((σ : ℝ) : ℂ) + I * (γ : ℝ)).re - R) + 1
      ≤ 3 + 1 / (|γ| - R) + (σ + |γ| + R) / (σ - R) := by
  have hcre : (((σ : ℝ) : ℂ) + I * (γ : ℝ)).re = σ := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  have hc1 : |γ| ≤ ‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ := abs_le_norm_center_sub_one
  have hcn : ‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ ≤ σ + |γ| := by
    have h := norm_center_le (σ := σ) (γ := γ)
    rwa [abs_of_pos (by linarith : (0 : ℝ) < σ)] at h
  have hden1 : (0 : ℝ) < |γ| - R := by linarith
  have hden2 : (0 : ℝ) < σ - R := by linarith
  rw [hcre]
  have hb1 : 1 / (‖((σ : ℝ) : ℂ) + I * (γ : ℝ) - 1‖ - R) ≤ 1 / (|γ| - R) :=
    one_div_le_one_div_of_le hden1 (by linarith)
  have hb2 : (‖((σ : ℝ) : ℂ) + I * (γ : ℝ)‖ + R) / (σ - R) ≤ (σ + |γ| + R) / (σ - R) := by
    gcongr
  linarith [hb1, hb2]

#print axioms inner_factor_le
