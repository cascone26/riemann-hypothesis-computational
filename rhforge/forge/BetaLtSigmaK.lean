import Mathlib.NumberTheory.LSeries.Nonvanishing

/-!
# §1455 — a ζ zero has `β < σK` (discharges §1454's `hβ`)

If `ζ(β+iγ) = 0` then `β < 1` (else `Re(β+iγ) = β ≥ 1` would force `ζ ≠ 0` by
`riemannZeta_ne_zero_of_one_le_re`), hence `β < 1 < 1 + 1/(2K')` for any `K' > 0`.  This discharges §1454's `hβ`
directly from the zero hypothesis.  EXIT0 / std-axioms.
-/

open Complex

theorem beta_lt_sigmaK {β γ K' : ℝ} (hK' : 0 < K')
    (hζ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0) : β < 1 + 1 / (2 * K') := by
  have hβ1 : β < 1 := by
    by_contra h
    push_neg at h
    have hre : (1 : ℝ) ≤ ((β : ℂ) + I * (γ : ℝ)).re := by
      simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
        Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
      exact h
    exact riemannZeta_ne_zero_of_one_le_re hre hζ
  have : (0 : ℝ) < 1 / (2 * K') := by positivity
  linarith

#print axioms beta_lt_sigmaK
