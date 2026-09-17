import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.NumberTheory.LSeries.Nonvanishing

/-!
# §1499 — a zero of `LFunction χ` has `β ≤ 1`, hence `β < σK` (the region `hβ`)

L-analogue of ζ §1455.  `LFunction χ` (χ ≠ 1) is nonvanishing on `Re ≥ 1`
(`DirichletCharacter.LFunction_ne_zero_of_one_le_re`), so any zero `β+iγ` has `β ≤ 1`, hence
`β < σK = 1 + 1/(2K)` for `K > 0`.  This discharges the `hβ` hypothesis of the region connector §1498.
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem beta_le_one_of_LFunction_zero {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {β γ : ℝ} (hL : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0) : β ≤ 1 := by
  by_contra h
  push_neg at h
  have hre : (1 : ℝ) ≤ ((β : ℂ) + I * (γ : ℝ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
    linarith
  exact LFunction_ne_zero_of_one_le_re χ (Or.inl hχ) hre hL

theorem beta_lt_sigmaK_of_LFunction_zero {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {β γ A B C : ℝ} (hK : 0 < 3 * A + 4 * B + C)
    (hL : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0) :
    β < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
  have hβ1 : β ≤ 1 := beta_le_one_of_LFunction_zero hχ hL
  have hpos : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
  linarith

#print axioms beta_lt_sigmaK_of_LFunction_zero
