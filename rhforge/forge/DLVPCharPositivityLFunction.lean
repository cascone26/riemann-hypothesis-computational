import DLVPCharLogDeriv
import LogDerivLFunctionEqTwist

/-!
# §1497 — DLVP `3-4-1` positivity in `LFunction` form (the region `hdlvp`)

Transporting §1468's positivity (stated with the naive `LSeries ↗χ`, `LSeries ↗χ²`) onto the entire `LFunction χ`,
`LFunction χ²` via §1472 (`logDeriv (LFunction χ) = logDeriv (LSeries ↗χ)` on `Re > 1`).  For `σ > 1`:

  `0 ≤ 3·Re(-logDeriv ζ σ) + 4·Re(-logDeriv (LFunction χ)(σ+it)) + Re(-logDeriv (LFunction χ²)(σ+2it))`.

This is exactly the `hdlvp` input the generic region composition §1336 consumes, with the three terms now in the
`LFunction` form matching the decomposition §1493 (`F`), the ζ pole (`T₀`), and the growth §1491 (`T₂`).  Needs only
`[NeZero N]` (the transport §1472 holds for every character).  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries DirichletCharacter
open scoped LSeries.notation

theorem dlvp_char_positivity_LFunction {N : ℕ} [NeZero N] (χ : DirichletCharacter ℂ N) {σ t : ℝ}
    (hσ : 1 < σ) :
    0 ≤ 3 * (-logDeriv riemannZeta ((σ : ℂ))).re
      + 4 * (-logDeriv (LFunction χ) ((σ : ℂ) + I * (t : ℝ))).re
      + (-logDeriv (LFunction (χ ^ 2)) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ))).re := by
  have hs1 : 1 < ((σ : ℂ) + I * (t : ℝ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have hs2 : 1 < ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ)).re := by
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]; simpa using hσ
  have h := dlvp_char_positivity_logDeriv (t := t) χ hσ
  rw [show logDeriv (LSeries (↗χ)) ((σ : ℂ) + I * (t : ℝ))
        = logDeriv (LFunction χ) ((σ : ℂ) + I * (t : ℝ)) from (logDeriv_LFunction_eq_twist χ hs1).symm,
    show logDeriv (LSeries (↗(χ ^ 2))) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ))
        = logDeriv (LFunction (χ ^ 2)) ((σ : ℂ) + I * ((2 * t : ℝ) : ℂ)) from
      (logDeriv_LFunction_eq_twist (χ ^ 2) hs2).symm] at h
  exact h

#print axioms dlvp_char_positivity_LFunction
