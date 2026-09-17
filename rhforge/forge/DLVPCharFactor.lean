import Mathlib.Analysis.SpecialFunctions.Complex.Circle

/-!
# §1465 — the character 3-4-1 positivity core `0 ≤ 3 + 4·Re w + Re(w²)` for `‖w‖ ≤ 1`

The abstract heart of the de la Vallée Poussin `3-4-1` inequality FOR DIRICHLET L-FUNCTIONS.  In the L-series
combination `ζ³·L(·,χ)⁴·L(·,χ²)` the per-term factor is `3 + 4Re(w) + Re(w²)` with `w = χ(n)·n^{-it}`, `‖w‖ ≤ 1`
(`= 1` for `n` coprime to the modulus, `= 0` otherwise).  Since `Re(w²) = Re(w)² - Im(w)²`,

  `3 + 4Re w + Re(w²) = 2(1+Re w)² + (1 - ‖w‖²) ≥ 0`.

This is the L-function generalization of the ζ per-term positivity (§1357) and the first brick of the
Dirichlet-L zero-free region.  EXIT0 / std-axioms.
-/

theorem dlvp_char_factor_nonneg {w : ℂ} (hw : ‖w‖ ≤ 1) : 0 ≤ 3 + 4 * w.re + (w ^ 2).re := by
  have hsq : (w ^ 2).re = w.re ^ 2 - w.im ^ 2 := by rw [sq, Complex.mul_re]; ring
  have hnorm : w.re ^ 2 + w.im ^ 2 ≤ 1 := by
    have h1 : ‖w‖ ^ 2 ≤ 1 := by nlinarith [norm_nonneg w, hw]
    have h2 : ‖w‖ ^ 2 = w.re ^ 2 + w.im ^ 2 := by
      have := Complex.sq_norm_sub_sq_im w; linarith
    linarith
  rw [hsq]
  nlinarith [sq_nonneg (1 + w.re), hnorm]

#print axioms dlvp_char_factor_nonneg
