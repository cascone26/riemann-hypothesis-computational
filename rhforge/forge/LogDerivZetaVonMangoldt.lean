import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# §1356 — `logDeriv ζ s = - L Λ s` on `Re s > 1` (von Mangoldt Dirichlet series)

Foundation of the DLVP positivity (the `3-4-1` inequality).  Mathlib's
`LSeries_vonMangoldt_eq_deriv_riemannZeta_div` gives `L ↗Λ s = - deriv ζ s / ζ s`; since
`logDeriv ζ s = deriv ζ s / ζ s`, this reads `logDeriv ζ s = - L ↗Λ s`.  The Dirichlet series
`L ↗Λ s = ∑ Λ(n) n^{-s}` has nonnegative coefficients, so `Re (-logDeriv ζ)` at `σ, σ+iγ, σ+2iγ`
combines through `3 + 4cosθ + cos2θ ≥ 0` into the classical DLVP positivity feeding §1355's `hdlvp`.
EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries
open scoped LSeries.notation

theorem logDeriv_riemannZeta_eq_neg_vonMangoldt {s : ℂ} (hs : 1 < s.re) :
    logDeriv riemannZeta s = - LSeries (↗Λ) s := by
  have h := LSeries_vonMangoldt_eq_deriv_riemannZeta_div hs
  rw [logDeriv_apply, h]
  ring

#print axioms logDeriv_riemannZeta_eq_neg_vonMangoldt
