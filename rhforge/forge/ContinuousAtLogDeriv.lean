import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.Analytic.Within
import Mathlib.Analysis.Complex.Basic

/-!
# §1344 — `logDeriv f` is continuous where `f` is analytic and nonzero

`logDeriv f = deriv f / f`.  Where `f` is analytic at `z`, `deriv f` is analytic (`AnalyticAt.deriv`) hence
continuous, and `f` is continuous; with `f z ≠ 0` the quotient is continuous at `z`.  This supplies the two
`ContinuousAt (logDeriv ·) z` hypotheses of §1342 at the DLVP test point `σ+iγ` (`Re > 1`), where both `ζ` and the
factorization factor `(∏ᶠ...)•g` are analytic and nonvanishing.  EXIT0 / std-axioms.
-/

open Complex

theorem continuousAt_logDeriv {f : ℂ → ℂ} {z : ℂ} (hf : AnalyticAt ℂ f z) (hfz : f z ≠ 0) :
    ContinuousAt (logDeriv f) z := by
  have heq : logDeriv f = fun x => deriv f x / f x := by funext x; rw [logDeriv_apply]
  rw [heq]
  exact (hf.deriv.continuousAt).div hf.continuousAt hfz

#print axioms continuousAt_logDeriv
