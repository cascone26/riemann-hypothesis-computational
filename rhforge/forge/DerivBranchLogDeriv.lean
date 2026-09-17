import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Analysis.SpecialFunctions.Complex.Analytic
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.Calculus.LogDeriv

/-!
# §1321 — A log-branch's derivative is the logarithmic derivative: `deriv f = logDeriv h`

If `f` is a holomorphic branch of `log h` (`exp∘f = h`) with `HasDerivAt f a z` and `HasDerivAt h b z`, then
`a = logDeriv h z`.  Proof: differentiate `exp∘f = h` (chain rule `HasDerivAt.cexp`): `exp(f z)·a = b`, and
`h z = exp(f z) ≠ 0`, so `a = b/h z = logDeriv h z`.  Together with §1320 (`Re f = log‖h‖`) this shows a branch
supplies BOTH inputs Borel–Carathéodory §1318 needs — reducing the whole ζ'/ζ estimate to branch EXISTENCE alone.
EXIT0 / std-axioms.
-/

open Complex

theorem deriv_branch_eq_logDeriv {f h : ℂ → ℂ} {z a b : ℂ}
    (hf : HasDerivAt f a z) (hh : HasDerivAt h b z) (hfh : ∀ w, Complex.exp (f w) = h w) :
    a = logDeriv h z := by
  have hexp : HasDerivAt (fun w => Complex.exp (f w)) (Complex.exp (f z) * a) z :=
    HasDerivAt.cexp hf
  rw [funext hfh] at hexp
  have hab : Complex.exp (f z) * a = b := hexp.unique hh
  have hez : Complex.exp (f z) ≠ 0 := Complex.exp_ne_zero _
  rw [logDeriv_apply, hh.deriv, ← hfh z, eq_div_iff hez]
  linear_combination hab

#print axioms deriv_branch_eq_logDeriv
