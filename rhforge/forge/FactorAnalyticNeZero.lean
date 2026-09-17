import Mathlib.Analysis.Meromorphic.FactorizedRational
import Mathlib.Analysis.Analytic.Constructions
import Mathlib.Analysis.Complex.Basic

/-!
# §1346 — The factorization factor `(∏ᶠ(·-u)^{d u})•g` is analytic and nonzero off the divisor support

At a point `z` with `d z = 0` (outside the divisor support — e.g. the DLVP test point `σ+iγ`, `Re > 1`, no zeros/poles),
the factorized rational `∏ᶠ u,(·-u)^{d u}` is analytic (`FactorizedRational.analyticAt`, `0 ≤ d z`) and nonzero
(`FactorizedRational.ne_zero`).  If `g` is analytic and nonzero at `z`, the product `(∏ᶠ...)•g` is analytic and
nonzero there.  This supplies the `AnalyticAt` (for §1344 ⟹ `ContinuousAt (logDeriv F)`) and nonvanishing facts for
the extract_zeros_poles factor `F` at the test point.  EXIT0 / std-axioms.
-/

open Complex

theorem analyticAt_and_ne_zero_factor {d : ℂ → ℤ} {g : ℂ → ℂ} {z : ℂ}
    (hdz : d z = 0) (hg : AnalyticAt ℂ g z) (hgz : g z ≠ 0) :
    AnalyticAt ℂ ((∏ᶠ u, (· - u) ^ (d u)) • g) z ∧ ((∏ᶠ u, (· - u) ^ (d u)) • g) z ≠ 0 := by
  have hFRa : AnalyticAt ℂ (∏ᶠ u, (· - u) ^ (d u)) z :=
    Function.FactorizedRational.analyticAt (le_of_eq hdz.symm)
  have hFRne : (∏ᶠ u, (· - u) ^ (d u)) z ≠ 0 := Function.FactorizedRational.ne_zero hdz
  have hmul : ((∏ᶠ u, (· - u) ^ (d u)) • g) = fun s => (∏ᶠ u, (· - u) ^ (d u)) s * g s := by
    funext s; rw [Pi.smul_apply', smul_eq_mul]
  refine ⟨?_, ?_⟩
  · rw [hmul]; exact hFRa.mul hg
  · rw [Pi.smul_apply', smul_eq_mul]; exact mul_ne_zero hFRne hgz

#print axioms analyticAt_and_ne_zero_factor
