import Mathlib.NumberTheory.LSeries.DirichletContinuation

/-!
# §1472 — `logDeriv (LFunction χ) = logDeriv (LSeries ↗χ)` on `Re s > 1`

The DLVP positivity §1468 is stated with the naive `LSeries (↗χ)`, but the Borel–Carathéodory factor bound
§1471 works with the *entire* function `LFunction χ` (`differentiable_LFunction`, `χ ≠ 1`).  On the open
half-plane `Re s > 1` the two agree (`LFunction_eq_LSeries`), so — being equal on a neighbourhood of any such
`s` — their logarithmic derivatives agree there:

  `logDeriv (LFunction χ) s = logDeriv (LSeries (↗χ)) s`   for `1 < Re s`.

This is the transport lemma that lets the entire-function BC machinery consume the `LSeries`-side positivity.
EXIT0 / std-axioms.
-/

open Complex DirichletCharacter
open scoped LSeries.notation

theorem logDeriv_LFunction_eq_twist {N : ℕ} [NeZero N] (χ : DirichletCharacter ℂ N) {s : ℂ}
    (hs : 1 < s.re) :
    logDeriv (LFunction χ) s = logDeriv (LSeries (↗χ)) s := by
  have hopen : IsOpen {z : ℂ | 1 < z.re} := isOpen_lt continuous_const Complex.continuous_re
  have hmem : s ∈ {z : ℂ | 1 < z.re} := hs
  have heq : Set.EqOn (LFunction χ) (LSeries (↗χ)) {z : ℂ | 1 < z.re} := fun z hz =>
    LFunction_eq_LSeries χ hz
  have hnhds : LFunction χ =ᶠ[nhds s] LSeries (↗χ) := heq.eventuallyEq_of_mem (hopen.mem_nhds hmem)
  rw [logDeriv_apply, logDeriv_apply, hnhds.deriv_eq, hnhds.eq_of_nhds]

#print axioms logDeriv_LFunction_eq_twist
