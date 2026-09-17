import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.Analysis.Complex.CauchyIntegral

open Complex DirichletCharacter

/-!
# §1523 — `LFunction χ` is `AnalyticOnNhd` on the open half-plane `{Re > 0}` (`χ ≠ 1`)

`LFunction χ` is entire (`differentiable_LFunction`), so it is `AnalyticOnNhd` on the open, convex (hence
preconnected) set `{s | 0 < s.re}` — the domain on which the analytic-continuation identity theorem (Track B)
identifies it with the Abel-summed twisted-series form.  Foundational input for the sub-`Re=1` growth extension.
EXIT0 / std-axioms.
-/

theorem LFunction_analyticOnNhd_re_pos {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1) :
    AnalyticOnNhd ℂ (LFunction χ) {s : ℂ | 0 < s.re} :=
  (differentiable_LFunction hχ).differentiableOn.analyticOnNhd
    (isOpen_lt continuous_const Complex.continuous_re)

theorem isPreconnected_re_pos : IsPreconnected {s : ℂ | 0 < s.re} :=
  (convex_halfSpace_re_gt 0).isPreconnected

#print axioms LFunction_analyticOnNhd_re_pos
