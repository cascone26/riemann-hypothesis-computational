import Mathlib.Analysis.SpecialFunctions.Complex.Analytic
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# §1477 — generic sphere factorization by division: `F = P·(F/P)` on a sphere where `P ≠ 0`

The generic (function-agnostic) component supplying the `hfac` hypothesis of the Borel–Carathéodory factor bound
§1471: wherever `P` is nonvanishing on `sphere c R`, taking `g := F/P` gives `F s = P s · g s` pointwise on the
sphere (`mul_div_cancel₀`).  Instantiated at `F = LFunction χ` and `P = ∏(s-ρ)^{m ρ}` (the half-disk zeros, §1476),
this is exactly the sphere-form factorization §1471 consumes — the analytic content then reduces to showing
`g = F/P` is holomorphic and nonvanishing on the inner disk (the generic §1383/§1346 bricks).  EXIT0 / std-axioms.
-/

open Complex Metric

theorem factor_div_on_sphere {F P : ℂ → ℂ} {c : ℂ} {R : ℝ}
    (hP : ∀ s ∈ sphere c R, P s ≠ 0) :
    ∀ s ∈ sphere c R, F s = P s * (F s / P s) := by
  intro s hs
  rw [mul_div_cancel₀ _ (hP s hs)]

#print axioms factor_div_on_sphere
