import Mathlib.Data.Complex.Basic
import Mathlib.Tactic

/-!
# §1798 — At a non-real node, a real linear form hits every complex value

Reusable core of the converse's interpolation step.  If `r₀` is non-real then `{1, r₀}` are ℝ-linearly
independent in ℂ, so the real-linear map `(a,b) ↦ a + b·r₀ : ℝ² → ℂ` is onto: for every target `w` there are
REAL `a, b` with `a + b r₀ = w`.  Explicitly `a = Re w − Im w·Re r₀/Im r₀`, `b = Im w/Im r₀`.  This is the
generalization of §1794 (the `w = i` case) and the engine that lets a real polynomial `P = Q·L` hit
`P(r₀) = i/Q(r₀)` while keeping real coefficients, completing the root-isolation for the Hermite converse.
EXIT0 / std-axioms.
-/

open Complex

/-- **Real linear surjectivity at a non-real node.**  For non-real `r₀` and any `w`, real `a,b` give
`a + b·r₀ = w`. -/
theorem real_linear_surj_at_nonreal (r₀ : ℂ) (hr : r₀.im ≠ 0) (w : ℂ) :
    ∃ a b : ℝ, (a : ℂ) + (b : ℂ) * r₀ = w := by
  refine ⟨w.re - w.im * r₀.re / r₀.im, w.im / r₀.im, ?_⟩
  apply Complex.ext <;>
    simp [Complex.add_re, Complex.add_im, Complex.mul_re, Complex.mul_im,
          Complex.ofReal_re, Complex.ofReal_im] <;>
    field_simp <;> ring
