import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Tactic

/-!
# §1794 — Hermite converse seed: a non-real root's conjugate pair can be made negative

The load-bearing crux of the converse (PSD Hankel ⟹ all roots real).  Over ℂ the Hankel form equals
`∑_{r∈roots} (P_x(r))²` with `P_x(t) = ∑ c_i t^i` (real coefficients `c_i`).  Real roots contribute
squares `≥ 0`, but a NON-REAL root `r` and its conjugate `r̄` contribute `(P_x(r))² + (P_x(r̄))²`, which
real coefficients can drive strictly negative.  Concretely, choosing `c₀ = -Re r / Im r`, `c₁ = 1/Im r`
makes `P_x(r) = i` and `P_x(r̄) = -i`, so the pair contributes `i² + (-i)² = -2 < 0`.

This is exactly why positivity of the (all-roots) Hankel form forbids non-real roots — the seed that, assembled
over the root multiset, yields the Hermite criterion's hard direction.  EXIT0 / std-axioms.
-/

open Complex

/-- **Conjugate-pair negativity.**  For a non-real `r`, real linear coefficients make the conjugate-pair
quadratic contribution equal `-2` (hence `< 0`). -/
theorem conj_pair_eval_neg (r : ℂ) (hr : r.im ≠ 0) :
    ∃ c₀ c₁ : ℝ,
      ((c₀ : ℂ) + (c₁ : ℂ) * r) ^ 2 + ((c₀ : ℂ) + (c₁ : ℂ) * (starRingEnd ℂ r)) ^ 2 = -2 := by
  refine ⟨-r.re / r.im, 1 / r.im, ?_⟩
  have hI : ((-r.re / r.im : ℝ) : ℂ) + ((1 / r.im : ℝ) : ℂ) * r = Complex.I := by
    apply Complex.ext <;>
      simp [Complex.add_re, Complex.add_im, Complex.mul_re, Complex.mul_im,
            Complex.ofReal_re, Complex.ofReal_im] <;>
      field_simp <;> ring
  have hIbar : ((-r.re / r.im : ℝ) : ℂ) + ((1 / r.im : ℝ) : ℂ) * (starRingEnd ℂ r) = -Complex.I := by
    apply Complex.ext <;>
      simp [Complex.add_re, Complex.add_im, Complex.mul_re, Complex.mul_im,
            Complex.ofReal_re, Complex.ofReal_im, Complex.conj_re, Complex.conj_im] <;>
      field_simp <;> ring
  rw [hI, hIbar, neg_sq, Complex.I_sq]
  norm_num
