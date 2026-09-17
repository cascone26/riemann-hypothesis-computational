import Mathlib

namespace CMKeyLemma

/-- **Geronimus simplicity core: a tridiagonal eigenvector is determined by its first component.**

For a real symmetric **tridiagonal (Jacobi)** matrix with diagonal `a` and off-diagonals `b` (all
`b i ≠ 0`), an eigenvector `v` for eigenvalue `μ` obeys the 3-term recurrence
`a i · v i + b i · v (i+1) + b (i-1) · v (i-1) = μ · v i` (with the `i=0` boundary having no `v (-1)`
term).  If the first component `v 0 = 0`, then `v = 0` identically.

Consequently the map `v ↦ v 0` is injective on each eigenspace, so **every eigenspace is at most
one-dimensional ⟹ every eigenvalue is simple** (the classical Geronimus/Jacobi result).  This reduces CM's
open "simple" (of the Weil form `QW_λ`) to the *structural* hypothesis "the Jacobi off-diagonals are
nonzero" — a cleaner named sufficient condition, extending the even-simple decomposition
(`ParityDefinite`, `EvenOfRealFourier`).  Machine-verified, axiom-clean. -/
theorem tridiag_eigvec_determined_by_first (a b v : ℕ → ℝ) (μ : ℝ)
    (hb : ∀ i, b i ≠ 0)
    (h0 : a 0 * v 0 + b 0 * v 1 = μ * v 0)
    (hrec : ∀ i, a (i + 1) * v (i + 1) + b (i + 1) * v (i + 2) + b i * v i = μ * v (i + 1))
    (hv0 : v 0 = 0) :
    ∀ i, v i = 0 := by
  -- prove the two-step invariant `v i = 0 ∧ v (i+1) = 0`
  have key : ∀ i, v i = 0 ∧ v (i + 1) = 0 := by
    intro i
    induction i with
    | zero =>
      refine ⟨hv0, ?_⟩
      have hb0 : b 0 * v 1 = 0 := by
        have := h0; rw [hv0, mul_zero, mul_zero, zero_add] at this; exact this
      exact (mul_eq_zero.mp hb0).resolve_left (hb 0)
    | succ k ih =>
      obtain ⟨hk, hk1⟩ := ih
      refine ⟨hk1, ?_⟩
      have hbk : b (k + 1) * v (k + 2) = 0 := by
        have := hrec k
        rw [hk, hk1, mul_zero, mul_zero, mul_zero, add_zero, zero_add] at this
        exact this
      exact (mul_eq_zero.mp hbk).resolve_left (hb (k + 1))
  intro i
  exact (key i).1

#print axioms tridiag_eigvec_determined_by_first

end CMKeyLemma
