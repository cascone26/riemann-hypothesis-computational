import DicksonSum

/-!
# §1813 — The Chebyshev-image polynomial `Qpoly` as a genuine `Polynomial`

The on-circle criterion (§1789/§1812) uses `Q(s) = v m + ∑_{k<m} v(m+1+k)·D_{k+1}(s)` as an evaluated
expression.  Here we package it as an actual `Polynomial ℂ` `Qpoly` (a Dickson combination) and prove its
evaluation matches, so its ROOT MULTISET `(Qpoly v m).roots` is available to feed the real-rootedness iffs
(§1806/§1811).  EXIT0 / std-axioms.
-/

open Polynomial Finset

/-- The half-degree Chebyshev-image polynomial as a Dickson combination. -/
noncomputable def Qpoly (v : ℕ → ℂ) (m : ℕ) : ℂ[X] :=
  C (v m) + ∑ k ∈ Finset.range m, C (v (m + 1 + k)) * dickson (k + 1)

/-- **`Qpoly` evaluates to the on-circle criterion's `Q`.** -/
theorem Qpoly_eval (v : ℕ → ℂ) (m : ℕ) (s : ℂ) :
    (Qpoly v m).eval s = v m + ∑ k ∈ Finset.range m, v (m + 1 + k) * (dickson (k + 1)).eval s := by
  simp only [Qpoly, Polynomial.eval_add, Polynomial.eval_C, Polynomial.eval_finset_sum,
    Polynomial.eval_mul]
