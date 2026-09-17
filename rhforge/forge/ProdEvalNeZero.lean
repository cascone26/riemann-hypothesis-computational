import Mathlib.Algebra.Polynomial.Roots
import Mathlib.Tactic

/-!
# §1801 — The isolating product is nonzero off its root set

General-case converse, step 2.  `Q = ∏_{s∈S}(X − C s)` evaluated at `r₀ ∉ S` is `∏_{s∈S}(r₀ − s) ≠ 0`
(product of nonzeros in the domain ℂ).  This guarantees §1798's linear factor `L` can be chosen with
`L(r₀) = i / Q(r₀)`, so `P = Q·L` satisfies `P(r₀) = i` while vanishing on all of `S`.  EXIT0 / std-axioms.
-/

open Polynomial

/-- **Isolating product nonzero at `r₀ ∉ S`.** -/
theorem prod_eval_ne_zero (S : Multiset ℂ) (r₀ : ℂ) (h : r₀ ∉ S) :
    ((S.map (fun s => X - C s)).prod).eval r₀ ≠ 0 := by
  simp only [Polynomial.eval_multiset_prod, Multiset.map_map]
  apply Multiset.prod_ne_zero
  simp only [Multiset.mem_map, Function.comp_apply, Polynomial.eval_sub, Polynomial.eval_X,
    Polynomial.eval_C, not_exists, not_and]
  intro s hs hz
  apply h
  rw [sub_eq_zero.mp hz]; exact hs
