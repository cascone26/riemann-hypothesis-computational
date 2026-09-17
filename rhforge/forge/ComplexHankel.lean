import Mathlib.Tactic
import Mathlib.Data.Complex.Basic

/-!
# §1795 — Complex-root Hankel form = sum of squared evaluations (over ℂ)

The object the Hermite converse acts on.  For a multiset of COMPLEX roots `R` (the full root set of a real
polynomial, closed under conjugation) and real coefficients `c : Fin (d+1) → ℝ`:

    ∑_{i,j} (∑_{r∈R} r^{i+j}) c_i c_j  =  ∑_{r∈R} (∑_i c_i r^i)²   (in ℂ).

The power sums `∑_{r∈R} r^{i+j}` are real (symmetric functions of a conj-closed set), so this is a real
quadratic form in `c`; the RHS exposes it as `∑ over roots of P_c(r)²`.  Paired with §1794 (a non-real root's
conjugate pair contributes `-2`), positivity of this form forbids non-real roots — the Hermite hard direction.
EXIT0 / std-axioms.
-/

open Finset

variable {d : ℕ}

/-- single-root complex identity -/
theorem single_root_hankel_C (a : ℂ) (c : Fin (d + 1) → ℝ) :
    ∑ i : Fin (d + 1), ∑ j : Fin (d + 1), a ^ ((i : ℕ) + (j : ℕ)) * (c i : ℂ) * (c j : ℂ)
      = (∑ i : Fin (d + 1), (c i : ℂ) * a ^ (i : ℕ)) ^ 2 := by
  rw [sq, Finset.sum_mul_sum]
  refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
  rw [pow_add]; ring

/-- **Complex-root Hankel form = sum of squared evaluations.** -/
theorem complex_hankel_eq_sum_sq (R : Multiset ℂ) (c : Fin (d + 1) → ℝ) :
    ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
        (R.map (fun r => r ^ ((i : ℕ) + (j : ℕ)))).sum * (c i : ℂ) * (c j : ℂ)
      = (R.map (fun r => (∑ i : Fin (d + 1), (c i : ℂ) * r ^ (i : ℕ)) ^ 2)).sum := by
  induction R using Multiset.induction with
  | empty => simp
  | cons a R ih =>
    simp only [Multiset.map_cons, Multiset.sum_cons, add_mul, Finset.sum_add_distrib]
    rw [single_root_hankel_C, ih]
