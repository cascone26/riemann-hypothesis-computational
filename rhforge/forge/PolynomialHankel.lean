import MultisetHankel
import Mathlib.Algebra.Polynomial.Roots

/-!
# §1793 — Polynomial-level Hankel PSD (the usable Hermite necessary condition)

Applies §1792 (`MultisetHankel`) to `p.roots`, the multiset of real roots of a real polynomial.  Since
`Polynomial.roots` collects only roots in the base field `ℝ`, the power-sum Hankel form built from `p.roots`
is *always* a sum of squares over those real roots, hence PSD.  This is the Hermite real-rootedness criterion's
necessary direction stated at the level of a concrete polynomial `p : ℝ[X]` — the form you actually evaluate.

    ∑_{i,j} (∑_{r∈p.roots} r^{i+j}) x_i x_j  =  ∑_{r∈p.roots} (∑_i x_i r^i)²  ≥ 0 .

EXIT0 / std-axioms.
-/

open Polynomial Finset

variable {d : ℕ}

/-- **Root-power-sum Hankel form of a real polynomial = sum of squares over its real roots.** -/
theorem polynomial_roots_hankel_eq_sum_sq (p : ℝ[X]) (x : Fin (d + 1) → ℝ) :
    ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
        (p.roots.map (fun r => r ^ ((i : ℕ) + (j : ℕ)))).sum * x i * x j
      = (p.roots.map (fun r => (∑ i : Fin (d + 1), x i * r ^ (i : ℕ)) ^ 2)).sum :=
  multiset_hankel_eq_sum_sq p.roots x

/-- **Hermite necessary condition (Polynomial level): the root-power-sum Hankel form is PSD.** -/
theorem polynomial_roots_hankel_nonneg (p : ℝ[X]) (x : Fin (d + 1) → ℝ) :
    0 ≤ ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
          (p.roots.map (fun r => r ^ ((i : ℕ) + (j : ℕ)))).sum * x i * x j :=
  multiset_hankel_nonneg p.roots x

/-- **Zeroth power sum = number of real roots (with multiplicity).**  The `(0,0)` Hankel datum `s_0` counts
the real roots — the top-left entry that, in the full Hermite criterion, equals the rank ceiling. -/
theorem polynomial_powersum_zero (p : ℝ[X]) :
    (p.roots.map (fun r => r ^ (0 : ℕ))).sum = Multiset.card p.roots := by
  simp

