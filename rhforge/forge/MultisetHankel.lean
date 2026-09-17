import Mathlib.Tactic

/-!
# §1792 — Multiset-level Hankel power-sum SOS (bridge to `Polynomial.roots`)

Lifts §1790's Hankel sum-of-squares identity from a `Fin n → ℝ` root family to a `Multiset ℝ` — which is
exactly what `Polynomial.roots` returns, so this is the connective tissue letting a concrete real polynomial's
roots feed the Hermite PSD criterion (§1790/§1791).

For `R : Multiset ℝ` (the roots) and `x : Fin (d+1) → ℝ`:

    ∑_{i,j} (∑_{r∈R} r^{i+j}) x_i x_j  =  ∑_{r∈R} (∑_i x_i r^i)²   ( ≥ 0 ).

EXIT0 / std-axioms.
-/

open Finset

variable {d : ℕ}

/-- **Single-root Hankel identity.**  For one real `a`, the rank-1 power-sum form is a perfect square. -/
theorem single_root_hankel (a : ℝ) (x : Fin (d + 1) → ℝ) :
    ∑ i : Fin (d + 1), ∑ j : Fin (d + 1), a ^ ((i : ℕ) + (j : ℕ)) * x i * x j
      = (∑ i : Fin (d + 1), x i * a ^ (i : ℕ)) ^ 2 := by
  rw [sq, Finset.sum_mul_sum]
  refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
  rw [pow_add]; ring

/-- **Multiset Hankel form = sum of squares over roots.** -/
theorem multiset_hankel_eq_sum_sq (R : Multiset ℝ) (x : Fin (d + 1) → ℝ) :
    ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
        (R.map (fun r => r ^ ((i : ℕ) + (j : ℕ)))).sum * x i * x j
      = (R.map (fun r => (∑ i : Fin (d + 1), x i * r ^ (i : ℕ)) ^ 2)).sum := by
  induction R using Multiset.induction with
  | empty => simp
  | cons a R ih =>
    simp only [Multiset.map_cons, Multiset.sum_cons, add_mul, Finset.sum_add_distrib]
    rw [single_root_hankel, ih]

/-- **Multiset Hankel form is PSD** (nonneg), directly usable on `Polynomial.roots`. -/
theorem multiset_hankel_nonneg (R : Multiset ℝ) (x : Fin (d + 1) → ℝ) :
    0 ≤ ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
          (R.map (fun r => r ^ ((i : ℕ) + (j : ℕ)))).sum * x i * x j := by
  rw [multiset_hankel_eq_sum_sq]
  refine Multiset.sum_nonneg ?_
  intro y hy
  simp only [Multiset.mem_map] at hy
  obtain ⟨r, _, rfl⟩ := hy
  exact sq_nonneg _
