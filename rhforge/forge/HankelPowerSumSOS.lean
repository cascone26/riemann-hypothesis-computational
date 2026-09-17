import Mathlib.Tactic

/-!
# §1790 — Hermite/Hankel power-sum form is a sum of squares (all-real ⟹ Hankel PSD)

Foundation brick of a real-rootedness module — the "Q has all real roots in [−2,2]" input that closes the
Carathéodory–Fejér / on-circle engine (§1781–1789).  Mathlib has NO Sturm / Bezoutian / Hankel-signature
theory (verified: the `Hermite` dir is Hermite *polynomials*, `Bezout` is the domain algebra), so this begins
that build from the cleanest atomic fact.

For a finite real family `r : Fin n → ℝ` (the roots) and a real coefficient vector `x : Fin (d+1) → ℝ`, the
quadratic form of the Newton power-sum Hankel matrix `H_{ij} = s_{i+j} = ∑_k r_k^{i+j}` equals a sum of squares:

    ∑_{i,j} (∑_k r_k^{i+j}) x_i x_j  =  ∑_k (∑_i x_i r_k^i)²   ( ≥ 0 ).

This is exactly the easy (necessary) direction of Hermite's real-root-counting criterion: when a polynomial
splits over ℝ, its power-sum Hankel matrix is positive semidefinite.  EXIT0 / std-axioms.
-/

open Finset

variable {n d : ℕ}

/-- **Hankel power-sum form = sum of squares.**  For a real root family `r` and coefficient vector `x`, the
quadratic form of the Newton power-sum Hankel matrix is `∑_k (∑_i x_i r_k^i)²`. -/
theorem hankel_powersum_eq_sum_sq (r : Fin n → ℝ) (x : Fin (d + 1) → ℝ) :
    ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
        (∑ k : Fin n, r k ^ ((i : ℕ) + (j : ℕ))) * x i * x j
      = ∑ k : Fin n, (∑ i : Fin (d + 1), x i * r k ^ (i : ℕ)) ^ 2 := by
  have key : ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
        (∑ k : Fin n, r k ^ ((i : ℕ) + (j : ℕ))) * x i * x j
      = ∑ i : Fin (d + 1), ∑ j : Fin (d + 1), ∑ k : Fin n,
          (x i * r k ^ (i : ℕ)) * (x j * r k ^ (j : ℕ)) := by
    refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
    rw [Finset.sum_mul, Finset.sum_mul]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [pow_add]; ring
  have rhs : ∑ k : Fin n, (∑ i : Fin (d + 1), x i * r k ^ (i : ℕ)) ^ 2
      = ∑ k : Fin n, ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
          (x i * r k ^ (i : ℕ)) * (x j * r k ^ (j : ℕ)) := by
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [sq, Finset.sum_mul_sum]
  rw [key, rhs]
  conv_lhs => enter [2, i]; rw [Finset.sum_comm]
  rw [Finset.sum_comm]

/-- **Hankel power-sum form is PSD.**  Immediate: the form is a sum of squares, hence nonnegative. -/
theorem hankel_powersum_form_nonneg (r : Fin n → ℝ) (x : Fin (d + 1) → ℝ) :
    0 ≤ ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
          (∑ k : Fin n, r k ^ ((i : ℕ) + (j : ℕ))) * x i * x j := by
  rw [hankel_powersum_eq_sum_sq]
  exact Finset.sum_nonneg (fun k _ => sq_nonneg _)
