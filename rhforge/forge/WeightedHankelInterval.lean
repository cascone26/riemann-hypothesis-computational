import Mathlib.Tactic

/-!
# §1791 — Weighted Hankel power-sum form + [−2,2] interval localization

Extends §1790 (HankelPowerSumSOS) with an arbitrary real weight `w : Fin n → ℝ`, then specializes to the
weight `w_k = 4 − r_k²` which is `≥ 0` exactly when `|r_k| ≤ 2`.  This is the INTERVAL half of the Hermite
real-rootedness criterion needed to close the on-circle engine (§1781–1789): "all roots real AND in [−2,2]".

    ∑_{i,j} (∑_k w_k r_k^{i+j}) x_i x_j  =  ∑_k w_k (∑_i x_i r_k^i)²  .

With `w_k ≥ 0` the RHS is a nonneg-weighted sum of squares, hence `≥ 0`; and `w_k = 4 − r_k² ≥ 0` under
`r_k² ≤ 4` gives the secondary Hankel form whose positivity localizes the roots to `[−2,2]`.
EXIT0 / std-axioms.
-/

open Finset

variable {n d : ℕ}

/-- **Weighted Hankel power-sum form = weighted sum of squares.** -/
theorem weighted_hankel_powersum_eq_sum_sq (w r : Fin n → ℝ) (x : Fin (d + 1) → ℝ) :
    ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
        (∑ k : Fin n, w k * r k ^ ((i : ℕ) + (j : ℕ))) * x i * x j
      = ∑ k : Fin n, w k * (∑ i : Fin (d + 1), x i * r k ^ (i : ℕ)) ^ 2 := by
  have key : ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
        (∑ k : Fin n, w k * r k ^ ((i : ℕ) + (j : ℕ))) * x i * x j
      = ∑ i : Fin (d + 1), ∑ j : Fin (d + 1), ∑ k : Fin n,
          w k * ((x i * r k ^ (i : ℕ)) * (x j * r k ^ (j : ℕ))) := by
    refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
    rw [Finset.sum_mul, Finset.sum_mul]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [pow_add]; ring
  have rhs : ∑ k : Fin n, w k * (∑ i : Fin (d + 1), x i * r k ^ (i : ℕ)) ^ 2
      = ∑ k : Fin n, ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
          w k * ((x i * r k ^ (i : ℕ)) * (x j * r k ^ (j : ℕ))) := by
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [sq, Finset.sum_mul_sum]
    simp_rw [Finset.mul_sum]
  rw [key, rhs]
  conv_lhs => enter [2, i]; rw [Finset.sum_comm]
  rw [Finset.sum_comm]

/-- **Nonneg-weight Hankel form is PSD.** -/
theorem weighted_hankel_form_nonneg (w r : Fin n → ℝ) (x : Fin (d + 1) → ℝ)
    (hw : ∀ k, 0 ≤ w k) :
    0 ≤ ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
          (∑ k : Fin n, w k * r k ^ ((i : ℕ) + (j : ℕ))) * x i * x j := by
  rw [weighted_hankel_powersum_eq_sum_sq]
  exact Finset.sum_nonneg (fun k _ => mul_nonneg (hw k) (sq_nonneg _))

/-- **[−2,2] interval localization.**  The secondary Hankel form with weight `4 − r_k²` is PSD whenever every
`|r_k| ≤ 2` (i.e. `r_k² ≤ 4`) — the Hermite-criterion witness that the roots lie in `[−2,2]`. -/
theorem hankel_interval_form_nonneg (r : Fin n → ℝ) (x : Fin (d + 1) → ℝ)
    (hr : ∀ k, r k ^ 2 ≤ 4) :
    0 ≤ ∑ i : Fin (d + 1), ∑ j : Fin (d + 1),
          (∑ k : Fin n, (4 - r k ^ 2) * r k ^ ((i : ℕ) + (j : ℕ))) * x i * x j :=
  weighted_hankel_form_nonneg (fun k => 4 - r k ^ 2) r x (fun k => by linarith [hr k])
