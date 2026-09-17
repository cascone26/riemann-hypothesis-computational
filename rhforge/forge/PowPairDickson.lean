import DicksonSum

/-!
# §1786 — Pair-to-Dickson: z^{m+k} + z^{m-k} = zᵐ·D_k(z+z⁻¹)

The atomic pairing step for the palindromic → `zᵐ·Q(z+z⁻¹)` factorization.  For `z ≠ 0` and `k ≤ m`, the
symmetric pair of powers about `zᵐ` collapses to `zᵐ` times a Dickson value:

    z^{m+k} + z^{m-k}  =  zᵐ · D_k(z+z⁻¹).

Summed over the palindromic pairs of a degree-`2m` polynomial, this yields `P(z) = zᵐ·Q(z+z⁻¹)` with
`Q = Σ cₖ D_k`.  Direct consequence of §1783 (`dickson_eval_sum`) plus `z^m·(z^k)⁻¹ = z^{m-k}` (needs `k ≤ m`).
EXIT0 / std-axioms.
-/

open Polynomial

variable {K : Type*} [Field K]

/-- **Pair-to-Dickson.**  For `z ≠ 0` and `k ≤ m`, `z^{m+k} + z^{m-k} = zᵐ · D_k(z+z⁻¹)`. -/
theorem pow_pair_eq_dickson (z : K) (hz : z ≠ 0) (m k : ℕ) (hk : k ≤ m) :
    z ^ (m + k) + z ^ (m - k) = z ^ m * (dickson k).eval (z + z⁻¹) := by
  have hzk : z ^ k ≠ 0 := pow_ne_zero k hz
  have hmk : z ^ (m - k) = z ^ m * (z ^ k)⁻¹ := by
    rw [eq_comm, mul_inv_eq_iff_eq_mul₀ hzk, ← pow_add, Nat.sub_add_cancel hk]
  rw [dickson_eval_sum z hz k, inv_pow, mul_add, ← pow_add, hmk]
