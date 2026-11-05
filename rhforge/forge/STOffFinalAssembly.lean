import HLGCRowSumVanishes
import STOffRealSymm
import STOffRealValued

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: THE MATRIX NULL-VECTOR CLAIM CLOSES — `S_T_off_real *ᵥ e_off = 0`

Step 4-g of §2391's 4-step plan — the literal closing theorem of the entire multi-cycle
matrix-bridge arc that opened at §2380. Assembles §2411-2412's `H_L_G_C_row_sum_vanishes`/`_conj`
against `S_T_off_sum`'s own definition (`∑_k(H_L_G_C(z_k)+H_L_G_C(z̄_k))`, `STOffPairedSum.lean`):
unfolds the matrix sum, swaps the `k`/`n` finite sums (`Finset.sum_comm`), distributes, and applies
the two row-vanishing facts termwise (giving `0+0` for every `k`). The resulting complex-valued
statement is then cast down to the REAL matrix `S_T_off_real` (§2355) via `S_T_off_sum_eq_ofReal_re`
(§2356) and injectivity of the ℝ→ℂ cast, landing the literal `Matrix.mulVec` null-vector statement:
`S_T_off_real *ᵥ e_off = 0`, i.e. Theorem 7.2's Eq (7.16)(ii) claim. Search-before-build: grepped
`S_T_off_sum_row_sum_vanishes`/`S_T_off_real_mulVec` — fresh. Fully general, no zeta/RH dependence
beyond `L≠0`/`∀k,η_k≠0`. -/
theorem S_T_off_sum_row_sum_vanishes (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (m : Fin (2 * (2 * K) + 1)) :
    ∑ n : Fin (2 * (2 * K) + 1), S_T_off_sum (2 * K) L K z m n * (e_off K L z n : ℂ) = 0 := by
  have hterm : ∀ n : Fin (2 * (2 * K) + 1),
      S_T_off_sum (2 * K) L K z m n * (e_off K L z n : ℂ)
        = ∑ k : Fin K, (H_L_G_C (2 * K) L (z k) m n * (e_off K L z n : ℂ)
            + H_L_G_C (2 * K) L (starRingEnd ℂ (z k)) m n * (e_off K L z n : ℂ)) := by
    intro n
    unfold S_T_off_sum
    rw [Matrix.sum_apply, Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro k _
    rw [Matrix.add_apply]
    ring
  rw [Finset.sum_congr rfl (fun n _ => hterm n), Finset.sum_comm]
  apply Finset.sum_eq_zero
  intro k _
  rw [Finset.sum_add_distrib, H_L_G_C_row_sum_vanishes K L hL z k (hz k) m,
    H_L_G_C_row_sum_vanishes_conj K L hL z k (hz k) m]
  ring

theorem S_T_off_real_mulVec_e_off_eq_zero (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (m : Fin (2 * (2 * K) + 1)) :
    (S_T_off_real (2 * K) L K z).mulVec (e_off K L z) m = 0 := by
  have hC := S_T_off_sum_row_sum_vanishes K L hL z hz m
  have hcast : (((S_T_off_real (2 * K) L K z).mulVec (e_off K L z) m : ℝ) : ℂ)
      = ∑ n : Fin (2 * (2 * K) + 1), S_T_off_sum (2 * K) L K z m n * (e_off K L z n : ℂ) := by
    unfold Matrix.mulVec dotProduct S_T_off_real
    push_cast
    apply Finset.sum_congr rfl
    intro n _
    rw [← S_T_off_sum_eq_ofReal_re]
  rw [hC] at hcast
  exact_mod_cast hcast

#print axioms S_T_off_sum_row_sum_vanishes
#print axioms S_T_off_real_mulVec_e_off_eq_zero
