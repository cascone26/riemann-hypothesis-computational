import SumIcoCharEqZero

open Finset

/-!
# §1514 — bounded character partial sums: `‖∑_{k ∈ range n} χ(↑k)‖ ≤ M` (`χ ≠ 1`, all `n`)

Strong induction: for `n < M` the sum of `n` terms each of norm `≤ 1` is `≤ n ≤ M`; for `n ≥ M` a full block
`Ico (n-M) n` sums to `0` (§1513), so `S(n) = S(n-M)` and the IH applies.  This uniform bound on the twisted
partial sums is the Dirichlet-test / Abel-summation input for the truncated `L(χ)` representation.  EXIT0 / std.
-/

theorem norm_partial_sum_char_le {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    (n : ℕ) : ‖∑ k ∈ Finset.range n, χ (k : ZMod M)‖ ≤ M := by
  have hMpos : 0 < M := Nat.pos_of_ne_zero (NeZero.ne M)
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    rcases lt_or_ge n M with hn | hn
    · calc ‖∑ k ∈ range n, χ (k : ZMod M)‖
          ≤ ∑ k ∈ range n, ‖χ (k : ZMod M)‖ := norm_sum_le _ _
        _ ≤ ∑ _k ∈ range n, (1 : ℝ) := Finset.sum_le_sum (fun k _ => χ.norm_le_one _)
        _ = (n : ℝ) := by simp
        _ ≤ (M : ℝ) := by exact_mod_cast hn.le
    · have hnM : (n - M) + M = n := Nat.sub_add_cancel hn
      have h0 : ∑ k ∈ Finset.Ico (n - M) n, χ (k : ZMod M) = 0 := by
        have h := sum_Ico_char_eq_zero hχ (n - M)
        rwa [hnM] at h
      have hsplit : ∑ k ∈ range n, χ (k : ZMod M) = ∑ k ∈ range (n - M), χ (k : ZMod M) := by
        rw [Finset.range_eq_Ico,
          ← Finset.sum_Ico_consecutive _ (Nat.zero_le (n - M)) (Nat.sub_le n M), h0,
          add_zero, ← Finset.range_eq_Ico]
      rw [hsplit]
      exact ih (n - M) (by omega)

#print axioms norm_partial_sum_char_le
