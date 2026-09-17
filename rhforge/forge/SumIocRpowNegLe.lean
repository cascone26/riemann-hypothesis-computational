import SumIcoSuccRpowNegLe

open Finset

/-!
# §1518 — Ioc form of the summed total-variation bound: `∑_{i∈Ioc m k} i^{-(σ+1)} ≤ m^{-σ}/σ`

Reindex `Ioc m k → Ico m k` by `i ↦ i-1` reduces to §1517.  The form the Abel telescoping sum (§1515, over
`Ioc m (n-1)`) directly consumes.  EXIT0 / std-axioms.
-/

theorem sum_Ioc_rpow_neg_le {σ : ℝ} (hσ : 0 < σ) {m k : ℕ} (hm : 1 ≤ m) (hmk : m ≤ k) :
    ∑ i ∈ Finset.Ioc m k, (i : ℝ) ^ (-(σ + 1)) ≤ (m : ℝ) ^ (-σ) / σ := by
  have hreindex : ∑ i ∈ Finset.Ioc m k, (i : ℝ) ^ (-(σ + 1))
      = ∑ i ∈ Finset.Ico m k, ((i : ℝ) + 1) ^ (-(σ + 1)) := by
    refine Finset.sum_nbij' (fun i => i - 1) (fun i => i + 1) ?_ ?_ ?_ ?_ ?_
    · intro a ha; rw [Finset.mem_Ioc] at ha; rw [Finset.mem_Ico]; omega
    · intro b hb; rw [Finset.mem_Ico] at hb; rw [Finset.mem_Ioc]; omega
    · intro a ha; rw [Finset.mem_Ioc] at ha; omega
    · intro b _; omega
    · intro a ha
      rw [Finset.mem_Ioc] at ha
      have hcast : ((a - 1 : ℕ) : ℝ) + 1 = (a : ℝ) := by
        rw [Nat.cast_sub (by omega), Nat.cast_one]; ring
      rw [hcast]
  rw [hreindex]
  exact sum_Ico_succ_rpow_neg_le hσ hm hmk

#print axioms sum_Ioc_rpow_neg_le
