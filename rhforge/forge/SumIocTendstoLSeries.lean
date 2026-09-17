import Mathlib.NumberTheory.LSeries.Dirichlet

set_option maxHeartbeats 1000000

/-!
# §1586 — the truncated sums converge to `L(Λ)`: `∑_{n≤N} Λ(n)n^{-s} → LSeries ↗Λ s` (`Re s > 1`)

The vonMangoldt L-series is summable for `Re s > 1` (`LSeriesSummable_vonMangoldt`); its partial sums over
`range (N+1)` tend to `LSeries ↗Λ s` (`HasSum.tendsto_sum_nat` + the `+1` shift), and they equal the `Ioc 0 N`
sums in the `Λ(n)·n^{-s}` form of §1584 (`term 0 = 0`; `term n = Λ(n)/n^s = Λ(n)·n^{-s}`).  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction Filter
open scoped LSeries.notation

theorem sum_Ioc_tendsto_LSeries {s : ℂ} (hs : 1 < s.re) :
    Tendsto (fun N : ℕ => ∑ n ∈ Finset.Ioc 0 N, (Λ n : ℂ) * (n : ℂ) ^ (-s))
      atTop (nhds (LSeries (↗Λ) s)) := by
  have hsum : LSeriesSummable (↗Λ) s := ArithmeticFunction.LSeriesSummable_vonMangoldt hs
  have hh : Tendsto (fun N : ℕ => ∑ n ∈ Finset.range (N + 1), LSeries.term (↗Λ) s n)
      atTop (nhds (LSeries (↗Λ) s)) :=
    (hsum.hasSum.tendsto_sum_nat).comp (tendsto_add_atTop_nat 1)
  apply hh.congr
  intro N
  rw [Nat.range_succ_eq_Icc_zero, ← Finset.add_sum_Ioc_eq_sum_Icc (Nat.zero_le N),
    LSeries.term_zero, zero_add]
  refine Finset.sum_congr rfl (fun n hn => ?_)
  rw [Finset.mem_Ioc] at hn
  rw [LSeries.term_of_ne_zero (by omega : n ≠ 0), Complex.cpow_neg, div_eq_mul_inv]

#print axioms sum_Ioc_tendsto_LSeries
