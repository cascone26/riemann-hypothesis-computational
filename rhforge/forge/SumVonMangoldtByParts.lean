import Mathlib.NumberTheory.Chebyshev
import Mathlib.Algebra.BigOperators.Module
import Mathlib.Analysis.SpecialFunctions.Pow.Complex

set_option maxHeartbeats 1000000

/-!
# §1578 — Abel step of the Mellin pillar: `∑_{n≤N} Λ(n)n^{-s} = N^{-s}ψ(N) - ∑ (Δn^{-s})·ψ(n)`

Opening brick of the PNT-ERROR arc.  Summation by parts (`Finset.sum_Ioc_by_parts`, the §1515 idiom) with
`f i = i^{-s}`, `g i = Λ(i)`: the partial sums of `g` are exactly Chebyshev's `ψ` (`psi_eq_sum_Icc`,
`range (N+1) = Icc 0 N`, `Λ 0 = 0`).  This expresses the truncated `-ζ'/ζ` Dirichlet sum through `ψ` — the
discrete half of the Mellin representation `-ζ'/ζ(s) = s∫₁^∞ ψ(t)t^{-s-1}dt`.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction Chebyshev Finset

theorem sum_vonMangoldt_cpow_by_parts {s : ℂ} {N : ℕ} (hN : 0 < N) :
    ∑ n ∈ Finset.Ioc 0 N, (Λ n : ℂ) * (n : ℂ) ^ (-s)
      = (N : ℂ) ^ (-s) * ((ψ (N : ℝ) : ℝ) : ℂ)
        - ∑ i ∈ Finset.Ioc 0 (N - 1),
            ((((i + 1 : ℕ) : ℂ)) ^ (-s) - (i : ℂ) ^ (-s)) * ((ψ (i : ℝ) : ℝ) : ℂ) := by
  have hrange : ∀ k : ℕ, ∑ i ∈ Finset.range (k + 1), ((Λ i : ℝ) : ℂ) = ((ψ (k : ℝ) : ℝ) : ℂ) := by
    intro k
    rw [Nat.range_succ_eq_Icc_zero, psi_eq_sum_Icc, Nat.floor_natCast]
    push_cast
    ring
  have hcomm : ∀ n : ℕ, (Λ n : ℂ) * (n : ℂ) ^ (-s) = (n : ℂ) ^ (-s) • ((Λ n : ℝ) : ℂ) := by
    intro n
    rw [smul_eq_mul, mul_comm]
  calc ∑ n ∈ Finset.Ioc 0 N, (Λ n : ℂ) * (n : ℂ) ^ (-s)
      = ∑ n ∈ Finset.Ioc 0 N, ((n : ℂ) ^ (-s)) • ((Λ n : ℝ) : ℂ) := by
        exact Finset.sum_congr rfl (fun n _ => hcomm n)
    _ = ((N : ℂ) ^ (-s)) • ∑ i ∈ Finset.range (N + 1), ((Λ i : ℝ) : ℂ)
        - (((0 + 1 : ℕ) : ℂ) ^ (-s)) • ∑ i ∈ Finset.range (0 + 1), ((Λ i : ℝ) : ℂ)
        - ∑ i ∈ Finset.Ioc 0 (N - 1),
            ((((i + 1 : ℕ) : ℂ)) ^ (-s) - (i : ℂ) ^ (-s)) • ∑ j ∈ Finset.range (i + 1), ((Λ j : ℝ) : ℂ) := by
        exact Finset.sum_Ioc_by_parts (fun i => (i : ℂ) ^ (-s)) (fun i => ((Λ i : ℝ) : ℂ)) hN
    _ = (N : ℂ) ^ (-s) * ((ψ (N : ℝ) : ℝ) : ℂ)
        - ∑ i ∈ Finset.Ioc 0 (N - 1),
            ((((i + 1 : ℕ) : ℂ)) ^ (-s) - (i : ℂ) ^ (-s)) * ((ψ (i : ℝ) : ℝ) : ℂ) := by
        rw [hrange N]
        have hzero : ∑ i ∈ Finset.range (0 + 1), ((Λ i : ℝ) : ℂ) = 0 := by
          simp [ArithmeticFunction.map_zero]
        rw [hzero, smul_zero, sub_zero, smul_eq_mul]
        congr 1
        exact Finset.sum_congr rfl (fun i _ => by rw [hrange i, smul_eq_mul])

#print axioms sum_vonMangoldt_cpow_by_parts
