import SumVonMangoldtByParts
import SumStripIntegrals

set_option maxHeartbeats 1000000

/-!
# §1584 ★★ — the TRUNCATED MELLIN representation: `∑_{n≤N} Λ(n)n^{-s} = N^{-s}ψ(N) + s·∫_1^N ψ(t)t^{-s-1} dt`

Seventh Mellin brick — the discrete half of the pillar is COMPLETE.  §1578 (Abel) writes the truncated
`-ζ'/ζ` sum through `ψ`; the difference sum re-indexes from `Ioc 0 (N-1)` to `Ico 1 N` (they coincide) and §1583
turns it into `-(-s∫_1^N) = +s∫_1^N`.  What remains for the full pillar is the `N → ∞` limit (both sides) and
`logDeriv` identification.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction Chebyshev

theorem truncated_mellin {s : ℂ} (hs : s ≠ 0) (hσ : 0 ≤ s.re) {N : ℕ} (hN : 0 < N) :
    ∑ n ∈ Finset.Ioc 0 N, (Λ n : ℂ) * (n : ℂ) ^ (-s)
      = (N : ℂ) ^ (-s) * ((ψ (N : ℝ) : ℝ) : ℂ)
        + s * ∫ t : ℝ in (1 : ℝ)..(N : ℝ), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1) := by
  rw [sum_vonMangoldt_cpow_by_parts hN]
  have hidx : Finset.Ioc 0 (N - 1) = Finset.Ico 1 N := by
    ext x
    simp only [Finset.mem_Ioc, Finset.mem_Ico]
    omega
  rw [hidx, sum_strip_integrals hs hσ hN]
  ring

#print axioms truncated_mellin
