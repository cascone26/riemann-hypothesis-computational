import StripDiffPsiEqIntegral
import PsiCpowIntervalIntegrable

set_option maxHeartbeats 1000000

/-!
# §1583 — the strips glue: `∑_{k=1}^{N-1} Δk^{-s}·ψ(k) = -s·∫_1^N ψ(t)t^{-s-1} dt`

Sixth Mellin brick.  Each Abel-difference term is a strip integral (§1581); the strips `[k, k+1]` for
`k ∈ Ico 1 N` tile `[1, N]`, so `intervalIntegral.sum_integral_adjacent_intervals_Ico` (with the §1582
integrability on each strip) glues them into a single integral.  With §1578 this gives the truncated Mellin
representation `∑_{n≤N} Λ(n)n^{-s} = N^{-s}ψ(N) + s∫_1^N ψ(t)t^{-s-1}dt`.  EXIT0 / std-axioms.
-/

open Complex Chebyshev

theorem sum_strip_integrals {s : ℂ} (hs : s ≠ 0) (hσ : 0 ≤ s.re) {N : ℕ} (hN : 0 < N) :
    ∑ k ∈ Finset.Ico 1 N, ((((k + 1 : ℕ) : ℂ)) ^ (-s) - (k : ℂ) ^ (-s)) * ((ψ (k : ℝ) : ℝ) : ℂ)
      = -s * ∫ t : ℝ in (1 : ℝ)..(N : ℝ), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1) := by
  have hstep : ∀ k ∈ Finset.Ico 1 N,
      ((((k + 1 : ℕ) : ℂ)) ^ (-s) - (k : ℂ) ^ (-s)) * ((ψ (k : ℝ) : ℝ) : ℂ)
      = -s * ∫ t : ℝ in ((k : ℕ) : ℝ)..(((k + 1 : ℕ)) : ℝ), ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1) := by
    intro k hk
    rw [Finset.mem_Ico] at hk
    have hcast : (((k + 1 : ℕ)) : ℝ) = (k : ℝ) + 1 := by push_cast; ring
    have h81 := strip_diff_psi_eq_integral hs hk.1
    rw [← hcast] at h81
    exact h81
  rw [Finset.sum_congr rfl hstep, ← Finset.mul_sum]
  congr 1
  have hint : ∀ k ∈ Finset.Ico 1 N, IntervalIntegrable
      (fun t : ℝ => ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1)) MeasureTheory.volume
      (((k : ℕ)) : ℝ) (((k + 1 : ℕ)) : ℝ) := by
    intro k hk
    rw [Finset.mem_Ico] at hk
    have hcast : (((k + 1 : ℕ)) : ℝ) = (k : ℝ) + 1 := by push_cast; ring
    have h82 := psi_cpow_intervalIntegrable hσ hk.1
    rw [← hcast] at h82
    exact h82
  have hintS : ∀ k ∈ Set.Ico 1 N, IntervalIntegrable
      (fun t : ℝ => ((ψ t : ℝ) : ℂ) * (t : ℂ) ^ (-s - 1)) MeasureTheory.volume
      (((k : ℕ)) : ℝ) (((k + 1 : ℕ)) : ℝ) := by
    intro k hk
    exact hint k (Finset.mem_Ico.mpr ⟨hk.1, hk.2⟩)
  have hadj := intervalIntegral.sum_integral_adjacent_intervals_Ico
    (a := fun k : ℕ => (k : ℝ)) hN hintS
  simpa using hadj

#print axioms sum_strip_integrals
