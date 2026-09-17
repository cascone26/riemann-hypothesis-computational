import NormSumIocTwistUnifRegion
import TendstoNatRpowNegBound
import Mathlib.NumberTheory.LSeries.Dirichlet
import Mathlib.Topology.UniformSpace.UniformConvergence

open Finset Filter Topology

/-! §1530 UniformCauchySeqOn of the twisted partial sums on {δ ≤ Re s, ‖s‖ ≤ R}. EXIT0 / std-axioms. -/

theorem uniformCauchySeqOn_twist_partial_sum {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M}
    (hχ : χ ≠ 1) {δ R : ℝ} (hδ : 0 < δ) :
    UniformCauchySeqOn (fun n s => ∑ k ∈ Finset.Icc (1 : ℕ) n, χ (k : ZMod M) * (k : ℂ) ^ (-s))
      atTop {s : ℂ | δ ≤ s.re ∧ ‖s‖ ≤ R} := by
  rw [Metric.uniformCauchySeqOn_iff]
  intro ε hε
  set C : ℝ := (M : ℝ) * (2 + R / δ) with hC
  have hb' : Tendsto (fun N : ℕ => (N : ℝ) ^ (-δ) * C) atTop (𝓝 0) := by
    have := tendsto_nat_rpow_neg_bound hδ C 0; simpa using this
  obtain ⟨N, hN⟩ := (Metric.tendsto_atTop.mp hb') ε hε
  refine ⟨max N 1, fun p hp q hq s hs => ?_⟩
  obtain ⟨hsδ, hsR⟩ := hs
  have hIcc : ∀ x : ℕ, Finset.Icc (1 : ℕ) x = Finset.Ioc 0 x := fun x => by
    ext y; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  have hMnn : (0 : ℝ) ≤ (M : ℝ) := Nat.cast_nonneg M
  have hCnn : (0 : ℝ) ≤ C := by rw [hC]; have : 0 ≤ R / δ := div_nonneg (le_trans (norm_nonneg s) hsR) hδ.le; positivity
  have hkey : ∀ a b : ℕ, 1 ≤ a → a ≤ b → N ≤ a →
      dist (∑ k ∈ Finset.Icc (1 : ℕ) b, χ (k : ZMod M) * (k : ℂ) ^ (-s))
        (∑ k ∈ Finset.Icc (1 : ℕ) a, χ (k : ZMod M) * (k : ℂ) ^ (-s)) < ε := by
    intro a b ha hab haN
    rcases eq_or_lt_of_le hab with h | h
    · subst h; simpa using hε
    · have hdiff : (∑ k ∈ Finset.Icc (1 : ℕ) b, χ (k : ZMod M) * (k : ℂ) ^ (-s))
          - (∑ k ∈ Finset.Icc (1 : ℕ) a, χ (k : ZMod M) * (k : ℂ) ^ (-s))
          = ∑ k ∈ Finset.Ioc a b, (k : ℂ) ^ (-s) • χ (k : ZMod M) := by
        rw [hIcc, hIcc, ← Finset.sum_Ioc_consecutive _ (Nat.zero_le a) hab, add_sub_cancel_left]
        apply Finset.sum_congr rfl; intro k _; rw [smul_eq_mul, mul_comm]
      rw [dist_eq_norm, hdiff]
      have haR : (1 : ℝ) ≤ (a : ℝ) := by exact_mod_cast ha
      have habR : (a : ℝ) ≤ (b : ℝ) := by exact_mod_cast hab
      have haσnn : (0 : ℝ) ≤ (a : ℝ) ^ (-δ) := Real.rpow_nonneg (by linarith) _
      have hu := norm_sum_Ioc_twist_unif_region hχ hδ hsδ hsR ha h
      have hbσ : (b : ℝ) ^ (-δ) ≤ (a : ℝ) ^ (-δ) :=
        Real.rpow_le_rpow_of_nonpos (by linarith) habR (by linarith)
      have ha1σ : ((a + 1 : ℕ) : ℝ) ^ (-δ) ≤ (a : ℝ) ^ (-δ) :=
        Real.rpow_le_rpow_of_nonpos (by linarith) (by push_cast; linarith) (by linarith)
      have hstep : ((b : ℝ) ^ (-δ) + ((a + 1 : ℕ) : ℝ) ^ (-δ) + R * (a : ℝ) ^ (-δ) / δ) * M
          ≤ (a : ℝ) ^ (-δ) * C := by
        have hCexp : (a : ℝ) ^ (-δ) * C
            = ((a : ℝ) ^ (-δ) + (a : ℝ) ^ (-δ) + R * (a : ℝ) ^ (-δ) / δ) * M := by
          rw [hC]; field_simp; ring
        rw [hCexp]
        gcongr
      have hlt : (a : ℝ) ^ (-δ) * C < ε := by
        have hthis := hN a haN
        rw [Real.dist_eq, sub_zero, abs_of_nonneg (by positivity)] at hthis
        exact hthis
      calc ‖∑ k ∈ Finset.Ioc a b, (k : ℂ) ^ (-s) • χ (k : ZMod M)‖
          ≤ ((b : ℝ) ^ (-δ) + ((a + 1 : ℕ) : ℝ) ^ (-δ) + R * (a : ℝ) ^ (-δ) / δ) * M := hu
        _ ≤ (a : ℝ) ^ (-δ) * C := hstep
        _ < ε := hlt
  rcases le_total p q with hpq | hpq
  · rw [dist_comm]; exact hkey p q (le_trans (le_max_right N 1) hp) hpq (le_trans (le_max_left N 1) hp)
  · exact hkey q p (le_trans (le_max_right N 1) hq) hpq (le_trans (le_max_left N 1) hq)

#print axioms uniformCauchySeqOn_twist_partial_sum
