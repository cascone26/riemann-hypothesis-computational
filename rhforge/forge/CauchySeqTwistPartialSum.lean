import NormSumIocTwistUniform
import TendstoNatRpowNegBound
import Mathlib.NumberTheory.LSeries.Dirichlet

open Finset Filter Topology

/-! §1526 the twisted partial sums are Cauchy at each s with Re s > 0. EXIT0 / std-axioms. -/

theorem cauchySeq_twist_partial_sum {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {s : ℂ} (hσ : 0 < s.re) :
    CauchySeq (fun n => ∑ k ∈ Finset.Icc (1 : ℕ) n, χ (k : ZMod M) * (k : ℂ) ^ (-s)) := by
  set C : ℝ := (M : ℝ) * (2 + ‖s‖ / s.re) with hC
  have hb' : Tendsto (fun N : ℕ => (N : ℝ) ^ (-s.re) * C) atTop (𝓝 0) := by
    have := tendsto_nat_rpow_neg_bound hσ C 0; simpa using this
  rw [Metric.cauchySeq_iff]
  intro ε hε
  obtain ⟨N, hN⟩ := (Metric.tendsto_atTop.mp hb') ε hε
  refine ⟨max N 1, fun p hp q hq => ?_⟩
  have hkey : ∀ a b : ℕ, 1 ≤ a → a ≤ b → N ≤ a →
      dist (∑ k ∈ Finset.Icc (1 : ℕ) b, χ (k : ZMod M) * (k : ℂ) ^ (-s))
        (∑ k ∈ Finset.Icc (1 : ℕ) a, χ (k : ZMod M) * (k : ℂ) ^ (-s)) < ε := by
    intro a b ha hab haN
    rcases eq_or_lt_of_le hab with h | h
    · subst h; simpa using hε
    · have hIcc : ∀ x : ℕ, Finset.Icc (1 : ℕ) x = Finset.Ioc 0 x := fun x => by
        ext y; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
      have hdiff : (∑ k ∈ Finset.Icc (1 : ℕ) b, χ (k : ZMod M) * (k : ℂ) ^ (-s))
          - (∑ k ∈ Finset.Icc (1 : ℕ) a, χ (k : ZMod M) * (k : ℂ) ^ (-s))
          = ∑ k ∈ Finset.Ioc a b, (k : ℂ) ^ (-s) • χ (k : ZMod M) := by
        rw [hIcc, hIcc, ← Finset.sum_Ioc_consecutive _ (Nat.zero_le a) hab, add_sub_cancel_left]
        apply Finset.sum_congr rfl; intro k _; rw [smul_eq_mul, mul_comm]
      rw [dist_eq_norm, hdiff]
      have hσC : (0 : ℝ) ≤ C := by rw [hC]; positivity
      have haR : (1 : ℝ) ≤ (a : ℝ) := by exact_mod_cast ha
      have habR : (a : ℝ) ≤ (b : ℝ) := by exact_mod_cast hab
      have ha1R : (1 : ℝ) ≤ ((a + 1 : ℕ) : ℝ) := by exact_mod_cast (by omega : 1 ≤ a + 1)
      have hMnn : (0 : ℝ) ≤ (M : ℝ) := Nat.cast_nonneg M
      have hu := norm_sum_Ioc_twist_uniform_le hχ hσ (by rintro rfl; simp at hσ) ha h
      -- bound the §1519 RHS by a^{-σ}*C
      have hbσ : (b : ℝ) ^ (-s.re) ≤ (a : ℝ) ^ (-s.re) :=
        Real.rpow_le_rpow_of_nonpos (by linarith) habR (by linarith)
      have ha1σ : ((a + 1 : ℕ) : ℝ) ^ (-s.re) ≤ (a : ℝ) ^ (-s.re) :=
        Real.rpow_le_rpow_of_nonpos (by linarith) (by push_cast; linarith) (by linarith)
      have hstep : (b : ℝ) ^ (-s.re) * M + ((a + 1 : ℕ) : ℝ) ^ (-s.re) * M
          + ‖s‖ * (a : ℝ) ^ (-s.re) / s.re * M ≤ (a : ℝ) ^ (-s.re) * C := by
        have hCexp : (a : ℝ) ^ (-s.re) * C
            = (a : ℝ) ^ (-s.re) * M + (a : ℝ) ^ (-s.re) * M
              + ‖s‖ * (a : ℝ) ^ (-s.re) / s.re * M := by rw [hC]; field_simp; ring
        rw [hCexp]
        have t1 : (b : ℝ) ^ (-s.re) * M ≤ (a : ℝ) ^ (-s.re) * M :=
          mul_le_mul_of_nonneg_right hbσ hMnn
        have t2 : ((a + 1 : ℕ) : ℝ) ^ (-s.re) * M ≤ (a : ℝ) ^ (-s.re) * M :=
          mul_le_mul_of_nonneg_right ha1σ hMnn
        linarith [t1, t2]
      have hlt : (a : ℝ) ^ (-s.re) * C < ε := by
        have hthis := hN a haN
        rw [Real.dist_eq, sub_zero, abs_of_nonneg (by positivity)] at hthis
        exact hthis
      calc ‖∑ k ∈ Finset.Ioc a b, (k : ℂ) ^ (-s) • χ (k : ZMod M)‖
          ≤ (b : ℝ) ^ (-s.re) * M + ((a + 1 : ℕ) : ℝ) ^ (-s.re) * M
            + ‖s‖ * (a : ℝ) ^ (-s.re) / s.re * M := hu
        _ ≤ (a : ℝ) ^ (-s.re) * C := hstep
        _ < ε := hlt
  rcases le_total p q with hpq | hpq
  · rw [dist_comm]; exact hkey p q (le_trans (le_max_right N 1) hp) hpq (le_trans (le_max_left N 1) hp)
  · exact hkey q p (le_trans (le_max_right N 1) hq) hpq (le_trans (le_max_left N 1) hq)

#print axioms cauchySeq_twist_partial_sum
