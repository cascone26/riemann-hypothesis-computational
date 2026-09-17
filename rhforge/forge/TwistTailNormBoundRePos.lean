import NormSumIocTwistUniform
import TendstoNatRpowNegBound
import Mathlib.NumberTheory.LSeries.Dirichlet

open Finset Filter Topology

/-! §1527 twisted-tail norm bound on Re s > 0 (conditional limit given). EXIT0 / std-axioms. -/

theorem twist_tail_norm_bound_re_pos {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {s : ℂ} (hσ : 0 < s.re) {N : ℕ} (hN : 1 ≤ N) {g : ℂ}
    (hg : Tendsto (fun n => ∑ k ∈ Finset.Icc (1 : ℕ) n, χ (k : ZMod M) * (k : ℂ) ^ (-s))
      atTop (𝓝 g)) :
    ‖g - ∑ k ∈ Finset.Icc (1 : ℕ) N, χ (k : ZMod M) * (k : ℂ) ^ (-s)‖
      ≤ ((N + 1 : ℕ) : ℝ) ^ (-s.re) * M + ‖s‖ * (N : ℝ) ^ (-s.re) / s.re * M := by
  have hs0 : s ≠ 0 := by rintro rfl; rw [Complex.zero_re] at hσ; exact lt_irrefl 0 hσ
  set P := ∑ k ∈ Finset.Icc (1 : ℕ) N, χ (k : ZMod M) * (k : ℂ) ^ (-s) with hP
  set D : ℝ := ((N + 1 : ℕ) : ℝ) ^ (-s.re) * M + ‖s‖ * (N : ℝ) ^ (-s.re) / s.re * M with hD
  have hDnn : 0 ≤ D := by rw [hD]; positivity
  have hIcc : ∀ x : ℕ, Finset.Icc (1 : ℕ) x = Finset.Ioc 0 x := fun x => by
    ext y; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  have htend : Tendsto (fun n => ∑ k ∈ Finset.Ioc N n, (k : ℂ) ^ (-s) • χ (k : ZMod M))
      atTop (𝓝 (g - P)) := by
    refine Tendsto.congr' ?_ (hg.sub_const P)
    filter_upwards [eventually_ge_atTop N] with n hn
    rw [hP, hIcc, hIcc, ← Finset.sum_Ioc_consecutive _ (Nat.zero_le N) hn, add_sub_cancel_left]
    apply Finset.sum_congr rfl; intro k _; rw [smul_eq_mul, mul_comm]
  have hnorm : Tendsto (fun n => ‖∑ k ∈ Finset.Ioc N n, (k : ℂ) ^ (-s) • χ (k : ZMod M)‖)
      atTop (𝓝 ‖g - P‖) := htend.norm
  have hb : Tendsto (fun n : ℕ => (n : ℝ) ^ (-s.re) * M + D) atTop (𝓝 D) := by
    have := tendsto_nat_rpow_neg_bound hσ (M : ℝ) D; simpa using this
  have hab : ∀ n, ‖∑ k ∈ Finset.Ioc N n, (k : ℂ) ^ (-s) • χ (k : ZMod M)‖
      ≤ (n : ℝ) ^ (-s.re) * M + D := by
    intro n
    by_cases hn : N < n
    · have h := norm_sum_Ioc_twist_uniform_le hχ hσ hs0 hN hn
      rw [hD]; linarith [h]
    · have hempty : Finset.Ioc N n = ∅ := by rw [Finset.Ioc_eq_empty]; omega
      rw [hempty, Finset.sum_empty, norm_zero]
      have : (0 : ℝ) ≤ (n : ℝ) ^ (-s.re) * M := by positivity
      linarith [hDnn]
  exact le_of_tendsto_of_tendsto' hnorm hb hab

#print axioms twist_tail_norm_bound_re_pos
