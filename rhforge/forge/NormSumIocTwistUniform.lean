import NormSumIocTwistByParts
import NormCpowNegSuccSub
import SumIocRpowNegLe

open Finset

/-!
# §1519 — uniform twisted block bound: `‖∑_{Ioc m n} i^{-s}·χ(↑i)‖ ≤ M·(n^{-σ}+(m+1)^{-σ}+‖s‖m^{-σ}/σ)`

Combining the Abel by-parts bound §1515, the difference bound §1516 (each `‖(i+1)^{-s}-i^{-s}‖≤‖s‖i^{-(σ+1)}`), and
the summed bound §1518 (`∑_{Ioc m (n-1)} i^{-(σ+1)}≤m^{-σ}/σ`) gives the uniform block estimate — independent of `n`
except through the vanishing `n^{-σ}` term.  This is the twisted-tail control for the truncated `L(χ)` rep and its
explicit growth bound.  For `0<σ`, `s≠0`, `1≤m<n`.  EXIT0 / std-axioms.
-/

theorem norm_sum_Ioc_twist_uniform_le {M : ℕ} [NeZero M] {χ : DirichletCharacter ℂ M} (hχ : χ ≠ 1)
    {s : ℂ} (hσ : 0 < s.re) (hs0 : s ≠ 0) {m n : ℕ} (hm : 1 ≤ m) (hmn : m < n) :
    ‖∑ i ∈ Finset.Ioc m n, (i : ℂ) ^ (-s) • χ (i : ZMod M)‖
      ≤ (n : ℝ) ^ (-s.re) * M + ((m + 1 : ℕ) : ℝ) ^ (-s.re) * M
        + ‖s‖ * (m : ℝ) ^ (-s.re) / s.re * M := by
  have hby := norm_sum_Ioc_twist_by_parts_le hχ (s := s) hmn
  have hnpos : (0 : ℝ) < (n : ℝ) := by exact_mod_cast (by omega : 0 < n)
  have hm1pos : (0 : ℝ) < ((m + 1 : ℕ) : ℝ) := by exact_mod_cast (by omega : 0 < m + 1)
  have hn : ‖(n : ℂ) ^ (-s)‖ = (n : ℝ) ^ (-s.re) := by
    rw [Complex.norm_natCast_cpow_of_pos (by omega : 0 < n), Complex.neg_re]
  have hmm : ‖((m + 1 : ℕ) : ℂ) ^ (-s)‖ = ((m + 1 : ℕ) : ℝ) ^ (-s.re) := by
    rw [Complex.norm_natCast_cpow_of_pos (by omega : 0 < m + 1), Complex.neg_re]
  rw [hn, hmm] at hby
  have htel : ∑ i ∈ Finset.Ioc m (n - 1), ‖((i + 1 : ℕ) : ℂ) ^ (-s) - (i : ℂ) ^ (-s)‖ * (M : ℝ)
      ≤ ‖s‖ * (m : ℝ) ^ (-s.re) / s.re * M := by
    have hstep : ∑ i ∈ Finset.Ioc m (n - 1), ‖((i + 1 : ℕ) : ℂ) ^ (-s) - (i : ℂ) ^ (-s)‖ * (M : ℝ)
        ≤ ∑ i ∈ Finset.Ioc m (n - 1), (‖s‖ * (i : ℝ) ^ (-(s.re + 1))) * (M : ℝ) := by
      refine Finset.sum_le_sum (fun i hi => ?_)
      rw [Finset.mem_Ioc] at hi
      exact mul_le_mul_of_nonneg_right (norm_cpow_neg_succ_sub_le hσ.le hs0 (by omega)) (by positivity)
    have hfact : ∑ i ∈ Finset.Ioc m (n - 1), (‖s‖ * (i : ℝ) ^ (-(s.re + 1))) * (M : ℝ)
        = ‖s‖ * (M : ℝ) * ∑ i ∈ Finset.Ioc m (n - 1), (i : ℝ) ^ (-(s.re + 1)) := by
      rw [Finset.mul_sum]; apply Finset.sum_congr rfl; intro i _; ring
    have hsum := sum_Ioc_rpow_neg_le hσ (m := m) (k := n - 1) hm (by omega)
    have hsM : ‖s‖ * (M : ℝ) * ∑ i ∈ Finset.Ioc m (n - 1), (i : ℝ) ^ (-(s.re + 1))
        ≤ ‖s‖ * (M : ℝ) * ((m : ℝ) ^ (-s.re) / s.re) := by
      apply mul_le_mul_of_nonneg_left hsum (by positivity)
    calc ∑ i ∈ Finset.Ioc m (n - 1), ‖((i + 1 : ℕ) : ℂ) ^ (-s) - (i : ℂ) ^ (-s)‖ * (M : ℝ)
        ≤ ‖s‖ * (M : ℝ) * ∑ i ∈ Finset.Ioc m (n - 1), (i : ℝ) ^ (-(s.re + 1)) := by
          rw [← hfact]; exact hstep
      _ ≤ ‖s‖ * (M : ℝ) * ((m : ℝ) ^ (-s.re) / s.re) := hsM
      _ = ‖s‖ * (m : ℝ) ^ (-s.re) / s.re * M := by ring
  linarith [hby, htel]

#print axioms norm_sum_Ioc_twist_uniform_le
