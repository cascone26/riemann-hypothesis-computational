import Psi1CDeriv
import BoundaryTendstoZero

set_option maxHeartbeats 1000000

/-!
# §1595 — smoothed boundary dies: `Ψ(N)·N^{-s-1} → 0` for `Re s > 1`

`‖Ψ(N)‖ = ‖∫_1^N ψ‖ ≤ (N-1)·ψ(N) ≤ (log 4+4)·N²` (`norm_integral_le_of_norm_le_const` + Chebyshev), so
`‖Ψ(N)N^{-s-1}‖ ≤ (log 4+4)·N^{1-σ} → 0` — the §1585 squeeze pattern.  Kills the boundary of §1594 as `N → ∞`.
EXIT0 / std-axioms.
-/

open Complex Chebyshev Filter

theorem psi1C_norm_le {N : ℕ} (hN : 1 ≤ N) :
    ‖psi1C (N : ℝ)‖ ≤ (Real.log 4 + 4) * (N : ℝ) ^ (2 : ℝ) := by
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hb : ∀ t ∈ Set.uIcc (1 : ℝ) (N : ℝ), ‖((ψ t : ℝ) : ℂ)‖ ≤ (Real.log 4 + 4) * (N : ℝ) := by
    intro t ht
    rw [Set.uIcc_of_le hNR, Set.mem_Icc] at ht
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (psi_nonneg t)]
    calc ψ t ≤ ψ (N : ℝ) := psi_mono ht.2
      _ ≤ (Real.log 4 + 4) * (N : ℝ) := psi_le_const_mul_self (by linarith)
  have h : ‖∫ t in (1:ℝ)..(N:ℝ), ((ψ t : ℝ) : ℂ)‖ ≤ (Real.log 4 + 4) * (N : ℝ) * |(N : ℝ) - 1| := by
    apply intervalIntegral.norm_integral_le_of_norm_le_const
    intro t ht
    rw [Set.uIoc_of_le hNR, Set.mem_Ioc] at ht
    exact hb t (by rw [Set.uIcc_of_le hNR]; exact Set.mem_Icc.mpr ⟨le_of_lt ht.1, ht.2⟩)
  calc ‖psi1C (N : ℝ)‖
      ≤ (Real.log 4 + 4) * (N : ℝ) * |(N : ℝ) - 1| := h
    _ ≤ (Real.log 4 + 4) * (N : ℝ) * (N : ℝ) := by
        have habs : |(N : ℝ) - 1| ≤ (N : ℝ) := by
          rw [abs_of_nonneg (by linarith)]; linarith
        have hc : (0 : ℝ) ≤ (Real.log 4 + 4) * (N : ℝ) := by
          have := Real.log_nonneg (by norm_num : (1:ℝ) ≤ 4)
          positivity
        exact mul_le_mul_of_nonneg_left habs hc
    _ = (Real.log 4 + 4) * (N : ℝ) ^ (2 : ℝ) := by
        rw [show ((2:ℝ)) = ((2:ℕ):ℝ) by norm_num, Real.rpow_natCast]
        ring

theorem psi1C_boundary_tendsto_zero {s : ℂ} (hs : 1 < s.re) :
    Tendsto (fun N : ℕ => psi1C (N : ℝ) * (N : ℂ) ^ (-s - 1)) atTop (nhds 0) := by
  apply squeeze_zero_norm' (a := fun N : ℕ => (Real.log 4 + 4) * (N : ℝ) ^ (1 - s.re))
  · filter_upwards [eventually_ge_atTop 1] with N hN
    have hN0 : 0 < N := hN
    have hNR : (0 : ℝ) < (N : ℝ) := by exact_mod_cast hN0
    have hNR1 : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
    rw [norm_mul, Complex.norm_natCast_cpow_of_pos hN0]
    have hexp : ((-s - 1 : ℂ)).re = -s.re - 1 := by
      simp [Complex.sub_re, Complex.neg_re]
    rw [hexp]
    calc ‖psi1C (N : ℝ)‖ * (N : ℝ) ^ (-s.re - 1)
        ≤ ((Real.log 4 + 4) * (N : ℝ) ^ (2 : ℝ)) * (N : ℝ) ^ (-s.re - 1) := by
          apply mul_le_mul_of_nonneg_right (psi1C_norm_le hN) (Real.rpow_nonneg hNR.le _)
      _ = (Real.log 4 + 4) * (N : ℝ) ^ (1 - s.re) := by
          rw [mul_assoc, ← Real.rpow_add hNR]
          congr 2
          ring
  · have h : Tendsto (fun N : ℕ => (N : ℝ) ^ (1 - s.re)) atTop (nhds 0) := by
      have h2 := (tendsto_rpow_neg_atTop (by linarith : (0 : ℝ) < s.re - 1)).comp
        tendsto_natCast_atTop_atTop
      have heq : (fun N : ℕ => (N : ℝ) ^ (1 - s.re)) = fun N : ℕ => (N : ℝ) ^ (-(s.re - 1)) := by
        funext N; congr 1; ring
      rw [heq]
      exact h2
    simpa using h.const_mul (Real.log 4 + 4)

#print axioms psi1C_boundary_tendsto_zero
