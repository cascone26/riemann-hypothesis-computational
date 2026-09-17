import Mathlib.NumberTheory.Chebyshev
import Mathlib.Analysis.SpecialFunctions.Pow.Complex
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics
import Mathlib.Order.Filter.AtTopBot.Archimedean

set_option maxHeartbeats 1000000

/-!
# §1585 — boundary term vanishes: `N^{-s}·ψ(N) → 0` for `Re s > 1`

`‖N^{-s}ψ(N)‖ = N^{-σ}ψ(N) ≤ (log 4 + 4)·N^{1-σ} → 0` (Chebyshev's bound `psi_le_const_mul_self` from Mathlib +
`tendsto_rpow_neg_atTop`).  Kills the boundary term of the truncated Mellin representation §1584 as `N → ∞`.
EXIT0 / std-axioms.
-/

open Complex Chebyshev Filter

theorem boundary_term_tendsto_zero {s : ℂ} (hs : 1 < s.re) :
    Tendsto (fun N : ℕ => (N : ℂ) ^ (-s) * ((ψ (N : ℝ) : ℝ) : ℂ)) atTop (nhds 0) := by
  apply squeeze_zero_norm' (a := fun N : ℕ => (Real.log 4 + 4) * (N : ℝ) ^ (1 - s.re))
  · filter_upwards [eventually_ge_atTop 1] with N hN
    have hN0 : 0 < N := hN
    have hNR : (0 : ℝ) < (N : ℝ) := by exact_mod_cast hN0
    rw [norm_mul, Complex.norm_natCast_cpow_of_pos hN0, Complex.norm_real, Real.norm_eq_abs,
      abs_of_nonneg (psi_nonneg _), Complex.neg_re]
    have hψ : ψ (N : ℝ) ≤ (Real.log 4 + 4) * (N : ℝ) :=
      psi_le_const_mul_self (by positivity)
    calc (N : ℝ) ^ (-s.re) * ψ (N : ℝ)
        ≤ (N : ℝ) ^ (-s.re) * ((Real.log 4 + 4) * (N : ℝ)) :=
          mul_le_mul_of_nonneg_left hψ (Real.rpow_nonneg hNR.le _)
      _ = (Real.log 4 + 4) * (N : ℝ) ^ (1 - s.re) := by
          rw [show (1 : ℝ) - s.re = -s.re + 1 by ring, Real.rpow_add hNR, Real.rpow_one]
          ring
  · have h : Tendsto (fun N : ℕ => (N : ℝ) ^ (1 - s.re)) atTop (nhds 0) := by
      have h2 := (tendsto_rpow_neg_atTop (by linarith : (0 : ℝ) < s.re - 1)).comp
        tendsto_natCast_atTop_atTop
      have heq : (fun N : ℕ => (N : ℝ) ^ (1 - s.re)) = fun N : ℕ => (N : ℝ) ^ (-(s.re - 1)) := by
        funext N; congr 1; ring
      rw [heq]
      exact h2
    simpa using h.const_mul (Real.log 4 + 4)

#print axioms boundary_term_tendsto_zero
