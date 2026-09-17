import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv

set_option maxHeartbeats 1000000

/-!
# §1554 — the finite-product correction is bounded: `Re(-logDeriv ∏(1-p^{-·}) c) ≤ ∑_{p|N} log p`

Second brick of the QUADRATIC (`χ²=1`) case.  §1553 split
`logDeriv(LFunctionTrivChar N) = logDeriv(∏) + logDeriv ζ`.  Here we bound the finite-product term:
`logDeriv(∏_{p|N}(1-p^{-·})) c = ∑_{p|N} logDeriv(1-p^{-·}) c` (`logDeriv_prod`), each `= log(p)·p^{-c}/(1-p^{-c})`
with `‖·‖ ≤ log p` (because `p^{-σ} ≤ 1/2 ≤ ‖1-p^{-c}‖` for `σ>1`, `p≥2`).  So
`Re(-logDeriv(∏) c) ≤ ‖logDeriv(∏) c‖ ≤ ∑_{p|N} log p` — a fixed constant of `N`.  EXIT0 / std.
-/

open Complex DirichletCharacter

theorem prod_correction_re_le {N : ℕ} [NeZero N] {c : ℂ} (hc : 1 < c.re) :
    (-logDeriv (fun w => ∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w))) c).re
      ≤ ∑ p ∈ N.primeFactors, Real.log p := by
  have hper : ∀ p ∈ N.primeFactors,
      ‖logDeriv (fun w : ℂ => (1 : ℂ) - (p : ℂ) ^ (-w)) c‖ ≤ Real.log p := by
    intro p hp
    have hpp := Nat.prime_of_mem_primeFactors hp
    have hp2 : 2 ≤ p := hpp.two_le
    have hp0 : 0 < p := by omega
    have hpC0 : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
    have hp1R : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (by omega : 1 < p)
    have hp2R : (2 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp2
    have hpRnn : (0 : ℝ) ≤ (p : ℝ) := by positivity
    have hlogp_pos : 0 < Real.log p := Real.log_pos hp1R
    have hnormpc : ‖(p : ℂ) ^ (-c)‖ = (p : ℝ) ^ (-c.re) := by
      rw [Complex.norm_natCast_cpow_of_pos hp0, Complex.neg_re]
    have hpc_half : (p : ℝ) ^ (-c.re) ≤ 1 / 2 := by
      calc (p : ℝ) ^ (-c.re) ≤ (p : ℝ) ^ (-(1 : ℝ)) := by
            apply Real.rpow_le_rpow_of_exponent_le (le_of_lt hp1R); linarith
        _ = 1 / (p : ℝ) := by rw [Real.rpow_neg_one, one_div]
        _ ≤ 1 / 2 := by rw [div_le_div_iff₀ (by linarith) (by norm_num)]; linarith
    have hpc_lt1 : ‖(p : ℂ) ^ (-c)‖ < 1 := by rw [hnormpc]; linarith
    have hden_ge : (p : ℝ) ^ (-c.re) ≤ ‖(1 : ℂ) - (p : ℂ) ^ (-c)‖ := by
      have htri : (1 : ℝ) - ‖(p : ℂ) ^ (-c)‖ ≤ ‖(1 : ℂ) - (p : ℂ) ^ (-c)‖ := by
        have := norm_sub_norm_le (1 : ℂ) ((p : ℂ) ^ (-c))
        simpa using this
      rw [hnormpc] at htri; linarith
    have hfc : (1 : ℂ) - (p : ℂ) ^ (-c) ≠ 0 := by
      intro h
      have h1 : (p : ℂ) ^ (-c) = 1 := by linear_combination -h
      rw [h1, norm_one] at hpc_lt1; exact lt_irrefl 1 hpc_lt1
    have hden_pos : 0 < ‖(1 : ℂ) - (p : ℂ) ^ (-c)‖ := norm_pos_iff.mpr hfc
    -- derivative of 1 - p^{-w}
    have hd1 : HasDerivAt (fun w : ℂ => -w) (-1) c := hasDerivAt_neg' c
    have hd2 : HasDerivAt (fun w : ℂ => (p : ℂ) ^ (-w))
        ((p : ℂ) ^ (-c) * Complex.log (p : ℂ) * (-1)) c := hd1.const_cpow (Or.inl hpC0)
    have hd3 : HasDerivAt (fun w : ℂ => (1 : ℂ) - (p : ℂ) ^ (-w))
        (0 - (p : ℂ) ^ (-c) * Complex.log (p : ℂ) * (-1)) c :=
      (hasDerivAt_const c (1 : ℂ)).sub hd2
    have hlogD : logDeriv (fun w : ℂ => (1 : ℂ) - (p : ℂ) ^ (-w)) c
        = (0 - (p : ℂ) ^ (-c) * Complex.log (p : ℂ) * (-1)) / ((1 : ℂ) - (p : ℂ) ^ (-c)) := by
      rw [logDeriv_apply, hd3.deriv]
    have hlogpn : ‖Complex.log (p : ℂ)‖ = Real.log p := by
      rw [show ((p : ℕ) : ℂ) = (((p : ℕ) : ℝ) : ℂ) by push_cast; ring,
        ← Complex.ofReal_log hpRnn, Complex.norm_real, Real.norm_eq_abs,
        abs_of_pos hlogp_pos]
    have hnum : ‖(0 : ℂ) - (p : ℂ) ^ (-c) * Complex.log (p : ℂ) * (-1)‖
        = Real.log p * (p : ℝ) ^ (-c.re) := by
      rw [zero_sub, norm_neg, norm_mul, norm_mul, norm_neg, norm_one, mul_one, hnormpc, hlogpn]
      ring
    rw [hlogD, norm_div, hnum, div_le_iff₀ hden_pos]
    exact mul_le_mul_of_nonneg_left hden_ge hlogp_pos.le
  -- assemble
  have hFval : logDeriv (fun w => ∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w))) c
      = ∑ p ∈ N.primeFactors, logDeriv (fun w : ℂ => (1 : ℂ) - (p : ℂ) ^ (-w)) c := by
    apply logDeriv_prod
    · intro p hp
      have hpp := Nat.prime_of_mem_primeFactors hp
      have hp2 : 2 ≤ p := hpp.two_le
      have hp0 : 0 < p := by omega
      have hp1R : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (by omega : 1 < p)
      intro h
      have hnormpc : ‖(p : ℂ) ^ (-c)‖ = (p : ℝ) ^ (-c.re) := by
        rw [Complex.norm_natCast_cpow_of_pos hp0, Complex.neg_re]
      have h1 : (p : ℂ) ^ (-c) = 1 := by linear_combination -h
      have : ‖(p : ℂ) ^ (-c)‖ = 1 := by rw [h1, norm_one]
      rw [hnormpc] at this
      have hlt : (p : ℝ) ^ (-c.re) < 1 := by
        apply Real.rpow_lt_one_of_one_lt_of_neg hp1R; linarith
      linarith
    · intro p hp
      have hpp := Nat.prime_of_mem_primeFactors hp
      have hp2 : 2 ≤ p := hpp.two_le
      have hpC0 : (p : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr (by omega)
      have hexp : DifferentiableAt ℂ (fun w : ℂ => -w) c := differentiableAt_id.neg
      exact (differentiableAt_const 1).sub (hexp.const_cpow (Or.inl hpC0))
  calc (-logDeriv (fun w => ∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w))) c).re
      ≤ ‖-logDeriv (fun w => ∏ p ∈ N.primeFactors, (1 - (p : ℂ) ^ (-w))) c‖ := Complex.re_le_norm _
    _ = ‖∑ p ∈ N.primeFactors, logDeriv (fun w : ℂ => (1 : ℂ) - (p : ℂ) ^ (-w)) c‖ := by
        rw [norm_neg, hFval]
    _ ≤ ∑ p ∈ N.primeFactors, ‖logDeriv (fun w : ℂ => (1 : ℂ) - (p : ℂ) ^ (-w)) c‖ :=
        norm_sum_le _ _
    _ ≤ ∑ p ∈ N.primeFactors, Real.log p := Finset.sum_le_sum hper

#print axioms prod_correction_re_le
