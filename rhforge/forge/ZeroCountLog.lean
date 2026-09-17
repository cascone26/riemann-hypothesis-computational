import ZeroCountFactor

set_option maxHeartbeats 1000000

open Complex Metric

/-!
# §1674 — the zero count in log form: `Σ m ≤ log(M/‖f c‖)/log 3`

The consumer-facing form of §1673 (`3^Σm ≤ M/‖f c‖`): take logs. This is the `O(log M)`
zero-count that multiplies `1/δ = log T/c` in the Landau bound — together the `O(log²T)`.
-/

theorem sum_mult_le_log_of_factor {c : ℂ} {R M : ℝ} {f g P : ℂ → ℂ} {S : Finset ℂ}
    {m : ℂ → ℕ} (hR : 0 < R) (hPdef : P = fun s => ∏ ρ ∈ S, (s - ρ) ^ (m ρ))
    (hS : ∀ ρ ∈ S, ρ ∈ closedBall c (R / 4))
    (hfac : ∀ s ∈ sphere c R, f s = P s * g s) (hfacc : f c = P c * g c)
    (hfsphere : ∀ s ∈ sphere c R, ‖f s‖ ≤ M)
    (hgdiff : DiffContOnCl ℂ g (ball c R)) (hfc0 : f c ≠ 0) :
    ((∑ ρ ∈ S, m ρ : ℕ) : ℝ) ≤ Real.log (M / ‖f c‖) / Real.log 3 := by
  have h3 := three_pow_sum_mult_le_of_factor hR hPdef hS hfac hfacc hfsphere hgdiff hfc0
  have hlog3 : (0 : ℝ) < Real.log 3 := Real.log_pos (by norm_num)
  have hpow : (0 : ℝ) < (3 : ℝ) ^ (∑ ρ ∈ S, m ρ) := by positivity
  have hmono := Real.log_le_log hpow h3
  rw [Real.log_pow] at hmono
  rw [le_div_iff₀ hlog3]
  exact hmono

#print axioms sum_mult_le_log_of_factor
