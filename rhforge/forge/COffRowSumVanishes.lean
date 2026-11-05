import HLGCEntryPartialFraction
import COffNullIdentityBoth
import STOffSumEntryClosedForm

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: THE ROW-SUM VANISHES at a single off-critical-line point (step 4 CLOSES)

Step 4-b of §2391's 4-step plan — the final assembly. For a fixed off-critical-line parameter
`z m0` (nonzero imaginary part), assembles: §2404's `H_L_G_C_entry_partial_fraction` (splits each
row entry into prefactor × two-term partial fraction), `Finset.sum_add_distrib`/`Finset.mul_sum`
(distributes the sum and factors the n-independent prefactors out), and §2402-2403's
`c_off_null_identity`/`c_off_null_identity_reflected` (both inner sums vanish) into the single
closing statement: the row-sum `∑_n H^G_{mn,L}(z_{m0})·e^{off}_n` is IDENTICALLY ZERO, for every row
`m` and every off-critical-line parameter with nonzero imaginary part. This is the scalar-level
content of Theorem 7.2's matrix claim `S^{(T,off)}·e^{off}=0` at a single `z_k`/`z̄_k` term — summing
this over `k` (via `S_T_off_sum_entry_closed_form`) and taking real parts (via `S_T_off_real`'s
definition) is the only remaining bookkeeping to reach the literal matrix statement, honestly not
attempted here (real remaining work, not a wall). Search-before-build: grepped
`c_off_row_sum_vanishes`/`row_sum_vanishes` — fresh. Fully general, no zeta/RH dependence beyond
`L≠0`. -/
theorem c_off_row_sum_vanishes (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) (m0 : Fin K)
    (hm0 : (z m0).im ≠ 0) (m : Fin (2 * (2 * K) + 1)) :
    ∑ n : Fin (2 * (2 * K) + 1),
      (4 * (Complex.sin (L * (z m0) / 2)) ^ 2
          * ((z m0) ^ 2 + (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)
                * (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))
          / (L ^ 2 * ((z m0) ^ 2 - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2)
              * ((z m0) ^ 2 - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2)))
        * (e_off K L z n : ℂ) = 0 := by
  have hcast : ∀ (n : Fin (2 * (2 * K) + 1)) (w : ℂ),
      (e_off K L z n : ℂ) * w = (Real.sqrt L : ℂ) * ((c_off K L z n : ℂ) * w) := by
    intro n w
    unfold e_off
    push_cast
    ring
  have hsum1 : ∑ n : Fin (2 * (2 * K) + 1),
      (e_off K L z n : ℂ)
        * ((z m0) - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ = 0 := by
    simp_rw [hcast]
    rw [← Finset.mul_sum, c_off_null_identity K L hL z m0 hm0, mul_zero]
  have hsum2 : ∑ n : Fin (2 * (2 * K) + 1),
      (e_off K L z n : ℂ)
        * ((z m0) + (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ = 0 := by
    simp_rw [hcast]
    rw [← Finset.mul_sum, c_off_null_identity_reflected K L hL z m0 hm0, mul_zero]
  have hterm : ∀ n : Fin (2 * (2 * K) + 1),
      (4 * (Complex.sin (L * (z m0) / 2)) ^ 2
          * ((z m0) ^ 2 + (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)
                * (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))
          / (L ^ 2 * ((z m0) ^ 2 - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2)
              * ((z m0) ^ 2 - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2)))
        * (e_off K L z n : ℂ)
      = (4 * (Complex.sin (L * (z m0) / 2)) ^ 2
            / (L ^ 2 * ((z m0) ^ 2 - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2))
          * (((z m0) + (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)) / 2))
          * ((e_off K L z n : ℂ)
              * ((z m0) - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹)
        + (4 * (Complex.sin (L * (z m0) / 2)) ^ 2
            / (L ^ 2 * ((z m0) ^ 2 - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2))
          * (((z m0) - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)) / 2))
          * ((e_off K L z n : ℂ)
              * ((z m0) + (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹) := by
    intro n
    have hzn : (z m0) ≠ (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) := by
      have := ne_of_im_ne_zero hm0 (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
      simpa using this
    have hzn' : (z m0) ≠ -(2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) := by
      have := ne_of_im_ne_zero hm0 (-(2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L))
      simpa using this
    rw [H_L_G_C_entry_partial_fraction L (z m0)
      (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
      (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) hzn hzn']
    simp only [div_eq_mul_inv, mul_inv]
    ring
  rw [Finset.sum_congr rfl (fun n _ => hterm n), Finset.sum_add_distrib, ← Finset.mul_sum,
    ← Finset.mul_sum, hsum1, hsum2]
  ring

#print axioms c_off_row_sum_vanishes
