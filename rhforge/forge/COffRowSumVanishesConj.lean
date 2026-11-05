import COffRowSumVanishes
import COffNullIdentityConj

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: THE ROW-SUM VANISHES at the CONJUGATE point too (conjugate branch closes)

Step 4-e of §2391's 4-step plan — mirrors §2405's `c_off_row_sum_vanishes` exactly at the conjugate
point `z̄_{m0}` in place of `z_{m0}`, using last cycle's `c_off_null_identity_conj`/`c_off_null_
identity_reflected_conj` in place of the plain versions. `H_L_G_C_entry_partial_fraction`
(§2404) is already fully generic in its evaluation point `z`, so the per-term algebra is a verbatim
copy; the only new input needed is `(starRingEnd ℂ (z m0)).im ≠ 0`, immediate from `(z m0).im ≠ 0`
via `Complex.conj_im`/`neg_ne_zero`. With both this and §2405 in hand, BOTH terms of Theorem 7.2's
Eq (7.2) paired sum `S^{(T,off)}=2∑_k(H^G_L(z_k)+H^G_L(z̄_k))` (`STOffPairedSum.lean`) now vanish
against `e_off` at every row — the scalar content of the FULL matrix null-vector claim is complete;
only the (purely mechanical) k-sum/real-part bookkeeping remains, honestly not attempted here.
Search-before-build: grepped `c_off_row_sum_vanishes_conj` — fresh. Fully general, no zeta/RH
dependence beyond `L≠0`. -/
theorem c_off_row_sum_vanishes_conj (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) (m0 : Fin K)
    (hm0 : (z m0).im ≠ 0) (m : Fin (2 * (2 * K) + 1)) :
    ∑ n : Fin (2 * (2 * K) + 1),
      (4 * (Complex.sin (L * (starRingEnd ℂ (z m0)) / 2)) ^ 2
          * ((starRingEnd ℂ (z m0)) ^ 2
                + (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)
                * (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))
          / (L ^ 2 * ((starRingEnd ℂ (z m0)) ^ 2
                - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2)
              * ((starRingEnd ℂ (z m0)) ^ 2
                - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2)))
        * (e_off K L z n : ℂ) = 0 := by
  have hm0conj : (starRingEnd ℂ (z m0)).im ≠ 0 := by
    simp only [Complex.conj_im]
    exact neg_ne_zero.mpr hm0
  have hcast : ∀ (n : Fin (2 * (2 * K) + 1)) (w : ℂ),
      (e_off K L z n : ℂ) * w = (Real.sqrt L : ℂ) * ((c_off K L z n : ℂ) * w) := by
    intro n w
    unfold e_off
    push_cast
    ring
  have hsum1 : ∑ n : Fin (2 * (2 * K) + 1),
      (e_off K L z n : ℂ)
        * ((starRingEnd ℂ (z m0))
            - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ = 0 := by
    simp_rw [hcast]
    rw [← Finset.mul_sum, c_off_null_identity_conj K L hL z m0 hm0, mul_zero]
  have hsum2 : ∑ n : Fin (2 * (2 * K) + 1),
      (e_off K L z n : ℂ)
        * ((starRingEnd ℂ (z m0))
            + (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ = 0 := by
    simp_rw [hcast]
    rw [← Finset.mul_sum, c_off_null_identity_reflected_conj K L hL z m0 hm0, mul_zero]
  have hterm : ∀ n : Fin (2 * (2 * K) + 1),
      (4 * (Complex.sin (L * (starRingEnd ℂ (z m0)) / 2)) ^ 2
          * ((starRingEnd ℂ (z m0)) ^ 2
                + (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)
                * (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))
          / (L ^ 2 * ((starRingEnd ℂ (z m0)) ^ 2
                - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2)
              * ((starRingEnd ℂ (z m0)) ^ 2
                - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2)))
        * (e_off K L z n : ℂ)
      = (4 * (Complex.sin (L * (starRingEnd ℂ (z m0)) / 2)) ^ 2
            / (L ^ 2 * ((starRingEnd ℂ (z m0)) ^ 2
                  - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2))
          * (((starRingEnd ℂ (z m0))
                + (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)) / 2))
          * ((e_off K L z n : ℂ)
              * ((starRingEnd ℂ (z m0))
                  - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹)
        + (4 * (Complex.sin (L * (starRingEnd ℂ (z m0)) / 2)) ^ 2
            / (L ^ 2 * ((starRingEnd ℂ (z m0)) ^ 2
                  - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) ^ 2))
          * (((starRingEnd ℂ (z m0))
                - (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)) / 2))
          * ((e_off K L z n : ℂ)
              * ((starRingEnd ℂ (z m0))
                  + (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹) := by
    intro n
    have hzn : (starRingEnd ℂ (z m0))
        ≠ (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) := by
      have := ne_of_im_ne_zero hm0conj (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
      simpa using this
    have hzn' : (starRingEnd ℂ (z m0))
        ≠ -(2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) := by
      have := ne_of_im_ne_zero hm0conj (-(2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L))
      simpa using this
    rw [H_L_G_C_entry_partial_fraction L (starRingEnd ℂ (z m0))
      (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
      (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L) hzn hzn']
    simp only [div_eq_mul_inv, mul_inv]
    ring
  rw [Finset.sum_congr rfl (fun n _ => hterm n), Finset.sum_add_distrib, ← Finset.mul_sum,
    ← Finset.mul_sum, hsum1, hsum2]
  ring

#print axioms c_off_row_sum_vanishes_conj
