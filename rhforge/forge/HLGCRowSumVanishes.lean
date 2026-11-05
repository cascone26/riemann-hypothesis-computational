import COffRowSumVanishesConj

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: `H_L_G_C`'s own row-sum vanishes (both `z_k` and `z̄_k`)

Step 4-f of §2391's 4-step plan: wraps §2405's `c_off_row_sum_vanishes` and §2410's
`c_off_row_sum_vanishes_conj` (both stated via the EXPLICIT closed-form expression) back onto the
actual matrix object `H_L_G_C` itself, via §2386's `H_L_G_C_entry_closed_form` bridge (supplying its
four pole-avoidance hypotheses uniformly via `ne_of_im_ne_zero`, exactly as
`S_T_off_sum_entry_closed_form`'s own proof does). This is the direct prerequisite for summing over
`k` against `S_T_off_sum`'s own `∑_k(H_L_G_C(z_k)+H_L_G_C(z̄_k))` definition — the last remaining
step toward the literal matrix statement. Search-before-build: grepped
`H_L_G_C_row_sum_vanishes` — fresh. Fully general, no zeta/RH dependence beyond `L≠0`. -/
theorem H_L_G_C_row_sum_vanishes (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) (m0 : Fin K)
    (hm0 : (z m0).im ≠ 0) (m : Fin (2 * (2 * K) + 1)) :
    ∑ n : Fin (2 * (2 * K) + 1), H_L_G_C (2 * K) L (z m0) m n * (e_off K L z n : ℂ) = 0 := by
  rw [← c_off_row_sum_vanishes K L hL z m0 hm0 m]
  apply Finset.sum_congr rfl
  intro n _
  congr 1
  refine H_L_G_C_entry_closed_form (2 * K) L hL (z m0) m n ?_ ?_ ?_ ?_
  · have := ne_of_im_ne_zero hm0 (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
    simpa using this
  · have := ne_of_im_ne_zero hm0 (-(2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L))
    simpa using this
  · have := ne_of_im_ne_zero hm0 (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
    simpa using this
  · have := ne_of_im_ne_zero hm0 (-(2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L))
    simpa using this

theorem H_L_G_C_row_sum_vanishes_conj (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) (m0 : Fin K)
    (hm0 : (z m0).im ≠ 0) (m : Fin (2 * (2 * K) + 1)) :
    ∑ n : Fin (2 * (2 * K) + 1),
      H_L_G_C (2 * K) L (starRingEnd ℂ (z m0)) m n * (e_off K L z n : ℂ) = 0 := by
  have hm0conj : (starRingEnd ℂ (z m0)).im ≠ 0 := by
    simp only [Complex.conj_im]
    exact neg_ne_zero.mpr hm0
  rw [← c_off_row_sum_vanishes_conj K L hL z m0 hm0 m]
  apply Finset.sum_congr rfl
  intro n _
  congr 1
  refine H_L_G_C_entry_closed_form (2 * K) L hL (starRingEnd ℂ (z m0)) m n ?_ ?_ ?_ ?_
  · have := ne_of_im_ne_zero hm0conj (2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
    simpa using this
  · have := ne_of_im_ne_zero hm0conj (-(2 * Real.pi * (((m : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L))
    simpa using this
  · have := ne_of_im_ne_zero hm0conj (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L)
    simpa using this
  · have := ne_of_im_ne_zero hm0conj (-(2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L))
    simpa using this

#print axioms H_L_G_C_row_sum_vanishes
#print axioms H_L_G_C_row_sum_vanishes_conj
