import QOffKNullIdentityConj
import COffNullIdentityBoth

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: `c_off` null identity at the CONJUGATE point, BOTH copies

Step 4-d of §2391's 4-step plan — mirrors §2402-2403's `c_off_null_identity`/`c_off_null_identity_
reflected` at the conjugate point `z̄_m`, using last cycle's `Q_off_K_null_identity_conj` in place
of `Q_off_K_null_identity`. The nodal-weight bridge (§2401's `c_off_eq_poly_mul_nodalWeight`) is
already fully generic in its evaluation point (it only concerns `c_off_n` vs. eval-at-NODE `ν_n`,
never referencing `z_m` at all), so the direct restatement carries over unchanged; the reflection
argument (§2392-2401's `freqReflect`/`c_off_freqReflect_invariant`/`freqNode_reflect_eq_neg`) is
likewise generic in the "evaluation point" position, so the reflected-copy proof is a verbatim
mirror. Together with §2405's `c_off_row_sum_vanishes`, these are the last two scalar ingredients
needed before mirroring the row-sum-vanishing theorem itself at `z̄_k`. Search-before-build: grepped
`c_off_null_identity_conj`/`c_off_null_identity_reflected_conj` — fresh. Fully general, no zeta/RH
dependence beyond `L≠0`/`η_m≠0`. -/
theorem c_off_null_identity_conj (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) (m : Fin K)
    (hm : (z m).im ≠ 0) :
    ∑ n : Fin (2 * (2 * K) + 1), (c_off K L z n : ℂ) *
      (starRingEnd ℂ (z m) - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ = 0 := by
  rw [← Q_off_K_null_identity_conj K L hL z m hm]
  apply Finset.sum_congr rfl
  intro n _
  rw [c_off_eq_poly_mul_nodalWeight]

theorem c_off_null_identity_reflected_conj (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) (m : Fin K)
    (hm : (z m).im ≠ 0) :
    ∑ n : Fin (2 * (2 * K) + 1), (c_off K L z n : ℂ) *
      (starRingEnd ℂ (z m) + (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ = 0 := by
  rw [← Function.Bijective.sum_comp (freqReflect_bijective K)
    (fun n : Fin (2 * (2 * K) + 1) => (c_off K L z n : ℂ) *
      (starRingEnd ℂ (z m) + (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹)]
  have heq : ∀ n : Fin (2 * (2 * K) + 1),
      (c_off K L z (freqReflect K n) : ℂ) *
        (starRingEnd ℂ (z m)
          + (2 * Real.pi * (((freqReflect K n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹
      = (c_off K L z n : ℂ) *
          (starRingEnd ℂ (z m)
            - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ := by
    intro n
    have hcast : (2 * Real.pi * (((freqReflect K n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)
        = -(2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) := by
      exact_mod_cast freqNode_reflect_eq_neg K L n
    rw [c_off_freqReflect_invariant, hcast]
    ring
  simp_rw [heq]
  exact c_off_null_identity_conj K L hL z m hm

#print axioms c_off_null_identity_conj
#print axioms c_off_null_identity_reflected_conj
