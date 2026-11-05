import COffNodalWeightBridge
import QOffKNullIdentity

open scoped MatrixOrder Matrix

/-!
# Off-critical-line lane: `c_off` null identity, BOTH copies — step 3 CLOSES

Step 3-f/g of §2391's 4-step plan — the actual closing pieces. First restates §2380-2383's
`Q_off_K_null_identity` (built via Theorem 5.2's abstract `barycentric_null_identity`, phrased over
`nodalWeight`) directly in terms of `c_off` using last cycle's
`c_off_eq_poly_mul_nodalWeight` bridge — this is exactly the "−ν_n" copy §2391's
`quadratic_partial_fraction_split` produces. Then, via the `freqReflect` bijection +
`c_off_freqReflect_invariant` + `freqNode_reflect_eq_neg` (all built this session), reindexes to
show the "+ν_n" copy ALSO vanishes — this was the genuinely open half of §2391's split, flagged as
"correctly NOT attempted" at the time. Together, both copies of the row-sum from §2391's partial
fraction split vanish, completing step 3 of the plan. Search-before-build: grepped
`c_off_null_identity` — fresh. Fully general, no zeta/RH dependence beyond `L≠0`/`η_m≠0`. -/
theorem c_off_null_identity (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) (m : Fin K)
    (hm : (z m).im ≠ 0) :
    ∑ n : Fin (2 * (2 * K) + 1), (c_off K L z n : ℂ) *
      (z m - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ = 0 := by
  rw [← Q_off_K_null_identity K L hL z m hm]
  apply Finset.sum_congr rfl
  intro n _
  rw [c_off_eq_poly_mul_nodalWeight]

theorem c_off_null_identity_reflected (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) (m : Fin K)
    (hm : (z m).im ≠ 0) :
    ∑ n : Fin (2 * (2 * K) + 1), (c_off K L z n : ℂ) *
      (z m + (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ = 0 := by
  rw [← Function.Bijective.sum_comp (freqReflect_bijective K)
    (fun n : Fin (2 * (2 * K) + 1) => (c_off K L z n : ℂ) *
      (z m + (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹)]
  have heq : ∀ n : Fin (2 * (2 * K) + 1),
      (c_off K L z (freqReflect K n) : ℂ) *
        (z m + (2 * Real.pi * (((freqReflect K n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹
      = (c_off K L z n : ℂ) *
          (z m - (2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ))⁻¹ := by
    intro n
    have hcast : (2 * Real.pi * (((freqReflect K n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ)
        = -(2 * Real.pi * (((n : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) / L : ℂ) := by
      exact_mod_cast freqNode_reflect_eq_neg K L n
    rw [c_off_freqReflect_invariant, hcast]
    ring
  simp_rw [heq]
  exact c_off_null_identity K L hL z m hm

#print axioms c_off_null_identity
#print axioms c_off_null_identity_reflected
