import CCMLoewner

/-!
# S2939 — CCM weight resonance set + Rayleigh/trace monotonicity (lane 2)

The CCM/Shi weight w_γ = 4 sin²(γL/2)/L vanishes EXACTLY on the resonance lattice γL ∈ 2πℤ:
`ccmW_eq_zero_iff` / `ccmW_pos_iff`.  This makes the no-resonance hypothesis of S2933 (`ccm_gram_posDef`)
an explicit arithmetic condition on the ordinates.  From the S2937 Loewner sandwich we also get the scalar
consequences: Rayleigh forms and traces are monotone in the weights, and the CCM ones are dominated by the
unweighted Cauchy Gram scaled by 4/L.
-/

open Matrix Finset

namespace CCMResonance

theorem ccmW_eq_zero_iff {L : ℝ} (hL : 0 < L) (γ : ℝ) :
    CCMWeight.ccmW L γ = 0 ↔ ∃ k : ℤ, γ * L = 2 * k * Real.pi := by
  unfold CCMWeight.ccmW
  rw [div_eq_zero_iff, mul_eq_zero, pow_eq_zero_iff (by norm_num), Real.sin_eq_zero_iff]
  constructor
  · rintro ((h | ⟨k, hk⟩) | h)
    · norm_num at h
    · exact ⟨k, by linarith⟩
    · exact absurd h hL.ne'
  · rintro ⟨k, hk⟩
    exact Or.inl (Or.inr ⟨k, by linarith⟩)

theorem ccmW_pos_iff {L : ℝ} (hL : 0 < L) (γ : ℝ) :
    0 < CCMWeight.ccmW L γ ↔ ∀ k : ℤ, γ * L ≠ 2 * k * Real.pi := by
  rw [lt_iff_le_and_ne, ne_comm, Ne, ccmW_eq_zero_iff hL]
  simp [CCMWeight.ccmW_nonneg hL γ]

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem rayleigh_mono (C : Matrix n m ℝ) {w w' : m → ℝ} (h : ∀ j, w j ≤ w' j) (x : n → ℝ) :
    x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x) ≤ x ⬝ᵥ ((C * diagonal w' * Cᵀ) *ᵥ x) := by
  have hpsd := CCMLoewner.gram_mono C h
  have := hpsd.dotProduct_mulVec_nonneg x
  simp only [star_trivial, Matrix.sub_mulVec, dotProduct_sub] at this
  linarith

theorem trace_mono (C : Matrix n m ℝ) {w w' : m → ℝ} (h : ∀ j, w j ≤ w' j) :
    (C * diagonal w * Cᵀ).trace ≤ (C * diagonal w' * Cᵀ).trace := by
  have := (CCMLoewner.gram_mono C h).trace_nonneg
  rw [Matrix.trace_sub] at this
  linarith

theorem ccm_rayleigh_le {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) (x : n → ℝ) :
    x ⬝ᵥ ((C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) *ᵥ x)
      ≤ 4 / L * (x ⬝ᵥ ((C * Cᵀ) *ᵥ x)) := by
  have h := rayleigh_mono C (w' := fun _ : m => 4 / L) (fun j => CCMWeight.ccmW_le hL (γ j)) x
  rw [CCMLoewner.const_gram, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul] at h
  exact h

theorem ccm_trace_le {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ).trace ≤ 4 / L * (C * Cᵀ).trace := by
  have h := trace_mono C (w' := fun _ : m => 4 / L) (fun j => CCMWeight.ccmW_le hL (γ j))
  rw [CCMLoewner.const_gram, Matrix.trace_smul, smul_eq_mul] at h
  exact h

end CCMResonance

#print axioms CCMResonance.ccmW_eq_zero_iff
#print axioms CCMResonance.ccmW_pos_iff
#print axioms CCMResonance.ccm_rayleigh_le
#print axioms CCMResonance.ccm_trace_le
