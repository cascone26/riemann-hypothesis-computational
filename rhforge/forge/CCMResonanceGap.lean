import CCMOscBounds

/-!
# S2956 — quantitative positivity from a resonance gap (lane 2)

If every ordinate stays δ-away from the resonance lattice in the sense |sin(γ_jL/2)| ≥ δ, then
w_j ≥ 4δ²/L and the zero side dominates a FIXED multiple of the Cauchy Gram:
* `ccmW_ge_of_gap` : δ ≤ |sin(γL/2)| ⟹ 4δ²/L ≤ w_γ;
* `ccm_gram_ge_gap` : (4δ²/L)·C·Cᵀ ⪯ G_ccm;
* `ccm_rayleigh_ge_gap` : (4δ²/L)·xᵀCCᵀx ≤ xᵀG_ccm x.
Together with S2937 this is the two-sided sandwich (4δ²/L)·CCᵀ ⪯ G_ccm ⪯ (4/L)·CCᵀ — the zero-side margin is
pinned to the Cauchy Gram up to the factor δ² (the resonance gap).
-/

open Matrix Finset

namespace CCMGap

theorem ccmW_ge_of_gap {L δ γ : ℝ} (hL : 0 < L) (hδ : 0 ≤ δ) (h : δ ≤ |Real.sin (γ * L / 2)|) :
    4 * δ ^ 2 / L ≤ CCMWeight.ccmW L γ := by
  unfold CCMWeight.ccmW
  apply div_le_div_of_nonneg_right _ hL.le
  have h2 : δ ^ 2 ≤ Real.sin (γ * L / 2) ^ 2 := by
    rw [← sq_abs (Real.sin _)]; exact pow_le_pow_left₀ hδ h 2
  linarith

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem ccm_gram_ge_gap {L δ : ℝ} (hL : 0 < L) (hδ : 0 ≤ δ) (C : Matrix n m ℝ) (γ : m → ℝ)
    (hgap : ∀ j, δ ≤ |Real.sin (γ j * L / 2)|) :
    (C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ - (4 * δ ^ 2 / L) • (C * Cᵀ)).PosSemidef := by
  rw [← CCMLoewner.const_gram]
  exact CCMLoewner.gram_mono C (fun j => ccmW_ge_of_gap hL hδ (hgap j))

theorem ccm_rayleigh_ge_gap {L δ : ℝ} (hL : 0 < L) (hδ : 0 ≤ δ) (C : Matrix n m ℝ) (γ : m → ℝ)
    (hgap : ∀ j, δ ≤ |Real.sin (γ j * L / 2)|) (x : n → ℝ) :
    4 * δ ^ 2 / L * (x ⬝ᵥ ((C * Cᵀ) *ᵥ x))
      ≤ x ⬝ᵥ ((C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) *ᵥ x) := by
  have h := CCMResonance.rayleigh_mono C (w := fun _ : m => 4 * δ ^ 2 / L)
    (fun j => ccmW_ge_of_gap hL hδ (hgap j)) x
  rw [CCMLoewner.const_gram, Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul] at h
  exact h

end CCMGap

#print axioms CCMGap.ccmW_ge_of_gap
#print axioms CCMGap.ccm_gram_ge_gap
#print axioms CCMGap.ccm_rayleigh_ge_gap
