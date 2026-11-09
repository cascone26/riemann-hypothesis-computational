import CCMResonanceGap

/-!
# S2958 — a small zero-side margin FORCES a near-resonant ordinate (lane 2)

Contrapositive of S2956: if some test vector x has Rayleigh ratio of the CCM zero side below (4δ²/L) relative to the
Cauchy Gram, i.e. xᵀG_ccm x < (4δ²/L)·xᵀCCᵀx, then some ordinate γ_j satisfies |sin(γ_jL/2)| < δ.
* `exists_near_resonant` — the structural reading of the numerically observed margin decay μ_λ ~ e^{−κλ²}:
  a margin of size μ (relative to CCᵀ) certifies an ordinate within resonance gap √(μL)/2.
* `exists_near_resonant_sqrt` — the same with δ = √(μL)/2 made explicit.
-/

open Matrix Finset

namespace CCMSmall

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem exists_near_resonant {L δ : ℝ} (hL : 0 < L) (hδ : 0 ≤ δ) (C : Matrix n m ℝ) (γ : m → ℝ)
    (x : n → ℝ)
    (hsmall : x ⬝ᵥ ((C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) *ᵥ x)
      < 4 * δ ^ 2 / L * (x ⬝ᵥ ((C * Cᵀ) *ᵥ x))) :
    ∃ j, |Real.sin (γ j * L / 2)| < δ := by
  by_contra h
  push Not at h
  exact absurd (CCMGap.ccm_rayleigh_ge_gap hL hδ C γ h x) (not_le.mpr hsmall)

theorem exists_near_resonant_sqrt {L μ : ℝ} (hL : 0 < L) (hμ : 0 ≤ μ) (C : Matrix n m ℝ)
    (γ : m → ℝ) (x : n → ℝ)
    (hsmall : x ⬝ᵥ ((C * diagonal (fun j => CCMWeight.ccmW L (γ j)) * Cᵀ) *ᵥ x)
      < μ * (x ⬝ᵥ ((C * Cᵀ) *ᵥ x))) :
    ∃ j, |Real.sin (γ j * L / 2)| < Real.sqrt (μ * L) / 2 := by
  have hδ : 0 ≤ Real.sqrt (μ * L) / 2 := by positivity
  refine exists_near_resonant hL hδ C γ x ?_
  have : 4 * (Real.sqrt (μ * L) / 2) ^ 2 / L = μ := by
    rw [div_pow, Real.sq_sqrt (by positivity)]; field_simp; ring
  rw [this]; exact hsmall

end CCMSmall

#print axioms CCMSmall.exists_near_resonant
#print axioms CCMSmall.exists_near_resonant_sqrt
