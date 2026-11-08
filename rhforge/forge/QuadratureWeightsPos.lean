/-
QuadratureWeightsPos.lean (HP RH center, LEDGER S2911, 2026-11-08)
Positivity half of S2868/S2870: if T = C A Cᵀ is positive semidefinite (as a quadratic form) and C is injective, then
A is positive semidefinite; in particular every diagonal entry of A is ≥ 0. With T = τ − λ_min I ⪰ 0 and A diagonal
(QuadratureDiagonal, S2896) this machine-verifies that the quadrature μ_N = Σ_j α_j δ_{t_j} has NONNEGATIVE weights.
-/
import Mathlib

open Matrix

namespace QuadPos

variable {m k : Type*} [Fintype m] [Fintype k] [DecidableEq m] [DecidableEq k]

/-- injective C ⟹ every v is Cᵀ u for some u (via a left inverse). -/
lemma exists_transpose_preimage (C : Matrix m k ℝ) (hinj : Function.Injective C.mulVec) (v : k → ℝ) :
    ∃ u : m → ℝ, Cᵀ *ᵥ u = v := by
  have hker : LinearMap.ker (Matrix.toLin' C) = ⊥ :=
    LinearMap.ker_eq_bot.mpr (by intro a b h; exact hinj (by simpa [Matrix.toLin'_apply] using h))
  obtain ⟨g, hg⟩ := LinearMap.exists_leftInverse_of_injective (Matrix.toLin' C) hker
  set L := LinearMap.toMatrix' g
  have hLC : L * C = 1 := by
    have : LinearMap.toMatrix' (g ∘ₗ Matrix.toLin' C) = LinearMap.toMatrix' LinearMap.id := by rw [hg]
    rwa [LinearMap.toMatrix'_comp, LinearMap.toMatrix'_id, LinearMap.toMatrix'_toLin'] at this
  refine ⟨Lᵀ *ᵥ v, ?_⟩
  rw [Matrix.mulVec_mulVec, ← transpose_mul, hLC, transpose_one, one_mulVec]

theorem psd_of_sandwich (C : Matrix m k ℝ) (A : Matrix k k ℝ) (hinj : Function.Injective C.mulVec)
    (hT : ∀ u : m → ℝ, 0 ≤ u ⬝ᵥ ((C * A * Cᵀ) *ᵥ u)) (v : k → ℝ) : 0 ≤ v ⬝ᵥ (A *ᵥ v) := by
  obtain ⟨u, rfl⟩ := exists_transpose_preimage C hinj v
  have h := hT u
  rw [← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec] at h
  nth_rewrite 1 [mulVec_transpose]
  exact h

theorem diag_nonneg (C : Matrix m k ℝ) (A : Matrix k k ℝ) (hinj : Function.Injective C.mulVec)
    (hT : ∀ u : m → ℝ, 0 ≤ u ⬝ᵥ ((C * A * Cᵀ) *ᵥ u)) (j : k) : 0 ≤ A j j := by
  have := psd_of_sandwich C A hinj hT (Pi.single j 1)
  simpa [mulVec, dotProduct, Pi.single_apply] using this

end QuadPos

#print axioms QuadPos.psd_of_sandwich
#print axioms QuadPos.diag_nonneg
