import Mathlib

open Matrix
open scoped BigOperators

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- Inverse of a diagonal matrix with nonvanishing entries. -/
theorem inv_diagonal_of_ne_zero (g : n → ℝ) (hg : ∀ j, g j ≠ 0) :
    (diagonal g)⁻¹ = diagonal (fun j => (g j)⁻¹) := by
  apply Matrix.inv_eq_right_inv
  rw [Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
  congr 1
  funext j
  exact mul_inv_cancel₀ (hg j)

/-- **Diagonal resolvent quadratic form.**

For a diagonal operator `diagonal g` on `ℝⁿ` with nonvanishing entries,
    `⟨η, (diagonal g)⁻¹ ξ⟩ = Σ_j η_j · (g_j)⁻¹ · ξ_j`.
Specializing `g_j = j − s` (with `s` off the spectrum) and `η_j ≡ 1` turns the CM spectral-realization
determinant (`detDprime_core`, cy380) `⟨η,(D−s)⁻¹ξ⟩` into `Σ_j (j−s)⁻¹ ξ_j = ξ̂(s)` — so the determinant's
zeros are the zeros of the finite Fourier transform of `ξ`.  Machine-verified, axiom-clean. -/
theorem dotProduct_inv_diagonal_mulVec (g ξ η : n → ℝ) (hg : ∀ j, g j ≠ 0) :
    η ⬝ᵥ ((diagonal g)⁻¹ *ᵥ ξ) = ∑ j, η j * ((g j)⁻¹ * ξ j) := by
  have hinv : (diagonal g)⁻¹ *ᵥ ξ = fun j => (g j)⁻¹ * ξ j := by
    rw [inv_diagonal_of_ne_zero g hg]
    funext j
    rw [Matrix.mulVec_diagonal]
  rw [hinv]
  simp only [dotProduct]

#print axioms dotProduct_inv_diagonal_mulVec

end CMKeyLemma
