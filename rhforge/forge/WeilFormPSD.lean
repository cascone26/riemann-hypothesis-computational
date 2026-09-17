import Mathlib

/-!
# Weil form positivity from on-line zeros (RH session 2026-08-29)

The classical Weil criterion's easy direction, and the object the Euler-product discriminator studies:
the compressed Weil form `∑_ρ v_ρ v_ρᵀ` — a sum of outer products of the REAL vectors `v_ρ = (F_j(γ_ρ))_j`
attached to on-line zeros `ρ = ½ + iγ_ρ` — is positive semidefinite (a real Gram matrix).
When zeros are on the line the `v_ρ` are real, so RH ⟹ Weil positivity.
-/

open Matrix

/-- A sum of outer products of real vectors is positive semidefinite (the on-line Weil form). -/
theorem weil_form_posSemidef {ι n : Type*} [Fintype ι] [Fintype n] [DecidableEq n]
    (v : ι → n → ℝ) :
    (∑ i, Matrix.vecMulVec (v i) (v i)).PosSemidef := by
  apply Matrix.posSemidef_sum Finset.univ
  intro i _
  simpa using Matrix.posSemidef_vecMulVec_self_star (v i)
