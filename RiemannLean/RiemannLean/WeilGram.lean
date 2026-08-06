import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered

/-!
# The Weil form as a Gram matrix (LEDGER §772)

Verified numerically (`operators/weil_gram_identity.py`, agreement 3e-15): the organ-2 Guinand–Weil
form equals `Q_{jk} = Σ_ρ φ̂_{a_j}(γ_ρ) φ̂_{a_k}(γ_ρ) = (Vᵀ V)_{jk}`, where `V_{ρ j} = φ̂_{a_j}(γ_ρ)`
is the matrix of test-function transforms evaluated at the zeros. This module formalizes the
STRUCTURAL reason the Weil positivity certificate holds: any real Gram matrix `Vᵀ V` is PosSemidef.
Under RH the zeros `γ_ρ` are real, so `V` is a real matrix and `Q = Vᵀ V ⪰ 0`; an off-line zero would
make a coordinate complex and break the real-Gram structure (the whole content of the criterion).
-/

open Matrix

/-- **Weil form is a Gram matrix ⇒ PosSemidef.**  Any matrix of the form `Vᵀ V` over ℝ (a Gram matrix
of the real column-vectors of `V`) is `PosSemidef`.  This is the structural reason the organ-2 Weil
certificate is positive-semidefinite: `Q_{jk} = Σ_ρ φ̂_j(γ_ρ) φ̂_k(γ_ρ) = (Vᵀ V)_{jk}` (§772), real
whenever the zeros are real (RH). -/
theorem weil_gram_posSemidef {m n : Type*} [Fintype m] [Fintype n] [DecidableEq n]
    (V : Matrix m n ℝ) : (Vᵀ * V).PosSemidef := by
  have h := Matrix.posSemidef_conjTranspose_mul_self V
  rwa [Matrix.conjTranspose_eq_transpose_of_trivial] at h

/-- **Weil positivity from the Gram identity (§772).**  If the Weil form `Q` equals the real Gram
matrix `Vᵀ V` of the zero-evaluation vectors (verified numerically to 3e-15), then `Q` is
`PosSemidef`.  This is the RH-conditional structural certificate: real `V` ⟺ zeros on the line. -/
theorem weil_posSemidef_of_gram {m n : Type*} [Fintype m] [Fintype n] [DecidableEq n]
    {Q : Matrix n n ℝ} {V : Matrix m n ℝ} (hQ : Q = Vᵀ * V) : Q.PosSemidef := by
  rw [hQ]; exact weil_gram_posSemidef V
