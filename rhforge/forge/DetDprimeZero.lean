import DprimeKernel

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **`D'` is singular — `det D' = 0` (the core of CM's `det0`).**

Since `D' = D − |Dξ⟩⟨η|` annihilates the nonzero vector `ξ` (`Dprime_mulVec_xi`, needs `⟨η,ξ⟩=1`),
`D'` has `0` as an eigenvalue, hence `det D' = 0`.  This is the reachable, precision-free core of CM's
`det0` factorization `Det(D'−s) = −s·Det(D''−s)`: it is exactly the statement that `0` is one of `D'`'s
eigenvalues, which supplies the `(−s)` factor.  (The full `det0` additionally identifies the *complementary*
determinant with the quotient operator `D''` on `E_N/ℂξ` — that needs an explicit-basis block-triangular
change of basis / quotient `LinearMap.det`, deferred as Mathlib-infra-heavy.)  Machine-verified, axiom-clean. -/
theorem det_Dprime_eq_zero (D : Matrix n n ℝ) (η ξ : n → ℝ)
    (hη : η ⬝ᵥ ξ = 1) (hξ0 : ξ ≠ 0) :
    (D - Matrix.vecMulVec (D *ᵥ ξ) η).det = 0 := by
  rw [← Matrix.exists_mulVec_eq_zero_iff]
  exact ⟨ξ, hξ0, Dprime_mulVec_xi D η ξ hη⟩

#print axioms det_Dprime_eq_zero

end CMKeyLemma
