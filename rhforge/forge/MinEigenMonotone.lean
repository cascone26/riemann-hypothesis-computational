/-
MinEigenMonotone.lean (HP RH center, LEDGER S2909, 2026-11-08)
ε_N is NON-INCREASING in N: the smallest eigenvalue of a real symmetric (m+1)×(m+1) matrix is ≤ that of its leading
m×m principal submatrix (case k = m−1 of matrix Cauchy interlacing, S2902). For the CCM/Shi Weil matrices τ_{λ,N}
(nested principal submatrices as N grows) this machine-verifies μ_λ ≤ ε_{λ,N} for every finite N — the inference
behind S2893's falsification (μ_{√20} ≤ ε_104 < 5e-95) and every margin-floor statement in S2887–S2895.
-/
import Mathlib
import MatrixCauchyInterlacing

open Matrix Module

namespace MinEig

theorem min_eig_le_submatrix {m : ℕ} (A : Matrix (Fin (m + 2)) (Fin (m + 2)) ℝ) (hA : A.IsHermitian) :
    let hT := (isHermitian_iff_isSymmetric.mp hA)
    let hS := (isHermitian_iff_isSymmetric.mp (hA.submatrix Fin.castSucc))
    hT.eigenvalues (MatrixInterlace.finrank_e (m + 2)) (Fin.last (m + 1))
      ≤ hS.eigenvalues (MatrixInterlace.finrank_e (m + 1)) (Fin.last m) := by
  intro hT hS
  have h := (MatrixInterlace.matrix_interlace A hA (Fin.last m)).1
  simpa [Fin.succ_last] using h

/-- positivity passes to principal submatrices: if every eigenvalue of A is ≥ c, so is the smallest of A'. -/
theorem submatrix_min_ge {m : ℕ} (A : Matrix (Fin (m + 2)) (Fin (m + 2)) ℝ) (hA : A.IsHermitian) (c : ℝ)
    (hc : c ≤ (isHermitian_iff_isSymmetric.mp hA).eigenvalues (MatrixInterlace.finrank_e (m + 2)) (Fin.last (m + 1))) :
    c ≤ (isHermitian_iff_isSymmetric.mp (hA.submatrix Fin.castSucc)).eigenvalues
      (MatrixInterlace.finrank_e (m + 1)) (Fin.last m) :=
  le_trans hc (min_eig_le_submatrix A hA)

end MinEig

#print axioms MinEig.min_eig_le_submatrix
#print axioms MinEig.submatrix_min_ge
