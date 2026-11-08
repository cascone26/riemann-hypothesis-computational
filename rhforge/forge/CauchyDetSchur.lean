/-
CauchyDetSchur.lean (HP RH center, LEDGER S2912, 2026-11-08) — Cauchy determinant formula, step 1 (Schur identity).
Absent from Mathlib (forge had only nonvanishing). For a, b with a + b' ≠ 0 etc., the Schur complement of the last
row/column of a Cauchy matrix is a diagonally-scaled Cauchy matrix:
  1/(x_i+y_j) − (1/(x_i+y_n))·(x_n+y_n)·(1/(x_n+y_j))
      = [(x_n−x_i)/(x_i+y_n)] · 1/(x_i+y_j) · [(y_n−y_j)/(x_n+y_j)].
Next: det_fromBlocks₂₂ + reindex ⟹ det C_{n+1} = (x_n+y_n)⁻¹ · Π u_i · Π v_j · det C_n ⟹ closed form by induction.
-/
import Mathlib

namespace CauchyDet

variable {K : Type*} [Field K]

theorem schur_entry (xi yj xn yn : K) (h1 : xi + yj ≠ 0) (h2 : xi + yn ≠ 0) (h3 : xn + yj ≠ 0)
    (h4 : xn + yn ≠ 0) :
    (xi + yj)⁻¹ - (xi + yn)⁻¹ * (xn + yn) * (xn + yj)⁻¹
      = ((xn - xi) / (xi + yn)) * (xi + yj)⁻¹ * ((yn - yj) / (xn + yj)) := by
  field_simp
  ring

/-- matrix form: the Schur complement is diag(u) · C · diag(v). -/
theorem schur_matrix {n : ℕ} (x y : Fin n → K) (xn yn : K)
    (h1 : ∀ i j, x i + y j ≠ 0) (h2 : ∀ i, x i + yn ≠ 0) (h3 : ∀ j, xn + y j ≠ 0) (h4 : xn + yn ≠ 0) :
    (Matrix.of fun i j => (x i + y j)⁻¹ - (x i + yn)⁻¹ * (xn + yn) * (xn + y j)⁻¹)
      = Matrix.diagonal (fun i => (xn - x i) / (x i + yn)) * Matrix.of (fun i j => (x i + y j)⁻¹)
          * Matrix.diagonal (fun j => (yn - y j) / (xn + y j)) := by
  ext i j
  simp only [Matrix.of_apply, Matrix.mul_diagonal, Matrix.diagonal_mul]
  rw [schur_entry _ _ _ _ (h1 i j) (h2 i) (h3 j) h4]

end CauchyDet

#print axioms CauchyDet.schur_entry
#print axioms CauchyDet.schur_matrix
