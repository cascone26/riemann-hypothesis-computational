import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Data.Real.StarOrdered
open Matrix
/-- **General LDL^T certificate** (the Arb→Lean lifting): any matrix `L · diag(d) · Lᴴ` with
`d ≥ 0` is `PosSemidef`.  An organ-2 rigorous Weil certificate is an exact-rational LDL^T
factorization `M = L·diag(d)·Lᴴ` with nonneg rational pivots `d`; this theorem turns ANY such
certificate into a machine-checked Lean proof of positive-semidefiniteness. -/
theorem ldlt_posSemidef {n : Type*} [Fintype n] [DecidableEq n]
    {d : n → ℝ} (hd : 0 ≤ d) (L : Matrix n n ℝ) :
    (L * Matrix.diagonal d * Lᴴ).PosSemidef :=
  (Matrix.PosSemidef.diagonal hd).mul_mul_conjTranspose_same L

#print axioms ldlt_posSemidef
