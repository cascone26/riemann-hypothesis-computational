import RiemannLean.Certificate
open Matrix
/-- **End-to-end machine-checked Weil-style positivity certificate.**  A concrete rational
SPD matrix `M` (Weil-form-representative) is `PosSemidef`, proven by exhibiting its exact LDL^T
factorization `M = L·diag(d)·Lᵀ` (nonneg rational pivots `d`) and applying `ldlt_posSemidef`.
This is the concrete instance of the organ-2 Arb → Lean certificate pipeline. -/
theorem weil_style_cert_posSemidef :
    ((!![(3 : ℝ)/2, (1 : ℝ)/4, (0 : ℝ); (1 : ℝ)/4, (5 : ℝ)/2, (1 : ℝ)/3; (0 : ℝ), (1 : ℝ)/3, (2 : ℝ)] : Matrix (Fin 3) (Fin 3) ℝ)).PosSemidef := by
  have hd : (0 : Fin 3 → ℝ) ≤ ![(3 : ℝ)/2, (59 : ℝ)/24, (346 : ℝ)/177] := by
    intro i; fin_cases i <;> norm_num
  have hL := ldlt_posSemidef hd (!![(1 : ℝ), (0 : ℝ), (0 : ℝ); (1 : ℝ)/6, (1 : ℝ), (0 : ℝ); (0 : ℝ), (8 : ℝ)/59, (1 : ℝ)] : Matrix (Fin 3) (Fin 3) ℝ)
  have heq : (!![(3 : ℝ)/2, (1 : ℝ)/4, (0 : ℝ); (1 : ℝ)/4, (5 : ℝ)/2, (1 : ℝ)/3; (0 : ℝ), (1 : ℝ)/3, (2 : ℝ)] : Matrix (Fin 3) (Fin 3) ℝ)
      = (!![(1 : ℝ), (0 : ℝ), (0 : ℝ); (1 : ℝ)/6, (1 : ℝ), (0 : ℝ); (0 : ℝ), (8 : ℝ)/59, (1 : ℝ)] : Matrix (Fin 3) (Fin 3) ℝ) * Matrix.diagonal ![(3 : ℝ)/2, (59 : ℝ)/24, (346 : ℝ)/177] * (!![(1 : ℝ), (0 : ℝ), (0 : ℝ); (1 : ℝ)/6, (1 : ℝ), (0 : ℝ); (0 : ℝ), (8 : ℝ)/59, (1 : ℝ)] : Matrix (Fin 3) (Fin 3) ℝ)ᴴ := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      norm_num [Matrix.mul_apply, Matrix.diagonal_apply, Matrix.conjTranspose_apply,
        Fin.sum_univ_three, Matrix.of_apply, Matrix.cons_val_zero, Matrix.cons_val_one,
        Matrix.head_cons, Matrix.cons_val', Matrix.empty_val', Matrix.cons_val_fin_one, star_trivial]
  rw [heq]; exact hL

#print axioms weil_style_cert_posSemidef
