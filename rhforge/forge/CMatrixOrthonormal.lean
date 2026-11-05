import ContrastMatrixPacked

open scoped MatrixOrder Matrix

lemma euclidean_equiv_dotProduct (k : ℕ) (v w : EuclideanSpace ℝ (Fin k)) :
    dotProduct (EuclideanSpace.equiv (Fin k) ℝ v) (EuclideanSpace.equiv (Fin k) ℝ w)
      = (inner ℝ v w : ℝ) := by
  rw [EuclideanSpace.inner_eq_star_dotProduct]
  unfold EuclideanSpace.equiv
  simp [PiLp.coe_continuousLinearEquiv, dotProduct_comm]

theorem C_matrix_dotProduct_mulVec_eq (n : ℕ) (hn : 1 ≤ n) (x y : Fin (n - 1) → ℝ) :
    dotProduct ((C_matrix n hn).mulVec x) ((C_matrix n hn).mulVec y)
      = dotProduct x y := by
  have hmv : (C_matrix n hn).mulVec = contrastEmbeddingPi n hn := by
    funext v
    exact LinearMap.toMatrix'_mulVec (contrastEmbeddingPi n hn) v
  rw [hmv]
  unfold contrastEmbeddingPi
  simp only [LinearMap.comp_apply, LinearEquiv.coe_toLinearMap, LinearIsometry.coe_toLinearMap]
  set X := (EuclideanSpace.equiv (Fin (n - 1)) ℝ).symm x with hX
  set Y := (EuclideanSpace.equiv (Fin (n - 1)) ℝ).symm y with hY
  have step1 := euclidean_equiv_dotProduct n (contrastEmbedding n hn X) (contrastEmbedding n hn Y)
  have step2 := (contrastEmbedding n hn).inner_map_map X Y
  have step3 := euclidean_equiv_dotProduct (n - 1) X Y
  rw [hX, hY, ContinuousLinearEquiv.apply_symm_apply, ContinuousLinearEquiv.apply_symm_apply] at step3
  calc (EuclideanSpace.equiv (Fin n) ℝ) (contrastEmbedding n hn X)
        ⬝ᵥ (EuclideanSpace.equiv (Fin n) ℝ) (contrastEmbedding n hn Y)
      = (inner ℝ (contrastEmbedding n hn X) (contrastEmbedding n hn Y) : ℝ) := step1
    _ = (inner ℝ X Y : ℝ) := step2
    _ = dotProduct x y := step3.symm

theorem C_matrix_transpose_mul_self (n : ℕ) (hn : 1 ≤ n) :
    (C_matrix n hn)ᵀ * (C_matrix n hn) = 1 := by
  have key : ∀ y : Fin (n - 1) → ℝ, ((C_matrix n hn)ᵀ * (C_matrix n hn)).mulVec y = y := by
    intro y
    funext i
    have hgen : ∀ x : Fin (n - 1) → ℝ,
        dotProduct x (((C_matrix n hn)ᵀ * (C_matrix n hn)).mulVec y) = dotProduct x y := by
      intro x
      rw [← Matrix.mulVec_mulVec, Matrix.dotProduct_mulVec, Matrix.vecMul_transpose,
          C_matrix_dotProduct_mulVec_eq]
    have hx := hgen (Pi.single i 1)
    rwa [single_dotProduct, single_dotProduct, one_mul, one_mul] at hx
  apply Matrix.ext
  intro i j
  have hk := key (Pi.single j 1)
  rw [Matrix.mulVec_single_one] at hk
  have hcol : ((C_matrix n hn)ᵀ * (C_matrix n hn)).col j i
      = (Pi.single j 1 : Fin (n - 1) → ℝ) i := by
    rw [hk]
  simpa [Matrix.col, Pi.single_apply, Matrix.one_apply, eq_comm] using hcol

#print axioms C_matrix_dotProduct_mulVec_eq
#print axioms C_matrix_transpose_mul_self
