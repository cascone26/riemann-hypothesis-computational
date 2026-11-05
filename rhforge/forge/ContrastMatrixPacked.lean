import ContrastMatrixIsometry

open scoped MatrixOrder Matrix RealInnerProductSpace

/-- **The Helmert-equivalent contrast matrix `C_{N,L}`, packaged as an honest `Matrix` (Shi Eq 4.5,
arXiv:2609.04908)**: the isometric embedding of §2332, made concrete as a `Matrix (Fin n) (Fin
(n-1)) ℝ` via `LinearMap.toMatrix'` after unwrapping `EuclideanSpace` through `EuclideanSpace.equiv`
on both sides. This is a genuine, usable instance of Eq (4.5)'s `C_{N,L}` — not merely an existence
statement. Fully general, no zeta/RH dependence. -/
noncomputable def contrastBasis (n : ℕ) (hn : 1 ≤ n) :
    OrthonormalBasis (Fin (n - 1)) ℝ (meanZeroSubmodule n) :=
  (stdOrthonormalBasis ℝ (meanZeroSubmodule n)).reindex (finCongr (meanZeroSubmodule_finrank n hn))

noncomputable def contrastEmbedding (n : ℕ) (hn : 1 ≤ n) :
    EuclideanSpace ℝ (Fin (n - 1)) →ₗᵢ[ℝ] EuclideanSpace ℝ (Fin n) :=
  (meanZeroSubmodule n).subtypeₗᵢ.comp (contrastBasis n hn).repr.symm.toLinearIsometry

noncomputable def contrastEmbeddingPi (n : ℕ) (hn : 1 ≤ n) :
    (Fin (n - 1) → ℝ) →ₗ[ℝ] (Fin n → ℝ) :=
  (EuclideanSpace.equiv (Fin n) ℝ).toLinearMap.comp
    ((contrastEmbedding n hn).toLinearMap.comp
      (EuclideanSpace.equiv (Fin (n - 1)) ℝ).symm.toLinearMap)

noncomputable def C_matrix (n : ℕ) (hn : 1 ≤ n) : Matrix (Fin n) (Fin (n - 1)) ℝ :=
  LinearMap.toMatrix' (contrastEmbeddingPi n hn)

/-- **`C_matrix`'s `mulVec` is injective** (Shi Eq 4.5's orthonormality condition, in matrix form):
follows since `contrastEmbeddingPi` is a composition of two bijective equivs around a genuine
isometry (`contrastEmbedding`, injective as all linear isometries are), and `LinearMap.toMatrix'`
composed with `mulVec` recovers the original linear map exactly. Fully general, no zeta/RH
dependence. -/
theorem C_matrix_mulVec_injective (n : ℕ) (hn : 1 ≤ n) :
    Function.Injective (C_matrix n hn).mulVec := by
  have hmv : (C_matrix n hn).mulVec = contrastEmbeddingPi n hn := by
    funext v
    exact LinearMap.toMatrix'_mulVec (contrastEmbeddingPi n hn) v
  rw [hmv]
  intro x y hxy
  have heq : contrastEmbeddingPi n hn x = contrastEmbeddingPi n hn y := hxy
  simp only [contrastEmbeddingPi, LinearMap.comp_apply] at heq
  exact (EuclideanSpace.equiv (Fin (n - 1)) ℝ).symm.injective
    ((contrastEmbedding n hn).injective ((EuclideanSpace.equiv (Fin n) ℝ).injective heq))

#print axioms C_matrix_mulVec_injective
