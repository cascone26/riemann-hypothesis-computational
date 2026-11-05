import ContrastMatrixExistence

open scoped MatrixOrder Matrix RealInnerProductSpace

/-- **An explicit isometric embedding realizing `C_{N,L}`'s existence (Shi Eq 4.5,
arXiv:2609.04908)**: for any `n ≥ 1`, there is a linear isometry from `EuclideanSpace ℝ (Fin (n-1))`
into `EuclideanSpace ℝ (Fin n)` whose range is exactly the zero-mean contrast subspace
`meanZeroSubmodule n`. Built via `stdOrthonormalBasis` on the submodule (reindexed to `Fin (n-1)`
using the dimension fact `meanZeroSubmodule_finrank`), composed with the submodule's isometric
inclusion into the ambient space. Matches the paper's own genericity remark that ANY orthonormal
basis of the contrast subspace gives a valid `C_{N,L}` — Helmert's specific closed form is one
choice among many, and existence alone suffices for Theorem 4.1/9.1's machinery. Fully general, no
zeta/RH dependence. -/
theorem exists_isometric_embedding_meanZeroSubmodule (n : ℕ) (hn : 1 ≤ n) :
    ∃ f : EuclideanSpace ℝ (Fin (n - 1)) →ₗᵢ[ℝ] EuclideanSpace ℝ (Fin n),
      LinearMap.range f.toLinearMap = meanZeroSubmodule n := by
  set b : OrthonormalBasis (Fin (Module.finrank ℝ (meanZeroSubmodule n))) ℝ (meanZeroSubmodule n) :=
    stdOrthonormalBasis ℝ (meanZeroSubmodule n) with hbdef
  set b' : OrthonormalBasis (Fin (n - 1)) ℝ (meanZeroSubmodule n) :=
    b.reindex (finCongr (meanZeroSubmodule_finrank n hn)) with hb'def
  set f : EuclideanSpace ℝ (Fin (n - 1)) →ₗᵢ[ℝ] EuclideanSpace ℝ (Fin n) :=
    (meanZeroSubmodule n).subtypeₗᵢ.comp b'.repr.symm.toLinearIsometry with hfdef
  refine ⟨f, ?_⟩
  have hcomp : f.toLinearMap
      = (meanZeroSubmodule n).subtype.comp b'.repr.symm.toLinearIsometry.toLinearMap := by
    rw [hfdef]
    show ((meanZeroSubmodule n).subtypeₗᵢ.comp b'.repr.symm.toLinearIsometry).toLinearMap = _
    rw [show ((meanZeroSubmodule n).subtypeₗᵢ.comp b'.repr.symm.toLinearIsometry).toLinearMap
        = (meanZeroSubmodule n).subtypeₗᵢ.toLinearMap.comp b'.repr.symm.toLinearIsometry.toLinearMap
        from rfl, Submodule.subtypeₗᵢ_toLinearMap]
  rw [hcomp, LinearMap.range_comp]
  have hsurj : LinearMap.range b'.repr.symm.toLinearIsometry.toLinearMap = ⊤ :=
    LinearMap.range_eq_top.mpr b'.repr.symm.surjective
  rw [hsurj, Submodule.map_top, Submodule.range_subtype]

#print axioms exists_isometric_embedding_meanZeroSubmodule
