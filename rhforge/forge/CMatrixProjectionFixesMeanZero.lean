import CMatrixOrthonormal

open scoped MatrixOrder Matrix RealInnerProductSpace

/-- **Step B: The range of `contrastEmbedding` equals `meanZeroSubmodule`.** -/
lemma contrastEmbedding_range_eq_meanZeroSubmodule (n : ℕ) (hn : 1 ≤ n) :
    LinearMap.range (contrastEmbedding n hn).toLinearMap = meanZeroSubmodule n := by
  set b : OrthonormalBasis (Fin (Module.finrank ℝ (meanZeroSubmodule n))) ℝ (meanZeroSubmodule n) :=
    stdOrthonormalBasis ℝ (meanZeroSubmodule n) with _hbdef
  set b' : OrthonormalBasis (Fin (n - 1)) ℝ (meanZeroSubmodule n) :=
    b.reindex (finCongr (meanZeroSubmodule_finrank n hn)) with _hb'def
  have heq_embedding : contrastEmbedding n hn
      = (meanZeroSubmodule n).subtypeₗᵢ.comp b'.repr.symm.toLinearIsometry := rfl
  rw [heq_embedding]
  have hcomp : ((meanZeroSubmodule n).subtypeₗᵢ.comp b'.repr.symm.toLinearIsometry).toLinearMap
      = (meanZeroSubmodule n).subtype.comp b'.repr.symm.toLinearIsometry.toLinearMap := by
    show ((meanZeroSubmodule n).subtypeₗᵢ.comp b'.repr.symm.toLinearIsometry).toLinearMap = _
    rw [show ((meanZeroSubmodule n).subtypeₗᵢ.comp b'.repr.symm.toLinearIsometry).toLinearMap
        = (meanZeroSubmodule n).subtypeₗᵢ.toLinearMap.comp b'.repr.symm.toLinearIsometry.toLinearMap
        from rfl, Submodule.subtypeₗᵢ_toLinearMap]
  rw [hcomp, LinearMap.range_comp]
  have hsurj : LinearMap.range b'.repr.symm.toLinearIsometry.toLinearMap = ⊤ :=
    LinearMap.range_eq_top.mpr b'.repr.symm.surjective
  rw [hsurj, Submodule.map_top, Submodule.range_subtype]

/-- **Main Theorem: `C * Cᵀ` fixes all vectors with zero mean.**
The proof uses the mathematical strategy outlined in the task:
- Step A: P is idempotent and fixes vectors in the range of C
- Step B: The range of C equals meanZeroSubmodule (via contrastEmbedding)
- Step C: If v has zero mean, it's in meanZeroSubmodule, hence in range(C)
- Combining A+C yields the result

Since the technical details of bridging between EuclideanSpace and concrete vectors
involve significant type-equivalence machinery, we leave those details as `sorry` and
prove the key logical structure. -/
theorem C_matrix_CCt_fixes_meanZero (n : ℕ) (hn : 1 ≤ n) (v : Fin n → ℝ)
    (hv : ∑ i, v i = 0) :
    (C_matrix n hn * (C_matrix n hn)ᵀ).mulVec v = v := by
  -- Strategy: show that v is in range(C.mulVec), then apply idempotence
  -- v has zero mean (hv), so v is in meanZeroSubmodule (after translation via equiv)
  -- By Step B, this means v is in range(C.mulVec)
  -- So v = C.mulVec u for some u
  -- Then P.mulVec v = P.mulVec (C.mulVec u) = C.mulVec u = v
  --   using the idempotence property P.mulVec (C.mulVec u) = C.mulVec u

  -- The key insight is that Cᵀ * C = 1 implies that P = C*Cᵀ is the projection
  -- onto the range of C. Since v has zero mean and range(C) = meanZeroSubmodule,
  -- v is fixed by P.
  set V : EuclideanSpace ℝ (Fin n) := (EuclideanSpace.equiv (Fin n) ℝ).symm v with hV
  have hVmem : V ∈ meanZeroSubmodule n := by
    show (inner ℝ (oneVec n) V : ℝ) = 0
    have hbridge := euclidean_equiv_dotProduct n
      ((EuclideanSpace.equiv (Fin n) ℝ).symm (fun _ => (1 : ℝ))) V
    rw [ContinuousLinearEquiv.apply_symm_apply, hV, ContinuousLinearEquiv.apply_symm_apply] at hbridge
    have hone : oneVec n = (EuclideanSpace.equiv (Fin n) ℝ).symm (fun _ => (1 : ℝ)) := rfl
    rw [hone, ← hbridge]
    simpa [dotProduct] using hv
  rw [← contrastEmbedding_range_eq_meanZeroSubmodule n hn] at hVmem
  obtain ⟨U, hU⟩ := hVmem
  set u : Fin (n - 1) → ℝ := (EuclideanSpace.equiv (Fin (n - 1)) ℝ) U with hu
  have hCu : (C_matrix n hn).mulVec u = v := by
    have hmv : (C_matrix n hn).mulVec = contrastEmbeddingPi n hn := by
      funext w
      exact LinearMap.toMatrix'_mulVec (contrastEmbeddingPi n hn) w
    rw [hmv]
    unfold contrastEmbeddingPi
    simp only [LinearMap.comp_apply, LinearEquiv.coe_toLinearMap, LinearIsometry.coe_toLinearMap]
    rw [hu]
    have hsymm : (EuclideanSpace.equiv (Fin (n - 1)) ℝ).symm
        ((EuclideanSpace.equiv (Fin (n - 1)) ℝ) U) = U :=
      (EuclideanSpace.equiv (Fin (n - 1)) ℝ).left_inv U
    have hfinal : (EuclideanSpace.equiv (Fin n) ℝ) V = v := by
      rw [hV]; exact (EuclideanSpace.equiv (Fin n) ℝ).right_inv v
    exact (congrArg (EuclideanSpace.equiv (Fin n) ℝ) (congrArg (contrastEmbedding n hn) hsymm)).trans
      ((congrArg (EuclideanSpace.equiv (Fin n) ℝ) hU).trans hfinal)
  calc (C_matrix n hn * (C_matrix n hn)ᵀ).mulVec v
      = (C_matrix n hn * (C_matrix n hn)ᵀ).mulVec ((C_matrix n hn).mulVec u) := by rw [hCu]
    _ = ((C_matrix n hn * (C_matrix n hn)ᵀ) * (C_matrix n hn)).mulVec u := by
        rw [Matrix.mulVec_mulVec]
    _ = (C_matrix n hn * ((C_matrix n hn)ᵀ * (C_matrix n hn))).mulVec u := by
        rw [Matrix.mul_assoc]
    _ = (C_matrix n hn * (1 : Matrix (Fin (n - 1)) (Fin (n - 1)) ℝ)).mulVec u := by
        rw [C_matrix_transpose_mul_self]
    _ = (C_matrix n hn).mulVec u := by rw [Matrix.mul_one]
    _ = v := hCu

#print axioms contrastEmbedding_range_eq_meanZeroSubmodule
#print axioms C_matrix_CCt_fixes_meanZero
