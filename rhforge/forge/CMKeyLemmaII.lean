import CMKeyLemmaI

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **CM 2025 (ccm_2511.22755) Lemma "key" part (ii), matrix core.**

With the setup of part (i) (`D`,`T` real symmetric, commutator `D·T − T·D = |β⟩⟨η| − |η⟩⟨β|`,
`T·ξ = 0`, `⟨η,ξ⟩ = 1`, `⟨β,ξ⟩ = 0`), CM's perturbed operator is `D' = D − |Dξ⟩⟨η|`.  Their
statement that `D'` descends to a **self-adjoint** operator on the Hilbert space carrying the
`T`-inner product `⟨f,g⟩_T = ⟨Tf,g⟩` is, at the matrix level, exactly the symmetry of `T·D'`:

    `(T · D')ᵀ = T · D'`.

Proof (CM §872): `T·D' = T·D + |β⟩⟨η|` (using part (i), `T·Dξ = −β`), whose transpose is
`D·T + |η⟩⟨β|`, and these agree by the commutator hypothesis.  Machine-verified, axiom-clean;
consumes `CMKeyLemma.T_mulVec_D_xi` (part i). -/
theorem T_mul_Dprime_isSymm (D T : Matrix n n ℝ) (β η ξ : n → ℝ)
    (hDsym : Dᵀ = D) (hTsym : Tᵀ = T)
    (hcomm : D * T - T * D = Matrix.vecMulVec β η - Matrix.vecMulVec η β)
    (hTxi : T *ᵥ ξ = 0) (hη : η ⬝ᵥ ξ = 1) (hβ : β ⬝ᵥ ξ = 0) :
    (T * (D - Matrix.vecMulVec (D *ᵥ ξ) η))ᵀ = T * (D - Matrix.vecMulVec (D *ᵥ ξ) η) := by
  have hTDξ : T *ᵥ (D *ᵥ ξ) = -β := CMKeyLemma.T_mulVec_D_xi D T β η ξ hcomm hTxi hη hβ
  -- T·D' = T·D + |β⟩⟨η|
  have hTDp : T * (D - Matrix.vecMulVec (D *ᵥ ξ) η) = T * D + Matrix.vecMulVec β η := by
    rw [Matrix.mul_sub, Matrix.mul_vecMulVec, hTDξ, Matrix.neg_vecMulVec, sub_neg_eq_add]
  -- commutator ⟹ T·D + |β⟩⟨η| = D·T + |η⟩⟨β|
  have key : T * D + Matrix.vecMulVec β η = D * T + Matrix.vecMulVec η β := by
    have h := sub_eq_iff_eq_add.mp hcomm
    rw [h]; abel
  rw [hTDp, Matrix.transpose_add, Matrix.transpose_mul, Matrix.transpose_vecMulVec, hDsym, hTsym]
  exact key.symm

#print axioms T_mul_Dprime_isSymm

end CMKeyLemma
