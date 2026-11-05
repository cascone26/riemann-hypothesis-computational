import Mathlib

open scoped MatrixOrder Matrix RealInnerProductSpace

/-- **The zero-mean contrast subspace `𝒱⁰_{N,L}` (Shi Eq 4.5, arXiv:2609.04908)**: vectors in
`EuclideanSpace ℝ (Fin n)` orthogonal to the all-ones vector, i.e. with coordinate-sum zero — the
subspace any valid `C_{N,L}` (Helmert or otherwise) must have as its column-span, per the paper's
own remark "Any other basis ... gives a congruent matrix pencil." Fully general, no zeta/RH
dependence. -/
noncomputable def oneVec (n : ℕ) : EuclideanSpace ℝ (Fin n) :=
  (EuclideanSpace.equiv (Fin n) ℝ).symm (fun _ => (1 : ℝ))

noncomputable def meanZeroFunctional (n : ℕ) : Module.Dual ℝ (EuclideanSpace ℝ (Fin n)) :=
  innerSL ℝ (oneVec n)

noncomputable def meanZeroSubmodule (n : ℕ) : Submodule ℝ (EuclideanSpace ℝ (Fin n)) :=
  LinearMap.ker (meanZeroFunctional n)

/-- The all-ones functional is nonzero (for `n ≥ 1`), since it evaluates to `n ≠ 0` on itself. -/
theorem meanZeroFunctional_ne_zero (n : ℕ) (hn : 1 ≤ n) :
    meanZeroFunctional n ≠ 0 := by
  intro h
  have heval : meanZeroFunctional n (oneVec n) = 0 := by rw [h]; rfl
  have hinner : ⟪oneVec n, oneVec n⟫ = 0 := by
    simpa [meanZeroFunctional, innerSL_apply_apply] using heval
  have hzero : oneVec n = 0 := inner_self_eq_zero.mp hinner
  have hcoord := congrFun (congrArg (EuclideanSpace.equiv (Fin n) ℝ) hzero) ⟨0, hn⟩
  simp [oneVec] at hcoord

/-- **The zero-mean subspace has dimension `n − 1`** (for `n ≥ 1`): via `Module.Dual.finrank_ker_
add_one_of_ne_zero`, the kernel of a nonzero functional on an `n`-dimensional space has dimension
`n − 1`. Fully general, no zeta/RH dependence. -/
theorem meanZeroSubmodule_finrank (n : ℕ) (hn : 1 ≤ n) :
    Module.finrank ℝ (meanZeroSubmodule n) = n - 1 := by
  unfold meanZeroSubmodule
  have h := Module.Dual.finrank_ker_add_one_of_ne_zero (meanZeroFunctional_ne_zero n hn)
  have hfin : Module.finrank ℝ (EuclideanSpace ℝ (Fin n)) = n :=
    finrank_euclideanSpace_fin
  rw [hfin] at h
  omega

#print axioms meanZeroFunctional_ne_zero
#print axioms meanZeroSubmodule_finrank
