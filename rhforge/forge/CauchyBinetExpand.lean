/-
CauchyBinetExpand.lean (HP RH center, LEDGER S2916, 2026-11-08)
Mathlib has no Cauchy–Binet formula. Step 1 of building it: the function-sum expansion
  det (A * B) = Σ_{p : m → n} det (A with columns p) * Π_i B (p i) i,
for A : m×n, B : n×m (rectangular), and the vanishing of non-injective terms.
Target (S2917+): sum over strictly-monotone p ⟹ Cauchy–Binet ⟹ TP of C·W·Cᵀ (S2859 band theorem).
-/
import Mathlib

open Matrix Finset Equiv

namespace CauchyBinet

variable {m n R : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n] [CommRing R]

/-- Function-sum expansion of det of a rectangular product. -/
theorem det_mul_expand (A : Matrix m n R) (B : Matrix n m R) :
    det (A * B) = ∑ p : m → n, det (A.submatrix id p) * ∏ i, B (p i) i := by
  calc det (A * B) = ∑ p : m → n, ∑ σ : Perm m,
        ((Perm.sign σ : ℤ) : R) * ∏ i, A (σ i) (p i) * B (p i) i := by
        simp only [det_apply', mul_apply, prod_univ_sum, mul_sum, Fintype.piFinset_univ]
        rw [Finset.sum_comm]
    _ = ∑ p : m → n, det (A.submatrix id p) * ∏ i, B (p i) i := by
        refine sum_congr rfl fun p _ => ?_
        rw [det_apply', Finset.sum_mul]
        refine sum_congr rfl fun σ _ => ?_
        simp only [submatrix_apply, id, prod_mul_distrib, mul_assoc]

/-- Non-injective column selections contribute zero. -/
theorem det_submatrix_eq_zero_of_not_inj (A : Matrix m n R) (p : m → n)
    (hp : ¬ Function.Injective p) : det (A.submatrix id p) = 0 := by
  obtain ⟨i, j, hij, hne⟩ : ∃ i j, p i = p j ∧ i ≠ j := by
    by_contra h
    exact hp fun i j hij => by_contra fun hne => h ⟨i, j, hij, hne⟩
  exact det_zero_of_column_eq hne (by intro k; simp [submatrix_apply, hij])

/-- Expansion restricted to injective selections. -/
theorem det_mul_expand_inj (A : Matrix m n R) (B : Matrix n m R) :
    det (A * B) = ∑ p ∈ univ.filter (fun p : m → n => Function.Injective p),
      det (A.submatrix id p) * ∏ i, B (p i) i := by
  rw [det_mul_expand]
  refine (sum_subset (filter_subset _ _) fun p _ hp => ?_).symm
  rw [det_submatrix_eq_zero_of_not_inj A p (by simpa using hp), zero_mul]

end CauchyBinet

#print axioms CauchyBinet.det_mul_expand
#print axioms CauchyBinet.det_submatrix_eq_zero_of_not_inj
#print axioms CauchyBinet.det_mul_expand_inj
