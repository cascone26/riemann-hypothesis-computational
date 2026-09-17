import Mathlib

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **The Weil τ-form is symmetric.**

CM's matrix `T` (`T_{ii}=a_i`, `T_{ij}=(b_i−b_j)/(c_i−c_j)` for `i≠j`, `c` injective) satisfies
`Tᵀ = T`: on the diagonal trivially; off-diagonal because swapping `i,j` flips the sign of both
numerator and denominator.  This *derives* the real-symmetry that the master theorem
`weil_spectral_realization` otherwise takes as a hypothesis — one fewer assumption.  Axiom-clean. -/
theorem T_isSymm_of_form (c a b : n → ℝ) (hc : Function.Injective c) (T : Matrix n n ℝ)
    (hdiag : ∀ i, T i i = a i)
    (hoff : ∀ i j, i ≠ j → T i j = (b i - b j) / (c i - c j)) :
    Tᵀ = T := by
  ext i j
  simp only [Matrix.transpose_apply]
  by_cases h : i = j
  · subst h; rfl
  · have hji : j ≠ i := fun e => h e.symm
    have hcij : c i - c j ≠ 0 := sub_ne_zero.mpr (fun e => h (hc e))
    have hcji : c j - c i ≠ 0 := sub_ne_zero.mpr (fun e => hji (hc e))
    rw [hoff j i hji, hoff i j h, div_eq_div_iff hcji hcij]
    ring

#print axioms T_isSymm_of_form

end CMKeyLemma
