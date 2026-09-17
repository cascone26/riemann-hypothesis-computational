import Mathlib

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **CM 2025 Lemma "basics" (i): the Weil matrix commutes with the grading `γ`.**

`γ : V_j ↦ V_{-j}` is the index-reversal involution `r` (`c(r j) = −c(j)`).  Commuting with `γ`
is `T(r i, r j) = T(i, j)`, i.e. `T.submatrix r r = T`.  For CM's τ-form with `a` even
(`a(r j)=a j`) and `b`, `c` odd (`b(r j)=−b j`, `c(r j)=−c j`), this holds: on the diagonal by
evenness of `a`; off-diagonal because the two sign flips in `(b_i−b_j)/(c_i−c_j)` cancel.

Together with `QD_of_loaded_form` and real-symmetry, this shows the concrete Weil matrix satisfies
every *structural* hypothesis of the abstract key lemma (`γ`-commuting + the QD commutator) — the only
remaining hypothesis is CM's open **even-simple**.  Unconditional, machine-verified, axiom-clean. -/
theorem T_submatrix_reversal (c a b : n → ℝ) (r : n → n) (T : Matrix n n ℝ)
    (hr : Function.Involutive r) (hc : Function.Injective c)
    (hae : ∀ j, a (r j) = a j) (hbo : ∀ j, b (r j) = -b j) (hco : ∀ j, c (r j) = -c j)
    (hdiag : ∀ i, T i i = a i)
    (hoff : ∀ i j, i ≠ j → T i j = (b i - b j) / (c i - c j)) :
    T.submatrix r r = T := by
  ext i j
  simp only [Matrix.submatrix_apply]
  by_cases h : i = j
  · subst h; rw [hdiag (r i), hdiag i, hae i]
  · have hrij : r i ≠ r j := fun e => h (hr.injective e)
    have hcij : c i - c j ≠ 0 := sub_ne_zero.mpr (fun e => h (hc e))
    have hcij2 : -c i - -c j ≠ 0 := by
      have : -c i - -c j = -(c i - c j) := by ring
      rw [this]; exact neg_ne_zero.mpr hcij
    rw [hoff (r i) (r j) hrij, hoff i j h, hbo, hbo, hco, hco, div_eq_div_iff hcij2 hcij]
    ring

#print axioms T_submatrix_reversal

end CMKeyLemma
