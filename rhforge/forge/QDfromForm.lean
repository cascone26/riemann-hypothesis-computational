import Mathlib

open Matrix

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **CM 2025 Lemma "basics" (ii): the QD commutator relation from the "loaded" matrix form.**

CM's Weil matrix `T = QW_λ^N` has the special form (their eq. `form`):
`T_{ii} = a_i` and `T_{ij} = (b_i − b_j)/(c_i − c_j)` for `i ≠ j`, where `c` are the (distinct)
eigenvalues of the scaling operator `D = diagonal c` (`c_i = i`).  For any such `T`,

    `D·T − T·D = |β⟩⟨η| − |η⟩⟨β|`,  with `β = b`, `η ≡ 1`.

This is the **concrete→abstract bridge**: it discharges the commutator hypothesis (`hcomm`) that the
entire abstract key-lemma chain (`T_mulVec_D_xi`, `T_mul_Dprime_isSymm`, `detDprime_core`, …) takes as
input — so those results apply to the *actual* Weil matrix.  Unconditional (needs only the form and
`c` injective).  Machine-verified, axiom-clean. -/
theorem QD_of_loaded_form (c a b : n → ℝ) (hc : Function.Injective c) (T : Matrix n n ℝ)
    (hdiag : ∀ i, T i i = a i)
    (hoff : ∀ i j, i ≠ j → T i j = (b i - b j) / (c i - c j)) :
    diagonal c * T - T * diagonal c
      = Matrix.vecMulVec b (fun _ => 1) - Matrix.vecMulVec (fun _ => 1) b := by
  ext i j
  simp only [Matrix.sub_apply, Matrix.diagonal_mul, Matrix.mul_diagonal, Matrix.vecMulVec_apply,
    mul_one, one_mul]
  by_cases h : i = j
  · subst h; ring
  · rw [hoff i j h]
    have hcij : c i - c j ≠ 0 := sub_ne_zero.mpr (fun e => h (hc e))
    field_simp

#print axioms QD_of_loaded_form

end CMKeyLemma
