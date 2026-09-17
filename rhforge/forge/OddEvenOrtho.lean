import Mathlib

open Matrix
open scoped BigOperators

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **Odd ⟂ even under the grading `γ` (CM's `⟨β,ξ⟩ = 0`).**

If `b` is `γ`-odd (`b(r j) = −b j`) and `ξ` is `γ`-even (`ξ(r j) = ξ j`) for the reversal involution
`r`, then `⟨b, ξ⟩ = 0`.  This is CM's "the two eigenspaces of `γ` are orthogonal", which discharges the
`hβ : β ⬝ᵥ ξ = 0` hypothesis of `T_mulVec_D_xi`/`T_mul_Dprime_isSymm` for the concrete Weil matrix
(where `β = b` is odd and the ground state `ξ` is even).  Machine-verified, axiom-clean. -/
theorem dotProduct_odd_even_eq_zero (r : n → n) (b ξ : n → ℝ)
    (hr : Function.Involutive r) (hb : ∀ j, b (r j) = -b j) (hξ : ∀ j, ξ (r j) = ξ j) :
    b ⬝ᵥ ξ = 0 := by
  have h1 : ∑ j, b (r j) * ξ (r j) = ∑ j, b j * ξ j := by
    have := Equiv.sum_comp (Function.Involutive.toPerm r hr) (fun j => b j * ξ j)
    simpa [Function.Involutive.coe_toPerm] using this
  have h2 : (∑ j, b (r j) * ξ (r j)) + (∑ j, b j * ξ j) = 0 := by
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_eq_zero
    intro j _
    rw [hb, hξ]; ring
  rw [h1] at h2
  have h3 : (∑ j, b j * ξ j) = 0 := by linarith
  exact h3

#print axioms dotProduct_odd_even_eq_zero

end CMKeyLemma
