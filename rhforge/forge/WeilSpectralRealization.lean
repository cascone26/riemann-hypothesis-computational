import QDfromForm
import OddEvenOrtho
import TSymmOfForm
import CMKeyLemmaII
import DetDprimeDiagonal

open Matrix
open scoped BigOperators

namespace CMKeyLemma

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **MASTER THEOREM — Connes–Moscovici's spectral realization for the concrete Weil matrix.**

Let `D = diagonal c` be the scaling operator (`c` the distinct eigenvalues, `c_j = j`) and `T` the
Weil quadratic form `QW_λ^N` in CM's τ-form (`T_{ii}=a_i`, `T_{ij}=(b_i−b_j)/(c_i−c_j)`), with `a` even
/ `b`,`c` odd under the reversal `r`.  Assume the **even-simple normalization** — CM's own open
missing-step, taken here as the *only* external hypothesis — namely `T·ξ = 0` (ξ the ground state)
with `ξ` even and `⟨η,ξ⟩ = 1` (`η ≡ 1`).  Then for `D' = D − |Dξ⟩⟨η|` and any `s` off the spectrum
(`c_j − s ≠ 0`):

  1. `(T·D')ᵀ = T·D'`  — `D''` is **self-adjoint** in the `T`-inner product;
  2. `det(D'−s·1) = −s·det(D−s·1)·Σ_j (c_j−s)⁻¹ ξ_j`  — its **spectrum is the zeros of `ξ̂(s)`**.

ALL structural hypotheses — real-symmetry (`T_isSymm_of_form`), the QD commutator (`QD_of_loaded_form`),
`⟨β,ξ⟩=0` (`dotProduct_odd_even_eq_zero`), `D` symmetric — are discharged internally from the τ-form and
the grading.  The ONLY external input is even-simple.  This is the complete, honest statement of CM's
construction for the real Weil matrix, machine-verified and axiom-clean; the single open hinge is
even-simple, exactly as CM leave it. -/
theorem weil_spectral_realization (c a b ξ : n → ℝ) (r : n → n) (s : ℝ) (T : Matrix n n ℝ)
    (hr : Function.Involutive r) (hc : Function.Injective c)
    (hdiag : ∀ i, T i i = a i) (hoff : ∀ i j, i ≠ j → T i j = (b i - b j) / (c i - c j))
    (hbo : ∀ j, b (r j) = -b j) (hξe : ∀ j, ξ (r j) = ξ j)
    (hTxi : T *ᵥ ξ = 0) (hη : (fun _ => (1 : ℝ)) ⬝ᵥ ξ = 1)
    (hg : ∀ j, c j - s ≠ 0) :
    (T * (diagonal c - Matrix.vecMulVec (diagonal c *ᵥ ξ) (fun _ => 1)))ᵀ
        = T * (diagonal c - Matrix.vecMulVec (diagonal c *ᵥ ξ) (fun _ => 1))
      ∧ (diagonal c - Matrix.vecMulVec (diagonal c *ᵥ ξ) (fun _ => 1)
            - s • (1 : Matrix n n ℝ)).det
        = -s * (diagonal c - s • (1 : Matrix n n ℝ)).det
            * (∑ j, (fun _ => (1 : ℝ)) j * ((c j - s)⁻¹ * ξ j)) := by
  have hTsym : Tᵀ = T := T_isSymm_of_form c a b hc T hdiag hoff
  have hcomm := QD_of_loaded_form c a b hc T hdiag hoff
  have hβ := dotProduct_odd_even_eq_zero r b ξ hr hbo hξe
  have hDsym : (diagonal c)ᵀ = diagonal c := Matrix.diagonal_transpose c
  refine ⟨?_, ?_⟩
  · exact T_mul_Dprime_isSymm (diagonal c) T b (fun _ => 1) ξ hDsym hTsym hcomm hTxi hη hβ
  · exact detDprime_diagonal c (fun _ => 1) ξ s hg hη

#print axioms weil_spectral_realization

end CMKeyLemma
