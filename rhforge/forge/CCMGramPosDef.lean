import CCMWeightGram

/-!
# S2933 — strict positivity of the CCM/Shi zero-side Gram (lane 2)

Upgrades S2931 (PSD) to POSITIVE DEFINITE: if no zero ordinate is a resonance of the window
(sin(γ_j L/2) ≠ 0 for all j) and the Cauchy columns separate vectors (Cᵀ injective, e.g. ≥ n distinct
ordinates with distinct shifts ν), then G = C·diag(w_γ)·Cᵀ is PosDef and every eigenvalue is > 0 —
i.e. the zero-side margin is strictly positive on that truncation.
* `ccmW_pos` : sin(γL/2) ≠ 0 ⇒ w_γ > 0.
* `ccm_gram_posDef`, `ccm_gram_eig_pos`.
Scope (Charisma): RH-conditional (γ real), finite truncation, zero side only.
-/

open Matrix Finset

namespace CCMWeight

theorem ccmW_pos {L γ : ℝ} (hL : 0 < L) (hs : Real.sin (γ * L / 2) ≠ 0) : 0 < ccmW L γ := by
  unfold ccmW
  have : 0 < Real.sin (γ * L / 2) ^ 2 := by positivity
  positivity

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem ccm_gram_posDef {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ)
    (hs : ∀ j, Real.sin (γ j * L / 2) ≠ 0) (hinj : Function.Injective (Cᵀ).mulVec) :
    (C * diagonal (fun j => ccmW L (γ j)) * Cᵀ).PosDef := by
  refine PosDef.of_dotProduct_mulVec_pos (ccm_gram_hermitian hL C γ) fun x hx => ?_
  simpa [star_trivial] using
    GramWeil.pos_of_pos C _ (fun j => ccmW_pos hL (hs j)) hinj x hx

theorem ccm_gram_eig_pos {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ)
    (hs : ∀ j, Real.sin (γ j * L / 2) ≠ 0) (hinj : Function.Injective (Cᵀ).mulVec) (i : n) :
    0 < (ccm_gram_posDef hL C γ hs hinj).1.eigenvalues i :=
  (ccm_gram_posDef hL C γ hs hinj).eigenvalues_pos i

end CCMWeight

#print axioms CCMWeight.ccm_gram_posDef
#print axioms CCMWeight.ccm_gram_eig_pos
