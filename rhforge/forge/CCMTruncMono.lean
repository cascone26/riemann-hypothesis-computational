import CCMResonance
import CCMGramPosDef

/-!
# S2941 — truncation monotonicity of the CCM zero side (lane 2)

Truncating the zero side to a finite set S of ordinates is the Gram with weights w·1_S.  Since 1_S·w ≤ 1_T·w for
S ⊆ T and w ≥ 0, S2937's Loewner monotonicity gives:
* `gram_trunc_mono` : G_S ⪯ G_T  (adding zeros can only increase the zero-side Gram);
* `ccm_trunc_rayleigh_mono` : xᵀG_S x ≤ xᵀG_T x for the CCM weights;
* `posDef_of_posDef_le` : if G_S ≻ 0 and G_S ⪯ G_T then G_T ≻ 0 — strict positivity certified on a finite
  truncation PROPAGATES to every larger truncation (the direction in which the Arb certificates are used).
-/

open Matrix Finset

namespace CCMTrunc

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

/-- Weights restricted to a finite set of ordinate indices. -/
def truncW (S : Finset m) (w : m → ℝ) : m → ℝ := fun j => if j ∈ S then w j else 0

theorem truncW_mono {S T : Finset m} (hST : S ⊆ T) {w : m → ℝ} (hw : ∀ j, 0 ≤ w j) (j : m) :
    truncW S w j ≤ truncW T w j := by
  unfold truncW
  by_cases hS : j ∈ S
  · simp [hS, hST hS]
  · by_cases hT : j ∈ T <;> simp [hS, hT, hw j]

theorem gram_trunc_mono (C : Matrix n m ℝ) {S T : Finset m} (hST : S ⊆ T) {w : m → ℝ}
    (hw : ∀ j, 0 ≤ w j) :
    (C * diagonal (truncW T w) * Cᵀ - C * diagonal (truncW S w) * Cᵀ).PosSemidef :=
  CCMLoewner.gram_mono C (truncW_mono hST hw)

theorem ccm_trunc_rayleigh_mono {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ)
    {S T : Finset m} (hST : S ⊆ T) (x : n → ℝ) :
    x ⬝ᵥ ((C * diagonal (truncW S (fun j => CCMWeight.ccmW L (γ j))) * Cᵀ) *ᵥ x)
      ≤ x ⬝ᵥ ((C * diagonal (truncW T (fun j => CCMWeight.ccmW L (γ j))) * Cᵀ) *ᵥ x) :=
  CCMResonance.rayleigh_mono C (truncW_mono hST (fun j => CCMWeight.ccmW_nonneg hL (γ j))) x

theorem posDef_of_posDef_le {A B : Matrix n n ℝ} (hA : A.PosDef) (hBA : (B - A).PosSemidef) :
    B.PosDef := by
  have := hA.add_posSemidef hBA
  simpa using this

/-- A strictly positive truncated CCM zero side stays strictly positive under any enlargement. -/
theorem ccm_trunc_posDef_mono {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ)
    {S T : Finset m} (hST : S ⊆ T)
    (hS : (C * diagonal (truncW S (fun j => CCMWeight.ccmW L (γ j))) * Cᵀ).PosDef) :
    (C * diagonal (truncW T (fun j => CCMWeight.ccmW L (γ j))) * Cᵀ).PosDef :=
  posDef_of_posDef_le hS (gram_trunc_mono C hST (fun j => CCMWeight.ccmW_nonneg hL (γ j)))

end CCMTrunc

#print axioms CCMTrunc.gram_trunc_mono
#print axioms CCMTrunc.ccm_trunc_rayleigh_mono
#print axioms CCMTrunc.ccm_trunc_posDef_mono
