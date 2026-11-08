import Mathlib
import GramWeilForm
import GramEigenBound

/-!
# S2931 — the CCM/Shi zero-side weights instantiated in the Gram framework (lane 2)

CCM τ ≡ Shi S_{N,L}: zero side S_mn = Σ_ρ w_ρ /((ν_m+γ)(ν_n+γ)) with w_γ = 4·sin²(γL/2)/L.
Here the weight is formalized and plugged into GramWeil / GramEig:
* `ccmW_nonneg`, `ccmW_le` : 0 ≤ w_γ ≤ 4/L (L > 0);  `ccmW_neg` : evenness in γ.
* `ccm_gram_posSemidef` : C·diag(w_γ)·Cᵀ is PSD for ANY real ordinates γ_j (finite truncation).
* `ccm_gram_eig_nonneg` : every eigenvalue of that Gram is ≥ 0.
* `ccm_form_le` : uniform upper bracket xᵀGx ≤ (4/L)·Σ_j ((Cᵀx)_j)².
Scope (Charisma): RH-conditional finite-truncation statement about the zero side; the arithmetic side
(W02 − W_R − W_p) and the identification τ = zero side are not formalized here.
-/

open Matrix Finset

namespace CCMWeight

noncomputable def ccmW (L γ : ℝ) : ℝ := 4 * Real.sin (γ * L / 2) ^ 2 / L

theorem ccmW_nonneg {L : ℝ} (hL : 0 < L) (γ : ℝ) : 0 ≤ ccmW L γ := by
  unfold ccmW; positivity

theorem ccmW_le {L : ℝ} (hL : 0 < L) (γ : ℝ) : ccmW L γ ≤ 4 / L := by
  unfold ccmW
  apply div_le_div_of_nonneg_right _ hL.le
  have h := Real.sin_sq_le_one (γ * L / 2)
  nlinarith

theorem ccmW_neg (L γ : ℝ) : ccmW L (-γ) = ccmW L γ := by
  unfold ccmW
  rw [show -γ * L / 2 = -(γ * L / 2) by ring, Real.sin_neg, neg_sq]

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem ccm_gram_posSemidef {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    (C * diagonal (fun j => ccmW L (γ j)) * Cᵀ).PosSemidef :=
  GramWeil.posSemidef_of_nonneg C _ (fun j => ccmW_nonneg hL (γ j))

theorem ccm_gram_hermitian {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) :
    (C * diagonal (fun j => ccmW L (γ j)) * Cᵀ).IsHermitian :=
  (ccm_gram_posSemidef hL C γ).1

theorem ccm_gram_eig_nonneg {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) (i : n) :
    0 ≤ (ccm_gram_hermitian hL C γ).eigenvalues i :=
  GramEig.gram_eig_nonneg C _ (fun j => ccmW_nonneg hL (γ j)) _ i

theorem ccm_form_le {L : ℝ} (hL : 0 < L) (C : Matrix n m ℝ) (γ : m → ℝ) (x : n → ℝ) :
    x ⬝ᵥ ((C * diagonal (fun j => ccmW L (γ j)) * Cᵀ) *ᵥ x)
      ≤ 4 / L * ∑ j, ((Cᵀ *ᵥ x) j) ^ 2 := by
  rw [GramWeil.quad_form, mul_sum]
  exact sum_le_sum fun j _ =>
    mul_le_mul_of_nonneg_right (ccmW_le hL (γ j)) (sq_nonneg _)

end CCMWeight

#print axioms CCMWeight.ccmW_neg
#print axioms CCMWeight.ccm_gram_eig_nonneg
#print axioms CCMWeight.ccm_form_le
