/-
GramEigenBound.lean (HP RH center, LEDGER S2929, 2026-11-08) — lane 2.
Quadratic-form lower bounds ⟹ eigenvalue lower bounds (the μ_λ = λ_min reading of the margin):
- `eig_ge_of_form_ge`: A Hermitian, ∀x, μ‖x‖² ≤ xᵀAx ⟹ every eigenvalue ≥ μ.
- `gram_eig_nonneg`: G = C·diag(w)·Cᵀ with w ≥ 0 ⟹ all eigenvalues ≥ 0.
- `gram_eig_stable`: certified margin μ at w ⟹ every eigenvalue of G(w') ≥ μ − Σ_j δ_j Σ_i C_ij² (|w'−w| ≤ δ).
Closes the loop S2923 (Weil form) + S2925 (truncation) + S2927 (perturbation) → spectral statement.
-/
import Mathlib
import GramPerturb

open Matrix Finset

namespace GramEig

variable {n m : Type*} [Fintype n] [Fintype m] [DecidableEq n] [DecidableEq m]

theorem eig_ge_of_form_ge (A : Matrix n n ℝ) (hA : A.IsHermitian) (μ : ℝ)
    (h : ∀ x : n → ℝ, μ * ∑ i, x i ^ 2 ≤ x ⬝ᵥ (A *ᵥ x)) (i : n) :
    μ ≤ hA.eigenvalues i := by
  rw [hA.eigenvalues_eq]
  set v : n → ℝ := ⇑(hA.eigenvectorBasis i) with hv
  have hn : ∑ j, v j ^ 2 = 1 := by
    have h1 : ‖hA.eigenvectorBasis i‖ = 1 := (hA.eigenvectorBasis).orthonormal.1 i
    rw [EuclideanSpace.norm_eq, Real.sqrt_eq_one] at h1
    simpa [hv, Real.norm_eq_abs, sq_abs] using h1
  have := h v
  rw [hn, mul_one] at this
  simpa [hv, star_trivial] using this

theorem gram_eig_nonneg (C : Matrix n m ℝ) (w : m → ℝ) (hw : ∀ j, 0 ≤ w j)
    (hG : (C * diagonal w * Cᵀ).IsHermitian) (i : n) :
    0 ≤ hG.eigenvalues i :=
  eig_ge_of_form_ge _ hG 0 (fun x => by
    rw [zero_mul]; exact GramWeil.psd_of_nonneg C w hw x) i

theorem gram_eig_stable (C : Matrix n m ℝ) (w w' δ : m → ℝ) (h : ∀ j, |w' j - w j| ≤ δ j)
    (μ : ℝ) (hμ : ∀ x : n → ℝ, μ * ∑ i, x i ^ 2 ≤ x ⬝ᵥ ((C * diagonal w * Cᵀ) *ᵥ x))
    (hG' : (C * diagonal w' * Cᵀ).IsHermitian) (i : n) :
    μ - ∑ j, δ j * ∑ k, C k j ^ 2 ≤ hG'.eigenvalues i :=
  eig_ge_of_form_ge _ hG' _ (GramPert.margin_stable C w w' δ h μ hμ) i

end GramEig

#print axioms GramEig.eig_ge_of_form_ge
#print axioms GramEig.gram_eig_nonneg
#print axioms GramEig.gram_eig_stable
