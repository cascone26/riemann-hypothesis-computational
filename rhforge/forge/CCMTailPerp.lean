import CCMTailPair

/-!
# S2980 — on 𝟙^⊥ the far zeros cost at most δ·n·‖x‖²: the CCM margin is decided by FINITELY many zeros (lane 2)

* `abs_sum_sq_le_card` — (Σ|x_i|)² ≤ n·Σ x_i²  (Cauchy–Schwarz).
* `near_const_quad_perp` — entries within δ of a constant, Σx = 0 ⇒ |xᵀMx| ≤ δ·n·Σx².
* `margin_transfer_perp` — if xᵀAx = xᵀGx + xᵀTx with xᵀGx ≥ m·Σx² and T near-constant (δ), then for Σx = 0:
  xᵀAx ≥ (m − δ·n)·Σx².
Instantiated with S2978 (δ = 11V²Σw/x⁴ from the ±γ-paired far zeros), the CCM margin on 𝟙^⊥ is bounded below by the
margin of the FINITE near-zero Gram minus an explicit tail term.
-/

open Finset Matrix

namespace CCMTailPerp

variable {n : Type*} [Fintype n]

theorem abs_sum_sq_le_card (x : n → ℝ) :
    (∑ i, |x i|) ^ 2 ≤ (Fintype.card n : ℝ) * ∑ i, x i ^ 2 := by
  have h := sq_sum_le_card_mul_sum_sq (s := (univ : Finset n)) (f := fun i => |x i|)
  simpa [sq_abs, Finset.card_univ] using h

theorem near_const_quad_perp (M : Matrix n n ℝ) (μ δ : ℝ) (hδ : 0 ≤ δ)
    (h : ∀ i k, |M i k - μ| ≤ δ) (x : n → ℝ) (hx : ∑ i, x i = 0) :
    |x ⬝ᵥ (M *ᵥ x)| ≤ δ * (Fintype.card n : ℝ) * ∑ i, x i ^ 2 := by
  have hq := CCMTailQuad.near_const_quad M μ δ h x
  rw [hx] at hq
  simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero, sub_zero] at hq
  have hcs := abs_sum_sq_le_card x
  calc |x ⬝ᵥ (M *ᵥ x)| ≤ δ * (∑ i, |x i|) ^ 2 := hq
    _ ≤ δ * ((Fintype.card n : ℝ) * ∑ i, x i ^ 2) := mul_le_mul_of_nonneg_left hcs hδ
    _ = δ * (Fintype.card n : ℝ) * ∑ i, x i ^ 2 := by ring

theorem margin_transfer_perp (A G T : Matrix n n ℝ) (m μ δ : ℝ) (hδ : 0 ≤ δ)
    (hT : ∀ i k, |T i k - μ| ≤ δ) (x : n → ℝ) (hx : ∑ i, x i = 0)
    (hsplit : x ⬝ᵥ (A *ᵥ x) = x ⬝ᵥ (G *ᵥ x) + x ⬝ᵥ (T *ᵥ x))
    (hG : m * ∑ i, x i ^ 2 ≤ x ⬝ᵥ (G *ᵥ x)) :
    (m - δ * (Fintype.card n : ℝ)) * ∑ i, x i ^ 2 ≤ x ⬝ᵥ (A *ᵥ x) := by
  have h := near_const_quad_perp T μ δ hδ hT x hx
  have h' := neg_abs_le (x ⬝ᵥ (T *ᵥ x))
  rw [hsplit]
  nlinarith

end CCMTailPerp

#print axioms CCMTailPerp.abs_sum_sq_le_card
#print axioms CCMTailPerp.near_const_quad_perp
#print axioms CCMTailPerp.margin_transfer_perp
