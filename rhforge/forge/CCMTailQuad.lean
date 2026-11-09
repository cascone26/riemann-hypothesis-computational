import CCMTailConst

/-!
# S2976 — quadratic-form control of the far-zero tail (lane 2)

* `near_const_quad` — if every entry of M is within δ of μ, then
  `|xᵀMx − μ(Σx)²| ≤ δ(Σ|x|)²`.
* `far_zero_tail_quad` — for the far-zero block T_ik = Σ_j w_j/((ν_i−x_j)(ν_k−x_j)) (w ≥ 0, |ν_i| ≤ V, |x_j| ≥ 2V):
  `|xᵀTx − μ(Σx)²| ≤ δ(Σ|x|)²` with `μ = Σ w_j/x_j²`, `δ = 10V·Σ w_j/|x_j|³`.
Reading (Report): combined with S2974's numerics (A − G_K ≈ μ_K𝟙𝟙ᵀ), the zeros beyond height T act on test vectors
only through the single functional x ↦ Σx (up to an O(V·Σw/|γ|³) error): on 𝟙^⊥ the tail is negligible, so the
CCM margin on 𝟙^⊥ is decided by finitely many zeros.
-/

open Finset Matrix

namespace CCMTailQuad

variable {n : Type*} [Fintype n]

theorem near_const_quad (M : Matrix n n ℝ) (μ δ : ℝ) (h : ∀ i k, |M i k - μ| ≤ δ) (x : n → ℝ) :
    |x ⬝ᵥ (M *ᵥ x) - μ * (∑ i, x i) ^ 2| ≤ δ * (∑ i, |x i|) ^ 2 := by
  have e1 : x ⬝ᵥ (M *ᵥ x) = ∑ i, ∑ k, x i * M i k * x k := by
    simp only [dotProduct, mulVec, Finset.mul_sum]
    apply Finset.sum_congr rfl; intro i _
    apply Finset.sum_congr rfl; intro k _
    ring
  have e2 : μ * (∑ i, x i) ^ 2 = ∑ i, ∑ k, x i * μ * x k := by
    rw [sq, Finset.sum_mul_sum, Finset.mul_sum]
    apply Finset.sum_congr rfl; intro i _
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl; intro k _
    ring
  have e3 : δ * (∑ i, |x i|) ^ 2 = ∑ i, ∑ k, |x i| * δ * |x k| := by
    rw [sq, Finset.sum_mul_sum, Finset.mul_sum]
    apply Finset.sum_congr rfl; intro i _
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl; intro k _
    ring
  rw [e1, e2, e3, ← Finset.sum_sub_distrib]
  refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum (fun i _ => ?_))
  rw [← Finset.sum_sub_distrib]
  refine le_trans (Finset.abs_sum_le_sum_abs _ _) (Finset.sum_le_sum (fun k _ => ?_))
  rw [show x i * M i k * x k - x i * μ * x k = x i * (M i k - μ) * x k by ring,
    abs_mul, abs_mul]
  have := h i k
  have hxi := abs_nonneg (x i)
  have hxk := abs_nonneg (x k)
  calc |x i| * |M i k - μ| * |x k| ≤ |x i| * δ * |x k| := by gcongr

theorem far_zero_tail_quad {m : Type*} [Fintype m] {V : ℝ} (ν : n → ℝ) (xz w : m → ℝ)
    (hν : ∀ i, |ν i| ≤ V) (hx : ∀ j, 2 * V ≤ |xz j|) (hx0 : ∀ j, 0 < |xz j|) (hw : ∀ j, 0 ≤ w j)
    (x : n → ℝ) :
    |x ⬝ᵥ ((fun i k => ∑ j, w j * (1 / ((ν i - xz j) * (ν k - xz j)))) *ᵥ x)
        - (∑ j, w j * (1 / xz j ^ 2)) * (∑ i, x i) ^ 2|
      ≤ (10 * V * ∑ j, w j / |xz j| ^ 3) * (∑ i, |x i|) ^ 2 :=
  near_const_quad _ _ _ (fun i k =>
    CCMTail.cauchy_tail_entry_near_const (ν i) (ν k) xz w (hν i) (hν k) hx hx0 hw) x

end CCMTailQuad

#print axioms CCMTailQuad.near_const_quad
#print axioms CCMTailQuad.far_zero_tail_quad
