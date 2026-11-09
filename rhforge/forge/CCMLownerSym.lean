import CCMLownerDeriv2

/-!
# S2990 — symmetrised second level: on 𝟙^⊥ a Löwner tail costs only D₂-sized terms (lane 2)

* `lowner_entry_sym` — from S2987 applied at both endpoints:
  |(g b − g a)/(π(a−b)) + (g'(a) + g'(b))/(2π)| ≤ D₂|b − a|/π.
* `quad_rowsum_zero` — x ⊥ 𝟙 ⇒ xᵀ(u𝟙ᵀ + 𝟙uᵀ)x = 0.
* `near_rowsum_quad_perp` — |M_ik − (u_i + u_k)| ≤ δ ∀ i k and Σx = 0 ⇒ |xᵀMx| ≤ δ(Σ|x|)².
With u_i = −g'(ν_i)/(2π), δ = 2V·D₂/π (|ν_i − ν_k| ≤ 2V) for the off-diagonal part: on 𝟙^⊥ the far-tail Löwner block costs
O(V·sup|g''|) instead of O(sup|g'|) — the sharp second level (S2985 showed the first level is sharp to 2–4%).
-/

open Set Finset Matrix

namespace CCMLownerSym

theorem lowner_entry_sym (g g' g'' : ℝ → ℝ) {V D₂ : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) {a b : ℝ} (ha : a ∈ Icc (-V) V) (hb : b ∈ Icc (-V) V)
    (hab : a ≠ b) :
    |(g b - g a) / (Real.pi * (a - b)) + (g' a + g' b) / (2 * Real.pi)| ≤ D₂ * |b - a| / Real.pi := by
  have h1 := CCMLownerDeriv2.lowner_entry_sub_rowconst g g' g'' hg hg' hD ha hb hab
  have h2 := CCMLownerDeriv2.lowner_entry_sub_rowconst g g' g'' hg hg' hD hb ha (Ne.symm hab)
  have e : (g a - g b) / (Real.pi * (b - a)) = (g b - g a) / (Real.pi * (a - b)) := by
    rw [show g a - g b = -(g b - g a) by ring, show b - a = -(a - b) by ring, mul_neg, neg_div_neg_eq]
  rw [e, abs_sub_comm a b] at h2
  have hpi := Real.pi_pos
  have e2 : (g b - g a) / (Real.pi * (a - b)) + (g' a + g' b) / (2 * Real.pi)
      = ((g b - g a) / (Real.pi * (a - b)) - (-(g' a) / Real.pi)
        + ((g b - g a) / (Real.pi * (a - b)) - (-(g' b) / Real.pi))) / 2 := by
    field_simp; ring
  rw [e2, abs_div, abs_two]
  have := abs_add_le ((g b - g a) / (Real.pi * (a - b)) - (-(g' a) / Real.pi))
    ((g b - g a) / (Real.pi * (a - b)) - (-(g' b) / Real.pi))
  linarith

variable {n : Type*} [Fintype n]

theorem quad_rowsum_zero (u x : n → ℝ) (hx : ∑ i, x i = 0) :
    x ⬝ᵥ ((vecMulVec u (fun _ => (1 : ℝ)) + vecMulVec (fun _ => (1 : ℝ)) u) *ᵥ x) = 0 := by
  have h1 : (fun _ => (1 : ℝ)) ⬝ᵥ x = 0 := by simpa [dotProduct] using hx
  rw [Matrix.add_mulVec, dotProduct_add, CCMMoment2.quad_vecMulVec, CCMMoment2.quad_vecMulVec, h1]
  ring

theorem near_rowsum_quad_perp (M : Matrix n n ℝ) (u : n → ℝ) (δ : ℝ)
    (h : ∀ i k, |M i k - (u i + u k)| ≤ δ) (x : n → ℝ) (hx : ∑ i, x i = 0) :
    |x ⬝ᵥ (M *ᵥ x)| ≤ δ * (∑ i, |x i|) ^ 2 := by
  set P : Matrix n n ℝ := vecMulVec u (fun _ => (1 : ℝ)) + vecMulVec (fun _ => (1 : ℝ)) u with hP
  have hPe : ∀ i k, P i k = u i + u k := by intro i k; simp [hP, vecMulVec_apply]
  have hq := CCMMoment2.near_matrix_quad M P δ (fun i k => by rw [hPe]; exact h i k) x
  rw [quad_rowsum_zero u x hx, sub_zero] at hq
  exact hq

end CCMLownerSym

#print axioms CCMLownerSym.lowner_entry_sym
#print axioms CCMLownerSym.quad_rowsum_zero
#print axioms CCMLownerSym.near_rowsum_quad_perp
