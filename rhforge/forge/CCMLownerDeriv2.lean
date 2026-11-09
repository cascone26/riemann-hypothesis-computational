import CCMLownerDeriv

/-!
# S2987 — second-level Löwner control via g'' (lane 2: sharp moment hierarchy, step 2)

S2985: Löwner off-diagonal entries are difference quotients, bounded (sharply, 2–4%) by sup|g'|/π.  Next level:
* `deriv_lipschitz` — |g''| ≤ D₂ on [−V,V] ⇒ |g'(x) − g'(a)| ≤ D₂|x − a|.
* `divdiff_sub_deriv` — |g b − g a − g'(a)(b − a)| ≤ D₂·(b − a)²  (MVT on the segment for h = g − g'(a)·x).
* `lowner_entry_sub_rowconst` — the Löwner entry (g b − g a)/(π(a − b)) equals −g'(a)/π up to D₂|b − a|/π.
Reading: off-diagonal Löwner tail = −(1/π)·g'(ν_i) (row-constant: u𝟙ᵀ, killed on 𝟙^⊥) + O(D₂|ν_i − ν_k|/π).
For the CCM far tail g_tail = πΣw·2ν/(ν²−x²): D₂ = sup|g''| = O(Σ w·V/x⁴) — the sharp second level.
-/

open Set

namespace CCMLownerDeriv2

theorem deriv_lipschitz (g' g'' : ℝ → ℝ) {V D₂ : ℝ}
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) {a x : ℝ} (ha : a ∈ Icc (-V) V) (hx : x ∈ Icc (-V) V) :
    |g' x - g' a| ≤ D₂ * |x - a| := by
  have h := (convex_Icc (-V) V).norm_image_sub_le_of_norm_hasDerivWithin_le hg'
    (fun y hy => by simpa [Real.norm_eq_abs] using hD y hy) ha hx
  simpa [Real.norm_eq_abs] using h

theorem divdiff_sub_deriv (g g' g'' : ℝ → ℝ) {V D₂ : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) {a b : ℝ} (ha : a ∈ Icc (-V) V) (hb : b ∈ Icc (-V) V) :
    |g b - g a - g' a * (b - a)| ≤ D₂ * (b - a) ^ 2 := by
  have hsub : uIcc a b ⊆ Icc (-V) V := uIcc_subset_Icc ha hb
  set h : ℝ → ℝ := fun x => g x - g' a * x with hh
  have hderiv : ∀ x ∈ uIcc a b, HasDerivWithinAt h (g' x - g' a) (uIcc a b) x := by
    intro x hx
    have h1 := ((hg x (hsub hx)).mono hsub)
    have h2 := (hasDerivWithinAt_id x (uIcc a b)).const_mul (g' a)
    have h3 := h1.sub h2
    rw [mul_one] at h3
    exact h3
  have hbound : ∀ x ∈ uIcc a b, ‖g' x - g' a‖ ≤ D₂ * |b - a| := by
    intro x hx
    rw [Real.norm_eq_abs]
    have hl := deriv_lipschitz g' g'' hg' hD ha (hsub hx)
    have hxa : |x - a| ≤ |b - a| := by
      rcases le_total a b with hab | hab
      · rw [uIcc_of_le hab] at hx
        rw [abs_of_nonneg (by linarith [hx.1]), abs_of_nonneg (by linarith)]; linarith [hx.2]
      · rw [uIcc_of_ge hab] at hx
        rw [abs_of_nonpos (by linarith [hx.2]), abs_of_nonpos (by linarith)]; linarith [hx.1]
    have hD0 : 0 ≤ D₂ := le_trans (abs_nonneg _) (hD a ha)
    exact le_trans hl (mul_le_mul_of_nonneg_left hxa hD0)
  have hmv := (convex_uIcc a b).norm_image_sub_le_of_norm_hasDerivWithin_le hderiv hbound
    left_mem_uIcc right_mem_uIcc
  simp only [hh, Real.norm_eq_abs] at hmv
  have e : g b - g' a * b - (g a - g' a * a) = g b - g a - g' a * (b - a) := by ring
  rw [e] at hmv
  calc |g b - g a - g' a * (b - a)| ≤ D₂ * |b - a| * |b - a| := hmv
    _ = D₂ * (b - a) ^ 2 := by rw [mul_assoc, ← sq, sq_abs]

theorem lowner_entry_sub_rowconst (g g' g'' : ℝ → ℝ) {V D₂ : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) {a b : ℝ} (ha : a ∈ Icc (-V) V) (hb : b ∈ Icc (-V) V)
    (hab : a ≠ b) :
    |(g b - g a) / (Real.pi * (a - b)) - (-(g' a) / Real.pi)| ≤ D₂ * |b - a| / Real.pi := by
  have hdd := divdiff_sub_deriv g g' g'' hg hg' hD ha hb
  have hne : b - a ≠ 0 := sub_ne_zero.mpr (Ne.symm hab)
  have hpos : 0 < |b - a| := abs_pos.mpr hne
  have hpi := Real.pi_pos
  have e : (g b - g a) / (Real.pi * (a - b)) - (-(g' a) / Real.pi)
      = -((g b - g a - g' a * (b - a)) / (b - a)) / Real.pi := by
    field_simp
    ring
  rw [e, abs_div, abs_neg, abs_div, abs_of_pos hpi, div_le_div_iff_of_pos_right hpi,
    div_le_iff₀ hpos]
  calc |g b - g a - g' a * (b - a)| ≤ D₂ * (b - a) ^ 2 := hdd
    _ = D₂ * |b - a| * |b - a| := by rw [mul_assoc, ← sq, sq_abs]

end CCMLownerDeriv2

#print axioms CCMLownerDeriv2.deriv_lipschitz
#print axioms CCMLownerDeriv2.divdiff_sub_deriv
#print axioms CCMLownerDeriv2.lowner_entry_sub_rowconst
