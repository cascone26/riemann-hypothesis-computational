import CCMTaylorHalf
import CCMLownerSym

/-!
# S2998 — symmetrised Löwner entry with Taylor's ½ at BOTH ends (lane 2: halves S2990's constant)

* `taylor_half_right` — a ≤ b ⇒ |g a − g b − g'(b)(a − b)| ≤ D₂(b − a)²/2  (S2996 applied to h(x) = g(−x)).
* `lowner_entry_sym_half` — |(g b − g a)/(π(a−b)) + (g'(a)+g'(b))/(2π)| ≤ D₂|b − a|/(2π).
Plugged into S2992/S2994: the 𝟙^⊥ far-tail bound constant halves (Schur version: ~20× → ~10× loose).
-/

open Set

namespace CCMSymHalf

theorem taylor_half_right (g g' g'' : ℝ → ℝ) {V D₂ : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) {a b : ℝ} (ha : a ∈ Icc (-V) V) (hb : b ∈ Icc (-V) V)
    (hab : a ≤ b) :
    |g a - g b - g' b * (a - b)| ≤ D₂ * (b - a) ^ 2 / 2 := by
  have hmaps : MapsTo (fun x : ℝ => -x) (Icc (-V) V) (Icc (-V) V) := by
    intro x hx; exact ⟨by linarith [hx.2], by linarith [hx.1]⟩
  have hneg : ∀ x ∈ Icc (-V) V, -x ∈ Icc (-V) V := fun x hx => hmaps hx
  set h : ℝ → ℝ := fun x => g (-x) with hh
  set h' : ℝ → ℝ := fun x => -g' (-x) with hh'
  set h'' : ℝ → ℝ := fun x => g'' (-x) with hh''
  have hdh : ∀ x ∈ Icc (-V) V, HasDerivWithinAt h (h' x) (Icc (-V) V) x := by
    intro x hx
    have := (hg (-x) (hneg x hx)).scomp x ((hasDerivAt_neg x).hasDerivWithinAt) hmaps
    exact this.congr_deriv (by simp [hh'])
  have hdh' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt h' (h'' x) (Icc (-V) V) x := by
    intro x hx
    have := ((hg' (-x) (hneg x hx)).scomp x ((hasDerivAt_neg x).hasDerivWithinAt) hmaps).neg
    exact this.congr_deriv (by simp [hh''])
  have hDh : ∀ x ∈ Icc (-V) V, |h'' x| ≤ D₂ := fun x hx => hD (-x) (hneg x hx)
  have key := CCMTaylorHalf.taylor_half h h' h'' hdh hdh' hDh (hneg b hb) (hneg a ha)
    (by linarith)
  simp only [hh, hh', neg_neg] at key
  have e : g a - g b - g' b * (a - b) = g a - g b - -g' b * (-a - -b) := by ring
  rw [e]
  calc _ ≤ D₂ * (-a - -b) ^ 2 / 2 := key
    _ = D₂ * (b - a) ^ 2 / 2 := by ring

theorem lowner_entry_sym_half (g g' g'' : ℝ → ℝ) {V D₂ : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) {a b : ℝ} (ha : a ∈ Icc (-V) V) (hb : b ∈ Icc (-V) V)
    (hab : a ≠ b) :
    |(g b - g a) / (Real.pi * (a - b)) + (g' a + g' b) / (2 * Real.pi)|
      ≤ D₂ * |b - a| / (2 * Real.pi) := by
  have hpi := Real.pi_pos
  -- both one-sided expansions, in a common form
  have hL : |g b - g a - g' a * (b - a)| ≤ D₂ * (b - a) ^ 2 / 2 ∧
      |g b - g a - g' b * (b - a)| ≤ D₂ * (b - a) ^ 2 / 2 := by
    rcases le_total a b with h | h
    · refine ⟨CCMTaylorHalf.taylor_half g g' g'' hg hg' hD ha hb h, ?_⟩
      have := taylor_half_right g g' g'' hg hg' hD ha hb h
      rw [show g b - g a - g' b * (b - a) = -(g a - g b - g' b * (a - b)) by ring, abs_neg]
      exact this
    · have h1 := CCMTaylorHalf.taylor_half g g' g'' hg hg' hD hb ha h
      have h2 := taylor_half_right g g' g'' hg hg' hD hb ha h
      constructor
      · rw [show (b - a) ^ 2 = (a - b) ^ 2 by ring]
        exact h2
      · rw [show g b - g a - g' b * (b - a) = -(g a - g b - g' b * (a - b)) by ring, abs_neg,
          show (b - a) ^ 2 = (a - b) ^ 2 by ring]
        exact h1
  obtain ⟨hL1, hL2⟩ := hL
  have hne : b - a ≠ 0 := sub_ne_zero.mpr (Ne.symm hab)
  have hpos : 0 < |b - a| := abs_pos.mpr hne
  have e : (g b - g a) / (Real.pi * (a - b)) + (g' a + g' b) / (2 * Real.pi)
      = -((g b - g a - g' a * (b - a)) + (g b - g a - g' b * (b - a))) / (2 * Real.pi * (b - a)) := by
    field_simp; ring
  rw [e, abs_div, abs_neg, abs_mul, abs_mul, abs_of_pos (by norm_num : (0:ℝ) < 2), abs_of_pos hpi,
    div_le_div_iff₀ (by positivity) (by positivity)]
  have hsum := abs_add_le (g b - g a - g' a * (b - a)) (g b - g a - g' b * (b - a))
  have hsq : (b - a) ^ 2 = |b - a| * |b - a| := by rw [← sq, sq_abs]
  rw [hsq] at hL1 hL2
  nlinarith [abs_nonneg (b - a)]

end CCMSymHalf

#print axioms CCMSymHalf.taylor_half_right
#print axioms CCMSymHalf.lowner_entry_sym_half
