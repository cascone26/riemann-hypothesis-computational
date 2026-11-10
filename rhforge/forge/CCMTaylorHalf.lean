import CCMSchurTest

/-!
# S2996 — Taylor's ½: |g b − g a − g'(a)(b−a)| ≤ D₂(b−a)²/2 (lane 2: removes S2987's factor-2 slack)

* `taylor_half` — a ≤ b in [−V,V], |g''| ≤ D₂ ⇒ |g b − g a − g'(a)(b − a)| ≤ D₂(b − a)²/2.
  Proof: φ±(x) = D₂(x−a)²/2 ∓ (g x − g a − g'(a)(x−a)) has φ±(a) = 0 and φ±' = D₂(x−a) ∓ (g'(x) − g'(a)) ≥ 0
  (S2987 `deriv_lipschitz`), so MVT (exists_hasDerivAt_eq_slope) gives φ±(b) ≥ 0.
-/

open Set

namespace CCMTaylorHalf

theorem half_aux (g g' g'' : ℝ → ℝ) {V D₂ : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) {a b : ℝ} (ha : a ∈ Icc (-V) V) (hb : b ∈ Icc (-V) V)
    (hab : a < b) (s : ℝ) (hs : s = 1 ∨ s = -1) :
    0 ≤ D₂ * (b - a) ^ 2 / 2 - s * (g b - g a - g' a * (b - a)) := by
  set φ : ℝ → ℝ := fun x => D₂ * (x - a) ^ 2 / 2 - s * (g x - g a - g' a * (x - a)) with hφ
  set φ' : ℝ → ℝ := fun x => D₂ * (x - a) - s * (g' x - g' a) with hφ'
  have hsub : Icc a b ⊆ Icc (-V) V := Icc_subset_Icc ha.1 hb.2
  have hcont : ContinuousOn φ (Icc a b) := by
    have hgc : ContinuousOn g (Icc a b) :=
      fun x hx => ((hg x (hsub hx)).continuousWithinAt).mono hsub
    rw [hφ]
    exact (continuousOn_const.mul ((continuousOn_id.sub continuousOn_const).pow 2)).div_const _
      |>.sub (continuousOn_const.mul ((hgc.sub continuousOn_const).sub
        (continuousOn_const.mul (continuousOn_id.sub continuousOn_const))))
  have hderiv : ∀ x ∈ Ioo a b, HasDerivAt φ (φ' x) x := by
    intro x hx
    have hxV : x ∈ Ioo (-V) V := ⟨lt_of_le_of_lt ha.1 hx.1, lt_of_lt_of_le hx.2 hb.2⟩
    have hnhds : Icc (-V) V ∈ nhds x := Icc_mem_nhds hxV.1 hxV.2
    have hgx : HasDerivAt g (g' x) x := (hg x (Ioo_subset_Icc_self hxV)).hasDerivAt hnhds
    have h1 : HasDerivAt (fun y => D₂ * (y - a) ^ 2 / 2) (D₂ * (x - a)) x := by
      have := ((((hasDerivAt_id' x).sub_const a).pow 2).const_mul D₂).div_const 2
      exact this.congr_deriv (by ring)
    have h2 : HasDerivAt (fun y => s * (g y - g a - g' a * (y - a))) (s * (g' x - g' a)) x := by
      have := ((hgx.sub_const (g a)).sub (((hasDerivAt_id' x).sub_const a).const_mul (g' a))).const_mul s
      exact this.congr_deriv (by ring)
    exact h1.sub h2
  obtain ⟨c, hc, hslope⟩ := exists_hasDerivAt_eq_slope φ φ' hab hcont hderiv
  have hcV : c ∈ Icc (-V) V := hsub (Ioo_subset_Icc_self hc)
  have hlip := CCMLownerDeriv2.deriv_lipschitz g' g'' hg' hD ha hcV
  have hca : 0 < c - a := by linarith [hc.1]
  rw [abs_of_pos hca] at hlip
  have hφ'c : 0 ≤ φ' c := by
    rw [hφ']; simp only
    rcases hs with rfl | rfl
    · have := (abs_le.mp hlip).2; linarith
    · have := (abs_le.mp hlip).1; linarith
  have hφa : φ a = 0 := by rw [hφ]; simp
  have hba : 0 < b - a := by linarith
  have : 0 ≤ φ b := by
    have h := hslope
    rw [hφa, sub_zero, eq_div_iff hba.ne'] at h
    nlinarith
  simpa [hφ] using this

theorem taylor_half (g g' g'' : ℝ → ℝ) {V D₂ : ℝ}
    (hg : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g (g' x) (Icc (-V) V) x)
    (hg' : ∀ x ∈ Icc (-V) V, HasDerivWithinAt g' (g'' x) (Icc (-V) V) x)
    (hD : ∀ x ∈ Icc (-V) V, |g'' x| ≤ D₂) {a b : ℝ} (ha : a ∈ Icc (-V) V) (hb : b ∈ Icc (-V) V)
    (hab : a ≤ b) :
    |g b - g a - g' a * (b - a)| ≤ D₂ * (b - a) ^ 2 / 2 := by
  rcases eq_or_lt_of_le hab with h | h
  · subst h; simp
  have hp := half_aux g g' g'' hg hg' hD ha hb h 1 (Or.inl rfl)
  have hm := half_aux g g' g'' hg hg' hD ha hb h (-1) (Or.inr rfl)
  rw [abs_le]; constructor <;> linarith

end CCMTaylorHalf

#print axioms CCMTaylorHalf.half_aux
#print axioms CCMTaylorHalf.taylor_half
