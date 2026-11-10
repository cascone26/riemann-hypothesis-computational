import Mathlib

/-!
# S3036 — the θ-integrals of the Poisson–Landau bound (lane 1, brick P4)

∫₀^{2π} max(−cos θ, 0)² dθ = π/2,  ∫₀^{2π} max(cos θ, 0) dθ = 2,  ∫₀^{2π} max(−cos θ, 0) dθ = 2.
With the Schwarz formula (S3034) Re(−H'(0)) = (2/R)·(2π)⁻¹∫(−cos θ)·log|G(Re^{iθ})| dθ and the left-half convexity bound
log|ζ| ≤ (R·max(−cos θ,0)/2)·log t + O(loglog t), the first integral gives (2/R)(2π)⁻¹(R/2)(π/2)·log t = (1/4)·log t;
the other two carry the O(log log t) terms.
-/

open Real intervalIntegral MeasureTheory

namespace CosHalf

lemma cos_nonneg_Icc1 {θ : ℝ} (h : θ ∈ Set.Icc 0 (π / 2)) : 0 ≤ cos θ :=
  cos_nonneg_of_mem_Icc ⟨by linarith [h.1, pi_pos], h.2⟩
lemma cos_nonpos_Icc2 {θ : ℝ} (h : θ ∈ Set.Icc (π / 2) (3 * π / 2)) : cos θ ≤ 0 :=
  cos_nonpos_of_pi_div_two_le_of_le h.1 (by linarith [h.2])
lemma cos_nonneg_Icc3 {θ : ℝ} (h : θ ∈ Set.Icc (3 * π / 2) (2 * π)) : 0 ≤ cos θ := by
  have : cos θ = cos (θ - 2 * π) := by rw [cos_sub_two_pi]
  rw [this]; exact cos_nonneg_of_mem_Icc ⟨by linarith [h.1], by linarith [h.2, pi_pos]⟩

theorem integral_negcos_pos_sq : ∫ θ in (0 : ℝ)..(2 * π), (max (-cos θ) 0) ^ 2 = π / 2 := by
  have hc : Continuous fun θ : ℝ => (max (-cos θ) 0) ^ 2 := by fun_prop
  have i1 : IntervalIntegrable (fun θ : ℝ => (max (-cos θ) 0) ^ 2) volume 0 (π / 2) := hc.intervalIntegrable _ _
  have i2 : IntervalIntegrable (fun θ : ℝ => (max (-cos θ) 0) ^ 2) volume (π / 2) (3 * π / 2) := hc.intervalIntegrable _ _
  have i3 : IntervalIntegrable (fun θ : ℝ => (max (-cos θ) 0) ^ 2) volume (3 * π / 2) (2 * π) := hc.intervalIntegrable _ _
  rw [← integral_add_adjacent_intervals (i1.trans i2) i3, ← integral_add_adjacent_intervals i1 i2]
  have e1 : ∫ θ in (0 : ℝ)..(π / 2), (max (-cos θ) 0) ^ 2 = 0 := by
    rw [integral_congr (g := fun _ => (0 : ℝ))]; · simp
    intro θ hθ; rw [Set.uIcc_of_le (by linarith [pi_pos])] at hθ
    simp [max_eq_right (by linarith [cos_nonneg_Icc1 hθ] : -cos θ ≤ 0)]
  have e3 : ∫ θ in (3 * π / 2)..(2 * π), (max (-cos θ) 0) ^ 2 = 0 := by
    rw [integral_congr (g := fun _ => (0 : ℝ))]; · simp
    intro θ hθ; rw [Set.uIcc_of_le (by linarith [pi_pos])] at hθ
    simp [max_eq_right (by linarith [cos_nonneg_Icc3 hθ] : -cos θ ≤ 0)]
  have e2 : ∫ θ in (π / 2)..(3 * π / 2), (max (-cos θ) 0) ^ 2 = π / 2 := by
    rw [integral_congr (g := fun θ => cos θ ^ 2)]
    · rw [integral_cos_sq]
      have h1 : sin (3 * π / 2) = -1 := by
        rw [show 3 * π / 2 = π / 2 + π by ring, sin_add_pi, sin_pi_div_two]
      have h2 : cos (3 * π / 2) = 0 := by
        rw [show 3 * π / 2 = π / 2 + π by ring, cos_add_pi, cos_pi_div_two, neg_zero]
      rw [h1, h2, sin_pi_div_two, cos_pi_div_two]; ring
    intro θ hθ; rw [Set.uIcc_of_le (by linarith [pi_pos])] at hθ
    simp only
    rw [max_eq_left (by linarith [cos_nonpos_Icc2 hθ] : (0 : ℝ) ≤ -cos θ)]; ring
  rw [e1, e2, e3]; ring

theorem integral_cos_pos : ∫ θ in (0 : ℝ)..(2 * π), max (cos θ) 0 = 2 := by
  have hc : Continuous fun θ : ℝ => max (cos θ) 0 := by fun_prop
  have i1 : IntervalIntegrable (fun θ : ℝ => max (cos θ) 0) volume 0 (π / 2) := hc.intervalIntegrable _ _
  have i2 : IntervalIntegrable (fun θ : ℝ => max (cos θ) 0) volume (π / 2) (3 * π / 2) := hc.intervalIntegrable _ _
  have i3 : IntervalIntegrable (fun θ : ℝ => max (cos θ) 0) volume (3 * π / 2) (2 * π) := hc.intervalIntegrable _ _
  rw [← integral_add_adjacent_intervals (i1.trans i2) i3, ← integral_add_adjacent_intervals i1 i2]
  have e1 : ∫ θ in (0 : ℝ)..(π / 2), max (cos θ) 0 = 1 := by
    rw [integral_congr (g := fun θ => cos θ)]
    · rw [integral_cos, sin_pi_div_two, sin_zero]; ring
    intro θ hθ; rw [Set.uIcc_of_le (by linarith [pi_pos])] at hθ
    exact max_eq_left (cos_nonneg_Icc1 hθ)
  have e2 : ∫ θ in (π / 2)..(3 * π / 2), max (cos θ) 0 = 0 := by
    rw [integral_congr (g := fun _ => (0 : ℝ))]; · simp
    intro θ hθ; rw [Set.uIcc_of_le (by linarith [pi_pos])] at hθ
    exact max_eq_right (cos_nonpos_Icc2 hθ)
  have e3 : ∫ θ in (3 * π / 2)..(2 * π), max (cos θ) 0 = 1 := by
    rw [integral_congr (g := fun θ => cos θ)]
    · rw [integral_cos, sin_two_pi]
      rw [show 3 * π / 2 = π / 2 + π by ring, sin_add_pi, sin_pi_div_two]; ring
    intro θ hθ; rw [Set.uIcc_of_le (by linarith [pi_pos])] at hθ
    exact max_eq_left (cos_nonneg_Icc3 hθ)
  rw [e1, e2, e3]; ring

theorem integral_negcos_pos : ∫ θ in (0 : ℝ)..(2 * π), max (-cos θ) 0 = 2 := by
  have hc : Continuous fun θ : ℝ => max (-cos θ) 0 := by fun_prop
  have i1 : IntervalIntegrable (fun θ : ℝ => max (-cos θ) 0) volume 0 (π / 2) := hc.intervalIntegrable _ _
  have i2 : IntervalIntegrable (fun θ : ℝ => max (-cos θ) 0) volume (π / 2) (3 * π / 2) := hc.intervalIntegrable _ _
  have i3 : IntervalIntegrable (fun θ : ℝ => max (-cos θ) 0) volume (3 * π / 2) (2 * π) := hc.intervalIntegrable _ _
  rw [← integral_add_adjacent_intervals (i1.trans i2) i3, ← integral_add_adjacent_intervals i1 i2]
  have e1 : ∫ θ in (0 : ℝ)..(π / 2), max (-cos θ) 0 = 0 := by
    rw [integral_congr (g := fun _ => (0 : ℝ))]; · simp
    intro θ hθ; rw [Set.uIcc_of_le (by linarith [pi_pos])] at hθ
    exact max_eq_right (by linarith [cos_nonneg_Icc1 hθ])
  have e2 : ∫ θ in (π / 2)..(3 * π / 2), max (-cos θ) 0 = 2 := by
    rw [integral_congr (g := fun θ => -cos θ)]
    · rw [intervalIntegral.integral_neg, integral_cos]
      rw [show 3 * π / 2 = π / 2 + π by ring, sin_add_pi, sin_pi_div_two]; ring
    intro θ hθ; rw [Set.uIcc_of_le (by linarith [pi_pos])] at hθ
    exact max_eq_left (by linarith [cos_nonpos_Icc2 hθ])
  have e3 : ∫ θ in (3 * π / 2)..(2 * π), max (-cos θ) 0 = 0 := by
    rw [integral_congr (g := fun _ => (0 : ℝ))]; · simp
    intro θ hθ; rw [Set.uIcc_of_le (by linarith [pi_pos])] at hθ
    exact max_eq_right (by linarith [cos_nonneg_Icc3 hθ])
  rw [e1, e2, e3]; ring

end CosHalf

#print axioms CosHalf.integral_negcos_pos_sq
#print axioms CosHalf.integral_cos_pos
#print axioms CosHalf.integral_negcos_pos
