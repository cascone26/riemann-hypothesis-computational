import GrowthConvex

/-!
# S3012 — clean convexity growth: leading coefficient 4 (lane 1)

`growth_convex_clean`: ∃ C ≥ 1, for 1 < σ, σ − 1 ≤ R ≤ 1, σ + R ≤ 2, 0 < R < R' ≤ |t|, |t| ≥ 3:
  (−ζ'/ζ(σ+it)).re ≤ 4·log(2+|t|) + (8/R)·(log(C+1) + log(log(2+|t|)+6) + log(3 + 1/(σ−1))).
(From S3010: Ms + 1 ≤ (C+1)X^p·Lg with X = 2+|t|, p = (1−σ+R)/2 ≤ R/2, Lg = log X + 6; (8/R)·p·log X ≤ 4·log X.)
Compare S2963: 8·log|t| + (8/R)(…). Leading coefficient HALVED — input for the quartic rung with κ → 4.
-/

open Complex

namespace GrowthConvexClean

theorem growth_convex_clean : ∃ C : ℝ, 1 ≤ C ∧ ∀ σ t R R' : ℝ, 1 < σ → σ - 1 ≤ R → R ≤ 1 →
    σ + R ≤ 2 → 0 < R → R < R' → R' ≤ |t| → 3 ≤ |t| →
    (-logDeriv riemannZeta ((σ : ℂ) + I * (t : ℂ))).re
      ≤ 4 * Real.log (2 + |t|) + 8 / R * (Real.log (C + 1) + Real.log (Real.log (2 + |t|) + 6)
          + Real.log (3 + 1 / (σ - 1))) := by
  obtain ⟨C, hC1, hC⟩ := GrowthConvex.growth_at_convex
  refine ⟨C, hC1, fun σ t R R' hσ hR0 hR1 hσR hR hRR' hR't ht => ?_⟩
  refine le_trans (hC σ t R R' hσ hR0 hR1 hσR hR hRR' hR't ht) ?_
  set X := 2 + |t| with hX
  set p := (1 - σ + R) / 2 with hp
  set Lg := Real.log X + 6 with hLg
  set D := 3 + 1 / (σ - 1) with hD
  have hX1 : 1 ≤ X := by rw [hX]; linarith [abs_nonneg t]
  have hX0 : 0 < X := by linarith
  have hp0 : 0 ≤ p := by rw [hp]; linarith
  have hpR : p ≤ R / 2 := by rw [hp]; linarith
  have hlogX : 0 ≤ Real.log X := Real.log_nonneg hX1
  have hLg6 : 6 ≤ Lg := by rw [hLg]; linarith
  have hD0 : 0 < D := by
    rw [hD]; have : 0 < 1 / (σ - 1) := by apply div_pos one_pos; linarith
    linarith
  have hXp : 1 ≤ X ^ p := Real.one_le_rpow hX1 hp0
  have hXp0 : 0 < X ^ p := by linarith
  have hM1 : C * X ^ p * Lg + 1 ≤ (C + 1) * (X ^ p * Lg) := by nlinarith
  have hMpos : 0 < C * X ^ p * Lg + 1 := by positivity
  have hlog : Real.log ((C * X ^ p * Lg + 1) * D)
      ≤ Real.log (C + 1) + p * Real.log X + Real.log Lg + Real.log D := by
    calc Real.log ((C * X ^ p * Lg + 1) * D)
        ≤ Real.log (((C + 1) * (X ^ p * Lg)) * D) :=
          Real.log_le_log (mul_pos hMpos hD0) (mul_le_mul_of_nonneg_right hM1 hD0.le)
      _ = Real.log (C + 1) + p * Real.log X + Real.log Lg + Real.log D := by
          rw [Real.log_mul (by positivity) hD0.ne', Real.log_mul (by positivity) (by positivity),
            Real.log_mul hXp0.ne' (by positivity), Real.log_rpow hX0]; ring
  have hsplit : 4 * Real.log ((C * X ^ p * Lg + 1) * D) / (R / 2)
      = 8 / R * Real.log ((C * X ^ p * Lg + 1) * D) := by
    rw [div_div_eq_mul_div]; ring
  rw [hsplit]
  have h8R : 0 ≤ 8 / R := by positivity
  have step1 := mul_le_mul_of_nonneg_left hlog h8R
  have step2 : 8 / R * (p * Real.log X) ≤ 4 * Real.log X := by
    have : 8 / R * (p * Real.log X) = Real.log X * (8 * p / R) := by ring
    rw [this]
    have hq : 8 * p / R ≤ 4 := by rw [div_le_iff₀ hR]; linarith
    nlinarith
  have hexp : 8 / R * (Real.log (C + 1) + p * Real.log X + Real.log Lg + Real.log D)
      = 8 / R * (p * Real.log X) + 8 / R * (Real.log (C + 1) + Real.log Lg + Real.log D) := by ring
  linarith

end GrowthConvexClean

#print axioms GrowthConvexClean.growth_convex_clean
