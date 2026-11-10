import ConvexSphereMs
import GrowthAtKGamma

/-!
# S3010 — growth input from the CONVEXITY sphere bound (lane 1)

Analogue of S2961 `growth_at_sharpN`, with Ms from S3008:
`growth_at_convex`: ∃ C ≥ 1, for 1 < σ, σ − 1 ≤ R ≤ 1, σ + R ≤ 2, 0 < R < R' ≤ |t|, |t| ≥ 3:
  (−ζ'/ζ(σ+it)).re ≤ 4·log((C(2+|t|)^{(1−σ+R)/2}(log(2+|t|)+6) + 1)(3 + 1/(σ−1)))/(R/2).
Leading behaviour: (8/R)·((1−σ+R)/2)·log|t| ≤ 4·log|t| (since 1−σ+R < R) — half of S2963's 8·log|t|.
-/

open Complex Metric

namespace GrowthConvex

theorem growth_at_convex : ∃ C : ℝ, 1 ≤ C ∧ ∀ σ t R R' : ℝ, 1 < σ → σ - 1 ≤ R → R ≤ 1 → σ + R ≤ 2 →
    0 < R → R < R' → R' ≤ |t| → 3 ≤ |t| →
    (-logDeriv riemannZeta ((σ : ℂ) + I * (t : ℂ))).re
      ≤ 4 * Real.log ((C * (2 + |t|) ^ ((1 - σ + R) / 2) * (Real.log (2 + |t|) + 6) + 1)
          * (3 + 1 / (σ - 1))) / (R / 2) := by
  obtain ⟨C, hC1, hC⟩ := ConvexSphereMs.convex_sphere_Ms
  refine ⟨C, hC1, fun σ t R R' hσ hR0 hR1 hσR hR hRR' hR't ht => ?_⟩
  set c : ℂ := (σ : ℂ) + I * (t : ℂ) with hc
  have hcre_eq : c.re = σ := by rw [hc]; simp [Complex.add_re, Complex.mul_re]
  have hcre : 1 < c.re := by rw [hcre_eq]; exact hσ
  have hball : (1 : ℂ) ∉ ball c R' := by rw [hc]; exact one_notMem_ball_test hR't
  have hc' : c = (σ : ℂ) + (t : ℂ) * I := by rw [hc]; ring
  have key := hgrowth_riemannZeta_scalar
    (Ms := C * (2 + |t|) ^ ((1 - σ + R) / 2) * (Real.log (2 + |t|) + 6))
    (B := 4 * Real.log ((C * (2 + |t|) ^ ((1 - σ + R) / 2) * (Real.log (2 + |t|) + 6) + 1)
      * (3 + 1 / (c.re - 1))) / (R / 2)) hR hRR' hcre hball
    (fun z hz => by
      rw [hc'] at hz
      exact hC σ t R z hσ hR0 hR1 hσR ht hz) le_rfl
  rw [hcre_eq] at key
  exact key

end GrowthConvex

#print axioms GrowthConvex.growth_at_convex
