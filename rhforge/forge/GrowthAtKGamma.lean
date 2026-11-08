/-
GrowthAtKGamma.lean (HP RH center, LEDGER S2927b→S2928, 2026-11-08) — lane 1.
Point-generic growth input for the P-generic region (S2926): at c = σ + it with σ > 1, R < 1, R < R' ≤ |t|,
  Re(−ζ'/ζ)(c) ≤ 4·log((Ms+1)(3 + 1/(σ−1)))/(R/2),   Ms = 2 + 1/(‖c−1‖−R) + (‖c‖+R)/(c.re−R).
Composes §1555 (`hgrowth_riemannZeta_scalar`), §1447 (`norm_riemannZeta_le_sphere_N1`), §1371 (`one_notMem_ball_test`),
`abs_le_norm_center_sub_one`. Instantiated at t = kγ (`growth_Tz`) — discharges S2926's `hgrowth` for every k ≥ 2,
in particular the 3iγ term that 5-8-4-1 needs and the 3-4-1 tower never had.
-/
import HgrowthRiemannZetaScalar
import ZetaSphereN1
import PoleAvoidGeom
import AbsLeNormCenterSubOne
import TrigPolyZetaRegion

open Complex Metric

namespace GrowthK

theorem growth_at {σ t R R' : ℝ} (hσ : 1 < σ) (hR : 0 < R) (hR1 : R < 1) (hRR' : R < R')
    (hR't : R' ≤ |t|) :
    (-logDeriv riemannZeta ((σ : ℂ) + I * (t : ℂ))).re
      ≤ 4 * Real.log ((2 + 1 / (‖(σ : ℂ) + I * (t : ℂ) - 1‖ - R)
          + (‖(σ : ℂ) + I * (t : ℂ)‖ + R) / (((σ : ℂ) + I * (t : ℂ)).re - R) + 1)
          * (3 + 1 / (((σ : ℂ) + I * (t : ℂ)).re - 1))) / (R / 2) := by
  set c : ℂ := (σ : ℂ) + I * (t : ℂ) with hc
  have hcre_eq : c.re = σ := by
    rw [hc]; simp [Complex.add_re, Complex.mul_re]
  have hcre : 1 < c.re := by rw [hcre_eq]; exact hσ
  have hcR : 0 < c.re - R := by rw [hcre_eq]; linarith
  have hc1 : R < ‖c - 1‖ :=
    lt_of_lt_of_le (lt_of_lt_of_le hRR' hR't) (by rw [hc]; exact abs_le_norm_center_sub_one)
  have hball : (1 : ℂ) ∉ ball c R' := by rw [hc]; exact one_notMem_ball_test hR't
  exact hgrowth_riemannZeta_scalar hR hRR' hcre hball
    (fun z hz => norm_riemannZeta_le_sphere_N1 hR hcR hc1 hz) le_rfl

/-- Growth at the k-th point σ + i·kγ, in the `Tz` form of S2926. -/
theorem growth_Tz {σ γ R R' : ℝ} (k : ℕ) (hσ : 1 < σ) (hR : 0 < R) (hR1 : R < 1)
    (hRR' : R < R') (hR't : R' ≤ |(k : ℝ) * γ|) :
    TrigPolyRegion.Tz σ γ k
      ≤ 4 * Real.log ((2 + 1 / (‖(σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ) - 1‖ - R)
          + (‖(σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ)‖ + R)
            / (((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ)).re - R) + 1)
          * (3 + 1 / (((σ : ℂ) + I * (((k : ℝ) * γ : ℝ) : ℂ)).re - 1))) / (R / 2) :=
  growth_at hσ hR hR1 hRR' hR't

end GrowthK

#print axioms GrowthK.growth_at
#print axioms GrowthK.growth_Tz
