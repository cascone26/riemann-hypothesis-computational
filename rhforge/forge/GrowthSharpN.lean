import GrowthAtKGamma
import ZetaSphereSharpN

/-!
# S2961 — growth input from the sharp N-term sphere bound (lane 1, wiring S2959)

`GrowthK.growth_at` (used by every 5-8-4-1 rung S2945–S2957) feeds Borel–Carathéodory with the N = 1 sphere
bound.  Here the same `hgrowth_riemannZeta_scalar` is fed the S2959 sharp bound instead, with all complex norms
replaced by real quantities (`‖c-1‖ ≥ |t|`, `‖c‖ ≤ σ+|t|`):

  `sharpMs σ t R N := (N^{1-a}-1)/(1-a) + 1 + N^{1-a}/(|t|-R) + (σ+|t|+R)·N^{-a}/a`,  `a := σ - R`,

* `growth_at_sharpN` — `(-ζ'/ζ(σ+it)).re ≤ 4·log((sharpMs+1)(3+1/(σ-1)))/(R/2)`  (needs `σ - R < 1`).
* `growth_Tz_sharpN` — the same at `σ + i·kγ`, in S2926's `Tz` form.
With `N ≍ |t|`, `sharpMs ≪ |t|^{1-a}` and `1 - a = 1 - σ + R < R`, so the leading term is `≤ 8·log|t|`
for EVERY radius R (vs 8/R·log|t| from N = 1).  Next brick: the clean scalar form with `N = ⌊|t|⌋`.
-/

open Complex Metric

namespace GrowthSharp

noncomputable def sharpMs (σ t R : ℝ) (N : ℕ) : ℝ :=
  ((N : ℝ) ^ (1 - (σ - R)) - 1) / (1 - (σ - R)) + 1 + (N : ℝ) ^ (1 - (σ - R)) / (|t| - R)
    + (σ + |t| + R) * ((N : ℝ) ^ (-(σ - R)) / (σ - R))

theorem sphere_le_sharpMs {σ t R R' : ℝ} {N : ℕ} (hN : 1 ≤ N) (hσ : 1 < σ) (hR : 0 < R)
    (hR1 : R < 1) (haR : σ - R < 1) (hRR' : R < R') (hR't : R' ≤ |t|) {z : ℂ}
    (hz : z ∈ sphere ((σ : ℂ) + I * (t : ℂ)) R) :
    ‖riemannZeta z‖ ≤ sharpMs σ t R N := by
  set c : ℂ := (σ : ℂ) + I * (t : ℂ) with hc
  have hcre : c.re = σ := by rw [hc]; simp [Complex.add_re, Complex.mul_re]
  have hcR : 0 < c.re - R := by rw [hcre]; linarith
  have ha1 : c.re - R < 1 := by rw [hcre]; exact haR
  have ht1 : |t| ≤ ‖c - 1‖ := by rw [hc]; exact abs_le_norm_center_sub_one
  have hc1 : R < ‖c - 1‖ := by linarith
  have h := norm_riemannZeta_le_sphere_sharpN hN hR hcR ha1 hc1 hz
  rw [hcre] at h
  have hNR : (1 : ℝ) ≤ (N : ℝ) := by exact_mod_cast hN
  have hNpos : (0 : ℝ) < (N : ℝ) := by linarith
  have hcn : ‖c‖ ≤ σ + |t| := by
    calc ‖c‖ ≤ ‖(σ : ℂ)‖ + ‖I * (t : ℂ)‖ := norm_add_le _ _
      _ = σ + |t| := by
        rw [norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Complex.norm_real,
          Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos (by linarith)]
  have htR : 0 < |t| - R := by linarith
  have hB : (N : ℝ) ^ (1 - (σ - R)) / (‖c - 1‖ - R) ≤ (N : ℝ) ^ (1 - (σ - R)) / (|t| - R) :=
    div_le_div_of_nonneg_left (Real.rpow_nonneg (le_of_lt hNpos) _) htR (by linarith)
  have hT : (‖c‖ + R) * ((N : ℝ) ^ (-(σ - R)) / (σ - R))
      ≤ (σ + |t| + R) * ((N : ℝ) ^ (-(σ - R)) / (σ - R)) := by
    have : 0 ≤ (N : ℝ) ^ (-(σ - R)) / (σ - R) :=
      div_nonneg (Real.rpow_nonneg (le_of_lt hNpos) _) (by linarith)
    exact mul_le_mul_of_nonneg_right (by linarith) this
  unfold sharpMs
  linarith

theorem growth_at_sharpN {σ t R R' : ℝ} {N : ℕ} (hN : 1 ≤ N) (hσ : 1 < σ) (hR : 0 < R)
    (hR1 : R < 1) (haR : σ - R < 1) (hRR' : R < R') (hR't : R' ≤ |t|) :
    (-logDeriv riemannZeta ((σ : ℂ) + I * (t : ℂ))).re
      ≤ 4 * Real.log ((sharpMs σ t R N + 1) * (3 + 1 / (σ - 1))) / (R / 2) := by
  set c : ℂ := (σ : ℂ) + I * (t : ℂ) with hc
  have hcre_eq : c.re = σ := by rw [hc]; simp [Complex.add_re, Complex.mul_re]
  have hcre : 1 < c.re := by rw [hcre_eq]; exact hσ
  have hball : (1 : ℂ) ∉ ball c R' := by rw [hc]; exact one_notMem_ball_test hR't
  have key := hgrowth_riemannZeta_scalar (B := 4 * Real.log ((sharpMs σ t R N + 1)
      * (3 + 1 / (c.re - 1))) / (R / 2)) hR hRR' hcre hball
    (fun z hz => sphere_le_sharpMs hN hσ hR hR1 haR hRR' hR't hz) le_rfl
  rw [hcre_eq] at key
  exact key

/-- Growth at the k-th point σ + i·kγ, in the `Tz` form of S2926. -/
theorem growth_Tz_sharpN {σ γ R R' : ℝ} {N : ℕ} (k : ℕ) (hN : 1 ≤ N) (hσ : 1 < σ)
    (hR : 0 < R) (hR1 : R < 1) (haR : σ - R < 1) (hRR' : R < R')
    (hR't : R' ≤ |(k : ℝ) * γ|) :
    TrigPolyRegion.Tz σ γ k
      ≤ 4 * Real.log ((sharpMs σ ((k : ℝ) * γ) R N + 1) * (3 + 1 / (σ - 1))) / (R / 2) :=
  growth_at_sharpN hN hσ hR hR1 haR hRR' hR't

end GrowthSharp

#print axioms GrowthSharp.sphere_le_sharpMs
#print axioms GrowthSharp.growth_at_sharpN
#print axioms GrowthSharp.growth_Tz_sharpN
