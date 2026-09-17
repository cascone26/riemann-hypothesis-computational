import PerronRemainderDiffAt
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic

set_option maxHeartbeats 1000000

/-!
# §1622 — vertical-segment interval-integrability of the remainder `R = dslope Ĥ 1`

Along any vertical segment `t ↦ σ + it` (`t ∈ [a,b]`) that stays in `R`'s holomorphy region (no `ζ`-zero, `≠ 0,-1`),
`t ↦ R(σ+it)` is continuous (§1621 gives `DifferentiableAt`, hence `ContinuousAt`, composed with the continuous
line map), so it is interval-integrable.  This is the concrete integrability input for splitting the Perron line
integral `∫ F = ∫ pole + ∫ R` and for shifting the contour left across the zero-free region.  EXIT0 / std-axioms.
-/

open Complex MeasureTheory

theorem intervalIntegrable_perron_remainder_vseg {x : ℝ} (hx : 0 < x) {σ a b : ℝ}
    (hne : ∀ t ∈ Set.uIcc a b,
      riemannZeta ((σ : ℂ) + (t : ℝ) * I) ≠ 0 ∧
        ((σ : ℂ) + (t : ℝ) * I) ≠ 0 ∧ ((σ : ℂ) + (t : ℝ) * I) ≠ -1) :
    IntervalIntegrable
      (fun t : ℝ =>
        dslope
          (Function.update
            (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
            1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((σ : ℂ) + (t : ℝ) * I))
      volume a b := by
  apply ContinuousOn.intervalIntegrable
  intro t ht
  have hL : ContinuousAt (fun t : ℝ => (σ : ℂ) + (t : ℝ) * I) t := by fun_prop
  obtain ⟨hζ, h0, hm1⟩ := hne t ht
  have hR := differentiableAt_perron_remainder hx (z := (σ : ℂ) + (t : ℝ) * I)
    (fun _ => hζ) h0 hm1
  have hcomp : ContinuousAt
      (fun t : ℝ =>
        dslope
          (Function.update
            (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
            1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((σ : ℂ) + (t : ℝ) * I)) t :=
    ContinuousAt.comp hR.continuousAt hL
  exact hcomp.continuousWithinAt

#print axioms intervalIntegrable_perron_remainder_vseg
