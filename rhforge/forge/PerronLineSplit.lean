import PerronPrincipalPart
import PerronRemainderVsegIntegrable
import PolePoleVsegIntegrable

set_option maxHeartbeats 1000000

/-!
# §1624 — the Perron line-integral split: `∫ F = ∫ (pole) + ∫ (remainder)`

Along a vertical segment `t ↦ σ + it` with `σ ≠ 1` (so `σ+it ≠ 1`) staying in the holomorphy region, the Perron
line integral of `F(s) = -ζ'/ζ(s)·x^{s+1}/(s(s+1))` splits as the integral of its explicit pole term
`(x²/2)/(s-1)` plus the integral of the holomorphic remainder `R = dslope Ĥ 1` — by the pointwise split §1614 plus
interval-additivity (§1622 + §1623 integrability).  Shifting `σ` left across the zero-free region moves this split's
pole-term integral through the residue while the remainder integral deforms freely.  EXIT0 / std-axioms.
-/

open Complex MeasureTheory

theorem perron_line_split {x : ℝ} (hx : 0 < x) {σ a b : ℝ} (hσ : σ ≠ 1)
    (hne : ∀ t ∈ Set.uIcc a b,
      riemannZeta ((σ : ℂ) + (t : ℝ) * I) ≠ 0 ∧
        ((σ : ℂ) + (t : ℝ) * I) ≠ 0 ∧ ((σ : ℂ) + (t : ℝ) * I) ≠ -1) :
    (∫ t : ℝ in a..b,
        -logDeriv riemannZeta ((σ : ℂ) + (t : ℝ) * I) *
            (x : ℂ) ^ (((σ : ℂ) + (t : ℝ) * I) + 1) /
          (((σ : ℂ) + (t : ℝ) * I) * (((σ : ℂ) + (t : ℝ) * I) + 1)))
      = (∫ t : ℝ in a..b, (x : ℂ) ^ (2 : ℂ) / 2 / (((σ : ℂ) + (t : ℝ) * I) - 1))
        + (∫ t : ℝ in a..b,
            dslope
              (Function.update
                (fun u => (u - 1) * (-logDeriv riemannZeta u * (x : ℂ) ^ (u + 1) / (u * (u + 1))))
                1 ((x : ℂ) ^ (2 : ℂ) / 2)) 1 ((σ : ℂ) + (t : ℝ) * I)) := by
  rw [← intervalIntegral.integral_add (intervalIntegrable_pole_vseg (x := x) hσ)
    (intervalIntegrable_perron_remainder_vseg hx hne)]
  apply intervalIntegral.integral_congr
  intro t _ht
  have hst1 : (σ : ℂ) + (t : ℝ) * I ≠ 1 := by
    intro h; apply hσ; have := congrArg Complex.re h; simpa using this
  have h14 := perron_dslope_eq (x := x) hst1
  simp only [Pi.add_apply]
  rw [h14]; ring

#print axioms perron_line_split
