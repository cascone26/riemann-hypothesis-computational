import RegionRiemannZetaWireup
import DivisorNonnegAll
import OneLeDivisorZeta
import DivisorFilterReLt
import ZetaZeroMemFilter

/-!
# §1409 — connector: discharge §1406's geometric hypotheses, reduce to the analytic bounds

Fixes the concrete zero-set `S = support(ball c R').filter p` (`p = (· ∈ closedBall c Rh)`, `c = σK+iγ`) and
`m = divisor` and discharges §1406's `hm` (§1387), `hm₀` (§1351), `hleft` (§1408), `hρ₀` (§1407) from the zero
data.  Reduces the quantitative region to the four analytic bounds `hdecomp` (§1400), `hgB` (§1402), `hpole`
(§1405), `hgrowth` (§1404) plus the geometry of the kept zero.  Conclusion: `β ≤ 1 - 1/(14K)`.  EXIT0 / std.
-/

open Complex Metric Classical

theorem region_riemannZeta_connect {γ β A B C : ℝ} {R' Rh : ℝ} {g : ℂ}
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hball : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hz₀ : (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hz₀ne : riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) ≠ 0)
    (hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) Rh)
    (hgB : -B ≤ g.re)
    (hdecomp : logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))
      = (∑ ρ ∈ (divisor_riemannZeta_ball_support_finite
            (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R').toFinset.filter
          (fun u => u ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) Rh),
          (MeromorphicOn.divisor riemannZeta
            (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R') ρ : ℂ)
          / (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ) - ρ)) + g)
    (hpole : (-logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ))).re
      ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : (-logDeriv riemannZeta
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))).re ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  refine region_riemannZeta_wireup (m := MeromorphicOn.divisor riemannZeta
    (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R'))
    hK hA hB hC hβ ?_ ?_ ?_ ?_ hgB hdecomp hpole hgrowth
  · exact zeta_zero_mem_filter hball hρ₀ball hζρ₀ hz₀ hz₀ne _ hp
  · exact fun ρ hρ => divisor_filter_support_re_lt hσK _ hρ
  · exact fun ρ _ => divisor_riemannZeta_nonneg_all hball ρ
  · exact one_le_divisor_riemannZeta_at_zero hball hρ₀ball hζρ₀ hz₀ hz₀ne

#print axioms region_riemannZeta_connect
