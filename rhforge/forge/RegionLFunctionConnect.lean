import RegionZpow
import DLVPCharPositivityLFunction
import DivisorLFunctionNonneg
import OneLeDivisorLFunction
import DivisorFilterReLtLFunction
import LFunctionZeroMemFilter

/-!
# §1498 — connector: the L(·,χ) quantitative zero-free region from the four analytic bounds

L-analogue of ζ §1409, but driven by the generic region master §1336 directly (its `hdlvp` supplied by §1497, the
`LFunction`-form DLVP positivity — the three `ζ³·L(χ)⁴·L(χ²)` terms map onto §1336's `T₀`/`F`/`T₂` slots).  Fixes the
concrete zero-set `S = support(ball c R').filter (·∈closedBall c Rh)`, `c = σK+iγ`, `m = divisor`, and discharges
the geometric hypotheses `hρ₀` (§1496), `hleft` (§1494), `hm` (§1479), `hm₀` (§1492) from the zero data.  Reduces
the region to the four analytic bounds — `hdecomp`/`hgB` (§1493), `hpole` (ζ §1405), `hgrowth` (§1491 for χ²) —
plus the kept-zero geometry.  Conclusion: `β ≤ 1 - 1/(14K)`.  For `χ ≠ 1`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter Metric Classical

theorem region_LFunction_connect {N : ℕ} [NeZero N] {χ : DirichletCharacter ℂ N} (hχ : χ ≠ 1)
    {γ β A B C : ℝ} {R' Rh : ℝ} {g : ℂ}
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hLρ₀ : LFunction χ ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hz₀ : (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hz₀ne : LFunction χ (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) ≠ 0)
    (hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) Rh)
    (hgB : -B ≤ g.re)
    (hdecomp : logDeriv (LFunction χ) (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))
      = (∑ ρ ∈ (divisor_LFunction_ball_support_finite hχ
            (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R').toFinset.filter
          (fun u => u ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) Rh),
          (MeromorphicOn.divisor (LFunction χ)
            (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R') ρ : ℂ)
          / (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ) - ρ)) + g)
    (hpole : (-logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ))).re
      ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : (-logDeriv (LFunction (χ ^ 2))
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))).re ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  refine region_of_zpow_decomp_and_bounds
    (m := MeromorphicOn.divisor (LFunction χ)
      (ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R'))
    hK hA hB hC hβ ?_ ?_ ?_ ?_ hgB hdecomp ?_ hpole hgrowth
  · exact LFunction_zero_mem_filter hχ hρ₀ball hLρ₀ hz₀ hz₀ne _ hp
  · exact fun ρ hρ => divisor_filter_support_re_lt_LFunction hχ hσK _ hρ
  · exact fun ρ _ => divisor_LFunction_nonneg_all hχ ρ
  · exact one_le_divisor_LFunction_at_zero hχ hρ₀ball hLρ₀ hz₀ hz₀ne
  · exact dlvp_char_positivity_LFunction χ hσK

#print axioms region_LFunction_connect
