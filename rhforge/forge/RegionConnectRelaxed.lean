import RegionRiemannZetaConnect

/-!
# §1411 — relaxation decoupler: `B` a FREE parameter above the (circular-looking) BC bound

§1409 wants the Borel–Carathéodory slot as `-B ≤ g.re`, but the concrete §1402 bound `BOUND` is evaluated at
`c = σK+iγ = 1+1/(2(3A+4B+C))+iγ`, which itself contains `B` — a genuine fixed-point circularity that blocks
setting `B := BOUND` directly.  This wrapper breaks it: take an **abstract** `bound : ℝ` with `-bound ≤ g.re`
(supplied by §1410 with `bound := BOUND`) and the scalar relaxation `bound ≤ B`, and discharge §1409's `hgB` via
`-B ≤ -bound ≤ g.re`.  Now `B` is a free parameter constrained only by the scalar inequality `BOUND ≤ B` — exactly
the shape the constant-choice fixed-point needs.  EXIT0 / std-axioms.
-/

open Complex Metric Classical

theorem region_riemannZeta_connect_relaxed {γ β A B C : ℝ} {R' Rh : ℝ} {g : ℂ} {bound : ℝ}
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
    (hgBex : -bound ≤ g.re) (hBub : bound ≤ B)
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
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) :=
  region_riemannZeta_connect hK hA hB hC hβ hball hρ₀ball hζρ₀ hz₀ hz₀ne hp
    (le_trans (neg_le_neg hBub) hgBex) hdecomp hpole hgrowth

#print axioms region_riemannZeta_connect_relaxed
