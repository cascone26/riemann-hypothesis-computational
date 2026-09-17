import RegionFromGrowthDisks
import ZetaBoundSphereCenter

/-!
# §1431 ★★ — quantitative ζ region, BOTH sphere existentials discharged

Discharges §1430's `2iγ`-disk sphere trio (`hζsphere₂`, `hsl₂`, `hMζ₂`) via §1398 at centre `σK+2iγ`
(`γ`-coordinate `2γ`), exactly as §1426 did for the `iγ` disk.  The region `β ≤ 1-1/(14K)` for real ζ now rests on
**purely scalar/geometry data on both disks** — no `Mζ` existential and no `g₀` anywhere: the two explicit BC
scalars (`hBexplicit`, `hBexplicit₂`), the pole value bound `hpole`, the kept-zero geometry, the two growth-geometry
constraints, and the constant inequalities.  Only the scalar constant-choice fixed-point remains.  EXIT0 / std.
-/

open Complex Metric

theorem region_riemannZeta_from_two_growth {γ β A B C R R' R'' R₂ R'₂ R''₂ : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hball : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hρ₀ball : ((β : ℂ) + I * (γ : ℝ))
      ∈ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R')
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hp : ((β : ℂ) + I * (γ : ℝ))
      ∈ closedBall (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) (R / 2))
    (hR''0 : 0 < R'') (hR''1 : R'' < 1) (hR''t : R'' < |γ|)
    (hincl : (1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 + R < R'')
    (hBexplicit : ∀ Ms : ℝ, (∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) R,
        ‖riemannZeta z‖ ≤ Ms) →
      4 * Real.log ((Ms + 1) / ‖riemannZeta
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))‖) / (R / 2) ≤ B)
    (hpole : (-logDeriv riemannZeta (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ))).re
      ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hR₂ : 0 < R₂) (hRR'₂ : R₂ < R'₂)
    (hball₂ : (1 : ℂ) ∉ ball (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)) R'₂)
    (hR''₂0 : 0 < R''₂) (hR''₂1 : R''₂ < 1) (hR''₂t : R''₂ < |2 * γ|)
    (hincl₂ : (1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 + R₂ < R''₂)
    (hBexplicit₂ : ∀ Ms₂ : ℝ, (∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
        + I * ((2 * γ : ℝ) : ℂ)) R₂, ‖riemannZeta z‖ ≤ Ms₂) →
      4 * Real.log ((Ms₂ + 1) / ‖riemannZeta
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))‖) / (R₂ / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  obtain ⟨Ms₂, hMs₂⟩ :=
    exists_norm_riemannZeta_le_sphere (σ := 1 + 1 / (2 * (3 * A + 4 * B + C))) (γ := 2 * γ) (R := R₂)
      hσK hR''₂0 hR''₂1 hR''₂t hincl₂
  have hs0mem₂ : (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)) + (R₂ : ℝ)
      ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)) R₂ := by
    simp only [mem_sphere_iff_norm, add_sub_cancel_left, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos hR₂]
  have hMs₂0 : 0 ≤ Ms₂ := (norm_nonneg _).trans (hMs₂ _ hs0mem₂)
  exact region_riemannZeta_from_growth_disks hR hRR' hK hA hB hC hβ hball hρ₀ball hζρ₀ hp
    hR''0 hR''1 hR''t hincl hBexplicit hpole hR₂ hRR'₂ hball₂ hMs₂ (by linarith) (by linarith)
    (hBexplicit₂ Ms₂ hMs₂)

#print axioms region_riemannZeta_from_two_growth
