import RegionExplicitScalar
import HfarHolds

/-!
# §1454 — quantitative ζ region with the geometry fixed (`R = R₂ = 1/4`, `R' = |γ|`, `R'₂ = |2γ|`)

Specializes §1453 to the concrete DLVP geometry: `R = R₂ = 1/4`, `R' = |γ|`, `R'₂ = |2γ|`.  `hfar` is discharged by
§1434 (`8 ≤ 7·K·(1/4)`, i.e. `K ≥ 32/7`), and the radius/ordering hypotheses reduce to `1 < |γ|`.  The region
`β ≤ 1-1/(14K)` for real ζ now follows from just: `ζ(β+iγ)=0`, `1 < |γ|`, `K ≥ 32/7`, `A,B,C ≥ 0`, `hβ`, the two
explicit log bounds `hBexp`/`hCexp` (at `R=1/4`), and `hA_pole`.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem region_riemannZeta_R14 {γ β A B C : ℝ}
    (hγ : 1 < |γ|) (hKbig : 32 / 7 ≤ 3 * A + 4 * B + C)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hBexp : 4 * Real.log ((3 + 1 / (|γ| - 1 / 4)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |γ| + 1 / 4)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 / 4))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (1 / 4 / 2) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    (hCexp : 4 * Real.log ((3 + 1 / (|2 * γ| - 1 / 4)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |2 * γ| + 1 / 4)
          / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1 / 4))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (1 / 4 / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hK : 0 < 3 * A + 4 * B + C := by linarith
  have hγ2 : (1 : ℝ) < |2 * γ| := by rw [abs_mul]; simp only [abs_two]; linarith [abs_nonneg γ]
  have hfar : (1 + 1 / (2 * (3 * A + 4 * B + C))) - (1 / 4 : ℝ) / 2
      ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := hfar_holds hK (by nlinarith [hKbig])
  exact region_riemannZeta_explicit_scalar (R := 1 / 4) (R' := |γ|) (R₂ := 1 / 4) (R'₂ := |2 * γ|)
    (by norm_num) (by linarith) (by norm_num) hK hA hB hC hβ hfar (le_refl _) hζρ₀ hBexp hA_pole
    (by norm_num) (by linarith) (by norm_num) (le_refl _) hCexp

#print axioms region_riemannZeta_R14
