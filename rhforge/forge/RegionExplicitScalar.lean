import RegionDischargeBoth
import HByscalarOfExplicit

/-!
# §1453 — quantitative ζ region from CLEAN explicit scalar bounds (no norms)

Discharges §1449's two `hByscalar` conditions via §1452, so the region `β ≤ 1-1/(14K)` for real ζ now follows from
the norm-free explicit inequalities

  `B ≥ 4·log((3+1/(|γ|-R)+(σK+|γ|+R)/(σK-R))·(3+1/(σK-1)))/(R/2)`   (iγ)
  `C ≥ 4·log((3+1/(|2γ|-R₂)+(σK+|2γ|+R₂)/(σK-R₂))·(3+1/(σK-1)))/(R₂/2)`   (2iγ)

together with `ζ(β+iγ)=0`, the constant/geometry scalars, and `hA_pole`.  `σK = 1+1/(2K)`, `K = 3A+4B+C`.  This is
the fully explicit scalar form — the last brick before the numeric constant choice.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem region_riemannZeta_explicit_scalar {γ β A B C R R' R₂ R'₂ : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hR1 : R < 1)
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hfar : (1 + 1 / (2 * (3 * A + 4 * B + C))) - R / 2 ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)))
    (hR'γ : R' ≤ |γ|)
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hBexp : 4 * Real.log ((3 + 1 / (|γ| - R)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |γ| + R) / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - R))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (R / 2) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    (hR₂ : 0 < R₂) (hRR'₂ : R₂ < R'₂) (hR₂1 : R₂ < 1) (hR'₂γ : R'₂ ≤ |2 * γ|)
    (hCexp : 4 * Real.log ((3 + 1 / (|2 * γ| - R₂)
        + ((1 + 1 / (2 * (3 * A + 4 * B + C))) + |2 * γ| + R₂) / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - R₂))
        * (3 + 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1))) / (R₂ / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  have hRγ : R < |γ| := lt_of_lt_of_le hRR' hR'γ
  have hRγ₂ : R₂ < |2 * γ| := lt_of_lt_of_le hRR'₂ hR'₂γ
  exact region_riemannZeta_discharge_both hR hRR' hR1 hK hA hB hC hβ hfar hR'γ hζρ₀
    (hByscalar_of_explicit hσK hR hRγ (by linarith) hBexp) hA_pole hR₂ hRR'₂ hR₂1 hR'₂γ
    (hByscalar_of_explicit hσK hR₂ hRγ₂ (by linarith) hCexp)

#print axioms region_riemannZeta_explicit_scalar
