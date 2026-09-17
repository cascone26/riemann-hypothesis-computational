import RegionDischargeIGamma
import ZetaSphereN1
import HBexplicitOfGe
import AbsLeNormCenterSubOne

/-!
# §1449 — discharge the 2iγ-disk sphere+BC hypotheses too (both disks now scalar)

Mirror of §1448 for the growth centre `σK+2iγ`: `Ms₂ := 2 + 1/(‖c₂-1‖-R₂) + (‖c₂‖+R₂)/(c₂.re-R₂)` (§1447),
`Mζ₂ := Ms₂+1`, discharging §1448's `hζsphere₂`/`hsl₂`/`hMζ₂`/`hBexplicit₂` (via §1447+§1444), `hc1₂` from §1446
at `γ→2γ`.  Both disks are now reduced to the two explicit scalar bounds `hByscalar` (iγ) and `hByscalar₂` (2iγ):
the quantitative region `β ≤ 1-1/(14K)` for real ζ now follows from PURELY scalar/geometry data — the last
mechanical wiring before the constant choice.  EXIT0 / std-axioms.
-/

open Complex Metric DirichletCharacter

theorem region_riemannZeta_discharge_both {γ β A B C R R' R₂ R'₂ : ℝ}
    (hR : 0 < R) (hRR' : R < R') (hR1 : R < 1)
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hfar : (1 + 1 / (2 * (3 * A + 4 * B + C))) - R / 2 ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)))
    (hR'γ : R' ≤ |γ|)
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hByscalar : 4 * Real.log ((2 + 1 / (‖(((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) - 1‖ - R)
        + (‖((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)‖ + R)
          / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)).re - R) + 1)
        * (3 + 1 / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)).re - 1)))
        / (R / 2) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    (hR₂ : 0 < R₂) (hRR'₂ : R₂ < R'₂) (hR₂1 : R₂ < 1) (hR'₂γ : R'₂ ≤ |2 * γ|)
    (hByscalar₂ : 4 * Real.log ((2 + 1 / (‖(((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
          + I * ((2 * γ : ℝ) : ℂ)) - 1‖ - R₂)
        + (‖((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)‖ + R₂)
          / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re - R₂) + 1)
        * (3 + 1 / ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ)).re - 1)))
        / (R₂ / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  set c₂ : ℂ := ((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ) with hc₂_def
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  have hcre_eq : c₂.re = 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    rw [hc₂_def]
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  have hcre : 1 < c₂.re := by rw [hcre_eq]; exact hσK
  have hcR : 0 < c₂.re - R₂ := by rw [hcre_eq]; linarith
  have hc1 : R₂ < ‖c₂ - 1‖ :=
    lt_of_lt_of_le (lt_of_lt_of_le hRR'₂ hR'₂γ) (by rw [hc₂_def]; exact abs_le_norm_center_sub_one)
  have hden1 : (0 : ℝ) < ‖c₂ - 1‖ - R₂ := by linarith
  set Ms₂ : ℝ := 2 + 1 / (‖c₂ - 1‖ - R₂) + (‖c₂‖ + R₂) / (c₂.re - R₂) with hMs₂_def
  have hMspos : 0 < Ms₂ + 1 := by rw [hMs₂_def]; positivity
  refine region_riemannZeta_discharge_iγ (Mζ₂ := Ms₂ + 1) (Ms₂ := Ms₂) hR hRR' hR1 hK hA hB hC hβ hfar
    hR'γ hζρ₀ hByscalar hA_pole hR₂ hRR'₂ hR'₂γ ?_ (by linarith) hMspos ?_
  · intro z hz
    exact norm_riemannZeta_le_sphere_N1 hR₂ hcR hc1 hz
  · exact hBexplicit_of_ge hcre hMspos hR₂ hByscalar₂

#print axioms region_riemannZeta_discharge_both
