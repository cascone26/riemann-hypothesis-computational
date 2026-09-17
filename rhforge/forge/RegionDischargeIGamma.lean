import RegionExplicitGeom
import ZetaSphereN1
import HBexplicitOfGe
import AbsLeNormCenterSubOne

/-!
# §1448 — discharge the iγ-disk sphere+BC hypotheses of §1437 (N=1 explicit `Ms`)

Sets `Ms := 2 + 1/(‖c-1‖-R) + (‖c‖+R)/(c.re-R)` (§1447 at `c = σK+iγ`), `Mζ := Ms+1`, and discharges §1437's
`hζsphere` (§1447), `hsl`/`hMζ` (arithmetic), and `hBexplicit` (§1444) — reducing the iγ side to a single scalar
bound `hByscalar : B ≥ 4·log((Ms+1)·(3+1/(σK-1)))/(R/2)`.  `hc1 : R < ‖c-1‖` is derived from `R < R' ≤ |γ| ≤ ‖c-1‖`
(§1446).  The 2iγ disk is passed through unchanged.  EXIT0 / std-axioms.
-/

open Complex Metric DirichletCharacter

theorem region_riemannZeta_discharge_iγ {γ β A B C R R' R₂ R'₂ Mζ₂ Ms₂ : ℝ}
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
    (hR₂ : 0 < R₂) (hRR'₂ : R₂ < R'₂) (hR'₂γ : R'₂ ≤ |2 * γ|)
    (hζsphere₂ : ∀ z ∈ sphere (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ)
        + I * ((2 * γ : ℝ) : ℂ)) R₂, ‖riemannZeta z‖ ≤ Ms₂)
    (hsl₂ : Ms₂ < Mζ₂) (hMζ₂ : 0 < Mζ₂)
    (hBexplicit₂ : 4 * Real.log (Mζ₂ / ‖riemannZeta
      (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * ((2 * γ : ℝ) : ℂ))‖) / (R₂ / 2) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  set c : ℂ := ((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ) with hc_def
  have hσK : (1 : ℝ) < 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    have : (0 : ℝ) < 1 / (2 * (3 * A + 4 * B + C)) := by positivity
    linarith
  have hcre_eq : c.re = 1 + 1 / (2 * (3 * A + 4 * B + C)) := by
    rw [hc_def]
    simp only [Complex.add_re, Complex.ofReal_re, Complex.mul_re, Complex.I_re, Complex.I_im,
      Complex.ofReal_im, zero_mul, mul_zero, sub_zero, add_zero]
  have hcre : 1 < c.re := by rw [hcre_eq]; exact hσK
  have hcR : 0 < c.re - R := by rw [hcre_eq]; linarith
  have hc1 : R < ‖c - 1‖ :=
    lt_of_lt_of_le (lt_of_lt_of_le hRR' hR'γ) (by rw [hc_def]; exact abs_le_norm_center_sub_one)
  have hden1 : (0 : ℝ) < ‖c - 1‖ - R := by linarith
  set Ms : ℝ := 2 + 1 / (‖c - 1‖ - R) + (‖c‖ + R) / (c.re - R) with hMs_def
  have hMspos : 0 < Ms + 1 := by rw [hMs_def]; positivity
  refine region_riemannZeta_explicit_geom (Mζ := Ms + 1) (Ms := Ms) hR hRR' hK hA hB hC hβ hfar hR'γ hζρ₀
    ?_ (by linarith) hMspos ?_ hA_pole hR₂ hRR'₂ hR'₂γ hζsphere₂ hsl₂ hMζ₂ hBexplicit₂
  · intro z hz
    exact norm_riemannZeta_le_sphere_N1 hR hcR hc1 hz
  · exact hBexplicit_of_ge hcre hMspos hR hByscalar

#print axioms region_riemannZeta_discharge_iγ
