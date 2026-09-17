import RegionR14
import HBexpRhsLe

/-!
# §1457 — quantitative ζ region from CLEAN linear-in-γ log bounds

Discharges §1454's giant `hBexp`/`hCexp` from the clean bounds `B ≥ 32·log((3|γ|+10)(3+2K))` and
`C ≥ 32·log((3|2γ|+10)(3+2K))` (via §1456 + `le_trans`).  The quantitative DLVP region `β ≤ 1-1/(14K)` for real ζ
now follows from: `ζ(β+iγ)=0`, `1<|γ|`, `K ≥ 32/7`, `A,B,C ≥ 0`, `hβ`, `hA_pole`, and these two linear-in-`|γ|`
log bounds — the shape the §1445 threshold closes.  EXIT0 / std-axioms.
-/

open Complex DirichletCharacter

theorem region_riemannZeta_R14_clean {γ β A B C : ℝ}
    (hγ : 1 < |γ|) (hKbig : 32 / 7 ≤ 3 * A + 4 * B + C)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hζρ₀ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hBclean : 32 * Real.log ((3 * |γ| + 10) * (3 + 2 * (3 * A + 4 * B + C))) ≤ B)
    (hA_pole : ‖logDeriv (LFunctionTrivChar₁ 1)
      ((((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ)) : ℂ))‖ ≤ A)
    (hCclean : 32 * Real.log ((3 * |2 * γ| + 10) * (3 + 2 * (3 * A + 4 * B + C))) ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  have hKhalf : (1 : ℝ) / 2 ≤ 3 * A + 4 * B + C := by linarith
  have hγ2 : (1 : ℝ) < |2 * γ| := by rw [abs_mul]; simp only [abs_two]; linarith [abs_nonneg γ]
  exact region_riemannZeta_R14 hγ hKbig hA hB hC hβ hζρ₀
    (le_trans (hBexp_rhs_le hKhalf hγ) hBclean) hA_pole
    (le_trans (hBexp_rhs_le hKhalf hγ2) hCclean)

#print axioms region_riemannZeta_R14_clean
