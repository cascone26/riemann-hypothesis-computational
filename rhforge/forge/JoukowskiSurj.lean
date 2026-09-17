import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Tactic

/-!
# §1782 — Joukowski surjectivity: every w ∈ [-2,2] is z + z⁻¹ for some z on the unit circle

Companion to §1781 (`JoukowskiOnCircle`).  §1781 gave `z + z⁻¹ ∈ [-2,2] real ⟺ |z|=1`; this brick gives
the SURJECTIVITY of the other direction: for every real `w ∈ [-2,2]` there is `z` with `|z|=1` and
`z + z⁻¹ = w`, via the explicit `z = (w + i√(4-w²))/2`.  Together §1781+§1782 establish that `z ↦ z+z⁻¹`
maps the unit circle ONTO the segment `[-2,2]` — the exact correspondence the Chebyshev/Schur-Cohn
route uses to translate "zeros on the circle" into "roots real in [-2,2]".  EXIT0 / std-axioms.
-/

open Complex

/-- **Joukowski surjectivity onto [-2,2].**  For real `w` with `-2 ≤ w ≤ 2`, there is a complex `z` on the
unit circle (`normSq z = 1`) with `z + z⁻¹ = w`. -/
theorem joukowski_surj (w : ℝ) (hw : -2 ≤ w) (hw2 : w ≤ 2) :
    ∃ z : ℂ, normSq z = 1 ∧ z + z⁻¹ = (w : ℂ) := by
  have h4 : (0 : ℝ) ≤ 4 - w ^ 2 := by nlinarith
  have hsqrt : Real.sqrt (4 - w ^ 2) ^ 2 = 4 - w ^ 2 := Real.sq_sqrt h4
  refine ⟨⟨w / 2, Real.sqrt (4 - w ^ 2) / 2⟩, ?_, ?_⟩
  · -- normSq = (w/2)² + (√(4-w²)/2)² = w²/4 + (4-w²)/4 = 1
    rw [normSq_apply]
    show w / 2 * (w / 2) + Real.sqrt (4 - w ^ 2) / 2 * (Real.sqrt (4 - w ^ 2) / 2) = 1
    nlinarith [hsqrt]
  · -- z + z⁻¹ = z + conj z (since |z|=1) = 2·Re(z) = w
    have hns : normSq (⟨w / 2, Real.sqrt (4 - w ^ 2) / 2⟩ : ℂ) = 1 := by
      rw [normSq_apply]
      show w / 2 * (w / 2) + Real.sqrt (4 - w ^ 2) / 2 * (Real.sqrt (4 - w ^ 2) / 2) = 1
      nlinarith [hsqrt]
    apply Complex.ext
    · rw [add_re, inv_re, hns]
      show w / 2 + w / 2 / 1 = (w : ℂ).re
      rw [Complex.ofReal_re]; ring
    · rw [add_im, inv_im, hns]
      show Real.sqrt (4 - w ^ 2) / 2 + -(Real.sqrt (4 - w ^ 2) / 2) / 1 = (w : ℂ).im
      rw [Complex.ofReal_im]; ring
