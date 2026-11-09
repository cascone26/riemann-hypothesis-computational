import ZetaEEdges

/-!
# S2991 — convexity edges for BOTH signs of t: |t| ≥ 2 (lane 1)

S2982 had t ≥ 2.  Negative heights:
* right edge: ζ(1 + it) with t ≤ −2 is ζ(1 − i|t|) — already S2977 `norm_zeta_one_sub_mul_I_le`;
* left edge: ζ(it) with t ≤ −2 equals conj ζ(i|t|) (Mathlib `riemannZeta_conj`), same norm.
* `norm_zeta_one_add_mul_I_le_abs` — 2 ≤ |t| ⇒ ‖ζ(1 + it)‖ ≤ log|t| + 5.
* `norm_zeta_mul_I_le_abs` — 2 ≤ |t| ⇒ ‖ζ(it)‖ ≤ √(|t|/2π)·(log|t| + 5).
Remaining edge gap for three-lines: |t| < 2 (compact, continuity of E = (s−1)ζ(s) on the two segments).
-/

open Complex

namespace ZetaEdgesUniform

theorem norm_zeta_one_add_mul_I_le_abs {t : ℝ} (ht : 2 ≤ |t|) :
    ‖riemannZeta (1 + (t : ℂ) * I)‖ ≤ Real.log |t| + 5 := by
  rcases le_or_gt 0 t with h | h
  · rw [abs_of_nonneg h] at ht ⊢
    exact ZetaEEdges.norm_zeta_one_add_mul_I_le ht
  · rw [abs_of_neg h] at ht ⊢
    have e : (1 + (t : ℂ) * I) = 1 - ((-t : ℝ) : ℂ) * I := by push_cast; ring
    rw [e]
    exact ZetaImagAxisBound.norm_zeta_one_sub_mul_I_le ht

theorem norm_zeta_mul_I_le_abs {t : ℝ} (ht : 2 ≤ |t|) :
    ‖riemannZeta ((t : ℂ) * I)‖ ≤ Real.sqrt (|t| / (2 * Real.pi)) * (Real.log |t| + 5) := by
  rcases le_or_gt 0 t with h | h
  · rw [abs_of_nonneg h] at ht ⊢
    exact ZetaEEdges.norm_zeta_mul_I_le ht
  · rw [abs_of_neg h] at ht ⊢
    have e : ((t : ℂ) * I) = (starRingEnd ℂ) (((-t : ℝ) : ℂ) * I) := by
      simp [map_mul, Complex.conj_ofReal, Complex.conj_I]
    rw [e, riemannZeta_conj, Complex.norm_conj]
    exact ZetaEEdges.norm_zeta_mul_I_le ht

end ZetaEdgesUniform

#print axioms ZetaEdgesUniform.norm_zeta_one_add_mul_I_le_abs
#print axioms ZetaEdgesUniform.norm_zeta_mul_I_le_abs
