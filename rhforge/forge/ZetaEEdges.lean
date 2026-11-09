import ZetaImagAxisBound
import ZetaEntireE

/-!
# S2982 — explicit edge bounds for E(s) = (s−1)ζ(s) on Re s = 1 and Re s = 0 (lane 1: PL edges)

* `norm_zeta_one_add_mul_I_le` — t ≥ 2: ‖ζ(1 + it)‖ ≤ log t + 5  (S2977's Re s = 1 bound at s = 1 + it).
* `norm_zeta_mul_I_le` — t ≥ 2: ‖ζ(it)‖ ≤ √(t/2π)·(log t + 5)  (S2977 squared form).
* `norm_zetaE_one_add_mul_I_le` — t ≥ 2: ‖E(1 + it)‖ ≤ t·(log t + 5).
* `norm_zetaE_mul_I_le` — t ≥ 2: ‖E(it)‖ ≤ (1 + t)·√(t/2π)·(log t + 5).
These are the two boundary inputs of `PhragmenLindelof.vertical_strip` for E on 0 ≤ Re s ≤ 1 (after normalising by a
zero-free factor of size |s|^{(3−s)/2}); the interpolated exponent is (1−σ)/2 + 1 for E, i.e. (1−σ)/2 for ζ.
-/

open Complex

namespace ZetaEEdges

theorem norm_zeta_one_add_mul_I_le {t : ℝ} (ht : 2 ≤ t) :
    ‖riemannZeta (1 + (t : ℂ) * I)‖ ≤ Real.log t + 5 := by
  have hN : 1 ≤ ⌊t⌋₊ := Nat.le_floor (by norm_num; linarith)
  have hs1 : (1 + (t : ℂ) * I).re = 1 := by simp
  have hsne : (1 + (t : ℂ) * I) ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp at this; linarith
  have h := ZetaImagAxisBound.norm_riemannZeta_le_re_one hN hs1 hsne
  have hNx : (⌊t⌋₊ : ℝ) ≤ t := Nat.floor_le (by linarith)
  have hNlt : t < (⌊t⌋₊ : ℝ) + 1 := Nat.lt_floor_add_one t
  have hNpos : (0 : ℝ) < (⌊t⌋₊ : ℝ) := by linarith
  have hlog : Real.log (⌊t⌋₊ : ℝ) ≤ Real.log t := Real.log_le_log hNpos hNx
  have hs1n : ‖(1 + (t : ℂ) * I) - 1‖ = t := by
    rw [show (1 + (t : ℂ) * I) - 1 = (t : ℂ) * I by ring, norm_mul, Complex.norm_I,
      mul_one, Complex.norm_real, Real.norm_of_nonneg (by linarith)]
  have hsn : ‖(1 + (t : ℂ) * I)‖ ≤ 1 + t := by
    calc ‖(1 + (t : ℂ) * I)‖ ≤ ‖(1 : ℂ)‖ + ‖(t : ℂ) * I‖ := norm_add_le _ _
      _ = 1 + t := by
        rw [norm_one, norm_mul, Complex.norm_I, mul_one, Complex.norm_real,
          Real.norm_of_nonneg (by linarith)]
  rw [hs1n] at h
  have h1t : 1 / t ≤ 1 := by rw [div_le_one (by linarith)]; linarith
  have h2 : ‖(1 + (t : ℂ) * I)‖ / (⌊t⌋₊ : ℝ) ≤ 3 := by
    rw [div_le_iff₀ hNpos]; linarith
  linarith

theorem norm_zeta_mul_I_le {t : ℝ} (ht : 2 ≤ t) :
    ‖riemannZeta ((t : ℂ) * I)‖ ≤ Real.sqrt (t / (2 * Real.pi)) * (Real.log t + 5) := by
  have h := ZetaImagAxisBound.norm_zeta_mul_I_sq_le_log ht
  have hlog : 0 ≤ Real.log t + 5 := by have := Real.log_nonneg (show (1:ℝ) ≤ t by linarith); linarith
  have hc : 0 ≤ t / (2 * Real.pi) := by positivity
  rw [← Real.sqrt_sq (norm_nonneg (riemannZeta ((t : ℂ) * I))),
    ← Real.sqrt_sq hlog, ← Real.sqrt_mul hc]
  exact Real.sqrt_le_sqrt h

theorem norm_zetaE_one_add_mul_I_le {t : ℝ} (ht : 2 ≤ t) :
    ‖ZetaEntire.zetaE (1 + (t : ℂ) * I)‖ ≤ t * (Real.log t + 5) := by
  have hne : (1 + (t : ℂ) * I) ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp at this; linarith
  rw [ZetaEntire.zetaE_eq hne, norm_mul,
    show (1 + (t : ℂ) * I) - 1 = (t : ℂ) * I by ring, norm_mul, Complex.norm_I, mul_one,
    Complex.norm_real, Real.norm_of_nonneg (by linarith)]
  exact mul_le_mul_of_nonneg_left (norm_zeta_one_add_mul_I_le ht) (by linarith)

theorem norm_zetaE_mul_I_le {t : ℝ} (ht : 2 ≤ t) :
    ‖ZetaEntire.zetaE ((t : ℂ) * I)‖
      ≤ (1 + t) * (Real.sqrt (t / (2 * Real.pi)) * (Real.log t + 5)) := by
  have hne : ((t : ℂ) * I) ≠ 1 := by
    intro h; have := congrArg Complex.re h; simp at this
  rw [ZetaEntire.zetaE_eq hne, norm_mul]
  have h1 : ‖(t : ℂ) * I - 1‖ ≤ 1 + t := by
    calc ‖(t : ℂ) * I - 1‖ ≤ ‖(t : ℂ) * I‖ + ‖(1 : ℂ)‖ := norm_sub_le _ _
      _ = 1 + t := by
        rw [norm_one, norm_mul, Complex.norm_I, mul_one, Complex.norm_real,
          Real.norm_of_nonneg (by linarith)]; ring
  exact mul_le_mul h1 (norm_zeta_mul_I_le ht) (norm_nonneg _) (by linarith)

end ZetaEEdges

#print axioms ZetaEEdges.norm_zeta_one_add_mul_I_le
#print axioms ZetaEEdges.norm_zeta_mul_I_le
#print axioms ZetaEEdges.norm_zetaE_one_add_mul_I_le
#print axioms ZetaEEdges.norm_zetaE_mul_I_le
