import Mathlib.Analysis.SpecialFunctions.Complex.Log
import Mathlib.Analysis.SpecialFunctions.Pow.Complex
import Mathlib.Analysis.SpecialFunctions.Pow.Real
/-!
Lean core of the §812 surgery-rigidity proposition (BARRIER_G_RESULTS.md):
an Euler factor `1 - n^(-s)` over a real base `n ≥ 2` can only vanish on the
imaginary axis `Re s = 0`. Hence multiplicative surgery on ℕ (deleting/adding
integer g-primes) contributes zeros and poles ONLY on `Re s = 0`, and every
sparse-surgery Beurling system inherits ζ's strip-zeros verbatim.
-/

/-- If an Euler-type factor `1 - n^(-s)` vanishes for a real base `n ≥ 2`,
then `s` lies on the imaginary axis. -/
theorem surgery_factor_zero_re_eq_zero {n : ℝ} (hn : 2 ≤ n) {s : ℂ}
    (h : 1 - (n : ℂ) ^ (-s) = 0) : s.re = 0 := by
  have hn0 : (0 : ℝ) < n := by linarith
  have h1 : (n : ℂ) ^ (-s) = 1 := by
    have := sub_eq_zero.mp h
    exact this.symm
  have hnorm : ‖(n : ℂ) ^ (-s)‖ = n ^ (-s).re := by
    simpa using Complex.norm_cpow_eq_rpow_re_of_pos hn0 (-s)
  have hone : n ^ (-s).re = 1 := by
    rw [← hnorm, h1, norm_one]
  have hb : (1 : ℝ) < n := by linarith
  have : (-s).re = 0 := by
    by_contra hne
    rcases lt_or_gt_of_ne hne with hlt | hgt
    · have := Real.rpow_lt_one_of_one_lt_of_neg hb hlt
      rw [hone] at this; exact lt_irrefl 1 this
    · have := Real.one_lt_rpow_iff_of_pos hn0 |>.mpr (Or.inl ⟨hb, hgt⟩)
      rw [hone] at this; exact lt_irrefl 1 this
  simpa using neg_eq_zero.mp (by simpa using this)

/-- Contrapositive form used by the rigidity proposition: off the imaginary
axis, surgery factors never vanish — so in any strip `Re s > σ₀ > 0` the
surgered zeta has exactly ζ's zeros. -/
theorem surgery_factor_ne_zero_of_re_ne {n : ℝ} (hn : 2 ≤ n) {s : ℂ}
    (hs : s.re ≠ 0) : 1 - (n : ℂ) ^ (-s) ≠ 0 :=
  fun h => hs (surgery_factor_zero_re_eq_zero hn h)

#print axioms surgery_factor_zero_re_eq_zero
#print axioms surgery_factor_ne_zero_of_re_ne
