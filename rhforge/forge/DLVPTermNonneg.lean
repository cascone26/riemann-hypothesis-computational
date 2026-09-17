import Mathlib.NumberTheory.LSeries.Dirichlet
import ComplexDLVP

/-!
# §1357 — termwise DLVP `3-4-1` nonnegativity for the von Mangoldt series

For every `n`, the `3-4-1` combination of the L-series *terms* of `↗Λ` at the three points
`σ, σ+iγ, σ+2iγ` is nonnegative.  For `n ≥ 1`, `term ↗Λ (σ+ib) n = a · (↑n)^{-ib}` with `a = Λ(n)·n^{-σ} ≥ 0`
real, so the combination is `a·(3 + 4 Re z + Re(z²))`, `z = (↑n)^{-iγ}`, `‖z‖ = 1`, `z² = (↑n)^{-2iγ}` —
nonnegative by `complex_dlvp_factor_nonneg` (§1116).  Summing over `n` (next brick) gives the DLVP positivity
feeding §1355's `hdlvp`.  EXIT0 / std-axioms.
-/

open Complex ArithmeticFunction LSeries
open scoped LSeries.notation

theorem dlvp_term_combination_nonneg {σ γ : ℝ} (n : ℕ) :
    0 ≤ 3 * (term (↗Λ) ((σ : ℂ)) n).re
      + 4 * (term (↗Λ) ((σ : ℂ) + I * (γ : ℝ)) n).re
      + (term (↗Λ) ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ)) n).re := by
  rcases eq_or_ne n 0 with hn0 | hn0
  · subst hn0; simp [LSeries.term]
  · have hn : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn0
    have hnpos : (0 : ℝ) < (n : ℝ) := by positivity
    set a : ℝ := (Λ n) * (n : ℝ) ^ (-σ) with ha_def
    set z : ℂ := (n : ℂ) ^ (-(I * (γ : ℝ))) with hz_def
    -- key: Re(term at σ+ib) = a · Re((↑n)^{-ib})
    have key : ∀ b : ℝ, (term (↗Λ) ((σ : ℂ) + I * (b : ℝ)) n).re
        = a * ((n : ℂ) ^ (-(I * (b : ℝ)))).re := by
      intro b
      rw [LSeries.term_of_ne_zero hn0]
      have hsplit : (n : ℂ) ^ ((σ : ℂ) + I * (b : ℝ))
          = (n : ℂ) ^ (σ : ℂ) * (n : ℂ) ^ (I * (b : ℝ)) := Complex.cpow_add _ _ hn
      rw [hsplit]
      have hreal : (n : ℂ) ^ (-(σ : ℂ)) = (((n : ℝ) ^ (-σ) : ℝ) : ℂ) := by
        rw [show ((n : ℂ)) = (((n : ℝ)) : ℂ) from by push_cast; ring,
          show (-(σ : ℂ)) = (((-σ) : ℝ) : ℂ) from by push_cast; ring,
          ← Complex.ofReal_cpow hnpos.le]
      have hcombine : ((Λ n : ℝ) : ℂ) / ((n : ℂ) ^ (σ : ℂ) * (n : ℂ) ^ (I * (b : ℝ)))
          = ((a : ℝ) : ℂ) * (n : ℂ) ^ (-(I * (b : ℝ))) := by
        rw [ha_def]
        rw [div_eq_mul_inv, mul_inv, ← Complex.cpow_neg, ← Complex.cpow_neg, hreal]
        push_cast [ha_def]
        ring
      rw [hcombine, Complex.re_ofReal_mul]
    -- the three points
    have e0 : (term (↗Λ) ((σ : ℂ)) n).re = a := by
      have h := key 0
      simp only [Complex.ofReal_zero, mul_zero, add_zero, neg_zero, Complex.cpow_zero,
        Complex.one_re, mul_one] at h
      exact h
    have e1 : (term (↗Λ) ((σ : ℂ) + I * (γ : ℝ)) n).re = a * z.re := key γ
    have hz2 : z ^ 2 = (n : ℂ) ^ (-(I * ((2 * γ : ℝ) : ℂ))) := by
      rw [sq, hz_def, ← Complex.cpow_add _ _ hn]
      congr 1; push_cast; ring
    have e2 : (term (↗Λ) ((σ : ℂ) + I * ((2 * γ : ℝ) : ℂ)) n).re = a * (z ^ 2).re := by
      have h := key (2 * γ)
      rw [hz2]
      exact h
    -- assemble
    have ha : 0 ≤ a := by
      rw [ha_def]; exact mul_nonneg vonMangoldt_nonneg (Real.rpow_nonneg hnpos.le _)
    have hznorm : ‖z‖ = 1 := by
      rw [hz_def, show (n : ℂ) = (((n : ℝ)) : ℂ) by push_cast; ring,
        Complex.norm_cpow_eq_rpow_re_of_pos hnpos]
      simp
    have hfac : 0 ≤ 3 + 4 * z.re + (z ^ 2).re := complex_dlvp_factor_nonneg z hznorm
    rw [e0, e1, e2]
    have : 3 * a + 4 * (a * z.re) + a * (z ^ 2).re = a * (3 + 4 * z.re + (z ^ 2).re) := by ring
    rw [this]
    exact mul_nonneg ha hfac

#print axioms dlvp_term_combination_nonneg
